#!/usr/bin/env python3
"""Profile wall time, aggregate process-tree RSS, and effective CPU cores."""

import argparse
import csv
import datetime as dt
import json
import os
import platform
import shlex
import subprocess
import sys
import time
from pathlib import Path


def process_snapshot():
    proc = subprocess.run(
        ["ps", "-axo", "pid=,ppid=,pgid=,rss=,%cpu=,command="],
        check=False,
        capture_output=True,
        text=True,
    )
    rows = []
    for line in proc.stdout.splitlines():
        fields = line.split(maxsplit=5)
        if len(fields) < 5:
            continue
        try:
            rows.append(
                {
                    "pid": int(fields[0]),
                    "ppid": int(fields[1]),
                    "pgid": int(fields[2]),
                    "rss_kb": float(fields[3]),
                    "cpu_pct": float(fields[4]),
                    "command": fields[5] if len(fields) == 6 else "",
                }
            )
        except ValueError:
            continue
    return rows


def descendant_pids(rows, root_pid):
    children = {}
    for row in rows:
        children.setdefault(row["ppid"], []).append(row["pid"])
    found = {root_pid}
    frontier = [root_pid]
    while frontier:
        parent = frontier.pop()
        for child in children.get(parent, []):
            if child not in found:
                found.add(child)
                frontier.append(child)
    # future/parallelly PSOCK workers can be re-parented to launchd on macOS.
    worker_tag = f"parallelly.parent={root_pid}."
    found.update(row["pid"] for row in rows if worker_tag in row.get("command", ""))
    return found


def hardware_info():
    info = {
        "hostname": platform.node(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "logical_cpu_count": os.cpu_count(),
    }
    for key in ("machdep.cpu.brand_string", "hw.memsize", "hw.physicalcpu", "hw.logicalcpu"):
        proc = subprocess.run(["sysctl", "-n", key], check=False, capture_output=True, text=True)
        if proc.returncode == 0:
            info[key] = proc.stdout.strip()
    return info


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--interval", type=float, default=0.25)
    parser.add_argument(
        "--configured-cores",
        type=float,
        default=None,
        help="Explicit algorithm/thread cap; omit when the source does not constrain it",
    )
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command and args.command[0] == "--" else args.command
    if not command:
        parser.error("a command is required after --")

    out_path = Path(args.out)
    samples_path = Path(args.samples)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    samples_path.parent.mkdir(parents=True, exist_ok=True)

    started_utc = dt.datetime.now(dt.timezone.utc).isoformat()
    start = time.perf_counter()
    child = subprocess.Popen(command, start_new_session=True)
    samples = []
    tracked_pids = {child.pid}
    while child.poll() is None:
        rows = process_snapshot()
        alive_pids = {row["pid"] for row in rows}
        # Preserve processes previously observed in the process tree. R PSOCK
        # workers may be re-parented to launchd while the root R session is
        # still alive; dropping them would undercount both RSS and CPU usage.
        tracked_pids.intersection_update(alive_pids)
        tracked_pids.update(descendant_pids(rows, child.pid))
        # The benchmark command starts in a new process group.  R's PSOCK
        # workers can be re-parented to launchd before the next 0.25-s sample,
        # but they retain this process-group ID.  Selecting by PGID therefore
        # captures the complete benchmark even after re-parenting.
        selected = [
            row
            for row in rows
            if row["pid"] in tracked_pids or row["pgid"] == child.pid
        ]
        samples.append(
            {
                "elapsed_seconds": time.perf_counter() - start,
                "process_count": len(selected),
                "aggregate_rss_mb": sum(row["rss_kb"] for row in selected) / 1024.0,
                "effective_cpu_cores": sum(row["cpu_pct"] for row in selected) / 100.0,
            }
        )
        time.sleep(args.interval)
    exit_code = child.wait()
    wall_seconds = time.perf_counter() - start
    finished_utc = dt.datetime.now(dt.timezone.utc).isoformat()

    peak_rss = max((x["aggregate_rss_mb"] for x in samples), default=0.0)
    peak_cores = max((x["effective_cpu_cores"] for x in samples), default=0.0)
    mean_cores = (
        sum(x["effective_cpu_cores"] for x in samples) / len(samples) if samples else 0.0
    )
    summary = {
        "label": args.label,
        "status": "complete" if exit_code == 0 else "failed",
        "exit_code": exit_code,
        "started_utc": started_utc,
        "finished_utc": finished_utc,
        "wall_seconds": wall_seconds,
        "wall_minutes": wall_seconds / 60.0,
        "peak_process_tree_rss_mb": peak_rss,
        "peak_process_tree_rss_gb": peak_rss / 1024.0,
        "configured_cores": args.configured_cores,
        "mean_effective_cpu_cores": mean_cores,
        "peak_effective_cpu_cores": peak_cores,
        "sample_interval_seconds": args.interval,
        "n_samples": len(samples),
        "command": shlex.join(command),
        "rss_definition": "sum of RSS across the isolated benchmark process group at each sample; shared pages may be counted in more than one process",
        "cpu_definition": "sum of ps percent-CPU across the isolated benchmark process group divided by 100",
        "hardware": hardware_info(),
    }
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
    with samples_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(samples[0].keys()) if samples else ["elapsed_seconds", "process_count", "aggregate_rss_mb", "effective_cpu_cores"])
        writer.writeheader()
        writer.writerows(samples)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
