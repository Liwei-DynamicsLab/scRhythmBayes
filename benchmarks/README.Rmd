# Benchmarks

This directory contains the benchmark analyses used to evaluate scRhythmBayes.

The benchmarks assess the main components of the workflow: cell-level circadian phase inference, gene-level rhythmic parameter estimation and differential rhythm-state module assignment. Simulated single-cell circadian datasets with known ground truth were used to evaluate phase recovery, rhythmic parameter accuracy and module classification.

The benchmark suite includes:

- seed-gene number robustness;
- no-ZT-prior ablation;
- sampling time-point number benchmark;
- cell-number benchmark;
- module assignment evaluation;
- runtime evaluation.

These analyses support the robustness tests described in the manuscript, including the effects of seed-gene support, Zeitgeber-time anchoring, temporal sampling density and cell number on scRhythmBayes performance.

Large generated outputs, including simulated objects and fitted model results, are not stored in this repository.