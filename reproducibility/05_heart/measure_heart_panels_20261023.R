# R measures actual heatmap/facet viewports; Python only audits those rectangles.
measure_heart_panels <- function(grob, prefix, width, height, body_pattern,
                                n_row, n_col, row_groups, column_groups) {
  probe <- paste0(prefix, ".layout-probe.pdf")
  grDevices::cairo_pdf(probe, width = width, height = height, family = "Arial")
  grid::grid.newpage()
  grid::grid.draw(grob)
  grid::grid.force()
  tree <- grid::grid.ls(viewports = TRUE, grobs = FALSE, print = FALSE)
  readr::write_tsv(data.frame(name = tree$name, path = tree$vpPath), paste0(prefix, ".viewports.tsv"))
  ix <- which(grepl(body_pattern, tree$name))
  stopifnot(length(ix) == n_row * n_col)
  bounds <- lapply(ix, function(i) {
    grid::upViewport(0)
    path <- sub("^ROOT::", "", paste(tree$vpPath[i], tree$name[i], sep = "::"))
    grid::downViewport(path)
    lo <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc"), valueOnly = TRUE)
    hi <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(1, "npc"), valueOnly = TRUE)
    unname(c(lo$x, lo$y, hi$x, hi$y) * 72)
  })
  grid::upViewport(0)
  grDevices::dev.off()
  by_y <- order(vapply(bounds, function(b) -(b[2] + b[4]) / 2, numeric(1)))
  rows <- split(by_y, rep(seq_len(n_row), each = n_col))
  ord <- unlist(lapply(rows, function(z) z[order(vapply(bounds[z], function(b) b[1], numeric(1)))]))
  panels <- lapply(seq_along(ord), function(i) list(id = paste0("p", i), bbox_pt = bounds[[ord[i]]]))
  manifest <- list(
    schema_version = 1L, backend = "r-grid-actual-viewports",
    figure = list(width_pt = width * 72, height_pt = height * 72),
    panels = panels, row_groups = unname(row_groups), column_groups = unname(column_groups)
  )
  jsonlite::write_json(manifest, paste0(prefix, ".alignment-layout.json"), auto_unbox = TRUE, pretty = TRUE)
  result <- system2(Sys.which("python3"), c(
    shQuote("external/plot_helpers/audit_panel_alignment.py"),
    shQuote(paste0(prefix, ".alignment-layout.json")),
    "--json-out", shQuote(paste0(prefix, ".alignment.json")),
    "--overlay-svg", shQuote(paste0(prefix, ".alignment.svg")),
    "--tolerance-pt", "1.5", "--gutter-tolerance-pt", "1.5", "--strict"
  ))
  stopifnot(result == 0L)
  unlink(probe)
  invisible(manifest)
}
