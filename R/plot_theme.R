#' scRhythmBayes ggplot theme
#'
#' This function provides the default ggplot2 theme used by scRhythmBayes
#' plotting functions.
#'
#' @param aspect.ratio A numeric value specifying the panel aspect ratio.
#' Defaults to \code{NULL}.
#' @param base_size A numeric value specifying the base font size. Defaults to
#' \code{12}.
#' @param base_family A character vector specifying the base font family.
#' Defaults to \code{""}.
#' @param ... Additional arguments passed to \code{\link[ggplot2]{theme}}.
#'
#' @return A ggplot2 theme object.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
#'   geom_point()
#'
#' p + theme_srb()
#' }
#'
#' @importFrom ggplot2 theme theme_bw element_blank element_text element_rect
#' @importFrom ggplot2 element_line margin
#' @export
theme_srb <- function(
    aspect.ratio = NULL,
    base_size = 12,
    base_family = "",
    ...
) {
  text_scale <- base_size / 12

  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      aspect.ratio = aspect.ratio,
      text = ggplot2::element_text(
        size = 12 * text_scale,
        color = "black",
        family = base_family
      ),
      plot.title = ggplot2::element_text(
        size = 14 * text_scale,
        color = "black",
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 6)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 12 * text_scale,
        color = "black",
        hjust = 0,
        margin = ggplot2::margin(b = 4)
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = "white"
      ),
      panel.background = ggplot2::element_rect(
        fill = "white",
        color = "white"
      ),
      panel.border = ggplot2::element_rect(
        fill = NA,
        color = "black",
        linewidth = 0.6
      ),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(
        size = 12 * text_scale,
        color = "black"
      ),
      axis.text = ggplot2::element_text(
        size = 11 * text_scale,
        color = "black"
      ),
      axis.ticks = ggplot2::element_line(
        color = "black",
        linewidth = 0.4
      ),
      strip.text = ggplot2::element_text(
        size = 12 * text_scale,
        color = "black",
        face = "bold",
        margin = ggplot2::margin(3, 3, 3, 3)
      ),
      strip.background = ggplot2::element_rect(
        fill = "grey95",
        color = "black",
        linewidth = 0.5
      ),
      legend.title = ggplot2::element_text(
        size = 11 * text_scale,
        color = "black",
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = 10 * text_scale,
        color = "black"
      ),
      legend.key = ggplot2::element_rect(
        fill = "transparent",
        color = "transparent"
      ),
      legend.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      ...
    )
}


#' scRhythmBayes void ggplot theme
#'
#' This function provides a blank ggplot2 theme used by circular or
#' annotation-only scRhythmBayes plotting functions.
#'
#' @param base_size A numeric value specifying the base font size. Defaults to
#' \code{12}.
#' @param base_family A character vector specifying the base font family.
#' Defaults to \code{""}.
#' @param ... Additional arguments passed to \code{\link[ggplot2]{theme}}.
#'
#' @return A ggplot2 theme object.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
#'   geom_point()
#'
#' p + theme_srb_void()
#' }
#'
#' @importFrom ggplot2 theme theme_void element_blank element_text element_rect
#' @importFrom ggplot2 margin
#' @export
theme_srb_void <- function(
    base_size = 12,
    base_family = "",
    ...
) {
  text_scale <- base_size / 12

  ggplot2::theme_void(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(
        size = 12 * text_scale,
        color = "black",
        family = base_family
      ),
      plot.title = ggplot2::element_text(
        size = 16 * text_scale,
        color = "black",
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 6)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 12 * text_scale,
        color = "black",
        hjust = 0.5,
        margin = ggplot2::margin(b = 4)
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = "white"
      ),
      panel.background = ggplot2::element_rect(
        fill = "white",
        color = "white"
      ),
      strip.text = ggplot2::element_text(
        size = 13 * text_scale,
        color = "black",
        face = "bold",
        margin = ggplot2::margin(3, 3, 3, 3)
      ),
      strip.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(
        size = 11 * text_scale,
        color = "black",
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = 10 * text_scale,
        color = "black"
      ),
      legend.key = ggplot2::element_rect(
        fill = "transparent",
        color = "transparent"
      ),
      legend.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(20, 30, 20, 20),
      ...
    )
}


#' scRhythmBayes phase palette
#'
#' This function returns a soft circadian color palette for Zeitgeber time
#' groups. The default palette is anchored to the six common sampling times
#' ZT02, ZT06, ZT10, ZT14, ZT18, and ZT22, and interpolated to 24 integer
#' hours.
#'
#' @param zt A numeric vector of Zeitgeber time values. If \code{NULL}, colors
#' for integer hours 0 to 23 are returned. Defaults to \code{NULL}.
#' @param matched Logical; if \code{TRUE}, colors are matched to \code{zt}. If
#' \code{FALSE}, the full 24-hour palette is returned. Defaults to \code{TRUE}.
#' @param names_as A character vector specifying the name format of returned
#' colors. Must be one of \code{"key"} or \code{"label"}. \code{"key"} uses
#' numeric hour names such as \code{"2"}; \code{"label"} uses names such as
#' \code{"ZT02"}. Defaults to \code{"key"}.
#'
#' @return A named character vector of colors.
#'
#' @examples
#' palette_srb_phase()
#' palette_srb_phase(c(2, 6, 10, 14, 18, 22))
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom stats setNames
#' @export
palette_srb_phase <- function(
    zt = NULL,
    matched = TRUE,
    names_as = c("key", "label")
) {
  names_as <- match.arg(names_as)

  anchor_cols <- c(
    "#D97B93",
    "#F08A72",
    "#F3B36A",
    "#A9D17A",
    "#7FC4A5",
    "#7E79BF"
  )

  pal_24 <- grDevices::colorRampPalette(
    c(anchor_cols, anchor_cols[1])
  )(25)[1:24]

  names(pal_24) <- as.character(0:23)

  if (!isTRUE(matched) || is.null(zt)) {
    if (names_as == "label") {
      names(pal_24) <- .srb_zt_label(as.numeric(names(pal_24)))
    }

    return(pal_24)
  }

  zt <- as.numeric(zt)
  zt <- zt[is.finite(zt)]

  if (length(zt) == 0) {
    stop("zt must contain at least one finite value.")
  }

  zt_key <- .srb_zt_key(zt)
  zt_lab <- .srb_zt_label(zt)

  if (all(abs(zt - round(zt)) < 1e-8)) {
    out <- pal_24[zt_key]
  } else {
    out <- grDevices::colorRampPalette(anchor_cols)(length(zt_key))
  }

  names(out) <- if (names_as == "label") {
    zt_lab
  } else {
    zt_key
  }

  out
}


.srb_zt_key <- function(x) {
  x <- as.numeric(x)

  ifelse(
    abs(x - round(x)) < 1e-8,
    as.character(as.integer(round(x)) %% 24),
    as.character(x %% 24)
  )
}


.srb_zt_label <- function(x) {
  x <- as.numeric(x)

  ifelse(
    abs(x - round(x)) < 1e-8,
    paste0("ZT", sprintf("%02d", as.integer(round(x)) %% 24)),
    paste0("ZT", x)
  )
}


.srb_match_zt_colors <- function(
    zt_levels,
    zt_colors = NULL
) {
  zt_levels <- as.numeric(zt_levels)
  zt_levels <- zt_levels[is.finite(zt_levels)]

  if (length(zt_levels) == 0) {
    stop("zt_levels must contain at least one finite value.")
  }

  zt_keys <- .srb_zt_key(zt_levels)
  zt_labels <- .srb_zt_label(zt_levels)

  if (is.null(zt_colors)) {
    zt_colors <- palette_srb_phase(
      zt = zt_levels,
      matched = TRUE,
      names_as = "key"
    )
  } else {
    if (is.null(names(zt_colors))) {
      if (length(zt_colors) != length(zt_keys)) {
        stop("zt_colors must be named or have the same length as zt_levels.")
      }

      names(zt_colors) <- zt_keys
    } else {
      if (all(zt_keys %in% names(zt_colors))) {
        zt_colors <- zt_colors[zt_keys]
      } else if (all(zt_labels %in% names(zt_colors))) {
        zt_colors <- zt_colors[zt_labels]
        names(zt_colors) <- zt_keys
      } else {
        stop("zt_colors must contain colors for all ZT levels.")
      }
    }
  }

  zt_colors
}
.srb_discrete_colors <- function(x) {
  x <- as.character(x)
  x <- unique(x)

  base_cols <- c(
    "#E64B35",
    "#4DBBD5",
    "#00A087",
    "#91D1C2",
    "#3C5488",
    "#F39B7F",
    "#8491B4",
    "#DC0000",
    "#7E6148",
    "#B09C85"
  )

  if (length(x) <= length(base_cols)) {
    cols <- base_cols[seq_along(x)]
  } else {
    cols <- grDevices::colorRampPalette(base_cols)(length(x))
  }

  stats::setNames(cols, x)
}
.srb_group_colors <- function(x) {
  x <- as.character(x)
  x <- unique(x)

  base_cols <- c(
    "#A8C7E6",
    "#F2A7A7",
    "#B7DDB2",
    "#F2D49B",
    "#CBB7E8",
    "#B8DDE1",
    "#E7C2A2",
    "#D6D6D6"
  )

  if (length(x) <= length(base_cols)) {
    cols <- base_cols[seq_along(x)]
  } else {
    cols <- grDevices::colorRampPalette(base_cols)(length(x))
  }

  stats::setNames(cols, x)
}
