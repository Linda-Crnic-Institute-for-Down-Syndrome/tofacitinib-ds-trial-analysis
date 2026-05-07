####################################################################
# Common helper functions                                          #
# Script authors: Matthew Galbraith & Micah Donovan                #
# version: 04_28_2026                                              #
####################################################################

# Setting and modifying default theme for ggplot -----
theme_set(theme_gray(base_size = 12, base_family = "Arial") +
            theme(
              panel.border = element_rect(colour="black", fill = "transparent"), 
              plot.title = element_text(face="bold", hjust = 0),
              axis.text = element_text(color="black", size = 14), 
              axis.text.x = element_text(angle = 0, hjust = 0.5),
              panel.background = element_blank(),
              panel.grid = element_blank(),
              plot.background = element_blank(),
              strip.background = element_blank(), # facet label borders
              legend.key=element_blank(), legend.background=element_blank() # remove grey bg from legend
            ))

# Custom function for heatmaps using tidyheatmap    ----
tidy_hm <- function(
    hm_dat,
    row,
    col,
    value,
    palette = "RdBu",
    palette_rev = NULL,
    scale_type = c("diverging", "continuous"),
    n_colors = 11,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    row_dend = TRUE,
    col_dend = TRUE,
    asterisk = FALSE,
    sig_val = NULL,
    title = "Heatmap title",
    subtitle = "",
    use_quantile = TRUE,
    quantile_val = 0.99,
    hm_max = NULL,
    col_size = 3,
    row_size = 3,
    row_font = 10,
    col_font = 10,
    rotate_text = FALSE
) {
  scale_type <- match.arg(scale_type)
  
  if (is.null(palette_rev)) {
    palette_rev <- if (scale_type == "diverging") TRUE else FALSE
  }
  
  row_var <- rlang::enquo(row)
  col_var <- rlang::enquo(col)
  val_var <- rlang::enquo(value)
  sig_var <- rlang::enquo(sig_val)
  
  vals <- dplyr::pull(hm_dat, !!val_var)
  vals_non_na <- vals[!is.na(vals)]
  
  if (length(vals_non_na) == 0) {
    stop("`value` contains only NA values.")
  }
  
  lim_source <- if (scale_type == "diverging") abs(vals_non_na) else vals_non_na
  
  hm_lim <- round(
    if (!is.null(hm_max)) {
      hm_max
    } else if (use_quantile) {
      as.numeric(stats::quantile(lim_source, probs = quantile_val, na.rm = TRUE))
    } else {
      max(lim_source, na.rm = TRUE)
    },
    1
  )
  
  if (!is.finite(hm_lim)) {
    stop("Could not determine heatmap scale limit.")
  }
  
  if (hm_lim == 0) {
    hm_lim <- 1
  }
  
  if (palette %in% scico:::scico_palette_names()) {
    pal_colors <- scico::scico(n_colors, palette = palette)
  } else if (palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    base_cols <- RColorBrewer::brewer.pal(max_n, palette)
    pal_colors <- grDevices::colorRampPalette(base_cols)(n_colors)
  } else {
    stop("Palette not recognized. Use a valid scico or RColorBrewer palette name.")
  }
  
  if (palette_rev) {
    pal_colors <- rev(pal_colors)
  }
  
  if (scale_type == "diverging") {
    breaks <- seq(-hm_lim, hm_lim, length.out = n_colors)
    legend_at <- c(-hm_lim, 0, hm_lim)
    hm_palette <- circlize::colorRamp2(breaks, pal_colors)
  } else {
    if (any(vals_non_na < 0)) {
      warning("Negative values detected with scale_type = 'continuous'. Values below 0 will be clipped by the color scale.")
    }
    
    breaks <- seq(0, hm_lim, length.out = n_colors)
    pal_colors[1] <- "#FFFFFF"
    legend_at <- c(0, hm_lim / 2, hm_lim)
    hm_palette <- circlize::colorRamp2(breaks, pal_colors)
  }
  
  ncol_hm <- hm_dat %>% dplyr::distinct(!!col_var) %>% nrow()
  nrow_hm <- hm_dat %>% dplyr::distinct(!!row_var) %>% nrow()
  
  hm <- hm_dat |>
    tidyHeatmap::heatmap(
      .row = !!row_var,
      .column = !!col_var,
      .value = !!val_var,
      palette_value = hm_palette,
      heatmap_legend_param = list(
        color_bar = "continuous",
        at = legend_at
      ),
      cluster_rows = cluster_rows,
      cluster_columns = cluster_cols,
      row_title = NULL,
      column_title = NULL,
      show_row_dend = row_dend,
      show_column_dend = col_dend,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = row_font),
      column_names_gp = grid::gpar(fontsize = col_font),
      column_names_rot = if (rotate_text) 45 else 0,
      border = TRUE,
      width = ncol_hm * grid::unit(col_size, "mm"),
      height = nrow_hm * grid::unit(row_size, "mm")
    )
  
  if (asterisk && !rlang::quo_is_null(sig_var)) {
    hm <- hm |>
      tidyHeatmap::layer_asterisk(!!sig_var < 0.05)
  }
  
  tidyHeatmap::wrap_heatmap(hm) +
    ggplot2::labs(title = title, subtitle = subtitle)
}
