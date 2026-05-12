# Common helper functions
mem_used <- function() lobstr::mem_used() %>% as.numeric() %>% R.utils::hsize()
obj_size <- function(x) object.size(x) %>% print(units = "auto")

# Setting and modifying default theme for plots
theme_set(theme_gray(base_size=12, base_family="Arial") +
            theme(
              panel.border=element_rect(colour="black", fill="transparent"),
              plot.title=element_text(face="bold", hjust=0),
              axis.text=element_text(color="black", size=14),
              axis.text.x=element_text(angle=0, hjust=0.5),
              axis.ticks = element_line(color = "black"), # make sure tick marks are black
              panel.background=element_blank(),
              panel.grid=element_blank(),
              plot.background=element_blank(),
              strip.background = element_blank(), # facet label borders
              legend.key=element_blank(), legend.background=element_blank() # remove grey bg from legend
            )
)
#

## Excel export function
export_excel <- function(named_list, filename = "") {
  wb <- openxlsx::createWorkbook()
  ## Loop through the list of split tables as well as their names
  ## and add each one as a sheet to the workbook
  Map(function(data, name){
    openxlsx::addWorksheet(wb, name)
    openxlsx::writeData(wb, name, data)
  }, named_list, names(named_list))
  ## Save workbook to working directory
  openxlsx::saveWorkbook(wb, file = here("results", paste0(out_file_prefix, filename, ".xlsx")), overwrite = TRUE)
  cat("Saved as:", here("results", paste0(out_file_prefix, filename, ".xlsx")))
} # end of function


# get standard ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# effectsize tidyHeatmap -----
hm_effsize <- function(hm_dat, title = "Heatmap title", subtitle = "") {
  hm_lim <- hm_dat |>
    pull(effsize) |>
    abs() |>
    max() |>
    round(2)
  breaks <- seq(-hm_lim, hm_lim, length.out = 11)
  hm_palette <- circlize::colorRamp2(
    breaks,
    RColorBrewer::brewer.pal(11, "RdBu") |> rev()
  )
  # plot heatmap
  hm_dat |>
    tidyHeatmap::heatmap(
      score_name_short,
      timepoint,
      effsize,
      palette_value = hm_palette,
      heatmap_legend_param = list(color_bar = "continuous", at = seq(-hm_lim, hm_lim, length.out = 5)),
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      row_title = NULL,
      show_column_names = TRUE,
      column_title = NULL,
      border = TRUE
    ) |>
    tidyHeatmap::layer_asterisk(BH_padj_endpoints < 0.1) |>
    tidyHeatmap::wrap_heatmap() +
    labs(
      title = title,
      subtitle = subtitle
    )
}
#
hm_effsize_p <- function(hm_dat, title = "Heatmap title", subtitle = "") {
  hm_lim <- hm_dat |>
    pull(effsize) |>
    abs() |>
    max() |>
    round(2)
  breaks <- seq(-hm_lim, hm_lim, length.out = 11)
  hm_palette <- circlize::colorRamp2(
    breaks,
    RColorBrewer::brewer.pal(11, "RdBu") |> rev()
  )
  # plot heatmap
  hm_dat |>
    tidyHeatmap::heatmap(
      score_name_short,
      timepoint,
      effsize,
      palette_value = hm_palette,
      heatmap_legend_param = list(color_bar = "continuous", at = seq(-hm_lim, hm_lim, length.out = 5)),
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      row_title = NULL,
      show_column_names = TRUE,
      column_title = NULL,
      border = TRUE
    ) |>
    tidyHeatmap::layer_asterisk(p < 0.05) |>
    tidyHeatmap::wrap_heatmap() +
    labs(
      title = title,
      subtitle = subtitle
    )
}
#


# effectsize tidyHeatmap stratified -----
hm_effsize_strat <- function(hm_dat, title = "Heatmap title", subtitle = "") {
  hm_lim <- hm_dat |>
    pull(effsize) |>
    abs() |>
    max() |>
    round(2)
  breaks <- seq(-hm_lim, hm_lim, length.out = 11)
  hm_palette <- circlize::colorRamp2(
    breaks,
    RColorBrewer::brewer.pal(11, "RdBu") |> rev()
  )
  # plot heatmap
  hm_dat |>
    tidyHeatmap::heatmap(
      score_name_short,
      level,
      effsize,
      palette_value = hm_palette,
      heatmap_legend_param = list(color_bar = "continuous", at = seq(-hm_lim, hm_lim, length.out = 5)),
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      row_title = NULL,
      show_column_names = TRUE,
      column_title = NULL,
      border = TRUE
    ) |>
    tidyHeatmap::layer_asterisk(p < 0.05) |>
    tidyHeatmap::wrap_heatmap() +
    labs(
      title = title,
      subtitle = subtitle
    )
}


# LMM tidyHeatmap stratified -----
hm_lmm_strat <- function(hm_dat, title = "Heatmap title", subtitle = "") {
  hm_lim <- hm_dat |>
    pull(estimate) |>
    abs() |>
    max() |>
    round(2)
  breaks <- seq(-hm_lim, hm_lim, length.out = 11)
  hm_palette <- circlize::colorRamp2(
    breaks,
    RColorBrewer::brewer.pal(11, "RdBu") |> rev()
  )
  # plot heatmap
  hm_dat |>
    tidyHeatmap::heatmap(
      score_name_short,
      level,
      estimate,
      palette_value = hm_palette,
      heatmap_legend_param = list(color_bar = "continuous", at = seq(-hm_lim, hm_lim, length.out = 5)),
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      row_title = NULL,
      show_column_names = TRUE,
      column_title = NULL,
      border = TRUE
    ) |>
    tidyHeatmap::layer_asterisk(BH_padj < 0.1) |>
    tidyHeatmap::wrap_heatmap() +
    labs(
      title = title,
      subtitle = subtitle
    )
}

