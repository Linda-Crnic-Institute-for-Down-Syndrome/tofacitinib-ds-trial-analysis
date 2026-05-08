# Common helper functions
# by Matthew Galbraith
mem_used <- function() lobstr::mem_used() %>% as.numeric() %>% R.utils::hsize()
obj_size <- function(x) object.size(x) %>% print(units = "auto")

## Setting and modifying default theme for plots
# theme_set(theme_gray(base_size=12, base_family="Arial") +
#             theme(
#               panel.border=element_rect(colour="black", fill="transparent"), 
#               plot.title=element_text(face="bold", hjust=0),
#               axis.text=element_text(color="black", size=14), 
#               axis.text.x=element_text(angle=0, hjust=0.5),
#               axis.ticks = element_line(color = "black"), # make sure tick marks are black
#               panel.background=element_blank(),
#               panel.grid=element_blank(),
#               plot.background=element_blank(),
#               strip.background = element_blank(), # facet label borders
#               legend.key=element_blank(), legend.background=element_blank() # remove grey bg from legend
#             )
# )
# pull from gist
suppressMessages(devtools::source_gist("https://gist.github.com/mattgalbraith/f082ed7d152729f4ae72383e564a70e8", filename = "ggplot_theme.R"))
# may need to add/update personal access token
# usethis::gh_token_help() # check token status
# usethis::git_sitrep()
# https://usethis.r-lib.org/articles/git-credentials.html
# Details in this gist: https://gist.github.com/mattgalbraith/0f9f2d75023be5355b693cb832b9abef
# usethis::create_github_token() # takes you to a pre-filled form to create a new PAT
# gitcreds::gitcreds_set() # register this token in the local Git credential store


## Density color function
getDenCols <- function(x, y, transform = TRUE) { # set to TRUE if using log2 transformation of data
  if(transform) {
    df <- data.frame(log2(x), log2(y))
  } else{
    df <- data.frame(x, y)
  }
  z <- grDevices::densCols(df, colramp = grDevices::colorRampPalette(c("black", "white")))
  df$dens <- grDevices::col2rgb(z)[1,] + 1L
  cols <-  grDevices::colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  return(df$dens)
} # End of function

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
    # arrange(-median_diff) |>
    # group_by(group) |>
    # mutate(lab_name = fct_inorder(lab_name)) |>
    tidyHeatmap::heatmap(
      lab_name,
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
    # tidyHeatmap::layer_asterisk(BH_padj < 0.1) |> 
    tidyHeatmap::layer_asterisk(p < 0.05) |>
    # tidyHeatmap::layer_point(p < 0.05) |> 
    tidyHeatmap::wrap_heatmap() +
    labs(
      title = title,
      subtitle = subtitle
    )
}
#



