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
