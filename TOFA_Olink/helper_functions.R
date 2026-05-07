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

# Excel export function -----
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

# labelled Volcano plot function for lms -----
volcano_plot_lab_lm <- function(res, title = "", 
                                subtitle = "down in Pos.                                                  up in Pos.",
                                y_lim = c(0, NA)){
  res <- res |> 
    mutate(
      color = if_else(qvalue < 0.1, "padj < 0.1", "All")
    )
  # get max for x-axis
  x_lim <- res |> 
    summarize(max = max(log2FoldChange, na.rm = TRUE), min = min(log2FoldChange, na.rm = TRUE)) |> 
    abs() |> 
    max() |> 
    ceiling()
  res |> 
    ggplot(aes(log2FoldChange, -log10(qvalue), color = color)) + 
    geom_hline(yintercept = -log10(0.1), linetype = 2) + 
    geom_vline(xintercept = 0, linetype = 2) + 
    geom_point() + 
    scale_color_manual(values = c("padj < 0.1" = "red", "All" = "black")) + 
    xlim(-x_lim, x_lim) +
    ylim(y_lim) +
    geom_text_repel(data = res |>  filter(qvalue<0.1, log2FoldChange>0) |>  slice_max(order_by = log2FoldChange, n = 5), aes(label = Assay), min.segment.length = 0, show.legend = FALSE, size = 3) +
    geom_text_repel(data = res |>  filter(qvalue<0.1, log2FoldChange<0) |>  slice_min(order_by = log2FoldChange, n = 5), aes(label = Assay), min.segment.length = 0, show.legend = FALSE, size = 3) +
    # geom_text_repel(data = res %>% filter(!is.na(qvalue)) %>% slice_min(order_by = qvalue, n = 3), aes(label = Assay), min.segment.length = 0, show.legend = FALSE, nudge_x = -0, nudge_y = 0.1) +
    theme(aspect.ratio=1.2) +
    labs(
      title = title,
      subtitle = subtitle
    )
} # end of function

# GSEA FUNCTIONS -----
# function to get combined pos and neg GSEA results -----
run_fgsea2 <- function(geneset, ranks, weighted = FALSE) {
  library("fgsea")
  # with gseaParam = 0, results are VERY similar to original GSEA # this seems to not be operating as expected as N^0 = 1, so all ranking stats would be 1
  weight = 0
  if(weighted) weight = 1
  # Run positive enrichment
  fgseaRes_POSITIVE <- fgseaMultilevel(
    geneset, 
    ranks, 
    minSize=15, 
    maxSize=500,
    gseaParam = weight,
    # nperm = 1000,
    eps = 0.0, # fgsea has a default lower bound eps=1e-10 for estimating P-values. If you need to estimate P-value more accurately, you can set the eps argument to zero
    scoreType = "pos"
  )
  # Run negative enrichment
  fgseaRes_NEGATIVE <- fgseaMultilevel(
    geneset,
    ranks,
    minSize=15,
    maxSize=500,
    gseaParam = 0,
    # nperm = 1000,
    eps = 0.0, # fgsea has a default lower bound eps=1e-10 for estimating P-values. If you need to estimate P-value more accurately, you can set the eps argument to zero
    scoreType = "neg"
  )
  # Combine positive and negative results + re-adjust pvals
  fgseaRes_POS_NEG <- inner_join(
    fgseaRes_POSITIVE %>% 
      as_tibble(),
    fgseaRes_NEGATIVE %>% 
      as_tibble(),
    by = c("pathway"),
    suffix = c("_POS", "_NEG")
  )
  fgseaRes_COMBINED <- bind_rows(
    fgseaRes_POS_NEG %>% filter(ES_POS > abs(ES_NEG)) %>% select(pathway) %>% inner_join(fgseaRes_POSITIVE),
    fgseaRes_POS_NEG %>% filter(ES_POS < abs(ES_NEG)) %>% select(pathway) %>% inner_join(fgseaRes_NEGATIVE)
  ) %>% 
    mutate(padj = p.adjust(pval, method = "BH"))%>% 
    arrange(padj, -abs(NES))
  return(fgseaRes_COMBINED)
} # end of function
