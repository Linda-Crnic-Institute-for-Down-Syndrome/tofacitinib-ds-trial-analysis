# Common helper functions for DESeq
#
# v0.3 -31125 by Matthew Galbraith
mem_used <- function() lobstr::mem_used() %>% as.numeric() %>% R.utils::hsize()
obj_size <- function(x) object.size(x) %>% print(units = "auto")

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

# get size factors function ----
get_size_fcts <- function(x) {
  name <- deparse(substitute(x))
  x |>  
    sizeFactors() |>  
    enframe(name = "Sampleid", value = "SizeFactor")
}

# get standard ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Functions to print results summaries ± independent filtering -------
# CONTRASTS VERSION
get_results_sum <- function(x, contrast = c(predictor, "denom", "num"), show_ind_filt_off=TRUE) {
  message("------------\nResults summary for ", contrast[2], " vs. ", contrast[3])
  message("Model formula: ", design(x), "\n------")
  if (show_ind_filt_off) {
    cat("independentFiltering OFF, Cooks cutoff ON")
    x %>%
      results(
        contrast,
        # name = name,
        cooksCutoff=TRUE,
        independentFiltering=FALSE
      ) %>%
      summary()
    message("----\n")
  } else {
    cat("independentFiltering ON, Cooks cutoff ON")
    x %>%
      results(
        contrast,
        # name = name,
        cooksCutoff=TRUE,
        independentFiltering=TRUE
      ) %>%
      summary()
    message("----\n")
  }
}
#

# Functions to generate final results tbls in our format (use with for loops to get all comparisons) ------
# CONTRASTS VERSION
get_results_tbl <- function(x, contrast = "", cooks = TRUE, ind_filt = TRUE, shrink_type = "apeglm") {
  #
  message("------\nGenerating DEseq2 results for ", paste0(contrast[1], ": ", contrast[2], " vs. ", contrast[3]))
  message("Model formula: ", design(x), "\n------")
  ## relevel variable of interest (and re-run nbinomWaldTest()) to ensure that
  # coefficient of interest is available (check with resultsNames()) for
  # lfcShrink() - otherwise need to use contrast argument which does not allow
  # use of apeglm shrinkage
  message(paste0("1. Re-leveling '", contrast[1], "' with '", contrast[3], "' as reference level"))
  # no longer hardcoded - see ?`SummarizedExperiment-class` for info on accessors for colData
  x[[contrast[1],]] <- x[[contrast[1],]] |> relevel(contrast[3])
  #
  message("2. Running negative binomial Wald test")
  x <- x %>% nbinomWaldTest(quiet=TRUE)
  ## Get results without LFC shrinkage (default)
  res <- x %>%
    results(
      contrast,
      cooksCutoff = cooks,
      independentFiltering = ind_filt
    )
  res_tbl <- res %>% # Could also use biobroom::tidy() on results or dds but has less useful colnames
    as.data.frame() %>%
    as_tibble(rownames = "Geneid")
  ## Apply LFC shrinkage
  message("3. Calculating log2 fold-change shrinkage")
  # NOTES ON LOG FOLD CHANGE SHRINKAGE
  # earlier DESeq2 versions (<1.16) carried out LFC shrinkage by
  # default(betaPrior=TRUE); more recent versions set betaPrior=FALSE and use
  # lfcShrink() in a separate step.
  # https://support.bioconductor.org/p/95695/
  # However: 
  # with betaPrior=TRUE, p-values are calculated for the shrunken LFC, while
  # betaPrior=FALSE + subsequent lfcShrink() calculates p-values based on
  # un-shrunken LFC and only shrinks them afterwards.
  # Also:
  # Difference between using lfcShrink() with the coef argument and using lfcShrink() with the contrast argument
  # lfcShrink(dds=dds, coef=2, res=res) OR lfcShrink(dds=dds, contrast=c("condition","B","A"), res=res)
  # https://support.bioconductor.org/p/98833/#98837 From Michael Love:
  # They are not identical. Using contrast is similar to what DESeq2 used to do:
  # it forms an expanded model matrix, treating all factor levels equally, and
  # averages over all distances between all pairs of factor levels to estimate the
  # prior. Using coef, it just looks at that column of the model matrix (so
  # usually that would be one level against the reference level) and estimates the
  # prior for that coefficient from the distribution of those MLE of coefficients.
  # I implemented both for lfcShrink, because 'contrast' provides backward support
  # (letting people get the same coefficient they obtained with previous
  # versions), while future types of shrinkage estimators will use the 'coef'
  # approach, which is much simpler.
  # my current recommendation would be to use the p-values from un-shrunken LFC
  # and then use the shrunken LFC for visualization or ranking of genes. This is
  # the table you get with default DESeq => results => lfcShrink. If you want to
  # be future-proof, I'd go with 'coef' with lfcShrink. All the methods I'm
  # planning on adding (ours and others) really just want to shrink one
  # coefficient at a time, not do the expanded model matrix thing.
  res_shrink <- x %>%
    lfcShrink(
      coef = paste(contrast[1], contrast[2], "vs", contrast[3], sep="_"), # REPLACE dashes or this fails
      type = shrink_type, # "apeglm" or "normal" or "ashr"; apeglim and ashr are better at preserving large LFCs
      parallel = FALSE,
      res = res
    )
  res_shrink_tbl <- res_shrink %>%
    as.data.frame() %>%
    as_tibble(rownames="Geneid")
  # Combine and mutate to get final results tbl
  message("4. Assembling final results table\n------")
  res_tbl %>%
    inner_join(res_shrink_tbl %>% select(Geneid, log2FoldChange), by = "Geneid") %>%
    dplyr::rename(
      log2FoldChange = log2FoldChange.x,
      log2FoldChange_adj = log2FoldChange.y
    ) %>%
    mutate(
      FoldChange = 2^log2FoldChange,
      FoldChange_adj = 2^log2FoldChange_adj,
      Model = design(x) |> paste(collapse = "")
    ) %>%
    dplyr::select(
      Geneid:baseMean,
      Model,
      FoldChange,
      log2FoldChange,
      FoldChange_adj,
      log2FoldChange_adj,
      pvalue,
      padj
    ) %>%
    arrange(padj) %>%
    inner_join(gene_anno, by="Geneid") %>%
    select(Gene_name = gene_name, chr, everything())
}

# labelled Volcano plot function ----
volcano_plot_lab <- function(res, 
                             labels = TRUE,
                             n_labels = 2,
                             title = "", 
                             subtitle = "",
                             y_lim = c(0, NA),
                             raster = FALSE
){res <- res %>% 
    mutate(
      color = if_else(padj < 0.1, "padj < 0.1", "All")
    )
  # Get max finite -log10(pval) and replace 0s if needed
  max_finite <- res %>%
    filter(padj > 0) %$%
    min(padj) %>%
    -log10(.)
  res <- res %>% mutate(
    shape = if_else(padj == 0, "infinite", "finite"),
    padj = if_else(padj == 0, 10^-(max_finite * 1.05), padj)
  )
  # get max for x-axis
  x_lim <- res %>% 
    summarize(max = max(log2(FoldChange_adj), na.rm = TRUE), min = min(log2(FoldChange_adj), na.rm = TRUE)) %>% 
    abs() %>% 
    max() %>% 
    ceiling()
  p <- res %>% 
    ggplot(aes(log2(FoldChange_adj), -log10(padj), color = color, shape = shape)) + 
    geom_hline(yintercept = -log10(0.1), linetype = 2) + 
    geom_vline(xintercept = 0, linetype = 2) + 
    geom_point() + 
    scale_color_manual(values = c("padj < 0.1" = "red", "All" = "black")) + 
    scale_shape_manual(values = c("infinite" = 2, "finite" = 16)) +
    guides(shape = "none") +
    # xlim(-x_lim, x_lim) +
    # ylim(y_lim) +
    theme(aspect.ratio=1.2) +
    labs(
      title = title,
      subtitle = subtitle
    )
  if(labels == TRUE) {
    p <- p +
      geom_text_repel(data = res %>% slice_min(order_by = padj, n = n_labels), aes(label = Gene_name), min.segment.length = 0, show.legend = FALSE) + # , nudge_y = -max_finite / 5
      geom_text_repel(data = res %>% slice_min(order_by = FoldChange_adj, n = n_labels) %>% filter(padj < 0.1), aes(label = Gene_name), min.segment.length = 0, show.legend = FALSE, nudge_x = -x_lim / 5, nudge_y = max_finite / 20, ylim = c(max_finite / 10, NA)) +
      geom_text_repel(data = res %>% slice_max(order_by = FoldChange_adj, n = n_labels) %>% filter(padj < 0.1), aes(label = Gene_name), min.segment.length = 0, show.legend = FALSE, nudge_x = x_lim / 5, nudge_y =  max_finite / 20, ylim = c(max_finite / 10, NA))
  }
  #
  if(raster) {
    # to rasterize all points:
    # rasterize(p, layers='Point', dpi = 600)
    p <- rasterize(p, layers='Point', dpi = 600, dev = "ragg_png")
    # otherwize can rasterize individual layers using rasterize(geom_sina())
  }
  #
  return(p)
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




