#' Title plot PCA
#'
#' @return
#' @export
#'
#' @examples
plot_PCA <- function(dds, save_Plot = FALSE, file_Dir = "./", file_Name = "All_Samples_PCA_removeBatch.pdf") {
  rld <- DESeq2::rlog(dds, blind = FALSE)
  data <- DESeq2::plotPCA(rld, returnData = TRUE, intgroup = "Condition")
  percentVar <- round(100 * attr(data, "percentVar"))
  pca <- ggpubr::ggscatter(data, "PC1", "PC2",
    color = "Condition",
    label = "name", size = 5, repel = T,
    palette = "jama", # 杂志jama的配色
    ellipse = TRUE, # 画椭圆
    mean.point = F,
    star.plot = F # 生成星图
  ) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_bw(base_size = 18, base_line_size = 1) +
    ggtitle("Principal Component Analysis")

  if (save_Plot) {
    ggplot2::ggsave(paste0(file_Dir, file_Name), width = 8, height = 6, plot = pca)
  }
  pca
}



# .plot_PCA <- function(dds,save_Plot = FALSE,file_Dir = './',file_Name = 'All_Samples_PCA_removeBatch.pdf'){
#   rld <- DESeq2::rlog(dds, blind = FALSE)
#   pheno <- SummarizedExperiment::colData(dds)
#   data_pca <- as.data.frame(t(as.data.frame(SummarizedExperiment::assay(rld))[,pheno$Sample]))
#   res.pca <- FactoMineR::PCA(data_pca, graph = FALSE)
#   pca <- factoextra::fviz_pca_ind(res.pca,
#                       geom = c("point", "text"),
#                       geom.ind = c("point", "text"), # show points only (nbut not "text")
#                       col.ind = as.character(pheno$Condition), # color by groups
#                       fill.ind = pheno$Condition,
#                       palette = "jco",
#                       mean.point = FALSE,pointsize = 3,
#                       pointshape = 20,
#                       repel = T,
#                       ellipse.type = "confidence",
#                       addEllipses = T, # Concentration ellipses
#                       legend.title = "Groups",
#                       ggtheme = theme_bw(),
#                       title = "Principal Component Analysis")
#   if (save_Plot){
#     ggplot2::ggsave(paste0(file_Dir,file_Name),width = 8,height = 6,plot = pca)
#   }
#   pca
# }


#' Title
#'
#' @param dds
#' @param fileDir
#' @param type
#'
#' @return
#' @export
#'
#' @examples
Plot_corr_heatmap <- function(dds, fileDir = "./", type = "pearson", savePlot = FALSE) {
  rld <- DESeq2::rlog(dds, blind = FALSE)
  cor_matrix <- stats::cor(SummarizedExperiment::assay(rld), method = type)
  cor_heatmap <- pheatmap::pheatmap(cor_matrix,
    show_rownames = T, show_colnames = T,
    main = "Sample correlation coefficient",
    fontsize = 15
  )
  if (savePlot == TRUE) {
    ggplot2::ggsave(paste0(fileDir, type, "_Samples_Corr.pdf"), width = 8, height = 6, plot = cor_heatmap)
  } else {
    return(cor_heatmap)
  }
}

#' Title
#'
#' @param Diff_res
#'
#' @return
#' @export
#'
#' @examples
getDEtab <- function(Diff_res) {
  DEdf <- Diff_res$Res
  DEdf$logP <- -log10(DEdf$padj)
  DEdf$Group <- "Not-significant"
  DEdf$Group[which((DEdf$padj < 0.05) & (-log2(2) >= DEdf$log2FoldChange & DEdf$log2FoldChange > -log2(5)))] <- "Down-regulated (FC:-2~-5)"
  DEdf$Group[which((DEdf$padj < 0.05) & (-log2(5) >= DEdf$log2FoldChange))] <- "Down-regulated (FC:< -5)"
  DEdf$Group[which((DEdf$padj < 0.05) & (log2(5) > DEdf$log2FoldChange & DEdf$log2FoldChange >= 1))] <- "Up-regulated (FC:2-5)"
  DEdf$Group[which((DEdf$padj < 0.05) & (DEdf$log2FoldChange >= log2(5)))] <- "Up-regulated (FC:>5)"
  down_genes <- subset(DEdf, DEdf$Group == "Down-regulated (FC:< -5)" | DEdf$Group == "Down-regulated (FC:-2~-5)") %>%
    as.data.frame() %>%
    dplyr::arrange(padj)
  up_genes <- subset(DEdf, DEdf$Group == "Up-regulated (FC:>5)" | DEdf$Group == "Up-regulated (FC:2-5)") %>%
    as.data.frame() %>%
    dplyr::arrange(padj)
  top10_label <- c(as.character(head(up_genes$Symbol, 5)), as.character(head(down_genes$Symbol, 5)))
  DEdf$Label <- ""
  DEdf$Label[match(top10_label, DEdf$Symbol)] <- top10_label
  DEdf$log2FC <- -DEdf$log2FoldChange
  DE_res <- list(DEgenes = DEdf, upGenes = up_genes, downGenes = down_genes)
  return(DE_res)
}



#' Title
#'
#' @param Diff_res
#' @param Dir
#'
#' @return
#' @export
#'
#' @examples
Plot_volcano <- function(Diff_res, Dir = "./", save_Plot = TRUE) {
  DEtab <- getDEtab(Diff_res = Diff_res)
  prefix <- paste(Diff_res$Group[1], Diff_res$Group[2], sep = "_VS_")
  filename <- paste(Dir, prefix, "_volcano.pdf")
  npg5 <- c("#F39B7FFF", "#3C5488FF", "#00A087FF", "#E64B35FF", "#4DBBD5FF")
  df <- data.frame(stats::na.omit(DEtab$DEgenes))
  volcano_plot <- ggpubr::ggscatter(df,
    x = "log2FoldChange", y = "logP", color = "Group", alpha = 0.8,
    palette = npg5, size = 1.8,
    xlab = "log2FoldChange", ylab = "-log10 (Adjust P-value)"
  ) +
    geom_hline(yintercept = 1.3, linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    ggpubr::theme_pubr(legend = "right", border = T) +
    ggplot2::annotate("text", x = -5, y = -2, label = "#mRNAs :", color = "black", size = 2.5) +
    ggplot2::annotate("text", x = -3.5, y = -2, label = table(df$Group)["Down-regulated (FC:< -5)"], color = "#3C5488FF", size = 2.5) +
    ggplot2::annotate("text", x = 3.5, y = -2, label = table(df$Group)["Up-regulated (FC:>5)"], color = "#E64B35FF", size = 2.5) +
    ggplot2::annotate("text", x = -2, y = -2, label = table(df$Group)["Down-regulated (FC:-2~-5)"], color = "#F39B7FFF", size = 2.5) +
    ggplot2::annotate("text", x = 2, y = -2, label = table(df$Group)["Up-regulated (FC:2-5)"], color = "#4DBBD5FF", size = 2.5) +
    ggpubr::font("xylab", size = 15, face = "bold") +
    ggpubr::font("title", size = 12, face = "bold") +
    ggpubr::font("xy.text", size = 12, face = "bold") +
    ggtitle(label = paste("Volcano plot : ", prefix, sep = ""))
  if (save_Plot == TRUE) {
    ggplot2::ggsave(filename = filename, plot = volcano_plot, width = 6, height = 4)
  }
  return(volcano_plot)
}


#' Title
#'
#' @param Diff_res
#' @param dds
#' @param Dir
#'
#' @return
#' @export
#'
#' @examples
Plot_heatmap <- function(Diff_res, dds, Dir = "./", save_Plot = TRUE, show_rownames = FALSE) {
  pheno <- SummarizedExperiment::colData(dds)
  logExpression <- as.data.frame(SummarizedExperiment::assay(rlog(dds, blind = FALSE)))
  prefix <- paste(Diff_res$Group[1], Diff_res$Group[2], sep = "_VS_")
  filename <- paste(Dir, prefix, "_DEgenes_Heatmap.pdf")
  tmp <- pheno[pheno$Condition %in% Diff_res$Group, ]
  group <- data.frame(row.names = rownames(tmp), Condition = tmp$Condition)
  upExpression <- logExpression[Diff_res$up_gene$SYMBOL, rownames(group)]
  downExpression <- logExpression[Diff_res$down_gene$SYMBOL, rownames(group)]
  if (nrow(upExpression) < 2 | nrow(downExpression) < 2) {
    print("DE genes number < 2,don't need to plot")
  } else {
    upExpression <- upExpression[order(rowMeans(upExpression), decreasing = TRUE), ]
    downExpression <- downExpression[order(rowMeans(downExpression), decreasing = TRUE), ]
    subexpression <- rbind(upExpression, downExpression)
    heat <- pheatmap::pheatmap(subexpression,
      show_rownames = show_rownames, annotation_col = group, fontsize_row = 10, cluster_cols = F,
      # color = grDevices::colorRampPalette(c(2, 'white', 3))(100),
      cluster_rows = FALSE, scale = "row",
      main = prefix
    )
    print("DE Genes Heatmap Plot Finish")
  }
  return(heat)
  if (save_Plot == TRUE) {
    ggplot2::ggsave(filename, plot = heat, width = nrow(group) * 0.6 + 0.5, height = nrow(subexpression) * 0.06 + 0.6)
  }
}



#' Title
#'
#' @param Diff_res
#' @param Dir
#'
#' @return
#' @export
#'
#' @examples
PlotKEGG <- function(Diff_res, Dir = "./", save_Plot = TRUE) {
  universeIDs <- rownames(Diff_res$res)
  prefix <- paste(Diff_res$Group[1], Diff_res$Group[2], sep = "_VS_")
  upfilename <- paste(Dir, prefix, "_KEGG_Pathway_Up.pdf")
  downfilename <- paste(Dir, prefix, "_KEGG_Pathway_Down.pdf")
  KEGG_up <- enricher(Diff_res$up_gene$SYMBOL,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = universeIDs,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    TERM2GENE,
    TERM2NAME = TERM2NAME
  )
  dotplot_KEGG_up <- dotplot(KEGG_up)
  dotplot_KEGG_up
  KEGG_down <- enricher(Diff_res$down_gene$SYMBOL,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = universeIDs,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    TERM2GENE,
    TERM2NAME = TERM2NAME
  )
  dotplot_KEGG_down <- dotplot(KEGG_down)
  dotplot_KEGG_down
  if (save_Plot == TRUE) {
    ggplot2::ggsave(filename = downfilename, plot = dotplot_KEGG_down, width = 8, height = 6)
    ggplot2::ggsave(filename = upfilename, plot = dotplot_KEGG_up, width = 8, height = 6)
  }
}


#' Title
#'
#' @param Diff_res
#' @param DIR
#' @param p_cutoff
#' @param version
#'
#' @return
#' @export
#'
#' @examples
Plot_DE_GOTerms <- function(Diff_res, DIR = "./", p_cutoff = 0.05, version) {
  prefix <- paste(Diff_res$Group[1], Diff_res$Group[2], sep = "_VS_")
  GO_tab <- getGO(Diff_res$up_gene$SYMBOL, version = version)$GO_df %>%
    dplyr::arrange(type, p.adjust) %>%
    dplyr::filter(p.adjust < p_cutoff)
  test <- data.frame(row.names = GO_tab$Description, `-log10(padjust)` = -log10(GO_tab$p.adjust))
  anno <- data.frame(row.names = GO_tab$Description, type = factor(GO_tab$type))
  # anno_color = list(type=c("Biological Process"="#EE0000FF",
  #                          "Cellular Component"="#631879FF","Molecular Function" = "#008280FF"))
  pheatmap::pheatmap(test,
    cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = anno,
    fontsize_row = 22,
    fontsize_col = 8,
    cellwidth = 40, cellheight = 30,
    show_colnames = F,
    show_rownames = TRUE,
    legend = TRUE,
    scale = "none",
    color = colorRampPalette(c("white", "firebrick3"), bias = 1)(100), main = paste("High expression in :", Diff_res$Group[1]),
    filename = paste(DIR, prefix, paste("_High_", Diff_res$Group[1]), "_GOTerms.pdf")
  )

  GO_tab2 <- getGO(Diff_res$down_gene$SYMBOL, version = version)$GO_df %>%
    arrange(type, p.adjust)
  test2 <- data.frame(row.names = GO_tab2$Description, `-log10(padjust)` = -log10(GO_tab2$p.adjust))
  anno <- data.frame(row.names = GO_tab2$Description, type = factor(GO_tab2$type))
  pheatmap::pheatmap(test2,
    cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = anno,
    fontsize_row = 22,
    fontsize_col = 8,
    cellwidth = 40, cellheight = 30,
    show_colnames = F,
    show_rownames = TRUE,
    legend = TRUE,
    scale = "none",
    color = colorRampPalette(c("white", "navy"), bias = 1)(100), main = paste("High expression in :", Diff_res$Group[2]),
    filename = paste(DIR, prefix, paste("_High_", Diff_res$Group[2]), "_GOTerms.pdf")
  )
}



#' Title
#'
#' @param dds
#' @param ID
#'
#' @return
#' @export
#'
#' @examples
Barplot_ID <- function(dds, ID) {
  NormData <- DESeq2::counts(dds, normalized = TRUE)
  pheno <- SummarizedExperiment::colData(dds)
  bardat <- NormData[ID, ]
  gene_dat <- data.frame(NormCounts = bardat, Condition = pheno$Condition)
  gene_dat$NormCounts <- as.numeric(gene_dat$NormCounts)
  ggpubr::ggboxplot(gene_dat,
    y = "NormCounts", x = "Condition", color = "Condition",
    add = "mean_se", palette = "aaas", bxp.errorbar = T
  ) +
    ggpubr::rremove("x.ticks") +
    ggplot2::ylab("Normalized Counts") + ggplot2::xlab("") +
    ggpubr::theme_pubr(border = F, x.text.angle = 45, legend = "right") +
    ggpubr::labs_pubr() +
    ggplot2::ggtitle(id2name[ID, "Gene_Name"])
}

#' Title
#'
#' @param df
#' @param n
#' @param wt
#'
#' @return
#' @export
#'
#' @examples
get_TopLabel <- function(df, n = 5) {
  Pos_topn <- df %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::top_n(n = n, log2FoldChange)
  Neg_topn <- df %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::top_n(n = -n, log2FoldChange)
  Label_n <- rbind(Pos_topn, Neg_topn)$SYMBOL
  return(Label_n)
}


#' Title
#'
#' @param Diff_res
#' @param n
#' @param Dir
#' @param save_Plot
#' @param p_cutoff
#' @param FC_cutoff
#' @param wt
#'
#' @return
#' @export
#'
#' @examples
Plot_enhancedVolcano <- function(
    Diff_res, n = 5, Dir = "./", save_Plot = FALSE,
    p_cutoff = 0.05, FC_cutoff = 1) {
  DF <- as.data.frame(Diff_res$Res)
  rownames(DF) <- Diff_res$Res$SYMBOL
  selectLab <- get_TopLabel(DF, n)
  prefix <- paste(Diff_res$Group[1], Diff_res$Group[2], sep = "_VS_")
  filename <- paste(Dir, prefix, "_enhanced_volcano.pdf")
  en_Volcano_Plot <- EnhancedVolcano::EnhancedVolcano(DF,
    lab = rownames(DF),
    selectLab = selectLab,
    x = "log2FoldChange",
    y = "padj",
    title = prefix, drawConnectors = T,
    pCutoff = p_cutoff,
    FCcutoff = FC_cutoff
  )
  if (save_Plot == TRUE) {
    ggplot2::ggsave(filename = filename, plot = en_Volcano_Plot, width = 12, height = 10)
  }
  return(en_Volcano_Plot)
}
