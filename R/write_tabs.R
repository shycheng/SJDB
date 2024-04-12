#' Title
#'
#' @param dds
#' @param Dir
#'
#' @return
#' @export
#'
#' @examples
ExptoXlSX <- function(dds, Dir = "./") {
  NormData <- BiocGenerics::counts(dds, normalized = TRUE)
  writexl::write_xlsx(as.data.frame(NormData),
    path = paste(Dir, "RNASeq_NormExpression.xlsx", sep = ""), col_names = T, format_headers = T
  )
}


#' Title save different expression genes to excel files
#'
#' @param Diff_res
#' @param Dir default dir path: "2.Results/R_outputs/Tables/"
#'
#' @return
#' @export
#'
#' @examples
DEtoXlSX <- function(Diff_res, Dir = "./") {
  prefix <- paste(Diff_res$Group[1], Diff_res$Group[2], sep = "_VS_")
  writexl::write_xlsx(
    x = list(
      Up_regulated = as.data.frame(Diff_res$up_gene),
      Down_regulated = as.data.frame(Diff_res$down_gene)
    ),
    path = paste(Dir, prefix, "_DEGenes.xlsx", sep = ""), col_names = T, format_headers = T
  )
  writexl::write_xlsx(as.data.frame(Diff_res$Res),
    path = paste(Dir, prefix, "_AllGenes.xlsx", sep = ""), col_names = T, format_headers = T
  )
}

#' Title
#'
#' @param Diff_res
#' @param Dir
#' @param version
#'
#' @return
#' @export
#'
#' @examples
GOtoXlSX <- function(Diff_res, Dir = "./", version = "V3") {
  prefix <- paste(Diff_res$Group[1], Diff_res$Group[2], sep = "_VS_")
  up_go <- getGO(
    genes = Diff_res$up_gene$SYMBOL, version = version,
    universeIDs = names(Diff_res$rankGeneList)
  )$GO_df
  down_go <- getGO(
    genes = Diff_res$down_gene$SYMBOL, version = version,
    universeIDs = names(Diff_res$rankGeneList)
  )$GO_df
  if (is.null(up_go)) {
    print("There are no enrichment results of up regulated genes")
  } else {
    GOList <- list(UP_GO = up_go)
  }
  if (is.null(down_go)) {
    print("There are no enrichment results of down regulated genes")
  } else {
    GOList[["DOWN_GO"]] <- down_go
  }
  if (length(GOList) > 0) {
    writexl::write_xlsx(GOList, path = paste(Dir, prefix, "_GOterms_enrichment.xlsx", sep = ""), col_names = T, format_headers = T)
  }
}



#' Title
#'
#' @param Diff_res
#' @param Dir
#' @param pvalueCutoff
#' @param version
#'
#' @return
#' @export
#'
#' @examples
GSEAtoXlSX <- function(Diff_res, Dir = "./",
                       pvalueCutoff = 0.5,
                       version = "V3") {
  if (version == "V3") {
    OrgDB <- "org.Sjaponicum.eg.db"
  } else if (version == "V4") {
    OrgDB <- "org.SjaponicumV4.eg.db"
  } else {
    print("Please provide genome version V3 or V4!")
  }


  plotDir <- Dir
  prefix <- paste(Diff_res$Group[1], Diff_res$Group[2], sep = "_VS_")
  kk <- clusterProfiler::GSEA(
    geneList = Diff_res$rankGeneList,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME
  )
  bp <- clusterProfiler::gseGO(
    geneList = Diff_res$rankGeneList,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500, ont = "BP", OrgDb = OrgDB, keyType = "GENEID"
  )
  cc <- clusterProfiler::gseGO(
    geneList = Diff_res$rankGeneList,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500, ont = "CC", OrgDb = OrgDB, keyType = "GENEID"
  )

  mf <- clusterProfiler::gseGO(
    geneList = Diff_res$rankGeneList,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500, ont = "MF", OrgDb = OrgDB, keyType = "GENEID"
  )
  if (is.null(kk)) {
    print("There are no enrichment results of gsea KEGG")
  } else {
    p_kk <- enrichplot::gseaplot2(kk, geneSetID = 1:5, pvalue_table = TRUE, ES_geom = "dot")
    ggsave(plot = p_kk, filename = paste(plotDir, prefix, "_GSEA_KEGGpathway.pdf", sep = ""), width = 9, height = 6)
    writexl::write_xlsx(x = as.data.frame(kk), path = paste(Dir, prefix, "_GSEA_KEGGpathway.xlsx", sep = ""), col_names = T, format_headers = T)
  }
  if (is.null(bp)) {
    print("There are no results:Biological Process")
  } else {
    p_bp <- enrichplot::gseaplot2(bp, geneSetID = 1:5, pvalue_table = TRUE, ES_geom = "line")
    ggsave(plot = p_bp, filename = paste(plotDir, prefix, "_GSEA_GOterm_BP.pdf", sep = ""), width = 9, height = 6)
    bp <- as.data.frame(bp)
  }
  if (is.null(cc)) {
    print("There are no results:Cellular Component")
  } else {
    p_cc <- enrichplot::gseaplot2(cc, geneSetID = 1:5, pvalue_table = TRUE, ES_geom = "line")
    ggsave(plot = p_cc, filename = paste(plotDir, prefix, "_GSEA_GOterm_CC.pdf", sep = ""), width = 9, height = 6)
    cc <- as.data.frame(cc)
  }
  if (is.null(mf)) {
    print("There are no results:Molecular Function")
  } else {
    p_mf <- enrichplot::gseaplot2(mf, geneSetID = 1:5, pvalue_table = TRUE, ES_geom = "line")
    ggsave(plot = p_mf, filename = paste(plotDir, prefix, "_GSEA_GOterm_mf.pdf", sep = ""), width = 9, height = 6)
    mf <- as.data.frame(mf)
  }
  GOLIST <- plyr::compact(list(MF = mf, BP = bp, CC = cc))
  if (length(GOLIST) == 0) {
    print("There are no enrichment terms")
  } else {
    writexl::write_xlsx(x = GOLIST, path = paste(Dir, prefix, "_GSEA_GOterm.xlsx", sep = ""), col_names = T, format_headers = T)
  }
}
