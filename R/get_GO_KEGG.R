#' Title get geneinfo
#'
#' @param geneid
#'
#' @return
#' @export
#'
#' @examples
getGeneInfo <- function(geneid){
  GeneInfo <- data.frame(ID=geneid,GeneName=id2name[geneid,"Gene_Name"],
                         Note=id2name[geneid,"Note"])
  return(GeneInfo)
}

#' Title
#'
#' @param genes
#' @param universeIDs
#' @param qCutoff
#'
#' @return
#' @export
#'
#' @examples
getGO <- function(genes,universeIDs,qCutoff=0.05){
  BP <- clusterProfiler::enrichGO(
    gene = genes,
    OrgDb = "org.Sjaponicum.eg.db",
    keyType = "GENEID",
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = universeIDs ,
    qvalueCutoff = qCutoff,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE
  )
  CC <- clusterProfiler::enrichGO(
    gene = genes,
    OrgDb = "org.Sjaponicum.eg.db",
    keyType = "GENEID",
    ont = "CC",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = universeIDs ,
    qvalueCutoff = qCutoff,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
  )
  MF <- clusterProfiler::enrichGO(
    gene = genes,
    OrgDb = "org.Sjaponicum.eg.db",
    keyType = "GENEID",
    ont = "MF",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = universeIDs ,
    qvalueCutoff = qCutoff,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE
  )

  if (is.null(CC) == FALSE) {
    if (dim(CC)[1] > 0 ) {
      CC_df <- base::as.data.frame(CC)
      CC_df$type <-"Cellular Component"
    }else{
      CC_df <- NULL
    }
  }
  if (is.null(MF) == FALSE) {
    if (dim(MF)[1] > 0 ) {
      MF_df <- base::as.data.frame(MF)
      MF_df$type <-"Molecular Function"
    }else{
      MF_df <- NULL
    }
  }
  if (is.null(BP) == FALSE) {
    if (dim(BP)[1] > 0 ) {
      BP_df <- base::as.data.frame(BP)
      BP_df$type <-"Biological Process"
    }else{
      BP_df <- NULL
    }
  }
  GO_df <- base::rbind(CC_df,BP_df,MF_df)
  GOList = list(BP =BP,CC = CC,MF = MF)
  if (is.null(BP) == TRUE){
    print("There are no significant enrichment GO Terms，try to increase p-value")
  }
  return(list(GO_df = GO_df,GOList = GOList))
}


#' Title
#'
#' @param geneid
#' @param universeIDs
#'
#' @return
#' @export
#'
#' @examples
getKEGG <- function(geneid,universeIDs){
  KEGG_tab <- clusterProfiler::enricher(geneid,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       universe = universeIDs,
                       minGSSize = 10,
                       maxGSSize = 500,
                       qvalueCutoff = 0.2,
                       TERM2GENE,
                       TERM2NAME = TERM2NAME)
  return(KEGG_tab)
}
