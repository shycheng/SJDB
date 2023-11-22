#' Title switch ID '-' or '_'
#'
#' @param ID
#'
#' @return
#' @export
#'
#' @examples
revID <- function(ID) {
  if (sum(stringr::str_count(ID, "-")) > 0) {
    id <- gsub("-", "_", ID)
  } else {
    id <- gsub("_", "-", ID)
  }
  return(id)
}


.check_phenoTab <- function(pheno) {
  required_cols <- c("Sample", "Sex", "Condition")
  if (all(required_cols %in% colnames(pheno))) {
    message("Pheno info looks good.")
  } else {
    missing_cols <- required_cols[!required_cols %in% colnames(pheno)]
    stop("Missing required columns in pheno data frame: ", paste(missing_cols, collapse = ", "))
  }
}


# 定义Assay类
setClass(
  "Assay",
  slots = c(
    counts = "data.frame",
    tpm = "data.frame",
    deseq2 = "data.frame"
  )
)

# 定义Assay的构造函数
Assay <- function(counts) {
  # 验证counts是否为data.frame
  if (!is(counts, "data.frame")) {
    stop("'counts' must be a data.frame")
  }
  # 创建Assay对象，其中tpm和deseq2是空的data.frames
  new("Assay", counts = counts, tpm = data.frame(), deseq2 = data.frame())
}


#' Title 定义一个名为Person的S4类，它有三个属性：Assay 类 , pheno 和 outDir
#'
#' @slot assay Assay.
#' @slot pheno data.frame.
#' @slot outDir character.
#' @slot dds DESeqDataSet.
#' @slot contrast_df data.frame.
#' @slot diffres list.
#'
#' @return
#' @export
#'
#' @examples
setClass(
  "RNASeqOBJ",
  slots = c(
    assay = "Assay", # 使用Assay类作为slot
    pheno = "data.frame",
    outDir = "character",
    dds = "ANY",
    contrast_df = "data.frame",
    diffres = "list"
  ),
  prototype = list(
    dds = NULL # 设置默认值为NULL
  )
)

.createDirectories <- function(outDir) {
  dirsToCreate <- c("Summary")
  for (dirName in dirsToCreate) {
    dirPath <- file.path(outDir, dirName)
    if (!dir.exists(dirPath)) {
      dir.create(dirPath, recursive = T)
    }
  }
}


#' Title create a RNAseq object for run RNAseq pipeline
#'
#' @param counts raw counts of RNA-seq data
#' @param pheno a dataframe including Sample, Condition and other information
#' @param outDir output path
#'
#' @return
#' @export
#'

RNASeqOBJ <- function(counts, pheno = data.frame(), outDir = "./") {
  # Check for pheno validity
  .check_phenoTab(pheno)

  # Create directories
  message(sprintf("Files will be saved to %s", outDir))
  .createDirectories(outDir)

  # Create Assay object
  # NOTE: You should define how the Assay function handles NA counts
  assay <- Assay(counts)

  # Create RNASeqOBJ object
  new("RNASeqOBJ",
    assay = assay,
    pheno = pheno,
    outDir = outDir,
    dds = NULL, # 允许dds在初始时为NULL
    contrast_df = data.frame(),
    diffres = list()
  )
}


.validRNASeqOBJ <- function(RNASeqOBJ = NULL) {
  stopifnot(!is.null(RNASeqOBJ))
  if (inherits(RNASeqOBJ, "RNASeqOBJ")) {
    return(RNASeqOBJ)
  } else {
    stop("Cannot validate RNASeqOBJ options are a valid RNASeqOBJ")
  }
}







#' Title TPM normalized of raw counts
#'
#' @param object
#' @param version
#' @param gtf_file
#'
#' @return
#' @export
#'
#' @examples
fill_TPM <- function(object = NULL,
                     version = "V3",
                     gtf_file = NULL) {
  .validRNASeqOBJ(object)

  # TPM 标准化
  ### 获取基因非重叠外显子长度
  if (version == "V3") {
    txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtf_file)
  } else if (version == "V4") {
    txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtf_file)
  } else {
    print("Error,please provide genome version V3 or V4!")
  }


  # then collect the exons per gene id
  exons.list.per.gene <- GenomicFeatures::exonsBy(txdb, by = "gene")

  # then for each gene, reduce all the exons to a set of non overlapping exons,
  # calculate their lengths (widths) and sum then
  exonic.gene.sizes <- sum(width(GenomicRanges::reduce(exons.list.per.gene)))

  ### TPM标准化公式
  r_tpm <- function(dfr, len) {
    dfr1 <- sweep(dfr, MARGIN = 1, (len / 10^4), `/`)
    scf <- colSums(dfr1) / (10^6)
    return(sweep(dfr1, 2, scf, `/`))
  }

  TPM <- r_tpm(object@assay@counts, exonic.gene.sizes[rownames(object@assay@counts)])
  object@assay@tpm <- TPM

  TPM$Annotation <- SJDB::id2name[rownames(object@assay@counts), "Note"]

  write.csv(TPM,
    file = file.path(object@outDir, "Summary", "Expression_TPM_addAnnotation.csv"),
    quote = F,
    row.names = T
  )

  object
}


#' Title plot_PCA using one of assay in RNAseq object
#'
#' @param object
#' @param assay
#' @param show
#'
#' @return
#' @export
#'
#' @examples
plot_PCA <- function(object = NULL,
                     assay = "tpm",
                     show = TRUE) {
  .validRNASeqOBJ(object)
  if (assay == "counts") {
    mat <- object@assay@counts
  } else if (assay == "tpm") {
    mat <- object@assay@tpm
  } else if (assay == "deseq2") {
    mat <- object@assay@deseq2
  } else {
    stop("'Please choose a assay to PCA'")
  }
  pca <- FactoMineR::PCA(t(mat), scale.unit = TRUE, ncp = 5)
  pca_df <- data.frame(pca$ind$coord)

  pca_df$condition <- object@pheno$Condition
  pca_df$Treat <- object@pheno$Treat
  pca_df$name <- rownames(pca_df)
  pca_plot <- ggpubr::ggscatter(pca_df, "Dim.1", "Dim.2", color = "condition", size = 5, label = "name", repel = T) +
    xlab(paste0("PC1: ", round(pca$eig[, 2][1], digits = 2), "% variance")) +
    ylab(paste0("PC2: ", round(pca$eig[, 2][2], digits = 2), "% variance")) +
    theme_bw(base_size = 18, base_line_size = 1) +
    ggtitle("PCA")

  ggplot2::ggsave(
    filename = file.path(object@outDir, "Summary", sprintf("%s_pca.pdf", assay)),
    plot = pca_plot,
    width = 12,
    height = 8
  )

  if (show) {
    pca_plot
  }
}

#' Title plot correlation coefficient heatmap of samples
#'
#' @param object
#' @param assay
#' @param method
#' @param cluster_rows
#' @param cluster_cols
#' @param show
#'
#' @return
#' @export
#'
#' @examples
plot_corHeatmap <- function(object = NULL,
                            assay = "tpm",
                            method = "pearson",
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            show = TRUE) {
  .validRNASeqOBJ(object)

  if (assay == "counts") {
    mat <- object@assay@counts
  } else if (assay == "tpm") {
    mat <- object@assay@tpm
  } else if (assay == "deseq2") {
    mat <- object@assay@deseq2
  } else {
    stop("'Please choose a assay to plot_corHeatmap'")
  }
  cor_matrix <- cor(mat,
    method = method
  )
  cor_heatmap <- pheatmap::pheatmap(cor_matrix,
    show_rownames = T, show_colnames = T,
    cluster_cols = cluster_cols,
    cluster_rows = cluster_rows,
    main = sprintf("Samples' %s correlation coefficient", method),
    fontsize = 15, display_numbers = F, fontsize_number = 12
  )

  ggsave(
    filename = file.path(object@outDir, "Summary", sprintf("%s_%s_corHeatmap.pdf", assay, method)),
    plot = cor_heatmap,
    width = 16,
    height = 16
  )

  if (show) {
    cor_heatmap
  }
}


#' Title
#'
#' @param object
#' @param design
#' @param force
#'
#' @return
#' @export
#'
#' @examples
run_DESeq2 <- function(object = NULL, design = ~Condition, force = FALSE) {
  .validRNASeqOBJ(object)
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("The DESeq2 package is required but not installed. Please install it using BiocManager::install('DESeq2').")
  }

  # 准备DESeqDataSet对象
  .prepareDESeqDataSet <- function(object, design) {
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = object@assay@counts,
      colData = object@pheno,
      design = design
    )
    dds <- dds[rowMeans(counts(dds)) >= 5, ]
    return(dds)
  }

  # 执行DESeq2分析
  .performDESeq2 <- function(dds) {
    dds <- DESeq2::DESeq(dds)
    return(dds)
  }

  if (is.null(object@dds) || force) {
    if (force) message("Re-running DESeq2 analysis.")
    dds <- .prepareDESeqDataSet(object, design)
    dds <- .performDESeq2(dds)
    object@dds <- dds
    object@assay@deseq2 <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
  } else {
    message("DESeq2 analysis has already been run. Use 'force = TRUE' to re-run.")
  }

  return(object)
}



#' Title
#'
#' @param object
#' @param control
#'
#' @return
#' @export
#'
#' @examples
createContrastDataFrame <- function(object = NULL,
                                    control = "dsGFP") {
  .validRNASeqOBJ(object)

  pheno <- object@pheno

  # 获取非控制组和控制组的条件
  treatmentConditions <- unique(pheno$Condition[pheno$Treat != control])
  controlConditions <- unique(pheno$Condition[pheno$Treat == control])

  # 如果控制组只有一个元素，但处理组有多个，重复控制组元素
  if (length(treatmentConditions) != length(controlConditions) && length(controlConditions) == 1) {
    controlConditions <- rep(controlConditions, length(treatmentConditions))
  }

  # 创建对比矩阵
  contrast_df <- data.frame(group.x = sort(treatmentConditions), group.y = sort(controlConditions))
  object@contrast_df <- contrast_df

  return(object)
}


#' Title
#'
#' @param object
#' @param p_cutoff
#' @param FC_cutoff
#'
#' @return
#' @export
#'
#' @examples
processContrasts <- function(object = NULL,
                             p_cutoff = 0.05,
                             FC_cutoff = log2(1.5)) {
  .validRNASeqOBJ(object)
  contrast_df <- object@contrast_df
  dds <- object@dds

  # 定义内部函数
  get_res <- function(i) {
    get_mRNA_diff(
      dds = dds,
      group = "Condition",
      x = contrast_df$group.x[i],
      y = contrast_df$group.y[i],
      p_cutoff = p_cutoff,
      FC_cutoff = FC_cutoff
    )
  }

  # 使用lapply遍历对比矩阵
  diff_res_list <- lapply(seq_len(nrow(contrast_df)), get_res)
  names(diff_res_list) <- paste(contrast_df$group.x, contrast_df$group.y, sep = "_vs_")

  object@diffres <- diff_res_list
  return(object)
}



run_RNASeq_pipeline <- function(object = NULL,
                                version = "V3") {
  .validRNASeqOBJ(object)


  # 在此方法中，使用lapply来迭代contrast_df中的行并处理每个对照

  run_RNASeq_pipeline <- function(root_dir, version = version,
                                  diff_res) {
    prefix <- paste(as.character(diff_res$Group)[1],
      as.character(diff_res$Group)[2],
      sep = "_vs_"
    ) # eg: C24M_vs_C12Mix
    out_dir <- file.path(root_dir, prefix, "/")
    if (dir.exists(out_dir) == FALSE) {
      dir.create(out_dir)
    }
    tryCatch(
      expr = {
        DEtoXlSX(Diff_res = diff_res, Dir = out_dir)
        GSEAtoXlSX(Diff_res = diff_res, Dir = out_dir, version = version)
        GOtoXlSX(Diff_res = diff_res, Dir = out_dir, version = version)
        Plot_DE_GOTerms(Diff_res = diff_res, DIR = out_dir, version = version)
        PlotKEGG(Diff_res = diff_res, Dir = out_dir, save_Plot = T)
        Plot_volcano(Diff_res = diff_res, Dir = out_dir, save_Plot = T)
        Plot_enhancedVolcano(Diff_res = diff_res, Dir = out_dir, save_Plot = T)
        Plot_heatmap(Diff_res = diff_res, dds = dds, Dir = out_dir, save_Plot = T)
      }, error = function(e) {
        # 当发生异常时，打印错误信息并返回 NULL
        cat("Error in run_RNASeq_pipeline for", prefix, ":", e$message, "\n")
        return(NULL)
      }, finally = print(paste0("Finish analysis:", prefix, "  Files were save to ", out_dir))
    )
  }
  lapply(object@diffres, run_RNASeq_pipeline, root_dir = object@outDir)
  invisible(object) # 修改为 invisible，通常情况下不需要打印整个对象
}
