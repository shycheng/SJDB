#' Title Plot gene expression bubbleplot / heatmap
#'
#' @param inp
#' @param inpGrp
#' @param inpPlt
#' @param inpsub1
#' @param inpH5
#' @param inpScl
#' @param inpRow
#' @param inpCol
#' @param inpcols
#' @param inpfsz
#' @param save
 #' @param shinycell_prefix
#'
#' @return
#' @export
#'
#' @examples
scBubbHeat <- function(shinycell_prefix,inp,inpH5, inpGrp, inpPlt = "Bubbleplot",
                       inpsub1=NULL,   inpScl=T, inpRow =T, inpCol=T,
                       inpcols = 'Blue-Yellow-Red',inpfsz = 'Small', save = FALSE){
  inpConf <- shinycell_prefix@conf
  inpMeta <- shinycell_prefix@meta
  # inpH5   <- shinycell_prefix@hdf5
  inpGene <- shinycell_prefix@gene
  cList   <- shinycell_prefix@cList
  sList   <- shinycell_prefix@sList

  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]}
  # Identify genes that are in our dataset
  geneList = inp[inp %in% names(inpGene)]

  # Prepare ggData
  h5file <- hdf5r::H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData = data.table::data.table()
  for(iGene in geneList){
    tmp = inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE]
    colnames(tmp) = c("sampleID", "sub")
    tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]]
    tmp$geneName = iGene
    tmp$val = h5data$read(args = list(inpGene[iGene], quote(expr=)))
    ggData = data.table::rbindlist(list(ggData, tmp))
  }
  h5file$close_all()

  # Aggregate
  ggData$val = expm1(ggData$val)
  ggData = ggData[, .(val = mean(val), prop = sum(val>0) / length(sampleID)),
                  by = c("geneName", "grpBy")]
  ggData$val = log1p(ggData$val)

  # Scale if required
  colRange = range(ggData$val)
  if(inpScl){
    ggData[, val:= scale(val), keyby = "geneName"]
    colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val))))
  }

  # hclust row/col if necessary
  ggMat = data.table::dcast.data.table(ggData, geneName~grpBy, value.var = "val")
  tmp = ggMat$geneName
  ggMat = as.matrix(ggMat[, -1])
  rownames(ggMat) = tmp
  ggMat <- na.omit(ggMat)
  ggData <- ggData[ggData$geneName %in% rownames(ggMat),]

  if(inpRow){
    hcRow = ggdendro::dendro_data(as.dendrogram(hclust(dist(ggMat))))
    ggRow = ggplot() + coord_flip() +
      geom_segment(data = hcRow$segments, aes(x=x,y=y,xend=xend,yend=yend)) +
      scale_y_continuous(breaks = rep(0, data.table::uniqueN(ggData$grpBy)),
                         labels = unique(ggData$grpBy), expand = c(0, 0)) +
      scale_x_continuous(breaks = seq_along(hcRow$labels$label),
                         labels = hcRow$labels$label, expand = c(0, 0.5)) +
      sctheme(base_size = sList[inpfsz]) +
      theme(axis.title = element_blank(), axis.line = element_blank(),
            axis.ticks = element_blank(), axis.text.y = element_blank(),
            axis.text.x = element_text(color="white", angle = 45, hjust = 1))
    ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label)
  } else {
    ggData$geneName = factor(ggData$geneName, levels = rev(geneList$gene))
  }
  if(inpCol){
    hcCol = dendro_data(as.dendrogram(hclust(dist(t(ggMat)))))
    ggCol = ggplot() +
      geom_segment(data = hcCol$segments, aes(x=x,y=y,xend=xend,yend=yend)) +
      scale_x_continuous(breaks = seq_along(hcCol$labels$label),
                         labels = hcCol$labels$label, expand = c(0.05, 0)) +
      scale_y_continuous(breaks = rep(0, data.table::uniqueN(ggData$geneName)),
                         labels = unique(ggData$geneName), expand=c(0,0)) +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      theme(axis.title = element_blank(), axis.line = element_blank(),
            axis.ticks = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(color = "white"))
    ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label)
  }

  # Actual plot according to plottype
  if(inpPlt == "Bubbleplot"){
    # Bubbleplot
    ggOut = ggplot(ggData, aes(grpBy, geneName, color = val, size = prop)) +
      geom_point() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      scale_x_discrete(expand = c(0.05, 0)) +
      scale_y_discrete(expand = c(0, 0.5)) +
      scale_size_continuous("proportion", range = c(0, 8),
                            limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) +
      scale_color_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) +
      guides(color = guide_colorbar(barwidth = 15)) +
      theme(axis.title = element_blank(), legend.box = "vertical")
  } else {
    # Heatmap
    ggOut = ggplot(ggData, aes(grpBy, geneName, fill = val)) +
      geom_tile() +
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +
      scale_x_discrete(expand = c(0.05, 0)) +
      scale_y_discrete(expand = c(0, 0.5)) +
      scale_fill_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) +
      guides(fill = guide_colorbar(barwidth = 15)) +
      theme(axis.title = element_blank())
  }

  # Final tidy
  ggLeg = g_legend(ggOut)
  ggOut = ggOut + theme(legend.position = "none")
  if(!save){
    if(inpRow & inpCol){ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),
                   layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
    } else if(inpRow){ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),
                   layout_matrix = rbind(c(1,3),c(2,NA)))
    } else if(inpCol){ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, ggCol, heights = c(1,7,2),
                   layout_matrix = rbind(c(3),c(1),c(2)))
    } else {ggOut =
      gridExtra::grid.arrange(ggOut, ggLeg, heights = c(7,2),
                   layout_matrix = rbind(c(1),c(2)))
    }
  } else {
    if(inpRow & inpCol){ggOut =
      gridExtra::arrangeGrob(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),
                  layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
    } else if(inpRow){ggOut =
      gridExtra::arrangeGrob(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),
                  layout_matrix = rbind(c(1,3),c(2,NA)))
    } else if(inpCol){ggOut =
      gridExtra::arrangeGrob(ggOut, ggLeg, ggCol, heights = c(1,7,2),
                  layout_matrix = rbind(c(3),c(1),c(2)))
    } else {ggOut =
      gridExtra::arrangeGrob(ggOut, ggLeg, heights = c(7,2),
                  layout_matrix = rbind(c(1),c(2)))
    }
  }
  return(ggOut)
}

#' Title a S4 object of shinycell_prefix
#'
#' @slot conf data.table.
#' @slot def list.
#' @slot gene integer.
#' @slot meta data.table.
#' @slot cList list.
#' @slot sList numeric.
#' @slot hdf5 character.
#'
#' @return
#' @export
#'
#' @examples
ShinyCell_prefix <- setClass(Class = "ShinyCell_prefix",slots = list(
  conf = "data.table",
  def = "list",
  gene = 'integer',
  meta = 'data.table',
  cList = 'list',
  sList = 'numeric',
  hdf5 = 'character'
))



#' Title Plot theme
#'
#' @param base_size
#' @param XYval
#' @param Xang
#' @param XjusH
#'
#' @return
#'
#' @examples
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){
  oupTheme = theme(
    text =             element_text(size = base_size, family = "Helvetica"),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line =   element_line(colour = "black"),
    axis.ticks =  element_line(colour = "black", size = base_size / 20),
    axis.title =  element_text(face = "bold"),
    axis.text =   element_text(size = base_size),
    axis.text.x = element_text(angle = Xang, hjust = XjusH),
    legend.position = "bottom",
    legend.key =      element_rect(colour = NA, fill = NA)
  )
  if(!XYval){
    oupTheme = oupTheme + theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  return(oupTheme)
}



#' Title Function to extract legend
#'
#' @param a.gplot
#'
#' @return
#'
#' @examples
g_legend <- function(a.gplot){
  tmp <- ggplot2::ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}




