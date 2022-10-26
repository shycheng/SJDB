 #' Title Plot gene expression bubbleplot / heatmap
#'
#' @param inpConf
#' @param inpMeta
#' @param inp
#' @param inpGrp
#' @param inpPlt
#' @param inpsub1
#' @param inpH5
#' @param inpGene
#' @param inpScl
#' @param inpRow
#' @param inpCol
#' @param inpcols
#' @param inpfsz
#' @param save
#'
#' @return
#' @export
#'
#' @examples
scBubbHeat <- function(shinycell_prefix,inp, inpGrp, inpPlt = "Bubbleplot",
                       inpsub1=NULL,   inpScl=T, inpRow =T, inpCol=T,
                       inpcols = 'Blue-Yellow-Red',inpfsz = 'Small', save = FALSE){
  inpConf <- shinycell_prefix@conf
  inpMeta <- shinycell_prefix@meta
  inpH5   <- shinycell_prefix@hdf5
  inpGene <- shinycell_prefix@gene
  cList   <- shinycell_prefix@cList
  sList   <- shinycell_prefix@sList

  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]}
  # Identify genes that are in our dataset
  geneList = inp[inp %in% names(inpGene)]

  # Prepare ggData
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  ggData = data.table()
  for(iGene in geneList){
    tmp = inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE]
    colnames(tmp) = c("sampleID", "sub")
    tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]]
    tmp$geneName = iGene
    tmp$val = h5data$read(args = list(inpGene[iGene], quote(expr=)))
    ggData = rbindlist(list(ggData, tmp))
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
  ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val")
  tmp = ggMat$geneName
  ggMat = as.matrix(ggMat[, -1])
  rownames(ggMat) = tmp
  ggMat <- na.omit(ggMat)
  ggData <- ggData[ggData$geneName %in% rownames(ggMat),]

  if(inpRow){
    hcRow = dendro_data(as.dendrogram(hclust(dist(ggMat))))
    ggRow = ggplot() + coord_flip() +
      geom_segment(data = hcRow$segments, aes(x=x,y=y,xend=xend,yend=yend)) +
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$grpBy)),
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
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$geneName)),
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
      grid.arrange(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),
                   layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
    } else if(inpRow){ggOut =
      grid.arrange(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),
                   layout_matrix = rbind(c(1,3),c(2,NA)))
    } else if(inpCol){ggOut =
      grid.arrange(ggOut, ggLeg, ggCol, heights = c(1,7,2),
                   layout_matrix = rbind(c(3),c(1),c(2)))
    } else {ggOut =
      grid.arrange(ggOut, ggLeg, heights = c(7,2),
                   layout_matrix = rbind(c(1),c(2)))
    }
  } else {
    if(inpRow & inpCol){ggOut =
      arrangeGrob(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),
                  layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))
    } else if(inpRow){ggOut =
      arrangeGrob(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),
                  layout_matrix = rbind(c(1,3),c(2,NA)))
    } else if(inpCol){ggOut =
      arrangeGrob(ggOut, ggLeg, ggCol, heights = c(1,7,2),
                  layout_matrix = rbind(c(3),c(1),c(2)))
    } else {ggOut =
      arrangeGrob(ggOut, ggLeg, heights = c(7,2),
                  layout_matrix = rbind(c(1),c(2)))
    }
  }
  return(ggOut)
}


# # create a prefix class
# sc1conf = readRDS("~/Projects/Rscripts_set/SJDB/inst/extdata/SingleCell_RDS/Shinydir/shiny_data/sc1conf.rds")
# sc1def  = readRDS("~/Projects/Rscripts_set/SJDB/inst/extdata/SingleCell_RDS/Shinydir/shiny_data/sc1def.rds")
# sc1gene = readRDS("~/Projects/Rscripts_set/SJDB/inst/extdata/SingleCell_RDS/Shinydir/shiny_data/sc1gene.rds")
# sc1meta = readRDS("~/Projects/Rscripts_set/SJDB/inst/extdata/SingleCell_RDS/Shinydir/shiny_data/sc1meta.rds")
#
# sc2conf = readRDS("~/Projects/Rscripts_set/SJDB/inst/extdata/SingleCell_RDS/Shinydir/shiny_data/sc2conf.rds")
# sc2def  = readRDS("~/Projects/Rscripts_set/SJDB/inst/extdata/SingleCell_RDS/Shinydir/shiny_data/sc2def.rds")
# sc2gene = readRDS("~/Projects/Rscripts_set/SJDB/inst/extdata/SingleCell_RDS/Shinydir/shiny_data/sc2gene.rds")
# sc2meta = readRDS("~/Projects/Rscripts_set/SJDB/inst/extdata/SingleCell_RDS/Shinydir/shiny_data/sc2meta.rds")
#
# ### Useful stuff
# # Colour palette
# cList = list(c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84",
#                "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"),
#              c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF",
#                "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)],
#              c("#FDE725","#AADC32","#5DC863","#27AD81","#21908C",
#                "#2C728E","#3B528B","#472D7B","#440154"))
# names(cList) = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")
#
#
# # Panel sizes
# pList = c("400px", "600px", "800px")
# names(pList) = c("Small", "Medium", "Large")
# pList2 = c("500px", "700px", "900px")
# names(pList2) = c("Small", "Medium", "Large")
# pList3 = c("600px", "800px", "1000px")
# names(pList3) = c("Small", "Medium", "Large")
# sList = c(18,24,30)
# names(sList) = c("Small", "Medium", "Large")
# lList = c(5,6,7)
# names(lList) = c("Small", "Medium", "Large")
#
#
# ShinyCell_prefix <- setClass(Class = "ShinyCell_prefix",slots = list(
#   conf = "data.table",
#   def = "list",
#   gene = 'integer',
#   meta = 'data.table',
#   cList = 'list',
#   sList = 'numeric',
#   hdf5 = 'character'
# ))
#
#
# ShinyCell_prefix_male <- ShinyCell_prefix(  conf = sc1conf,
#                                             def  = sc1def ,
#                                             gene = sc1gene,
#                                             meta = sc1meta,
#                                             cList = cList,
#                                             sList = sList,
#                                             hdf5 = "~/Projects/Rscripts_set/SJDB/inst/extdata/SingleCell_RDS/Shinydir/shiny_data/sc1gexpr.h5")
# ShinyCell_prefix_female <- ShinyCell_prefix(conf = sc2conf,
#                                             def = sc2def ,
#                                             gene = sc2gene,
#                                             meta = sc2meta,
#                                             cList = cList,
#                                             sList = sList,
#                                             hdf5 = "~/Projects/Rscripts_set/SJDB/inst/extdata/SingleCell_RDS/Shinydir/shiny_data/sc2gexpr.h5")
#
#













