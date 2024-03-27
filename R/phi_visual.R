#' @title phiPlot
#' @description plot mean psi vs. phi scatter plot for each exon, colored by the confidence
#' interval of phi
#' @param data the dataframe of phi, which should include columns for mean psi, phi, phi confidence
#' interval and the exon name
#' @param psi_range Only exons with mean psi within this range will be annotated
#' @param phi_conf_upr Only exons with phi confidence interval lower than this threshold
#' will be annotated
#' @param psi_col,phi_col,phi_conf_col,annot_col The colnames to indicate the mean psi,phi,phi confidence
#' interval and exon annotation within the input data
#' @param annot_size the text size for exon annotation
#' @import ggplot2 RColorBrewer
#' @importFrom latex2exp TeX
#' @importFrom ggrepel geom_text_repel
#' @export
phiPlot = function(data,psi_range = c(0.3,0.7),phi_conf_upr = 0.2,
                   psi_col = "mean_psi",phi_col = "phi",
                   phi_conf_col = "phi_conf",annot_col = "symbol",
                   annot_size = 4){
  data_mid <- data[data[,psi_col] > psi_range[1] & data[,psi_col] < psi_range[2] &
                     data[,phi_conf_col] < phi_conf_upr,]
  out = ggplot()+
    geom_point(data = data,
               aes_string(x = psi_col,y = phi_col,color = phi_conf_col),size = 2)+
    scale_color_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) +
    ggrepel::geom_text_repel(data=data_mid, aes_string(x = psi_col,y = phi_col,label = annot_col),
                             size =annot_size)+
    xlab(TeX("$\\bar{\\psi}$"))+ylab(TeX("$\\phi$"))+
    labs(color=TeX('$CI\\,of\\,\\phi$'))+
    theme_classic()+
    theme(text = element_text(size = 15))
  return(out)
}


#' @title plot psi distribution
#' @description plot psi distribution for a cell group into histogram
#' @param spliceOb a Splice object
#' @param gene the name of specified gene
#' @param site the selected site to be plot
#' @param cells the cells to be included in the psi plot
#' @param alpha the parameter alpha for fitted beta distribution for psi
#' @param beta the parameter beta for fitted beta distribution for psi
#' @param exprs_thresh the minimum gene expression required for psi calculation
#' @param binwidth the bin size for the histogram ranging from (0,1]
#' @param low,high,mid Colors for the bins to represent the mean gene count in each bin
#' @param midpoint the gene count thresh to set the middle color
#' @param ... aesthetic parameters for ggplot
#' @return a ggplot object of psi distribution histogram
#' @import ggplot2
#' @importFrom latex2exp TeX
#' @importFrom dplyr group_by summarise
#' @importFrom stats dbeta median na.omit
psiHist.base <- function(spliceOb, gene, site, cells = "all",
                         alpha = NULL,beta = NULL,
                         exprs_thresh = 5,
                         binwidth = 0.05,low = "steelblue",
                         high = "coral",mid = "whitesmoke",
                         midpoint = NULL,binsize = 8){
  count_list = getMetasites.base(spliceOb,gene)
  gene_count_table = count_list@cellGeneCount
  exon_count_table = count_list@cellSiteCount
  metasite = count_list@sites

  if(cells[1] != "all"){
    allcells = rownames(gene_count_table)
    cells = cells[cells %in% allcells]
    gene_count_table = gene_count_table[cells,]
    exon_count_table = exon_count_table[cells,]
  }

  if(site %in% colnames(gene_count_table)){
    group_exon = exon_count_table[,site]
    group_gene = gene_count_table[,site]
  }
  else{
    id = sapply(metasite,function(x) site %in% x)
    if(sum(id) == 0){
      stop("The exon doesn't exist in this gene")
    }
    else{
      group_exon = exon_count_table[,id]
      group_gene = gene_count_table[,id]
    }
  }

  psi = unlist(group_exon/group_gene)

  plot_data = na.omit(as.data.frame(cbind(group_gene,group_exon,psi)))
  colnames(plot_data) = c("gene","exon","psi")
  plot_data = plot_data[plot_data$gene >= exprs_thresh,]
  plot_data$psi_bin = cut(plot_data$psi,seq(0,1,binwidth),include.lowest = TRUE)

  plot_data = plot_data %>%
    group_by(psi_bin) %>%
    summarise(gene = mean(gene),count = n())
  plot_data$freq = plot_data$count/sum(plot_data$count)
  plot_data$freq = plot_data$freq/binwidth

  plot_data$psi_bin = as.numeric(plot_data$psi_bin)
  plot_data$psi_bin = plot_data$psi_bin/max(plot_data$psi_bin)
  plot_data$psi_bin = plot_data$psi_bin-min(plot_data$psi_bin)/2

  if(is.null(midpoint)){
    midpoint = median(plot_data$gene)
  }
  out = ggplot()+
    geom_segment(data = na.omit(plot_data),
                 aes(x = psi_bin,xend = psi_bin,y = 0,
                     yend = freq,color = gene),size = binsize)+
    scale_color_gradient2(low = low,high = high,
                          mid = mid,midpoint = midpoint,name = "gene count")+
    #scale_color_gradientn(name = "gene count",colors = rev(brewer.pal(11, "RdYlBu")))+
    xlab(TeX("$\\psi$"))+ylab("Density")+
    theme_classic()+
    theme(text = element_text(size = 15))

  if(!is.null(alpha) & !is.null(beta)){
    group_density = dbeta(seq(0,1,0.01),alpha,beta)
    out = out+geom_line(aes(x = seq(0,1,0.01),y = group_density))
  }

  return(out)
}

#' @title generic psiHist function definition
#' @param object the Splice or Seurat object
#' @inheritParams psiHist.base
#' @param ... parameters for psiHist.base
#' @export
setGeneric("psiHist",
           function(object,gene,site,...) standardGeneric("psiHist"))

#' @title generic psiHist function for Splice object
#' @description  plot psi distribution for a cell group into histogram for a Splice object
#' @param object The Splice object
#' @inheritParams psiHist.base
#' @return a ggplot object of psi distribution histogram
#' @export
setMethod("psiHist",
          signature(object = "Splice",gene = "character",site = "character"),
          function(object, gene, site,...){
            psiHist.base(spliceOb = object,gene = gene,
                         site = site,...)
          }
)


#' @title psiCellPlot.base
#' @description plot psi value in cell embeddings
#' @param spliceOb a Splice object
#' @param gene the name of specified gene
#' @param sites the selected exon to be plot
#' @param cell_embedding A dataframe of cell embedding coordinates
#' @param dims Two dimensions to use to plot cell embedding
#' @param cells the cells to be included in the psi plot
#' @param exprs_thresh the minimum gene expression required for psi calculation
#' @param ... parameters for gene_exons_psi.base to calculate psi value
#' @return a ggplot object of psi in cell embedding
#' @import ggplot2
#' @importFrom latex2exp TeX
#' @importFrom stats median na.omit
#' @export
psiCellPlot.base = function(spliceOb,gene,sites,cell_embedding,dims = c(1,2),
                            cells = "all",exprs_thresh = 10,...){
  cell_embedding = as.data.frame(cell_embedding)
  dims = colnames(cell_embedding)[dims]

  psi = geneSitesPsi.base(spliceOb = spliceOb,gene = gene,
                          cells = cells,exprs_thresh = exprs_thresh,...)

  if(sum(colnames(psi) %in% sites) == 0){
    stop("The selected sites don't exist, maybe be filtered out due to low expression!")
  }

  psi = psi[,colnames(psi) %in% sites,drop=FALSE]

  cell_embedding = cbind(cell_embedding,psi[rownames(cell_embedding),,drop=FALSE])
  #return(cell_embedding)

  out_plot = lapply(colnames(psi),function(x){
    out = ggplot()+
      geom_point(data = cell_embedding,
                 aes_string(x = dims[1],y = dims[2]),
                 color = "grey",alpha = 1,size = 1)+
      geom_point(data = na.omit(cell_embedding),
                 aes_string(x = dims[1],y = dims[2],color = x),
                 alpha = 1,size = 2)+
      #scale_color_gradient2(low = "steelblue",high = "coral",
      #                      mid = "whitesmoke",midpoint = median(psi[,x],na.rm = TRUE))+
      scale_color_viridis_c(name = TeX("$\\psi$"))+
      theme_classic()+
      xlab(dims[1])+ylab(dims[2])+
      theme(text = element_text(size = 15))
    return(out)
  })
  return(out_plot)
}


#' @title generic psiCellPlot function definition
#' @param object the Splice or Seurat object
#' @inheritParams psiCellPlot.base
#' @param ... parameters for psiCellPlot.base
#' @export
setGeneric("psiCellPlot",
           function(object,gene,exons,cell_embedding,...)
             standardGeneric("psiCellPlot"))

