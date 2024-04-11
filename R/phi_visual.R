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
  count_list = getMetaSites(spliceOb,gene,cells,verbose = FALSE)
  gene_count_table = count_list@cellGeneCount
  exon_count_table = count_list@cellSiteCount
  metasite = count_list@sites

  if(cells[1] != "all"){
    allcells = rownames(gene_count_table)
    cells = cells[cells %in% allcells]
    gene_count_table = gene_count_table[cells,,drop = FALSE]
    exon_count_table = exon_count_table[cells,,drop = FALSE]
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

  plot_data = plot_data %>% mutate(psi_bin = as.numeric(psi_bin)) %>%
    mutate(psi_bin = binwidth*(psi_bin-1) + binwidth/2)

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


#' @title psiCellPlot
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
psiCellPlot = function(spliceOb,gene,sites,cell_embedding,dims = c(1,2),
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

#' @title isoform2bases
#' @description transform an isoform sequence into a vector of base this isoform covered
#' @param isoform The isoform sequence, should be formatted as "s1,e1|s2,e2|...", s means the start site of an exon
#' while e means the end site.
#' @param count The count for each base
#' @param sep The character to split the start and end sites of an exon.
#' @param split The character to split the exon in the isoform sequence.
#' @return a vector of bases
isoform2bases <- function(isoform,count = 1,
                          sep = ",",split = "|"){
  exons = unlist(strsplit(isoform,split = split,fixed = TRUE))

  bases = lapply(exons,function(x){
    x = unlist(strsplit(x,split = sep,fixed = TRUE))
    left = as.numeric(x[1])
    right = as.numeric(x[2])
    return(c(left:right))
  })
  bases = rep(unlist(bases),count)
  return(bases)
}

#' @title junction_count
#' @description collect the splice junction in an isoform sequence
#' @inheritParams isoform2bases
#' @param polyA A flag to indicate if the read has a polyA tail
#' @return A matrix recording the junction count, the first column is the junction
#' identity and the second column is the count
junction_count <- function(isoform,count = 1,polyA = 0,
                           sep = ",",split = "|"){
  junctions = unlist(strsplit(isoform,split = sep,fixed = TRUE))
  end = junctions[length(junctions)]
  if(length(junctions) > 2){
    junctions = junctions[2:(length(junctions)-1)]
    junc_count = lapply(junctions,function(x){
      x = unlist(strsplit(x,split = split,fixed = TRUE))
      return(x)
    })
    junc_count = do.call(rbind,junc_count)
    junc_count = cbind(junc_count,count)
  }
  else{
    junc_count = NULL
  }
  if(polyA > 0.5){
    end_count = c(end,"polyA",count)
  }
  else{
    end_count = NULL
  }
  junc_count = rbind(junc_count,end_count)
  return(junc_count)
}

#' @title sashimi_plot_data
#' @description transform the isoform count data to the format for sahimi plot
#' @inheritParams junction_count
#' @return A list recording the base count and the junction count.
sashimi_plot_data <- function(isoforms,counts,polyA = NULL,
                              sep = ",",split = "|"){
  if(length(isoforms) != length(counts)){
    stop("The size of isoforms and their counts don't match!")
  }
  if(!is.null(polyA)){
    if(length(polyA) != length(isoforms)){
      stop("The size of isoforms and their polyA don't match!")
    }
  }
  else{
    polyA = rep(0,length(isoforms))
  }

  iso_count = as.data.frame(cbind(isoforms,counts))
  colnames(iso_count) = c("isoform","count")
  iso_count$count = as.numeric(iso_count$count)
  iso_count = iso_count %>%
    group_by(isoform) %>%
    summarise(count = sum(count))
  bases = lapply(1:nrow(iso_count),function(i){
    isoform = iso_count$isoform[i]
    count = iso_count$count[i]
    sub_base = isoform2bases(isoform = isoform,count = count,
                             sep = sep,split = split)
    return(sub_base)
  })
  bases = as.data.frame(table(unlist(bases)))
  colnames(bases) = c("chr","coverage")
  bases$chr = as.numeric(as.character(bases$chr))

  junctions = lapply(1:length(isoforms),function(i){
    sub_junc = junction_count(isoform = isoforms[i],count = counts[i],
                              polyA = polyA[i],sep = sep,split = split)
    return(sub_junc)
  })
  junctions = as.data.frame(do.call(rbind,junctions))
  colnames(junctions) = c("start","end","count")
  junctions$count = as.numeric(junctions$count)
  junctions = junctions %>%
    group_by(start,end) %>%
    summarise(count = sum(count),.groups = "drop")

  return(list(bases,junctions))
}

#' @title range_scale
#' @description scale the data into a fixed range
#' @param data a numeric vector
#' @param lwr,upr the lower/upper bound of the range
#' @return a numeric vector recording the scaled data
range_scale = function(data,lwr = 1,upr = 5){
  diff = max(data)-min(data)
  unit = (upr-lwr)/diff
  scale_data = (data-min(data))*unit+lwr
  return(scale_data)
}


#' @title sashimi_plot
#' @description generate the sashimi plot for the isoform
#' @param coverage The base count
#' @param junction The junction count
#' @param filter_ratio Threshold to filter out lowly covered base
#' @param color_id The prdefined color set to use, currently there are two sets
#' @param region The xlab annotation for the sashimi plot
#' @inheritParams range_scale
#' @param color_set a vector to record user defined color set to use, should include two colors, one for the
#' base coverage and the other for the junction connection.
#' @import ggplot2
#' @importFrom geomtextpath geom_textcurve
#' @return a ggplot object.
sashimi_plot <- function(coverage,junction,filter_ratio = 20,
                         color_id = 1,region = NULL,lwr = 1,upr = 5,
                         color_set = NULL){
  rownames(coverage) = coverage$chr

  junction$y_start = coverage[junction$start,"coverage"]
  junction$y_end = coverage[junction$end,"coverage"]
  junction = na.omit(junction)
  junction$start = as.numeric(junction$start)
  junction$end = as.numeric(junction$end)

  filter = max(coverage$coverage)/filter_ratio
  junction_filter = junction[junction$count > filter &
                               junction$y_start > filter &
                               junction$y_end > filter,]

  if(is.null(region)){
    region = "chr"
  }
  if(is.null(color_set)){
    color_set = list(c("steelblue","LightSkyBlue"),c("coral","LightSalmon"))
  }

  sashimi = ggplot()+
    geom_segment(data = coverage,aes(x = chr,y = coverage,xend = chr,yend = 0),
                 color = color_set[[color_id]][1],alpha = 0.5)+
    geom_textcurve(data = junction_filter,aes(x = start,y = y_start,xend = end,yend = y_end,label = round(count),
                                              linewidth = range_scale(log(count,10),lwr = lwr,upr = upr)),
                   linecolour = color_set[[color_id]][2],lineend = "butt",curvature = -0.1,textcolour = "black",
                   alpha = 0.75)+
    theme_classic()+
    xlab(region)+ylab("reads coverage")+ylim(c(0,max(coverage$coverage)*1.3))+
    theme(text = element_text(size = 15))+
    theme(legend.position="none")

  return(sashimi)
}

#' @title sashimi
#' @description generate the sashimi plot for a gene from the Splice object
#' @param spliceOb The input splice object
#' @param gene The target gene
#' @param cells The target cell population
#' @param ... Omitted parameters for aesthetic  parameters in sashimi_plot
#' @return a ggplot object.
sashimi <- function(spliceOb,gene,cells,...){
  data = getIsoform(spliceOb,gene,cells)

  sashimi_data = sashimi_plot_data(data$isoform,data$count)
  out_plot = sashimi_plot(sashimi_data[[1]],sashimi_data[[2]],...)

  return(out_plot)
}
