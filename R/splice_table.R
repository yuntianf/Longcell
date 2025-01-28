#' @include spliceOb.R


#' @title isoform2sites
#' @description  split the isoform string to a vector of splicing sites
#' @details split the isoform string to a vector of splicing sites
#' @param isoform A string representing the isoform
#' @param split The chracter to split each exon within the isoform
#' @param sep The chracter to split each the start and end sites of an exon
#' @return A string vector, each element is a splicing site
isoform2sites = function(isoform,split = "|",sep = ","){
  split = paste(c("\\",split,"|","\\",sep),collapse = "")
  sites = unlist(strsplit(isoform,split = split))

  return(sites)
}

#' @title sites2vec
#' @description  convert the splicing sites vector to a binary vector
#' @details convert the splicing sites vector to a binary vector, the spliced-in sites will be set as 1,
#' the spliced-out sites will be set as 0, the depleted sites which may be due to technical truncations
#' are set as NA.
#' @param start A number recording the start position for a read
#' @param in_site A string vector recording the start position for each exon for a read
#' @param out_site A string vector recording the end position for each exon for a read
#' @param end A number recording the end position for a read
#' @param polyA A bool variable to indicate if the read has a polyA
#' @param end A number recording the end position for a read
#' @param sites_index A string vector recording all sites across all reads in this gene
#' @param strand The strand of this gene
#' @return A binary vector indicating the splicing sites existence for a read
sites2vec = function(start,in_site,out_site,end,polyA,sites_index,strand){
  if(length(in_site) != length(out_site)){
    stop("The length of in and out splicing sites should match!")
  }

  if(!is.na(in_site[1]) & !is.na(out_site[1])){
    sites_vec = as.data.frame(rbind(cbind(in_site,"in"),cbind(out_site,"out")))
    colnames(sites_vec) = c("sites","state")
    sites_vec = sites_vec %>% mutate(sites = as.numeric(sites),status = 1)

    sites_index = left_join(sites_index,sites_vec,by = c("sites","state"))
    sites_index$status[is.na(sites_index$status)] = 0
  }
  else{
    sites_index$status = 0
  }


  if(strand == "+"){
    sites_index[sites_index$sites <= start,"status"] = NA
    if(!polyA){
      sites_index[sites_index$sites >= end,"status"] = NA
    }
  }else if(strand == "-"){
    sites_index[sites_index$sites >= end,"status"] = NA
    if(!polyA){
      sites_index[sites_index$sites <= start,"status"] = NA
    }
  }
  return(sites_index$status)
}


#' @title sites2matrix
#' @description  convert the splicing sites vector for all reads to a binary matrix
#' @details This function applies sites2vec to all reads and concat all vector results to a matrix
#' @param data A dataframe recording the splicing sites for a read, including start, middle splicing sites, end
#' and polyA
#' @param in_sites_index A string vector recording unique start position for all middle exons
#' @param out_sites_index A string vector recording unique end position for all middle exons
#' @param strand The strand of this gene
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom magrittr %>%
#' @return A binary matrix indicating the splicing sites existence for all reads from a gene
sites2matrix = function(data,in_sites_index,out_sites_index,strand){
  sites_index = as.data.frame(rbind(cbind(in_sites_index,"in"),cbind(out_sites_index,"out")))
  colnames(sites_index) = c("sites","state")
  sites_index = sites_index %>% mutate(sites = as.numeric(sites)) %>% arrange(sites,state)

  out = lapply(1:nrow(data),function(i){
    x = data[i,]
    vec = sites2vec(x$start,unlist(x$in_site),unlist(x$out_site),
                    x$end,x$polyA,sites_index,strand)
    return(vec)
  })

  sites = paste(sites_index$sites,sites_index$state,sep = "_")
  out = as.data.frame(do.call(rbind,out))
  colnames(out) = sites
  return(out)
}

#' @title sites_coef
#' @description  judge if two sites can be merged
#' @details If two exons are co-existing or mutually exclusive, they have the same information,
#' and can considered as one exon in downstream analysis
#' @param a A numerical vector for the existence of an exon in each isoform
#' @param b A numerical vector for the existence of an exon in each isoform
#' @param count The read count
#' @importFrom stats na.omit
#' @return TRUE or FALSE, indicating if the two exons should be merged
sites_coef <- function(a,b,count){
  if(length(a) != length(b) | length(a) != length(count)){
    stop("Length of two input vectors and the count should be the same!")
  }

  non_na = sum(count[!is.na(a) & !is.na(b)])
  one_na = sum(count[!is.na(a) & is.na(b)]) + sum(count[is.na(a) & !is.na(b)])

  flag = !is.na(a) & !is.na(b)
  sub_a = a[flag]
  sub_b = b[flag]
  sub_count = count[flag]

  coexist = sum(sub_count[sub_a == sub_b])
  exclusion = sum(sub_count[sub_a != sub_b])
  coef = (coexist-exclusion)/(coexist+exclusion)

  if(non_na < one_na){
    return(0.01*coef/abs(coef))
  }

  return(coef)
}

#' @title sites_coef_matrix
#' @description  judge if each site in the gene can be merged with other
#' @details Iteratively check if each pair of exons in the gene co-existing or
#' mutually exclusive
#' @param sites_table The site count table output from sites2matrix()
#' @return A square matrix of logical values, indicating if each two exons should be merged
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr select
#' @importFrom magrittr %>%
sites_coef_matrix <- function(sites_table,count){
  if(nrow(sites_table) != length(count)){
    stop("The size of reads and their count should match!")
  }
  if(sum(count) < 10){
    warning("Too few reads for this gene, the estimation of coefficiency between each splicing sites may not be accurate!")
  }
  sites = colnames(sites_table)

  coef_matrix = lapply(1:length(sites),function(i){
    temp = lapply(i:length(sites),function(j){
      coef = sites_coef(sites_table[,i],sites_table[,j],count)
      return(c(i,j,coef))
    })
    temp = as.data.frame(do.call(rbind,temp))
    return(temp)
  })
  coef_matrix = as.data.frame(do.call(rbind,coef_matrix))
  colnames(coef_matrix) = c("site1","site2","coef")
  # coef_matrix$coef = abs(coef_matrix$coef)

  coef_matrix = tidyr::pivot_wider(coef_matrix,names_from = "site2",values_from = "coef")
  coef_matrix = as.data.frame(coef_matrix)

  rownames(coef_matrix) = sites[coef_matrix$site1]
  coef_matrix = coef_matrix %>% dplyr::select(-site1)
  colnames(coef_matrix) = rownames(coef_matrix)

  coef_matrix[lower.tri(coef_matrix)] <- t(coef_matrix)[lower.tri(coef_matrix)]
  diag(coef_matrix) = 1

  coef_matrix[is.na(coef_matrix)] = 0
  return(coef_matrix)
}

#' @title moduleExtract
#' @description  extract exons to be merged according to the exon merge flag table
#' @details sites which can be merged according to the exon merge flag table are extracted
#' together into a group. All different merge groups are collected as a list
#' @param coef_matrix A square matrix of logical values, indicating if each two sites should be merged
#' @param eps The threshold for two sites to be merged
#' @importFrom dbscan dbscan
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @return A list of vectors with variable length,each vector represents a group of sites to be merged
moduleExtract <- function(coef_matrix,eps = 0.05){
  if(nrow(coef_matrix) != ncol(coef_matrix)){
    stop("The similarity matrix should be symmetric!")
  }
  if(is.null(rownames(coef_matrix)) | is.null(colnames(coef_matrix))){
    warning("There is no annotation of each row and col of the data, will use the id")
    rownames(coef_matrix) = 1:nrow(coef_matrix)
    colnames(coef_matrix) = 1:ncol(coef_matrix)
  }
  dis_matrix = 1- abs(coef_matrix)
  cluster = dbscan::dbscan(as.dist(dis_matrix),eps = eps,minPts = 1)$cluster
  cluster = as.data.frame(cbind(rownames(coef_matrix),cluster))
  colnames(cluster) = c("node","cluster")

  cluster = cluster %>% group_by(cluster) %>% summarise(node = list(node))
  return(cluster$node)
}

#' @title siteMerge
#' @description  Merge sites within the same group
#' @details Merge sites within the same group
#' @param site_matrix A binary matrix to indicate if the sites exist
#' @return A vector for the meta splicing site.
siteMerge = function(site_matrix){
  if(is.null(nrow(site_matrix))){
    return(site_matrix)
  }
  merge = round(rowMeans(site_matrix,na.rm = TRUE))
  return(merge)
}

#' @title metaSpliceSite
#' @description  Merge sites into meta sites based on their co-existence or if they are mutually exclusive
#' @details Merge sites into meta sites based on their co-existence or if they are mutually exclusive
#' @inheritParams sites_coef_matrix
#' @inheritParams moduleExtract
#' @return A list with two elements, the first one is the count matrix of meta splicing sites,
#' the second one is the component of the meta splicing sites
metaSpliceSite = function(sites_table,count,eps = 0.05){
  if(is.na(eps)){
    meta = sites_table
    cluster = as.list(colnames(sites_table))
  }
  else{
  coef = sites_coef_matrix(sites_table,count)
  cluster = moduleExtract(coef,eps = eps)
  meta = lapply(cluster,function(x){
    sub_sites_table = sites_table[,x]
    if(length(x) > 1){
      sub_coef = as.matrix(coef[x,x])
      relation = sub_coef %*% sub_coef

      inverse_id = relation[1,] < 0
      if(sum(inverse_id) > 0){
        sub_sites_table[,inverse_id] = 1-sub_sites_table[,inverse_id]
      }
    }
    return(siteMerge(sub_sites_table))
  })
  meta = as.data.frame(do.call(cbind,meta))
  colnames(meta) = paste("metasite",1:length(cluster),sep = "_")
  }
  out = list(meta,cluster)
  names(out) = c("meta_splice_site_table","meta_splice_site")
  return(out)
}

#' @title spliceTable
#' @description  Merge sites into meta sites for a gene
#' @details Merge sites into meta sites for a gene
#' @param data A dataframe recording the isoform for a gene
#' @param strand The strand of the target gene
#' @param eps The threshold for two sites to be merged
#' @param cell_col The name of the column recording the cell barcode
#' @param iso_col The name of the column recording the isoform string
#' @param count_col The name of the column recording the isoform count in each cell
#' @param polyA_col The name of the column recording the polyA existence for a read
#' @inheritParams isoform2sites
#' @importFrom dplyr mutate_at
#' @importFrom dplyr group_by_at
#' @importFrom dplyr summarise
#' @importFrom dplyr filter
#' @return A metaSite object for a gene
spliceTable = function(data,strand,eps = 0.05, cell_col = "cell", iso_col = "isoform",
                       count_col = "count",polyA_col = "polyA", split = "|",sep = ","){
  iso_uniq = unique(data[,iso_col])

  iso_sites_index = lapply(iso_uniq,function(x){
    sites = isoform2sites(x)
    len = length(sites)
    start = sites[1]
    end = sites[len]
    if(len >= 4){
      in_sites = sites[seq(3,(len-1),2)]
      out_sites = sites[seq(2,(len-2),2)]
    }
    else{
      in_sites = NA
      out_sites = NA
    }
    return(list(x,start,in_sites,out_sites,end))
  })
  iso_sites_index = as.data.frame(do.call(rbind,iso_sites_index))
  if(length(iso_sites_index) == 0 || nrow(iso_sites_index) == 0){
    return(NULL)
  }
  colnames(iso_sites_index) = c("isoform","start","in_site","out_site","end")
  iso_sites_index = iso_sites_index %>% mutate_at(c("isoform","start","end"),unlist)

  data = inner_join(data, iso_sites_index, by = "isoform")
  data[,"polyA"] = ifelse(data[,"polyA"] > 0.5,TRUE,FALSE)
  data = data %>% group_by_at(c("cell","isoform","polyA")) %>%
    summarise(start = unique(start),
              in_site = unique(in_site),
              out_site = unique(out_site),
              end = unique(end),
              count = sum(count),.groups = "drop")
  data = data %>% filter(count > 0)

  in_sites_index = sort(as.numeric(unique(unlist(data$in_site))))
  out_sites_index = sort(as.numeric(unique(unlist(data$out_site))))

  if(length(in_sites_index) == 0 & length(out_sites_index) == 0){
    return(NULL)
  }

  site_vec = sites2matrix(data,in_sites_index,out_sites_index,strand)

  #return(list(data,site_vec))

  metasite = metaSpliceSite(site_vec,unlist(data[,count_col]),eps)

  metasite[[1]] = as.data.frame(cbind(data[,c(cell_col,count_col)],metasite[[1]]))

  return(metasite)
}

#' @title Search if an exon is within an exon group.
#' @description Search if an exon is within an exon group.
#' @details Some exons can be merged into a meta-exon due to co-existence or mutually exclusiveness,
#' when such exons are required to be observed, we need to search them in the meta-exon.
#' @param exon The specified exon to be observed.
#' @param exons_group A vector of strings, each one is a meta-exon
#' @param sep The separator between each exon in the meta-exon
#' @return A numerical value representing the id of the meta-exon which incorporated the specified
#' exon.
exon_search <- function(exon,exons_group,sep = "|"){
  exons_group = strsplit(exons_group,split = sep,fixed = TRUE)
  id = sapply(exons_group,function(x){
    flag = exon %in% x
    return(flag)
  })
  if(sum(id) == 0){
    stop("This exon doesn't exist!")
  }
  else if(sum(id) > 1){
    warning("The exon is duplicated!")
  }
  return(which(id)[1])
}

#' @title build an exon count table for a gene
#' @description build an exon count table for a gene, each row is an isoform in a cell
#' and each column is a meta-exon
#' @details build an exon count table for a gene, each row is an isoform in a cell
#' and each column is a meta-exon
#' @param spliceOb the input splice object
#' @param gene name of the target gene
#' @param gene_bed gene bed annotation
#' @param gtf gtf annotation
#' @param cells cells to use in building exon count table
#' @param exons exon to use in building exon count table
#' @param exon_len_thresh the minimum length of an exon to be preserved in exon count table
#' @param sep The seperation character between exons in the isoform string
#' @return an exon count table for isoforms in each cell
extractExonTable <- function(spliceOb,gene,gene_bed = NULL,
                                  cells = "all",exons = "all",
                                  exon_len_thresh = 10,sep = "|"){
  isoform = getIsoform(spliceOb,genes = gene,cells = cells)
  isoform = isoform %>% filter(count > 0)

  if(length(isoform) == 0 || nrow(isoform) == 0){
    return(NULL)
  }

  cell_id = isoform$cell
  isoform_table = exonTable(isoform$exon,isoform$count,isoform$polyA)

  row_filter = rowSums(isoform_table,na.rm = TRUE) > 0
  col_filter = colSums(isoform_table,na.rm = TRUE) > 0
  isoform_table = isoform_table[row_filter,col_filter]
  cell_id = cell_id[row_filter]

  if(length(isoform_table) == 0 || nrow(isoform_table) == 0){
    return(NULL)
  }

  isoform_table = exonTableMerge(isoform_table,sep = sep)
  gene_count = isoform_table[[2]]
  isoform_table = isoform_table[[1]]


  if(exons[1] == "all"){
    if(!is.null(gene_bed)){
      exons = exon_len_filter(gene = gene,gene_bed = gene_bed,
                              thresh = exon_len_thresh,exons = colnames(isoform_table))
    }
    else{
      exons = colnames(isoform_table)
    }
  }
  else{
    exons_id = sapply(exons,function(i){
      exon_search(exon = i,exons_group = colnames(isoform_table))
    })
    exons = colnames(isoform_table)[unique(unlist(exons_id))]
  }

  isoform_table = isoform_table %>% dplyr::select_at(exons)
  isoform_table = as.data.frame(cbind(cell_id,isoform_table,gene_count))
  return(isoform_table)
}

#' @title cellGeneSiteCount
#' @description summarise the site count and gene count for each single cell
#' @details Due to the truncation, the sites are counted in a censoring way, only
#' sites between two sites or between a site and the polyA will be treated as spliced out,
#' otherwise it is censor and won't be counted.
#' @param splice_site_table a dataframe of count for splicing sites in each isoform
#' @param cell_col The name of the column recording the cell barcode
#' @param count_col The name of the column recording the isoform count in each cell
#' @return a list of site count table and gene count table for each cell
#' @importFrom dplyr group_by_at
#' @importFrom dplyr summarise_at
#' @importFrom dplyr select
#' @importFrom dplyr one_of
#' @importFrom stats na.omit
#' @return A list including two dataframes, the fisrt one records the gene count for each site in each cell,
#' while the second one record the spliced-in site count.
cellGeneSiteCount <- function(splice_site_table,cell_col = "cell",count_col = "count"){
  cell_gene_count = splice_site_table %>% group_by_at(cell_col) %>%
    summarise_at(setdiff(colnames(splice_site_table),c(cell_col,count_col)),
                 ~sum(get(count_col)*!is.na(.)))

  cell_site_count = splice_site_table %>% group_by_at(cell_col) %>%
    summarise_at(setdiff(colnames(splice_site_table),c(cell_col,count_col)),
                 ~sum(get(count_col)*.,na.rm = TRUE))

  cell_gene_count = as.data.frame(cell_gene_count)
  rownames(cell_gene_count) = cell_gene_count[,cell_col]
  cell_gene_count = cell_gene_count %>% dplyr::select(-one_of(cell_col))

  cell_site_count = as.data.frame(cell_site_count)
  rownames(cell_site_count) = cell_site_count[,cell_col]
  cell_site_count = cell_site_count %>% dplyr::select(-one_of(cell_col))
  return(list(cell_gene_count,cell_site_count))
}


#' @title geneSiteTable_df
#' @description merge sites into meta sites and summarise the site count and gene count for each single cell
#' @details An integration function for spliceTable() and cellGeneSiteCount()
#' @inheritParams spliceTable
#' @inheritParams cellGeneSiteCount
#' @return A list including two dataframes, the fisrt one records the gene count for each site in each cell,
#' while the second one record the spliced-in site count.
#' @export
geneSiteTable_df = function(data,strand,eps = 0.05, cell_col = "cell", iso_col = "isoform",
                               count_col = "count",polyA_col = "polyA", split = "|",sep = ","){
  splice_table = spliceTable(data = data,strand = strand, eps = eps,cell_col = cell_col, iso_col = iso_col,
                             count_col = count_col,polyA_col = polyA_col, split = split,sep = sep)
  if(is.null(splice_table)){
    return(NULL)
  }
  cell_gene_site_count = cellGeneSiteCount(splice_table[[1]],cell_col = cell_col,count_col = count_col)
  cell_gene_site_count[[3]] = splice_table[[2]]
  return(cell_gene_site_count)
}

#' @title geneSiteTable.base
#' @description An interface function to apply geneSiteTable_df for target genes in a Splice object
#' @details An interface function to apply geneSiteTable_df for target genes in a Splice object
#' @param spliceOb A Splice object
#' @param gene_bed The gene bed annotation, should include the gene name and the gene strand.
#' @param genes The target genes to generate meta sites, default for all genes
#' @param eps The threshold for two sites to be merged
#' @param bed_gene_col The name of the column to record the gene name in the gene bed
#' @param bed_strand_col The name of the column to record the gene strand in the gene bed
#' @param overwrite The flag to decide if overwrite the computed metaSite object
#' @param verbose The flag to indicate if print the information
#' @return A Splice object with meta sites stored the in the meta_sites slot.
#' @export
geneSiteTable.base = function(spliceOb,gene_bed,genes= "all",eps = 0.05,
                              bed_gene_col = "gene",bed_strand_col = "strand",
                              overwrite = FALSE,
                              verbose = TRUE){
  if(genes[1] == "all"){
    genes = spliceOb@genes$gene
  }
  diff = setdiff(genes,spliceOb@genes$gene)
  if(length(diff) > 0){
    warning(paste(head(diff,10),collapse = ","),"... don't exsit in the splice object, will be ignored!")
    genes = genes[genes %in% spliceOb@genes$gene]
  }

  sub_bed = unique(gene_bed %>%
                   filter_at(bed_gene_col, ~.%in%genes) %>%
                   dplyr::select_at(c(bed_gene_col,bed_strand_col)))
  sub_bed = as.data.frame(sub_bed)

  data = getIsoform(spliceOb,genes = sub_bed[,bed_gene_col])

  temp = future_lapply(1:nrow(sub_bed),function(i){
    gene_i = sub_bed[i,bed_gene_col]
    if(verbose){
      print(gene_i)
    }

    exist = TRUE
    tryCatch({
      cache = getMetaSites(spliceOb,gene_i)
    },error = function(e) {
      exist <<- FALSE
    })

    if(exist & !overwrite){
      return(cache)
    }

    strand_i = sub_bed[i,bed_strand_col]
    sub_data = as.data.frame(data %>% filter(gene == gene_i))
    gene_i_site_table = geneSiteTable_df(sub_data,strand_i,eps = eps)
    if(is.null(gene_i_site_table)){
      return(NULL)
    }
    #return(gene_i_site_table)
    meta_site_ob = new("metaSite",
                       sites = gene_i_site_table[[3]],
                       cellGeneCount = gene_i_site_table[[1]],
                       cellSiteCount = gene_i_site_table[[2]])
    return(meta_site_ob)
  },future.seed=TRUE, future.globals=TRUE,future.packages = c("Longcell"))

  names(temp) = sub_bed[,bed_gene_col]
  temp<-temp[!sapply(temp,is.null)]
  spliceOb@meta_sites = temp
  return(spliceOb)
}


#' @title generic geneSiteTable function definition
#' @param object the Splice object
#' @inheritParams geneSiteTable.base
#' @param ... other possible parameters for genes_groups_GLRT.base
#' @export
setGeneric("geneSiteTable",
           function(object,gene_bed,genes,eps,...) standardGeneric("geneSiteTable"))

#' @title generic gene_exons_table function for Splice object
#' @param object a Splice object
#' @inheritParams geneSiteTable.base
#' @param ... parameters for geneSiteTable.base
#' @export
setMethod("geneSiteTable",
          signature(object = "Splice",gene_bed = "list",genes = "character"),
          function(object,gene_bed,genes = "all", ...){
            geneSiteTable.base(spliceOb = object,gene_bed = gene_bed,
                               genes= genes,...)
          }
)
