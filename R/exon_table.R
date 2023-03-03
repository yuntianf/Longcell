#' @include spliceOb.R

#' @title assign canonical exons to self annotated sub-exons
#' @description  assign canonical exons to self annotated sub-exons
#' @param start,end The start and end site of the canonical exon
#' @param bed self annotated bed for a gene
#' @param exon_id A vector of string to name the self annotated exon, should have the same length
#' as the self-annotated exons
#' @param strand_col The colname to indicate the strand in the self-annotated gene bed
#' @param start_col,end_col The colnames of start and end sites in the self-annotated gene bed
#' @return A vector of sub-exons contained by the canonical exon.
exon_corres <- function(start,end,bed,exon_id = NULL,
                        strand_col = "strand",start_col = "start",end_col = "end"){
  if(is.null(exon_id)){
    if(unique(bed[,strand_col]) == "+"){
      exon_id = 1:nrow(bed)
    }
    else{
      exon_id = nrow(bed):1
    }
  }
  else{
    if(length(exon_id) != nrow(bed)){
      stop("The exon id annotation should have the same length of gene bed!")
    }
  }

  exons = exon_id[which(bed[,start_col] >= start & bed[,end_col] <= end)]
  return(exons)
}

#' @title assign canonical isoform to self annotated sub-exons
#' @description  assign canonical isoform (represented by multi canonical exons) to
#' self annotated sub-exons
#' @param starts,ends The start and end sites of the canonical exons within the isoform
#' @param bed the gene bed annotation
#' @param exon_id A vector of string to name the self annotated exon, should have the same length
#' as the self-annotated exons
#' @param strand_col The colname to indicate the strand in the self-annotated gene bed
#' @param start_col,end_col The colnames of start and end sites in the self-annotated gene bed
#' @return A string of sub-exons within the isoform, seperated by "|".
transcript_corres <- function(starts,ends,bed,exon_id = NULL,
                              strand_col = "strand",start_col = "start",end_col = "end"){
  if(length(starts) != length(ends)){
    stop("The start and end point of exons don't correspond!")
  }

  exons <- lapply(1:length(starts),function(i){
    exon = exon_corres(starts[i],ends[i],bed,
                       exon_id,strand_col,start_col,end_col)
    return(exon)
  })
  exons <- unlist(exons)
  names(exons) <- NULL

  if(unique(bed[,strand_col]) == "+"){
    exons = sort(exons)
  }
  else{
    exons = sort(exons,decreasing = TRUE)
  }
  exons = paste(exons,collapse = "|")
  return(exons)
}

#' @title assign gtf exon annotations to self annotated sub-exons
#' @description  assign gtf exon annotations for a gene to self annotated sub-exons
#' @param gtf a dataframe of gtf annotation
#' @param bed a dataframe of self annotated gene bed
#' @param gene the gene to be be matched
#' @param gtf_gene_col the colname in the gtf which indicates the gene id
#' @param bed_gene_col the colname in the bed which indicates the gene id
#' @param transcript_col the colname in the gtf which indicates the name of transcripts
#' @param exon_id A vector of string to name the self annotated exon, should have the same length
#' as the self-annotated exons
#' @param strand_col The colname to indicate the strand in the self-annotated gene bed
#' @param gtf_start_col,gtf_end_col The colnames for start and end sites of the canonical exons
#' in the gtf annotation
#' @param bed_start_col,bed_end_col The colnames for start and end sites of the canonical exons
#' in the bed annotation
#' @return A list of sub-exon string named by isoforms in the gtf.
gtf_bed_corres <- function(gtf,bed,gene,
                           gtf_gene_col = "gene_id",bed_gene_col = "gene",
                           transcript_col = "transname",
                           exon_id = NULL,strand_col = "strand",
                           gtf_start_col = "start",gtf_end_col = "end",
                           bed_start_col = "start",bed_end_col = "end"){
  gene_gtf = gtf[gtf[,gtf_gene_col] == gene,]
  gene_bed = bed[bed[,bed_gene_col] == gene,]

  transcript_exons <- gene_gtf %>%
    group_by(across(all_of(transcript_col))) %>%
    summarise(exons = transcript_corres(starts = !!sym(gtf_start_col),
                                        ends = !!sym(gtf_end_col),
                                        bed = gene_bed,
                                        exon_id = NULL,
                                        strand_col = strand_col,
                                        start_col = bed_start_col,
                                        end_col = bed_end_col))

  return(transcript_exons)
}

#' @title Build exon table from isoforms
#' @description  Count the exon incorporation for each isoform
#' @details Input isoforms and return a table of exon count for each isoform in cencoring way
#' @param exons A vector of string representing isoform, should be a sequence of exon id.
#' @param count A numerical vector representing the count for each isoform
#' @param polyA A numerical vector representing the polyA existence for each isoform
#' @param split A charater used to split the exon ids in the isoform
#' @importFrom stats na.omit
#' @return A data frame for the exon count in each isoform,
#' each row is an isoform, each col is an exon.
exonTable <- function(exons,count,polyA, split = "|"){
  exon_uniq = unique(exons)

  exon_vec_dic <- lapply(exon_uniq,function(x){
    exon_seq = as.numeric(unlist(strsplit(x,split = split,fixed = TRUE)))
    na_id = which(is.na(exon_seq))
    if(sum(is.na(exon_seq)) > 0){
      return(NULL)
    }
    exon_seq = na.omit(exon_seq)
    return(exon_seq)
  })
  names(exon_vec_dic) <- exon_uniq

  exon_dic = sort(unique(unlist(exon_vec_dic)))
  if(is.null(exon_dic)){
    return(NULL)
  }

  exon_vec <- lapply(1:length(exons),function(i){
    vec = exon_vec_dic[[exons[i]]]
    exon_count = count[[i]]

    out = rep(0,length(exon_dic))
    names(out) <- exon_dic
    if(length(vec) == 0){
      out[1:length(out)] <- NA
    }
    else{
      out[as.character(vec)] = exon_count
      out[exon_dic < min(vec)] = NA
      if(polyA[i] < 0.5){
        out[exon_dic > max(vec)] = NA
      }
    }
    return(out)
  })

  exon_vec <- as.data.frame(do.call(rbind,exon_vec))
  return(exon_vec)
}

#' @title Build exon table from annotated isoforms
#' @description  Count the exon incorporation for each annotated isoform
#' @details Input annotated isoforms and return a table of exon count for each isoform
#' @param gene_bed A dataframe of a bed annotation.
#' @param gtf A dataframe of gtf annotation.
#' @param gene The target gene to build annotated exon table
#' @param split A charater used to split the exon ids in the isoform
#' @return A data frame for the exon count in each annotated isoform,
#' each row is an isoform, each col is an exon.
annoExonTable <- function(gene_bed,gtf,gene,split = "|"){
  isoforms = gtf_bed_corres(gtf = gtf,bed = gene_bed,gene = gene)
  isoform_table = exonTable(exons = isoforms$exons,count = rep(1,nrow(isoforms)),
                            polyA = rep(TRUE,nrow(isoforms)),split = split)
  isoform_table[is.na(isoform_table)] = 0
  return(isoform_table)
}


#' @title judge if two exons can be merged
#' @description  judge if two exons can be merged
#' @details If two exons are co-existing or mutually exclusive, they have the same information,
#' and can considered as one exon in downstream analysis
#' @param a A numerical vector for the existence of an exon in each isoform
#' @param b A numerical vector for the existence of an exon in each isoform
#' @importFrom stats na.omit
#' @return TRUE or FALSE, indicating if the two exons should be merged
ifMerge <- function(a,b){
  if(length(a) != length(b)){
    stop("Length of two input vectors should be the same!")
  }

  exon_corres = sum(na.omit(xor(a,b)))
  if(exon_corres == 0 | exon_corres == length(na.omit(a))){
    na_corres = xor(is.na(a),is.na(b))
    if(sum(na_corres) == 0){
      return(TRUE)
    }
    else{
      return(NA)
    }
  }
  else{
    return(FALSE)
  }
}

#' @title judge if exons can be merged
#' @description  judge if each exon in the gene can be merged with other
#' @details Iteratively check if each pair of exons in the gene co-existing or
#' mutually exclusive
#' @param exon_table The exon count table output from exonTable()
#' @return A square matrix of logical values, indicating if each two exons should be merged
#' @import tidyr
exonTableMergeFlag <- function(exon_table){
  exons = colnames(exon_table)

  merge_table <- lapply(1:ncol(exon_table),function(i){
    merge_vec <- lapply(i:ncol(exon_table),function(j){
      return(c(i,j,ifMerge(exon_table[,i],exon_table[,j])))
    })
    merge_vec = do.call(rbind,merge_vec)
    return(merge_vec)
  })

  merge_table = as.data.frame(do.call(rbind,merge_table))
  colnames(merge_table) = c("exon1","exon2","if_merge")

  merge_matrix = tidyr::pivot_wider(merge_table,names_from = "exon2",values_from = "if_merge")
  merge_matrix = as.data.frame(merge_matrix)

  rownames(merge_matrix) = exons[merge_matrix$exon1]
  merge_matrix = merge_matrix[,-1]
  colnames(merge_matrix) = rownames(merge_matrix)

  merge_matrix[lower.tri(merge_matrix)] <- t(merge_matrix)[lower.tri(merge_matrix)]
  diag(merge_matrix) = 1
  return(merge_matrix)
}


#' @title combine the table for exon merge
#' @description  combine merge table from data(a) and annotation(b)
#' @details Exons are merged according to their co-existence or mutual exclusiveness in isoforms.
#' Isoforms in data can be different in annotations, this function will choose the intersection for
#' the merge flag
#' @param a exon merge flag table from data
#' @param b exon merge flag table from annotation
#' @return A square matrix of logical values, indicating if each two exons should be merged
dataAnnoCombine <- function(a,b){
  exons = colnames(a)
  if(length(setdiff(exons,colnames(b))) > 0){
    stop("The exons in the data and annotation don't correspond")
  }

  b = b[exons,exons]
  if(length(a) != length(b)){
    stop("The size of the matrix should be the same")
  }
  a = as.matrix(a)
  b = as.matrix(b)

  merge_flag = sapply(1:length(a),function(i){
    if(!is.na(a[i])){
      return(a[i])
    }
    else if(!is.na(b[i])){
      return(b[i])
    }
    else{
      stop("The value from annotation should be definite!")
    }
  })

  merge_flag = matrix(merge_flag,nrow(a),ncol(a))
  rownames(merge_flag) = rownames(a)
  colnames(merge_flag) = colnames(a)
  return(merge_flag)
}


#' @title collect exons to be merged into a list
#' @description  extract exons to be merged according to the exon merge flag table
#' @details exons which can be merged according to the exon merge flag table are extracted
#' together into a group. All different merge groups are collected as a list
#' @param data A square matrix of logical values, indicating if each two exons should be merged
#' @return A list of vectors with variable length,each vector represents a group of exons to be merged
moduleExtract <- function(data){
  if(nrow(data) != ncol(data)){
    stop("The similarity matrix should be symmetric!")
  }
  if(is.null(rownames(data)) | is.null(colnames(data))){
    warning("There is no annotation of each row and col of the data, will use the id")
    rownames(data) = 1:nrow(data)
    colnames(data) = 1:ncol(data)
  }

  data = data[order(rowSums(data),decreasing = TRUE),
              order(rowSums(data),decreasing = TRUE)]
  exons = rownames(data)

  merge_list = list()
  while(length(data) > 1){
    if(sum(data) == nrow(data)){
      merge_list = c(merge_list,as.list(rownames(data)))
      break
    }
    id = which(data[1,] == 1)
    subdata = data[id,id]
    while(sum(subdata == 0) > 0){
      filter = which(rowSums(subdata) == min(rowSums(subdata)))[1]
      id = id[-filter]
      subdata = subdata[-filter,-filter]
    }
    merge_list = c(merge_list,list(rownames(data)[id]))
    data = data[-id,-id]
  }
  merge_list <- c(merge_list,as.list(setdiff(exons,unlist(merge_list))))
  return(merge_list)
}

#' @title combine the merge flag from data and annotation
#' @description combine the merge flag from data and annotation
#' @details Input the exon merge flag matrix generated from data and annotation and output
#' exon groups to be merged
#' @param dataExonMergeMatrix exon merge flag matrix generated from data
#' @param annoExonMergeMatrix exon merge flag matrix generated from annotation
#' @return A list of vectors with variable length,each vector represents a group of exons to be merged
exonMatrixMerge <- function(dataExonMergeMatrix,annoExonMergeMatrix = NULL){
  if(is.null(annoExonMergeMatrix)){
    annoExonMergeMatrix = matrix(0,nrow(dataExonMergeMatrix),ncol(dataExonMergeMatrix))
    annoExonMergeMatrix = as.data.frame(annoExonMergeMatrix)
    rownames(annoExonMergeMatrix) = rownames(dataExonMergeMatrix)
    colnames(annoExonMergeMatrix) = colnames(dataExonMergeMatrix)
  }

  exons = rownames(dataExonMergeMatrix)
  annoExonMergeMatrix = annoExonMergeMatrix[exons,exons]
  if(sum(is.na(annoExonMergeMatrix)) > 0){
    stop("The exons in the data and annotation don't correspond!")
  }

  if(nrow(dataExonMergeMatrix) != nrow(annoExonMergeMatrix) |
     ncol(dataExonMergeMatrix) != ncol(annoExonMergeMatrix)){
    stop("The number of exons in data and annotations should be the same!")
  }

  exon_merge_matrix = dataAnnoCombine(dataExonMergeMatrix,annoExonMergeMatrix)

  merge_list = moduleExtract(exon_merge_matrix)
  return(merge_list)
}

#' @title merge two exons
#' @description merge the existence of two exons in each isoform
#' @details Input the count vectors for two exons in each isoform, cencoring in one exon can
#' be recovered by the information from the other
#' @param a count vector for exon a
#' @param b count vector for exon b
#' @return a count vector for merged exon
exonExonMerge <- function(a,b){
  if(length(a) != length(b)){
    stop("The size of the matrix should be the same")
  }

  if(sum(xor(a,b),na.rm = TRUE) > 0){
    b = !b
  }

  merge_flag = sapply(1:length(a),function(i){
    if(is.na(a[i])){
      return(b[i])
    }
    if(is.na(b[i])){
      return(a[i])
    }
    if(!is.na(a[i]) & !is.na(b[i])){
      if(a[i] == b[i]){
        return(a[i])
      }
      else{
        stop("There exist collision between exons to be merged!")
      }
    }
  })

  return(merge_flag)
}

#' @title merge multi-exons
#' @description merge the existence of multi-exons in each isoform
#' @details Input a dataframe of count for mutil-exons in each isoform, mering them by pair iteratively
#' @param data a dataframe of count for mutil-exons in each isoform
#' @return a count vector for merged exon
multiExonsMerge <- function(data){
  if(is.null(ncol(data)) || ncol(data) == 1){
    return(data)
  }

  merged_exon = rep(NA,nrow(data))
  for(i in 1:ncol(data)){
    merged_exon = exonExonMerge(merged_exon,data[,i])
  }
  return(merged_exon)
}

#' @title merge exons in an exon table according to to be merged list
#' @description merge exons in an exon table according to to be merged list
#' @details Input a dataframe of count for mutil-exons in each isoform and a list
#' indicating which exons should be merged together. Exons in the same group will be
#' merged into a meta-exon
#' @param exon_table a dataframe of count for mutil-exons in each isoform
#' @param merge_list A list of vectors with variable length,each vector represents a group of exons to be merged
#' @param sep The seperation character between exons in the isoform string
#' @return a count vector for merged exon
exonTableMergeList <- function(exon_table,merge_list,sep = "|"){
  gene_count = apply(exon_table,1,function(x){
    return(max(x,na.rm = TRUE))
  })
  exon_table = exon_table %>% mutate_all(as.logical)
  merged_exon_table <- lapply(merge_list,function(x){
    temp = multiExonsMerge(exon_table[,x])
  })

  merged_exon_table = as.data.frame(do.call(cbind,merged_exon_table))

  merged_exons = sapply(merge_list,function(x){
    exons = paste(x,collapse = sep)
  })

  colnames(merged_exon_table) = merged_exons
  merged_exon_table$gene_count = gene_count
  return(merged_exon_table)
}

#' @title merge exons
#' @description merge exons according to their relationship, co-existing or mutual exclusive
#' @details Input a dataframe of count for mutil-exons in each isoform and the bed
#' and gtf annotation to identify the relationship between exon pairs. Exon sets which are
#' agreed to be merged by both data and annotation are merged
#' @param exon_table a dataframe of count for mutil-exons in each isoform
#' @param gene_bed gene bed annotation
#' @param gtf gtf annotation
#' @param gene the gene to refer annotation
#' @param sep The seperation character between exons in the isoform string
#' @return a count matrix for merged exon
exonTableMerge <- function(exon_table,
                           gene_bed = NULL,gtf = NULL,
                           gene = NULL,sep = "|"){
  if(is.null(ncol(exon_table)) || ncol(exon_table) == 1){
    exon_table = as.data.frame(exon_table)
    gene_count = apply(exon_table,1,max)
    gene_count[is.na(gene_count)] = 0
    merged_table = as.data.frame(cbind(exon_table,gene_count))
    colnames(merged_table) = c("exon","gene_count")
    return(merged_table)
  }
  exon_merge_flag = exonTableMergeFlag(exon_table)

  if(!is.null(gene_bed) & !is.null(gtf) & !is.null(gene)){
    anno_exon_table = annoExonTable(gene_bed = gene_bed,gtf = gtf,
                                    gene = gene,split = sep)
    anno_merge_flag = exonTableMergeFlag(anno_exon_table)
  }
  else{
    anno_merge_flag = NULL
  }
  merge_list = exonMatrixMerge(exon_merge_flag,anno_merge_flag)
  merged_table = exonTableMergeList(exon_table,merge_list)

  return(merged_table)
}

#' @title filter out too short exons
#' @description filter out too short exons which can be confused with insertions or deletions
#' @details Some sub-exons are too short and their existence can be confused by insertions
#' or deletions, which can misled the downstream alternative splicing analysis.
#' @param gene name of the target gene
#' @param gene_bed gene bed annotation
#' @param thresh the minimum length for the exon to be
#' @param exons selected exons, exon filtering will be within this selection
#' @param exon_id the vector of names for each exon
#' @param gene_id the name of gene column in the bed file
#' @param len the name of exon length column in the bed file
#' @param strand the name of strand column in the bed file
#' @param sep The seperation character between exons in the isoform string
#' @return a vector of preserved exons
exon_len_filter <- function(gene,gene_bed,thresh = 10,exons = NULL,
                            exon_id = NULL,gene_id = "gene",
                            len = "width",strand = "strand",sep = "|"){
  gene_bed = gene_bed[gene_bed[,gene_id] == gene,]

  if(is.null(exon_id)){
    if(unique(gene_bed[,strand]) == "+"){
      exon_id = 1:nrow(gene_bed)
    }
    else{
      exon_id = nrow(gene_bed):1
    }
  }

  if(is.null(exons)){
    exons = exon_id
  }

  filtered = sapply(exons,function(x){
    e = unlist(strsplit(x,split = sep,fixed = TRUE))
    exons_len = sum(gene_bed[exon_id %in% e,len])
    if(exons_len >= thresh){
      return(TRUE)
    }
    else{
      return(FALSE)
    }
  })
  return(exons[filtered])
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
extractExonTable <- function(spliceOb,gene,gene_bed = NULL,gtf = NULL,
                                  cells = "all",exons = "all",
                                  exon_len_thresh = 10,sep = "|"){
  isoform = getIsoform(spliceOb,genes = gene,cells = cells)

  if(length(isoform) == 0 || nrow(isoform) == 0){
    return(NULL)
  }

  isoform_table = suppressWarnings(exonTable(isoform$exons,isoform$count,isoform$polyA))
  isoform_table = cbind(1:nrow(isoform),isoform_table)
  colnames(isoform_table)[1] = "cell_id"

  row_filter = rowSums(isoform_table,na.rm = TRUE)-isoform_table[,1]
  col_filter = colSums(isoform_table,na.rm = TRUE)
  isoform_table = isoform_table[row_filter > 0,col_filter > 0]

  if(length(isoform_table) == 0 || nrow(isoform_table) == 0){
    return(NULL)
  }
  all_exons = setdiff(colnames(isoform_table),"cell_id")

  isoform_table = cbind(isoform_table[,"cell_id"],
                  suppressWarnings(exonTableMerge(isoform_table[,all_exons],
                                       gene_bed = gene_bed,gtf = gtf,
                                       gene = gene,sep = sep)))
  colnames(isoform_table)[1] = "cell_id"
  if(length(all_exons) == 1){
    colnames(isoform_table)[2] = all_exons
  }

  if(exons[1] == "all"){
    exons = setdiff(colnames(isoform_table),c("cell_id","gene_count"))
    if(!is.null(gene_bed)){
      exons = exon_len_filter(gene = gene,gene_bed = gene_bed,
                              thresh = exon_len_thresh,exons = exons)
    }
  }
  else{
    if(length(setdiff(exons,colnames(isoform_table))) > 0){
      warning("Some designated exons don't exsit or be merged, will use the intersection!")
    }
    exons_id = sapply(exons,function(i){
      exon_search(exon = i,exons_group = colnames(isoform_table))
    })
    exons = colnames(isoform_table)[unique(unlist(exons_id))]
  }

  isoform_table = isoform_table[,colnames(isoform_table) %in% c("cell_id",exons,"gene_count")]
  isoform_table$cell_id = isoform[isoform_table$cell_id,"cell"]
  return(isoform_table)
}

#' @title count exons in each cell
#' @description summarise the exon count and gene count for each single cell
#' @details Due to the truncation, the exons are counted in a censoring way, only
#' exons between two exons or between an exon and polyA will be treated as spliced out,
#' otherwise it is censor and won't be counted.
#' @param exon_table a dataframe of count for mutil-exons in each isoform
#' @param cell_col the name of the cell column in the exon table
#' @return a list of exon count table and gene count table for each cell
#' @import dplyr
#' @importFrom stats na.omit
cellGeneExonCount <- function(exon_table,cell_col = "cell_id"){
  exon_table <- exon_table %>% mutate_at(vars(-all_of(cell_col)),as.numeric)

  gene_table = apply(exon_table,1,function(x){
    x[2:(ncol(exon_table)-1)][!is.na(x[2:(ncol(exon_table)-1)])] = x[ncol(exon_table)]
    return(x[1:(ncol(exon_table)-1)])
  })
  gene_table = as.data.frame(t(gene_table))
  gene_table <- gene_table %>% mutate_at(vars(-all_of(cell_col)),as.numeric)

  exons_table = apply(exon_table,1,function(x){
    x[2:(ncol(exon_table) - 1)] = as.numeric(x[2:(ncol(exon_table) - 1)])*as.numeric(x[ncol(exon_table)])
    return(x[1:(ncol(exon_table)-1)])
  })
  exons_table = as.data.frame(t(exons_table))
  exons_table <- exons_table %>% mutate_at(vars(-all_of(cell_col)),as.numeric)

  cell_gene_count = gene_table %>%
    group_by(across(all_of(cell_col))) %>%
    summarise_all(list(~sum(na.omit(.))))
  cell_exon_count = exons_table %>%
    group_by(across(all_of(cell_col))) %>%
    summarise_all(list(~sum(na.omit(.))))

  return(list(cell_gene_count,cell_exon_count))
}

#' @title count exons in each cell from a Splice object
#' @description summarise the exon count and gene count for each single cell from a Splice object
#' @details Due to the truncation, the exons are counted in a censoring way, only
#' exons between two exons or between an exon and polyA will be treated as spliced out,
#' otherwise it is censor and won't be counted.
#' @param spliceOb the Splice object
#' @param gene name of the target gene
#' @param gene_bed gene bed annotation
#' @param gtf gtf annotation
#' @param cells cells to use in building exon count table
#' @param exons exon to use in building exon count table
#' @param exon_len_thresh the minimum length of an exon to be preserved in exon count table
#' @param sep The seperation character between exons in the isoform string
#' @return a list of exon count table and gene count table for each cell
#' @export
gene_exons_table.base <- function(spliceOb,gene,gene_bed = NULL,gtf = NULL,
                                  cells = "all",exons = "all",
                                  exon_len_thresh = 10,sep = "|"){
  exon_table = extractExonTable(spliceOb,gene,gene_bed = gene_bed,gtf = gtf,
                                cells = cells,exons = exons,
                                exon_len_thresh = exon_len_thresh,sep = sep)
  if(is.null(exon_table) || nrow(exon_table) == 0){
    return(NULL)
  }

  count_list = cellGeneExonCount(exon_table)
  return(count_list)
}

#' @title generic gene_exons_table function definition
#' @param object the Splice or Seurat object
#' @inheritParams gene_exons_table.base
#' @param ... other possible parameters for genes_groups_GLRT.base
#' @export
setGeneric("gene_exons_table",
           function(object,gene,...) standardGeneric("gene_exons_table"))

#' @title generic gene_exons_table function for Splice object
#' @param object a Splice object
#' @inheritParams gene_exons_table.base
#' @param ... parameters for gene_exons_table.base
#' @export
setMethod("gene_exons_table",
          signature(object = "Splice",gene = "character"),
          function(object,gene,gene_bed = NULL,gtf = NULL,...){
            gene_exons_table.base(object,gene = gene,
                                  gene_bed = gene_bed,gtf = gtf,...)
          }
)

#' @title generic gene_exons_table function for Seurat object
#' @param object a Seurat object with a Splice object in misc
#' @param slot the slot name for the Splice object embedded in the Seurat object
#' @inheritParams gene_exons_table.base
#' @param ... parameters for gene_exons_table.base
#' @import Seurat
#' @export
setMethod("gene_exons_table",
          signature(object = "Seurat",gene = "character"),
          function(object,gene,slot = "splice",
                   gene_bed = NULL,gtf = NULL,...){
            spliceOb = getSplice(object,slot = slot)
            gene_exons_table.base(spliceOb,gene = gene,
                                  gene_bed = gene_bed,gtf = gtf,...)
          }
)
