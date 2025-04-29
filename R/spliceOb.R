#' @import methods

#' @title The metaSite Class
#' @description The metaSite class is an object to store the meta splicing sites information. It includes three slots, the gene
#' count and spliced-in count per cell for each meta site the orginal splicing sites that each meta site corresponds to.
#' @slot cellGeneCount A cell by meta site matrix, recording the gene expression count.
#' @slot cellSiteCount A cell by meta site matrix, recording the spliced-in count.
#' @slot sites A list with length as the number of meta splicing sites, Each element records the corresponding original splicing sites.
#' @name metaSite-class
#' @rdname metaSite-class
#' @concept objects
#' @exportClass metaSite
#'
setClass("metaSite",
         representation(sites = "list",
                        cellGeneCount = "data.frame",
                        cellSiteCount = "data.frame"),
         prototype(sites = NULL,cellGeneCount = NULL,cellSiteCount = NULL))

setMethod("show",
          signature(object="metaSite"),
          function(object) {
            if(length(object@sites) > 0){
              cat("metaSite object with ", nrow(object@cellGeneCount)," cells and ",
                  ncol(object@cellGeneCount), "meta splicing sites, which are merged from ",
                  length(unlist(object@sites)), "original splicing sites.\n")
              cat("The head of the gene count matrix for each meta site is:\n")
              print(head(object@cellGeneCount))
              cat("The head of the spliced in count matrix for each meta site is:\n")
              print(head(object@cellSiteCount))
            }
            else{
              cat("Empty metaSite object.")
            }
          })




#' @title The Splice Class
#' @description The Splice class is an intermediate data storage class that stores the isoforms and other
#' related information needed for performing downstream analyses: including highly variable
#' exons detection and differential alternative splicing analysis
#' @slot cells The vector to store cell names
#' @slot genes The vector to store gene names.
#' @slot isoforms The vector to store isoform names.
#' @slot count A data frame to store the expression, polyA and other related information
#' for each isoform in each cell
#'
#' @name Splice-class
#' @rdname Splice-class
#' @concept objects
#' @exportClass Splice
#'
setClass("Splice",
         representation(cells = "data.frame",
                        genes = "data.frame",
                        isoforms = "character",
                        count = "data.frame",
                        meta_sites = "list"),
         prototype(cells = NULL, genes = NULL,
                   isoforms = NULL,count = NULL,
                   meta_sites = NULL))

setMethod("show",
          signature(object="Splice"),
          function(object) {
            cat("Splice object with ", nrow(object@cells),
                " cells and ", nrow(object@genes), " genes\n")
            cat("The top 10 cells are:", paste(head(object@cells, 10)$cell, collapse = ","), "\n")
            cat("The top 10 genes are:", paste(head(object@genes, 10)$gene, collapse = ","), "\n")
          })

#' @title createSplice
#'
#' @description Create a Splice object from the output file from saveExonList().
#' @details Create a Splice object from the output file from saveExonList().
#'
#' @param path The path to store the output files from saveExonList().
#' @param cell The vector of cells to be preserved in the Splice object
#' @param gene The vector of genes to be preserved in the Splice object
#' @param cell_col The colname to for the cell id in the count file
#' @param gene_col The colname to for the gene id in the count file
#' @param iso_col The colname to for the isoform id in the count file
#' @importFrom utils read.table
#'
createSplice <- function(path,cell = "all", gene = "all",
                               cell_col = "cell",gene_col = "gene",iso_col = "exons"){
  cells = scan(paste(path,"barcodes.tsv",sep = "/"),character(), quote = "")
  genes = scan(paste(path,"features.tsv",sep = "/"),character(), quote = "")
  isoforms = scan(paste(path,"isoforms.tsv",sep = "/"),character(), quote = "")
  count = read.table(paste(path,"isoform_count.tsv",sep = "/"),header = TRUE)

  if(cell[1] != "all"){
    cell_id = which(cells %in% cell)
    cells = cells[cell_id]
    count = count[count[,cell_col] %in% cell_id,]
  }
  if(gene[1] != "all"){
    gene_id = which(genes %in% gene)
    genes = genes[gene_id]
    count = count[count[,gene_col] %in% gene_id,]
  }
  cells = cells[which(1:length(cells) %in% count[,cell_col])]
  genes = genes[which(1:length(genes) %in% count[,gene_col])]
  isoforms = isoforms[which(1:length(isoforms) %in% count[,iso_col])]
  count[,cell_col] = as.numeric(as.factor(count[,cell_col]))
  count[,gene_col] = as.numeric(as.factor(count[,gene_col]))
  count[,iso_col] = as.numeric(as.factor(count[,iso_col]))

  splice_object = new("Splice",
               cells = cells,genes = genes,
               isoforms = isoforms, count = count,
               meta_sites = list())
  return(object)
}

#' @title creatSplice_from_df
#'
#' @description Create a Splice object from the output file from a dataframe.
#' @details Create a Splice object from the output file from dataframe, the dataframe should
#' contain at least four columns, including cell, gene, isoform and umi count.
#'
#' @param df The input dataframe.
#' @inheritParams createSplice
#' @importFrom dplyr mutate_at group_by summarise
#' @importFrom magrittr %>%
#' @return A Splice object
#' @export
#'
creatSplice_from_df <- function(df,cell = "all",gene = "all",
                                cell_col = "cell",gene_col = "gene",
                                iso_col = "isoform",exprs_col = "count"){
  df = as.data.frame(df)
  nece = c(cell_col,gene_col,iso_col,exprs_col)
  meta = setdiff(colnames(df),nece)
  df = df[,c(nece,meta)]
  colnames(df)[1:length(nece)] = c("cell","gene","isoform","count")

  if(cell[1] != "all"){
    df = df[df$cell %in% cell,]
  }
  if(gene[1] != "all"){
    df = df[df$gene %in% gene,]
  }

  cells = df %>% group_by(cell) %>% summarise(nUMI = sum(count),nGene = length(unique(gene)))
  genes = df %>% group_by(gene) %>% summarise(count = sum(count),nCell = length(unique(cell)))
  isoforms = names(table(df$isoform))

  df = df %>% mutate_at(c("cell","gene","isoform"),~as.numeric(as.factor(.)))

  splice_object = new("Splice",
                      cells = cells,genes = genes,
                      isoforms = isoforms, count = df,
                      meta_sites = list())
  return(splice_object)
}

#' @title getMetaSites
#' @description Extract the metaSites object for a gene for a subset of cells from the splice object
#' @param spliceOb the Splice object
#' @param gene name of the target gene
#' @param cells cells to use in building exon count table
#' @param verbose The flag to indicate if output message
#' @return a metaSite object
#' @export
getMetaSites = function(spliceOb,gene,cells = "all",verbose = TRUE){
  if(is.null(spliceOb@meta_sites[[gene]])){
    stop("The set of meta splicing sites for gene ", gene,
         " hasn't been built yet, please run geneSiteTable() first!")
  }

  meta_sites = spliceOb@meta_sites[[gene]]
  if(cells[1] == "all"){
    cells = rownames(meta_sites@cellGeneCount)
  }
  else{
    diff = setdiff(cells,rownames(meta_sites@cellGeneCount))
    if(length(diff) > 0 & verbose){
      warning(paste(head(diff,10),collapse = ","),"... don't have enough expression for this gene ",
              gene," will be ignored!")
    }
    cells = cells[cells %in% rownames(meta_sites@cellGeneCount)]
  }

  meta_sites@cellGeneCount = meta_sites@cellGeneCount[cells,,drop=FALSE]
  meta_sites@cellSiteCount = meta_sites@cellSiteCount[cells,,drop=FALSE]

  return(meta_sites)
}

#' @title extract isoforms
#' @description  extract isoforms for multigenes from selected cells
#' @param spliceob A Splice object
#' @param genes A vector of genes to extract isoforms
#' @param cells A vector of cells to choose to extract isoforms. Default to choose
#' all cells
#' @return A dataframe which store the expression and related information for each
#' isoform.
getIsoform.base <- function(spliceob,genes, cells = "all"){
  isoforms = spliceob@count

  if(class(spliceob@cells) == "character"){
    cells_pool = spliceob@cells
  }
  else{
    cells_pool = spliceob@cells$cell
  }

  if(class(spliceob@genes) == "character"){
    genes_pool = spliceob@genes
  }
  else{
    genes_pool = spliceob@genes$gene
  }

  if(cells[1] == "all"){
    cell_id = 1:length(cells_pool)
  }
  else{
    cell_id = which(cells_pool %in% cells)
  }

  gene_id = which(genes_pool %in% genes)

  select_iso = isoforms[isoforms$gene %in% gene_id &
                          isoforms$cell %in% cell_id,]

  select_iso$gene = genes_pool[select_iso$gene]
  select_iso$cell = cells_pool[select_iso$cell]
  select_iso$isoform = spliceob@isoforms[select_iso$isoform]

  return(select_iso)
}

#' @title generic getIsoform function definition
#' @param object the Splice object
#' @inheritParams getIsoform.base
#' @param ... other possible parameters for getIsoform.base
#' @export
setGeneric("getIsoform",
           function(object,genes,...) standardGeneric("getIsoform"))

#' @title generic getIsoform function for Splice object
#' @description  extract isoforms for multigenes from selected cells for Splice object
#' @param object a Splice object
#' @inheritParams getIsoform.base
#' @param ... parameters for getIsoform.base
#' @return A dataframe which store the expression and related information for each
#' isoform.
#'
#' @export
setMethod("getIsoform",
          signature(object = "Splice",genes = "character"),
          function(object,genes,...){
            getIsoform.base(object,genes = genes,...)
          }
)

#' @title metasite_retrieve
#' @description  recover the metasites to the corresponding original sites for multiple genes
#' @param spliceob A Splice object
#' @param genes A vector of genes to extract isoforms
#' @param sites A vector of meta splicing sites, usually has format as "metasite_id"
#' @return A vector of string, each element is the pasted splicing sites for a meta site.
metasite_retrieve = function(splice_ob,genes,sites){
  if(length(genes) != length(sites)){
    stop("The length of genes and sites should correspond!")
  }
  n = length(genes)

  out = sapply(1:n,function(i){
    gene_i = genes[i]
    sites_i = sites[i]
    print(gene_i)

    site_i_id = as.numeric(gsub("metasite_","",sites_i))

    meta = getMetaSites(splice_ob,gene_i)
    result = meta@sites[[site_i_id]]
    result = paste(result,collapse = ",")
    return(result)
  })
  return(out)
}
