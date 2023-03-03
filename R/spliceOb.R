#' @import methods

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
         representation(cells = "vector",
                        genes = "vector",
                        isoforms = "vector",
                        count = "data.frame"),
         prototype(cells = NULL, genes = NULL,isoforms = NULL,count = NULL))

setMethod("show",
          signature(object="Splice"),
          function(object) {
            cat("Splice object with ", length(object@cells),
                " cells and ", length(object@genes), "genes")
          })

#' @title Create a Splice object
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
#' @export
#'
createSpliceObject <- function(path,cell = "all", gene = "all",
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

  object = new("Splice",
               cells = cells,genes = genes,
               isoforms = isoforms, count = count)
  return(object)
}

#' @title Create a Seurat object with Splice object
#'
#' @description Create a Seurat and Splice object from the output file from saveExonList().
#' @details The Seraut object is created to store the gene expression information,
#' and the Splice object for isoform information. The Splice object is embedded in the
#' misc slot in Seurat object
#'
#' @param path The path to store the output files from saveExonList()
#' @param project The name of the project for this object
#' @param assay The assay name for the count matrix
#' @param min.cells The minimum number of cells to express a gene
#' @param min.features The minimum number of genes cell to express
#' @param names.field For the initial identity class for each cell, choose this field from the cell's name.
#' @param names.delim For the initial identity class for each cell, choose this delimiter from the cell's column name
#' @param meta.data Additional cell-level metadata to add to the Seurat object.
#' @param ... other parameters for Seurat::CreateSeuratObject
#' @import Seurat
#' @importFrom Matrix sparseMatrix
#' @importFrom utils read.table
#' @export
#'
createSpliceSeurat <- function(path,project = "SeuratProject", assay = "RNA",
                               min.cells = 3, min.features = 20, names.field = 1,
                               names.delim = "_", meta.data = NULL,...){
  cells = scan(paste(path,"barcodes.tsv",sep = "/"),character(), quote = "")
  genes = scan(paste(path,"features.tsv",sep = "/"),character(), quote = "")
  count = read.table(paste(path,"matrix.mtx",sep = "/"))
  colnames(count) = c("cell","gene","count")

  count$cell = as.factor(count$cell)
  count$gene = as.factor(count$gene)
  count_matrix = Matrix::sparseMatrix(j=as.numeric(count$cell),
                              i=as.numeric(count$gene),
                              x=as.numeric(count$count),
                              dimnames=list(genes,cells))

  object = Seurat::CreateSeuratObject(counts = count_matrix,...)

  spliceOb = createSpliceObject(path,
                                cell = colnames(object),
                                gene = rownames(object))
  Seurat::Misc(object,slot = "splice") = spliceOb
  return(object)
}

#' @title extract the Splice object from the Seurat object
#' @description  extract the Splice object from the Seurat object
#' @param object A seurat object with a Splice object embedded in
#' @param slot the slot name to store Splice object in Seurat object
#' @return A Splice object.
#' @import Seurat
#' @export
getSplice <- function(object,slot = "splice"){
  splice = Misc(object)[[slot]]
  return(splice)
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

  if(cells[1] == "all"){
    cell_id = 1:length(spliceob@cells)
  }
  else{
    cell_id = which(spliceob@cells %in% cells)
  }
  gene_id = which(spliceob@genes %in% genes)

  select_iso = isoforms[isoforms$gene %in% gene_id &
                          isoforms$cell %in% cell_id,]

  select_iso$gene = spliceob@genes[select_iso$gene]
  select_iso$cell = spliceob@cells[select_iso$cell]
  select_iso$exons = spliceob@isoforms[select_iso$exons]

  return(select_iso)
}

#' @title generic getIsoform function definition
#' @param object the Splice or Seurat object
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

#' @title generic getIsoform function for Seurat object
#' @description  extract isoforms for multigenes from selected cells for Seurat object
#' @param object a Seurat object
#' @param slot the slot name to store Splice object in Seurat object
#' @inheritParams getIsoform.base
#' @param ... parameters for getIsoform.base
#' @import Seurat
#' @return A dataframe which store the expression and related information for each
#' isoform.
#'
#' @export
setMethod("getIsoform",
          signature(object = "Seurat",genes = "character"),
          function(object,genes, slot = "splice",...){
            spliceOb = getSplice(object,slot = slot)
            getIsoform.base(spliceOb,genes = genes,...)
          }
)
