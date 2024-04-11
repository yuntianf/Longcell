#' @include spliceOb.R
#' @include splice_table.R


#' @title geneSitesPsi.base
#' @description compute psi for sites in a gene from a Splice object
#' @param spliceOb the Splice object
#' @param gene name of the target gene
#' @param cells cells to use in building exon count table
#' @param exprs_thresh The minimum threshold for the gene expression to estimate the psi
#' @importFrom dplyr select_at
#' @importFrom magrittr %>%
#' @return a dataframe contains the psi estimation for each exon in each cell
#' @export
geneSitesPsi.base <- function(spliceOb,gene,
                                cells = "all",
                                exprs_thresh = 10,
                                ...){
  meta_sites = getMetaSites(spliceOb,gene,cells)

  cell_gene_count = meta_sites@cellGeneCount
  cell_site_count = meta_sites@cellSiteCount

  cell_site_psi = cell_site_count/cell_gene_count
  cell_site_psi[cell_gene_count < exprs_thresh] = NA

  row_filter = rowSums(is.na(cell_site_psi)) < ncol(cell_site_psi)
  col_filter = colnames(cell_site_psi)[colSums(is.na(cell_site_psi)) < nrow(cell_site_psi)]

  cell_site_psi = cell_site_psi[row_filter,] %>% select_at(col_filter)

  return(cell_site_psi)
}

#' @title generic genes_exons_psi function definition
#' @param object the Splice or Seurat object
#' @inheritParams  geneSitesPsi.base
#' @param ... other possible parameters for geneSitesPsi.base
#' @export
setGeneric("geneSitesPsi",
           function(object,gene,...)
             standardGeneric("geneSitesPsi"))

#' @title generic geneSitesPsi function for Splice object
#' @param object A Splice object
#' @inheritParams geneSitesPsi.base
#' @param ... other possible parameters for geneSitesPsi.base
#' @export
setMethod("geneSitesPsi",
          signature(object = "Splice",gene="character"),
          function(object,gene,...){
            geneSitesPsi.base(spliceOb = object,gene = gene,...)
          }
)


#' @title phi_from_psi
#' @description calculate the alternative splicing heterogeneity for an exon based on psi estimation
#' @param x a vector of exon psi across cells
#' @importFrom stats var
#' @return a numeric value for exon phi
phi_from_psi <- function(x){
  return(var(x)/(mean(x)*(1-mean(x))))
}

#' @title phi_from_beta
#' @description phi estimation from beta distribution
#' @param alpha alpha parameter for beta distribution
#' @param beta beta parameter for beta distribution
#' @return a numeric value for exon phi
phi_from_beta <- function(alpha,beta){
  return(1/(alpha+beta+1))
}


#' @title phi_conf_bootstrap
#' @description using bootstrap to calculate the confidence interval for exon phi
#' @param x a vector of exon psi across cells
#' @param iters the number of sampling iteration in bootstrap
#' @importFrom stats na.omit quantile
#' @return a numeric vector of the confidence interval
phi_conf_bootstrap <- function(x,iters = 100){
  conf <- sapply(1:iters,function(i){
    sample_x <- sample(x,length(x),replace = TRUE)
    sample_phi <- phi(sample_x)
    return(sample_phi)
  })
  out <- na.omit(conf)
  return(c(mean(out),quantile(out,c(0.025,0.975),na.rm = TRUE)))
}

#' @title phiNabla
#' @description The Nabla for phi estimation from alpha and beta
#' @inheritParams phi_from_beta
#' @return a numerical vector for the first derivative for alpha and beta
phiNabla = function(alpha,beta){
  return(rep(-1/(alpha+beta+1)^2,2))
}

#' @title phi_conf_delta
#' @description using delta method to estimate the confidence interval for exon phi
#' @inheritParams phi_from_beta
#' @param covmat The covariance matrix for alpha and beta
#' @param a The range of the confidence interval
#' @return a numeric vector of the confidence interval
phi_conf_delta <- function(alpha,beta,covmat,a = 0.05){
  phi_v = phiNabla(alpha,beta) %*% covmat %*% phiNabla(alpha,beta)
  CI = phi_from_beta(alpha,beta)+qnorm(c(a/2,(1-a/2)),0,sd = sqrt(phi_v))
  return(CI)
}


#' @title geneSitesPhi_from_psi
#' @description compute phi for exons in a gene from a Splice object
#' @param spliceOb the Splice object
#' @param gene name of the target gene
#' @param cells cells to use in building exon count table
#' @param cell_count_thresh the minimum number of cells expressing this gene
#' @param exprs_thresh the minimum gene expression to calculate the exon psi (for unbiased estimation)
#' @param iters iteration for bootstrap to generate the confidence interval of phi
#' @param ... Other parameters for geneSitesPsi.base
#' @importFrom dplyr mutate_at
#' @importFrom stats na.omit
#' @return a dataframe recording the phi estimation
geneSitesPhi_from_psi <- function(spliceOb,gene,cells = "all",
                           cell_count_thresh = 30,exprs_thresh = 10,
                           iters = 100,...){
  cell_site_psi = geneSitesPsi.base(spliceOb = spliceOb,gene = gene,
                                    cells = cells,exprs_thresh = exprs_thresh,...)
  if(is.null(cell_site_psi) || length(cell_site_psi) == 0){
    return(NULL)
  }
  sites = colnames(cell_site_psi)[which(colSums(!is.na(cell_site_psi)) >= cell_count_thresh)]

  if(length(sites) == 0){
    return(NULL)
  }
  else{
    sites_phi <- lapply(sites,function(i){
      psi = na.omit(unlist(cell_site_psi[,i]))
      conf = phi_conf(psi,iters = iters)
      return(c(i,mean(psi),conf,length(psi)))
    })
    sites_phi <- as.data.frame(do.call(rbind,sites_phi))
    colnames(sites_phi) <- c("meta_sites","mean_psi","phi","phi_lwr","phi_upr","count")
    sites_phi = sites_phi %>% mutate_at(setdiff(colnames(sites_phi),"meta_sites"),as.numeric)

    return(sites_phi)
  }
}

#' @title geneSitesPhi_from_beta
#' @description compute phi for exons in a gene from a Splice object
#' @param spliceOb the Splice object
#' @param gene name of the target gene
#' @param cells cells to use in building exon count table
#' @param cell_count_thresh the minimum number of cells expressing this gene
#' @param exon_len_thresh the minimum length of an exon to be preserved in exon count table
#' @param exprs_thresh the minimum gene expression to calculate the exon psi (for unbiased estimation)
#' @param CI_level The alpha for the confidence interval of phi estimation
#' @importFrom stats na.omit
#' @importFrom dplyr mutate_at
#' @return a dataframe recording the phi estimation
geneSitesPhi_from_beta = function(spliceOb,gene,cells = "all",
                                  mean_psi_lwr = 0.01,mean_psi_upr = 0.99,
                                  cell_count_thresh = 20,exprs_thresh = 5,
                                  CI_level = 0.05){
  meta_sites = getMetaSites(spliceOb,gene,cells)

  cell_gene_count = meta_sites@cellGeneCount
  cell_site_count = meta_sites@cellSiteCount

  sites = colnames(cell_site_count)

  sites_phi = lapply(sites,function(i){
    exon_count = cell_site_count[,i]
    gene_count = cell_gene_count[,i]

    exon_count = round(exon_count[gene_count >= exprs_thresh])
    gene_count = round(gene_count[gene_count >= exprs_thresh])

    if(length(gene_count) >= cell_count_thresh){
      psi = mean(na.omit(exon_count/gene_count))
      if(psi <= mean_psi_lwr | psi >= mean_psi_upr){
        return(NULL)
      }
      else{
        fit = mle_bb(exon_count,gene_count)
        if(fit@details$convergence == 1){
          return(NULL)
        }
        a = fit@coef[1]
        b = fit@coef[2]

        psi = a/(a+b)
        phi = phi_from_beta(a,b)
        phi_conf = phi_conf_delta(alpha = a,beta = b,covmat = fit@vcov,a = CI_level)
        return(c(i,psi,phi,phi_conf,length(gene_count)))
      }
    }
  })
  sites_phi = as.data.frame(do.call(rbind,sites_phi))
  colnames(sites_phi) <- c("meta_sites","mean_psi","phi","phi_lwr","phi_upr","count")
  sites_phi = sites_phi %>% mutate_at(setdiff(colnames(sites_phi),"meta_sites"),as.numeric)

  return(sites_phi)
}


#' @title compute phi for exons in multi-genes
#' @description compute phi for exons in a gene from a Splice object
#' @inheritParams geneSitesPhi_from_psi
#' @inheritParams geneSitesPhi_from_beta
#' @param psi_lwr,psi_upr The lower and upper bound to preserve as alternative spliced exons
#' @param phi_conf_thresh The maximum length of confidence interval to filter out inconfident
#' phi estimation
#' @param method The method to use for phi estimation
#' @param cores The number of cores for parallization
#' @param mode The mode of parallization
#' @importFrom stats na.omit
#' @importFrom parallel detectCores
#' @importFrom future sequential
#' @importFrom future multisession
#' @importFrom future multicore
#' @importFrom future.apply future_lapply
#' @importFrom dplyr mutate filter
#' @return a dataframe recording the phi estimation
#' @export
genesSitesPhi.base <- function(spliceOb,genes,
                                 cells = "all",
                                 cell_count_thresh = 20,
                                 exprs_thresh = 10,iters = 100,
                                 psi_lwr = 0.1, psi_upr = 0.9,
                                 phi_conf_thresh = 0.2,
                                 method = match.arg("beta","psi"),
                                 CI_level = 0.05,
                                 cores = 4,
                                 mode = c("sequential","multisession","multicore"),...){
  core_avai = parallel::detectCores()
  cores = ifelse(cores > core_avai,core_avai,cores)

  #   mode = match.arg(mode,c("sequential","multisession","multicore"))

  #   if(cores > 1 & mode == "sequential"){
  #     warning("Set more than one core to use, will use multisession mode instead!")
  #     mode = "multisession"
  #   }
  #   if(mode == "sequential"){
  #     plan(strategy = mode)
  #   }
  #   else{
  #     cat("The job will be paralleled on ",cores," cores.\n")
  #     plan(strategy = mode,workers = cores)
  #   }

  phi <- parallel::mclapply(genes,function(x){
    error <- try({
      if(method == "psi"){
        temp = geneSitesPhi_from_psi(spliceOb,x,
                            cells = cells,
                            cell_count_thresh = cell_count_thresh,
                            exprs_thresh = exprs_thresh,
                            iters = iters,...)
      }
      else if(method == "beta"){
        temp = geneSitesPhi_from_beta(spliceOb,x,cells,
                                      mean_psi_lwr = psi_lwr,mean_psi_upr = psi_upr,
                                      cell_count_thresh,exprs_thresh = exprs_thresh,
                                      CI_level = CI_level)
      }

      if(!is.null(temp)){
        temp = cbind(x,temp)
      }
      return(temp)
    })
    if(class(error) == "try-error"){
      print(x)
      return(NULL)
    }
  },mc.cores = getOption("mc.cores", cores))
  phi <- as.data.frame(do.call(rbind,phi))
  phi = na.omit(phi)

  if(length(phi) > 0){
    colnames(phi)[1] <- "gene"
    phi = phi %>% mutate(phi_conf = phi_upr-phi_lwr) %>%
                  filter(mean_psi > psi_lwr,mean_psi < psi_upr,phi_conf<phi_conf_thresh)
  }
  return(phi)
}

#' @title generic genesSitesPhi function definition
#' @param object the Splice or Seurat object
#' @param genes vector of names of the target gene
#' @param ... other possible parameters
#' @export
setGeneric("genesSitesPhi",
           function(object,genes,...)
             standardGeneric("genesSitesPhi"))

#' @title generic genesSitesPhi function for Splice object
#' @param object the Splice object
#' @inheritParams genesSitesPhi.base
#' @param ... parameters for genesSitesPhi.base
#' @export
setMethod("genesSitesPhi",
          signature(object = "Splice",genes = "character"),
          function(object,genes,...){
            genesSitesPhi.base.base(object,genes = genes,...)
          }
)



#' @title generic HighExprsGene function definition
#' @param object the Splice object
#' @param thresh the threshold for minimum total gene expression
#' @param ... other possible parameters
#' @export
setGeneric("HighExprsGene",
           function(object,thresh,...) standardGeneric("HighExprsGene"))

#' @title generic HighExprsGene function for Seurat object
#' @description choose genes with high total expression for phi calculation
#' @param object the Seurat object
#' @param slot the slot name RNA expression assay in the Seurat object
#' @param thresh the minimum total gene expression
#' @importFrom dplyr filter
#' @export
setMethod("HighExprsGene",
          signature(object = "Splice",thresh = "numeric"),
          function(object,thresh = 100){
            exprs = object@gene %>% filter(count > thresh)
            return(exprs$gene)
          }
)
