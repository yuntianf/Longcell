#' @include spliceOb.R
#' @include exon_table.R


#' @title negative log likelihood for exon psi
#' @description compute log likelihood for exon count and gene count given a
#' beta-binomial distribution
#' @param exon_count a vector of exon count
#' @param n a vector of gene count with the same length as exon_count
#' @param alpha alpha for beta distribution
#' @param beta beta for beta distribution
#' @return a numerical value as negative log likelihood
#' @importFrom extraDistr dbbinom
ll <- function(exon_count,n,alpha, beta) {
  temp = sapply(1:length(n), function(i) {
    -suppressWarnings(extraDistr::dbbinom(exon_count[i], n[i], alpha, beta, log = TRUE))
  })
  return(sum(temp))
}

#' @title mle for beta-binomial
#' @description mle for beta-binomial with multiple start points to avoid local optimization
#' @param exon_count a vector of exon count
#' @param gene_count a vector of gene count with the same length as exon_count
#' @param start a list of start points, each line should be a pair of alpha and beta
#' @return a list of result for mle
#' @importFrom bbmle mle2 summary
mle_multi_start <- function(exon_count,gene_count,start){
    test <- lapply(start,function(x){
        x = as.numeric(x)
        error <- try({
          H <- suppressWarnings(bbmle::mle2(ll,start = list(alpha = x[1],beta = x[2]),
                    data = list(exon_count = exon_count,n = gene_count)))
            })
        if(class(error) == "try-error"){
            return(NULL)
        }
        return(H)
    })
    test[vapply(test,is.null,logical(1L))] <- NULL
    if(length(test) == 0){
        return(NULL)
    }
    else{
        log_lik <- sapply(test,function(x){
            return(x@min)
        })
        p_value <- sapply(test,function(x){
          test_summary = bbmle::summary(x)
          p_value = colMeans(test_summary@coef)["Pr(z)"]
        })
        id = order(p_value,log_lik)[1]
        return(test[[id]])
    }
}


#' @title Generalized likelihood ratio test.
#' @description Using beta-binomial to model the exon count and do
#' generalized likelihood ratio test for two cell populations to test if they have the
#' same psi distribution.
#' @param exon_count1 a vector of exon count for cell population 1
#' @param gene_count1 a vector of gene count for cell population 1
#' @param exon_count2 a vector of exon count for cell population 2
#' @param gene_count2 a vector of gene count for cell population 2
#' @param iters the number of iterations for bootstrap
#' @param psi_num the sample size from two beta distributions to calculate the earth-moving distance
#' @param p_thresh the threshold for the p-value
#' @importFrom stats pchisq
#' @return a vector of p-value and parameters for the psi distributions in two
#' cell populations if they have significant difference
GLRT <- function(exon_count1,gene_count1,exon_count2,gene_count2,
                 iters = 100,psi_num = 500,p_thresh = 0.05){
    if(length(exon_count1) != length(gene_count1) |
       length(exon_count2) != length(gene_count2)){
        stop("The length of exon count and gene count should be the same!")
    }

    start = list(c(0.5,0.5),c(2,2),c(1,3),c(3,1))
    H1_1 <- mle_multi_start(exon_count1,gene_count1,start = start)
    H1_2 <- mle_multi_start(exon_count2,gene_count2,start = start)
    H0 <- mle_multi_start(c(exon_count1,exon_count2),
                          c(gene_count1,gene_count2),
                          start = start)

    if(is.null(H1_1) | is.null(H1_2) | is.null(H0)){
        print("optimization failed!")
        return(NULL)
    }

    chi = 2*(H0@min-H1_1@min-H1_2@min)
    p= pchisq(chi,2,lower.tail = F)

    #if(p <= p_thresh){
    #    dis = beta_dis(exon_count1,gene_count1,exon_count2,gene_count2,
    #                   iters = iters,psi_num = psi_num)
    #}
    #else{
    #    dis = rep(NA,3)
    #}
    return(c(p,H1_1@coef,H1_2@coef))
}

#' @title Generalized likelihood ratio test for a gene.
#' @description an interface function to run GLRT for specified exons in a gene
#' for specified two cell groups from a Splice object.
#' @param spliceOb a Splice object
#' @param gene the name of specified gene
#' @param group1 cells for group 1
#' @param group2 cells for group 2
#' @param gene_bed a dataframe of gene bed annotation
#' @param gtf a dataframe of gtf annotation
#' @param exons specified exons to do GLRT
#' @param exon_len_thresh the minimum length of exons for GLRT
#' @param iters the number of iterations for bootstrap
#' @param psi_num the sample size from two beta distributions to calculate
#' the earth-moving distance
#' @param group_cell_thresh the minimum cell number for each group for GLRT
#' @param verbose whether to print out the iteration information
#' @param ... parameters specified for gene_exons_table.base()
#' @return a dataframe of the test result. Each line is an exon with p-values
#' and parameters for the psi distributions in two cell populations if they
#' have significant difference
gene_GLRT.base <- function(spliceOb,gene,group1,group2,
                      gene_bed = NULL,gtf = NULL,exons = "all",
                      exon_len_thresh = 10,
                      iters = 100,psi_num = 500,
                      group_cell_thresh = 30,
                      verbose = FALSE,...){
  count_list = gene_exons_table.base(spliceOb,gene,
                                     gene_bed = gene_bed,gtf = gtf,
                                     cells = c(group1,group2),exons = exons,
                                     exon_len_thresh = exon_len_thresh,...)

  gene_count_table = count_list[[1]]
  exon_count_table = count_list[[2]]

  gene_count_table1 = gene_count_table[gene_count_table$cell_id %in% group1,]
  gene_count_table2 = gene_count_table[gene_count_table$cell_id %in% group2,]
  exon_count_table1 = exon_count_table[exon_count_table$cell_id %in% group1,]
  exon_count_table2 = exon_count_table[exon_count_table$cell_id %in% group2,]

  if(nrow(gene_count_table1) < group_cell_thresh |
     nrow(gene_count_table2) < group_cell_thresh){
    return(NULL)
  }
  exons = setdiff(colnames(exon_count_table),c("cell_id","gene_count"))

  test = lapply(exons,function(i){
    exon_count1 = round(unlist(exon_count_table1[,i]))
    gene_count1 = round(unlist(gene_count_table1[,i]))
    exon_count1 = exon_count1[gene_count1 > 0]
    gene_count1 = gene_count1[gene_count1 > 0]

    exon_count2 = round(unlist(exon_count_table2[,i]))
    gene_count2 = round(unlist(gene_count_table2[,i]))
    exon_count2 = exon_count2[gene_count2 > 0]
    gene_count2 = gene_count2[gene_count2 > 0]

    if(length(gene_count1) > group_cell_thresh &
       length(gene_count2) > group_cell_thresh){
      psi = c(exon_count1,exon_count2)/c(gene_count1,gene_count2)
      if(verbose){
        cat("--- exon ",i,":","psi=",mean(psi),"\n")
      }
      if(mean(psi) >= 0.05 & mean(psi) <= 0.95){
        out = suppressWarnings(GLRT(exon_count1,gene_count1,exon_count2,gene_count2,
                   iters = iters,psi_num = psi_num))
        if(length(out) > 0){
          return(c(i,out))
        }
      }
    }
    return(NULL)
  })

  test = as.data.frame(do.call(rbind,test))
  if(length(test) > 0 && nrow(test) > 0){
    colnames(test) = c("exon","p","alpha1","beta1","alpha2","beta2")
    #                   "dis_left","dis","dis_right")
    return(test)
  }

  return(test)
}

#' @title generic gene_GLRT function definition
#' @param object the Splice or Seurat object
#' @inheritParams gene_GLRT.base
#' @param ... other possible parameters for gene_GLRT.base
#' @export
setGeneric("gene_GLRT",
           function(object,gene,group1,group2,...)
             standardGeneric("gene_GLRT"))

#' @title generic gene_GLRT function for Splice object
#' @description  an interface function to run GLRT for specified exons in a gene
#' for specified two cell groups from a Splice object.
#' @param object the Splice object
#' @inheritParams gene_GLRT.base
#' @return a dataframe of the test result. Each line is an exon with p-values
#' and parameters for the psi distributions in two cell populations if they
#' have significant difference
#' @export
setMethod("gene_GLRT",
          signature(object = "Splice",gene = "character",
                    group1 = "vector",group2 = "vector"),
          function(object,gene,group1,group2,...){
            gene_GLRT.base(spliceOb = object,gene = gene,
                                  group1 = group1,group2 = group2,...)
          }
)

#' @title generic gene_GLRT function for Seurat object
#' @description  an interface function to run GLRT for specified exons in a gene
#' for specified two cell groups from a Seurat object.
#' @param object the Seurat object
#' @param group_col if group_col is NULL, group1 and group2 should be cell names in
#' each group. If group_col is specified as a colname in Seurat cell meta data,
#' group1 and group2 should be group names.
#' @param group1 group name or cell name
#' @param group2 group name or cell name
#' @param slot the slot name to store Splice object in Seurat object
#' @inheritParams gene_GLRT.base
#' @import Seurat
#' @return a dataframe of the test result. Each line is an exon with p-values
#' and parameters for the psi distributions in two cell populations if they
#' have significant difference
#' @export
setMethod("gene_GLRT",
          signature(object = "Seurat",gene = "character",
                    group1 = "character",group2 = "character"),
          function(object,gene,group1,group2,
                   slot = "splice",group_col = NULL,...){
            spliceOb = getSplice(object,slot = slot)
            if(is.null(group_col)){
              warning("The group_col is not specified, you need to input cell name vectors for group1 and group2!")
              gene_GLRT.base(spliceOb,gene,group1,group2,...)
            }
            else{
              meta = object@meta.data
              cell_group1 = rownames(meta[meta[,group_col] %in% group1,])
              cell_group2 = rownames(meta[meta[,group_col] %in% group2,])
              gene_GLRT.base(spliceOb,gene = gene,
                             group1 = cell_group1,group2 = cell_group2,...)
            }
          }
)

#' @title calculate mean value for a beta distribution
#' @description calculate mean value for a beta distribution
#' @param alpha alpha for the beta distribution
#' @param beta beta for the beta distribution
#' @return a numeric value as the mean
beta_mean <- function(alpha,beta){
  return(alpha/(alpha+beta))
}

#' @title calculate variance for a beta distribution
#' @description calculate variance for a beta distribution
#' @param alpha alpha for the beta distribution
#' @param beta beta for the beta distribution
#' @return a numeric value as the variance
beta_var <- function(alpha,beta){
  return(alpha*beta/((alpha+beta)^2*(alpha+beta+1)))
}

#' @title distance between two psi distribution modeled by beta.
#' @description for significant difference of alternative splicing, this function calculates
#' the earth-moving distance between the psi distribution of two cell populations and its confidence interval
#' through bootstrap.
#' @param exon_count1 a vector of exon count for cell population 1
#' @param gene_count1 a vector of gene count for cell population 1
#' @param exon_count2 a vector of exon count for cell population 2
#' @param gene_count2 a vector of gene count for cell population 2
#' @param iters the number of iterations for bootstrap
#' @param psi_num the sample size from two beta distributions to calculate the earth-moving distance
#' @return a vector of the distance and its confidence interval
#' @importFrom transport wasserstein1d
#' @importFrom stats quantile rbeta
beta_dis <- function(exon_count1,gene_count1,exon_count2,gene_count2,
                     iters = 100, psi_num = 500){
  paras <- sapply(1:iters,function(x){
    id1 <- sample(1:length(gene_count1),length(gene_count1),replace = TRUE)
    id2 <- sample(1:length(gene_count2),length(gene_count2),replace = TRUE)

    start = list(c(0.5,0.5),c(2,2),c(1,3),c(3,1))

    H1_1 <- mle_multi_start(exon_count1[id1],gene_count1[id1],start = start)
    H1_2 <- mle_multi_start(exon_count2[id2],gene_count2[id2],start = start)

    if(is.null(H1_1) | is.null(H1_2)){
      return(NULL)
    }

    else{
      psi_sample1 <- rbeta(psi_num,H1_1@coef[1],H1_1@coef[2])
      psi_sample2 <- rbeta(psi_num,H1_2@coef[1],H1_2@coef[2])
      dis <- transport::wasserstein1d(psi_sample1,psi_sample2)
      return(dis)
    }
  })

  paras_conf <- quantile(unlist(paras),c(0.025,0.975))

  return(c(paras_conf[1],mean(unlist(paras)),paras_conf[2]))
}

#' @title Generalized likelihood ratio test for multiples gene.
#' @description an interface function to run GLRT for all exons in multiple genes
#' for multiple cell group pairs from a Splice object.
#' @param genes a vector of gene names to do GLRT
#' @param group1s a named list of cell names, each element should be cells belonging to a group
#' @param group2s a named list of cell names, each element should be cells belonging to a group
#' @param q_thresh the threshold for q value by fdr control
#' @param cores The number of cores used for parallization
#' @param verbose if print the process information
#' @param ... other parameters for gene_GLRT.base
#' @inheritParams gene_GLRT.base
#' @import dplyr
#' @importFrom stats p.adjust
#' @importFrom parallel mclapply detectCores
#' @return a dataframe of the test result. Each line is an exon with p-values
#' and parameters for the psi distributions in two cell populations if they
#' have significant difference
genes_groups_GLRT.base <- function(spliceOb,genes,group1s,group2s,
                              gene_bed = NULL,gtf = NULL,
                              exon_len_thresh = 10,
                              iters = 100,psi_num = 500,
                              group_cell_thresh = 30,q_thresh = 0.05,
                              cores = 4,verbose = FALSE,...){
  core_avai = parallel::detectCores()
  cores = ifelse(cores > core_avai,core_avai,cores)
  cat("The job will be paralleled on ",cores," cores.\n")
  #results <- lapply(genes,function(i){
  results <- parallel::mclapply(genes,function(i){
        group1 <- lapply(1:length(group1s),function(x){
            group2 <- lapply(1:length(group2s),function(y){
                error <- try({
                  test = gene_GLRT.base(spliceOb,gene = i,
                                        group1 = group1s[[x]],group2 = group2s[[y]],
                                        gene_bed = gene_bed,gtf = gtf,
                                        exon_len_thresh = exon_len_thresh,
                                        iters = iters,psi_num = psi_num,
                                        group_cell_thresh = group_cell_thresh,
                                        verbose = verbose,...)
                  },silent = TRUE)

                group_name1 = names(group1s)[x]
                group_name2 = names(group2s)[y]

                if(is.null(group_name1) || is.null(group_name2)){
                  warnings("The name for the cell group are not provided, will set with null!\n")
                }

                if(class(error) == "try-error"){
                  cat("Something wrong happened for GLRT of gene ",i,
                      " for groups ",group_name1,"-",group_name2,"!\n")
                  return(NULL)
                }

                if(length(test) >0 && nrow(test) >0){
                  test = cbind(group_name1,group_name2,test)
                  colnames(test)[c(1,2)] = c("group1","group2")
                }
                return(test)
            })
            group2 <- do.call(rbind,group2)
        })
        group1 <- as.data.frame(do.call(rbind,group1))

        if(length(group1) > 0){
          group1 = cbind(i,group1)
          colnames(group1)[1] = "gene"
        }
        return(group1)
    }, mc.cores = getOption("mc.cores", cores))
  #})
    results <- as.data.frame(do.call(rbind,results))
    if(length(results) >0 && nrow(results) > 0){
      results <- results %>% mutate_at(vars(-all_of(c("gene","exon","group1","group2"))),
                                       as.numeric)
      results$q = p.adjust(p = results$p,method = "fdr",n = nrow(results))
      results = results[results$q <= q_thresh,]

      results$mean_diff <- beta_mean(results$alpha2,results$beta2) -
        beta_mean(results$alpha1,results$beta1)
      results$var_diff <- beta_var(results$alpha2,results$beta2) -
        beta_var(results$alpha1,results$beta1)
    }
    return(results)
}

#' @title generic genes_groups_GLRT function definition
#' @param object the Splice or Seurat object
#' @inheritParams genes_groups_GLRT.base
#' @param ... other possible parameters for genes_groups_GLRT.base
#' @export
setGeneric("genes_groups_GLRT",
           function(object,genes,group1s,group2s,...) standardGeneric("genes_groups_GLRT"))

#' @title generic genes_groups_GLRT function for Splice object
#' @description  an interface function to run GLRT for all exons in multiple genes
#' for multiple cell group pairs from a Splice object.
#' @param object the Splice object
#' @inheritParams genes_groups_GLRT.base
#' @return a dataframe of the test result. Each line is an exon with p-values
#' and parameters for the psi distributions in two cell populations if they
#' have significant difference
#' @export
setMethod("genes_groups_GLRT",
          signature(object = "Splice",genes = "character",
                    group1s = "character",group2s = "character"),
          function(object,genes,group1s,group2s,...){
            genes_groups_GLRT.base(spliceOb = object,genes = genes,
                                  group1s = group1s,group2s = group2s,...)
          }
)

#' @title generic genes_groups_GLRT function for Seurat object
#' @description  an interface function to run GLRT for all exons in multiple genes
#' for multiple cell group pairs from a Splice object.
#' @param object the Seurat object
#' @param group_col if group_col is NULL, group1s and group2s should be cell names in
#' each group. If group_col is specified as a colname in Seurat cell meta data,
#' group1s and group2s should be group names.
#' @param group1s group names or cell names
#' @param group2s group names or cell names
#' @param slot the slot name to store Splice object in Seurat object
#' @inheritParams genes_groups_GLRT.base
#' @param ... parameters for genes_groups_GLRT.base
#' @return a dataframe of the test result. Each line is an exon with p-values
#' and parameters for the psi distributions in two cell populations if they
#' have significant difference
#' @export
setMethod("genes_groups_GLRT",
          signature(object = "Seurat",genes = "character",
                    group1s = "character",group2s = "character"),
          function(object,genes,group1s,group2s,
                   group_col = NULL,slot = "splice",...){
            spliceOb = getSplice(object,slot = slot)
            if(is.null(group_col)){
              warning("The group_col is not specified, you need to input cell name vectors for group1 and group2!")
              out = genes_groups_GLRT.base(spliceOb = spliceOb,genes = genes,
                                     group1s = group1s,group2s = group2s,...)
            }
            else{
              meta = object@meta.data
              cell_group1 = lapply(group1s,function(x){
                rownames(meta[meta[,group_col] %in% x,])
              })
              names(cell_group1) = group1s
              cell_group2 = lapply(group2s,function(x){
                rownames(meta[meta[,group_col] %in% x,])
              })
              names(cell_group2) = group2s
              out = genes_groups_GLRT.base(spliceOb,genes = genes,
                                    group1s = cell_group1,group2s = cell_group2,...)
            }
            return(out)
          }
)


#' @title plot for GLRT significant results
#' @description  plot the scatter plot for mean change versus variance change of exon psi
#' between two cell groups
#' @param data the GLRT result, a dataframe output from genes_groups_GLRT()
#' @param gene_col The name of column in the data which stores the gene name
#' @param mean_diff The name of column in the data which stores the mean difference
#' of the psi distributions
#' @param var_diff The name of column in the data which stores the variance difference
#' of the psi distributions
#' @param q The name of column in the data which stores the q value
#' @param q_thresh the minimum of q to show the gene name in the scatter plot
#' @param q_mid The scatter plot is colored by q value, q_mid is set for the middle color point
#' @param pt.size The point size in the scatter
#' @param text.size The size of the gene name annotation on the point
#' @param low,mid,high The color to indicate the q value
#' @import ggplot2
#' @importFrom latex2exp TeX
#' @return a ggplot object
#' @export
GLRT_sig_plot = function(data,gene_col = "gene",
                         mean_diff = "mean_diff",var_diff = "var_diff",q = "q",
                         q_thresh = 3,q_mid = NULL,pt.size = 2,text.size = 3,
                         low = "steelblue",mid = "whitesmoke",high = "tomato"){
  data$log_q = -log(data[,q],10)
  if(is.null(q_mid)){
    q_mid = (max(data$log_q)+min(data$log_q))/2
  }

  out = ggplot()+
    geom_point(data = data,aes_string(x = mean_diff,y=var_diff,color ="log_q"),size = pt.size)+
    scale_color_gradient2(low = low,mid = mid,high = high,midpoint = q_mid)+
    geom_text(data = data[data$log_q > q_thresh,],
              aes_string(x = mean_diff,y = var_diff,label = gene_col),
              nudge_x = 0, nudge_y = 0, check_overlap = T,size =text.size)+
    xlab(latex2exp::TeX("mean change of $\\psi$"))+ylab(latex2exp::TeX("var change of $\\psi$"))+
    theme_classic()+
    theme(text = element_text(size = 15))

  return(out)
}
