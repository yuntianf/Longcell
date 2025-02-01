#' @include spliceOb.R
#' @include splice_table.R
#' @include beta_binom_mle.R
#'

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

    H1_1 <- mle_bb(exon_count1,gene_count1)
    H1_2 <- mle_bb(exon_count2,gene_count2)
    H0 <- mle_bb(c(exon_count1,exon_count2),
                          c(gene_count1,gene_count2))

    if(is.null(H1_1) | is.null(H1_2) | is.null(H0)){
        print("optimization failed!")
        return(NULL)
    }

    chi = 2*(H0@min-H1_1@min-H1_2@min)
    p= pchisq(chi,2,lower.tail = F)

    emd = beta_emd(H1_1@coef[1],H1_1@coef[2],H1_2@coef[1],H1_2@coef[2])

    #if(p <= p_thresh){
    #    dis = beta_dis(exon_count1,gene_count1,exon_count2,gene_count2,
    #                   iters = iters,psi_num = psi_num)
    #}
    #else{
    #    dis = rep(NA,3)
    #}
    return(c(p,H1_1@coef,H1_2@coef,emd))
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
#' @param ... parameters specified for geneSiteTable()
#' @return a dataframe of the test result. Each line is an exon with p-values
#' and parameters for the psi distributions in two cell populations if they
#' have significant difference
gene_GLRT.base <- function(spliceOb,gene,group1,group2,
                           exons = "all",
                           iters = 100,psi_num = 500,
                           group_cell_thresh = 30,
                           verbose = FALSE,...){
  if(is.null(spliceOb@meta_sites[[gene]])){
    warning("The splice site table is not pre-computed,
            please run geneSiteTable for the target gene first!")
    return(NULL)
  }
  count_list = spliceOb@meta_sites[[gene]]

  gene_count_table = count_list@cellGeneCount
  exon_count_table = count_list@cellSiteCount

  group1 = group1[group1 %in% rownames(gene_count_table)]
  group2 = group2[group2 %in% rownames(gene_count_table)]

  gene_count_table1 = gene_count_table[group1,,drop=FALSE]
  gene_count_table2 = gene_count_table[group2,,drop=FALSE]
  exon_count_table1 = exon_count_table[group1,,drop=FALSE]
  exon_count_table2 = exon_count_table[group2,,drop=FALSE]

  if(nrow(gene_count_table1) < group_cell_thresh ||
     nrow(gene_count_table2) < group_cell_thresh){
    if(verbose){
      cat("The cell number for this gene is not enough, will skip!\n")
    }
    return(NULL)
  }
  if(exons == "all"){
    exons = colnames(exon_count_table)
  }
  else if(is.numeric(exons)){
    exons = colnames(exon_count_table)[exons]
  }
  else{
    exons = exons[exons %in% colnames(exon_count_table)]
  }

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
    colnames(test) = c("site","p","alpha1","beta1","alpha2","beta2","emd")
    #                   "dis_left","dis","dis_right")
    return(test)
  }

  return(test)
}

#' @title generic gene_GLRT function definition
#' @param object the Splice object
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

    H1_1 <- mle_bb(exon_count1[id1],gene_count1[id1])
    H1_2 <- mle_bb(exon_count2[id2],gene_count2[id2])

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

beta_emd = function(alpha1,beta1,alpha2,beta2){
  cdf1 <- function(x) pbeta(x, shape1 = alpha1, shape2 = beta1)
  cdf2 <- function(x) pbeta(x, shape1 = alpha2, shape2 = beta2)

  emd <- integrate(function(x) abs(cdf1(x) - cdf2(x)), lower = 0, upper = 1)$value
  return(emd)
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
#' @importFrom future sequential
#' @importFrom future multisession
#' @importFrom future multicore
#' @importFrom future.apply future_lapply
#' @return a dataframe of the test result. Each line is an exon with p-values
#' and parameters for the psi distributions in two cell populations if they
#' have significant difference
genes_groups_GLRT.base <- function(spliceOb,genes,group1s,group2s,
                              iters = 100,psi_num = 500,
                              group_cell_thresh = 30,q_thresh = 0.05,
                              filter = TRUE,verbose = FALSE,...){
  if(is.null(names(group1s))){
    names(group1s) = 1:length(group1s)
  }
  if(is.null(names(group2s))){
    names(group2s) = 1:length(group2s)
  }
  #results <- lapply(genes,function(i){
  results <- future_lapply(genes,function(i){
        group1 <- lapply(1:length(group1s),function(x){
            group2 <- lapply(1:length(group2s),function(y){
                error <- try({
                  test = gene_GLRT.base(spliceOb,gene = i,
                                        group1 = group1s[[x]],
                                        group2 = group2s[[y]],
                                        iters = iters,psi_num = psi_num,
                                        group_cell_thresh = group_cell_thresh,
                                        verbose = verbose,...)
                  },silent = TRUE)

                group_name1 = names(group1s)[x]
                group_name2 = names(group2s)[y]

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
    }, future.seed = TRUE,future.packages = c("Longcell"))
  #})
    results <- as.data.frame(do.call(rbind,results))
    if(length(results) >0 && nrow(results) > 0){
      results <- results %>% mutate_at(vars(-all_of(c("gene","site","group1","group2"))),
                                       as.numeric)
      if(filter == TRUE){
        results$q = p.adjust(p = results$p,method = "fdr",n = nrow(results))
        results = results[results$q <= q_thresh,]

        results$mean_diff <- beta_mean(results$alpha2,results$beta2) -
          beta_mean(results$alpha1,results$beta1)
        results$var_diff <- beta_var(results$alpha2,results$beta2) -
          beta_var(results$alpha1,results$beta1)
      }
    }

    return(results)
}

#' @title generic genes_groups_GLRT function definition
#' @param object the Splice object
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
                    group1s = "list",group2s = "list"),
          function(object,genes,group1s,group2s,...){
            genes_groups_GLRT.base(spliceOb = object,genes = genes,
                                  group1s = group1s,group2s = group2s,...)
          }
)


#' @title genes_multigroups_GLRT
#' @description  Do the generalized likelihood ratio test for multiples gene across multiple cell groups.
#' @param spliceOb the Splice object
#' @param genes target genes to be tested
#' @param meta a dataframe recording the meta data for each, should have at lease two columns,
#' one indicate the cell barcode and the other one indicate the cell group information.
#' @param cell_col,group_col The column names in the meta to record the cell/group info.
#' @param filter A bool flag to indicate if insignificant results should be filtered out.
#' @param q_thresh The threshold for the FDR controled q value.
#' @param ... Other parameters used in gene_GLRT.base.
#' @return a dataframe of the test result. Each line is an exon with p-values
#' and parameters for the psi distributions in two cell populations if they
#' have significant difference
#' @export
genes_multigroups_GLRT = function(spliceOb,genes,meta,cell_col = "cell",group_col = "group",
                                  filter = TRUE,q_thresh = 0.05,...){
  group_uniq = unique(meta[,group_col])

  results = lapply(group_uniq,function(x){
    group1 = meta %>% filter_at(group_col,~.== x)
    group2 = meta %>% filter_at(group_col,~.!= x)

    sub = future_lapply(genes,function(i){
      error <- try({
        test = gene_GLRT.base(spliceOb,gene = i,
                              group1 = unlist(group1[,cell_col]),
                              group2 = unlist(group2[,cell_col]),...)
      },silent = FALSE)
      # print(error)
      if(class(error) == "try-error"){
        cat("Something wrong happened for GLRT of gene ",i,
            " for group ",x,"!\n")
        return(NULL)
      }

      if(length(test) >0 && nrow(test) >0){
        test = cbind(i,test)
      }
      return(test)
    },future.seed = TRUE,future.packages = c("Longcell"))

    sub = as.data.frame(do.call(rbind,sub))

    if(length(sub) >0 && nrow(sub) >0){
      sub = cbind(x,sub)
      colnames(sub)[1:2] = c("group","gene")
    }

    return(sub)
  })

  results = do.call(rbind,results)

  if(length(results) >0 && nrow(results) > 0){
    results <- results %>% mutate_at(vars(-all_of(c("gene","site","group"))),
                                     as.numeric)
    if(filter == TRUE){
      results$q = p.adjust(p = results$p,method = "fdr",n = nrow(results))
      results = results[results$q <= q_thresh,]

      results$mean_diff <- beta_mean(results$alpha2,results$beta2) -
        beta_mean(results$alpha1,results$beta1)
      results$var_diff <- beta_var(results$alpha2,results$beta2) -
        beta_var(results$alpha1,results$beta1)
    }
  }

  return(results)
}

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
#' @param mean_diff_thresh the minimum of mean difference of psi to show the gene name in the scatter plot
#' @param var_diff_thresh the minimum of variance difference of psi to show the gene name in the scatter plot
#' @param q_thresh the minimum of -log(q,10) to show the gene name in the scatter plot
#' @import ggplot2
#' @importFrom latex2exp TeX
#' @importFrom ggrepel geom_text_repel
#' @return a ggplot object
#' @export
GLRT_sig_plot = function(data,gene_col = "gene",color_col = NULL,
                         mean_diff = "mean_diff",var_diff = "var_diff",q = "q",
                         mean_diff_thresh = 0.15,var_diff_thresh = 0.1,log_q_thresh = 3){
  data = as.data.frame(data)
  data$log_q = -log(data[,q],10)

  out = ggplot(mapping = aes(x = !!sym(mean_diff),y = !!sym(var_diff)))

  if(is.null(color_col)){
    out = out+geom_point(data = data,aes(size = log_q),color = "grey",alpha = 0.75)
    color_row = 0
  }
  else{
    out = out+geom_point(data = data,aes(color = !!sym(color_col),size = log_q),alpha = 0.75)
    color_row = length(unique(data[,color_col]))
  }


  out = out+
    ggrepel::geom_text_repel(data = data %>% filter(abs(mean_diff) > mean_diff_thresh|abs(var_diff) > var_diff_thresh,
                                                    log_q > log_q_thresh),
                             aes(label = !!sym(gene_col)),color = "black")+
    theme_classic()+
    xlab(TeX("mean change of $\\psi$"))+
    ylab(TeX("var change of $\\psi$"))+
    theme(text = element_text(size = 15))+
    theme(
      legend.position = "top",           # Move all legends to the top initially
      legend.box = "horizontal"          # Arrange legends horizontally
    ) +
    guides(
      color = guide_legend(order = 1,nrow = color_row,override.aes = list(size = 5)),   # Set color legend order to be first
      size = guide_legend(order = 2,nrow = 4)    # Set size legend order to be second
    )


  return(out)
}
