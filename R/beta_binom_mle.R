#' @title ll
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

#' @title mle_module
#' @description MLE estimation for a beta-binomial distribution with different number of total experiments
#' @param exon_count a vector of exon count
#' @param gene_count a vector of gene count with the same length as exon_count
#' @param start a list of start points, each line should be a pair of alpha and beta
#' @return a mle2 object
#' @importFrom bbmle mle2
mle_module = function(exon_count,gene_count,start = NULL){
  if(is.null(start)){
    start = bbmoment(exon_count,gene_count)
  }
  error <- try({
    H <- suppressWarnings(bbmle::mle2(ll,start = list(alpha = start[1],beta = start[2]),
                                      data = list(exon_count = exon_count,n = gene_count)))
  })
  if(class(error) == "try-error"){
    return(NULL)
  }
  return(H)
}


#' @title mle_bb
#' @description mle for beta-binomial with multiple start points to avoid local optimization
#' @inheritParams mle_module
#' @return a list of result for mle
#' @importFrom bbmle mle2 summary
#' @export
mle_bb <- function(exon_count,gene_count,start = NULL){
  H = mle_module(exon_count,gene_count,start)
  if(!is.null(H)){
    return(H)
  }
  else{
    print("The MLE based on current starting points doesn't converge, will try other starting points!")
    start = list(c(0.5,0.5),c(2,2),c(1,3),c(3,1))

    test = lapply(start,function(x){
      x = as.numeric(x)
      H = mle_module(exon_count,gene_count,x)
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
}

#' @title bbmoment
#' @description moment estimation for alpha and beta parameters in beta-binomial
#' @param success a vector for the number of success
#' @param total a vector for the number of total trials
#' @return a vector include the moment estimation for alpha and beta
bbmoment = function(success,total){
  psi = na.omit(success/total)
  m = mean(psi)
  v = var(psi)
  alpha = m*(m*(1-m)/v-1)
  beta = alpha*(1-m)/m

  return(c(alpha,beta))
}
