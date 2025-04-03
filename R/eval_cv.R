#' Function to compute the Population Average Value under SMART for cross-validated DTR.
#'
#' Given a A,D,Y dataset, this function estimates the Population Average Value of an DTR derived from the same data by CV.
#'
#' @param Acv A matrix of the actual treatment. Each row represents a unit-level randomized sequential assignments.
#' @param Dhatcv A list of DTR matrix derived from each trail of cross-validation.
#' @param Ycv A vector of the final outcome variable of interest for each sample.
#' @param indcv A vector of integers (between 1 and number of folds inclusive) indicating which testing set does each sample belong to.
#' @param centered If \code{TRUE}, the outcome variables would be centered before processing. This minimizes
#' the variance of the estimator. Default is \code{FALSE}.
#' @return A list that contains the following items: \item{pav}{The estimated population value of the estimated DTR.},
#'\item{sd}{The standard deviation of the estimated population value.}
#' @examples
#' Acv <- matrix(rep(0,1000*2),nrow = 1000,ncol= 2)
#' for(i in 1:2){
#'  n1=round(1000/2)
#'  n0=1000-n1
#'  ind=sample(1:1000,size=n1)
#'  Acv[ind,i] = 1}
#' Dhatcv <- list()
#' for (i in 1:5) {
#' D <- matrix(0, nrow = 1000, ncol = 2)
#' num_ones <- 500 + (i - 1)
#' D[1:num_ones, 2] <- 1
#' Dhatcv [[i]] <- D
#' }
#' Ycv =  rnorm(1000,3,1.5)
#' indcv <- sample(rep(1:5, each = 200))
#' pav_est <- cvPAV_dtr(Acv,Dhatcv,Ycv,indcv,centered = FALSE)
#' pav_est$pav
#' pav_est$sd
#' @export cvPAV_dtr
cvPAV_dtr <- function (Acv,Dhatcv,Ycv,indcv,centered = FALSE){
  t = ncol(Acv)
  n = nrow(Acv)
  nfolds = max(indcv)
  if (centered) {
    Ycv = Ycv - mean(Ycv)
  }
  M_Dhat = lapply(Dhatcv, function(x){
    M_d = matrix(nrow = nrow(x),ncol = 2^2)
    for (i in 1:nrow(x)) {
      M_d[i,] = dtrMatrix(x[i,],s=2)
    }
    return(M_d)
  })
  M_A = matrix(nrow = n,ncol = 2^t)
  for(i in 1:n) {
    # M[i,] = trtMatrix(A[i,],D[i,]) #follow dtr
    M_A[i,] = trtMatrix(Acv[i,],M_Dhat[[1]][i,],s=1) #follow randomization
  }
  n_A = colSums(M_A)
  p = n_A/n
  pavfold = c()
  Sf = 0
  #St = c()
  pav_ztr= matrix(nrow=nfolds,ncol = 2^t)
  #pav_true =  matrix(nrow=nfolds,ncol = 2^t)
  for (k in 1:nfolds) {
    a = Acv[indcv==k,]
    y = Ycv[indcv==k]
    dhat = Dhatcv[[k]][indcv==k,]
    output = PAV_dtr(a,dhat,y)
    pavfold = c(pavfold, output$pav)
    st = output$sd^2/nfolds
    Sf = Sf+st
    #St = c(St,output$St)
    pav = PAV_dtr(Acv,Dhatcv[[k]],Ycv)
    pav_ztr[k,]= pav$sav
  }
  cov_E = matrix(nrow = 2^t,ncol = 2^t)
  for(i in 1:2^t){
    for(j in 1:2^t){
      if(i ==j){cov_E[i,j] = var(pav_ztr[,i])}
      else {cov_E[i,j] = cov(pav_ztr[,i],pav_ztr[,j])}
    }
  }
  SF2 = var(pavfold)
  varexp = Sf + sum(cov_E)  - (nfolds-1)/nfolds*min(SF2,Sf+sum(cov_E))
  return(list(pav=mean(pavfold),sd=sqrt(max(varexp,0))))
}


#' Function to compute the Population Average Prescriptive Effect under SMART for cross-validated DTR.
#'
#' Given a A,D,Y dataset, this function estimates the  Population Average Prescriptive Effect of an DTR derived from the same data by CV.
#'
#' @param Acv A matrix of the actual treatment. Each row represents a unit-level randomized sequential assignments.
#' @param Dhatcv A list of DTR matrix derived from each trail of cross-validation.
#' @param Ycv A vector of the final outcome variable of interest for each sample.
#' @param indcv A vector of integers (between 1 and number of folds inclusive) indicating which testing set does each sample belong to.
#' @param centered If \code{TRUE}, the outcome variables would be centered before processing. This minimizes
#' the variance of the estimator. Default is \code{FALSE}.
#' @return A list that contains the following items: \item{pape}{The estimated pape of the estimated DTR.},
#'\item{sd}{The standard deviation of the estimated pape.}
#' @examples
#' Acv <- matrix(rep(0,1000*2),nrow = 1000,ncol= 2)
#' for(i in 1:2){
#'  n1=round(1000/2)
#'  n0=1000-n1
#'  ind=sample(1:1000,size=n1)
#'  Acv[ind,i] = 1}
#' Dhatcv <- list()
#' for (i in 1:5) {
#' D <- matrix(0, nrow = 1000, ncol = 2)
#' num_ones <- 500 + (i - 1)
#' D[1:num_ones, 2] <- 1
#' Dhatcv [[i]] <- D
#' }
#' Ycv =  rnorm(1000,3,1.5)
#' indcv <- sample(rep(1:5, each = 200))
#' pape_est <- cvPAPE_dtr(Acv,Dhatcv,Ycv,indcv,centered = FALSE)
#' pape_est$pape
#' pape_est$sd
#' @export cvPAPE_dtr

cvPAPE_dtr <- function (Acv,Dhatcv,Ycv,indcv,centered = FALSE){
  t = ncol(Acv)
  n = nrow(Acv)
  nfolds = max(indcv)
  if (centered) {
    Ycv = Ycv - mean(Ycv)
  }
  M_Dhat = lapply(Dhatcv, function(x){
    M_d = matrix(nrow = nrow(x),ncol = 2^2)
    for (i in 1:nrow(x)) {
      M_d[i,] = dtrMatrix(x[i,],s=2)
    }
    return(M_d)
  })
  M_A = matrix(nrow = n,ncol = 2^t)
  for(i in 1:n) {
    M_A[i,] = trtMatrix(Acv[i,],M_Dhat[[1]][i,],s=1)
  }
  n_A = colSums(M_A)
  p = n_A/n
  papefold = c()
  covfold = c()
  sav_ztr = matrix(nrow = nfolds,ncol = 2^t)
  rpav_ztr = matrix(nrow = nfolds,ncol = 2^t)
  Sf = 0
  for (k in 1:nfolds) {
    a = Acv[indcv==k,]
    y = Ycv[indcv==k]
    dhat = Dhatcv[[k]][indcv==k,]
    output = PAPE_dtr(a,dhat,y)
    papefold = c(papefold, output$pape)
    st = output$st/nfolds
    Sf = Sf + st
    pape =  PAPE_dtr(Acv,Dhatcv[[k]],Ycv)
    covfold = c(covfold,pape$cov)
    sav_ztr[k,] = (n-1)/n*pape$sav
    rpav_ztr[k,] = (n-1)/n*pape$rpav
  }
  mF = n/nfolds
  SF2 = var(papefold)
  cov_E = matrix(nrow =2^t, ncol = 2^t)
  for(i in 1:2^t){
    for(j in 1:2^t){
      if(i ==j){
        cov_E[i,j]=var(sav_ztr[,i])-cov(rpav_ztr[,i],sav_ztr[,j])-
          cov(sav_ztr[,i],rpav_ztr[,j])+var(rpav_ztr[,i])
      }
      else{
        cov_E[i,j]=cov(sav_ztr[,i],sav_ztr[,j])-cov(rpav_ztr[,i],sav_ztr[,j])-
          cov(sav_ztr[,i],rpav_ztr[,j])+cov(rpav_ztr[,i],rpav_ztr[,j])
      }
    }}
  var_ztr = mF^2/(mF-1)^2*(Sf+mean(covfold)+sum(cov_E))
  varexp = var_ztr-(nfolds-1)/nfolds*min(SF2,var_ztr)
  return(list(pape=mean(papefold),sd=sqrt(max(varexp,0))))
}


#' Function to compute the Local Population Average Value under SMART for cross-validated DTR.
#'
#' Given a A,D,Y dataset, this function estimates the Population Average Value of an DTR derived from the same data by CV up to stage s.
#'
#' @param Acv A matrix of the actual treatment. Each row represents a unit-level randomized sequential assignments.
#' @param Dhatcv A list of DTR matrix derived from each trail of cross-validation.
#' @param Ycv A vector of the final outcome variable of interest for each sample.
#' @param indcv A vector of integers (between 1 and number of folds inclusive) indicating which testing set does each sample belong to.
#' @param s The stage up to which DTR should be evaluated..
#' @param centered If \code{TRUE}, the outcome variables would be centered before processing. This minimizes
#' the variance of the estimator. Default is \code{FALSE}.
#' @return A list that contains the following items: \item{stage}{The stage up to which DTR has been evaluated}，
#' \item{lpav}{The estimated local population value of the estimated DTR.},
#' \item{sd}{The standard deviation of the estimated local population value.}
#' @examples
#' Acv <- matrix(rep(0,1000*2),nrow = 1000,ncol= 2)
#' for(i in 1:2){
#'  n1=round(1000/2)
#'  n0=1000-n1
#'  ind=sample(1:1000,size=n1)
#'  Acv[ind,i] = 1}
#' Dhatcv <- list()
#' for (i in 1:5) {
#' D <- matrix(0, nrow = 1000, ncol = 2)
#' num_ones <- 500 + (i - 1)
#' D[1:num_ones, 2] <- 1
#' Dhatcv [[i]] <- D
#' }
#' Ycv =  rnorm(1000,3,1.5)
#' indcv <- sample(rep(1:5, each = 200))
#' lpav_est <- cvLocal_PAV_dtr(Acv,Dhatcv,Ycv,indcv,2,centered = FALSE)
#' lpav_est$lpav
#' lpav_est$sd
#' @export cvLocal_PAV_dtr
cvLocal_PAV_dtr <- function(Acv,Dhatcv,Ycv,indcv,s,centered=FALSE){
  t = ncol(Acv)
  n = nrow(Acv)
  nfolds = max(indcv)
  if (centered) {
    Ycv = Ycv - mean(Ycv)
  }
  M_Dhat = lapply(Dhatcv, function(x){
    M_d = matrix(nrow = nrow(x),ncol = 2^t)
    for (i in 1:nrow(x)) {
      M_d[i,] = dtrMatrix(x[i,],s=t)
    }
    return(M_d)
  })
  M_A = matrix(nrow = n,ncol = 2^t)
  for(i in 1:n) {
    M_A[i,] = trtMatrix(Acv[i,],M_Dhat[[1]][i,],s=1)
  }
  n_A = colSums(M_A)
  p = n_A/n
  lpavfold = c()
  Sf = 0
  covfold = c()
  lsavloo1est_ztr= matrix(nrow=nfolds,ncol = 2^t)
  for (k in 1:nfolds) {
    a = Acv[indcv==k,]
    y = Ycv[indcv==k]
    dhat = Dhatcv[[k]][indcv==k,]
    output = Local_PAV_dtr(a,dhat,y,s)
    lpavfold = c(lpavfold, output$lpav)
    st = output$St/nfolds
    st = ifelse(is.na(st),0,st)
    Sf = Sf + st

    lpav = Local_PAV_dtr(Acv,Dhatcv[[k]],Ycv,s)
    covfold = c(covfold,lpav$cov)
    m_D = M_Dhat[[k]]
    m_S = matrix(nrow = n,ncol = 2^t)
    m_D_s <- matrix(nrow = n,ncol = 2^(s-1))
    for(i in 1:n) {
      m_S[i,] = trtMatrix(Acv[i,],Dhatcv[[k]][i,],s=s)
      m_D_s[i,] = dtrMatrix(Dhatcv[[k]][i,],s=(s-1))
    }
    m_pd_s = t(apply(m_D_s,1,function(x){
      M <- c()
      for(i in 1:length(x)){M <- c(M,rep(x[i],2^(t-s+1)))}
      return(M)}))
    mdSum <- colSums(m_D)
    mpdSum <- colSums(m_pd_s)
    m_L = matrix(nrow = n,ncol = 2^t)
    p_loo_1 <- matrix(nrow = n,ncol = 2^t)
    for(i in 1:n){
      p_loo_1[i,] <- (mdSum - m_D[i,])/(mpdSum-m_pd_s[i,])
      p_loo_1[i,] <- ifelse(is.na(p_loo_1[i,]),0,p_loo_1[i,])
      m_L[i,] = m_S[i,]*Ycv[i]*p_loo_1[i,]/p
    }
    loo1est = colMeans(p_loo_1)
    lsav = colSums(sweep(m_S, MARGIN=2,1/p, `*`) * Ycv)/n
    lsavloo1est_ztr[k,] = lsav*loo1est
  }
  cov_E = matrix(nrow = 2^t,ncol = 2^t)
  for(i in 1:2^t){
    for(j in 1:2^t){
      if(i ==j){cov_E[i,j]=var(lsavloo1est_ztr[,i])}
      else{cov_E[i,j]= cov(lsavloo1est_ztr[,i],lsavloo1est_ztr[,j])}
    }
  }
  SF2 = var(lpavfold)
  varexp = Sf + mean(covfold) + sum(cov_E) - (nfolds-1)/nfolds*min(SF2,Sf + mean(covfold) + sum(cov_E))
  list <- list("stage" = s, "lpav" = mean(lpavfold), "sd"= sqrt(max(varexp,0)))
  return(list)
}

#' Function to compute the Local Population Average Prescriptive Effect under SMART for cross-validated DTR.
#'
#' Given a A,D,Y data set, this function estimates the Population Average Prescriptive of an DTR derived from the same data by CV up to stage s.
#'
#' @param Acv A matrix of the actual treatment. Each row represents a unit-level randomized sequential assignments.
#' @param Dhatcv A list of DTR matrix derived from each trail of cross-validation.
#' @param Ycv A vector of the final outcome variable of interest for each sample.
#' @param indcv A vector of integers (between 1 and number of folds inclusive) indicating which testing set does each sample belong to.
#' @param s The stage up to which DTR should be evaluated..
#' @param centered If \code{TRUE}, the outcome variables would be centered before processing. This minimizes
#' the variance of the estimator. Default is \code{FALSE}.
#' @return A list that contains the following items: \item{stage}{The stage up to which DTR has been evaluated}，
#' \item{lpape}{The estimated local pape of the estimated DTR.},
#' \item{sd}{The standard deviation of the estimated local pape.}
#' @examples
#' Acv <- matrix(rep(0,1000*2),nrow = 1000,ncol= 2)
#' for(i in 1:2){
#'  n1=round(1000/2)
#'  n0=1000-n1
#'  ind=sample(1:1000,size=n1)
#'  Acv[ind,i] = 1}
#' Dhatcv <- list()
#' for (i in 1:5) {
#' D <- matrix(0, nrow = 1000, ncol = 2)
#' num_ones <- 500 + (i - 1)
#' D[1:num_ones, 2] <- 1
#' Dhatcv [[i]] <- D
#' }
#' Ycv =  rnorm(1000,3,1.5)
#' indcv <- sample(rep(1:5, each = 200))
#' lpape_est <- cvLocal_PAPE_dtr(Acv,Dhatcv,Ycv,indcv,2,centered = FALSE)
#' lpape_est$lpape
#' lpape_est$sd
#' @export cvLocal_PAPE_dtr
cvLocal_PAPE_dtr <- function(Acv,Dhatcv,Ycv,indcv,s,centered=FALSE){
  t = ncol(Acv)
  n = nrow(Acv)
  nfolds = max(indcv)
  if (centered) {
    Ycv = Ycv - mean(Ycv)
  }
  M_Dhat = lapply(Dhatcv, function(x){
    M_d = matrix(nrow = nrow(x),ncol = 2^t)
    for (i in 1:nrow(x)) {
      M_d[i,] = dtrMatrix(x[i,],s=t)
    }
    return(M_d)
  })
  M_A = matrix(nrow = n,ncol = 2^t)
  for(i in 1:n) {
    # M[i,] = trtMatrix(A[i,],D[i,]) #follow dtr
    M_A[i,] = trtMatrix(Acv[i,],M_Dhat[[1]][i,],s=1) #follow randomization
  }
  n_A = colSums(M_A)
  p = n_A/n
  lpapefold = c()
  sdfold = c()
  Sf = 0
  covfold = c()
  lsavloo1est_ztr= matrix(nrow=nfolds,ncol = 2^t)
  sateloo2est_ztr= matrix(nrow=nfolds,ncol = 2^t)
  for (k in 1:nfolds) {
    a = Acv[indcv==k,]
    y = Ycv[indcv==k]
    dhat = Dhatcv[[k]][indcv==k,]
    output = Local_PAPE_dtr(a,dhat,y,s)
    lpapefold = c(lpapefold, output$lpape)
    st = output$St/nfolds
    st = ifelse(is.na(st),0,st)
    Sf = Sf + st

    lpape =  Local_PAPE_dtr(Acv,Dhatcv[[k]],Ycv,s)
    covfold = c(covfold, lpape$cov)
    #sdfold = c(sdfold,output$sd^2)
    #m_A = matrix(nrow = m,ncol = 2^t) #follow randomization
    m_D =M_Dhat[[k]]#dtr metrix
    m_S = matrix(nrow = n,ncol = 2^t) #local randomization
    m_D_s <- matrix(nrow = n,ncol = 2^(s-1)) #construct matrix for loo esitmation p(d_T|d_s)
    for(i in 1:n) {
      # m_A[i,] = trtMatrix(a[i,],dhat[i,],s=1) # [(1,1),(1,0),(0,1),(0,0)] A1A2
      # m_D[i,] = dtrMatrix(dhat[i,],s=t) #D1D1
      m_S[i,] = trtMatrix(Acv[i,],Dhatcv[[k]][i,],s=s) #D1A2
      m_D_s[i,] = dtrMatrix(Dhatcv[[k]][i,],s=(s-1)) # denominator move one stage ahead P(a2|a1)
    }
    m_pd_s = t(apply(m_D_s,1,function(x){
      M <- c()
      for(i in 1:length(x)){M <- c(M,rep(x[i],2^(t-s+1)))}
      return(M)}))
    mdSum <- colSums(m_D)
    #maSum <- colSums(M_A)
    mpdSum <- colSums(m_pd_s)
    m_L = matrix(nrow = n,ncol = 2^t)
    m_L2 = matrix(nrow = n,ncol = 2^t)
    p_loo_1 <- matrix(nrow = n,ncol = 2^t) #p(a2|a1)
    p_loo_2 <- matrix(nrow = n,ncol = 2^t) # p(a1a2)
    for(i in 1:n){
      p_loo_1[i,] <- (mdSum - m_D[i,])/(mpdSum-m_pd_s[i,])
      p_loo_1[i,] <- ifelse(is.na(p_loo_1[i,]),0,p_loo_1[i,])
      m_L[i,] = m_S[i,]*Ycv[i]*p_loo_1[i,]/p
      p_loo_2[i,] <- (mdSum - m_D[i,])/(n-1)
      m_L2[i,] = M_A[i,]*Ycv[i]*p_loo_2[i,]/p
    }
    #covariance
    loo1est = colMeans(p_loo_1) #p(a2|a1)
    loo2est = colMeans(p_loo_2) #p(\bat a2)
    lsav = colSums(sweep(m_S, MARGIN=2,1/p, `*`) * Ycv)/n #E[d(a1)Y(a1a2)]
    sate = colSums(sweep(M_A, MARGIN=2,1/p, `*`) * Ycv)/n  #E[Y(a1a2)]
    lsavloo1est_ztr[k,] = lsav*loo1est
    sateloo2est_ztr[k,] = sate*loo2est
  }
  cov_E = matrix(nrow = 2^t,ncol = 2^t)
  for(i in 1:2^t){
    for(j in 1:2^t){
      if(i ==j){
        cov_E[i,j]=var(lsavloo1est_ztr[,i])-cov(lsavloo1est_ztr[,i],sateloo2est_ztr[,j])-
          cov(sateloo2est_ztr[,i],lsavloo1est_ztr[,j])+
          var(sateloo2est_ztr[,i])
      }
      else{
        cov_E[i,j]= cov(lsavloo1est_ztr[,i],lsavloo1est_ztr[,j])-
          cov(lsavloo1est_ztr[,i],sateloo2est_ztr[,j])-
          cov(sateloo2est_ztr[,i],lsavloo1est_ztr[,j])+
          cov(sateloo2est_ztr[,i],sateloo2est_ztr[,j])
      }
    }
  }
  SF2 = var(lpapefold)
  varexp = Sf + mean(covfold)+sum(cov_E) -
    (nfolds-1)/nfolds*min(SF2,Sf + mean(covfold)+sum(cov_E))
  list <- list("stage" = s, "lpape" = mean(lpapefold), "sd"= sqrt(max(varexp,0)))
  return(list)
}
