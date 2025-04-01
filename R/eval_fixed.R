#' Function to compute the Population Average Value under SMART.
#'
#' Given a A,D,Y dataset, this function estimates the Population Average Value of a given DTR.
#'
#' @param A A matrix of the actual treatment. Each row represents a unit-level randomized sequential assignments.
#' @param D A matrix of the unit-level dynamic treatment that would have been assigned by the
#' individualized dynamic treatment rule.
#' @param y A vector of the final outcome variable of interest for each sample.
#' @return A list that contains the following items: \item{sav}{sample average value for each possible treatment}
#' \item{pav}{The estimated population value of a DTR.}, \item{sd}{The standard deviation of the estimated population value.}
#' @examples
#' A <- matrix(rep(0,100*2),nrow = 100,ncol= 2)
#' for(i in 1:2){
#'  n1=round(100/2)
#'  n0=100-n1
#'  ind=sample(1:100,size=n1)
#'  A[ind,i] = 1}
#' D <- matrix(rep(0,100*2),nrow = 100,ncol= 2)
#' D[1:50,2] = 1
#' Y =  rnorm(100,3,1.5)
#' pav_est <- PAV_dtr(A,D,Y)
#' pav_est$pav
#' pav_est$sd
#' @export PAV_dtr

PAV_dtr <- function(A,D,y){
  t = ncol(A)
  n = nrow(A)
  M = matrix(nrow = n,ncol = 2^t)
  M_A = matrix(nrow = n,ncol = 2^t)
  M_D <- matrix(nrow = n,ncol = 2^t)
  for(i in 1:n) {
    M[i,] = trtMatrix(A[i,],D[i,])
    M_A[i,] = trtMatrix(A[i,],D[i,],s=1)
    M_D[i,] = dtrMatrix(D[i,],s=t)
  }
  p = colSums(M_A)/n
  sav = colSums(sweep(M, MARGIN=2,1/p, `*`)*y)/n
  Sf = numeric(2^t)
  for(i in 1:2^t){
    Sf[i] = var((M_D[,i]*y)[M_A[,i]==1])
  }
  Sf = ifelse(is.na(Sf),0,Sf)
  varexp = sum(Sf/(p*n))
  dtr_list <- list("sav"=sav,"pav" = sum(sav),"sd"=sqrt(max(varexp,0)))
  return(dtr_list)
}


#' Function to compute the Population Average Prescriptive Effect under SMART.
#'
#' Given a A,D,Y dataset, this function estimates Population Average Prescriptive Effect of a given DTR.
#'
#' @param A A matrix of the actual treatment. Each row represents a unit-level randomized sequential assignments.
#' @param D A matrix of the unit-level dynamic treatment that would have been assigned by the
#' individualized dynamic treatment rule.
#' @param y A vector of the final outcome variable of interest for each sample.
#' @return A list that contains the following items:
#'  \item{pape}{The estimated population average prescriptvie effect of a given DTR.},
#'  \item{sd}{The standard deviation of the estimated pape.},
#'  \item{st}{E[S_{t}^2]/n_t (for internal use).},
#'  \item{cov}{covariance term in the variance formula (for internal use).},
#'  \item{sate}{vector of sample average value of each possible assignment under actual treatment(for internal use).},
#'  \item{sav}{vector of sample average value of each possible assignment under dtr (for internal use).},
#'  \item{rpav}{vector of ramdomized average value of each possible assignment under actual treatment (for internal use).},
#'  \item{sape}{vector of sample pape of each possible assignment under dtr (for internal use).},
#' @examples
#' A <- matrix(rep(0,100*2),nrow = 100,ncol= 2)
#' for(i in 1:2){
#'  n1=round(100/2)
#'  n0=100-n1
#'  ind=sample(1:100,size=n1)
#'  A[ind,i] = 1}
#' D <- matrix(rep(0,100*2),nrow = 100,ncol= 2)
#' D[1:50,2] = 1
#' Y =  rnorm(100,3,1.5)
#' pape_est <- PAPE_dtr(A,D,Y)
#' pape_est$pape
#' pape_est$sd
#' @export PAPE_dtr

PAPE_dtr <- function(A,D,y){
  t = ncol(A)
  n = nrow(A)

  M = matrix(nrow = n,ncol = 2^t)
  M_A = matrix(nrow = n,ncol = 2^t)
  M_D <- matrix(nrow = n,ncol = 2^t)
  for(i in 1:n) {
    M[i,] = trtMatrix(A[i,],D[i,])
    M_A[i,] = trtMatrix(A[i,],D[i,],s=1)
    M_D[i,] = dtrMatrix(D[i,],s=t)
  }
  p = colSums(M_A)/n
  p_d = colSums(M_D)/n

  M_cov <- matrix(nrow =2^t, ncol = 2^t)
  for(i in 1:(2^t)){
    for(j in 1:(2^t)){
      M_cov[i,j] = sum((M_A[,i]*M_D[,j] - M_A[,i]*p_d[j])/p[i]*y)/(n-1)
    }
  }

  sate = colSums(sweep(M_A, MARGIN=2,1/p, `*`) * y)/n
  sav = colSums(sweep(M, MARGIN=2,1/p, `*`)*y)/n
  sape = diag(M_cov)
  rpav = sate*p_d

  S_t = apply((M-sweep(M_A,MARGIN=2,p_d, `*`))*y,2,function(x) var(x[which(x!= 0)]))
  S_t = ifelse(is.na(S_t),0,S_t)

  cov1 <-c()
  for(i in 1:2^t){
    cov1 <- c(cov1,(sape[i]^2+2*(n-1)*(2*p_d[i]-1)*sape[i]*sate[i]-n*p_d[i]*(1-p_d[i])*sate[i]^2)/n^2)
  }
  cov2 <- c()
  for(i in 1:(2^t-1)){
    for(j in (i+1):(2^t)){
      cov2 <- c(cov2,(M_cov[i,j]*M_cov[j,i]+n*p_d[i]*p_d[j]*sate[i]*sate[j]+(n-1)*(p_d[i]*M_cov[i,j]*sate[j]+p_d[j]*M_cov[j,i]*sate[i]+p_d[i]*sate[i]*sape[j]+p_d[j]*sate[j]*sape[i]))/n^2)
    }
  }
  cov = sum(cov1)+ 2*sum(cov2)
  St = sum(S_t/(p*n))
  varexp = (n/(n-1))^2*(St+cov)
  dtr_list <- list("pape" = sum(sape),"sd"=sqrt(max(varexp,0)),"st"= St,"cov"= cov,"sate" = sate,"sav"=sav,"rpav"=rpav,"sape" = sape)
  return(dtr_list)
}

#' Function to compute the Local Population Average Value under SMART.
#'
#' Given a A,D,Y dataset and stage s, this function estimates the Population Average Value of a given DTR up stage s.
#'
#' @param A A matrix of the actual treatment. Each row represents a unit-level randomized sequential assignments.
#' @param D A matrix of the unit-level dynamic treatment that would have been assigned by the
#' individualized dynamic treatment rule.
#' @param Y A vector of the final outcome variable of interest for each sample.
#' @param s The stage up to which DTR should be evaluated..
#' @return A list that contains the following items: \item{stage}{The stage up to which DTR has been evaluated}，
#' \item{lpav}{The estimated local population value of given DTR.}, \item{sd}{The standard deviation of the estimated local population value.}，
#' \item{cov}}{covariance term in the variance formula (for internal use).},\item{St}{E[S_{t}^2]/n_t (for internal use).}
#' @examples
#' A <- matrix(rep(0,100*2),nrow = 100,ncol= 2)
#' for(i in 1:2){
#'  n1=round(100/2)
#'  n0=100-n1
#'  ind=sample(1:100,size=n1)
#'  A[ind,i] = 1}
#' D <- matrix(rep(0,100*2),nrow = 100,ncol= 2)
#' D[1:50,2] = 1
#' Y =  rnorm(100,3,1.5)
#' lpav_est <- Local_PAV_dtr(A,D,Y,2)
#' lpav_est$lpav
#' lpav_est$sd
#' @export Local_PAV_dtr
Local_PAV_dtr <- function(A,D,Y,s){
  s = s
  t = ncol(A)
  n = nrow(A)
  m_A = matrix(nrow = n,ncol = 2^t)
  m_D <- matrix(nrow = n,ncol = 2^t)
  m_S = matrix(nrow = n,ncol = 2^t)
  m_D_s <- matrix(nrow = n,ncol = 2^(s-1))
  for(i in 1:n) {
    m_A[i,] = trtMatrix(A[i,],D[i,],s=1)
    m_D[i,] = dtrMatrix(D[i,],s=t)
    m_S[i,] = trtMatrix(A[i,],D[i,],s=s)
    m_D_s[i,] = dtrMatrix(D[i,],s=(s-1))
  }

  p = colSums(m_A)/n
  p_d = colSums(m_D)/n

  m_pd_s = t(apply(m_D_s,1,function(x){
    M <- c()
    for(i in 1:length(x)){
      M <- c(M,rep(x[i],2^(t-s+1)))
    }
    return(M)
  }))

  mdSum <- colSums(m_D)
  maSum <- colSums(m_A)
  mpdSum <- colSums(m_pd_s)

  m_L = matrix(nrow = n,ncol = 2^t)
  p_loo_1 <- matrix(nrow = n,ncol = 2^t)
  for(i in 1:n){
    p_loo_1[i,] <- (mdSum - m_D[i,])/(mpdSum-m_pd_s[i,])
    p_loo_1[i,] <- ifelse(is.na(p_loo_1[i,]),0,p_loo_1[i,])
    m_L[i,] = m_S[i,]*Y[i]*p_loo_1[i,]/p
  }
  lpav = sum(colSums(m_L)/n)

  S_t = numeric(2^t)
  m_St = sweep(m_L,MARGIN=2,p, `*`)
  for(i in 1:2^t){
    S_t[i] = var(m_St[,i][m_A[,i]==1])
  }
  S_t = ifelse(is.na(S_t),0,S_t)

  loo1est = colMeans(p_loo_1)
  lsav = colSums(sweep(m_S, MARGIN=2,1/p, `*`) * Y)/n
  lsavSum <- colSums(m_L)
  cov_t1 <- matrix(nrow =2^t, ncol = 2^t)
  for(k in 1:(2^t)){
    for(j in 1:(2^t)){
      t1 <- c()
      for(i in 1:n){
        t1[i] <- (lsavSum[j]*m_L[i,k] - m_L[i,k]*m_L[i,j])/(n-1)
      }
      cov_t1[k,j] <- mean(t1)
    }
  }

  cov <- matrix(nrow = 2^t,ncol = 2^t)
  for(k in 1:2^t){
    for(j in 1:2^t){
      cov[k,j] = (cov_t1[k,j]-loo1est[j]*loo1est[k]*lsav[j]*lsav[k])/n
    }
  }
  St = sum(S_t/(p*n))
  cov_sum = sum(cov,na.rm = T)
  varexp = St + cov_sum
  list <- list("stage" = s, "lpav"= lpav, "sd" = sqrt(max(varexp,0)),"cov" = cov_sum,"St" = St)
  return(list)
}

#' Function to compute the Local Population Average Prescriptive Effect under SMART.
#'
#' Given a A,D,Y dataset and stage s, this function estimates the Population Average Prescriptive Effect of a given DTR up stage s.
#'
#' @param A A matrix of the actual treatment. Each row represents a unit-level randomized sequential assignments.
#' @param D A matrix of the unit-level dynamic treatment that would have been assigned by the
#' individualized dynamic treatment rule.
#' @param Y A vector of the final outcome variable of interest for each sample.
#' @param s The stage up to which DTR should be evaluated..
#' @return A list that contains the following items: \item{stage}{The stage up to which DTR has been evaluated}，
#' \item{lpape}{The estimated local population average prescriptive effect of given DTR.},
#' \item{sd}{The standard deviation of the estimated local population average prescriptive effect .}，
#' \item{cov}}{covariance term in the variance formula (for internal use).},\item{St}{E[S_{t}^2]/n_t (for internal use).}
#' @examples
#' A <- matrix(rep(0,100*2),nrow = 100,ncol= 2)
#' for(i in 1:2){
#'  n1=round(100/2)
#'  n0=100-n1
#'  ind=sample(1:100,size=n1)
#'  A[ind,i] = 1}
#' D <- matrix(rep(0,100*2),nrow = 100,ncol= 2)
#' D[1:50,2] = 1
#' Y =  rnorm(100,3,1.5)
#' lpape_est <- Local_PAPE_dtr(A,D,Y,2)
#' lpape_est$lpape
#' lpape_est$sd
#' @export Local_PAPE_dtr

Local_PAPE_dtr <- function(A,D,Y,s){
  t = ncol(A)
  n = nrow(A)
  m_A = matrix(nrow = n,ncol = 2^t)
  m_D <- matrix(nrow = n,ncol = 2^t)
  m_S = matrix(nrow = n,ncol = 2^t)
  m_D_s <- matrix(nrow = n,ncol = 2^(s-1))
  for(i in 1:n) {
    m_A[i,] = trtMatrix(A[i,],D[i,],s=1)
    m_D[i,] = dtrMatrix(D[i,],s=t)
    m_S[i,] = trtMatrix(A[i,],D[i,],s=s)
    m_D_s[i,] = dtrMatrix(D[i,],s=(s-1))
  }
  p = colSums(m_A)/n
  p_d = colSums(m_D)/n

  m_pd_s = t(apply(m_D_s,1,function(x){
    M <- c()
    for(i in 1:length(x)){
      M <- c(M,rep(x[i],2^(t-s+1)))
    }
    return(M)
  }))

  mdSum <- colSums(m_D)
  maSum <- colSums(m_A)
  mpdSum <- colSums(m_pd_s)

  m_L = matrix(nrow = n,ncol = 2^t)
  m_L2 = matrix(nrow = n,ncol = 2^t)
  p_loo_1 <- matrix(nrow = n,ncol = 2^t)
  p_loo_2 <- matrix(nrow = n,ncol = 2^t)

  for(i in 1:n){
    p_loo_1[i,] <- (mdSum - m_D[i,])/(mpdSum-m_pd_s[i,])
    p_loo_1[i,] <- ifelse(is.na(p_loo_1[i,]),0,p_loo_1[i,])
    m_L[i,] = m_S[i,]*Y[i]*p_loo_1[i,]/p
    p_loo_2[i,] <- (mdSum - m_D[i,])/(n-1)
    m_L2[i,] = m_A[i,]*Y[i]*p_loo_2[i,]/p
  }
  lpav = sum(colSums(m_L)/n)
  rpav = sum(colSums(m_L2)/n)
  lpape = lpav - rpav

  S_t = apply(sweep((m_L-m_L2),MARGIN=2,p, `*`),2,function(x)var(x[which(x!= 0)]))
  S_t = ifelse(is.na(S_t),0,S_t)

  loo1est = colMeans(p_loo_1)
  lsav = colSums(sweep(m_S, MARGIN=2,1/p, `*`) * Y)/n
  sate = colSums(sweep(m_A, MARGIN=2,1/p, `*`) * Y)/n

  m_cov <- matrix(nrow =2^t, ncol = 2^t)
  for(i in 1:(2^t)){
    for(j in 1:(2^t)){
      m_cov[i,j] = sum(m_S[,i]*m_D[,j]/p[i]*Y)/n
    }
  }

  m_cov2 <- matrix(nrow =2^t, ncol = 2^t)
  for(i in 1:(2^t)){
    for(j in 1:(2^t)){
      m_cov2[i,j] = sum((m_A[,i]*m_D[,j]-m_A[,i]*p_d[j])/p[i]*Y)/(n-1)
    }
  }

  m_cov3 <- matrix(nrow =2^t, ncol = 2^t)
  for(i in 1:(2^t)){
    for(j in 1:(2^t)){
      m_cov3[i,j] = sum(m_A[,i]*m_D[,j]/p[i]*Y)/n
    }
  }

  lsavSum <- colSums(m_L)
  loo1Sum <- colSums(p_loo_1)
  ySum <- colSums(m_A*Y/p)
  crosSum <- colSums(m_D*(m_A*Y)/p)

  cov_t1 <- matrix(nrow =2^t, ncol = 2^t)
  cov_t3 <- matrix(nrow =2^t, ncol = 2^t)
  cov_t4 <- matrix(nrow =2^t, ncol = 2^t)

  for(k in 1:(2^t)){
    for(j in 1:(2^t)){
      t1 <- c()
      t3 <- c()
      t4 <- c()
      for(i in 1:n){
        t1[i] <- (lsavSum[j]*m_L[i,k] - m_L[i,k]*m_L[i,j])/(n-1)
        t3[i] <-  (loo1Sum[j]*(m_A[i,k]*Y[i]) - p_loo_1[i,j]*(m_A[i,k]*Y[i]))/(n-1)/p[k]
        t4[i] <- (mdSum[j]*ySum[j]*p_loo_1[i,k] - ySum[j]*m_D[i,j]*p_loo_1[i,k]- mdSum[j]*(m_A[i,j]*Y[i])/p[j]*p_loo_1[i,k] -
                    crosSum[j]*p_loo_1[i,k] + 2*m_D[i,j]*(m_A[i,j]*Y[i])/p[j]*p_loo_1[i,k])/((n-1)*(n-2))
      }
      cov_t1[k,j] <- mean(t1)
      cov_t3[k,j] <- mean(t3)
      cov_t4[k,j] <- mean(t4)
    }
  }

  cov <- matrix(nrow = 2^t,ncol = 2^t)
  for(k in 1:2^t){
    for(j in 1:2^t){
      cov[k,j] = cov_t1[k,j]-loo1est[j]*loo1est[k]*lsav[j]*lsav[k] -
        1/(n-1)*(cov_t3[k,j]*m_cov[k,j]+cov_t3[j,k]*m_cov[j,k]) -
        (n-2)/(n-1)*(lsav[j]*cov_t4[j,k]+lsav[k]*cov_t4[k,j]) +
        p_d[k]*loo1est[j]*lsav[j]*sate[k]+p_d[j]*loo1est[k]*lsav[k]*sate[j]+
        1/(n-1)^2*(#p_d[j]*sate[k]*m_cov2[k,j]+p_d[k]*sate[j]*m_cov2[j,k] + m_cov2[j,k]*m_cov2[k,j]+
          m_cov3[k,j]*m_cov3[j,k]-p_d[j]*p_d[k]*sate[j]*sate[k]+
            (n-2)*(p_d[j]*sate[k]*m_cov2[j,k]+p_d[k]*sate[j]*m_cov2[k,j]))
      if(j == k){
        cov[k,j] = cov[k,k]/n + (n-2)*sate[k]*sate[k]*p_d[k]*(1-p_d[k])/(n-1)^2
      }
      else cov[k,j] = cov[k,j]/n - (n-2)*sate[j]*sate[k]*p_d[j]*p_d[k]/(n-1)^2
    }
  }
  cov_sum = sum(cov,na.rm = T)
  varexp = sum(S_t/(p*n)) + cov_sum
  list <- list("stage" = s, "lpape" = lpape,"sd"= sqrt(max(varexp,0)),"cov"= cov_sum, "St" = sum(S_t/(p*n)))
  return(list)
}


