#' Indicator vector showing if Treatment align with DTR
#' @param A_i A vector of actual treatment.
#' @param D_i A vector of DTR.
#' @param s stage to start randomization (when not to consider D)
#' If s=NULL, if Assinment follow DTR till the end; if s=1 Indicator of actual treatment from the start.
#' if s = 2, Indicator of follwing DTR up to the first stage and then actual treatment from stage 2.
#' @return A matrix of treatment indicator.
#'
#' @export trtMatrix
#'
# Helper Functions

trtMatrix <- function(A_i,D_i, s = NULL){
  t = length(A_i)
  if(is.null(s)) s = t + 1
  M_i = c(1)
  for (i in 0:(t-1)){
    for(j in (2^i):(2^(i+1)-1)){
      if(i < (s-1)){
        M_i[2*j] = M_i[j]*A_i[i+1]*D_i[i+1]
        M_i[2*j+1] = M_i[j]*(1-A_i[i+1])*(1-D_i[i+1])
      }
      else{
        M_i[2*j] = M_i[j]*A_i[i+1]
        M_i[2*j+1] = M_i[j]*(1-A_i[i+1])
      }
    }
  }
  return(M_i[2^t:(2^(t+1)-1)])
}

#' indicator of DTR
#'
#' @param D_i A vector of DTR.
#' @param s time to start randomization
#'
#' @return DTR vector up to time s
#'
#' @export dtrMatrix
#'
dtrMatrix <- function(D_i,s){
  s = s
  t = length(D_i)
  M_i = c(1)
  for (i in 0:(t-1)){
    for(j in (2^i):(2^(i+1)-1)){
      M_i[2*j] = M_i[j]*D_i[i+1]
      M_i[2*j+1] = M_i[j]*(1-D_i[i+1])
    }
  }
  return(M_i[(2^s):(2^(s+1)-1)])
}
