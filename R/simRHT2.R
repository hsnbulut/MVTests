#' @title Monte Carlo Simulation to obtain d and q constants for RHT2 function
#'
#' @description
#' Monte Carlo Simulation to obtain d and q constants for RHT2 function
#'
#' @details
#' \code{simRHT2} function computes d and q constants to construct an approximate 
#' F distribution of robust Hotelling T^2 statistic in high dimensional data. 
#' These constants are used in \code{RHT2} function.
#' For more detailed information, you can see the study by Bulut (2021).
#' 
#' @importFrom rrcov CovControlMrcd CovMrcd
#' @importFrom mvtnorm rmvnorm
#' @param n the sample size
#' @param p the number of variables
#' @param nrep the number of iteration. The default value is 500.
#' 
#' @export
#'
#' @return a list with 2 elements:
#' \item{q}{The q value}
#' \item{d}{The d value}
#' @references Bulut, H (2021). A robust Hotelling test statistic for one sample case in highdimensional data,
#' Communication in Statistics: Theory and Methods.
#' @author Hasan BULUT <hasan.bulut@omu.edu.tr>




simRHT2<-function(n,p,nrep=500){
  mu<-rep(0,p)
  sigma<-diag(p)
  control1<-CovControlMrcd(maxcsteps=100)  
  
  robHT2<-function(mu,sigma,n){
    return(n*t(mu)%*%solve(sigma)%*%mu)
  }
  
  T2r<-rep(0,nrep)
  for (i in 1:nrep){
    data<-rmvnorm(n,mean = mu,sigma = sigma) 
    mrcd<-CovMrcd(x=data,control=control1)
    T2r[i]<-robHT2(mu = as.matrix(mrcd@center),
                   sigma = as.matrix(mrcd@cov),
                   n=n)
  }
  
  muT2<-mean(T2r)
  varT2<-var(T2r)
  q<-((varT2/muT2^2)*p/2-1)^(-1)*(p+2)+4
  d<-muT2*(q-2)/q
  
  return(list(q=q,d=d))
}

