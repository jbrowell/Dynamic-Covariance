## Functions ####
#
# Parametric covariance functions, utilities for creating
# new covariance functions, and tools for working with them.

## To do:
#
# - Return Inf for bad params for use in numeric optimisation 
# - Add tests on imputs and return helpful errors when incoherent
# - Add facility for nugget effect
# - Add non-seperable functions. See Gneiting2007a for generic form
# - Add generic periodic fucntions for selection of kernels
# - Add facility for Lagrangian covariance functions
#


#' Powered Exponential Covariance Function
#'
#' Functional form of the Powered Exponential Covariance Function
#' 
#' @author Jethro Browell, \email{jethro.browell@@glasgow.ac.uk}
#' @param r Vector or matrix of separation distances
#' @param params A list of parameters with default \code{list(sigma=1,theta=1,gamma=1)}.
#' Parameters may be supplied as vectors of length equal to \code{length(r)}.
#' @param return_param_limits Boolean. If \code{TRUE} parameter limits are returned.
#' @param optim_bound Boolean. If \code{TRUE} parameters out side of limits do not result
#' in errors and instead a value of \code{Inf} is returned for these elements. For use in numerical
#' methods/optimisation.
#' @details Function that returns the value of the Powered Exponential Covariance function.
#' @return A vector of matrix the same size as \code{r} containing corresponding values of
#' the Powered Exponential Covariance Function, unless \code{return_param_limits=T} in which
#' case parameter limits are returned instead.
#' @keywords Covariance Function
#' @export
PowExp <- function(r=NULL,params=list(sigma=1,theta=1,gamma=1),return_param_limits=F,optim_bound=F){
  
  if(return_param_limits){
    return(list(lower=c(0,.Machine$double.eps,.Machine$double.eps),
                upper=c(Inf,Inf,2)))
  }else{
    
    if(is.null(r)){stop("r missing with no default.")}
    
    if(!optim_bound){
      if(any(params[[1]] < 0)){stop("sigma<0")}
      if(any(params[[2]] <=0)){stop("theta<=0")}
      if(any(params[[3]] <= 0 | params[[3]] > 2)){stop("gamma <= 0 | gamma > 2")}
    }
    
    output <- (params[[1]]^2)*exp(-(params[[2]]*r)^params[[3]])
    
    if(optim_bound){
      output[params[[1]] < 0] <- Inf
      output[params[[2]] <=0] <- Inf
      output[params[[3]] <= 0 | params[[3]] > 2] <- Inf
    }
    
    return(output)
    
  }
}



#' Whittle-Matern Covariance Function
#'
#' Functional form of the Whittle-Matern covariance function.
#' 
#' @author Jethro Browell, \email{jethro.browell@@glasgow.ac.uk}
#' @param r Vector or matrix of separation distances
#' @param params A list of parameters with default \code{list(sigma=1,theta=1,nu=1)}.
#' Parameters may be supplied as vectors of length equal to \code{length(r)}.
#' @details Function that returns the value of the Whittle-Matern covariance function.
#' @return A vector of matrix the same size as \code{r} containing corresponding values of
#' the Whittle-Matern covariance function.
#' @keywords Covariance Function
#' @export
Whittle_Matern <- function(r,params=list(sigma=1,theta=1,nu=1)){
  
  if(params[[1]] < 0){stop("sigma<0")}
  if(params[[2]] <= 0){stop("theta<=0")}
  if(params[[3]] <= 0){stop("nu<=0")}
  
  return((params[[1]]^2) * (2^(1-params[[3]]))/gamma(params[[3]]) * (params[[2]]*r)*besselK(params[[2]]*r,params[[3]]))
}


#' Cauchy Covariance Function
#'
#' Functional form of the Cauchy covariance function.
#' 
#' @author Jethro Browell, \email{jethro.browell@@glasgow.ac.uk}
#' @param r Vector or matrix of separation distances
#' @param params A list of parameters with default \code{list(sigma=1,gamma=1,nu=1,theta=1)}.
#' Parameters may be supplied as vectors of length equal to \code{length(r)}.
#' @details Function that returns the value of the Cauchy covariance function.
#' @return A vector of matrix the same size as \code{r} containing corresponding values of
#' the Cauchy covariance function.
#' @keywords Covariance Function
#' @export
Cauchy <- function(r,params=list(sigma=1,gamma=1,nu=1,theta=1)){
  
  if(params[[1]] < 0){stop("sigma<0")}
  if(params[[2]] <= 0 | gamma > 2){stop("gamma <= 0 | gamma > 2")}
  if(params[[3]] <= 0){stop("nu<=0")}
  if(params[[4]] <= 0){stop("theta<=0")}
  
  
  return(params[[1]]^2 * (1+(params[[4]]*r)^params[[2]])^(-params[[3]]))
}



#' Spherical Covariance Function
#'
#' Functional form of the Spherical covariance function.
#' 
#' @author Jethro Browell, \email{jethro.browell@@glasgow.ac.uk}
#' @param r Vector or matrix of separation distances
#' @param params A list of parameters with default \code{list(sigma=1,theta=1)}.
#' Parameters may be supplied as vectors of length equal to \code{length(r)}.
#' @details Function that returns the value of the Spherical covariance function.
#' @return A vector of matrix the same size as \code{r} containing corresponding values of
#' the Spherical covariance function.
#' @keywords Covariance Function
#' @export
Spherical <- function(r,params=list(sigma=1,theta=1)){
  
  if(params[[1]] < 0){stop("sigma<0")}
  if(params[[2]] <= 0){stop("theta<=0")}
  
  suppressWarnings(
    output <- params[[1]]^2 * (1 - (2/pi)*((r/params[[2]])*sqrt(1-(r/params[[2]])^2) + asin(r/params[[2]])))
  )
  output[r>params[[2]]] <- 0
  
  return(output)
  
}


#' Construction of a seperable covariance function
#'
#' A seperable covariance function, the product of two individual covariance functions.
#' 
#' @author Jethro Browell, \email{jethro.browell@@glasgow.ac.uk}
#' @param params A vector of parameter values, a concationation of parameters for 
#' \code{f1} and \code{f2}.
#' @param r1 Vector of matrix of separation distances for first separation variable
#' @param r2 Vector of matrix of separation distances for second separation variable
#' @param f1 Covariance function for first separation variable
#' @param f2 Covariance function for second separation variable
#' @param n_params1 Number of parameters to be passed to \code{f1}. Remainder are passed
#' to \code{f2}.
#' @details Function that returns the value of the Whittle-Matern covariance function.
#' @return A vector of matrix the same size as \code{r} containing corresponding values of
#' the Whittle-Matern covariance function.
#' @keywords Covariance Function
#' @export
Seperable2 <- function(params, # vector of parameter values
                       r1,r2, # vectors or matrices of separation in domain 1 and 2
                       f1,f2, # covariance function for domains 1 and 2
                       n_params1 # number of parameters passed to f1. Remainder passed to f2
){
  
  # Kronecker product
  f1(r1,params[1:n_params1]) %x% f1(r2,params[(n_params1+1):length(params)])
  
}

Seperable_func <- function(f1,f2, # covariance function for domains 1 and 2
                           n_params1 # number of parameters passed to f1. Remainder passed to f2
){
  f <- function(params,r1,r2){
    Seperable2(params=params,r1=r1,r2=r2,f1=f1,f2=f2,n_params1 = n_params1)
  } 
}
