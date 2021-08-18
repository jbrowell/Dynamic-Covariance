## Functions ####
#
# Parametric covariance functions, utilities for creating
# new covariance functions, and tools for working with them.

## To do:
#
# - Add tests on imputs and return helpful errors when incoherent
# - Add facility for nugget effect
# - Add non-seperable functions. See Gneiting2007a for generic form
# - Add generic periodic fucntions for selection of kernels
# - Add facility for Lagrangian covariance functions
#
# - Allow parameters to vary with covariate...
## Add return_param_limits option? (Not used in optim in the end...)
## Return Inf for bad params for use in numeric optimisation 

## Powered Exponential
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



## Whittle-Matern
Whittle_Matern <- function(r,params=list(sigma=1,theta=1,nu=1)){
  
  if(params[[1]] < 0){stop("sigma<0")}
  if(params[[2]] <= 0){stop("theta<=0")}
  if(params[[3]] <= 0){stop("nu<=0")}
  
  return((params[[1]]^2) * (2^(1-params[[3]]))/gamma(params[[3]]) * (params[[2]]*r)*besselK(params[[2]]*r,params[[3]]))
}


## Cauchy
Cauchy <- function(r,params=list(sigma=1,gamma=1,nu=1,theta=1)){
  
  if(params[[1]] < 0){stop("sigma<0")}
  if(params[[2]] <= 0 | gamma > 2){stop("gamma <= 0 | gamma > 2")}
  if(params[[3]] <= 0){stop("nu<=0")}
  if(params[[4]] <= 0){stop("theta<=0")}
  
  
  return(params[[1]]^2 * (1+(params[[4]]*r)^params[[2]])^(-params[[3]]))
}



## Spherical
Spherical <- function(r,params=list(sigma=1,theta=1)){
  
  if(params[[1]] < 0){stop("sigma<0")}
  if(params[[2]] <= 0){stop("theta<=0")}
  
  suppressWarnings(
    output <- params[[1]]^2 * (1 - (2/pi)*((r/params[[2]])*sqrt(1-(r/params[[2]])^2) + asin(r/params[[2]])))
  )
  output[r>params[[2]]] <- 0
  
  return(output)
  
}


## A seperable covariance function, the product of f1 and f2.
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
