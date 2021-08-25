#' Objective functions for covariance functions
#'
#' Function to evaluate a selected objective function for a given covariance
#' function and empirical covairance function.
#' 
#' @author Jethro Browell, \email{jethro.browell@@glasgow.ac.uk}
#' @param params A list of parameters to be passed to specified covairance function \code{cov_func}
#' @param R Separation matrix or vector
#' @param Emp_Cov Empirical covariance matrix against which covariance function is evaluated
#' @param cov_func A covariance function
#' @param loss Chosen loss function. Options are:
#' \itemize{
#'  \item{WLS}{Weighted Leased Squares - weigthing by correlation.}
#'  \item{LS}{Least Squares}
#' }
#' @param ... Additional arguments passed to \code{cov_func}
#' @details Function that returns the value of the chosen loss function for given
#' empirical covariance matrix and covairance function for use in numerical methods.
#' @keywords Covariance Function
#' @export
gac_obj <- function(params,R,Emp_Cov,cov_func,loss="WLS",...){
  
  # Calculate parametric covairance matrix from supplied parameters
  Par_cov <- cov_func(r = R,params,...)
  
  if(loss=="WLS"){
    # Weighted least squares of covariance only (VARIANCE excluded!)
    mean(((Emp_Cov-Par_cov)/abs(1-cov2cor(Par_cov)))[R!=0]^2)
    
  }else if(loss=="LS"){
    # Least squares (including variance)
    mean((Emp_Cov-Par_cov))^2
    
  }else{
    stop("Loss not recognised.")
  }
  
}


#
# Need to add smoothness parameter - ridge penalty
#

#' Generalised Additive Covariance Functions
#'
#' Function to fit a covairance function with generalised additive-type models
#' for parameters.
#' 
#' @author Jethro Browell, \email{jethro.browell@@glasgow.ac.uk}
#' @param params R Separation matrix
#' @param X Some container of covariates, e.g. a list with elements same size as \code{R}
#' @param Emp_CovEmpirical covaraiance to fit model to
#' @param cov_func A parametric covarianve function
#' @param param_eqns A list of equations for each parameter of \code{cov_func}
#' @param param_init Parameter values to initialise optimisaion
#' @param loss Chosen loss function. Options are:
#' \itemize{
#'  \item{WLS}{Weighted Leased Squares - weigthing by correlation.}
#'  \item{LS}{Least Squares}
#' }
#' @details Fits models for generalised additive covariance functions. Work in progress!
#' @return Returns an object of class \code{gac} for which methods will be writtern...
#' @keywords Covariance Function
#' @export
gac <- function(R,
                X,
                Emp_Cov,
                cov_func,
                param_eqns,
                param_init=NULL,
                loss="WLS"){
  
  ## Create modeling table of expanded basis from param_equations
  modelling_table <- data.frame(y=c(Emp_Cov),
                                r=c(R))
  design_mat <- list()
  for(i in names(X)){
    modelling_table[[i]] <- c(X[[i]])
  }
  for(i in 1:length(param_eqns)){
    
    design_mat[[i]] <- model.matrix(param_eqns[[i]],data = modelling_table)
    
  }
  
  ## Create objective function for model parameters ~ need some penalty on smoothness?
  internal_gac_obj <- function(gac_coef,design_mat,R,Emp_Cov,cov_func){
    
    # Calculate parametric covairance matrix from supplied parameters and equations/design matrix
    n_gac_coef <- unlist(lapply(design_mat,ncol))
    params <- list()
    for(i in 1:length(design_mat)){
      params[[i]] <-  matrix(design_mat[[i]] %*% gac_coef[
        sum(n_gac_coef[0:(i-1)])+1:n_gac_coef[i]],
        ncol = ncol(R))
    }
    
    gac_obj(params = params,R = R,Emp_Cov = Emp_Cov,cov_func = cov_func,loss = loss, optim_bound=T)
    
  }
  
  
  ## Estimate parameters
  
  # Prep initial values
  if(is.null(param_init)){
    param_init <- rep(1,length(param_eqns))
  }
  gac_coef_init <- c()
  for(i in 1:length(design_mat)){
    gac_coef_init <- c(gac_coef_init,param_init[i],rep(0,ncol(design_mat[[i]])-1))
  }
  
  # Check first evaluation...
  temp_test <- try(internal_gac_obj(gac_coef = gac_coef_init,design_mat = design_mat,
                                    R = R,Emp_Cov = Emp_Cov,cov_func = cov_func))
  if(class(temp_test)=="try.error"){stop("First evaluation of objective function (with param_init) failed.")}
  rm(temp_test)
  
  
  # Perform optimisation...
  Fit1 <- optim(par = gac_coef_init,
                fn = internal_gac_obj,
                design_mat = design_mat,
                R = R,
                Emp_Cov = Emp_Cov,
                cov_func = cov_func,
                # loss = loss,
                method = "BFGS")
  
  
  
  
  ## Return as "gac" class
  param_est <- list()
  gac_coef <- list()
  n_gac_coef <- unlist(lapply(design_mat,ncol))
  for(i in 1:length(design_mat)){
    gac_coef[[i]] <- Fit1$par[sum(n_gac_coef[0:(i-1)])+1:n_gac_coef[i]]
    
    param_est[[i]] <-  matrix(design_mat[[i]] %*% gac_coef[[i]],ncol = ncol(R))
  }
  
  output <- list(call=match.call(),
                 R=R,X=X,
                 Cov_Est = cov_func(r=R,params = param_est),
                 param_est=param_est,
                 gac_coef=gac_coef)
  
  class(output) <- "gac"
  return(output)
  
}
