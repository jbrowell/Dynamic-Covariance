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
#'  \item{WLS}{Weighted Leased Squares - weigthing by cross-correlation}
#'  \item{WLSf}{Weighted Leased Squares - weigthing by all correlations}
#'  \item{LS}{Least Squares}
#' }
#' @param pen Additional penalty (numeric), e.g. on smoothness, added to loss
#' @param ... Additional arguments passed to \code{cov_func}
#' @details Function that returns the value of the chosen loss function for given
#' empirical covariance matrix and covairance function for use in numerical methods.
#' @keywords Covariance Function
#' @export
gac_obj <- function(params,R,Emp_Cov,cov_func,loss="WLS",pen=0,...){
  
  # Calculate parametric covairance matrix from supplied parameters
  Par_cov <- cov_func(r = R,params,...)
  
  if(loss=="WLS"){
    # Weighted least squares of covariance only (VARIANCE excluded!)
    mean(((Emp_Cov-Par_cov)/abs(1-cov2cor(Par_cov)))[R!=0]^2) + pen
    
  }else if(loss=="WLSf"){
    # Weighted least squares
    mean((abs(cov2cor(Par_cov))*(Emp_Cov-Par_cov))^2) + pen
    
  }else if(loss=="LS"){
    # Least squares (including variance)
    mean((Emp_Cov-Par_cov)^2) + pen
    
  }else if(loss=="AL"){
    # Absolute Loss
    mean(abs(Emp_Cov-Par_cov)) + pen
  }else if(loss=="WALf"){
    # Weighted Absolute
    mean(abs(cov2cor(Par_cov))*abs(Emp_Cov-Par_cov)) + pen
    
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
#' @param X Some container of covariates, e.g. a list with elements same size as \code{R} or
#' the same length as \code{data}
#' @param Emp_Cov Empirical covaraiance to fit model a static model to.
#' @param data Dataset to fit dynamic covariance model to. \code{ncol(data)} must equal
#' \code{nrow(R)}
#' @param cov_func A parametric covarianve function.
#' @param param_eqns A list of equations for each parameter of \code{cov_func}
#' @param param_init Parameter values to initialise optimisaion (constant/intercept parameters only)
#' @param loss Chosen loss function. Options are:
#' \itemize{
#'  \item{WLS}{Weighted Leased Squares - weigthing by correlation}
#'  \item{WLSf}{Weighted Leased Squares - weigthing by all correlations}
#'  \item{LS}{Least Squares}
#' }
#' @param smoothness_param Parameter to control smoothness, is the weight that multiples
#' squared second derivative of smooth splines.
#' @details Fits models for generalised additive covariance functions. Work in progress!
#' @return Returns an object of class \code{gac} for which methods will be written...
#' @keywords Covariance Function
#' @export
gac <- function(R,
                X=list(),
                Emp_Cov=NULL,
                data=NULL,
                cov_func,
                param_eqns,
                param_init=rep(1,length(param_eqns)),
                loss="WLS",
                smoothness_param=0){
  
  if(is.null(Emp_Cov) & is.null(data)){stop("One of Emp_Cov or data must be supplied.")}
  if(!is.null(Emp_Cov)){
    if(nrow(R)!=nrow(Emp_Cov) | ncol(R)!=ncol(Emp_Cov)){stop("Size of R and Emp_Cov don't match.")}
  }
  if(!is.null(data)){
    
    if(!is.null(Emp_Cov)){warning("data supplied so Emp_Cov will not used.")}
    Emp_Cov <- cov(data)
    
    if(nrow(R)!=ncol(R) | ncol(R)!=ncol(data)){stop("Size of R and Emp_Cov don't match.")}
  }
  
  ## Create modeling table of expanded basis from param_equations
  modelling_table <- data.frame(y=c(Emp_Cov),
                                r=c(R))
  
  design_mat <- list()
  pen_mat <- list()
  dyn_covariates <- c()
  gam_prefits <- list()
  for(i in names(X)){
    # if(all(dim(X[[i]]==dim(R)))){
    if(length(X[[i]])==length(R)){
      # Static component
      modelling_table[[i]] <- c(X[[i]])
    }else if(!is.null(data)){
      # Dynamic component
      if(length(X[[i]])!=nrow(data)){stop(paste("Length of covarite",i,"does not match data"))}
      modelling_table[[i]] <- quantile(X[[i]],probs=seq(0,1,length.out=length(R)))  #rep(NA,length(R))
      dyn_covariates <- c(dyn_covariates,i)
    }
  }
  
  for(i in 1:length(param_eqns)){
    
    gam_prefits[[i]] <- gam(update(param_eqns[[i]],"y~."),data = modelling_table,fit = F)
    
    colnames(gam_prefits[[i]]$X) <- gam_prefits[[i]]$term.names
    
    # Design matrix (for single row of "data")
    design_mat[[i]] <- gam_prefits[[i]]$X
    
    # Penalty matrix
    pen_mat[[i]] <- gam_prefits[[i]]$S
    
  }
  
  
  ## Construct design mat for dynamic components
  if(!is.null(data)){
    
    design_dyn_mat <- list()  
    for(i in 1:length(param_eqns)){
      
      design_dyn_mat[[i]] <- data.frame(id=1:nrow(data))
      
      # Cycle through covariates
      for(k in dyn_covariates){
        
        design_dyn_mat[[i]][,k] <- X[[k]]
        
        # Cycle through each basis of covariate k (STRING MATCHING NOT ROBUST!!!)
        for(k2 in colnames(design_mat[[i]])[
          grep(pattern = paste0("(",k,")"),colnames(design_mat[[i]]))] # linear terms ~ "k"
        ){
          
          design_dyn_mat[[i]][,k2] <- approx(x=modelling_table[[k]],
                                             y=gam_prefits[[i]]$X[,k2],
                                             xout = design_dyn_mat[[i]][,k])$y
        }
      }
    }
  }
  
  
  ## Create objective function for model parameters ~ need some penalty on smoothness?
  internal_gac_obj <- function(gac_coef){#,design_mat,R,Emp_Cov,cov_func){
    
    
    n_gac_coef <- unlist(lapply(design_mat,ncol))
    
    ## Smoothness penalty?
    penalty <- 0
    for(i in 1:length(design_mat)){
      if(length(pen_mat[[i]])==0){next}
      
      pen_dim <- cumsum(c(0,sapply(pen_mat[[i]],ncol)))
      
      temp_coef <- gac_coef[sum(n_gac_coef[0:(i-1)])+1:n_gac_coef[i]][colnames(design_mat[[i]])!="(Intercept)"]
      for(j in 1:length(pen_mat[[i]])){
        temp_coef2 <- temp_coef[1:ncol(pen_mat[[i]][[j]])+pen_dim[j]]
        penalty <- penalty + t(temp_coef2) %*% pen_mat[[i]][[j]] %*% temp_coef2
      }
    }
    
    ## Main loss
    if(is.null(data)){
      # Calculate parametric covairance matrix from supplied parameters and equations/design matrix
      
      params <- list()
      for(i in 1:length(design_mat)){
        params[[i]] <-  matrix(design_mat[[i]] %*% gac_coef[
          sum(n_gac_coef[0:(i-1)])+1:n_gac_coef[i]],
          ncol = ncol(R))
        
      }
      
      # Value to return:
      gac_obj(params = params,R = R,Emp_Cov = Emp_Cov,cov_func = cov_func,loss = loss,
              optim_bound=T,pen=penalty*smoothness_param)
      
    }else{
      
      obj_vlaue <- rep(NA,nrow(data))
      E_after_j <- matrix(0,nrow(R),ncol(R))
      for(j in 1:nrow(data)){
        params <- list()
        for(i in 1:length(design_mat)){
          
          ## Compute covariate/value of basis for current row [j] of data
          
          ### < This does not work as intended! >
          # for(k in dyn_covariates[dyn_covariates %in% colnames(design_mat[[i]])]){
          #   
          #   design_mat[[i]][,k] <- approx(x=modelling_table[[k]],
          #                                 y=gam_prefits[[i]]$X[,k],
          #                                 xout = rep(X[[k]][j],length(R)))$y
          #   
          # }
          ### < This does not work as intended! >
          
          # Cycle through covariates
          for(k in dyn_covariates){
            
            # Cycle through each basis of covariate k (STRING MATCHING NOT ROBUST!!!)
            for(k2 in colnames(design_mat[[i]])[
              c(grep(pattern = paste0("(",k,")"),colnames(design_mat[[i]])), # smooths ~ s()
                which(colnames(design_mat[[i]])==k))] # linear terms ~ "k"
            ){
              # ## SLOW! Construct this data matrix (design mat for dynamic) once outside of this function and look-up
              # design_mat[[i]][,k2] <- approx(x=modelling_table[[k]],
              #                                y=gam_prefits[[i]]$X[,k2],
              #                                xout = rep(X[[k]][j],length(R)))$y
              
              design_mat[[i]][,k2] <- rep(design_dyn_mat[[i]][j,k2],length(R))
              
            }
            
          }
          
          
          
          
          params[[i]] <-  matrix(design_mat[[i]] %*% gac_coef[
            sum(n_gac_coef[0:(i-1)])+1:n_gac_coef[i]],
            ncol = ncol(R))
          
        }
        
        Par_cov <- cov_func(r=R,params = params,optim_bound=T)
        
        if(loss=="LS"){
          E_after_j <- E_after_j + data[j,] %*% t(data[j,]) - Par_cov
          
        }else if(loss=="WLS"){
          E_after_j <- E_after_j + (data[j,] %*% t(data[j,]) - Par_cov)/abs(1-cov2cor(Par_cov))
          E_after_j[(R==0)] <- 0
        }else if(loss=="WLSf"){
          E_after_j <- E_after_j + (data[j,] %*% t(data[j,]) - Par_cov)*abs(cov2cor(Par_cov))
        }else{
          stop("Loss not recognised. Options: LS, WLS, WLSf")
        }
      }
      
      ## Least Squares
      obj_vlaue <- mean((E_after_j/j)^2)
      
      # Value to return:
      print(obj_vlaue)
      obj_vlaue + penalty*smoothness_param
    }
    
  }
  # debug(internal_gac_obj)
  
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
  temp_test <- try(internal_gac_obj(gac_coef = gac_coef_init))#,design_mat = design_mat,
  # R = R,Emp_Cov = Emp_Cov,cov_func = cov_func))
  if(class(temp_test)[1]=="try.error"){stop("First evaluation of objective function (with param_init) failed.")}
  rm(temp_test)
  
  
  # Perform optimisation...
  Fit1 <- optim(par = gac_coef_init,
                fn = internal_gac_obj,
                method = "BFGS"
                # method = "Nelder-Mead",hessian = F
  )
  
  
  
  
  ## Return as "gac" class
  param_est <- list()
  gac_coef <- list()
  n_gac_coef <- unlist(lapply(design_mat,ncol))
  penalty <- 0
  for(i in 1:length(design_mat)){
    gac_coef[[i]] <- Fit1$par[sum(n_gac_coef[0:(i-1)])+1:n_gac_coef[i]]
    
    param_est[[i]] <-  matrix(design_mat[[i]] %*% gac_coef[[i]],ncol = ncol(R))
    
    ## Smoothness penalty?
    if(length(pen_mat[[i]])==0){next}
    pen_dim <- cumsum(c(0,sapply(pen_mat[[i]],ncol)))
    temp_coef <- gac_coef[[i]][colnames(design_mat[[i]])!="(Intercept)"]
    for(j in 1:length(pen_mat[[i]])){
      temp_coef2 <- temp_coef[1:ncol(pen_mat[[i]][[j]])+pen_dim[j]]
      penalty <- penalty + t(temp_coef2) %*% pen_mat[[i]][[j]] %*% temp_coef2
    }
  }
  
  output <- list(call=match.call(),
                 R=R,X=X,
                 Cov_Est = cov_func(r=R,params = param_est),
                 modelling_table = modelling_table,
                 gam_prefits = gam_prefits,
                 param_est=param_est,
                 gac_coef=gac_coef,
                 loss_final = Fit1$value,
                 penalty = as.numeric(penalty))
  
  class(output) <- "gac"
  return(output)
  
}
