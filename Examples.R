rm(list=ls())
require(rstudioapi)
require(data.table)
require(plot3D)
require(plot3Drgl)
require(mvnfast)
require(splines)
setwd(dirname(getActiveDocumentContext()$path))

source("CovFunctions.R")

## Notes for the future ####

# Try nls()... e.g. formula = cov ~ cov_function(...)

# End goal: estimate and sample mvn. This script focuses on estimating
# parametric contrivance functions.

# Alternative, estimate (sparse) precision matrix e.g. "glasso" then
# sample using "sparseMVN::rmvn.sparse()"...


## Form single-variate symmetric matrix, sample and fit  ####

r <- seq(0,3,by=0.1)
R <- as.matrix(dist(r))
Cov_R <- PowExp(R,params = c(sqrt(2),1.5,0.8))
# Cov_R <- Spherical(R)

image(t(Cov_R))
surf3D(matrix(r,length(r),length(r),byrow = F),
       matrix(r,length(r),length(r),byrow = T),
       Cov_R,
       colvar = Cov_R, colkey = F, facets = F,bty="f",
       xlab="Lead-time",ylab="Lead-time",zlab="Covariance",
       zlim=c(0,1),theta = -10,phi = 10)
# plotrgl()

## Visualisations
# r <- seq(0,3,by=0.01)
# plot(r,PowExp(r),type="l",ylim = c(0,1.1))
# lines(r,Whittle_Matern(r),col=2)
# lines(r,Cauchy(r),col=3)
# lines(r,Spherical(r,theta = 2),col=4)

## Sample

data_sim <- mvnfast::rmvn(n = 2^11,mu=rep(0,ncol(Cov_R)),sigma = Cov_R)
Cov_R_sim <- cov(data_sim)
image(t(Cov_R_sim))
image(t(Cov_R_sim-Cov_R))

obj <- function(params,R,Emp_Cov,cov_func,loss="WLS",...){
  
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

obj(c(1,1,1),R,Cov_R_sim,cov_func=PowExp)

Fit1 <- optim(par=c(1,1,1),
              obj,
              method = "L-BFGS-B",
              lower=c(0,0,0),
              upper = c(Inf,Inf,2),
              R=R,Emp_Cov = Cov_R_sim,
              cov_func=Spherical)

Cov_R_fit <- PowExp(R,params =  Fit1$par)

plot(c(R),c(Cov_R_sim),pch=16,col=rgb(0,1,0,alpha = .1))
points(c(R),c(Cov_R),pch=16)
points(c(R),c(Cov_R_fit),pch=16,col=2)



# Remove everything apart from functions
rm(list = setdiff(ls(), lsf.str()))



## Example with changing parameter ####
require(mgcv)
require(Matrix)

r <- seq(0,1,length.out=100)
R <- as.matrix(dist(r))
Z <- r %*% t(r) # NB: Cov is no longer a function of separation only... 

image(t(Z))

# True Covariance
Cov_R <- as.matrix(nearPD(PowExp(R,params = list(sigma=sqrt(2),theta=2+1/(0.2+Z),gamm=1)))$mat)
image(t(Cov_R))

# Empirical from simulation
data_sim <- mvnfast::rmvn(n = 2^11,mu=rep(0,ncol(Cov_R)),sigma = Cov_R)
Cov_R_sim <- cov(data_sim)
image(t(Cov_R_sim))



modelling_table <- data.frame(y=c(Cov_R_sim),
                              r=c(R),
                              x1=c(Z))

plot(x=modelling_table$r,
     y=modelling_table$y,
     col=rgb(1-modelling_table$x1,0,modelling_table$x1))


## Function to fit model with gam-like expressions for cov function parameters
# Generalised Additive Covariance
#
# Smoothness parameter - ridge penalty
#
gac <- function(R, # Separation matrix
                X, # Some container of covariates... list with elements same size as R
                Emp_Cov, # Empirical covaraiance to fit model to alt: = cov(data)
                cov_func, # parametric covarianve function
                param_eqns, # list of gam expressions for parameters
                param_init=NULL, # parameter values to initialise optimisaiont
                loss="WLS" # Loss function for optimisation
){
  
  # ## For dev:
  # R <- R
  # X <- list(x1=Z)
  # Emp_Cov <- Cov_R_sim
  # cov_func <- PowExp
  # param_eqns <- list(~1,
  #                    ~bs(x1,df=5,intercept = F),
  #                    ~1)
  # ##
  
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
  obj_gac <- function(gac_coef,design_mat,R,Emp_Cov,cov_func,loss){
    
    # Calculate parametric covairance matrix from supplied parameters and equations/design matrix
    n_gac_coef <- unlist(lapply(design_mat,ncol))
    params <- list()
    for(i in 1:length(design_mat)){
      params[[i]] <-  matrix(design_mat[[i]] %*% gac_coef[
        sum(n_gac_coef[0:(i-1)])+1:n_gac_coef[i]],
        ncol = ncol(R))
    }
    obj(params = params,R = R,Emp_Cov = Emp_Cov,cov_func = cov_func,loss=loss,optim_bound=T)
    
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
  temp_test <- try(obj_gac(gac_coef = gac_coef_init,design_mat = design_mat,
                           R = R,Emp_Cov = Emp_Cov,cov_func = cov_func))
  if(class(temp_test)=="try.error"){stop("First evaluation of objective function (with param_init) failed.")}
  rm(temp_test)
  
  
  # Perform optimisation...
  Fit1 <- optim(par=gac_coef_init,
                obj_gac,
                ## Can impose box constraints, but not for anything but trivial models
                # method = "L-BFGS-B",
                # lower=cov_func(return_param_limits=T)$lower,
                # upper = cov_func(return_param_limits=T)$upper,
                method="BFGS",
                design_mat = design_mat,
                R=R,
                loss=loss,
                Emp_Cov = Emp_Cov,
                cov_func=cov_func)
  
  
  
  
  ## Return as "gac" class
  param_est <- list()
  gac_coef <- list()
  n_gac_coef <- unlist(lapply(design_mat,ncol))
  for(i in 1:length(design_mat)){
    gac_coef[[i]] <- Fit1$par[sum(n_gac_coef[0:(i-1)])+1:n_gac_coef[i]]
    
    param_est[[i]] <-  matrix(design_mat[[i]] %*% gac_coef[[i]],ncol = ncol(R))
  }
  
  output <- list(call=match.call(),
                 Cov_Est = cov_func(r=R,params = param_est),
                 param_est=param_est,
                 gac_coef=gac_coef)
  
  class(output) <- "gac"
  return(output)
  
}



## Form two-variable symmetric matrices, sample and fit ####

r1 <- seq(0,3,length.out = 5)
r2 <- seq(0,3,length.out = 12)

R1 <- as.matrix(dist(r1))
R2 <- as.matrix(dist(r2))

## Create separable two-domain covariance function
Cov_Func2 <- Seperable_func(f1=PowExp,f2=Cauchy,n_params1 = 3)

Cov2_R <- Cov_Func2(rep(1,7),R1,R2)

image(t(Cov2_R))

sepfun <- Cov_Func2(rep(1,7),r1,t(r2))
surf3D(matrix(r1,length(r1),length(r2),byrow = F),
       matrix(r2,length(r1),length(r2),byrow = T),
       sepfun,
       colvar = sepfun, colkey = F, facets = F,bty="f",
       xlab="Spatial Separation",ylab="Temporal Separation",zlab="Covariance",
       zlim=c(0,1),theta = -10,phi = 10)
# plotrgl()


## Sample
data_sim2 <- mvnfast::rmvn(n = 2^12,mu=rep(0,ncol(Cov2_R)),sigma = Cov2_R)
Cov2_R_sim <- cov(data_sim2)
image(t(Cov2_R_sim))
image(t(Cov2_R_sim-Cov2_R))

obj2 <- function(params,R1,R2,Emp_Cov,cov_func,loss="WLS"){
  
  Par_cov <- cov_func(params,R1,R2)
  
  if(loss=="WLS"){
    # Weighted least squares of covariance only (VARIANCE excluded!)
    mean(((Emp_Cov-Par_cov)/abs(1-cov2cor(Par_cov)))[(R1 %x% R2) !=0]^2)
    
  }else if(loss=="LS"){
    # Least squares (including variance)
    mean((Emp_Cov-Par_cov))^2
    
  }else{
    stop("Loss not recognised.")
  }
  
}

obj2(rep(1,7),R1,R2,Cov2_R_sim,cov_func=Cov_Func2)

Fit2 <- optim(par=rep(1,7),
              obj2,
              method = "L-BFGS-B",
              lower=c(0,0,0,0,0,0,0),
              upper = c(Inf,Inf,2,Inf,2,Inf,Inf),
              R1=R1,R2=R2,Emp_Cov = Cov2_R_sim,
              cov_func=Cov_Func2)

Cov2_R_fit <- Cov_Func2(params =  Fit2$par,R1,R2)

image(t(Cov2_R_fit))

scatter3D(R1 %x% matrix(1,nrow(R2),nrow(R2)),
          matrix(1,nrow(R1),nrow(R1)) %x% R2,
          Cov2_R_sim,
          col = rgb(0,1,0,alpha = 0.2),cex=0.5,
          xlab="Spatial Separation",ylab="Temporal Separation",zlab="Covariance",
          zlim=c(0,1),theta = -10,phi = 10)

scatter3D(R1 %x% matrix(1,nrow(R2),nrow(R2)),
          matrix(1,nrow(R1),nrow(R1)) %x% R2,
          Cov2_R,
          col = 2,add = T)

plotrgl()

# Save as html!!! :)
# saveWidget(widget = rglwidget(x = scene3d()),file="test.html")

# Remove everything apart from functions
rm(list = setdiff(ls(), lsf.str()))




## Fit to GSP Group Data ####
load("data/GSPG_9_2_GAM-Grid_cGPD.rda")
load("data/GSPG_Centroids.rda")

## Get Time-difference matrix
tt <- unique(as.numeric(gsub(pattern = "[A-P]",replacement = "",names(UDATA)[-1])))
R2 <- abs(matrix(tt,length(tt),length(tt),byrow = T) - matrix(tt,length(tt),length(tt),byrow = F))
rm(tt)

## Get spatial separation matrix
# h <- data.table(GSPG=gsub(pattern = "[0-9\\.]",replacement = "",names(UDATA)[-1]))
setkey(Centro,"GSPG")
# h <- merge(h,Centro,by="GSPG",all.x=T)
R1 <- as.matrix(dist(Centro[,.(lon,lat)],diag=T))
R1 <- R1/median(R1)


## Empirical Cov
Cov_e <- cov(qnorm(as.matrix(UDATA[,-1])),use = "pairwise.complete.obs")

image(t(Cov_e))

### Single GSPG/Temporal ####
G=6
surf3D(matrix(R2[1,],nrow(R2),nrow(R2)),
       matrix(R2[1,],nrow(R2),nrow(R2),byrow = T),
       Cov_e[1:nrow(R2)+nrow(R2)*(G-1),1:nrow(R2)+nrow(R2)*(G-1)],
       colvar = Cov_e[1:nrow(R2)+nrow(R2)*(G-1),1:nrow(R2)+nrow(R2)*(G-1)],
       colkey = F, facets = F,bty="f",
       xlab="Temporal Separation",ylab="Temporal Separation",zlab="Covariance",
       theta = -10,phi = 10)
plotrgl()


plot(c(R2),c(Cov_e[1:nrow(R2)+nrow(R2)*(G-1),1:nrow(R2)+nrow(R2)*(G-1)]))
image(t(Cov_e[1:nrow(R2)+nrow(R2)*(G-1),1:nrow(R2)+nrow(R2)*(G-1)]))


### Spatio-temporal cov function ####
# what a mess!
scatter3D(R1 %x% matrix(1,nrow(R2),nrow(R2)),
          matrix(1,nrow(R1),nrow(R1)) %x% R2,
          Cov_e,
          col = rgb(0,1,0,alpha = 0.2),cex=0.5,
          xlab="Spatial Separation",ylab="Temporal Separation",zlab="Covariance",
          zlim=c(0,1),theta = -10,phi = 10)

plotrgl()

