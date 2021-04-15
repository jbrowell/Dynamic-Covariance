rm(list=ls())
require(rstudioapi)
require(data.table)
require(plot3D)
require(plot3Drgl)
require(mvnfast)
setwd(dirname(getActiveDocumentContext()$path))


## Notes for the future ####

# End goal: estimate and sample mvn. This script focuses on estiamting
# parametric covaraince functions.

# Alternative, estiamte (sparse) precision matrix e.g. "glasso" then
# sample using "sparseMVN::rmvn.sparse()"...


## Functions ####


## Powered Exponential
PowExp <- function(r,params=list(sigma=1,theta=1,gamma=1)){
  
  if(params[[1]] < 0){stop("sigma<0")}
  if(params[[2]] <=0){stop("theta<=0")}
  if(params[[3]] <= 0 | params[[3]] > 2){stop("gamma <= 0 | gamma > 2")}
  
  return((params[[1]]^2)*exp(-(params[[2]]*r)^params[[3]]))
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


## Seperable
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



## Form single-variate symetric matrix, sample and fit  ####

r <- seq(0,3,by=0.1)
R <- as.matrix(dist(r))
# Cov_R <- PowExp(R,params = c(sqrt(2),1.5,0.8))
Cov_R <- Spherical(R)

image(t(R))
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

data_sim <- mvnfast::rmvn(n = 2^10,mu=rep(0,ncol(Cov_R)),sigma = Cov_R)
Cov_R_sim <- cov(data_sim)
image(t(Cov_R_sim))
image(t(Cov_R_sim-Cov_R))

obj <- function(params,R,Emp_Cov,cov_func,loss="WLS"){
  
  Par_cov <- cov_func(r = R,params)
  
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

obj(c(1,1),R,Cov_R_sim,cov_func=Spherical)

Fit1 <- optim(par=c(1,1),
              obj,
              method = "L-BFGS-B",
              lower=c(0,0,0),
              upper = c(Inf,Inf,2),
              R=R,Emp_Cov = Cov_R_sim,
              cov_func=Spherical)

Cov_R_fit <- Spherical(R,params =  Fit1$par)

plot(c(R),c(Cov_R_sim),pch=16,col=rgb(0,1,0,alpha = .1))
points(c(R),c(Cov_R),pch=16)
points(c(R),c(Cov_R_fit),pch=16,col=2)











## Form two-variable symetric matrics, sample and fit ####
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




