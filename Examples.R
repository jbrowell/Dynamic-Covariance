rm(list=ls())
require(rstudioapi)
require(data.table)
require(plot3D)
require(plot3Drgl)
require(mvnfast)
require(splines)

require(roxygen2)
require(devtools)

setwd(dirname(getActiveDocumentContext()$path))


# Update package documentation
document(pkg = ".")
# Install from local repository
install(".")
# Load Package
require(gac)


##  faster VS score, half the time on sR... ####
vs_sample_quick <- function (y, dat, w = NULL, p = 0.5) {
  
  d <- length(y)
  
  out <- 0
  for (i in seq_len(d)) {
    
    for (j in seq_len(i)){
      
      vdat <- mean(abs(dat[i, ] - dat[j, ])^p)
      vy <- abs(y[i] - y[j])^p
      
      if (is.null(w)) {
        out <- out + (vy - vdat)^2
      } else {
        out <- out + w[i, j] * (vy - vdat)^2  
        
      }
      
    }
  }
  
  return(2*out)
}





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


gac_obj(c(1,1,1),R,Cov_R_sim,cov_func=PowExp,loss="WLS")

Fit1 <- optim(par=c(1,1,1),
              gac_obj,
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
require(Matrix)

r <- seq(0,1,length.out=24)
R <- as.matrix(dist(r))
Z <- r %*% t(r) # NB: Cov is no longer a function of separation only... 

image(t(Z))

# True Covariance
Cov_R <- as.matrix(nearPD(PowExp(R,params = list(sigma=1,theta=2+1/(.1+sqrt(Z)),gamm=1)))$mat)
image(t(Cov_R))

# Empirical from simulation
data_sim <- mvnfast::rmvn(n = 720,
                          mu=rep(0,ncol(Cov_R)),sigma = Cov_R)
Cov_R_sim <- cov(data_sim)
image(t(Cov_R_sim))



modelling_table <- data.frame(y=c(Cov_R_sim),
                              r=c(R),
                              x1=c(Z))

plot(x=modelling_table$r,
     y=modelling_table$y,
     col=rgb(1-modelling_table$x1,0,modelling_table$x1))



test_static_fit <- gac(R = R,
                       X = list(x1=Z),
                       Emp_Cov = Cov_R_sim,
                       cov_func = PowExp,
                       param_eqns = list(~1,
                                         ~1,
                                         ~1),
                       loss="WLS")



test_fit <- gac(R = R,
                X = list(x1=Z),
                Emp_Cov = Cov_R_sim,
                cov_func = PowExp,
                param_eqns = list(~1,
                                  ~bs(x1,df=5,intercept = F),
                                  ~1),
                loss="WLS")



image(t(test_fit$Cov_Est))


modelling_table$y_est <- c(test_fit$Cov_Est)
plot(x=modelling_table$r,
     y=modelling_table$y_est,
     col=rgb(1-modelling_table$x1,0,modelling_table$x1))




## Evaluate with log score, variogram score and energy score
require(scoringRules)

# Realisations
actuals <- mvnfast::rmvn(n = 100,mu=rep(0,ncol(Cov_R)),sigma = Cov_R)

cov_mats <- list(Name=c("True","Empirical","Static","GAC"),
                 mat=list(Cov_R,
                          nearPD(Cov_R_sim)$mat,
                          nearPD(test_static_fit$Cov_Est)$mat,
                          nearPD(test_fit$Cov_Est)$mat))
scores <- data.table(index=rep(1:nrow(actuals),length(cov_mats$Name)),
                     Name=rep(cov_mats$Name,each=nrow(actuals)))
# es=NULL,vs=NULL,`Log Score`=NULL)

for(cov_i in 1:length(cov_mats$Name)){
  for(i in 1:nrow(actuals)){
    
    # Draw sample trajectories
    traj <- mvnfast::rmvn(n = 250,mu=rep(0,ncol(Cov_R)),sigma = cov_mats$mat[[cov_i]])
    
    # es and vs
    scores[index==i & Name==cov_mats$Name[[cov_i]], es := es_sample(y=actuals[i,],dat = t(traj))]
    scores[index==i & Name==cov_mats$Name[[cov_i]], vs_0_5 := vs_sample_quick(y=actuals[i,],dat = t(traj), p = 0.5)]
    scores[index==i & Name==cov_mats$Name[[cov_i]], vs_1 := vs_sample_quick(y=actuals[i,],dat = t(traj), p = 1)]
    
    # log score
    scores[index==i & Name==cov_mats$Name[[cov_i]], ls := -log(mvnfast::dmvn(X = actuals[i,],mu = rep(0,ncol(actuals)),sigma = cov_mats$mat[[cov_i]]))]
    
  }
}


scores[,mean(es),by="Name"]
scores[,mean(vs_0_5),by="Name"]
scores[,mean(vs_1),by="Name"]
scores[,mean(ls),by="Name"]






## Notes on data from Ciaran:

# Large RDAs of data, and RMD to load




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

