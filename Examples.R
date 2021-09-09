rm(list=ls())
require(rstudioapi)
require(data.table)
require(plot3D)
require(plot3Drgl)
require(mvnfast)
require(mgcv)
require(Matrix)
require(scoringRules)
require(xtable)

require(roxygen2)
require(devtools)

setwd(dirname(getActiveDocumentContext()$path))
source("R/CovarianceFunctions.R")
source("R/Optimisation.R")

# remove.packages("gac")
# # Update package documentation
# document(pkg = ".")
# # Install from local repository
# install(".")
# # Load Package
# require(gac)


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

## Entropy Score
cov_entropy <- function(R_true,R_est){
  RR <- solve(R_true) %*% R_est
  sum(RR[diag(ncol(RR))==1]) - log(det(RR)) - ncol(RR)
}



## Notes for the future ####

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
       zlim=c(0,2),theta = -10,phi = 10)
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


# Evaluate fit to empirical covariance for given parameters (used later to
# optimise parameters)
gac_obj(c(1,1,1),R,Cov_R_sim,cov_func=PowExp,loss="WLS")

# Find parameters to minimise gac_obj
Fit1 <- optim(par=c(1,1,1),
              gac_obj,
              method = "L-BFGS-B",
              lower=c(0,0,0),
              upper = c(Inf,Inf,2),
              R=R,Emp_Cov = Cov_R_sim,
              cov_func=Spherical)

# Calculate estimated covariance matrix from result of optim
Cov_R_fit <- PowExp(R,params =  Fit1$par)

# Plot covariance functions (simulated data, true, fit)
plot(c(R),c(Cov_R_sim),pch=16,col=rgb(0,1,0,alpha = .1))
points(c(R),c(Cov_R),pch=16)
points(c(R),c(Cov_R_fit),pch=16,col=2)

# Remove everything apart from functions
rm(list = setdiff(ls(), lsf.str()))





## Dynamic Isotropic Example ####


r <- seq(0,1,length.out=6)
R <- as.matrix(dist(r))

N <- 5000
x <- runif(N)
# theta_fn <- function(x){x^2+1}
# theta_fn <- function(x){1/(1+exp(-10*(x-0.5)))+0.5}
theta_fn <- function(x){sin(2*pi*x)+2}

# True Cov
image(t(PowExp(R,params = list(sigma=1,theta=min(theta_fn(x)),gamm=1))))
image(t(PowExp(R,params = list(sigma=1,theta=max(theta_fn(x)),gamm=1))))

# Empirical from simulation
Z <- matrix(NA,N,ncol(R))
for(i in 1:N){
  Z[i,] <- mvnfast::rmvn(n = 1,
                         mu=rep(0,ncol(R)),sigma = PowExp(R,params = list(sigma=1,theta=theta_fn(x[i]),gamma=1)))
}

image(t(cov(Z)))

image(t(cov(Z[x<0.5,])))
image(t(cov(Z[x>0.5,])))


PowExp_1param <- function(params,...){PowExp(params = list(sigma=1,theta=params[[1]],gamma=1),...)}

iso_ex_static_fit <- gac(R = R,
                         Emp_Cov = cov(Z),
                         cov_func = PowExp_1param,
                         param_eqns = list(~1),
                         loss="WLS")
iso_ex_static_fit$gac_coef


iso_ex_linear_fit <- gac(R = R,
                         X = list(x1=x),
                         Emp_Cov = NULL,
                         data = Z,
                         cov_func = PowExp_1param,
                         param_eqns = list(~x1),
                         loss="WLS")


iso_ex_smooth_fit <- gac(R = R,
                         X = list(x1=x),
                         Emp_Cov = NULL,
                         data = Z,
                         cov_func = PowExp_1param,
                         param_eqns = list(~s(x1,bs="cr",k=5)),
                         loss="WLS",smoothness_param = 0)

# Plot Basis Functions
matplot(x=iso_ex_smooth_fit$modelling_table$x1,
        y=iso_ex_smooth_fit$gam_prefits[[1]]$X,type="l")

# Plot Smooth Fit
matplot(x=iso_ex_smooth_fit$modelling_table$x1,
        y=iso_ex_smooth_fit$gam_prefits[[1]]$X * 
          matrix(rep(iso_ex_smooth_fit$gac_coef[[1]],each=nrow(iso_ex_smooth_fit$gam_prefits[[1]]$X)),
                 ncol=ncol(iso_ex_smooth_fit$gam_prefits[[1]]$X)),
        type="l")
lines(iso_ex_smooth_fit$modelling_table$x1,
      iso_ex_smooth_fit$gam_prefits[[1]]$X %*% iso_ex_smooth_fit$gac_coef[[1]],col=3)


## Plot different estimates of theta_fn:
plot(x[order(x)],theta_fn(x)[order(x)],type="l",ylim=c(1,3))
lines(iso_ex_linear_fit$modelling_table$x1,
      iso_ex_linear_fit$gam_prefits[[1]]$X %*% iso_ex_linear_fit$gac_coef[[1]],col=2)
lines(iso_ex_smooth_fit$modelling_table$x1,
      iso_ex_smooth_fit$gam_prefits[[1]]$X %*% iso_ex_smooth_fit$gac_coef[[1]],col=3)


## Predict function for fits
theta_fn_linear_pred <- approxfun(x = iso_ex_linear_fit$modelling_table$x1,
                                  y = iso_ex_linear_fit$gam_prefits[[1]]$X %*% iso_ex_linear_fit$gac_coef[[1]],
                                  rule = 2)

theta_fn_smooth_pred <- approxfun(x = iso_ex_smooth_fit$modelling_table$x1,
                                  y = iso_ex_smooth_fit$gam_prefits[[1]]$X %*% iso_ex_smooth_fit$gac_coef[[1]],
                                  rule = 2)


## Evaluation of fits


# Realisations
N_oos <- 5000
x_oos <- runif(N_oos)

# Empirical from simulation
Z_oos <- matrix(NA,N_oos,ncol(R))
for(i in 1:N_oos){
  Z_oos[i,] <- mvnfast::rmvn(n = 1,
                             mu=rep(0,ncol(R)),sigma = PowExp(R,params = list(sigma=1,theta=theta_fn(x_oos[i]),gamma=1)))
}

cov_mats <- list(Name=c("True","Static Empirical","GAC-Linear","GAC-CR"))
scores <- data.table(index=rep(1:nrow(Z_oos),length(cov_mats$Name)),
                     Name=rep(cov_mats$Name,each=nrow(Z_oos)))
# es=NULL,vs=NULL,`Log Score`=NULL)

for(cov_i in 1:length(cov_mats$Name)){
  for(i in 1:nrow(Z_oos)){
    
    if(cov_mats$Name[cov_i]=="True"){
      cov_temp <- PowExp(R,params = list(sigma=1,theta=theta_fn(x_oos[i]),gamma=1))  
    }else if(cov_mats$Name[cov_i]=="Static Empirical"){
      cov_temp <- cov(Z)
    }else if(cov_mats$Name[cov_i]=="GAC-Linear"){
      cov_temp <- PowExp(R,params = list(sigma=1,theta=theta_fn_linear_pred(x_oos[i]),gamma=1))  
    }else if(cov_mats$Name[cov_i]=="GAC-CR"){
      cov_temp <- PowExp(R,params = list(sigma=1,theta=theta_fn_smooth_pred(x_oos[i]),gamma=1))  
    }else{
      stop("cov_i not recognised.")
    }
    
    
    # Draw sample trajectories
    traj <- mvnfast::rmvn(n = 1000,
                          mu=rep(0,ncol(R)),sigma = cov_temp)
    
    # es and vs
    scores[index==i & Name==cov_mats$Name[[cov_i]], es := es_sample(y=Z_oos[i,],dat = t(traj))]
    scores[index==i & Name==cov_mats$Name[[cov_i]], vs_0_5 := vs_sample_quick(y=Z_oos[i,],dat = t(traj), p = 0.5)]
    scores[index==i & Name==cov_mats$Name[[cov_i]], vs_1 := vs_sample_quick(y=Z_oos[i,],dat = t(traj), p = 1)]
    
    # log score
    scores[index==i & Name==cov_mats$Name[[cov_i]], ls := -log(mvnfast::dmvn(X = Z_oos[i,],mu = rep(0,ncol(R)),sigma = cov_temp))]
    
    # Entropy
    scores[index==i & Name==cov_mats$Name[[cov_i]], entropy := 
             cov_entropy(R_true = PowExp(R,params = list(sigma=1,theta=theta_fn(x_oos[i]),gamma=1)),
                         R_est = cov_temp)]
    
  }
}
rm(cov_temp)

print(setorder(scores[,mean(es),by="Name"],V1))
print(setorder(scores[,mean(vs_0_5),by="Name"],V1))
print(setorder(scores[,mean(vs_1),by="Name"],V1))
print(setorder(scores[,mean(ls),by="Name"],V1))
print(setorder(scores[,mean(entropy),by="Name"],V1))


print(xtable(scores[,.(`Energy`=mean(es),
                       `Log`=mean(ls),
                       `VS-0.5`=mean(vs_0_5),
                       `VS-1`=mean(vs_1),
                       Entropy = round(mean(entropy),3)),
                    by="Name"][order(-Log),],digits = 3,
             caption = c("Results of simulation experiment for example \\ref{sec:iso_dyn_cov}: Isotropic dynamic covariance."),
             label = c("tab:iso_dyn_cov")),
      caption.placement = "top",
      include.rownames=F)






# Remove everything apart from functions
rm(list = setdiff(ls(), lsf.str()))







## Anisotropic Example ####


r <- seq(0,1,length.out=24)
R <- as.matrix(dist(r))

# Z <- r %*% t(r) # NB: Cov is no longer a function of separation only... 
# Cov_R <- as.matrix(nearPD(PowExp(R,params = list(sigma=1,theta=2+1/(.1+sqrt(Z)),gamm=1)))$mat)

N <- length(r)
Z <- matrix(rep(1:N,N) + rep(1:N,each=N) - 1, N, N)/(2*N-1)
Cov_R <- as.matrix(nearPD(PowExp(R,params = list(sigma=1,theta=1+1/Z,gamm=1)))$mat)


image(t(Z))
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
                       loss="LS")
test_static_fit$gac_coef

test_fit <- gac(R = R,
                X = list(x1=Z),
                Emp_Cov = Cov_R_sim,
                cov_func = PowExp,
                param_eqns = list(~1,
                                  ~s(x1,bs="bs"),
                                  ~1),
                loss="WLSf",
                smoothness_param = 1e-4)

image(t(test_fit$Cov_Est))


modelling_table$y_est <- c(test_fit$Cov_Est)
plot(x=modelling_table$r,
     y=modelling_table$y_est,
     col=rgb(1-modelling_table$x1,0,modelling_table$x1))





image(Cov_R_sim-nearPD(Cov_R_sim)$mat)
image(test_fit$Cov_Est-nearPD(test_fit$Cov_Est)$mat)

## Evaluate with log score, variogram score and energy score


# Realisations
actuals <- mvnfast::rmvn(n = 1000,mu=rep(0,ncol(Cov_R)),sigma = Cov_R)

cov_mats <- list(Name=c("True","Empirical","Constant","GAC"),
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
    traj <- mvnfast::rmvn(n = 1000,mu=rep(0,ncol(Cov_R)),sigma = cov_mats$mat[[cov_i]])
    
    # es and vs
    scores[index==i & Name==cov_mats$Name[[cov_i]], es := es_sample(y=actuals[i,],dat = t(traj))]
    scores[index==i & Name==cov_mats$Name[[cov_i]], vs_0_5 := vs_sample_quick(y=actuals[i,],dat = t(traj), p = 0.5)]
    scores[index==i & Name==cov_mats$Name[[cov_i]], vs_1 := vs_sample_quick(y=actuals[i,],dat = t(traj), p = 1)]
    
    # log score
    scores[index==i & Name==cov_mats$Name[[cov_i]], ls := -log(mvnfast::dmvn(X = actuals[i,],mu = rep(0,ncol(actuals)),sigma = cov_mats$mat[[cov_i]]))]
    
    # Entropy
    scores[index==i & Name==cov_mats$Name[[cov_i]], entropy := 
             cov_entropy(R_true = Cov_R,R_est = cov_mats$mat[[cov_i]])]
    
  }
}


print(setorder(scores[,mean(es),by="Name"],V1))
print(setorder(scores[,mean(vs_0_5),by="Name"],V1))
print(setorder(scores[,mean(vs_1),by="Name"],V1))
print(setorder(scores[,mean(ls),by="Name"],V1))


print(xtable(scores[,.(`Energy`=mean(es),
                       `Log`=mean(ls),
                       `VS-0.5`=mean(vs_0_5),
                       `VS-1`=mean(vs_1),
                       Entropy = round(mean(entropy),3)),
                    by="Name"][order(-Log),],digits = 3,
             caption = c("Results of simulation experiment for example \\ref{sec:aniso_cov}: Isotropic dynamic covariance."),
             label = c("tab:aniso_cov")),
      caption.placement = "top",
      include.rownames=F)



# Remove everything apart from functions
rm(list = setdiff(ls(), lsf.str()))






## Wind and SOlar ####
## Notes on data from Ciaran:

# Large RDAs of data, and RMD to load



load("data/windsolar_fc.rda")

zone_dat[,g_val:=qnorm(u_val)]
zone_dat[g_val%in%c(-Inf,Inf),g_val:=NA]

gobs <- dcast(zone_dat,
              kfold+issueTime~id+lead_time,
              value.var = c("g_val"),
              drop = TRUE)
setorder(gobs,issueTime)
gobs[,1:10]

WindSolar_Cov <- cor(gobs[,-c(1,2)],use = "pairwise.complete.obs")
# image(WindSolar_Cov)

col6 <- colorRampPalette(c("blue","cyan","yellow","red"))
lattice::levelplot(WindSolar_Cov,xlab="node id", ylab="node id",
                   col.regions=col6(600), cuts=100, at=seq(-0,1,0.01),
                   scales=list(x=list(at = seq(0,48,12),rot=45),y=list(rot=45),tck=0.3,cex=0.1))



## Wind - Temporal ####
# issue time dependency
gobs <- dcast(zone_dat[id=="wind_scotland",],
              kfold+issueTime~id+lead_time,
              value.var = c("g_val"),
              drop = TRUE)
setorder(gobs,issueTime)
gobs[,1:10]


WindScot_Cov <- cor(gobs[hour(issueTime)==12,][,-c(1,2)],use = "pairwise.complete.obs")
WindScot_Cov <- cor(gobs[hour(issueTime)==0,][,-c(1,2)],use = "pairwise.complete.obs")

col6 <- colorRampPalette(c("blue","cyan","yellow","red"))
lattice::levelplot(WindScot_Cov,,xlab = "lead time [hours]",ylab = "lead time [hours]",
                   col.regions=col6(600), cuts=100, at=seq(-.1,1,0.01),
                   scales=list(x=list(at = seq(1,97,12),lab=seq(0,96,12)/2,rot=45),
                               y=list(at = seq(1,97,12),lab=seq(0,96,12)/2,rot=45),tck=0.3,cex=0.75))

r <- seq(0,48,by=0.5)
R <- as.matrix(dist(r))


# surf3D(matrix(r,length(r),length(r),byrow = F),
#        matrix(r,length(r),length(r),byrow = T),
#        WindScot_Cov,
#        colvar = WindScot_Cov, colkey = F, facets = F,bty="f",
#        xlab="Lead-time",ylab="Lead-time",zlab="Covariance",
#        zlim=c(0,1),theta = -10,phi = 10)
# plotrgl()

## For correlation, need to fix sigma = 1, so covert to correlation function
PowExp_cor <- function(params,...){
  PowExp(params = list(sigma=1,theta=params[[1]],gamma=params[[2]]),...)
}


## Static fit
WindScot_static_fit <- gac(R = R,
                           X = list(),
                           Emp_Cov = WindScot_Cov,
                           cov_func = PowExp_cor,
                           param_eqns = list(#~1,
                             ~1,
                             ~1),
                           loss="WLS")


lattice::levelplot(WindScot_static_fit$Cov_Est,xlab="node id", ylab="node id",
                   col.regions=col6(600), cuts=100, at=seq(-0.1,1,0.01),
                   scales=list(x=list(rot=45),y=list(rot=45),tck=0.3,cex=0.1))


## GAC fit

# "Dist along diagonal"
N <- length(r)
Z <- matrix(rep(1:N,N) + rep(1:N,each=N) - 1, N, N)/(2*N-1)

N <- matrix(rep(r,each=length(r)),length(r),length(r))
Z1 = cos(2*pi*N/12) + cos(2*pi*t(N)/12)
# image(Z1)

WindScot_gac_fit <- gac(R = R,
                        X = list(x1=Z,
                                 x2=Z1),
                        Emp_Cov = WindScot_Cov,
                        cov_func = PowExp_cor,
                        param_eqns = list(#~1,
                          ~s(x1,bs="bs",k=20),
                          # ~s(x1,bs="bs",k=20)+x2,
                          # ~s(x2,k=10),
                          ~1),
                        loss="WLS",smoothness_param = 0)


lattice::levelplot(WindScot_gac_fit$Cov_Est,xlab="node id", ylab="node id",
                   col.regions=col6(600), cuts=100, at=seq(0,1,0.01),
                   scales=list(x=list(rot=45),y=list(rot=45),tck=0.3,cex=0.1))




modelling_table <- data.frame(y=c(WindScot_Cov),
                              r=c(R),
                              x1=c(Z))

plot(x=modelling_table$r,
     y=modelling_table$y,
     col="grey",pch=".")


modelling_table$y_est <- c(WindScot_static_fit$Cov_Est)
points(x=modelling_table$r,
       y=modelling_table$y_est,
       col="red")


modelling_table$y_est <- c(WindScot_gac_fit$Cov_Est)
points(x=modelling_table$r,
       y=modelling_table$y_est,
       col="blue")






scatter3D(modelling_table$r,
          modelling_table$x1,
          modelling_table$y,
          col = "red",xlab="time diff",ylab="dist_along diag",zlab="covariance")

plotrgl()


plot1 <- plotly::plot_ly(x=modelling_table$r,
                         y=modelling_table$x1,
                         z=modelling_table$y, type="scatter3d", mode="markers", color = I("grey"),
                         marker = list(size = 1,opacity = 0.5),
                         name = "empirical")

modelling_table$y_est <- c(WindScot_static_fit$Cov_Est)
plot2 <- plotly::plot_ly(x=modelling_table$r,
                         y=modelling_table$x1,
                         z=modelling_table$y_est, type="scatter3d", mode="markers", color = I("red"),
                         marker = list(size = 1,opacity = 0.5),
                         name = "static")



modelling_table$y_est <- c(WindScot_gac_fit$Cov_Est)
plot3 <- plotly::plot_ly(x=modelling_table$r,
                         y=modelling_table$x1,
                         z=modelling_table$y_est, type="scatter3d", mode="markers", color = I("blue"),
                         marker = list(size = 1,opacity = 0.5),
                         name = "gac")

plot <- plotly::subplot(plot1,plot2,plot3,nrows = 1)
plotly::layout(plot, scene = list(xaxis = list(title = 'dt [hours]'), yaxis = list(title = 'dist along diag [-]'), zaxis = list(title = 'covariance [-]')))

invisible(gc())
















## Additive 2D-fourier term
PowExp_cor_spec <- function(params,...){
  PowExp(params = list(sigma=1,theta=params[[1]],gamma=params[[2]]),...) + params[[3]]*(2+Z1)
}

WindScot_gac_fit2 <- gac(R = R,
                         X = list(x1=c(Z),
                                  x2=c(Z1)),
                         Emp_Cov = WindScot_Cov,
                         cov_func = PowExp_cor_spec,
                         param_eqns = list(#~1,
                           ~s(x1,bs="bs",k=20),
                           ~1,
                           ~1),
                         loss="WLS")

lattice::levelplot(WindScot_gac_fit2$Cov_Est,xlab="node id", ylab="node id",
                   col.regions=col6(600), cuts=100, at=seq(-0.1,1,0.01),
                   scales=list(x=list(rot=45),y=list(rot=45),tck=0.3,cex=0.1))






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

