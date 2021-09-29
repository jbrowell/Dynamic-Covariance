

## Preamble ####

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
require(boot)
require(progress)
require(latex2exp)


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



## From here you can skip ahead to any of the examples...





## Form single-variate symmetric matrix, sample and fit  ####

# This example is not included in the paper but outlines the procedure used
# in the "gac()" function.

## Define data generating process
r <- seq(0,3,by=0.1)
R <- as.matrix(dist(r))
Cov_R <- PowExp(R,params = c(sqrt(2),1.5,0.8))

## Visualise
image(t(Cov_R))
surf3D(matrix(r,length(r),length(r),byrow = F),
       matrix(r,length(r),length(r),byrow = T),
       Cov_R,
       colvar = Cov_R, colkey = F, facets = F,bty="f",
       xlab="Lead-time",ylab="Lead-time",zlab="Covariance",
       zlim=c(0,2),theta = -10,phi = 10)
# plotrgl()


## Sample training data from N(0,Cov_R) 

data_sim <- mvnfast::rmvn(n = 2^11,mu=rep(0,ncol(Cov_R)),sigma = Cov_R)
Cov_R_sim <- cov(data_sim)
image(t(Cov_R_sim))
image(t(Cov_R_sim-Cov_R))


## Single evaluation of objective function for model fit
gac_obj(c(1,1,1),R,Cov_R_sim,cov_func=PowExp,loss="WLS")

## Find parameters to minimise gac_obj
Fit1 <- optim(par=c(1,1,1),
              gac_obj,
              method = "L-BFGS-B",
              lower=c(0,0,0),
              upper = c(Inf,Inf,2),
              R=R,Emp_Cov = Cov_R_sim,
              cov_func=Spherical)

## Calculate estimated covariance matrix from result of optim
Cov_R_fit <- PowExp(R,params =  Fit1$par)

## Plot covariance functions (simulated data, true, fit)
plot(c(R),c(Cov_R_sim),pch=16,col=rgb(0,1,0,alpha = .1))
points(c(R),c(Cov_R),pch=16)
points(c(R),c(Cov_R_fit),pch=16,col=2)
legend("top",c("Data","Truth","Fit"),pch=16,col=c(rgb(0,1,0),1,2))

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

## True Cov
image(t(PowExp(R,params = list(sigma=1,theta=min(theta_fn(x)),gamm=1))))
image(t(PowExp(R,params = list(sigma=1,theta=max(theta_fn(x)),gamm=1))))

## Empirical cov from simulation
Z <- matrix(NA,N,ncol(R))
for(i in 1:N){
  Z[i,] <- mvnfast::rmvn(n = 1,
                         mu=rep(0,ncol(R)),sigma = PowExp(R,params = list(sigma=1,theta=theta_fn(x[i]),gamma=1)))
}

image(t(cov(Z)))
image(t(cov(Z[x<0.5,])))
image(t(cov(Z[x>0.5,])))

## Define covairance function (only single gac model to estimate)
PowExp_1param <- function(params,...){PowExp(params = list(sigma=1,theta=params[[1]],gamma=1),...)}

## Static model (constant theta)
iso_ex_static_fit <- gac(R = R,
                         Emp_Cov = cov(Z),
                         cov_func = PowExp_1param,
                         param_eqns = list(~1),
                         loss="WLS")
iso_ex_static_fit$gac_coef

## Simple linear fit for theta
iso_ex_linear_fit <- gac(R = R,
                         X = list(x1=x),
                         Emp_Cov = NULL,
                         data = Z,
                         cov_func = PowExp_1param,
                         param_eqns = list(~x1),
                         loss="WLS")

## GAC fit for theta
iso_ex_smooth_fit <- gac(R = R,
                         X = list(x1=x),
                         Emp_Cov = NULL,
                         data = Z,
                         cov_func = PowExp_1param,
                         param_eqns = list(~s(x1,bs="cr",k=5)),
                         loss="WLS",smoothness_param = 5e-5)



## Predict function for fits
theta_fn_linear_pred <- approxfun(x = iso_ex_linear_fit$modelling_table$x1,
                                  y = iso_ex_linear_fit$gam_prefits[[1]]$X %*% iso_ex_linear_fit$gac_coef[[1]],
                                  rule = 2)

theta_fn_smooth_pred <- approxfun(x = iso_ex_smooth_fit$modelling_table$x1,
                                  y = iso_ex_smooth_fit$gam_prefits[[1]]$X %*% iso_ex_smooth_fit$gac_coef[[1]],
                                  rule = 2)


## Evaluation of fits


## Realisations
N_oos <- N
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
  
  print(cov_mats$Name[cov_i])
  pb <- progress_bar$new(total = nrow(Z_oos))
  
  for(i in 1:nrow(Z_oos)){
    
    pb$tick()
    
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
    
    ## Quick scores to calculate
    
    # log score
    scores[index==i & Name==cov_mats$Name[[cov_i]], ls := -log(mvnfast::dmvn(X = Z_oos[i,],mu = rep(0,ncol(R)),sigma = cov_temp))]
    
    # Entropy
    scores[index==i & Name==cov_mats$Name[[cov_i]], entropy := 
             cov_entropy(R_true = PowExp(R,params = list(sigma=1,theta=theta_fn(x_oos[i]),gamma=1)),
                         R_est = cov_temp)]
    
    ## Slow scores to calculate (comment out to save time...)
    
    # Draw sample trajectories
    traj <- mvnfast::rmvn(n = 1000,
                          mu=rep(0,ncol(R)),sigma = cov_temp)
    
    # es and vs
    scores[index==i & Name==cov_mats$Name[[cov_i]], es := es_sample(y=Z_oos[i,],dat = t(traj))]
    scores[index==i & Name==cov_mats$Name[[cov_i]], vs_0_5 := vs_sample_quick(y=Z_oos[i,],dat = t(traj), p = 0.5)]
    scores[index==i & Name==cov_mats$Name[[cov_i]], vs_1 := vs_sample_quick(y=Z_oos[i,],dat = t(traj), p = 1)]
    
    
    
  }
}
rm(cov_temp)



## Save/load results used in paper
# save.image(file="PSCC22_plots/iso_example.Rda")
# load("PSCC22_plots/iso_example.Rda")
##

print(setorder(scores[,mean(es),by="Name"],V1))
print(setorder(scores[,mean(vs_0_5),by="Name"],V1))
print(setorder(scores[,mean(vs_1),by="Name"],V1))
print(setorder(scores[,mean(ls),by="Name"],V1))
print(setorder(scores[,mean(entropy),by="Name"],V1))


print(xtable(scores[,.(`Energy`=signif(mean(es),4),
                       `Log`=signif(mean(ls),4),
                       `VS-0.5`=signif(mean(vs_0_5),4),
                       `VS-1`=signif(mean(vs_1),4),
                       KL = signif(mean(entropy),4)),
                    by="Name"][order(-Log),],digits = 3,
             caption = c("Results of simulation experiment for example \\ref{sec:iso_dyn_cov}: Isotropic dynamic covariance.  Italics indicate that the corresponding skill score relative to the GAC-CR model are not significantly different from zero."),
             label = c("tab:iso_dyn_cov")),
      caption.placement = "top",table.placement="",
      include.rownames=F)



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

setEPS(); postscript("PSCC22_plots/iso_theta_fits.eps",width = 6,height = 5)

plot(x[order(x)],theta_fn(x)[order(x)],type="l", lwd=2,
     ylim=c(1,4),axes=F,
     xlab=TeX("$x_t$"),ylab=TeX("$\\hat{\\theta}$"))
axis(1);axis(2); grid()
lines(iso_ex_linear_fit$modelling_table$x1,
      iso_ex_linear_fit$gam_prefits[[1]]$X %*% iso_ex_linear_fit$gac_coef[[1]],col=2,lty=2,lwd=2)
lines(iso_ex_smooth_fit$modelling_table$x1,
      iso_ex_smooth_fit$gam_prefits[[1]]$X %*% iso_ex_smooth_fit$gac_coef[[1]],col=3,lty=3,lwd=2)
legend("top",c("True","GAC-Linear","GAC-CR"),col=1:3,lty=1:3,lwd=2)

dev.off()


## Bootstrap skill score
scores_boot <- data.table()
ref <- "GAC-CR"
for(m in scores[Name!=ref,unique(Name)]){
  for(s in names(scores)[-c(1:2)]){
    
    boot_data <- boot(data = cbind(scores[Name==m,get(s)],scores[Name==ref,get(s)]),
                      statistic = function(data, i){
                        100-100*mean(data[i,1])/mean(data[i,2])
                      },
                      R=1000)
    
    
    scores_boot <- rbind(scores_boot,
                         data.table(Name=m,Score=s,
                                    Mean=mean(boot_data$t),SD=sd(boot_data$t),
                                    ci_L=NA,
                                    ci_R=NA))
    
  }
}; rm(boot_data,s,m) 

scores_boot[,ci_L := Mean - SD*qt(0.025,999,lower.tail = F)]
scores_boot[,ci_R := Mean + SD*qt(0.025,999,lower.tail = F)]

scores_boot[ci_L>0,Diff_from_zero:="+"]
scores_boot[ci_R<0,Diff_from_zero:="-"]
scores_boot[ci_L<=0 & ci_R>=0,Diff_from_zero:="0"]


scores_boot





# Remove everything apart from functions
rm(list = setdiff(ls(), lsf.str()))







## Non-stationary Example ####


r <- seq(0,1,length.out=51)
R <- as.matrix(dist(r))
N <- 5000

# Z <- r %*% t(r) # NB: Cov is no longer a function of separation only... 
# Cov_R <- as.matrix(nearPD(PowExp(R,params = list(sigma=1,theta=2+1/(.1+sqrt(Z)),gamm=1)))$mat)

Z <- matrix(rep(r,length(r)) + rep(r,each=length(r)), length(r), length(r)) + 1
Cov_R <- as.matrix(nearPD(PowExp(R,params = list(sigma=1,theta=5/(Z),gamm=0.8)))$mat)

image(t(Z))
image(t(Cov_R))


nonstat_cor_plot <- function(data,filename=NULL,...){
  
  col6 <- colorRampPalette(c("blue","cyan","yellow","red"))
  
  h <- lattice::levelplot(data,
                          xlab = list(label=TeX("$l_1$"),cex=1.2),
                          ylab = list(label=TeX("$l_2$"),cex=1.2),
                          col.regions=col6(600), cuts=100, at=seq(-.1,1.05,0.01),
                          scales=list(x=list(at = seq(1,ncol(data),length.out=6),
                                             lab=0:5/5,rot=45),
                                      y=list(at = seq(1,ncol(data),length.out=6),
                                             lab=0:5/5,rot=45),
                                      tck=0.3,cex=1.2),...)
  
  if(!is.null(filename)){
    setEPS(); postscript(filename)
    print(h)
    dev.off()
  }
  
  print(h)
  
  
}

nonstat_cor_plot(Cov_R, filename = "PSCC22_plots/Non-stationary_example.eps")


# Empirical from simulation
data_sim <- mvnfast::rmvn(n = N,
                          mu=rep(0,ncol(Cov_R)),sigma = Cov_R)
Cov_R_sim <- cov(data_sim)

nonstat_cor_plot(Cov_R_sim, filename = "PSCC22_plots/Non-stationary_empirical.eps")



modelling_table <- data.frame(y=c(Cov_R_sim),
                              r=c(R),
                              x1=c(Z))

plot(x=modelling_table$r,
     y=modelling_table$y,
     col=rgb(1-modelling_table$x1/max(modelling_table$x1),0,modelling_table$x1/max(modelling_table$x1)))



test_static_fit <- gac(R = R,
                       X = list(x1=Z),
                       Emp_Cov = Cov_R_sim,
                       cov_func = PowExp,
                       param_eqns = list(~1,
                                         ~1,
                                         ~1),
                       loss="WLSf")
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

nonstat_cor_plot(test_fit$Cov_Est,filename = "PSCC22_plots/Non-stationary_gac.eps")


modelling_table$y_est <- c(test_fit$Cov_Est)
plot(x=modelling_table$r,
     y=modelling_table$y_est,
     col=rgb(1-modelling_table$x1/max(modelling_table$x1),0,modelling_table$x1/max(modelling_table$x1)))





image(Cov_R_sim-nearPD(Cov_R_sim)$mat)
image(test_fit$Cov_Est-nearPD(test_fit$Cov_Est)$mat)

## Evaluate with log score, variogram score and energy score


# Realisations
actuals <- mvnfast::rmvn(n = N,mu=rep(0,ncol(Cov_R)),sigma = Cov_R)

cov_mats <- list(Name=c("True","Empirical","Stationary","GAC"),
                 mat=list(Cov_R,
                          nearPD(Cov_R_sim)$mat,
                          nearPD(test_static_fit$Cov_Est)$mat,
                          nearPD(test_fit$Cov_Est)$mat))
scores <- data.table(index=rep(1:nrow(actuals),length(cov_mats$Name)),
                     Name=rep(cov_mats$Name,each=nrow(actuals)))
# es=NULL,vs=NULL,`Log Score`=NULL)

for(cov_i in 1:length(cov_mats$Name)){
  
  print(cov_mats$Name[cov_i])
  pb <- progress_bar$new(total = nrow(actuals))
  
  for(i in 1:nrow(actuals)){
    
    pb$tick()
    
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


## Save/load results used in paper
# save.image(file="PSCC22_plots/nonstat_example.Rda")
# load("PSCC22_plots/nonstat_example.Rda")
##


nonstat_cor_plot(Cov_R, filename = "PSCC22_plots/Non-stationary_example.eps")
nonstat_cor_plot(Cov_R_sim, filename = "PSCC22_plots/Non-stationary_empirical.eps")
nonstat_cor_plot(test_fit$Cov_Est,filename = "PSCC22_plots/Non-stationary_gac.eps")

print(setorder(scores[,mean(es),by="Name"],V1))
print(setorder(scores[,mean(vs_0_5),by="Name"],V1))
print(setorder(scores[,mean(vs_1),by="Name"],V1))
print(setorder(scores[,mean(ls),by="Name"],V1))
print(setorder(scores[,mean(entropy),by="Name"],V1))


print(xtable(scores[,.(`Energy`=signif(mean(es),4),
                       `Log`=signif(mean(ls),4),
                       `VS-0.5`=signif(mean(vs_0_5),4),
                       `VS-1`=signif(mean(vs_1),4),
                       KL = signif(mean(entropy),4)),
                    by="Name"][order(-Log),],digits = c(NA,NA,3,2,1,1,3),
             caption = c("Results of simulation experiment for example \\ref{sec:non-stationary_cov}: Isotropic dynamic covariance."),
             label = c("tab:non-stationary_cov")),
      caption.placement = "top",table.placement="",
      include.rownames=F)




## Bootstrap skill score
scores_boot <- data.table()
ref <- "GAC"
for(m in scores[Name!=ref,unique(Name)]){
  for(s in names(scores)[-c(1:2)]){
    
    boot_data <- boot(data = cbind(scores[Name==m,get(s)],scores[Name==ref,get(s)]),
                      statistic = function(data, i){
                        100-100*mean(data[i,1])/mean(data[i,2])
                      },
                      R=1000)
    
    
    scores_boot <- rbind(scores_boot,
                         data.table(Name=m,Score=s,
                                    Mean=mean(boot_data$t),SD=sd(boot_data$t),
                                    ci_L=NA,
                                    ci_R=NA))
    
  }
}; rm(boot_data,s,m) 

scores_boot[,ci_L := Mean - SD*qt(0.025,999,lower.tail = F)]
scores_boot[,ci_R := Mean + SD*qt(0.025,999,lower.tail = F)]

scores_boot[ci_L>0,Diff_from_zero:="+"]
scores_boot[ci_R<0,Diff_from_zero:="-"]
scores_boot[ci_L<=0 & ci_R>=0,Diff_from_zero:="0"]


scores_boot




# Remove everything apart from functions
rm(list = setdiff(ls(), lsf.str()))













## Wind - Temporal ####

# Load wind and solar data

load("data/windsolar_fc.rda")

zone_dat[,g_val:=qnorm(u_val)]
zone_dat[g_val%in%c(-Inf,Inf),g_val:=NA]



# issue time dependency, neet to 1) not to mix, or 2) model explicitly
gobs <- dcast(zone_dat[id=="wind_scotland",],
              kfold+issueTime~id+lead_time,
              value.var = c("g_val"),
              drop = TRUE)
setorder(gobs,issueTime)
gobs[,1:10]


# Convert to correlation matrix
WindScot_Cov <- cov(gobs[hour(issueTime)==0 & kfold!="Test",][,-c(1,2)],use = "pairwise.complete.obs")
WindScot_K <- diag(sqrt(1/diag(WindScot_Cov)))
WindScot_K_inv <- solve(WindScot_K)
WindScot_Cor <- WindScot_K %*% WindScot_Cov %*% WindScot_K
diag(WindScot_Cor) <- 1


col6 <- colorRampPalette(c("blue","cyan","yellow","red"))

r <- seq(0,48,by=0.5)
R <- as.matrix(dist(r))


## For correlation, need to fix sigma = 1, so covert to correlation function
PowExp_cor <- function(params,...){
  PowExp(params = list(sigma=1,theta=params[[1]],gamma=params[[2]]),...)
}


## Static fit
WindScot_static_fit <- gac(R = R,
                           X = list(),
                           Emp_Cov = WindScot_Cor,
                           cov_func = PowExp_cor,
                           param_eqns = list(~1,
                                             ~1),
                           loss="WLS")



## GAC fit

# Z = "Dist along diagonal"
N <- length(r)
Z <- matrix(rep(1:N,N) + rep(1:N,each=N) - 1, N, N)/(2*N-1)

N <- matrix(rep(r,each=length(r)),length(r),length(r))

WindScot_gac_fit <- gac(R = R,
                        X = list(x1=Z),
                        Emp_Cov = WindScot_Cor,
                        cov_func = PowExp_cor,
                        param_eqns = list(~s(x1,bs="cr",k=20),
                                          ~1),
                        loss="WLS",smoothness_param = 0.1)



## Plots
wind_cor_plot <- function(data,filename=NULL,...){
  h <- lattice::levelplot(data,
                          xlab = list(label="Lead-time [hours]",cex=1.2),
                          ylab = list(label="Lead-time [hours]",cex=1.2),
                          col.regions=col6(600), cuts=100, at=seq(-.1,1.05,0.01),
                          scales=list(x=list(at = seq(1,97,12),lab=seq(0,96,12)/2,rot=45),
                                      y=list(at = seq(1,97,12),lab=seq(0,96,12)/2,rot=45),tck=0.3,cex=1.2),
                          ...)
  if(!is.null(filename)){
    setEPS(); postscript(filename)
    print(h)
    dev.off()
  }
  
  print(h)
  
}

lattice::lattice.options(
  layout.heights=list(bottom.padding=list(x=1), top.padding=list(x=-1)),
  layout.widths=list(left.padding=list(x=-0.5), right.padding=list(x=-1))
)


wind_cor_plot(WindScot_Cor,colorkey=F,filename = "PSCC22_plots/Wind_Emp.eps")

wind_cor_plot(WindScot_static_fit$Cov_Est,colorkey=F,filename = "PSCC22_plots/Wind_Const.eps")
# colorkey = list(space = "bottom",width=0.8,labels=list(cex=1.2)))

wind_cor_plot(WindScot_gac_fit$Cov_Est,colorkey=F,filename = "PSCC22_plots/Wind_gac.eps")


# modelling_table <- data.frame(y=c(WindScot_Cor),
#                               r=c(R),
#                               x1=c(Z),
#                               y_est_static=c(WindScot_static_fit$Cov_Est),
#                               y_est_gac=c(WindScot_gac_fit$Cov_Est))
# 
# plot(x=modelling_table$r,
#      y=modelling_table$y,
#      col="grey",pch=".")
# 
# points(x=modelling_table$r,
#        y=modelling_table$y_est_static,
#        col="red")
# 
# points(x=modelling_table$r,
#        y=modelling_table$y_est_gac,
#        col="blue")
# 
# 
# scatter3D(modelling_table$r,
#           modelling_table$x1,
#           modelling_table$y,
#           col = "red",xlab="time diff",ylab="dist_along diag",zlab="covariance")
# 
# # plotrgl()
# 
# plot1 <- plotly::plot_ly(x=modelling_table$r,
#                          y=modelling_table$x1,
#                          z=modelling_table$y, type="scatter3d", mode="markers", color = I("grey"),
#                          marker = list(size = 1,opacity = 0.5),
#                          name = "empirical")
# 
# modelling_table$y_est <- c(WindScot_static_fit$Cov_Est)
# plot2 <- plotly::plot_ly(x=modelling_table$r,
#                          y=modelling_table$x1,
#                          z=modelling_table$y_est, type="scatter3d", mode="markers", color = I("red"),
#                          marker = list(size = 1,opacity = 0.5),
#                          name = "static")
# 
# 
# modelling_table$y_est <- c(WindScot_gac_fit$Cov_Est)
# plot3 <- plotly::plot_ly(x=modelling_table$r,
#                          y=modelling_table$x1,
#                          z=modelling_table$y_est, type="scatter3d", mode="markers", color = I("blue"),
#                          marker = list(size = 1,opacity = 0.5),
#                          name = "gac")
# 
# plot <- plotly::subplot(plot1,plot2,plot3,nrows = 1)
# plotly::layout(plot, scene = list(xaxis = list(title = 'dt [hours]'), yaxis = list(title = 'dist along diag [-]'), zaxis = list(title = 'covariance [-]')))
# invisible(gc())


## Evaluation

actuals <- as.matrix(na.omit(gobs[hour(issueTime)==0 & kfold=="Test",][,-c(1,2)]))

cov_mats <- list(Name=c("Empirical","Constant","GAC"),
                 mat=list(WindScot_K_inv %*% WindScot_Cor %*% WindScot_K_inv,
                          WindScot_K_inv %*% WindScot_static_fit$Cov_Est %*% WindScot_K_inv,
                          WindScot_K_inv %*% WindScot_gac_fit$Cov_Est %*% WindScot_K_inv))

scores <- data.table(index=rep(1:nrow(actuals),length(cov_mats$Name)),
                     Name=rep(cov_mats$Name,each=nrow(actuals)))

for(cov_i in 1:length(cov_mats$Name)){
  
  print(cov_mats$Name[cov_i])
  pb <- progress_bar$new(total = nrow(actuals))
  
  cov_mats$mat[[cov_i]] <- nearPD(cov_mats$mat[[cov_i]])$mat
  
  for(i in 1:nrow(actuals)){
    
    pb$tick()
    
    # Draw sample trajectories
    traj <- mvnfast::rmvn(n = 2500,mu=rep(0,ncol(cov_mats$mat[[cov_i]])),sigma = cov_mats$mat[[cov_i]])
    
    # es and vs
    scores[index==i & Name==cov_mats$Name[[cov_i]], es := es_sample(y=actuals[i,],dat = t(traj))]
    scores[index==i & Name==cov_mats$Name[[cov_i]], vs_0_5 := vs_sample_quick(y=actuals[i,],dat = t(traj), p = 0.5)]
    scores[index==i & Name==cov_mats$Name[[cov_i]], vs_1 := vs_sample_quick(y=actuals[i,],dat = t(traj), p = 1)]
    
    # log score
    scores[index==i & Name==cov_mats$Name[[cov_i]], ls := -log(mvnfast::dmvn(X = actuals[i,],mu = rep(0,ncol(actuals)),sigma = cov_mats$mat[[cov_i]]))]
    
  }
}


## Save/load results used in paper
# save.image(file="PSCC22_plots/wind_example.Rda")
# load("PSCC22_plots/wind_example.Rda")
##



print(setorder(scores[,mean(es),by="Name"],V1))
print(setorder(scores[,mean(vs_0_5),by="Name"],V1))
print(setorder(scores[,mean(vs_1),by="Name"],V1))
print(setorder(scores[,mean(ls),by="Name"],V1))


print(xtable(scores[,.(`Energy`=signif(mean(es),4),
                       `Log`=signif(mean(ls),4),
                       `VS-0.5`=signif(mean(vs_0_5),4),
                       `VS-1`=signif(mean(vs_1),4)),
                    by="Name"][order(-Log),],digits = c(NA,0,3,2,0,0),
             caption = c("Results for temporal wind power forecasting.  Italics  indicate  that  the  corresponding  skill  score  relative  to  the GAC  model are not significantly different from zero."),
             label = c("tab:wind_cov")),
      caption.placement = "top",table.placement="",
      include.rownames=F)


## Bootstrap skill score
scores_boot <- data.table()
ref <- "GAC"
for(m in scores[Name!=ref,unique(Name)]){
  for(s in names(scores)[-c(1:2)]){
    
    boot_data <- boot(data = cbind(scores[Name==m,get(s)],scores[Name==ref,get(s)]),
                      statistic = function(data, i){
                        100-100*mean(data[i,1])/mean(data[i,2])
                      },
                      R=1000)
    
    
    scores_boot <- rbind(scores_boot,
                         data.table(Name=m,Score=s,
                                    Mean=mean(boot_data$t),SD=sd(boot_data$t)))
    
  }
}; rm(boot_data,s,m) 

scores_boot[,ci_L := Mean - SD*qt(0.025,999,lower.tail = F)]
scores_boot[,ci_R := Mean + SD*qt(0.025,999,lower.tail = F)]

scores_boot[ci_L>0,Diff_from_zero:="+"]
scores_boot[ci_R<0,Diff_from_zero:="-"]
scores_boot[ci_L<=0 & ci_R>=0,Diff_from_zero:="0"]

scores_boot


# Remove everything apart from functions
rm(list = setdiff(ls(), lsf.str()))

