# Libraries
library(sbw)
library(ggplot2)
library(foreign)
library(gtable)
library(gridExtra)
library(devtools)

# Source files
source('implied weighting function repl.R')


## Lalonde data

# Lalonde Experimental Data - Full
exp.full = read.dta('nsw_dw.dta')


# Lalonde Observational Data
# read treatment group - same as experimental
obs.treated = read.table("http://users.nber.org/~rdehejia/data/nswre74_treated.txt")
# read control group
obs.control = read.table("http://users.nber.org/~rdehejia/data/psid_controls.txt")
# smaller group
#obs.control = read.table("http://users.nber.org/~rdehejia/data/cps3_controls.txt")


colnames(obs.treated) = colnames(exp.full)[-1]
colnames(obs.control) = colnames(exp.full)[-1]

# Full dataset
obs.full = rbind(obs.treated, obs.control)

X = as.matrix(obs.full[,2:9]) # covariates
Z = as.vector(obs.full[,1])   # treatment indicator
y = as.vector(obs.full[,10])  # outcome

# Data frame with covariates and treatment indicator
df.sample = cbind(X,Z)

N = nrow(df.sample) # full sample size
nc = sum(1-Z) # control group size 
nt = sum(Z)   # treatment group size
k = ncol(X)   # no. of covariates


########################################################################################
#########################################################################################

### MRI weights

w_trial_MRI2 = imp.weights(data = df.sample, treatment = 'Z', intercept = TRUE, 
                           WLS.weights = NULL, method = 'MRI', 
                           estimand = "ATE", target_val = NULL)

### URI weights

w_trial_URI2 = imp.weights(data = df.sample, treatment = 'Z', intercept = TRUE, 
                           WLS.weights = NULL, method = 'URI', 
                           estimand = "ATE", target_val = NULL)


###########################################################################
###########################################################################
## Plot of ASMD and TASMD

w = w_trial_MRI2

asmd_before_weighting = rep(0,k)
asmd_after_weighting = rep(0,k)

tasmd_before_weighting_control = rep(0,k)
tasmd_after_weighting_control = rep(0,k)

tasmd_before_weighting_treat = rep(0,k)
tasmd_after_weighting_treat = rep(0,k)


for(j in 1:k)
{
  V = (var(X[Z==1,j]) + var(X[Z==0,j]))/2
  V_target = var(X[,j])
  
  asmd_before_weighting[j] = abs(mean(X[Z==1,j]) - mean(X[Z==0,j]))/sqrt(V)
  asmd_after_weighting[j] = abs( sum(X[Z==1,j]*w[Z==1]) - sum(X[Z==0,j]*w[Z==0]))/sqrt(V)
  
  tasmd_before_weighting_treat[j] = abs(mean(X[Z==1,j]) - mean(X[,j]))/sqrt(V_target)
  tasmd_after_weighting_treat[j] = abs(sum(X[Z==1,j]*w[Z==1]) - mean(X[,j]))/sqrt(V_target)
  
  tasmd_before_weighting_control[j] = abs(mean(X[Z==0,j]) - mean(X[,j]))/sqrt(V_target)
  tasmd_after_weighting_control[j] = abs(sum(X[Z==0,j]*w[Z==0]) - mean(X[,j]))/sqrt(V_target)
  
}


ASMD.matrix.MRI = round(cbind(asmd_before_weighting, asmd_after_weighting),5)
colnames(ASMD.matrix.MRI) = c('Unweighted', 'MRI Weights')
rownames(ASMD.matrix.MRI) = colnames(X)

TASMD.matrix_treat.MRI = round(cbind(tasmd_before_weighting_treat, tasmd_after_weighting_treat),5)
colnames(TASMD.matrix_treat.MRI) = c('Unweighted', 'MRI Weights')
rownames(TASMD.matrix_treat.MRI) = colnames(X)

TASMD.matrix_control.MRI = round(cbind(tasmd_before_weighting_control, tasmd_after_weighting_control),5)
colnames(TASMD.matrix_control.MRI) = c('Unweighted', 'MRI Weights')
rownames(TASMD.matrix_control.MRI) = colnames(X)


#########################################################################

w = w_trial_URI2

asmd_before_weighting = rep(0,k)
asmd_after_weighting = rep(0,k)

tasmd_before_weighting_control = rep(0,k)
tasmd_after_weighting_control = rep(0,k)

tasmd_before_weighting_treat = rep(0,k)
tasmd_after_weighting_treat = rep(0,k)


for(j in 1:k)
{
  V = (var(X[Z==1,j]) + var(X[Z==0,j]))/2
  V_target = var(X[,j])
  
  asmd_before_weighting[j] = abs(mean(X[Z==1,j]) - mean(X[Z==0,j]))/sqrt(V)
  asmd_after_weighting[j] = abs( sum(X[Z==1,j]*w[Z==1]) - sum(X[Z==0,j]*w[Z==0]))/sqrt(V)
  
  tasmd_before_weighting_treat[j] = abs(mean(X[Z==1,j]) - mean(X[,j]))/sqrt(V_target)
  tasmd_after_weighting_treat[j] = abs(sum(X[Z==1,j]*w[Z==1]) - mean(X[,j]))/sqrt(V_target)
  
  tasmd_before_weighting_control[j] = abs(mean(X[Z==0,j]) - mean(X[,j]))/sqrt(V_target)
  tasmd_after_weighting_control[j] = abs(sum(X[Z==0,j]*w[Z==0]) - mean(X[,j]))/sqrt(V_target)
  
}


ASMD.matrix.URI = round(cbind(asmd_before_weighting, asmd_after_weighting),5)
colnames(ASMD.matrix.URI) = c('Unweighted', 'URI Weights')
rownames(ASMD.matrix.URI) = colnames(X)

TASMD.matrix_treat.URI = round(cbind(tasmd_before_weighting_treat, tasmd_after_weighting_treat),5)
colnames(TASMD.matrix_treat.URI) = c('Unweighted', 'URI Weights')
rownames(TASMD.matrix_treat.URI) = colnames(X)

TASMD.matrix_control.URI = round(cbind(tasmd_before_weighting_control, tasmd_after_weighting_control),5)
colnames(TASMD.matrix_control.URI) = c('Unweighted', 'URI Weights')
rownames(TASMD.matrix_control.URI) = colnames(X)


##########################################################################
##########################################################################

## Absolute SIC for MRI and URI

X_treat = as.matrix(X[Z==1,]) # matrix of covariates in the treatment group
X_control = as.matrix(X[Z==0,]) # matrix of covariates in the control group

# Full-sample design matrix (with intercept)
X_tilde = as.matrix(cbind(rep(1,N), X, Z))

# Hat matrices
H_treat = X_treat %*% solve(t(X_treat)%*%X_treat)%*%t(X_treat) # treatment
H_control = X_control %*% solve(t(X_control)%*%X_control)%*%t(X_control) # control

H_tilde = X_tilde %*% solve(t(X_tilde)%*%X_tilde)%*%t(X_tilde) # full

h_MRI = c(diag(H_treat), diag(H_control))

## MRI fits and errors
fit_treat = lm(y[Z==1]~X_treat)
fit_control = lm(y[Z==0]~X_control)

# Residuals under MRI and URI
e_MRI = rep(0,N)
e_MRI[1:nt] = fit_treat$residuals
e_MRI[(nt+1):N] = fit_control$residuals

fit_full = lm(y~X_tilde[,-1])
e_URI = fit_full$residuals


# Sample influence curves under MRI and URI
a_SIC_MRI = rep(0,N)
a_SIC_URI = rep(0,N)

for(i in 1:N)
{
  a_SIC_MRI[i] = (nt*Z[i] + (1-Z[i])*nc)*abs(e_MRI[i]*w_trial_MRI2[i])/(1-h_MRI[i])
  a_SIC_URI[i] = (N-1)*abs(e_URI[i]*w_trial_URI2[i])/(1-H_tilde[i,i])
}


## Scaled SIC
a_SIC_MRI_std = a_SIC_MRI/max(a_SIC_MRI)
a_SIC_URI_std = a_SIC_URI/max(a_SIC_URI)

################################################################################################################
################################################################################################################
################################################################################################################
pdf(file = "Plots_combined_repl.pdf",
    width=6.5,
    height=9)
layout.mat <- rbind(c(1, 1, 4, 4), 
                    c(1, 1, 4, 4),
                    c(1, 1, 4, 4),
                    c(1, 1, 5, 5),
                    c(2, 2, 5, 5),
                    c(2, 2, 5, 5),
                    c(2, 2, 6, 6),
                    c(2, 2, 6, 6),
                    c(3, 3, 6, 6),
                    c(3, 3, 7, 7),
                    c(3, 3, 7, 7),
                    c(3, 3, 7, 7),
                    c(8, 9, 16, 16),
                    c(8, 9, 16, 16),
                    c(8, 9, 16, 16),
                    c(10, 11, 16, 16),
                    c(10, 11, 16, 16),
                    c(10, 11, 16, 16),
                    c(12, 13, 16, 16),
                    c(12, 13, 16, 16),
                    c(12, 13, 16, 16),
                    c(14, 15, 16, 16),
                    c(14, 15, 16, 16),
                    c(14, 15, 16, 16))
layout(layout.mat)
par(#mfrow=c(3,1), 
  mar=c(2.5, 5.5, 2, 0.5),
  mgp=c(1.5, 0.5, 0))
## ASMD
plot(y=1:nrow(ASMD.matrix.URI), 
     x=rev(ASMD.matrix.URI[,2]),
     pch=3,
     xlim=c(0, 3),
     ylim=c(0, nrow(ASMD.matrix.URI)+0.5),
     yaxt="n",
     ylab="",
     xlab="ASMD",
     cex=1.5,
     cex.axis=0.8,
     cex.lab=0.8)
title(main = "Balance and representativeness", outer=TRUE, adj=0.0125, line=-1, cex.main=1, font.main=1)
points(y=1:nrow(ASMD.matrix.URI), 
       x=rev(ASMD.matrix.URI[,1]),
       pch=4,
       cex=1.5)
points(y=1:nrow(ASMD.matrix.MRI), 
       x=rev(ASMD.matrix.MRI[,2]),
       pch=20,
       cex=1.5)
legend(x=1.9,
       y=nrow(ASMD.matrix.URI)+0.65,
       legend=c("Before weighting",
                "After URI weighting",
                "After MRI weighting"),
       pch=c(4, 3, 20),
       bty="n",
       cex=0.8) 
axis(side=2,
     las=2,
     at=1:nrow(ASMD.matrix.URI),
     labels=rev(rownames(ASMD.matrix.URI)),
     cex.axis=0.8,
     mgp=c(3, 1.25, 0))
mtext(text = "Covariate",
      side = 2,
      line = 4.5,
      cex=0.55)

## TASMD - treatment
plot(y=1:nrow(TASMD.matrix_treat.URI), 
     x=rev(TASMD.matrix_treat.URI[,2]),
     pch=3,
     xlim=c(0, 3),
     ylim=c(0, nrow(TASMD.matrix_treat.URI)+0.5),
     yaxt="n",
     ylab="",
     xlab="TASMD (treatment group)",
     cex=1.5,
     cex.axis=0.8,
     cex.lab=0.8)
points(y=1:nrow(TASMD.matrix_treat.URI), 
       x=rev(TASMD.matrix_treat.URI[,1]),
       pch=4,
       cex=1.5)
points(y=1:nrow(TASMD.matrix_treat.MRI), 
       x=rev(TASMD.matrix_treat.MRI[,2]),
       pch=20,
       cex=1.5)
legend(x=1.9,
       y=nrow(TASMD.matrix_treat.URI)+0.65,
       legend=c("Before weighting",
                "After URI weighting",
                "After MRI weighting"),
       pch=c(4, 3, 20),
       bty="n",
       cex=0.8) 
axis(side=2,
     las=2,
     at=1:nrow(ASMD.matrix.URI),
     labels=rev(rownames(ASMD.matrix.URI)),
     cex.axis=0.8,
     mgp=c(3, 1.25, 0))
mtext(text = "Covariate",
      side = 2,
      line = 4.5,
      cex=0.55)

## TASMD - control
plot(y=1:nrow(TASMD.matrix_control.URI), 
     x=rev(TASMD.matrix_control.URI[,2]),
     pch=3,
     xlim=c(0, 3),
     ylim=c(0, nrow(TASMD.matrix_control.URI)+0.5),
     yaxt="n",
     ylab="",
     xlab="TASMD (control group)",
     cex=1.5,
     cex.axis=0.8,
     cex.lab=0.8)
points(y=1:nrow(TASMD.matrix_control.URI), 
       x=rev(TASMD.matrix_control.URI[,1]),
       pch=4,
       cex=1.5)
points(y=1:nrow(TASMD.matrix_control.MRI), 
       x=rev(TASMD.matrix_control.MRI[,2]),
       pch=20,
       cex=1.5)
legend(x=1.9,
       y=nrow(TASMD.matrix_control.URI)+0.65,
       legend=c("Before weighting",
                "After URI weighting",
                "After MRI weighting"),
       pch=c(4, 3, 20),
       bty="n",
       cex=0.8) 
axis(side=2,
     las=2,
     at=1:nrow(ASMD.matrix.URI),
     labels=rev(rownames(ASMD.matrix.URI)),
     cex.axis=0.8,
     mgp=c(3, 1.25, 0))
mtext(text = "Covariate",
      side = 2,
      line = 4.5,
      cex=0.55)

### Distribution of the weights
### URI weights - treatment group
xmin <- min(min(w_trial_URI2[Z==1]),
            min(w_trial_MRI2[Z==1]),
            min(w_trial_URI2[Z==0]),
            min(w_trial_MRI2[Z==0]))
xmax <- max(max(w_trial_URI2[Z==1]),
            max(w_trial_MRI2[Z==1]),
            max(w_trial_URI2[Z==0]),
            max(w_trial_MRI2[Z==0]))
ess_t.URI = (sum(abs(w_trial_URI2[Z==1]))^2)/sum(w_trial_URI2[Z==1]^2)
ess_c.URI = (sum(abs(w_trial_URI2[Z==0]))^2)/sum(w_trial_URI2[Z==0]^2)
SS_t.URI = length(w_trial_URI2[Z==1])
SS_c.URI = length(w_trial_URI2[Z==0])

ess_t.MRI = (sum(abs(w_trial_MRI2[Z==1]))^2)/sum(w_trial_MRI2[Z==1]^2)
ess_c.MRI = (sum(abs(w_trial_MRI2[Z==0]))^2)/sum(w_trial_MRI2[Z==0]^2)
SS_t.MRI = length(w_trial_MRI2[Z==1])
SS_c.MRI = length(w_trial_MRI2[Z==0])

par(#mfrow=c(4,1), 
  mar=c(2.5, 3, 3, 1),
  mgp=c(1.33, 0.5, 0))
plot(density(w_trial_URI2[Z==1]),
     ylab="Density",
     xlab="URI weights",
     #xlim=c(xmin, xmax),
     yaxt="n",
     cex.lab=0.8,
     cex.axis=0.8,
     main=""
)
title(main="URI Weights: Treatment Group",
      cex.main=0.9,
      font.main=1,
      line=0.5)
title(main = "Dispersion and effective sample size", outer=TRUE, adj=0.7125, line=-1, cex.main=1, font.main=1)
axis(side=2,
     at=c(0, 800),
     cex.axis=0.8)
rug(w_trial_URI2[Z==1], 
    col = "black", 
    ticksize=0.1,
    lwd = 0.3, 
    side = 1)
legend(x=0.006,
       y=800,
       legend=paste0("Effective sample size = ", 
                     round(ess_t.URI, 1)),
       bty="n",
       cex=0.8)
legend(x=0.006,
       y=575,
       legend=paste0("Sample size = ", 
                     round(SS_t.URI, 1)),
       bty="n",
       cex=0.8)
plot(density(w_trial_MRI2[Z==1]),
     ylab="Density",
     xlab="MRI weights",
     main="",
     #xlim=c(xmin, xmax)
     yaxt="n",
     cex.lab=0.8,
     cex.axis=0.8)
title(main="MRI Weights: Treatment Group",
      cex.main=0.9,
      font.main=1,
      line=0.5)
rug(w_trial_MRI2[Z==1], 
    col = "black", 
    ticksize=0.1,
    lwd = 0.3, 
    side = 1)
axis(side=2, at=c(0, 20),
     cex.axis=0.8)
legend(x=0.10,
       y=27,
       legend=paste0("Effective sample size = ", 
                     round(ess_t.MRI, 1)),
       bty="n",
       cex=0.8)
legend(x=0.10,
       y=19,
       legend=paste0("Sample size = ", 
                     round(SS_t.MRI, 1)),
       bty="n",
       cex=0.8)
plot(density(w_trial_URI2[Z==0]),
     ylab="Density",
     xlab="URI weights",
     main="",
     #xlim=c(xmin, xmax)
     yaxt="n",
     cex.lab=0.8,
     cex.axis=0.8)
title(main="URI Weights: Control Group",
      cex.main=0.9,
      font.main=1,
      line=0.5)
axis(side=2, at=c(0, 500),
     cex.axis=0.8)
rug(w_trial_URI2[Z==0], 
    col = "black", 
    ticksize=0.1,
    lwd = 0.3, 
    side = 1)
legend(x=0,
       y=650,
       legend=paste0("Effective sample size = ", 
                     round(ess_c.URI, 1)),
       bty="n",
       cex=0.8)
legend(x=0,
       y=450,
       legend=paste0("         Sample size = ", 
                     round(SS_c.URI, 1)),
       bty="n",
       cex=0.8)
plot(density(w_trial_MRI2[Z==0]),
     ylab="Density",
     xlab="MRI weights",
     main="",
     #xlim=c(xmin, xmax),
     yaxt="n",
     cex.lab=0.8,
     cex.axis=0.8)
title(main="MRI Weights: Control Group",
      cex.main=0.9,
      font.main=1,
      line=0.5)
axis(side=2, at=c(0, 8000),
     cex.axis=0.8)
rug(w_trial_MRI2[Z==0], 
    col = "black", 
    ticksize=0.1,
    lwd = 0.3, 
    side = 1)
legend(x=4e-04,
       y=9000,
       legend=paste0("Effective sample size = ", 
                     round(ess_c.MRI, 1)),
       bty="n",
       cex=0.8)
legend(x=4e-04,
       y=6500,
       legend=paste0("      Sample size = ", 
                     round(SS_c.MRI, 1)),
       bty="n",
       cex=0.8)

### Positive/negative weights plots
## Age
## URI - treatment group
set.seed(1234)
black.URI <- cbind(df.sample[, c("black", "Z")], w_trial_URI2)
black.MRI <- cbind(df.sample[, c("black", "Z")], w_trial_MRI2)
re75.URI <- cbind(df.sample[, c("re75", "Z")], w_trial_URI2)
re75.MRI <- cbind(df.sample[, c("re75", "Z")], w_trial_MRI2)

par(mar=c(2.5, 0.5, 3.5, 0.5),
    mgp=c(1.5, 0.5, 0))
plot(x=black.URI[Z == 1, "black"] + rnorm(n=nrow(black.URI[Z == 1, ]), mean=0, sd=0.01),
     y=ifelse(black.URI[Z == 1, "w_trial_URI2"] < 0, 1, 2) + rnorm(n=nrow((black.URI[Z == 1,])), mean=0, sd=0.05), 
     pch=20,
     col=ifelse(black.URI[Z == 1, "w_trial_URI2"] < 0, "red", "black"),
     ylim=c(0.5, 2.5),
     xlab="Black (0/1)",
     yaxt="n",
     cex=300*abs(black.URI[Z == 1,"w_trial_URI2"])/sum(abs(black.URI[Z == 1,"w_trial_URI2"])),
     ylab="",
     main="",
     frame.plot = TRUE,
     cex.axis=0.8,
     cex.lab=0.8)
title(main="URI Weights: Treatment Group",
      cex.main=0.9,
      font.main=1,
      line=0.5)
title(main="Negative weights and extrapolation",
      outer=TRUE,
      adj=0.0125,
      line=-35.5, cex.main=1, font.main=1)
abline(v=weighted.mean(black.URI[Z == 1,"black"], black.URI[Z == 1, "w_trial_URI2"]))
points(x=mean(df.sample[, "black"]), y=1.5, pch="*", cex=4, lwd=0.5)

plot(x=re75.URI[Z == 1, "re75"],
     y=ifelse(re75.URI[Z == 1, "w_trial_URI2"] < 0, 1, 2) + rnorm(n=nrow((re75.URI[Z == 1,])), mean=0, sd=0.05), 
     pch=20,
     col=ifelse(re75.URI[Z == 1, "w_trial_URI2"] < 0, "red", "black"),
     ylim=c(0.5, 2.5),
     xlab="Earnings '75",
     yaxt="n",
     cex=300*abs(re75.URI[Z == 1,"w_trial_URI2"])/sum(abs(re75.URI[Z == 1,"w_trial_URI2"])),
     ylab="",
     main="",
     frame.plot = TRUE,
     font.main=1,
     cex.main=0.9,
     cex.axis=0.8,
     cex.lab=0.8)
title(main="URI Weights: Treatment Group",
      cex.main=0.9,
      font.main=1,
      line=0.5)
abline(v=weighted.mean(re75.URI[Z == 1,"re75"], re75.URI[Z == 1, "w_trial_URI2"]))
points(x=mean(df.sample[, "re75"]), y=1.5, pch="*", cex=4, lwd=0.5)

plot(x=black.MRI[Z == 1, "black"] + rnorm(n=nrow(black.MRI[Z == 1, ]), mean=0, sd=0.01),
     y=ifelse(black.MRI[Z == 1, "w_trial_MRI2"] < 0, 1, 2) + rnorm(n=nrow((black.MRI[Z == 1,])), mean=0, sd=0.05), 
     pch=20,
     col=ifelse(re75.MRI[Z == 1, "w_trial_MRI2"] < 0, "red", "black"),
     ylim=c(0.5, 2.5),
     xlab="Black (0/1)",
     yaxt="n",
     cex=70*abs(black.MRI[Z == 1,"w_trial_MRI2"])/sum(abs(black.MRI[Z == 1,"w_trial_MRI2"])),
     ylab="",
     main="",
     frame.plot = TRUE,
     font.main=1,
     cex.main=0.9,
     cex.axis=0.8,
     cex.lab=0.8)
title(main="MRI Weights: Treatment Group",
      cex.main=0.9,
      font.main=1,
      line=0.5)
abline(v=weighted.mean(black.MRI[Z == 1,"black"], black.MRI[Z == 1, "w_trial_MRI2"]))
points(x=mean(df.sample[, "black"]), y=1.5, pch="*", cex=4, lwd=0.5)

plot(x=re75.MRI[Z == 1, "re75"],
     y=ifelse(re75.MRI[Z == 1, "w_trial_MRI2"] < 0, 1, 2) + rnorm(n=nrow((re75.MRI[Z == 1,])), mean=0, sd=0.05), 
     pch=20,
     col=ifelse(re75.MRI[Z == 1, "w_trial_MRI2"] < 0, "red", "black"),
     ylim=c(0.5, 2.5),
     xlab="Earnings '75",
     yaxt="n",
     cex=70*abs(re75.MRI[Z == 1,"w_trial_MRI2"])/sum(abs(re75.MRI[Z == 1,"w_trial_MRI2"])),
     ylab="",
     main="",
     frame.plot = TRUE,
     font.main=1,
     cex.main=0.9,
     cex.axis=0.8,
     cex.lab=0.8)
title(main="MRI Weights: Treatment Group",
      cex.main=0.9,
      font.main=1,
      line=0.5)
abline(v=weighted.mean(re75.MRI[Z == 1,"re75"], re75.MRI[Z == 1, "w_trial_MRI2"]))
points(x=mean(df.sample[, "re75"]), y=1.5, pch="*", cex=4, lwd=0.5)

plot(x=black.URI[Z == 0, "black"] + rnorm(n=nrow(black.URI[Z == 0, ]), mean=0, sd=0.01),
     y=ifelse(black.URI[Z == 0, "w_trial_URI2"] < 0, 1, 2) + rnorm(n=nrow((black.URI[Z == 0,])), mean=0, sd=0.05), 
     pch=20,
     col=ifelse(re75.URI[Z == 0, "w_trial_URI2"] < 0, "red", "black"),
     ylim=c(0.5, 2.5),
     xlab="Black (0/1)",
     yaxt="n",
     cex=2000*abs(black.URI[Z == 0,"w_trial_URI2"])/sum(abs(black.URI[Z == 0,"w_trial_URI2"])),
     ylab="",
     main="",
     frame.plot = TRUE,
     font.main=1,
     cex.main=0.9,
     cex.axis=0.8,
     cex.lab=0.8)
title(main="URI Weights: Control Group",
      cex.main=0.9,
      font.main=1,
      line=0.5)
abline(v=weighted.mean(black.URI[Z == 0,"black"], black.URI[Z == 0, "w_trial_URI2"]))
points(x=mean(df.sample[, "black"]), y=1.5, pch="*", cex=4, lwd=0.5)

plot(x=re75.URI[Z == 0, "re75"],
     y=ifelse(re75.URI[Z == 0, "w_trial_URI2"] < 0, 1, 2) + rnorm(n=nrow((re75.URI[Z == 0,])), mean=0, sd=0.05), 
     pch=20,
     col=ifelse(re75.URI[Z == 0, "w_trial_URI2"] < 0, "red", "black"),
     ylim=c(0.5, 2.5),
     xlab="Earnings '75",
     yaxt="n",
     cex=2000*abs(re75.URI[Z == 0,"w_trial_URI2"])/sum(abs(re75.URI[Z == 0,"w_trial_URI2"])),
     ylab="",
     main="",
     frame.plot = TRUE,
     font.main=1,
     cex.main=0.9,
     cex.axis=0.8,
     cex.lab=0.8)
title(main="URI Weights: Control Group",
      cex.main=0.9,
      font.main=1,
      line=0.5)
abline(v=weighted.mean(re75.URI[Z == 0,"re75"], re75.URI[Z == 0, "w_trial_URI2"]))
points(x=mean(df.sample[,"re75"]), y=1.5, pch="*", cex=4, lwd=0.5)

plot(x=black.MRI[Z == 0, "black"] + rnorm(n=nrow(black.MRI[Z == 0, ]), mean=0, sd=0.01),
     y=ifelse(black.MRI[Z == 0, "w_trial_MRI2"] < 0, 1, 2) + rnorm(n=nrow((black.MRI[Z == 0,])), mean=0, sd=0.05), 
     pch=20,
     col=ifelse(re75.MRI[Z == 0, "w_trial_MRI2"] < 0, "red", "black"),
     ylim=c(0.5, 2.5),
     xlab="Black (0/1)",
     yaxt="n",
     cex=2000*abs(black.MRI[Z == 0,"w_trial_MRI2"])/sum(abs(black.MRI[Z == 0,"w_trial_MRI2"])),
     ylab="",
     main="",
     frame.plot = TRUE,
     font.main=1,
     cex.main=0.9,
     cex.axis=0.8,
     cex.lab=0.8)
title(main="MRI Weights: Control Group",
      cex.main=0.9,
      font.main=1,
      line=0.5)
abline(v=weighted.mean(black.MRI[Z == 0,"black"], black.MRI[Z == 0, "w_trial_MRI2"]))
points(x=mean(df.sample[, "black"]), y=1.5, pch="*", cex=4, lwd=0.5)

plot(x=re75.MRI[Z == 0, "re75"],
     y=ifelse(re75.MRI[Z == 0, "w_trial_MRI2"] < 0, 1, 2) + rnorm(n=nrow((re75.MRI[Z == 0,])), mean=0, sd=0.05), 
     pch=20,
     col=ifelse(re75.MRI[Z == 0, "w_trial_MRI2"] < 0, "red", "black"),
     ylim=c(0.5, 2.5),
     xlab="Earnings '75",
     yaxt="n",
     cex=2000*abs(re75.MRI[Z == 0,"w_trial_MRI2"])/sum(abs(re75.MRI[Z == 0,"w_trial_MRI2"])),
     ylab="",
     main="",
     frame.plot = TRUE,
     font.main=1,
     cex.main=0.9,
     cex.axis=0.8,
     cex.lab=0.8)
title(main="MRI Weights: Control Group",
      cex.main=0.9,
      font.main=1,
      line=0.5)
abline(v=weighted.mean(re75.MRI[Z == 0,"re75"], re75.MRI[Z == 0, "w_trial_MRI2"]))
points(x=mean(df.sample[, "re75"]), y=1.5, pch="*", cex=4, lwd=0.5)


# Sample influence curves under MRI and URI

# Plot of SIC
par(mar=c(2.5, 3, 2, 2),
    mgp=c(1.5, 1, 0))
matplot(1:N, 
        cbind(-a_SIC_MRI_std,a_SIC_URI_std), 
        type = c("h","h"),
        lty = c("solid","solid"),
        xlab = "Index", 
        ylab = "Scaled SIC", 
        cex.lab = 0.8,
        col = grey(c(0,0.6)), 
        axes = FALSE)
title(main="Influence",
      outer=TRUE,
      adj=0.555,
      line=-35.5, cex.main=1, font.main=1)
axis(side = 1, 
     at = c(1,seq(200,2500,200),N),
     cex.axis = 0.8,
     mgp=c(0, 0.5, 0),
     cex.axis=0.7)
axis(side = 2,
     ylim = range(a_SIC_URI_std), 
     at = seq(0,1,0.2), 
     cex.axis = 0.7,
     pos=-50)
axis(side = 4, ylim = -range(a_SIC_MRI_std), 
     at = -seq(0,1,0.2), 
     labels = seq(0,1,0.2), 
     cex.axis = 0.7,
     pos=2750)
legend(x=2200,
       y=1,
       legend=c("URI", "MRI"), 
       lty=c(1,1), 
       pch=c(NA, NA), 
       col=grey(c(0.6,0)), 
       cex = 0.6,
       bty="n")
par(xpd=NA)
abline(v=-500, lwd=0.1)
abline(h=1.2, lwd=0.1)
dev.off()
