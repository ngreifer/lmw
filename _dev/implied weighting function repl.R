library(scatterplot3d)

# Function to compute the vector of URI/MRI weights for two treatment groups.
imp.weights = function(data, treatment, outcome = NULL, intercept = TRUE, WLS.weights = NULL, method = 'MRI', 
                       estimand = "ATE", target_val = NULL)
{
  data = as.matrix(data) # convert the data.frame into matrix 
  Z = data[,treatment] # vector of treatment indicator in the whole sample
  
  # matrix of covariates in the whole sample
  if(is.null(outcome) == TRUE)
  {
    X = subset(data, select = -which(colnames(data) == treatment)) 
  } else {
    X = subset(data, select = -which(colnames(data) %in% c(treatment, outcome)))
  }
  
  n = nrow(data) # sample size
  nt = sum(Z) # treatment group size
  nc = n - nt # control group size
  k = ncol(X) # number of baseline covariates
  
  
  # matrix of covariates in the treatment and control group.
  Xt = X[Z==1,]
  Xc = X[Z==0,]
  
  # create the design matrices in the treatment group, control group and overall sample
  if(intercept == TRUE)
  {
    X_tilde = cbind(int = rep(1,n), X)
    X_tilde_t = X_tilde[Z==1,]
    X_tilde_c = X_tilde[Z==0,]
  }
  
  if(intercept == FALSE)
  {
    X_tilde = X
    X_tilde_t = Xt
    X_tilde_c = Xc
  }
  
  if(method == 'URI') 
  {
    # fit a linear regression model of Y on X_tilde and Z in the whole sample
    # no matter what the estimand is, the estimator is tau_hat (coefficient of Z) 
    G = cbind(X_tilde, Z)
    
    # incorporating the WLS weights
    if(is.null(WLS.weights) == TRUE)
    {
      D = diag((1/n), n) 
    } else{
      D = diag(WLS.weights)
    }
    
    MAT = t(D) %*% G %*% solve(t(G)%*% D %*% G) 
    w_tilde = MAT[,ncol(MAT)]
    w = (2*Z - 1)* w_tilde # final weights (sign changed for control group)
    
  }
  
  if(method == 'MRI')
  {
    w_c = rep(0,nc) # initialize vector of treatment weights
    w_t = rep(0,nt) # initialize vector of control weights
    if(estimand == 'ATE')
    {
      x_star = colMeans(X_tilde)
    }
    if(estimand == 'ATT')
    {
      x_star = colMeans(X_tilde[Z==1,])
    }
    if(estimand == 'ATC')
    {
      x_star = colMeans(X_tilde[Z==0,])
    }
    if(estimand == 'CATE')
    {
      if(intercept == TRUE)
      {
        x_star = c(1,target_val)
      }
      if(intercept == FALSE)
      {
        x_star = target_val
      }
    }  
    
    # incorporating the WLS weights
    if(is.null(WLS.weights) == TRUE)
    {
      Dt = diag((1/nt), nt)
      Dc = diag((1/nc), nc)
    } else{
      Dt = diag(WLS.weights[Z==1])
      Dc = diag(WLS.weights[Z==0])
    }
    
    # calculation of weights
    for(i in 1:nt)
    {
      w_t[i] = Dt[i,i] * (t(X_tilde_t[i,]) %*% solve(t(X_tilde_t)%*% Dt %*% X_tilde_t) %*% as.matrix(x_star))
      
    }
    for(i in 1:nc)
    {
      w_c[i] = Dc[i,i]* (t(X_tilde_c[i,]) %*% solve(t(X_tilde_c)%*% Dc %*% X_tilde_c) %*% as.matrix(x_star))
    }
    
    w = rep(0,n) # vector of weights for the full sample
    w[Z==1] = w_t
    w[Z==0] = w_c
    
    if(estimand == 'ATT') # no model fitted in treatment group, i.e. we use uniform weights in the treatment group
    {
      w[Z==1] = 1/nt
    }
    if(estimand == 'ATC') # no model fitted in control group, i.e. we use uniform weights in the control group
    {
      w[Z==0] = 1/nc
    }
  }
  
  # return a list of objects
  # return the vector of weights
  
  # calculate the total ESS = treatment ESS + control ESS
  #ess.weights = ess.gen(weights.treat = w[Z==1], weights.control = w[Z==0])
  
  # wave plot of the weights
  return(w)
  
  
}




### Examples

#df.sample = cbind(X,Z,y)
#w_trial_MRI = imp.weights(data = df.sample, treatment = 'Z', outcome = 'y', intercept = TRUE, WLS.weights = NULL, method = 'MRI', 
#                      estimand = "ATE", target_val = NULL)

#w_trial_MRI2 = imp.weights(data = df.sample, treatment = 'Z', outcome = 'y', intercept = TRUE, WLS.weights = NULL, method = 'MRI', 
#                         estimand = "CATE", target_val = colMeans(X))

#w_trial_MRI3 = imp.weights(data = df.sample, treatment = 'Z', outcome = 'y', intercept = TRUE, WLS.weights = NULL, method = 'MRI', 
#                          estimand = "ATT", target_val = NULL)

#w_trial_MRI3 = imp.weights(data = df.sample, treatment = 'Z', outcome = 'y', intercept = TRUE, WLS.weights = NULL, method = 'MRI', 
#                          estimand = "ATT", target_val = NULL)

# example on wls
#xx = rpois(nc,10)
#ww = rep(1/nt, n)
#ww[Z==0] = xx/sum(xx)
#w_trial_MRI4 = imp.weights(data = df.sample, treatment = 'Z', outcome = 'y', intercept = TRUE, WLS.weights = ww, method = 'MRI', 
#                          estimand = "ATT", target_val = NULL)



#w_trial_URI.1 = imp.weights(data = df.sample, treatment = 'Z', outcome = 'y', intercept = TRUE, WLS.weights = NULL, method = 'URI', 
#                         estimand = "ATE", target_val = NULL)

#w_trial_URI.2 = imp.weights(data = df.sample, treatment = 'Z', outcome = 'y', intercept = TRUE, WLS.weights = rep((1/n),n), method = 'URI', 
#                           estimand = "ATT", target_val = NULL)



