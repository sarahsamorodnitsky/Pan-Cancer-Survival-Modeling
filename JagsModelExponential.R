# Script for exponential survival model fit in JAGS 
# Author: Sarah Samorodnitsky and Eric Lock

library(rjags)
library(doParallel)
library(foreach)

load("FSCG.rda")

# For scaling
CenterByMean = function(vec) {
  # centering the variable
  new_vec = vec - mean(na.omit(vec))
  return(new_vec)
}

# Generate 5 fold training sets
Generate5FoldTrainingSet = function(X, Y, cancer_types, n_vec) {
  # Generates 5 folds of training values for each cancer type
  # non-randomly because of the seed
  
  InitManyLists = function(n) {
    new_list = list()
    for (i in 1:n) {
      new_list[[i]] = list()
    }
    return(new_list)
  }
  
  set.seed(1) # setting a seed so the same training data is produced each time
  
  # subsetting Y into training set
  prop = 0.2 # proportion to go into test set
  n = length(cancer_types)
  Training_obs = InitManyLists(n)
  names(Training_obs) = cancer_types
  for (k in 1:n) { # for each cancer type
    current_n = n_vec[k] # size of current cancer type
    reorder_ind = sample(1:current_n, size = current_n, replace = F) # reorder the observations
    current_test_n = round(prop*current_n) # how many observations should go in the test set
    fold_ind = seq(1, current_n, by = current_test_n) # the indices of the start of each test set
    all_ind = 1:current_n
    
    test_folds = lapply(1:length(fold_ind), function(i) { # the indices of observations for the test set
      if (i < length(fold_ind)) {
        test_ind = seq(fold_ind[i], fold_ind[i + 1] - 1) # get a sequence of indices for the test set
        reorder_ind[test_ind] # select from the reordered observations
      } else {
        test_ind = seq(fold_ind[i], current_n) # get a sequence of indices for the test set
        reorder_ind[test_ind] # select from the reordered observations
      }
    })
    
    if (length(fold_ind) > 5) { # if the set of obs was not evenly split into 5 groups
      extra_obs = test_folds[[6]]
      for (i in 1:length(extra_obs)) {
        test_folds[[i %% 5]] = c(test_folds[[i %% 5]], extra_obs[i]) # add an extra obs to another fold in a circular fashion
      }
    }
    test_folds = test_folds[-6] # get rid of the extra set of obs
    
    training_folds = lapply(1:length(test_folds), function(i) { # indices for observations in the training set
      all_ind[!(all_ind %in% test_folds[[i]])]
    })
    
    Training_obs[[k]] = training_folds
  }
  return(list(Training_obs = Training_obs))
}

# Function for calculating the posterior likelihood
PosteriorLikelihood = function(betas, X, Y, Last_Contact, Training_obs) {
  # posteriors is the result from the model computation
  # X and Y are the datasets, not in the form of an intercept and centered age already
  logSum <- function(l) {
    # given l = log(v), where v is a vector,
    # computes log(sum(v))
    return(max(l) + log(sum(exp(l-max(l)))))
  }
  
  CenterByMean = function(vec) {
    # centering the variable
    new_vec = vec - mean(na.omit(vec))
    return(new_vec)
  }
  
  LogExponentialLikelihood = function(x_ij, y_ij, y_ijc, betas_it) {
    # For testing
    x_ij = List_XY[[k]][i, -c(y_col, yc_col)]
    y_ij = List_XY[[k]][i, y_col] # for each patient in the test set
    y_ijc = List_XY[[k]][i, yc_col] # ith observation, censor time
    betas_it = betas[[k]][[j]]
    
    # calculate the likelihood for just one set of posteriors
    rate_ij = 1/exp(sum(x_ij * betas_it))
    
    if (!is.na(y_ij)) { # y_ij is not censored
      post_ij = dexp(y_ij, rate = rate_ij, log = TRUE)
      
    } else { # if y_ij is censored, use censoring time
      post_ij = log(pexp(y_ijc, rate = rate_ij, lower.tail = FALSE)) 
    }
    
    return(post_ij)
  }
  
  # Setting up X 
  n = length(X)
  # Centering age
  for (i in 1:n) {
    if (is.null(nrow(X[[i]]))) {
      current = X[[i]]
      current.2 = CenterByMean(current)
      X[[i]] = current.2
    } else {
      current = X[[i]]
      current[,1] = CenterByMean(current[,1])
      X[[i]] = current
    }
  }
  
  # Adding an intercept
  for (k in 1:n) {
    if (is.null(nrow(X[[k]]))) {
      len_X = length(X[[k]])
      intercept = rep(1, len_X)
      X[[k]] = cbind(intercept, X[[k]])
    } else {
      len_X = nrow(X[[k]])
      intercept = rep(1, len_X)
      X[[k]] = cbind(intercept, X[[k]])
    }
  }
  
  # Extracting the necessary information
  n_vec = sapply(X, nrow) # number of patients per cancer type
  
  # gathering test observation indices
  Test_obs = lapply(seq(length(cancer_types_27)), 
                    function(i) seq(n_vec[i])[!(seq(n_vec[i]) %in% Training_obs[[i]])]) 
  
  # Obtaining test observations and their times of last contact, etc
  Y_test = lapply(seq(length(cancer_types_27)), function(i) Y[[i]][Test_obs[[i]]]) 
  X_test = lapply(seq(length(cancer_types_27)), function(i) X[[i]][Test_obs[[i]], ])
  Last_Contact_test = lapply(seq(length(cancer_types_27)), function(i) Last_Contact[[i]][Test_obs[[i]]])
  
  # combine test data into one list
  List_XY = lapply(seq(length(cancer_types_27)), function(i) cbind(X_test[[i]], Y_test[[i]], 
                                                                   Last_Contact_test[[i]]))
  
  y_col = ncol(List_XY[[1]]) - 1 # the column in List_XY that contains the survival time
  yc_col = ncol(List_XY[[1]]) # the column in List_XY that contains the censor time
  
  posterior.vec = list()
  
  for (j in 1:1000) { 
    
    # for each cancer type
    # for each row in that cancer type
    # apply Likelihood with the corresponding cancer type's posterior parameters
    # at the jth iteration
    posterior.vec[[j]] = unlist(lapply(seq(length(cancer_types_27)), 
                                       function(k) sapply(seq(nrow(List_XY[[k]])), # k is for each cancer type
                                                          function(i) LogExponentialLikelihood(List_XY[[k]][i, -c(y_col, yc_col)], 
                                                                                               List_XY[[k]][i, y_col],  # i is for each patient in the test set
                                                                                               List_XY[[k]][i, yc_col], # ith observation, censor time
                                                                                               betas[[k]][[j]]))))  # kth cancer type, jth iteration posterior
  }
  
  #EFL: organize log-likelihoods into M by N matrix, where M is number of test samples and N is number of iterations
  logLike.mat <- matrix(nrow = length(unlist(Test_obs)), ncol= length(posterior.vec))
  for(j in 1:length(posterior.vec)){
    logLike.mat[,j] <- unlist(posterior.vec[[j]])
  }
  #EFL: compute joint log-likelihoods, at each iteration, by summing the logs
  logLike.joint <- colSums(logLike.mat)
  #EFL: Compute log of full posterior predictive likelihood, using logSum
  PL_Survival<- -log(length(posterior.vec))+logSum(logLike.joint)
  
  
  return(PL_Survival)
}

Last_Contact_27 = Last_Contact

# 5-fold cross validation
cl <- makeCluster(12)
registerDoParallel(cl)
PL.exp = foreach(cv_iter = 1:5, .combine = c, .packages = c("rjags")) %dopar%{
  n_vec = sapply(F27.50.3, nrow)
  F5Training_obs = Generate5FoldTrainingSet(F27.50.3, S27.50.3, cancer_types_27, n_vec) # the training sets
  
  Training_obs_cv_iter = sapply(F5Training_obs[[1]], '[', cv_iter) # cv_iter'th training fold
  
  # subsetting the data for the training fold
  # '_S' indicates that it is a subset of the full data 
  F27.50.3_S = lapply(1:length(F27.50.3), function(k) F27.50.3[[k]][Training_obs_cv_iter[[k]], ]); names(F27.50.3_S) = cancer_types_27 
  S27.50.3_S = lapply(1:length(S27.50.3), function(k) S27.50.3[[k]][Training_obs_cv_iter[[k]]]) 
  Last_Contact_27_S = lapply(1:length(Last_Contact_27), function(k) Last_Contact_27[[k]][Training_obs_cv_iter[[k]]])
  
  # Adding cancer types to data
  names(S27.50.3_S) = cancer_types_27
  names(Last_Contact_27_S) = cancer_types_27
  
  
  # Centering age in the design matrix for each inner matrix
  F27.50.3_S = lapply(F27.50.3_S, function(mat) {
    centered_age = CenterByMean(mat[,1])
    cbind(centered_age, mat[,2:ncol(mat)])
  })
  
  # Creating a list of matrices of censor and survival time combined together
  SLC.27.50_S = lapply(1:length(S27.50.3_S), function(i) {
    tcens = S27.50.3_S[[i]] # the tcens times
    
    # replacing the NAs in the tcens times (the censored observations) with the censor times 
    # (obs with tcens events don't have censor times)
    tcens[is.na(tcens)] = Last_Contact_27_S[[i]][is.na(S27.50.3_S[[i]])] 
    any(is.na(tcens))
    
    # vector for whether the obs is censored or not
    censored = as.numeric(is.na(S27.50.3_S[[i]])) # TRUE if censored, FALSE otherwise
    
    # returning tcens and censored together
    cbind(tcens, censored)
  })
  
  # Adding a column to each inner matrix that indicates which cancer type it is
  # and including the survival/censor time and the indicator for censoring
  F27.50.3_S.extra_col = lapply(1:length(F27.50.3_S), function(i) {
    cancer = rep(i, nrow(F27.50.3_S[[i]]))
    cbind(F27.50.3_S[[i]], cancer, SLC.27.50_S[[i]])
  })
  
  # Concatenating each inner matrix together to make one big matrix
  # Then adding an intercept
  F27.50.3_S.mat = do.call(rbind, F27.50.3_S.extra_col)
  F27.50.3_S.mat = cbind(intercept = 1, F27.50.3_S.mat)
  F27.50.3_S.mat = as.data.frame(F27.50.3_S.mat)
  
  # Setting up the data for JAGS
  tcga.data_S <- with(F27.50.3_S.mat,
                      list("X" = F27.50.3_S.mat[, 1:52],
                           "is.censored" = (censored==1),
                           "tcens" = tcens,
                           "t" = ifelse(censored==1, tcens, NA),
                           "cancer" = as.numeric(cancer), 
                           "shape"=1) #set shape parameter to 1 for Weibull (gives exponential distribution)
  )
  
  # setting up the model
  jags.model.full.cv_iter <- "model {
  for(i in 1:length(t)) { # for each observation
  is.censored[i] ~ dinterval(t[i], tcens[i])
  t[i] ~ dweib(shape, 1/lambda[i])
  log(lambda[i]) <- inprod(X[i,1:52], beta[cancer[i],])
  }
  
  # Priors
  for (i in 1:27) { # for each cancer type
  for (j in 1:52) { # for each covariate including the intercept
  beta[i, j] ~ dnorm(betatilde[j], g[j])
  }
  }
  
  tau2 <- 1/(10000*10000)
  for (j in 1:52) { # for each covariate
  g[j] ~ dgamma(0.01, 0.01)
  betatilde[j] ~ dnorm(0, tau2)
  }
  
}"  
  
  # Load and compile the model
  mod.full.cv_iter <- jags.model(textConnection(jags.model.full.cv_iter),
                                 data = tcga.data_S) 
  
  # Burn-in
  update(mod.full.cv_iter, 10000) 
  
  # Generate 10000 more iterations
  s.cv_iter <- jags.samples(mod.full.cv_iter, c('beta', 'betatilde', 'g'), 
                            n.iter = 10000, thin = 10) 
  
  iters = 10000
  # Save the posteriors for betas in format compatible with posterior likelihood function
  betas = lapply(1:length(F27.50.3_S), function(type) { # for each cancer type
    lapply(1:(iters/10), function(iter) { # for each iteration (there are 1000 because the chain was thinned every 10 samples)
      type_iter = s.cv_iter$beta[type, , iter, 1]
      names(type_iter) = c("intercept", colnames(F27.50.3_S[[1]]))
      type_iter
    })
  })
  names(betas) = cancer_types_27 # name the betas by the cancer type
  
  PosteriorLikelihood(betas = betas, X = F27.50.3, Y = S27.50.3, Last_Contact = Last_Contact_27, Training_obs = Training_obs_cv_iter)
  }
stopCluster(cl)

mean(PL.exp) 

