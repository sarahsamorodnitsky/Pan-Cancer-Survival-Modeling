# Forward selection for iteratively adding genes to model and determing which 
# maximizes the log-posterior likelihood. 
# Author: Sarah Samorodnitsky

load("FSCG.rda")

library(MASS)
library(truncnorm)
library(EnvStats)
library(doParallel)
library(foreach)

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

# The model which returns the set of training observations
ModelOnTrainingData = function(Full_Input_C, Survival_Input_C, Last_Contact, sigma2_start, genes,
                               beta_tilde_start, lambda2_start, iters, p, cv_iter, F5Training_Obs) {
  # Survival_Dist = "normal" if survival times generated from truncated normal dist or
  # "log-normal" if survival times generated from truncated log-normal distribution
  # p = the number of covariates: age + 55 gene mutation status + the intercept (user-supplied)
  
  ListToMatrix = function(K) {
    # Useful helper function to transform lists into matrices
    # returns a 3-column matrix
    if (length(K[[1]]) == 1) { # then K = Y, which is a list where each entry is one number
      return(matrix(unlist(K)))
    } else { # then K = X where X is a list of vectors
      matrix = matrix(nrow = length(K), ncol = p)
      for (i in 1:length(K)) {
        matrix[i,] = K[[i]]
      }
      return(matrix)
    }
  }
  
  CenterByMean = function(vec) {
    new_vec = vec - mean(na.omit(vec))
    return(new_vec)
  }
  
  InitManyLists = function(n) {
    new_list = list()
    for (i in 1:n) {
      new_list[[i]] = list()
    }
    return(new_list)
  }
  
  # For testing
  # x_i = X[[j]]; y_i = Y[[j]]
  # sigma2_i = sigma2[i]
  # lambda2_i = lambda2[i,]
  # current_beta_tilde = beta_tilde[i,]
  
  # Functions to help calculate posteriors
  Beta_Posterior_Helper = function(x_i, y_i, sigma2_i, lambda2_i, current_beta_tilde) {
    # returns the mean and the variance of the multivariate normal posterior for each beta_i
    # the mean should be a vector of length 2 because there are 2 covariates
    # the variance should be a 2x2 covariance matrix because 2 covariates
    x_i = as.matrix(x_i); y_i = as.matrix(y_i)
    B = solve((1/sigma2_i)*t(x_i)%*%x_i + diag((1/lambda2_i))) # not sure if its right to diag(1/lambda2)
    b = (1/sigma2_i)*t(x_i) %*% y_i + (1/lambda2_i) * current_beta_tilde
    mean = B%*%b; var = B
    return(list(mean = mean, var = var))
  }
  
  Lambda2_Posterior_Helper = function(current_betas, current_beta_tilde, I) {
    # calculates posterior variance for each coefficent across all 4 cancer types
    tot = rep(0, p)
    for (i in 1:I) {
      tot = tot + (current_betas[[i]] - current_beta_tilde)^2
    }
    return(tot)
  }
  
  Sigma2_Posterior_Rate_Helper = function(X, Y, current_betas, I) {
    # returns the rate for the posterior of sigma2
    tot = 0
    for (i in 1:I) {
      x_i = X[[i]]; y_i = Y[[i]]; current_betas_i = current_betas[[i]]
      x_i = as.matrix(x_i); y_i = as.matrix(y_i)
      tot = tot + sum((y_i - x_i %*% current_betas_i)^2)
      # if (is.na(tot)) {print(i)}
    }
    return(tot)
  }
  
  Beta_Tilde_Posterior_Helper = function(current_betas_mean, current_lambda2) {
    post_beta_tilde_mean = (I*tau2*current_betas_mean)/(current_lambda2 + I*tau2)
    post_beta_tilde_var = (current_lambda2*tau2)/(current_lambda2 + I*tau2)
    return(list(mean = post_beta_tilde_mean, var = post_beta_tilde_var))
  }
  
  X = Full_Input_C
  Y = Survival_Input_C
  
  n = length(X)
  # Centering age
  for (i in 1:n) {
    current = X[[i]]
    current[,1] = CenterByMean(current[,1])
    X[[i]] = current
  }
  
  # Adding an intercept
  for (k in 1:n) {
    len_X = nrow(X[[k]])
    intercept = rep(1, len_X)
    X[[k]] = cbind(intercept, X[[k]])
  }
  
  # Determining which observations in each cancer type are censored 
  Censored = list()
  for (j in 1:n) {
    current_rows = which(is.na(Y[[j]]))
    Censored[[j]] = current_rows
  }
  
  n_vec =  c(unlist(lapply(X, nrow))) # the number of observations per group
  I = length(n_vec) # the number of cancer types
  
  # creating other important quantities
  Identity = diag(p)
  tau2 = 10000^2
  total_obs = sum(n_vec)
  
  ### initializing the objects to store the posterior samples
  ### first initialize for log-normal posterior samples
  betas = InitManyLists(n) # list to store posterior betas, four inner lists because four cancer types
  names(betas) = names(X)
  beta_tilde = matrix(ncol = p, nrow = iters + 1) # list to include posterior mean for age covariate and mutation status covariate
  beta_tilde[1,] = beta_tilde_start
  lambda2 = matrix(ncol = p, nrow = iters + 1)
  lambda2[1,] = lambda2_start
  colnames(lambda2) = colnames(beta_tilde)
  sigma2 = c()
  sigma2[1] = sigma2_start
  
  if (ncol(X[[1]]) == 2) { # just intercept and age
    colnames(beta_tilde) = c("Intercept","Age")
  } else {
    colnames(beta_tilde) = c("Intercept","Age", paste("Mut Status", genes))
  }
  
  # Creating a training set of observations
  current_train_ind = sapply(F5Training_obs[[1]], '[', cv_iter)
  X_S = lapply(1:length(X), function(k) X[[k]][current_train_ind[[k]], ]); names(X_S) = cancer_types_27
  Y_S = lapply(1:length(Y), function(k) Y[[k]][current_train_ind[[k]]])
  Last_Contact_S = lapply(1:length(Last_Contact), function(k) Last_Contact[[k]][current_train_ind[[k]]])
  Censored_S = lapply(1:length(Y), function(k) which(is.na(Y_S[[k]])))
  
  total_obs2 = sum(sapply(Y_S, length)) # calculating the total number of observations after subsetting
  
  # initializing Y so that each censored observation is replaced with last contact time
  for (k in 1:n) {
    Y_S[[k]][Censored_S[[k]]] = Last_Contact_S[[k]][Censored_S[[k]]]
  }
  
  # taking the log of the survival times for Y1
  for (k in 1:n) {
    Y_S[[k]] = log(Y_S[[k]])
  }
  
  for (i in 1:iters) {
    
    # posterior for the betas
    for (j in 1:I) { # I is the number of groups/cancer types
      beta_post_params = Beta_Posterior_Helper(X_S[[j]], Y_S[[j]], sigma2[i], lambda2[i,], beta_tilde[i,])
      betas[[j]][[i]] = mvrnorm(1, mu = beta_post_params$mean, Sigma = beta_post_params$var) # posterior beta vector for jth group
    }
    
    # posterior for lambda^2
    current_betas = sapply(betas, `[`, i) # retrieve the current betas for the current iteration
    current_beta_tilde = beta_tilde[i,]
    W = Lambda2_Posterior_Helper(current_betas, current_beta_tilde, I)
    post_lambda2_shape = (I/2) + 0.01 # change 1 to 0.01
    post_lambda2_rate = 0.01 + (0.5)*W # change 1 to 0.01
    lambda2[i+1, ] = 1/rgamma(p, shape = post_lambda2_shape, rate = post_lambda2_rate) 
    
    # posterior for beta tilde
    current_betas_mean = apply(ListToMatrix(current_betas), 2, mean)
    current_lambda2 = lambda2[i,]
    post_beta_tilde = Beta_Tilde_Posterior_Helper(current_betas_mean, current_lambda2)
    beta_tilde[i+1, ] = rnorm(p, post_beta_tilde$mean, sqrt(post_beta_tilde$var))
    
    
    # posterior for sigma^2
    post_sigma2_shape = (total_obs2/2) + 0.01 # change 1 to 0.01
    post_sigma2_rate = 0.5*Sigma2_Posterior_Rate_Helper(X_S, Y_S, current_betas, I) + 0.01 # change 1 to 0.01
    sigma2[i+1] = 1/rgamma(1, post_sigma2_shape, rate = post_sigma2_rate)
    
    
    # Generate Ys for censored observations
    for (k in 1:n) {
      n_gens = length(Censored_S[[k]])
      Censored_Obs = X_S[[k]][Censored_S[[k]],]
      Censored_Last_Contact = Last_Contact_S[[k]][Censored_S[[k]]]
      Censor_Lower_Bound = Last_Contact_S[[k]][Censored_S[[k]]]

      Mu_Survival = Censored_Obs %*% t(t(current_betas[[k]]))
      random_survival = rtruncnorm(n_gens, a = log(abs(Censor_Lower_Bound)), 
                                   mean = Mu_Survival, 
                                   sd = rep(sqrt(sigma2[i]), n_gens))
      Y_S[[k]][Censored_S[[k]]] = random_survival
    }
  }
  
  
  return(list(betas = betas, beta_tilde = beta_tilde, lambda2 = lambda2, sigma2 = sigma2, 
              current_train_ind = current_train_ind))
}

# To subset out a specific set of genes
SubsetNGenes = function(F27.50.2, Genes) {
  # Subsets entire dataset to just a selected set of genes
  n = length(F27.50.2)
  F27.N.2 = F27.50.2 # note the difference here - N instead of 50
  for (i in 1:n) {
    current_type = F27.50.2[[i]]
    current_genes = colnames(current_type)
    current_available = current_genes %in% Genes
    current_subset = current_type[, current_available]
    current_subset = cbind(current_type[,1:2], current_subset)
    colnames(current_subset) = c("Age", "LastContact", current_genes[current_genes %in% Genes])
    F27.N.2[[i]] = current_subset
  }
  return(F27.N.2)
}

# Function for calculating the posterior likelihood
PosteriorLikelihood = function(posteriors, X, Y) {
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
  
  LogLikelihood = function(x_ij, y_ij, y_ijc, betas_it, sigma2_it) {
    # calculate the likelihood for just one set of posteriors
    mu_ij = sum(x_ij * betas_it) 
    
    LogNormalLikelihood = function(y, mu, sig2) {
      # Calculates log density of normal 
      # nl = (1/sqrt(2*pi*sig2))*exp(-(1/(2*sig2))*(y - mu)^2)
      lnl = (-0.5)*log(2*pi*sig2) - (1/(2*sig2))*(y - mu)^2
      return(lnl)
    }
    
    if (!is.na(y_ij)) { # y_ij is not censored
      post_ij = LogNormalLikelihood(log(y_ij), mu = mu_ij, sig2 = sigma2_it) # dnorm(log(y_ij), mean = mu_ij1, sd = sqrt(sigma21_it)) # sample from log-normal
      
    } else { # if y_ij is censored, use censoring time
      post_ij = log(1 - pnorm((log(y_ijc) - mu_ij)/sqrt(sigma2_it), mean = 0, sd = 1)) 
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
  Training_obs_cv_iter = posteriors$current_train_ind
  betas = posteriors$betas
  sigma2 = posteriors$sigma2
  Test_obs = lapply(seq(length(cancer_types_27)), function(i) seq(n_vec[i])[!(seq(n_vec[i]) %in% Training_obs_cv_iter[[i]])]) # gathering test observation indices
  
  
  # Obtaining test observations and their times of last contact, etc
  Y_test = lapply(seq(length(cancer_types_27)), function(i) Y[[i]][Test_obs[[i]]]) 
  X_test = lapply(seq(length(cancer_types_27)), function(i) X[[i]][Test_obs[[i]], ])
  Last_Contact_test = lapply(seq(length(cancer_types_27)), function(i) Last_Contact[[i]][Test_obs[[i]]])
  # combine test data into one list
  List_XY = lapply(seq(length(cancer_types_27)), function(i) cbind(X_test[[i]], Y_test[[i]], Last_Contact_test[[i]]))
  
  y_col = ncol(List_XY[[1]]) - 1 # the column in List_XY that contains the survival time
  yc_col = ncol(List_XY[[1]]) # the column in List_XY that contains the censor time
  
  posterior.vec = list()
  ind_to_use = seq(10000, 20000, by = 10)
  
  for (j in ind_to_use) { ##EFL: changed index 'i' to 'j' here
    
    # for each cancer type
    # for each row in that cancer type
    # apply Likelihood with the corresponding cancer type's posterior parameters
    # at the jth iteration
    # takes 7 seconds to run
    #EFL: also changed here 'i' to 'j' on left-hand side below
    posterior.vec[[j]] = unlist(lapply(seq(length(cancer_types_27)), function(k) sapply(seq(nrow(List_XY[[k]])), # for each cancer type
                                                                                            function(i) LogLikelihood(List_XY[[k]][i, -c(y_col, yc_col)], List_XY[[k]][i, y_col],  # for each patient in the test set
                                                                                                                      List_XY[[k]][i, yc_col], # ith observation, censor time
                                                                                                                      betas[[k]][[j]],  # kth cancer type, jth iteration posterior
                                                                                                                      sigma2[j]))))
    
  }
  
  # calculating the posterior likelihood
  #EFL: I've commented out the code below  -- see my email
  #PL = 0
  #for (i in 1:length(posterior.vec)) {
  #  iter_i = posterior.vec[[i]] # the posterior likelihood values for the ith iteration
  ##EFL: Added line below to define iter_i.vec
  #  iter_i.vec = unlist(iter_i)
  
  # log-likelihood: sum of the logs
  # ll_i = sum(sapply(iter_i, function(j) sum(log(j))))
  #  ll_i = logSum(iter_i.vec)
  
  #  PL = PL + ll_i
  #}
  # PL_Survival = PL/length(posterior.vec)
  
  #EFL: here is revised code to aggrate the log-likelihoods
  #EFL: organize log-likelihoods into M by N matrix, where M is number of test samples and N is number of iterations
  logLike.mat <- matrix(nrow = length(unlist(Test_obs)), ncol= length(ind_to_use))
  for (j in ind_to_use) { # corresponds to every 10th sample 
    k = which(j == ind_to_use) # need to explicitly state the column number because the values of j do not correspond
    logLike.mat[,k] <- unlist(posterior.vec[[j]])
  }
  #EFL: compute joint log-likelihoods, at each iteration, by summing the logs
  logLike.joint <- colSums(logLike.mat)
  #EFL: Compute log of full posterior predictive likelihood, using logSum
  PL_Survival<- -log(length(posterior.vec))+logSum(logLike.joint)
  
  
  return(PL_Survival)
}

Genes = NewGenes

# To encourage continuous check-ups on the forward selection procedure, 
# and to avoid the hassle of the required importation when using parallel loops
# in R, each step of the forward selection procedure was done in this way: 

# Step one: no genes in the model
### Copy and paste from here...
CurrentCovariates = NULL ### Change this everytime a new gene is added
NumberGenes = length(CurrentCovariates) # number genes currently in model
PossibleNewCovariates = Genes[!(Genes %in% CurrentCovariates)]
n_try = length(PossibleNewCovariates)

p = 2 + NumberGenes + 1 # the number of betas in the model: intercept + age = 2, + number of genes already in model, + 1 for possible new gene
beta_tilde_start = rep(0, p)
sigma2_start = .01
lambda2_start = rep(1, p)
iters = 20000

start <- Sys.time()
cl <- makeCluster(12)
registerDoParallel(cl)
PL.p = foreach(g = 1:n_try, .packages = c("MASS", "truncnorm", "EnvStats"),
               .export = c("SubsetNGenes", "PreProcessing", "Generate5FoldTrainingSet", "ModelOnTrainingData",
                           "PosteriorLikelihood")) %dopar% {
    gene = PossibleNewCovariates[g]
    vec = c()
    for (cv_iter in 1:5) {
      GenesToTry = sort(c(CurrentCovariates, gene))
      F27.p.3 = SubsetNGenes(F27.50.3, Genes = GenesToTry)
      S27.p.3 = S27.50.3; Last_Contact.p = Last_Contact
      
      # Generate the sets of training obs
      n_vec =  c(unlist(lapply(F27.p.3, nrow)))
      F5Training_obs = Generate5FoldTrainingSet(F27.p.3, S27.p.3, cancer_types_27, n_vec)
      Model.p = ModelOnTrainingData(F27.p.3, S27.p.3, Last_Contact.p, sigma2_start,
                                    genes = GenesToTry, beta_tilde_start, lambda2_start, iters, p, cv_iter, F5Training_obs)
      test = PosteriorLikelihood(Model.p, F27.p.3, S27.p.3) # current posterior likelihood
      vec[cv_iter] = test
    }
    vec
}
stopCluster(cl)
end <- Sys.time()
end - start
  
names(PL.p) = sapply(PossibleNewCovariates, function(i) paste(i, CurrentCovariates))
PL.ave = sapply(PL.p, mean)
PL.ave[PL.ave == max(PL.ave)]

### ... to here each time you want to add another gene to the model. 

# Forward selection without TP53
GenesNoTP53 = Genes[!(Genes %in% "TP53")]
Genes = GenesNoTP53

### Copy and paste from here...
CurrentCovariates = NULL ### Change this every time a gene is added
NumberGenes = length(CurrentCovariates) # number genes currently in model
PossibleNewCovariates = Genes[!(Genes %in% CurrentCovariates)]
n_try = length(PossibleNewCovariates)

p = 2 + NumberGenes + 1 # the number of betas in the model: intercept + age = 2, + number of genes already in model, + 1 for possible new gene
beta_tilde_start = rep(0, p)
sigma2_start = 1
lambda2_start = rep(1, p)
iters = 20000

start <- Sys.time()
cl <- makeCluster(12)
registerDoParallel(cl)
PL.p = foreach(g = 1:n_try, .packages = c("MASS", "truncnorm", "EnvStats"),
                .export = c("SubsetNGenes", "PreProcessing", "Generate5FoldTrainingSet", "ModelOnTrainingData",
                            "PosteriorLikelihood")) %dopar% {
                              gene = PossibleNewCovariates[g]
                              vec = c()
                              for (cv_iter in 1:5) {
                                F27.p.2 = SubsetNGenes(F27.50.2, Genes = GenesToTry)
                                FS27.p = PreProcessing(F27.p.2, S27.50.2, cancer_types_27)
                                F27.p.3 = FS27.p$Full; S27.p.3 = FS27.p$Survival; Last_Contact.p = FS27.p$LastContact
                                
                                # Generate the sets of training obs
                                n_vec =  c(unlist(lapply(F27.p.3, nrow)))
                                F5Training_obs = Generate5FoldTrainingSet(F27.p.3, S27.p.3, cancer_types_27, n_vec)
                                Model.p = ModelOnTrainingData(F27.p.3, S27.p.3, Last_Contact.p, sigma2_start,
                                                              genes = GenesToTry, beta_tilde_start, lambda2_start, iters, p, cv_iter, F5Training_obs)
                                test = PosteriorLikelihood(Model.p, F27.p.3, S27.p.3) # current posterior likelihood
                                vec[cv_iter] = test
                              }
                              vec
                            }
stopCluster(cl)
end <- Sys.time()
end - start

names(PL.p) = sapply(PossibleNewCovariates, function(i) paste(i, CurrentCovariates))
PL.ave = sapply(PL.p, mean)
PL.ave[PL.ave == max(PL.ave)]
### ... to here each time you want to add another gene to the model. 



