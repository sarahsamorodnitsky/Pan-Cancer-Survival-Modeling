# Model validation for Gibbs sampler
# Author: Sarah Samorodnitsky

library(truncnorm)
library(MASS)
library(EnvStats)
library(svMisc)

clinical_data = read.csv("TCGA-CDR.csv")
cancer_types = as.character(unique(clinical_data$type))

# Assess coverage rates
# Simulation to evaluate coverage
SimulateCoverage = function(iters) {
  # Use to evaluate coverage of function written above for data
  GenerateData = function(true_betas, true_sigma2, n_vec) {
    # Generate data in the same format as the TCGA data so as to use the function written above
    # survival outcomes include 50% censored observations
    Design = list()
    Response = list()
    Censored = list()
    tot = length(true_betas)
    for (i in 1:tot) { # 33 groups
      n = n_vec[i]
      beta_coefs = as.matrix(true_betas[[i]])
      X = cbind(rep(1,n), matrix(rnorm(n*2, 0, 1), nrow = n, ncol = 2)) # intercept, age, mutation status
      Y = matrix(rnorm(n, mean = X %*% beta_coefs, sd = sqrt(true_sigma2)), nrow = n, ncol = 1) # generating survival times
      C = matrix(rnorm(n, mean = X %*% beta_coefs, sd = sqrt(true_sigma2)), nrow = n, ncol = 1) # generating time of censoring
      
      # if censor time comes before survival time, then subject is censored and we don't know how long they survived until
      num_censored = sum(Y > C)
      num_not_censored = length(Y) - num_censored
      which_to_censor = Y > C # the observations which are censored
      which_to_not_censor = Y <= C # the observations which are not censored

      Y[which_to_censor,] = rep(NA, num_censored) # censor all observations whose survival times exceeded last contact times
      C[which_to_not_censor,] = rep(NA, num_not_censored) # ignore all censored times for observations whose survival times were less than their censor times

      Design[[i]] = X
      Response[[i]] = Y
      Censored[[i]] = C
    }
    return(list(Design = Design, Response = Response, Censored = Censored))
  }
  
  GetLastContact = function(C2) {
    # From the list of censor times (times of censoring and NAs for uncensored obs), gets only the times of last contact for censored observations
    # without the NAs
    n = length(C2)
    CensorTimes = list() # list of times of last contact for observations that are censored
    for (i in 1:n) {
      CensorTimes[[i]] = C2[[i]][which(!is.na(C2[[i]]))] # of censored times (NA for those uncensored) get the ones who are censored (not NAs)
    } 
    return(CensorTimes)
  }
  
  Check_if_in = function(true_param, interval) {
    # checks if the true parameter is contained in the credible interval
    check = true_param >= interval[1] & true_param <= interval[2]
    return(as.numeric(check))
  }
  
  InitManyLists = function(n) {
    new_list = list()
    for (i in 1:n) {
      new_list[[i]] = list()
    }
    return(new_list)
  }
  
  CredibleInterval = function(sorted_parameter, iters) {
    # Returns credible interval for sorted parameter
    # if sorted_parameter is a matrix, then it returns a list of credible intervals (for the betas)
    # otherwise returns a single credible interval, for sigma^2 or lambda^2 or beta_tilde
    cutoff_2.5 = iters*0.025 # cut offs for credible interval
    cutoff_97.5 = iters*0.975 
    if (is.vector(sorted_parameter)) { # if the parameter is sigma2 or beta_tilde or lambda2
      interval = c(sorted_parameter[cutoff_2.5 + 1], 
                   sorted_parameter[cutoff_97.5 - 1]) # cutting off upper and lower 2.5%tile
    } else { # if the parameter is the betas
      interval = list()
      p = ncol(sorted_parameter)
      for (i in 1:p) {
        current_interval = c(sorted_parameter[cutoff_2.5 + 1, i], 
                             sorted_parameter[cutoff_97.5 - 1, i])
        interval[[i]] = current_interval
      }
    }
    return(interval)
  }
  
  PosteriorBetasReformat = function(betas, cancer_types) {
    # lapply(test, function(x) apply(matrix(unlist(x), ncol = 3, byrow = TRUE)[1001:2000,], 2, sort)) another way of doing this w/o a loop
    n = length(betas)
    Reformatted_Betas = list()
    for (i in 1:n) {
      betas_current = apply(matrix(unlist(betas[i]), 
                                   ncol = 3, byrow = TRUE)[1001:2000,], 2, sort)
      colnames(betas_current) = c("Intrcpt", "Age", "Mut")
      Reformatted_Betas[[i]] = betas_current
    }
    names(Reformatted_Betas) = cancer_types
    return(Reformatted_Betas)
  }
  
  PosteriorBetasIntervals = function(betas, cancer_types, iters) {
    Intervals = list()
    n = length(betas)
    for (i in 1:n) {
      current = CredibleInterval(betas[[i]], iters)
      Intervals[[i]] = current
    }
    names(Intervals) = cancer_types
    return(Intervals)
  }
  
  Model_Modified = function(SimulatedFull, SimulatedSurvival, SimulatedCensored, 
                            sigma2_start, beta_tilde_start, 
                            lambda2_start, Survival_Dist) {
    # modified version of original Model to exclude all data rearrangements 
    
    X2 = SimulatedFull
    Y2 = SimulatedSurvival
    C2 = SimulatedCensored
    
    ListToMatrix = function(K) {
      # Useful helper function to transform lists into matrices
      # returns a 3-column matrix
      if (length(K[[1]]) == 1) { # then K = Y, which is a list where each entry is one number
        return(matrix(unlist(K)))
      } else { # then K = X where X is a list of vectors
        matrix = matrix(nrow = length(K), ncol = 3)
        for (i in 1:length(K)) {
          matrix[i,] = K[[i]]
        }
        return(matrix)
      }
    }
    
    # Functions to help calculate posteriors
    Beta_Posterior_Helper = function(x_i, y_i, sigma2, lambda2, current_beta_tilde) {
      # returns the mean and the variance of the multivariate normal posterior for each beta_i
      # the mean should be a vector of length 2 because there are 2 covariates
      # the variance should be a 2x2 covariance matrix because 2 covariates
      B = solve((1/sigma2)*t(x_i)%*%x_i + diag((1/lambda2))) # not sure if its right to diag(1/lambda2)
      b = (1/sigma2)*t(x_i) %*% y_i + (1/lambda2) * current_beta_tilde
      mean = B%*%b; var = B
      return(list(mean = mean, var = var))
    }
    
    Lambda2_Posterior_Helper = function(current_betas, current_beta_tilde, I) {
      tot = c(0, 0, 0)
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
        tot = tot + sum((y_i - x_i %*% current_betas_i)^2)
      }
      return(tot)
    }
    
    Beta_Tilde_Posterior_Helper = function(current_betas_mean, current_lambda2) {
      post_beta_tilde_mean = (I*tau2*current_betas_mean)/(current_lambda2 + I*tau2)
      post_beta_tilde_var = (current_lambda2*tau2)/(current_lambda2 + I*tau2)
      return(list(mean = post_beta_tilde_mean, var = post_beta_tilde_var))
    }
    
    n_vec =  c(unlist(lapply(X2, nrow)))
    I = length(n_vec) # the number of cancer types
    p = 3 # the number of covariates: age and TP53 mutation status + the intercept
    
    # creating other important quantities
    Identity = diag(p)
    tau2 = 1
    total_obs = sum(n_vec)
    n = length(X2)
    
    num_censored = sapply(sapply(Y2, is.na), sum) # the number of censored observations in each group
    which_censored = sapply(Y2, is.na) # True/false list
    last_contact = GetLastContact(C2) # gets the censored times
    
    # initializing the objects to store the posterior samples
    betas = InitManyLists(I) # list to store posterior betas, four inner lists because four cancer types
    names(betas) = names(X2)
    beta_tilde = matrix(ncol = 3, nrow = 2001) # list to include posterior mean for age covariate and mutation status covariate
    beta_tilde[1,] = beta_tilde_start
    colnames(beta_tilde) = c("Intercept","Age", "Mut_Status")
    lambda2 = matrix(ncol = 3, nrow = 2001)
    lambda2[1,] = lambda2_start
    colnames(lambda2) = colnames(beta_tilde)
    sigma2 = c()
    sigma2[1] = sigma2_start
    
    # initialize censored observations to have their last_contact time be their survival time
    # within loop, generate their survival times randomly
    for (k in 1:n) {
      Y2[[k]][which_censored[[k]]] = last_contact[[k]]
    }
    
    for (i in 1:2000) {
      # posterior for the betas
      for (j in 1:I) { # I is the number of groups/cancer types
        beta_post_params = Beta_Posterior_Helper(X2[[j]], Y2[[j]], sigma2[i], lambda2[i,], beta_tilde[i,])
        betas[[j]][[i]] = mvrnorm(1, mu = beta_post_params$mean, Sigma = beta_post_params$var) # posterior beta vector for jth group
      }
      
      # posterior for lambda^2
      current_betas = sapply(betas, `[`, i) # retrieve the current betas for the current group
      current_beta_tilde = beta_tilde[i,]
      W = Lambda2_Posterior_Helper(current_betas, current_beta_tilde, I)
      post_lambda2_shape = (I)/2 + 1
      post_lambda2_rate = 1 + (0.5)*W
      lambda2[i+1, 1] = 1/rgamma(1, shape = post_lambda2_shape, rate = post_lambda2_rate[1]) # for the intercept
      lambda2[i+1, 2] = 1/rgamma(1, shape = post_lambda2_shape, rate = post_lambda2_rate[2]) # for age
      lambda2[i+1, 3] = 1/rgamma(1, shape = post_lambda2_shape, rate = post_lambda2_rate[3]) # for TP53 mutation status
      
      # posterior for beta tilde
      current_betas_mean = apply(ListToMatrix(current_betas), 2, mean) # mean(ListToMatrix(current_betas))
      current_lambda2 = lambda2[i,]
      post_beta_tilde = Beta_Tilde_Posterior_Helper(current_betas_mean, current_lambda2)
      beta_tilde[i+1, 1] = rnorm(1, post_beta_tilde$mean[1], sqrt(post_beta_tilde$var[1]))
      beta_tilde[i+1, 2] = rnorm(1, post_beta_tilde$mean[2], sqrt(post_beta_tilde$var[2]))
      beta_tilde[i+1, 3] = rnorm(1, post_beta_tilde$mean[3], sqrt(post_beta_tilde$var[3]))
      
      # posterior for sigma^2
      post_sigma2_shape = (total_obs/2) + 1
      post_sigma2_rate = 0.5*Sigma2_Posterior_Rate_Helper(X2, Y2, current_betas, I) + 1
      sigma2[i+1] = 1/rgamma(1, post_sigma2_shape, rate = post_sigma2_rate)
      
      # Generate Ys for censored observations
      for (k in 1:n) { # n needs to be the number of cancer types
        n_gens = num_censored[[k]]
        Censored_Obs = X2[[k]][which_censored[[k]],]
        Censored_Last_Contact = last_contact[[k]]
        Censor_Lower_Bound = Censored_Last_Contact
        if (Survival_Dist == "normal") {
          if (nrow(Censored_Obs) != 0) { # only do this if there are censored observations
            Mu_Survival = current_betas[[k]][1] +
              current_betas[[k]][2]*Censored_Obs[,2] +
              current_betas[[k]][3]*Censored_Obs[,3]
            random_survival = rtruncnorm(n_gens, a = Censor_Lower_Bound,
                                         mean = Mu_Survival, sd = rep(sqrt(sigma2[i]), n_gens))
          }
        }
        Y2[[k]][which_censored[[k]]] = random_survival
      }
    }
    return(list(betas = betas, beta_tilde = beta_tilde, lambda2 = lambda2, sigma2 = sigma2))
  }
  
  n_params = 106 # the total number of estimate parameters (ie all the sigmas, lambdas, betas, beta_tildes)
  Coverage = rep(0, n_params) # proportion for each parameter
  
  for (i in 1:iters) {
    progress(i/(iters/100))
    
    true_beta_tilde = rnorm(3, 0, 1) # generate from the prior
    true_sigma2 = 1/rgamma(1, 1, 1) # generate from the prior
    true_lambda2 = 1/rgamma(3, 1, 1) # generate from the prior
    true_betas_matrix = mvrnorm(length(n_vec), true_beta_tilde, diag(true_lambda2));
    true_betas = split(true_betas_matrix, row(true_betas_matrix))
    SimulatedData = GenerateData(true_betas, true_sigma2, n_vec)
    SimulatedFull = SimulatedData$Design
    SimulatedSurvival = SimulatedData$Response
    SimulatedCensored = SimulatedData$Censored
    
    # Generating starting values
    sigma2_start = 1
    lambda2_start = rep(1, 3)
    beta_tilde_start = rnorm(3, 0, 10)
    Survival_Dist = "normal"
    
    Posteriors = Model_Modified(SimulatedFull, SimulatedSurvival, SimulatedCensored,
                                sigma2_start = sigma2_start,
                                beta_tilde_start = beta_tilde_start,
                                lambda2_start = lambda2_start,
                                "normal")
    
    beta_tilde = apply(Posteriors$beta_tilde[1001:2000,], 2, sort)
    sigma2 = sort(Posteriors$sigma2[1001:2000])
    lambda2 = apply(Posteriors$lambda2[1001:2000,], 2, sort)
    betas = PosteriorBetasReformat(Posteriors$betas, cancer_types)
    
    # Creating credible intervals
    beta_tilde_interval = CredibleInterval(beta_tilde, iters)
    sigma2_interval = CredibleInterval(sigma2, iters)
    lambda2_interval = CredibleInterval(lambda2, iters)
    beta_intervals = PosteriorBetasIntervals(betas, cancer_types, iters)
    
    Intervals = matrix(unlist(list(beta_tilde_interval, sigma2_interval, 
                                   lambda2_interval, beta_intervals)), ncol = 2, byrow = T)
    TrueValues = unlist(list(true_beta_tilde, true_sigma2,
                             true_lambda2, unlist(unname(true_betas))))
    
    for (j in 1:n_params) {
      In = Check_if_in(TrueValues[j], Intervals[j,])
      Coverage[j] = Coverage[j] + In
    }
  }
  names(Coverage) = c("beta_tilde_intercept", "beta_tilde_age", "beta_tilde_mutation",
                      "sigma2", "lambda2_intercept", "lambda2_age", "lambda2_mutation",
                      unlist(lapply(cancer_types, paste, c("intercept", "age", "mutation"))))
  Coverage = Coverage/iters
  return(Coverage)
}

n_vec = sample(20:500, size = 33)
TestCoverage = SimulateCoverage(1000)


