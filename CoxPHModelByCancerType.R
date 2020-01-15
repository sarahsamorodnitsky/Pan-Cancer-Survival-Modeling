# Setting up a Cox Proportional Hazards model for each cancer type to see 
# if it works. 
# Author: Sarah Samorodnitsky & Eric Lock

load("FSCG.rda")

library(survival)

# to store models in 
list_of_models = list()

# for each cancer type
for (i in 1:length(cancer_types_27)) {
  
  # save the survival data for that cancer type
  survival_censor_combined = S27.50.3[[i]]
  
  # replace censored times with times-of-last-contact
  survival_censor_combined[is.na(survival_censor_combined)] = Last_Contact[[i]][!is.na(Last_Contact[[i]])]
  
  # store whether or not an observation was censored
  event = as.numeric(!is.na(S27.50.3[[i]])) # the observations who are dead (=1) or alive (=0) (if NA -> censored and alive, if !NA -> not censored and dead)

  vars_in_list = lapply(seq(ncol(F27.50.3[[i]])), function(col) F27.50.3[[i]][, col]) # separating all the columns into entries in a list
  names(vars_in_list) = colnames(F27.50.3[[i]])
  
  # combine all the data together
  data_i = c(list(survival_censor_combined, event), vars_in_list)
  names(data_i)[1:2] = c("Time", "Event")
  
  # write up formula for model function call
  f_i = as.formula(paste("Surv(Time, Event) ~ ", paste(names(data_i)[-c(1,2)], collapse= "+")))
  mod_i = coxph(f_i, data = data_i)
  
  list_of_models[[i]] = mod_i
}

