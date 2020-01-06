# Code for importing data and calculating mutation rates across cancer types
# Using TCGA2STAT library, imports data through Broad Firehose pipeline in ready-to-use
# format for statistical analysis. 
# Author: Sarah Samorodnitsky

library(TCGA2STAT)

load("CancerDataForFindingGenes.rda") # Use this if having trouble with TCGA2STAT package
clinical_data = read.csv("TCGA-CDR.csv", header = T) # loads in TCGA clinical data
cancer_types = levels(as.factor(clinical_data$type)) # cancer types available in TCGA dataset

# For loading in the data through the TCGA2STAT pipeline
GetTCGAData = function(clinical_data, cancer_types) {
  n = length(cancer_types)
  mutation_list = list()
  for (i in 1:n) {
    current_type = cancer_types[i]
    mut.current = getTCGA(disease = current_type, data.type="Mutation", type="somatic")$dat
    mutation_list[[i]] = mut.current
  }
  names(mutation_list) = cancer_types
  return(mutation_list)
}

# For matching the bar codes of the patients in the TCGA clinical dataset and those imported through TCGA2STAT
MatchBarCodes = function(CancerData, clinical_data, cancer_types) {
  # Returns the mutation data so that the patient barcodes match the order they are in the clinical data set
  # CancerData_S = list of matrices; mutation data 
  # clinical_data = dataframe; survival, censor time, age data
  # cancer_types = character vector; contains names of all the cancer types (should be of length 28)
  
  # For some cancer types, there will be more survival data than there is available mutation data
  # In that case, need to remove observations for which we have survival data but don't have mutation data.
  # Then need to remember to remove the observations for which we don't have mutation data, which is why we update Y
  
  n = length(CancerData)
  X_new = list() 
  Y_new = list() 
  
  for (i in 1:n) {
    current_type = cancer_types[i]
    current_Y = clinical_data[clinical_data$type == current_type, ] # subset the clinical data
    current_X = CancerData[[i]] # subset mutation data
    ord_ind = match(current_Y$bcr_patient_barcode, colnames(current_X)) # the indices to reorder the mutation data so the bar codes match the bar codes in the survival data, this equals the number of observations in the survival data 
    new_current_X = current_X[, na.omit(ord_ind)] # reordering the columns of X (each column is a patient)
    new_current_Y = current_Y[!is.na(ord_ind), ]
    
    if (i == 19) { # for MESO 
      new_current_Y = current_Y
    }
    
    X_new[[i]] = new_current_X
    Y_new[[i]] = new_current_Y
  }
  
  names(X_new) = cancer_types
  names(Y_new) = cancer_types
  return(list(X_new = X_new, Y_new = Y_new))
}

# Calculating the rate of mutation by taking the average mutation rate for each gene across all 27 cancer types being considered
MutationRateByCancer = function(CancerData, cancer_types) {
  n = length(CancerData)
  
  # first collect all the genes we have data on
  all_genes = c()
  for (i in 1:n) {
    type = CancerData[[i]]
    genes.i = rownames(type) # the genes for which data is available for the ith cancer type
    
    genes_to_add = genes.i[!(genes.i %in% all_genes)] # genes in current cancer type which have not been saved yet
    all_genes = c(all_genes, genes_to_add)
  }
  
  # put those genes into a dataframe
  all_genes = sort(all_genes)
  rates = data.frame(gene = all_genes)
  
  # now go through each cancer type and find the rate of mutation for that gene by cancer type
  # add a column to the dataset for that cancer type
  # for any gene that a cancer type has no data on, add a zero
  # make sure the names of the genes for each cancer type align with the order in the dataframe rates
  for (i in 1:n) {
      type = CancerData[[i]]
      rates.i = apply(type, 1, function(row) sum(row)/length(row))
      genes.i = rownames(type) # current genes
      
      missing_genes.i = rates$gene[!(rates$gene %in% genes.i)] # the genes for which there is no data in this cancer type
      rates_to_add.i = rep(0, length(missing_genes.i))
      names(rates_to_add.i) = missing_genes.i
      
      # combining the data for available genes with the added genes
      rates.i = c(rates.i, rates_to_add.i) # update rates vector so it is the same length as the rates dataframe
      rates.i.ord = rates.i[order(factor(names(rates.i), levels = as.character(rates$gene)))]
      rates = cbind(rates, rates.i.ord)
  }
  
  rownames(rates) = rates$gene
  colnames(rates) = c("gene", cancer_types) # naming the columns appropriately
  rates = rates[,-1]
  
  # now I have a matrix where each row is a gene and each column is a cancer type and
  # each entry is the mutation rate for that cancer type and that gene
  
  # subset the rates only for the cancer types of interest
  # rates = rates[,!(colnames(rates) %in% c("PCPG", "PRAD", "TGCT", "THCA", "THYM"))]
  
  # now take the average mutation rate for each gene for each cancer type
  avg_rates = apply(rates, 1, mean)
  return(avg_rates)
}

# Select the genes of interest from the reordered datasets
SelectGenesFromOrigData = function(CancerData_RO_S, clinical_data_list_S, NewGenes, cancer_types_27){ 
  # Using the reordered data, subset the genes to be considered 
  
  n = length(CancerData_RO_S)
  CancerData_RO_S2 = list()
  
  for (i in 1:n) {
    current_type = CancerData_RO_S[[i]]
    avail_genes = rownames(current_type) %in% NewGenes # the genes for which there is data for this cancer type
    current_type_RO = current_type[avail_genes, ]
    
    num_missing = length(NewGenes) - sum(avail_genes)
    if (num_missing > 0) { # then we need to add rows of zeros for the missing genes
      which_missing = NewGenes[!(NewGenes %in% rownames(current_type_RO))]
      mat_to_add = matrix(0, nrow = num_missing, ncol = ncol(current_type_RO))
      rownames(mat_to_add) = which_missing
      current_type_RO = rbind(current_type_RO, mat_to_add)
    }
    
    # Now reorder the rows so the genes are in alphabetical order
    ord_ind = order(rownames(current_type_RO))
    current_type_RO = current_type_RO[ord_ind, ]
    current_type_RO_t = t(current_type_RO) # transpose the matrix so the columns are the genes and the rows are the patient bar codes
    CancerData_RO_S2[[i]] = current_type_RO_t
  }
  
  names(CancerData_RO_S2) = cancer_types_27
  return(CancerData_RO_S2)
}

# Adds age and last contact column to the list of mutation data
# Will give warning because it is changing the way NAs are inputted in the dataset to NAs in R
AddAgeCensorToX = function(CancerData_SelectGenes, clinical_data_list_S, cancer_types_28, NewGenes) {
  n = length(CancerData_SelectGenes)
  
  X = list()
  Y = list()
  
  for (i in 1:n) {
    current_X = CancerData_SelectGenes[[i]]
    current_Y = clinical_data_list_S[[i]]
    age = as.numeric(as.character(current_Y$age_at_initial_pathologic_diagnosis))
    
    current_X_new = cbind(age, as.numeric(as.character(current_Y$last_contact_days_to)), current_X) # concatenating age, censor time, and mutation status
    colnames(current_X_new) = c("Age", "LastContact", NewGenes)
    # rownames(current_X_new) = NULL
    
    current_Y_new = as.numeric(as.character(current_Y$death_days_to))
    names(current_Y_new) = current_Y$bcr_patient_barcode
    
    X[[i]] = current_X_new
    Y[[i]] = current_Y_new
  }
  
  names(X) = cancer_types_28
  names(Y) = cancer_types_28
  
  return(list(X = X, Y = Y))
}

# Removes extraneous variables from the mutation, survival, and censor data i.e.
# any variables missing both a censor time & survival time, ones with a survival time
# that is <= 0, or have a last contact time <= 0
PreProcessing = function(X, Y, cancer_types) {
  # Separates the censoring times from the rest of the data matrix
  # Putting the data in usable form, saving the last contact column separately
  # Returning data matrix, survival data, and last contact vector
  
  n = length(X)
  Last_Contact = list()
  
  for (i in 1:n) {
    current = X[[i]]
    lc = current[,2]
    Last_Contact[[i]] = lc
    current = current[,-2]
    X[[i]] = current
    
    # getting rid values that have NA for both last contact and for survival
    if (is.null(nrow(X[[i]]))) {
      current_length = length(X[[i]])
    } else {
      current_length = nrow(X[[i]])
    }
    
    to_delete = c() # saving entries to delete because they are NAs
    j = 1 # index of entries to delete
    for (k in 1:current_length) {
      last_contact_missing = is.na(Last_Contact[[i]][k])
      NA_in_survival = is.na(Y[[i]][k])
      if (last_contact_missing & NA_in_survival) {
        to_delete[j] = k
        j = j + 1
      }
    }
    if (!is.null(to_delete)) {
      if (is.null(nrow(X[[i]]))) {
        Last_Contact[[i]] = Last_Contact[[i]][-to_delete]
        X[[i]] = X[[i]][-to_delete]
        Y[[i]] = Y[[i]][-to_delete]
      } else {
        Last_Contact[[i]] = Last_Contact[[i]][-to_delete]
        X[[i]] = X[[i]][-to_delete,]
        Y[[i]] = Y[[i]][-to_delete]
      }
    }
  }
  names(Last_Contact) = cancer_types
  
  # Getting rid of any observations who have a last contact value that is less than or equal to 0
  for (k in 1:n) {
    current_lc = Last_Contact[[k]]
    current_data = X[[k]]
    current_survival = Y[[k]]
    any_negative_lc = current_lc <= 0 # which values are negative, this will be mostly falses
    any_negative_lc[is.na(any_negative_lc)] = rep(FALSE, sum(is.na(any_negative_lc))) 
    
    if (is.null(nrow(current_data))) {
      current_data = current_data[!any_negative_lc]
      current_survival = current_survival[!any_negative_lc]
      current_lc = current_lc[!any_negative_lc]
    } else {
      current_data = current_data[!any_negative_lc, ]
      current_survival = current_survival[!any_negative_lc]
      current_lc = current_lc[!any_negative_lc]
    }
    
    X[[k]] = current_data # saving the updated variables back into the overall data
    Y[[k]] = current_survival
    Last_Contact[[k]] = current_lc
  }
  
  # Getting rid of observations whose survival time is 0
  for (k in 1:n) {
    current_lc = Last_Contact[[k]]
    current_data = X[[k]]
    current_survival = Y[[k]]
    
    any_negative_survival = current_survival <= 0 
    any_negative_survival[is.na(any_negative_survival)] = rep(FALSE, sum(is.na(any_negative_survival))) 
    
    current_data = current_data[!any_negative_survival, ]
    current_survival = current_survival[!any_negative_survival]
    current_lc = current_lc[!any_negative_survival]
    
    X[[k]] = current_data # saving the updated variables back into the overall data
    Y[[k]] = current_survival
    Last_Contact[[k]] = current_lc
  }
  
  return(list(Full = X, Survival = Y, LastContact = Last_Contact))
}

# Downloading the data
CancerData = GetTCGAData(clinical_data, cancer_types)
names(CancerData) = cancer_types 

# Removing any observations with NAs for age
clinical_data$age_at_initial_pathologic_diagnosis = as.numeric(as.character(clinical_data$age_at_initial_pathologic_diagnosis)) # will cause NAs
clinical_data = clinical_data[!is.na(clinical_data$age_at_initial_pathologic_diagnosis), ]

# matching the patient barcodes between the mutation data from TCGA2STAT and 
# TCGA clinical data 
CancerDataAndClinical = MatchBarCodes(CancerData, clinical_data, cancer_types) 

# The reordered data (data with bar codes matched)
CancerData_RO = CancerDataAndClinical$X_new
clinical_data_list = CancerDataAndClinical$Y_new 

# Now removing the cancer types not being considered (c("MESO", "PCPG", "PRAD", "TGCT", "THCA", "THYM"))
CancerData_RO_S = CancerData_RO[!(names(CancerData_RO) %in% c("MESO", "PCPG", "PRAD", "TGCT", "THCA", "THYM"))]
clinical_data_list_S = clinical_data_list[!(names(clinical_data_list) %in% c("MESO", "PCPG", "PRAD", "TGCT", "THCA", "THYM"))]
cancer_types_27 = cancer_types[!(cancer_types %in% c("MESO", "PCPG", "PRAD", "TGCT", "THCA", "THYM"))]

# Check that the barcodes are identical between the imported mutation data and the survival data
all.match = c()
for (i in 1:length(CancerData_RO_S)) {
  all.match[i] = identical(colnames(CancerData_RO_S[[i]]), as.character(clinical_data_list_S[[i]]$bcr_patient_barcode))
}
all.match

# Calculating the average mutation rate across all cancer types
AvgMutRate = MutationRateByCancer(CancerData_RO_S, cancer_types_27)
AvgMutRate = AvgMutRate[!(names(AvgMutRate) %in% c("Unknown", "."))] # Getting rid of the gene with the name "Unknown"
AvgMutRate = sort(AvgMutRate, decreasing = T)[1:50]
NewGenes = sort(names(AvgMutRate))

# Number of samples
sum(sapply(CancerData_RO_S, ncol))

# Selecting the genes of interest from the full dataset
CancerData_SelectGenes = SelectGenesFromOrigData(CancerData_RO_S, clinical_data_list_S, NewGenes, cancer_types_27)

# Check that each cancer type has the same number of columns/genes (should be length(NewGenes) = 50)
sapply(CancerData_SelectGenes, ncol)

# Check that there are no NAs in the age column
any.NAs = c()
for (i in 1:length(F27.50.2)) {
  any.NAs[i] = any(is.na(F27.50.2[[i]][,1]))
}
any.NAs

# Add the time of last contact to the mutation dataset
FS27.50 = AddAgeCensorToX(CancerData_SelectGenes, clinical_data_list_S, cancer_types_27, NewGenes) # will give warnings
F27.50.2 = FS27.50$X
S27.50.2 = FS27.50$Y

# Processing the data so that extraneous observations are removed
FS27.50.2 = PreProcessing(F27.50.2, S27.50.2, cancer_types_27)
F27.50.3 = FS27.50.2$Full
S27.50.3 = FS27.50.2$Survival
Last_Contact = FS27.50.2$LastContact

n.vec = sapply(F27.50.3, nrow)

# Checking for any NAs in the mutation data 
check = sapply(F27.50.3, function(i) apply(i, 2, function(k) any(is.na(k))))
any(check)

# Check for negative survival or last contact times
any.neg.surv = c()
any.neg.lc = c()
for (i in 1:length(F27.50.3)) {
  any.neg.surv[i] = any(na.omit(S27.50.3[[i]]) <= 0)
  any.neg.lc[i] = any(na.omit(Last_Contact[[i]]) <= 0)
}
any(any.neg.surv)
any(any.neg.lc)

# Check the data for NAs in both last contact and survival
same.length = c()
NA.lc.surv = c()
for (i in 1:length(S27.50.3)) {
  same.length[i] = length(S27.50.3[[i]]) == length(Last_Contact[[i]])
  for (j in 1:length(S27.50.3[[i]])) {
    NA.lc.surv[j] = is.na(S27.50.3[[i]][j]) & is.na(Last_Contact[[i]][j])
  }
}
same.length
any(NA.lc.surv)

# Check that the data used in the analysis matches the original clinical data
# even after processing. 

surv.time.matches = c()
cens.time.matches = c()
age.matches = c()
for (i in 1:length(cancer_types_27)) {
  current.cans = F27.50.3[[i]]
  surv.i = c()
  cens.i = c()
  ages.i = c()
  for (j in 1:nrow(F27.50.3[[i]])) {
    current.id = rownames(current.cans)[j]
    current.clin = clinical_data[clinical_data$bcr_patient_barcode == current.id,] 
    current.surv = current.clin$death_days_to
    current.cens = current.clin$last_contact_days_to
    current.age = current.clin$age_at_initial_pathologic_diagnosis
    surv.saved = S27.50.3[[i]][j]
    cens.saved = Last_Contact[[i]][j]
    age.saved = F27.50.3[[i]][,1][j]
    surv.i[j] = if (!is.na(surv.saved)) surv.saved == current.surv else (is.na(surv.saved) & current.surv == "#N/A")
    cens.i[j] = if (!is.na(cens.saved)) cens.saved == current.cens else (is.na(cens.saved) & current.cens == "#N/A")
    ages.i[j] = age.saved == current.age
  }
  surv.time.matches[i] = all(surv.i)
  cens.time.matches[i] = all(cens.i)
  age.matches[i] = all(ages.i)
}
all(surv.time.matches)
all(cens.time.matches)
all(age.matches)

# Save the mutation data, survival data, names of the 27 cancer types, and the genes selected in an RDA file called "FSCG.rda"
save(F27.50.3, S27.50.3, Last_Contact, cancer_types_27, NewGenes, file = "FSCG.rda", version = 2)


