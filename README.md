# Pan-Cancer-Survival-Modeling

This is the repository with code used in the Pan-Cancer Survival Modeling project, titled "A Pan-Cancer and Polygenic Bayesian Hierarchical Model for the Effect of Somatic Mutations on Survival." In this repository, one will find code for (1) data importation and gene set selection, (2) the four models we tried and the calculation of the cross-validated out-of-sample log-posterior likelihoods to compare them, (3) the forward selection procedure we employed, (4) the survival plots we created, and (5) the simulation we used to validate our in-house Gibbs sampler.

# Description of each file 
(1) ImportData.R: In this file, one can find code on importing data through the TCGA2STAT pipeline and how we selected genes for our study. One will need to save the data from this code in order to run the rest of the analysis. At the bottom of this script there is a save() statement and at the top of each additional analysis script there is a load() statement to facilitate this. 

(2) The four models we tried were the log-normal, normal, exponential, and Weibull. The corresponding scripts for these models are as follows: PL5FCVLogNorm50Genes.R, PL5FCVNorm50Genes.R, JagsModelExponential.R, JagsModelWeibull.R. Based on the results in these, we used the log-normal model to continue our analysis. 

(3) To do forward selection, first we ran the forward selection procedure without any genes. This served as our baseline. The code for this can be found in ForwardSelectionNoGenes.R. Once this has been run, then one can move onto comparing single-gene models, two-gene models, and so on. 

In ForwardSelection_5FoldCV.R, one can find the code to run the forward selection procedure. Running parallel loops in R with the foreach package requires the user to import everything in the outside environment into the parallel environment which is a bit tricky. To avoid this, we chose not write the entire forward selection procedure in a function. This also encourages the user to come back and check every once in awhile where the procedure is at and see how the results are progressing. For this reason, in order to run the one-gene model and the two-gene model and so on, copy-and-paste the chunk of code between the ### comments. Every line with a ### between these two comments must be changed at each step of the forward selection procedure. 

(5) Survival_Curves_ACC.R: In this script, one can find the code we used to create the survival curves for a patient with ACC depending on various ages and combinations of mutation statuses. The script imports the posteriors from fitting the final model, the .rda file of which can be found in this repository. 

(6) SimulateCoverage.R: Here, one can find the simulation we used to assess the coverage rates of our in-house Gibbs sampler. This script generates 1000 datasets based on fixed true values of the parameters, calculates the posteriors for the parameters using our model, calculates the credible intervals for each, and then checks if the true value of each parameter used to generate the data is contained in the credible interval. 



