library(scales)

load('FSCG.rda')
load('~/Dropbox/PanCancerOutcome/Data/FinalModelFromForwardSelectionFAT4_TP53.rda')

SubsetNGenes = function(F33.55.2, Genes) {
  # Subsets entire dataset to just a selected set of genes
  n = length(F33.55.2)
  F33.N.2 = F33.55.2 # note the difference here - N instead of 55
  for (i in 1:n) {
    current_type = F33.55.2[[i]]
    current_genes = colnames(current_type)
    current_available = current_genes %in% Genes
    current_subset = current_type[, current_available]
    current_subset = cbind(current_type[,1:2], current_subset)
    colnames(current_subset) = c("Age", "LastContact", current_genes[current_genes %in% Genes])
    F33.N.2[[i]] = current_subset
  }
  return(F33.N.2)
}

F27.3.3 = SubsetNGenes(F27.50.3, c("FAT4","TP53"))
S27.3.3 = S27.50.3

AgeMeans <- c() ##mean age for each cancer type (to adjust for centering)
for(i in 1:length(F27.3.3)){
  AgeMeans[i] <- mean(F27.3.3[[i]][,1])
}
age.vec <- 25+5*c(1:12)

Time <- seq(1,3650/2, length.out=100) #5-year survival

genes_mutated = matrix(c(NA, "FAT4", NA, "TP53"), ncol = 2, byrow = T)
Mat <- matrix(unlist(FinalModelPosteriors$betas[[i]]), ncol=5,byrow=TRUE)[10001:20000,]
#j indexes age, k MUC5B status, l PIK3CA, and m TP53
age.vec.30.50.80 = c(30, 50, 80)

i = which(cancer_types_27 == "ACC")
Mat <- matrix(unlist(FinalModelPosteriors$betas[[i]]), ncol=4,byrow=TRUE)[10001:20000,]

png('ACC_figs_v2.png',width=900,height=700)
par(mfrow = c(2,2))
colors = c("black", "orange", "deepskyblue")
 # for each age potential age (of 30, 50, 80)
for(Fig in 1:4){
    if(Fig==1){ k=1; l=1} 
    if(Fig==2){ k=2; l=1}
    if(Fig==3){ k=1; l=2}
    if(Fig==4){ k=2; l=2}
    for(j in 1:length(age.vec.30.50.80)){
      age = age.vec.30.50.80[j]
      X <- as.matrix(c(1,age-AgeMeans[i],k-1,l-1))
      mus <- Mat%*%X
      sds <- sqrt(unlist(FinalModelPosteriors$sigma2))[10001:20000]
      SurvMat <- matrix(ncol=length(Time),nrow=10000)
      for(t in 1:length(Time)){
        SurvMat[,t] <- plnorm(Time[t], meanlog=mus, sdlog = sds,lower = FALSE)
      } #save figure
      # png(paste0('surv_',i,'_',age.vec[j],'_',k-1,'_',l-1,'_',m-1,'.png'),width=400,height=400)
      current_status = c(k,l)
      mutations = sapply(1:nrow(genes_mutated), function(i) genes_mutated[i, ][current_status[i]])
      mutations[is.na(mutations)] = rep(' ', sum(is.na(mutations)))
      
      if (k ==1 & l == 1) {
        title = paste("No Mutations")
      } else {
        title = paste("Mutation at: ", paste(mutations, collapse = ' '), sep = '', collapse = '')
      }
      if (j == 1) {
        plot(Time/365,colMeans(SurvMat), type = 'l', ylim=c(0,1), ylab = "Probability of survival", xlab = 'Time (years)', 
             main = title,lwd=3, col = colors[j], cex.main=1.25, cex.lab=1.5, cex.axis=1.0)
      } else {
        lines(Time/365, colMeans(SurvMat), ylim=c(0,1), lwd = 3, col = colors[j])
      }
      points(Time/365,apply(SurvMat,2,quantile,probs=c(0.975)), type = 'l', lty=2, col = alpha(colors[j], 0.3), lwd = 2)
      points(Time/365,apply(SurvMat,2,quantile,probs=c(0.025)), type = 'l', lty=2, col = alpha(colors[j], 0.3), lwd = 2)
      # dev.off()
  }
}
dev.off()

