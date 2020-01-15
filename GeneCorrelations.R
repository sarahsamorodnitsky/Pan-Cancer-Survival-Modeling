# Calculating the correlations between genes across all cancer types, performing a 
# correlation test to assess the significance of each correlation, and creating a 
# plot to show which correlations are significant and their magnitude of correlation. 
# Author: Sarah Samorodnitsky and Eric Lock

load("FSCG.rda")

library(corrplot)

# Removing age so that it is not included in correlation calculations
F27.50.3.noage = lapply(1:length(F27.50.3), function(i) F27.50.3[[i]][,-1])

# Combining the mutation data from each of the cancers together into one large matrix
# MutationListMatrix = matrix(unlist(F27.50.3.noage), ncol = ncol(F27.50.3.noage[[1]]), byrow = T)
MutationListMatrix = do.call(rbind, F27.50.3.noage)

# Check that the mutation data was combined correctly
n.vec = sapply(F27.50.3.noage, nrow)
sum(n.vec) == nrow(MutationListMatrix) # total numbers of patients align

k = 1 # index for MutationListMatrix (1, ..., 5698)
combined = c()
for (i in 1:length(n.vec)) { # index for cancer type
  same.i = c()
  for (j in 1:n.vec[i]) { # index for F27.50.3.noage
    same.i[j] = all(MutationListMatrix[k, ] == F27.50.3.noage[[i]][j, ])
    k = k + 1
  }
  combined[i] = all(same.i)
}
all(combined) 

# Naming the MutationListMatrix 
colnames(MutationListMatrix) = colnames(F27.50.3.noage[[1]])

# Creating the correlation matrix
GeneCorrelations <- cor(MutationListMatrix)
rownames(GeneCorrelations) = colnames(F27.50.3.noage[[1]])
colnames(GeneCorrelations) = colnames(F27.50.3.noage[[1]])

# Calculating significance

cor.mtest <- function(mat) {
  # Performs Fisher's exact test on every gene mutation vector
  # Returns a matrix of p-values (only the upper half of the matrix is filled in)
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      geneA = mat[, i]
      geneB = mat[, j]
      tmp = fisher.test(x=as.factor(geneA),y=as.factor(geneB))
      p.mat[i, j] = tmp$p.value
    }
  }
  colnames(p.mat) <- colnames(mat)
  rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

# Matrix of p-values
p.mat <- cor.mtest(MutationListMatrix)

# Number of significant correlations
n.compare = sum(!is.na(p.mat) & p.mat != 0) # number of unique comparisons made
n.sig = sum(p.mat[!is.na(p.mat) & p.mat != 0] <= 0.05)
n.sig/n.compare # percentage of correlation sthat were significant


# Plotting
col1 = colorRampPalette(c("medium blue", "white", "red")) # setting the colors
corrplot(GeneCorrelations, type="upper", order="original", 
         p.mat = p.mat, sig.level = 0.05,
         col = col1(20),
         tl.col="black",
         insig = "blank")
