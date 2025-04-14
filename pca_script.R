library(vcfR)
library(vegan)
library(ggplot2)
library(ggpubr)

#reading & transposing data

snps_filter <- vcfR::read.vcfR("raw.maf5.minq40.miss60.mdp3.recode.vcf", convertNA  = TRUE)

snps_num_f <- vcfR::extract.gt(snps_filter, 
                              element = "GT",
                              IDtoRowNames  = F,
                              as.numeric = T,
                              convertNA = T,
                              return.alleles = F)


snps_num_t_f <- t(snps_num_f) 

# COunting and removing NAs 

find_NAs <- function(x){
  NAs_TF <- is.na(x)
  i_NA <- which(NAs_TF == TRUE)
  N_NA <- length(i_NA)
  
  cat("Results:",N_NA, "NAs present\n.")
  return(i_NA)
}

N_rows <- nrow(snps_num_t_f)
N_NA   <- rep(x = 0, times = N_rows)
N_SNPs <- ncol(snps_num_t_f)

for(i in 1:N_rows){
  i_NA <- find_NAs(snps_num_t_f[i,]) 
  N_NA_i <- length(i_NA)
  N_NA[i] <- N_NA_i
}  

percent_NA <- N_NA/N_SNPs*100
print(percent_NA)

# happened to not use as the % is less than 50 everywhere
# i_NA_50percent <- which(percent_NA > 50) 
# snps_num_t_f1 <- snps_num_t_f[-i_NA_50percent, ]

invar_omit <- function(x){
  cat("Dataframe of dim",dim(x), "processed...\n")
  sds <- apply(x, 2, sd, na.rm = TRUE)
  i_var0 <- which(sds == 0)
  cat(length(i_var0),"columns removed\n")
  if(length(i_var0) > 0){
    x <- x[, -i_var0]
  }
  return(x)                      
}

# Removing invariant SNPs

snps_f_no_invar <- invar_omit(snps_num_t_f) 

N_col <- ncol(snps_f_no_invar)

for(i in 1:N_col){
  # get the current column
  column_i <- snps_f_noNAs[, i]
  # get the mean of the current column
  mean_i <- mean(column_i, na.rm = TRUE)
  # get the NAs in the current column
  NAs_i <- which(is.na(column_i))
  # record the number of NAs
  N_NAs <- length(NAs_i)
  # replace the NAs in the current column
  column_i[NAs_i] <- mean_i
  # replace the original column with the
  ## updated columns
  snps_f_noNAs[, i] <- column_i
}

# check
all(is.finite(snps_f_noNAs))

# scaling

SNPs_scaled_f <- scale(snps_f_noNAs)

# PCA

# if doing all over
# saveRDS(SNPs_scaled_f, file = "matrix_filtered.rds")
# SNPs_scaled_f <- readRDS("matrix_filtered.rds")

pca_scaled_f <- prcomp(SNPs_scaled_f)

# Graph for PCs' relative importance

screeplot(pca_scaled_f, 
          ylab  = "Relative importance",
          main = "SNPs Data Analysis")

summary_out_scaled_f <- summary(pca_scaled_f)

PCA_variation <- function(pca_summary, PCs = 2){
  var_explained <- pca_summary$importance[2,1:PCs]*100
  var_explained <- round(var_explained,1)
  return(var_explained)
}

var_out <- PCA_variation(summary_out_scaled_f,PCs = 10)

N_columns <- ncol(SNPs_scaled_f)
barplot(var_out,
        main = "Percent variation Scree plot",
        ylab = "Percent variation explained")
abline(h = 1/N_columns*100, col = 2, lwd = 2)

biplot(pca_scaled_f)

# Relative importance in numbers

eigenvectors <- pca_scaled_f$rotation
eigenvectors

eigenvalues <- pca_scaled_f$sdev^2
eigenvalues
total_variance <- sum(eigenvalues)
variance_explained <- eigenvalues / total_variance * 100
variance_explained

# Making populations vector

popullations <- read.table("maf5q40miss60dp3.list", header = FALSE)
popul_vector <- popullations$V2

# PCA graph, colored by populations 

pca_scores <- vegan::scores(pca_scaled_f)

pca_scores2 <- data.frame(popul_vector,
                          pca_scores)

ggpubr::ggscatter(data = pca_scores2,
                  y = "PC2",
                  x = "PC1",
                  color = "popul_vector",
                  shape = "popul_vector",
                  xlab = "PC1",
                  ylab = "PC2",
                  main = "PCA scores")


ggpubr::ggscatter(data = pca_scores2,                   
                 y = "PC2",                   
                 x = "PC1",                   
                 color = "popul_vector",                   
                 shape = "popul_vector",                   
                 xlab = "PC1",                   
                 ylab = "PC2",                   
                 main = "PCA scores") +
        geom_text(aes(label = rownames(pca_scores2), color = popul_vector), 
                  vjust = -0.5,   
                  hjust = 0.5,    
                  size = 3,  
                  check_overlap = TRUE)  

