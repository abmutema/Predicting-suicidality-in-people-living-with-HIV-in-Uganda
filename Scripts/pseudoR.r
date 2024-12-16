
# Load required libraries

library(stats)

# Load the PRS values
prs_data <- read.table("UG.best", header = TRUE)

# Load the covariates
covariate_data <- read.table("mhs.cov", header = TRUE)

# Load the principal components
pc_data <- read.table("mhs.eigenvec", header = FALSE)

# Assign column names: FID, IID, and PC1, PC2, ..., PCn
num_pcs <- ncol(pc_data) - 2  # Number of PCs
colnames(pc_data) <- c("FID", "IID", paste0("PC", 1:num_pcs))

# Load the phenotype data from mhs.fam
fam_data <- read.table("mhs.fam", header = FALSE)

# Assign column names to mhs.fam
colnames(fam_data) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype")

# Ensure Phenotype is binary 
fam_data$Phenotype <- ifelse(fam_data$Phenotype == 2, 1, 0)

# Merge all data sets based on FID and IID
merged_data <- merge(prs_data, covariate_data, by = "FID")

merged_data <- merge(merged_data, pc_data, by = "FID")

merged_data <- merge(merged_data, fam_data[, c("FID", "IID", "Phenotype")], by = c("FID", "IID"))

# Extract variables for modeling
phenotype <- merged_data$Phenotype  
prs <- merged_data$PRS            
covariates <- merged_data[, c("Age", "MDD","studsite")]  
pcs <- merged_data[, grep("PC", names(merged_data))] 

# Combine covariates and PCs
covariate_matrix <- cbind(covariates, pcs)

# Fit logistic regression models
# Null model (without PRS)
null_model <- glm(phenotype ~ ., data = covariate_matrix, family = binomial)

# Full model (with PRS added)
full_model <- glm(phenotype ~ . + prs, data = cbind(covariate_matrix, prs = prs), family = binomial)

# Compute Nagelkerke's pseudo-R2 
ll_null <- logLik(null_model)

ll_full <- logLik(full_model)

n <- length(phenotype)

nagelkerke_r2 <- (1 - exp((2 / n) * (ll_null - ll_full))) / (1 - exp((2 / n) * ll_null))

cat("Nagelkerke pseudo-R2:", nagelkerke_r2, "\n")


