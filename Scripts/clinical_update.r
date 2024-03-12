
###
#load required packages
library(dplyr)
library(tidyr)

# Load the data
df_hiv <- read.delim("./Data/edtcp_clinical.txt")

#Load the plink formated file with participant IDs
pheno <- read.table("./Data/Genotype/edtcp.fam")

#rename the columns of the phenotype data file 
names(pheno)<-c('FID', 'IID', 'PID', 'MID', 'SEX', 'P')

#select required columns
edtcp_clinical <- select(df_hiv, FID,IID, sex,phenotype)

#Subset baseline data for patients with matching IID in genotype data
edtcp <- subset(edtcp_clinical, IID %in% pheno$IID)

#Rearrange the data in the order provided in the genotype file
edtcp_pheno <- edtcp[order(match(edtcp$IID, pheno$IID)), ]

#Check for duplicated rows based on IID
duplicated_ID <- edtcp_pheno[duplicated(edtcp_pheno$IID) 
                    | duplicated(edtcp_pheno$IID, fromLast = TRUE), ]

#Select FID,IID,sex and phenotype columns
new_phenotype <- subset(edtcp_pheno, select = c(IID,sex,phenotype))

#now drop duplicated rows from the data using distinct
new_pheno <- distinct(new_phenotype)

#combine the ordered data frame and the GWAS.fam data frame into a 
#single data frame by merging the data frames based on IID
pheno_data <-pheno %>%
  left_join(new_pheno, by = "IID") %>%
  select(FID, IID, PID, MID,SEX= sex, P=phenotype)

# Replace "NAN" values with 0 in column 5 and -9 in column 6
pheno_data$SEX[is.na(pheno_data$SEX)] <- 0
pheno_data$P[is.na(pheno_data$P)] <- -9
#save the output 
write.table(pheno_data, file = "./Data/Genotype/clinical.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")




