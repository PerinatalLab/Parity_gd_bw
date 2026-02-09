library(tableone)
library(data.table)
library(dplyr)

########### Loading data needed ###########
pheno = fread(snakemake@input[[2]],drop=c("FID")) # pheno file
covar = fread(snakemake@input[[1]], drop= c("FID","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG")) # covar file
outcome = colnames(pheno)[2]
colnames(pheno) = c("IID","outcome")


########### Creating table ###########

dat = inner_join(pheno, covar, by = "IID") # merge 

if (outcome == "gd") {
	myVars = c("outcome", "maternalage","KJONN","genotyping_chip") # Vector of variables to summarize
	catVars = c("KJONN","genotyping_chip") # Vector of categorical variables 
	nonnormalVar = c("outcome", "maternalage") # Summarizing nonnormal variables (median)
	
	tab3 = CreateTableOne(vars = myVars, strata = "parity" , data = dat, factorVars = catVars) # Create a TableOne object, parity summary
	tab3Mat = print(tab3, nonnormal = nonnormalVar, formatOptions = list(big.mark = ","), howAllLevels = TRUE,quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

} else if (outcome == "parity") {
	myVars = c("outcome", "maternalage","genotyping_chip") # Vector of variables to summarize
        catVars = c("genotyping_chip") # Vector of categorical variables 
        nonnormalVar = c("outcome", "maternalage") # Summarizing nonnormal variables (median)

        tab3 = CreateTableOne(vars = myVars, strata = "outcome" , data = dat, factorVars = catVars) # Create a TableOne object, parity summary
        tab3Mat = print(tab3, nonnormal = nonnormalVar, formatOptions = list(big.mark = ","), howAllLevels = TRUE,quote = TRUE, noSpaces = TRUE, printToggle = FALSE)

} else {
	myVars = c("outcome", "maternalage","SVLEN_DG","KJONN","genotyping_chip") # Vector of variables to summarize
	catVars = c("KJONN","genotyping_chip") # Vector of categorical variables 
	nonnormalVar = c("outcome", "SVLEN_DG","maternalage") # Summarizing nonnormal variables (median)
	
	tab3 = CreateTableOne(vars = myVars, strata = "parity" , data = dat, factorVars = catVars) # Create a TableOne object, parity summary
	tab3Mat = print(tab3, nonnormal = nonnormalVar, formatOptions = list(big.mark = ","), howAllLevels = TRUE,quote = TRUE, noSpaces = TRUE, printToggle = FALSE)
}




########### Saving table ###########
fwrite(tab3Mat, snakemake@output[[1]])
