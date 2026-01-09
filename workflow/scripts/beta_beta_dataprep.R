########### r-packages needed ###########
library(dplyr)
library(data.table)




########### Loading data  ###########
gwas1 = fread(snakemake@input[[1]], select = c("ID","BETA","SE","PVALUE"))
gwas2 = fread(snakemake@input[[2]], select = c("ID","BETA","SE","PVALUE"))
topsnp1 = fread(snakemake@input[[3]])[grep("5e-08", cut_off, ignore.case = TRUE)]
topsnp2 = fread(snakemake@input[[4]])[grep("5e-08", cut_off, ignore.case = TRUE)]




########### Format and merge topsnps from the two parity groups ###########
gwas1$CImin = gwas1$BETA - 1.96*gwas1$SE; gwas2$CImin = gwas2$BETA - 1.96*gwas2$SE
gwas1$CImax = gwas1$BETA + 1.96*gwas1$SE; gwas2$CImax = gwas2$BETA + 1.96*gwas2$SE

topsnp_merged = full_join(topsnp1, topsnp2, by="ID")[,c("ID")]
dat = left_join(topsnp_merged, gwas1, by="ID")
colnames(dat) <- c("ID", "BETA0","SE0","PVALUE0", "CImin0", "CImax0")
dat = left_join(dat, gwas2, by="ID")




########### Saving data  ###########
fwrite(dat, snakemake@output[[1]], sep = "\t")
