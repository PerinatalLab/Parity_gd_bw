########### r-packages needed ###########
library(dplyr)
library(data.table)




########### Loading data  ###########
PGS = fread(snakemake@input[[1]], drop = c("FID"))
covar = fread(snakemake@input[[2]], drop = c("FID"))
pheno = fread(snakemake@input[[3]], drop = c("FID"))




########### Merge data and only keep individuals with genotyping_chip "Illumina_GlobalScreeningArray_MD_v.3.0" ###########
# to only keep only the samples newly genotyped from moba

colnames(pheno) = c("IID","outcome")

dat = inner_join(PGS, covar, by="IID")
dat = inner_join(dat, pheno, by="IID")

dat = dat %>% filter(genotyping_chip == "Illumina_GlobalScreeningArray_MD_v.3.0")

cat("Number of rows in data frame:", nrow(dat), "\n")
print(head(dat))




########### Linear models ###########
if (grepl("clean_phenotypes/gd", snakemake@input[[2]])) {
	mod = lm(as.numeric(outcome) ~ scale(as.numeric(SCORESUM),scale = F),data = dat )# + KJONN + poly(maternalage,2) + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG, data = dat)
	cat("Whole population","\n");print(summary(mod))
	coeff = as.data.frame(summary(mod)$coefficients); coeff$lm = "lm"; coeff$r2 = summary(mod)$adj.r.squared; coeff$n = nrow(dat)


	mod_p0 = lm(as.numeric(outcome) ~ scale(as.numeric(SCORESUM),scale = F), data = filter(dat, parity ==0)) #+ KJONN + poly(maternalage,2) + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG, data = filter(dat, parity ==0))
	cat("Parity 0","\n"); print(summary(mod_p0))
	coeff_p0 = as.data.frame(summary(mod_p0)$coefficients); coeff_p0$lm = "p0"; coeff_p0$r2 = summary(mod_p0)$adj.r.squared; coeff_p0$n = nrow(filter(dat, parity ==0))


	mod_p1 = lm(as.numeric(outcome) ~ scale(as.numeric(SCORESUM),scale = F), data = filter(dat, parity ==1)) #+ KJONN + poly(maternalage,2) + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG, data = filter(dat, parity ==1))
	cat("Parity > 0","\n"); print(summary(mod_p1))
	coeff_p1 = as.data.frame(summary(mod_p1)$coefficients); coeff_p1$lm = "p1"; coeff_p1$r2 = summary(mod_p1)$adj.r.squared; coeff_p1$n = nrow(filter(dat, parity ==1))


	mod_gxe = lm(as.numeric(outcome) ~ scale(as.numeric(SCORESUM),scale = F) * as.factor(parity),data = dat) #+ KJONN + poly(maternalage,2) + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG, data = dat)
	cat("GxE","\n");print(summary(mod_gxe))
	coeff_gxe = as.data.frame(summary(mod_gxe)$coefficients); coeff_gxe$lm = "gxe";  coeff_gxe$r2 = summary(mod_gxe)$adj.r.squared; coeff_gxe$n = nrow(dat)

} else if (grepl("zscore", snakemake@input[[2]])) {
        mod = lm(as.numeric(outcome) ~ scale(as.numeric(SCORESUM),scale = F), data=dat) #+ sex + pregnancy_duration + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG, data = dat)
        cat("Whole population","\n");print(summary(mod))
	coeff = as.data.frame(summary(mod)$coefficients); coeff$lm = "lm"; coeff$r2 = summary(mod)$adj.r.squared; coeff$n = nrow(dat)


        mod_p0 = lm(as.numeric(outcome) ~ scale(as.numeric(SCORESUM),scale = F), data = filter(dat, parity ==0))# + sex + pregnancy_duration + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG, data = filter(dat, parity ==0))
        cat("Parity 0","\n");print(summary(mod_p0))
        coeff_p0 = as.data.frame(summary(mod_p0)$coefficients); coeff_p0$lm = "p0"; coeff_p0$r2 = summary(mod_p0)$adj.r.squared; coeff_p0$n = nrow(filter(dat, parity ==0))


        mod_p1 = lm(as.numeric(outcome) ~ scale(as.numeric(SCORESUM),scale = F), data = filter(dat, parity ==1))# + sex + pregnancy_duration + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG, data = filter(dat, parity ==1))
        cat("Parity > 0","\n");print(summary(mod_p1))
        coeff_p1 = as.data.frame(summary(mod_p1)$coefficients); coeff_p1$lm = "p1"; coeff_p1$r2 = summary(mod_p1)$adj.r.squared; coeff_p1$n = nrow(filter(dat, parity ==1))


        mod_gxe = lm(as.numeric(outcome) ~ scale(as.numeric(SCORESUM),scale = F) * as.factor(parity), data = dat) #+ sex + pregnancy_duration + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG, data = dat)
        cat("GxE","\n");print(summary(mod_gxe))
        coeff_gxe = as.data.frame(summary(mod_gxe)$coefficients); coeff_gxe$lm = "gxe"; coeff_gxe$r2 = summary(mod_gxe)$adj.r.squared; coeff_gxe$n = nrow(dat)
}




########### Save data ###########
dat = rbind(coeff,coeff_p0,coeff_p1,coeff_gxe)
dat$names = rownames(dat)
fwrite(dat, snakemake@output[[1]], sep = "\t")
