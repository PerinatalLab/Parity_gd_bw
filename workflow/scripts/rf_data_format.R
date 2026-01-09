########### r-packages needed ###########
library(data.table)
library(dplyr)
library(ranger)
library(ggplot2)


########### Loading data ###########
dosage = fread(snakemake@input[[1]])  # Matrix or data frame of genetic variables such as SNPs coded as 0-1-2
outcome = fread(snakemake@input[[2]])     # Numeric vector of the outcome/phenotype
covariates = fread(snakemake@input[[3]])            # Optional data frame containing potentially confounding variables to be adjusted for
filename = snakemake@input[[1]]


########### Prep data fromat ###########

#### SNPS as columns and observations as rows in dosage ####
colnames(dosage)[1]="chr"
dosage$chr = paste("CHROM",dosage$chr, sep="")
dosage$chom_pos = paste(dosage$chr, dosage$POS, sep="POS")
dosage = dosage %>% filter(ifelse(POS == 114296770 & ALT =="C",FALSE,TRUE))
names = dosage[,dosage$chom_pos]
dosage = dosage %>% select(-chr,-POS,-REF,-ALT,-chom_pos)
dosage = round(t(dosage))
colnames(dosage) = names
roworder = rownames(dosage)
dosage = tibble::rownames_to_column(as.data.frame(dosage), var = "IID")


#### Only one ID colum ###
outcome = outcome %>% slice(match(roworder, IID)) %>% select(-FID)
covariates = as.data.frame(covariates %>% slice(match(roworder, IID)) %>% select(-FID))


#### Combinde dosage, outcome and covariate ###
dat1 = inner_join(outcome, dosage, by="IID")
dat = inner_join(dat1, covariates, by ="IID")
rm(dat1)


#### Parameter values numeric or factor, exclude parity ###
dat = dat %>% mutate(across(-c(genotyping_chip,IID), as.numeric)) 
if (any(colnames(dat) == 'parity')) {
	dat = dat %>% select(-parity)}
dat$genotyping_chip = as.factor(dat$genotyping_chip)


if (grepl("dosage/bw",filename)) {
	### Residualized birth weight zscores ###
        lmm = lm(bw_zscore ~ pregnancy_duration + sex  + genotyping_chip + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG, dat = dat)
        dat$bw_res = mean(dat$bw_zscore)+lmm$residuals

        dat2 = dat %>%
                select(-bw_zscore, -genotyping_chip, -sex, -pregnancy_duration,-PC1_AVG,-PC2_AVG,-PC3_AVG,-PC4_AVG,-PC5_AVG,-PC6_AVG,-PC7_AVG,-PC8_AVG,-PC9_AVG,-PC10_AVG) %>%
                rename("bw_zscore"="bw_res") # genetic features only
        dat = dat %>% select(-bw_res) # genetic, gestational duration, sex, chip and the ten first genotype principal components

} else {
	### Residualized gestational durations ###
	lmm = lm(gd ~ KJONN + genotyping_chip + poly(maternalage,2) + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG, dat = dat)
	dat$gd_res = mean(dat$gd)+lmm$residuals

	dat2 = dat %>% 
		select(-gd, -KJONN, -genotyping_chip, -maternalage,-PC1_AVG,-PC2_AVG,-PC3_AVG,-PC4_AVG,-PC5_AVG,-PC6_AVG,-PC7_AVG,-PC8_AVG,-PC9_AVG,-PC10_AVG) %>% 
		rename("gd"="gd_res") # genetic features only
	dat = dat %>% select(-gd_res) # genetic, maternal age, sex, chip and the ten first genotype principal components
}

fwrite(dat, snakemake@output[[1]], sep=",") 
fwrite(dat2, snakemake@output[[2]], sep=",") 

