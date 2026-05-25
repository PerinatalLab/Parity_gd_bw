########### r-packages needed ###########
library(data.table)
library(dplyr)
library(ranger)
library(ggplot2)


########### Loading data ###########
dosage = fread(snakemake@input[[1]])  # Matrix or data frame of genetic variables such as SNPs coded as 0-1 (mt, mnt, pt)
outcome = fread(snakemake@input[[2]])     # Numeric vector of the outcome/phenotype
covariates = fread(snakemake@input[[3]])            # Optional data frame containing potentially confounding variables to be adjusted for
filename = snakemake@input[[1]]


########### Prep data fromat ###########

#### SNPS as columns and observations as rows in dosage ####
format_haps= function(hap){
  variants = t(hap)[1,]  
  ids= names(hap)[2:ncol(hap)]
  hap= as.data.frame(t(hap[, 2:ncol(hap)]))
  names(hap)= variants
  hap$IID <- sub("^.*:", "", ids)
  return(hap)
}

dosage = format_haps(dosage)
print(nrow(dosage))
roworder = rownames(dosage)

#### Only one ID colum ###
outcome = outcome  %>% select(-FID)
covariates = covariates %>% select(-FID)


#### Combinde dosage, outcome and covariate ###
dat1 = inner_join(outcome, dosage, by="IID")
dat = inner_join(dat1, covariates, by ="IID")
rm(dat1)
print(nrow(dat))

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

print(nrow(dat))

fwrite(dat, snakemake@output[[1]], sep=",") 
fwrite(dat2, snakemake@output[[2]], sep=",") 
