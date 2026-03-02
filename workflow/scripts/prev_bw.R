library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(grid)



mfr = fread(snakemake@input[[1]])
parental_ids = fread(snakemake@input[[2]])
mfr_time_of_death = fread(snakemake@input[[3]])
sentrix_m = fread(snakemake@input[[4]])
sentrix_c = fread(snakemake@input[[5]])
qc_geno = fread(snakemake@input[[6]])


########### Child or mother ########### 
if (grepl("mother", snakemake@input[[7]])) {
  motherorchild = "mother"
  cat("Preparing the data for maternal genome and the mothers previous bw")
} else if (grepl("child", snakemake@input[[7]])) {
  motherorchild = "child"
  cat("Preparing the data for fetal genome and the mothers previous bw")
} 





########### CLEANING MFR ###########
cat("Initial number of pregnanices in medical birth registry (mfr): ", nrow(mfr),"\n")

#### No birth resulting in twins, triplet etc. ####
# FLERFODSEL =  The birth is a plural birth (twins, triplet etc.), if FLERFODSEL = 1 the delivery was a multiple birth
mfr = mfr %>% filter(is.na(FLERFODSEL))
cat("Number of pregnanices left after removing multiple deliveries (twins, etc): ", nrow(mfr),"\n")


#### Remove pregnancies treated with Asisted reproductive technology (ART) ####
# ART = Assisted Reproductive Technology treatment (IVF, ICSI or other ART treatments), if not NA the woman was treated with ART
mfr = mfr %>% filter(is.na(ART))
cat("Number of pregnanices left after removing woman treated with ART: ", nrow(mfr),"\n")


#### Unreasonable gestational age ####
# SVLEN_DG = Length of gestation in days (best estimate)
mfr = mfr %>% filter(SVLEN_DG >= 259 & SVLEN_DG <300 & !is.na(SVLEN_DG))
cat("Number of pregnanices left after excluding pregnancies with missing gestational length, is shorter than 154 days or is longer than 309 days: ", nrow(mfr),"\n")


#### Remove children born with malformations ####
# MISD = Congenital malformations, if MISD = 1 the child was born with malformations
mfr = mfr %>% filter(is.na(MISD))
cat("Number of pregnanices left after removing children born with malformations: ", nrow(mfr),"\n")


#### Remove children dying the same year as they was born ####
mfr = left_join(mfr,mfr_time_of_death, by="PREG_ID_1724")
mfr = mfr %>% filter(DODKAT != "Non-live birth")
cat("Number of pregnanices left after removing stillbirths: ", nrow(mfr),"\n")




########### Add maternal and paternal IDs to mfr via parental_ids ###########
# Matching by PREG_ID_1724
mfr = inner_join(mfr, parental_ids, by="PREG_ID_1724")
cat("Adding a maternal and paternal ID to each pregnancy from ''parental_ID_to_PREG_ID.csv''.", "\n")


#### Remove mothers with uncertain parity  ####
# Here we remove mothers which  has several children as the same parity, and where the parity do not match the year of delivery (e.g. parity 0 was delivered later than a parity 1)
mfr = mfr %>% group_by(M_ID_1724,PARITET_5) %>%
  filter(n() == 1 | PARITET_5 == 4) %>%
  ungroup() %>% group_by(M_ID_1724) %>% arrange(PARITET_5) %>%
  mutate(diff_parity = FAAR - dplyr::lag(FAAR)) %>%
  filter(!any(diff_parity<0,na.rm = T))
cat("Pregnencies left after removing mothers which  has several children as the same parity, and where the parity do not match the year of delivery: ", nrow(mfr), "\n")


#### Calculating the maternal age when delivering child ####
mfr = mfr %>% mutate(maternalage = FAAR - MOR_FAAR)


#### Birth weigth z-score #####
mfr = mfr %>% filter(!is.na(VEKT),VEKT < mean(mfr$VEKT,na.rm=TRUE)+5*sd(mfr$VEKT,na.rm=TRUE), VEKT > mean(mfr$VEKT,na.rm=TRUE)-5*sd(mfr$VEKT,na.rm=TRUE))

mfr$bw_zscore = (mfr$VEKT - mean(mfr$VEKT))/sd(mfr$VEKT)
cat("Number of pregnanices left after birth weight cleaning: ", nrow(mfr),"\n")
mfr %>% group_by(PARITET_5) %>% summarize(var(VEKT), sd(VEKT))


#### Previous delivery, mothers with only two pregnancies #### 
mfr = mfr %>% group_by(M_ID_1724) %>% filter(as.numeric(PARITET_5)<4) %>% filter(n() >= 2) 
table(mfr%>% summarise(n = n()) %>% pull(n))

mfr = mfr %>% group_by(M_ID_1724) %>% arrange(factor(PARITET_5, levels = c(0,1,2,3)) ) %>% slice_head(n = 2) %>%
  filter(PARITET_5[1] == "0",PARITET_5[2] == "1") %>% mutate(bw1 = bw_zscore[1], bw2 = bw_zscore[2])%>% 
  mutate(gd1 = SVLEN_DG[1], gd2 = SVLEN_DG[2]) %>% 
  mutate(ma1 = maternalage[1], ma2 = maternalage[2])%>% 
  mutate(sex1 = KJONN[1], sex2 = KJONN[2])

mfr = mfr %>% group_by(M_ID_1724) %>% arrange(factor(PARITET_5, levels = c(0,1,2,3)) ) %>% filter(row_number() == 2)
summary(mfr$bw1-mfr$bw2)


#### With genetic data #### 
if (motherorchild == "mother") {
  #### Mohter population #### 
  sentrix_m$M_ID_1724 = trimws(sentrix_m$M_ID_1724)
  sentrix_m = sentrix_m %>% select("M_ID_1724","SENTRIX_ID")
  mfr = inner_join(mfr, sentrix_m, by =c("M_ID_1724"))
  mfr = inner_join(mfr,qc_geno, by =c("SENTRIX_ID" = "SENTRIXID"))
} else if (motherorchild == "child") {
  #### Child population#### 
  sentrix_c$M_ID_1724 = trimws(sentrix_c$M_ID_1724)
  sentrix_c = sentrix_c %>% select("PREG_ID_1724","SENTRIX_ID")
  mfr = inner_join(mfr, sentrix_c, by =c("PREG_ID_1724"))
  mfr = inner_join(mfr, qc_geno, by =c("SENTRIX_ID" = "SENTRIXID"))
}


#### Dosage file #### 
dat_dosageall = fread(snakemake@input[[7]])
dat_dosage0 = fread(snakemake@input[[8]])
dat_dosage1 = fread(snakemake@input[[9]])

for_headerall = fread(snakemake@input[[10]])
for_header0 = fread(snakemake@input[[11]])
for_header1 = fread(snakemake@input[[12]])

listdosage = list(dat_dosage0,dat_dosage1,dat_dosageall)
listheader = list(for_header0,for_header1,for_headerall)

for (i in 1:length(listdosage)) {
  dat_dosage =  listdosage[[i]]
  for_header = listheader[[i]]
  dat_dosage$RSID = for_header %>% arrange(CHROM,GENPOS) %>% pull(nearestGene)
  colnam = dat_dosage$RSID
  dat_dosage = dat_dosage %>% select(-c("#CHROM","POS","REF","ALT","RSID"))
  dat_dosage = t(dat_dosage)
  colnames(dat_dosage) = colnam
  dat_dosage = as.data.frame(dat_dosage)
  dat_dosage$IID = row.names(dat_dosage)
  
  if (i == 1) {
    dosage = dat_dosage
    
  } else {
    unique_cols = setdiff(names(dat_dosage), names(dosage))
    unique = dat_dosage[, c("IID", unique_cols)]
    if (is.null(nrow(unique))) {
      cat("No new genetic variants to add","\n")
    }else {
      dosage = merge(dosage, unique, by = "IID", all = TRUE)}
  }
  
}

keep = as.data.frame(mfr$SENTRIX_ID)
colnames(keep) ="IID"
dat_dosage = inner_join(keep, dosage, by="IID")


#### Merge covar and phenotype, remove FID #### 
covar = fread(snakemake@input[[13]])
covar = covar %>% select(-c("FID","parity","sex","pregnancy_duration"))

mfr= mfr %>% ungroup() %>% select("SENTRIX_ID","bw1","bw2","gd1","gd2","ma1","ma2","sex1","sex2")

covar = inner_join(covar,mfr, by=c("IID" ="SENTRIX_ID"))



#### Linear models with and without prev bw #### 
## lm models, per snp
results = data.frame()
for (i in 1:(ncol(dat_dosage)-1)){
  cat("Regressions for: ",colnames(dat_dosage)[i+1],"\n")
  dat_snp = dat_dosage[,c(1,i+1)]
  rsid = colnames(dat_dosage)[i+1]
  covar.i = inner_join(covar, dat_snp, by="IID")
  
  formula = as.formula(
    paste(
      "bw2 ~ sex2 + genotyping_chip + ma2 + SVLEN_DG",
      "PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG +",
      "PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG + bw1 +",
      rsid
    )
  )
  
  fit = lm(formula, data = covar.i)
  
  variablesofinterest = c(paste(rsid), "bw1")
  coefficients = summary(fit)$coefficients
  tosave = as.data.frame(cbind(coefficients[variablesofinterest, c("Estimate")],confint(fit)[variablesofinterest,]))
  colnames(tosave) = c("Estimate","clmin","clmax")
  tosave$rsid = rsid
  tosave$n = length(fit$residuals)
  tosave$model = "With previous bw"
  results = rbind(results,tosave)
  
  formula = as.formula(
    paste(
      "bw2 ~ sex2 + genotyping_chip + ma2 + SVLEN_DG",
      "PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG +",
      "PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG +",
      rsid
    )
  )
  
  fit = lm(formula, data = covar.i)
  
  variablesofinterest = c(paste(rsid), "sex2")
  coefficients = summary(fit)$coefficients
  tosave = as.data.frame(cbind(coefficients[variablesofinterest, c("Estimate")],confint(fit)[variablesofinterest,]))
  colnames(tosave) = c("Estimate","clmin","clmax")
  tosave$rsid = rsid
  tosave$n = length(fit$residuals)
  tosave$model = "No previous bw"
  results = rbind(results,tosave)
}


#### Data for plotting #### 
toplot = results[grepl("^rs", rownames(results)), ]

x= toplot[toplot$model == "With previous bw",]
names(x) = paste0(names(x), "1")
y= toplot[toplot$model == "No previous bw",]
names(y) = paste0(names(y), "2")

df =cbind(x,y)

fwrite(df, snakemake@output[[1]], sep = "\t")


#### Prep data for first vs second pregnancy ####
if (motherorchild == "mother") {
  results = data.frame()
  for (i in 1:(ncol(dat_dosage)-1)){
    cat("Regressions for: ",colnames(dat_dosage)[i+1],"\n")
    dat_snp = dat_dosage[,c(1,i+1)]
    rsid = colnames(dat_dosage)[i+1]
    covar.i = inner_join(covar, dat_snp, by="IID")

    formula = as.formula(
      paste(
        "bw1 ~ sex1 + genotyping_chip + ma1 + SVLEN_DG",
        "PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG +",
        "PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG  +",
        rsid
      )
    )

    fit = lm(formula, data = covar.i,)

    variablesofinterest = c(paste(rsid), "sex1")
    coefficients = summary(fit)$coefficients
    tosave = as.data.frame(cbind(coefficients[variablesofinterest, c("Estimate")],confint(fit)[variablesofinterest,]))
    colnames(tosave) = c("Estimate","clmin","clmax")
    tosave$rsid = rsid
    tosave$n = length(fit$residuals)
    tosave$model = "Parity = zero"
    results = rbind(results,tosave)

    formula = as.formula(
      paste(
        "bw2 ~ sex2 + genotyping_chip + ma2 + SVLEN_DG",
        "PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG +",
        "PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG +",
        rsid
      )
    )

    fit = lm(formula, data = covar.i)

    variablesofinterest = c(paste(rsid), "sex2")
    coefficients = summary(fit)$coefficients
    tosave = as.data.frame(cbind(coefficients[variablesofinterest, c("Estimate")],confint(fit)[variablesofinterest,]))
    colnames(tosave) = c("Estimate","clmin","clmax")
    tosave$rsid = rsid
    tosave$n = length(fit$residuals)
    tosave$model = "Parity = one"
    results = rbind(results,tosave)

  }



  toplot = results[!grepl("sex", rownames(results)), ]

  x= toplot[toplot$model == "Parity = zero",]
  names(x) = paste0(names(x), "1")
  y= toplot[toplot$model == "Parity = one",]
  names(y) = paste0(names(y), "2")

  df.2 =cbind(x,y)

  fwrite(df.2, snakemake@output[[2]], sep = "\t")
  cat("An additional file was saved for plotting first vs second pregnancy SNP effects in mothers with two pregnancies in MoBa")
} else {
  df.2 = data.frame()
  fwrite(df.2, snakemake@output[[2]], sep = "\t")
}

