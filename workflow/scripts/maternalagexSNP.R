########### r-packages needed ##########
library(dplyr)
library(data.table)


########### Loading data ###########
dosage = fread(snakemake@input[[1]])
for_header = fread(snakemake@input[[2]],header = F)

covar = fread(snakemake@input[[3]])
covar = covar %>% select(-c("FID"))
pheno = fread(snakemake@input[[4]])
covar= inner_join(covar,pheno, by=c("IID"))


########### Format dosage ###########
colnames(dosage)[1]="chr"
dosage$chr = paste("CHROM",dosage$chr, sep="")
dosage$chom_pos = paste(dosage$chr, dosage$POS, sep="POS")
dosage = dosage %>% filter(ifelse(POS == 114296770 & ALT =="C",FALSE,TRUE))
names = as.data.frame(dosage[,dosage$chom_pos])
colnames(names) = "V1"
names = left_join(names,for_header, by= "V1")
dosage = dosage %>% select(-chr,-POS,-REF,-ALT,-chom_pos)
dosage = round(t(dosage))
colnames(dosage) = names[,"V2"]

dosage = as.data.frame(dosage)
colnames(dosage) = gsub("\\s+\\(", "(", colnames(dosage))
dosage$IID = row.names(dosage)

keep = as.data.frame(covar$IID)
colnames(keep) ="IID"
dat_dosage = inner_join(keep, dosage, by="IID")


########### Regressions  ###########

if (grepl("bw_zscore", snakemake@input[[3]])) {
  #### Linear models for bw #### 
  results = data.frame()
  cat("BW","\n")
  for (i in 1:(ncol(dat_dosage)-1)){
    cat("Regressions for: ",colnames(dat_dosage)[i+1],"\n")
    dat_snp = dat_dosage[,c(1,i+1)]
    locus = colnames(dat_dosage)[i+1]
    covar.i = inner_join(covar, dat_snp, by="IID")
    formula = as.formula(
      paste(
        "bw_zscore ~ KJONN+ genotyping_chip + maternalage + SVLEN_DG +",
        "PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG +",
        "PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG +`",
        locus,"`",sep=""
      )
    )
    
    fit = lm(formula, data = covar.i)
    
    if (grepl("NR2F2-AS1|SLC2A1-AS1|HLA-C|NKX6-3|KCNQ1|HMGA2|IGF1R|SLC6A2|ZBTB7B|LOC339593|KLHL29|LCORL|PDE10A|PTCH1",locus) & grepl("child", snakemake@input[[3]])){
      variablesofinterest = c(paste("`",locus,"`",sep=""), "maternalage")

    } else if (!grepl("NR2F2-AS1|SLC2A1-AS1|HLA-C|NKX6-3|KCNQ1|HMGA2|IGF1R|SLC6A2|ZBTB7B|LOC339593|KLHL29|LCORL|PDE10A|PTCH1",locus) & grepl("child", snakemake@input[[3]]))  {
      variablesofinterest = c(paste(locus), "maternalage")}  
    
    if (grepl("PLCE1|LINC00485|GLIS3",locus) & grepl("mother", snakemake@input[[3]])){
      variablesofinterest = c(paste("`",locus,"`",sep=""), "maternalage")
    } else if (!grepl("PLCE1|LINC00485|GLIS3",locus) & grepl("mother", snakemake@input[[3]])){
      variablesofinterest = c(paste(locus), "maternalage")}   
    
    coefficients = summary(fit)$coefficients

    tosave = as.data.frame(coefficients[variablesofinterest, c("Estimate", "Std. Error", "Pr(>|t|)")])
    tosave$locus = locus
    tosave$n = length(fit$residuals)
    tosave$model = "Main model"
    results = rbind(results,tosave)
    formula = as.formula(
      paste(
        "bw_zscore ~ KJONN + genotyping_chip + SVLEN_DG +",
        "PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG +",
        "PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG +maternalage*`",
        locus,"`",
        "+parity*`",
        locus,"`",sep=""
      )
    )
    
    fit = lm(formula, data = covar.i)

    if (grepl("NR2F2-AS1|SLC2A1-AS1|HLA-C|NKX6-3|KCNQ1|HMGA2|IGF1R|SLC6A2|ZBTB7B|LOC339593|KLHL29|LCORL|PDE10A|PTCH1",locus) & grepl("child", snakemake@input[[3]])){
      variablesofinterest = c(paste("`",locus,"`",sep=""), "maternalage",paste("maternalage:`",locus, "`",sep=""),paste("`",locus,"`:parity",sep=""))
    } else if (!grepl("NR2F2-AS1|SLC2A1-AS1|HLA-C|NKX6-3|KCNQ1|HMGA2|IGF1R|SLC6A2|ZBTB7B|LOC339593|KLHL29|LCORL|PDE10A|PTCH1",locus) & grepl("child", snakemake@input[[3]]))  {
      variablesofinterest = c(paste(locus), "maternalage",paste("maternalage:",locus, sep=""),paste(locus,":parity",sep=""))}  
    
    if (grepl("PLCE1|LINC00485|GLIS3",locus) & grepl("mother", snakemake@input[[3]])){
      variablesofinterest = c(paste("`",locus,"`",sep=""), "maternalage",paste("maternalage:`",locus, "`",sep=""),paste("`",locus,"`:parity",sep=""))
    } else if (!grepl("PLCE1|LINC00485|GLIS3",locus) & grepl("mother", snakemake@input[[3]])){
      variablesofinterest = c(paste(locus), "maternalage",paste("maternalage:",locus, sep=""),paste(locus,":parity",sep=""))}   
    
    coefficients = summary(fit)$coefficients
    tosave = as.data.frame(coefficients[variablesofinterest, c("Estimate", "Std. Error", "Pr(>|t|)")])
    tosave$locus = locus
    tosave$n = length(fit$residuals)
    tosave$model = "Parity x Maternal age"
    results = rbind(results,tosave)
    
  }
  
  
} else {
  #### interaction #### 
  results = data.frame()
  cat("GD","\n")
  for (i in 1:(ncol(dat_dosage)-1)){
    cat("Regressions for: ",colnames(dat_dosage)[i+1],"\n")
    dat_snp = dat_dosage[,c(1,i+1)]
    
    locus = colnames(dat_dosage)[i+1]
    covar.i = inner_join(covar, dat_snp, by="IID")
    formula = as.formula(
      paste(
        "gd ~ KJONN+ genotyping_chip + maternalage +",
        "PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG +",
        "PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG +`",
        locus,"`",sep=""
      )
    )
    
    fit = lm(formula, data = covar.i)
    
    if (grepl("EBF1|HLA-DQA1", locus) ){
    variablesofinterest = c(paste("`",locus,"`",sep=""), "maternalage")} else {
    variablesofinterest = c(paste(locus), "maternalage")}
    
    coefficients = summary(fit)$coefficients
    tosave = as.data.frame(coefficients[variablesofinterest, c("Estimate", "Std. Error", "Pr(>|t|)")])
    tosave$locus = locus
    tosave$n = length(fit$residuals)
    tosave$model = "Main model"
    results = rbind(results,tosave)
    formula = as.formula(
      paste(
        "gd ~ KJONN + genotyping_chip + maternalage +",
        "PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG +",
        "PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG +maternalage*`",
        locus,"`",
        "+parity*`",
        locus,"`",sep=""
      )
    )
    
    fit = lm(formula, data = covar.i)
    
    if (grepl("EBF1|HLA-DQA1", locus) ){
      variablesofinterest = c(paste("`",locus,"`",sep=""), "maternalage",paste("maternalage:`",locus, "`",sep=""),paste("`",locus,"`:parity",sep=""))
      } else {
        variablesofinterest = c(paste(locus), "maternalage",paste("maternalage:",locus, sep=""),paste(locus,":parity",sep=""))}
    
    coefficients = summary(fit)$coefficients
    tosave = as.data.frame(coefficients[variablesofinterest, c("Estimate", "Std. Error", "Pr(>|t|)")])
    tosave$locus = locus
    tosave$n = length(fit$residuals)
    tosave$model = "Parity x Maternal age"
    results = rbind(results,tosave)
    
  }
}

cat("rows result ",ncol(dat_dosage),"\n")

#### Data for saving #### 
x = results[grepl(paste0(c(colnames(dosage),"EBF1","KCNQ1","HMGA2","IGF1R","SLC6A2","ZBTB7B",
                            "LOC339593","KLHL29","LCORL","PDE10A","PTCH1","PLCE1","LINC00485","GLIS3"),
                          collapse = "|"),rownames(results)) & results$model == "Main model", ]
names(x) = paste0(names(x), " Main model")
cat("x",dim(x),"\n")

y = results[grepl(paste0(c(colnames(dosage),"EBF1","KCNQ1","HMGA2","IGF1R","SLC6A2","ZBTB7B",
                            "LOC339593","KLHL29","LCORL","PDE10A","PTCH1","PLCE1","LINC00485","GLIS3"),collapse = "|"),rownames(results)) & results$model == "Parity x Maternal age" & !grepl(":",rownames(results)), ]
rownames(y) = substr(rownames(y), 1, nchar(rownames(y)) - 1)
names(y) = paste0(names(y), "Parity x Maternal age")
cat("y",dim(y),"\n")


interactions_ma = results[grepl(paste0("maternalage:"),rownames(results)) ,]
names(interactions_ma) = paste0(names(interactions_ma), "_maternalage_interaction")

interactions_p = results[grepl(paste0(":parity"),rownames(results)),]
names(interactions_p) = paste0(names(interactions_p), "_parity_interaction")

df = cbind(x,y,interactions_ma,interactions_p)




#### Save data #### 
fwrite(df, snakemake@output[[1]], sep = ",")
#fwrite(df, "/mnt/scratch/karin/for_review/resources/forfigs/gd/mother_Xmaternalage.csv", sep = ",")

