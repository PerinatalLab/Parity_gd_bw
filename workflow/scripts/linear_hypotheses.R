library(data.table)
library(dplyr)
library(tidyr)

if (grepl("gd/child-parity1",snakemake@input[[1]])) {
      ## No  SNP to analyse
      pp =  data.frame()
      fwrite(pp, snakemake@output[[1]], sep=",")
      fwrite(pp, snakemake@output[[2]], sep=",")
    } else {




format_haps= function(hap){
#variants= paste(hap$chr, hap$pos, hap$ref, hap$eff, sep =':')
variants = t(hap)[1,]  
ids= names(hap)[2:ncol(hap)]
hap= as.data.frame(t(hap[, 2:ncol(hap)]))
names(hap)= variants

if (grepl("mother",snakemake@input[[1]])) {
  hap$IID <- sub("^.*:", "", ids)
  
  
} else if (grepl("child",snakemake@input[[1]])) {
  hap$IID <- sub(":.*$", "", ids)
  
} 

return(hap)
}


h1= fread(snakemake@input[[1]])
h2= fread(snakemake@input[[2]])
h3= fread(snakemake@input[[3]])

h1= format_haps(h1)
h2= format_haps(h2)
h3= format_haps(h3)

if (grepl("bw_zscore", snakemake@input[[1]])) {
  pheno= fread(snakemake@input[[4]], select = c("IID","bw_zscore"))
} else if (grepl("haplotypes/gd", snakemake@input[[1]])) {
  pheno= fread(snakemake@input[[4]], select = c("IID","gd"))
}


covar_m= fread(snakemake@input[[5]])
covar_c= fread(snakemake@input[[6]])

if (grepl("mother",snakemake@input[[1]])) {
  covar_m= fread(snakemake@input[[5]], drop = c("FID"))
  covar_c= fread(snakemake@input[[6]], select = c("IID","parity"))
  colnames(covar_c) = c("IID","parity_c")
  duos = fread(snakemake@input[[7]])

  dat = inner_join(duos,covar_m, by = c("Mother"="IID"))
  dat = inner_join(dat,covar_c, by = c("Child"="IID"))
  dat = dat[dat$parity == dat$parity_c,]
  covar = dat %>% select(-c("parity_c","Child"))
  colnames(covar)[colnames(covar) == "Mother"] = "IID"

} else if (grepl("child",snakemake@input[[1]])) {
  covar = fread(snakemake@input[[6]])

}


pheno = inner_join(pheno, covar, by = c("IID"))

if (grepl("parity0",snakemake@input[[1]])) {
  pheno = pheno %>% filter(parity == 0)
  
} else if (grepl("parity1",snakemake@input[[1]])) {
  pheno = pheno %>% filter(parity == 1)
  
} 


print(nrow(pheno))
write( paste('snp', 'n', 'beta_h1', 'se_h1', 'pvalue_h1', 'beta_h2', 'se_h2', 'pvalue_h2', 'beta_h3', 'se_h3', 'pvalue_h3', 'beta_h4', 'se_h4', 'pvalue_h4', sep= '\t'), snakemake@output[[1]], append= T)

results_list= lapply(names(h1)[1:(length(names(h1))-1)], function(snp) {

print(snp)

h1_temp= h1[, c('IID', snp)]
h2_temp= h2[, c('IID', snp)]
h3_temp= h3[, c('IID', snp)]

names(h1_temp)= c('IID', 'h1')
names(h2_temp)= c('IID', 'h2')
names(h3_temp)= c('IID', 'h3')

d= inner_join(pheno, h1_temp, by= 'IID') %>% inner_join(., h2_temp, by= 'IID') %>% inner_join(., h3_temp, by= 'IID')

if (grepl("bw_zscore", snakemake@input[[1]])) {
  m1= lm(bw_zscore~ h1 + h2 + h3 + KJONN + maternalage + SVLEN_DG + genotyping_chip  + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG, data = d)
  
} else if (grepl("haplotypes/gd", snakemake@input[[1]])) {
  m1= lm(gd ~ h1 + h2 + h3 + KJONN + maternalage + genotyping_chip  + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG, data = d)
}

n= length(resid(m1))
coefs= summary(m1)$coefficients[2:5,]
beta_h1= coefs[1,1]
se_h1= coefs[1,2]
pvalue_h1= coefs[1,4]
beta_h2= coefs[2,1]
se_h2= coefs[2,2]
pvalue_h2= coefs[2,4]
beta_h3= coefs[3,1]
se_h3= coefs[3,2]
pvalue_h3= coefs[3,4]
beta_h4= coefs[4,1]
se_h4= coefs[4,2]
pvalue_h4= coefs[4,4]

results= paste(snp, n, beta_h1, se_h1, pvalue_h1, beta_h2, se_h2, pvalue_h2, beta_h3, se_h3, pvalue_h3, beta_h4, se_h4, pvalue_h4, sep= '\t')
write(results, file= snakemake@output[[1]], append=TRUE)

ci = data.frame(confint(m1))
ci$term = row.names(ci)
names(ci) = c('lo95', 'up95', 'term')
vars <- c('h1','h2','h3')
aa <- coefs[vars, , drop=FALSE]  # keep as data.frame
aa <- cbind(term = rownames(aa), aa)
aa <- merge(aa, ci, by='term')
aa$snp = snp


rm(h1_temp, h2_temp, h3_temp, d, m1, coefs)
gc()


aa
}

)

toplot <- bind_rows(results_list)

fwrite(toplot, snakemake@output[[2]], sep=",")
}
