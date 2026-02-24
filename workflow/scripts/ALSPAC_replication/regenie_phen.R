library(optparse)
library(data.table)
library(dplyr)

parsed_opts <- list(
  make_option("--phen",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--ID",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--pca",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--traits",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--cohorts",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--phe",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--cov",
              action = "store",
              default = NA,
              type = 'character')
  
)
if (interactive()) {
  opt <- list(phen="/mnt/archive/alspac/pheno/B4346/alspac-B4346-phenotype.tsv",
              ID="results/replication/gd_Mother_parity0/ID.tsv",
              pca="results/replication/gd_Mother_parity0/ID_unrelated.eigenvec",
              traits="gd",
              cohorts="parity0",
              phe="results/replication/gd_Mother_parity0/phe.tsv",
              cov="results/replication/gd_Mother_parity0/cov.tsv"
  )
} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
  
}

phen <- fread(opt$phen) %>%
   mutate(PARITET=case_when(
    parity==0 ~ 0,
    parity>0  ~ 1,
    TRUE ~ NA)
  )%>%
  mutate(maternal_age2=maternal_age_at_delivery^2) %>%
  mutate(cidB4346= cidb4346) %>%
  mutate(GA=gestation_length*7,birth_weight=birthweight_from_notifications_or_clinical_records_grams) %>%
  select(cidB4346,PARITET,maternal_age_at_delivery,maternal_age2,sex_assigned_at_birth,GA,birth_weight) %>%
  filter(sex_assigned_at_birth==1 | sex_assigned_at_birth==2)

ID <- fread(opt$ID)
pca <- fread(opt$pca) %>%
  rename(FID = `#FID`)


traits <- opt$traits
cohorts <- opt$cohorts
  
if (traits=="gd" & cohorts=="parity0"){
  df <- pca %>%
    inner_join(ID,by=c("FID"="chipID")) %>%
    inner_join(phen,by="cidB4346") %>%
    filter(PARITET==0) %>%
    mutate(FID=imputedID,IID=imputedID,phe=GA) %>%
    select(FID,IID,phe,maternal_age_at_delivery,sex_assigned_at_birth,paste0("PC",1:10,sep="")) %>%
    filter(if_all(everything(), ~ !is.na(.))) %>%
    distinct()
    

  
}else if(traits=="gd" & cohorts=="parity1"){
  df <- pca %>%
    inner_join(ID,by=c("FID"="chipID")) %>%
    inner_join(phen,by="cidB4346") %>%
    filter(PARITET==1) %>%
    mutate(FID=imputedID,IID=imputedID,phe=GA) %>%
    select(FID,IID,phe,maternal_age_at_delivery,sex_assigned_at_birth,paste0("PC",1:10,sep="")) %>%
    filter(if_all(everything(), ~ !is.na(.))) %>%
    distinct()
  
  
}else if(traits=="gd" & cohorts=="parityall"){
  df <- pca %>%
    inner_join(ID,by=c("FID"="chipID")) %>%
    inner_join(phen,by="cidB4346") %>%
    mutate(FID=imputedID,IID=imputedID,phe=GA) %>%
    select(FID,IID,phe,PARITET,maternal_age_at_delivery,sex_assigned_at_birth,paste0("PC",1:10,sep="")) %>%
    filter(if_all(everything(), ~ !is.na(.))) %>%
    distinct() %>%
    filter(!duplicated(IID))
  
  
}else if(traits=="bw_zscore" & cohorts=="parity0"){
  df <- pca %>%
    inner_join(ID,by=c("FID"="chipID")) %>%
    inner_join(phen,by="cidB4346") %>%
    filter(PARITET==0) %>%
    mutate(bw_zscore=(birth_weight - mean(birth_weight, na.rm = TRUE)) / sd(birth_weight, na.rm = TRUE)) %>%
    mutate(FID=imputedID,IID=imputedID,phe=bw_zscore) %>%
    select(FID,IID,phe,maternal_age_at_delivery,sex_assigned_at_birth,GA,paste0("PC",1:10,sep="")) %>%
    filter(if_all(everything(), ~ !is.na(.))) %>%
    distinct()
  
}else if(traits=="bw_zscore" & cohorts=="parity1"){
  df <- pca %>%
    inner_join(ID,by=c("FID"="chipID")) %>%
    inner_join(phen,by="cidB4346") %>%
    filter(PARITET==1) %>%
    mutate(bw_zscore=(birth_weight - mean(birth_weight, na.rm = TRUE)) / sd(birth_weight, na.rm = TRUE)) %>%
    mutate(FID=imputedID,IID=imputedID,phe=bw_zscore) %>%
    select(FID,IID,phe,maternal_age_at_delivery,sex_assigned_at_birth,GA,paste0("PC",1:10,sep="")) %>%
    filter(if_all(everything(), ~ !is.na(.))) %>%
    distinct()
  
}else if(traits=="bw_zscore" & cohorts=="parityall"){
  df <- pca %>%
    inner_join(ID,by=c("FID"="chipID")) %>%
    inner_join(phen,by="cidB4346") %>%
    mutate(bw_zscore=(birth_weight - mean(birth_weight, na.rm = TRUE)) / sd(birth_weight, na.rm = TRUE)) %>%
    mutate(FID=imputedID,IID=imputedID,phe=bw_zscore) %>%
    select(FID,IID,phe,PARITET,maternal_age_at_delivery,sex_assigned_at_birth,GA,paste0("PC",1:10,sep="")) %>%
    filter(if_all(everything(), ~ !is.na(.))) %>%
    distinct() %>%
    filter(!duplicated(IID))
  
}else {
  print("Wrong trait name or parity name")
  print(traits)
  print(cohorts)
}

phe <- df %>%
  select(FID,IID,phe)
cov <- df %>%
  select(-phe) %>%
  arrange(PARITET) #make sure 0 start first

write.table(phe,opt$phe,col.names = T,row.names=F,quote=F,sep="\t")
write.table(cov,opt$cov,col.names = T,row.names=F,quote=F,sep="\t")
