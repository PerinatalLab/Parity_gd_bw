library(optparse)
library(data.table)
library(dplyr)

parsed_opts <- list(
  make_option("--phen",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--linkfile",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--genoID",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--outdir",
              action = "store",
              default = NA,
              type = 'character')
  
)
if (interactive()) {
  opt <- list(phen="/mnt/archive/alspac/pheno/B4346/alspac-B4346-phenotype.tsv",
              linkfile="/mnt/archive/alspac/pheno/B4346/alspac-B4346-linkage.tsv",
              genoID="/mnt/archive/alspac/alspac/B4346/imputed/data/swapped.sample",
              outdir="results/replication/gd_mother_parity0/ID.txt"
  )
} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
  
}

phen <- fread(opt$phen) %>%
  mutate(iflivein1y=ifelse(stillbirths==2 & outcome_b_livebirth == 1 & live_births == 3,TRUE,FALSE)) %>%
  mutate(ifhealth=ifelse(death_with_congenital_defect==-2 & termination_for_fetal_problem==7,TRUE,FALSE)) %>%
  mutate(ifnottwins=ifelse(pregnancy_size_summary == 1 & duplicate_pregnancy ==2,TRUE,FALSE)) %>%
  mutate(ifnotART=ifelse(b6s_in_vitro_fertilisation_conception != 1,TRUE,FALSE)) %>%
  mutate(ifspontaneous=ifelse(c4i_how_labour_started!=2 & c4i_how_labour_started!=3 & c4i_how_labour_started!=4 & c4i_how_labour_started!=-9999 ,TRUE,FALSE)) %>%
  mutate(ifnotCS=ifelse(c6a_method_of_delivery != 3 & c6c_delivery_by_caesarean_section != 1 & c6c_delivery_by_caesarean_section !=2,TRUE,FALSE)) %>%
  mutate(PARITET=case_when(
    parity==0 ~ 0,
    parity>0  ~ 1,
    TRUE ~ NA)
  )%>%
  mutate(maternal_age2=maternal_age_at_delivery^2) %>%
  mutate(cidB4346= cidb4346) %>%
  mutate(GA=gestation_length*7,birth_weight=birthweight_from_notifications_or_clinical_records_grams) %>%
  select(cidB4346,PARITET,maternal_age_at_delivery,maternal_age2,sex_assigned_at_birth,GA,birth_weight,iflivein1y,ifhealth,ifnottwins,ifnotART,ifspontaneous,ifnotCS) 

genoID<- fread(opt$genoID) %>%
  select(ID_2) 

linkfile <- fread(opt$linkfile)
linkfile$imputed_m <- paste0(linkfile$gi_hrc_g0m_g1,"M")
linkfile$imputed_f <- paste0(linkfile$gi_hrc_g0m_g1,"A")
linkfile$chip_m <- paste0(linkfile$gwa_660_g0m,"M")
linkfile$chip_f <- linkfile$gwa_550_g1

getID <- function(traits,cohorts,genome){
              if (traits=="gd" & cohorts=="parity0" & genome=="Mother"){
                ID <- phen %>%
                  filter(iflivein1y==T & ifhealth==T & ifnottwins==T & ifnotART==T & ifspontaneous==T & ifnotCS==T) %>%
                  filter(GA>154 & GA<308 & PARITET==0) %>%
                  inner_join(linkfile,by="cidB4346") %>%
                  inner_join(genoID,by=c("imputed_m"="ID_2")) %>%
                  mutate(imputedID=imputed_m,chipID=chip_m) %>%
                  select(cidB4346,imputedID,chipID) 
                
              }else if(traits=="gd" & cohorts=="parity1" & genome=="Mother"){
                ID <- phen %>%
                  filter(iflivein1y==T & ifhealth==T & ifnottwins==T & ifnotART==T & ifspontaneous==T & ifnotCS==T) %>%
                  filter(GA>154 & GA<308 & PARITET==1) %>%
                  inner_join(linkfile,by="cidB4346") %>%
                  inner_join(genoID,by=c("imputed_m"="ID_2")) %>%
                  mutate(imputedID=imputed_m,chipID=chip_m) %>%
                  select(cidB4346,imputedID,chipID)
                
              }else if(traits=="gd" & cohorts=="parityall" & genome=="Mother"){
                ID <- phen %>%
                  filter(iflivein1y==T & ifhealth==T & ifnottwins==T & ifnotART==T & ifspontaneous==T & ifnotCS==T) %>%
                  filter(GA>154 & GA<308) %>%
                  inner_join(linkfile,by="cidB4346",relationship = "many-to-many") %>%
                  inner_join(genoID,by=c("imputed_m"="ID_2")) %>%
                  mutate(imputedID=imputed_m,chipID=chip_m) %>%
                  select(cidB4346,imputedID,chipID)
                
              }else if(traits=="gd" & cohorts=="parity0" & genome=="Child"){
                ID <- phen %>%
                  filter(iflivein1y==T & ifhealth==T & ifnottwins==T & ifnotART==T & ifspontaneous==T & ifnotCS==T) %>%
                  filter(GA>154 & GA<308 & PARITET==0) %>%
                  inner_join(linkfile,by="cidB4346") %>%
                  inner_join(genoID,by=c("imputed_f"="ID_2")) %>%
                  mutate(imputedID=imputed_f,chipID=chip_f) %>%
                  select(cidB4346,imputedID,chipID)
              }else if(traits=="gd" & cohorts=="parity1" & genome=="Child"){
                ID <- phen %>%
                  filter(iflivein1y==T & ifhealth==T & ifnottwins==T & ifnotART==T & ifspontaneous==T & ifnotCS==T) %>%
                  filter(GA>154 & GA<308 & PARITET==1) %>%
                  inner_join(linkfile,by="cidB4346") %>%
                  inner_join(genoID,by=c("imputed_f"="ID_2")) %>%
                  mutate(imputedID=imputed_f,chipID=chip_f) %>%
                  select(cidB4346,imputedID,chipID) 
              }else if(traits=="gd" & cohorts=="parityall" & genome=="Child"){
                ID <- phen %>%
                  filter(iflivein1y==T & ifhealth==T & ifnottwins==T & ifnotART==T & ifspontaneous==T & ifnotCS==T) %>%
                  filter(GA>154 & GA<308) %>%
                  inner_join(linkfile,by="cidB4346",relationship = "many-to-many") %>%
                  inner_join(genoID,by=c("imputed_f"="ID_2")) %>%
                  mutate(imputedID=imputed_f,chipID=chip_f) %>%
                  select(cidB4346,imputedID,chipID) 
              }else if(traits=="bw_zscore" & cohorts=="parity0" & genome=="Mother"){
                ID <- phen %>%
                  filter(iflivein1y==T & ifhealth==T & ifnottwins==T & ifnotART==T) %>%
                  filter(GA>259 & GA<300 & birth_weight>0 & PARITET==0) %>%
                  filter(
                    abs(birth_weight - mean(birth_weight, na.rm = TRUE)) <= 5 * sd(birth_weight, na.rm = TRUE)
                  ) %>%
                  inner_join(linkfile,by="cidB4346") %>%
                  inner_join(genoID,by=c("imputed_m"="ID_2")) %>%
                  mutate(imputedID=imputed_m,chipID=chip_m) %>%
                  select(cidB4346,imputedID,chipID) 
                
              }else if(traits=="bw_zscore" & cohorts=="parity1" & genome=="Mother"){
                ID <- phen %>%
                  filter(iflivein1y==T & ifhealth==T & ifnottwins==T & ifnotART==T) %>%
                  filter(GA>259 & GA<300 & birth_weight>0 & PARITET==1) %>%
                  filter(
                    abs(birth_weight - mean(birth_weight, na.rm = TRUE)) <= 5 * sd(birth_weight, na.rm = TRUE)
                  ) %>%
                  inner_join(linkfile,by="cidB4346") %>%
                  inner_join(genoID,by=c("imputed_m"="ID_2")) %>%
                  mutate(imputedID=imputed_m,chipID=chip_m) %>%
                  select(cidB4346,imputedID,chipID)
                
              }else if(traits=="bw_zscore" & cohorts=="parityall" & genome=="Mother"){
                ID <- phen %>%
                  filter(iflivein1y==T & ifhealth==T & ifnottwins==T & ifnotART==T) %>%
                  filter(GA>259 & GA<300 & birth_weight>0) %>%
                  filter(
                    abs(birth_weight - mean(birth_weight, na.rm = TRUE)) <= 5 * sd(birth_weight, na.rm = TRUE)
                  ) %>%
                  inner_join(linkfile,by="cidB4346",relationship = "many-to-many") %>%
                  inner_join(genoID,by=c("imputed_m"="ID_2")) %>%
                  mutate(imputedID=imputed_m,chipID=chip_m) %>%
                  select(cidB4346,imputedID,chipID)
                
              }else if(traits=="bw_zscore" & cohorts=="parity0" & genome=="Child"){
                ID <- phen %>%
                  filter(iflivein1y==T & ifhealth==T & ifnottwins==T & ifnotART==T) %>%
                  filter(GA>259 & GA<300 & birth_weight>0 & PARITET==0) %>%
                  filter(
                    abs(birth_weight - mean(birth_weight, na.rm = TRUE)) <= 5 * sd(birth_weight, na.rm = TRUE)
                  ) %>%
                  inner_join(linkfile,by="cidB4346") %>%
                  inner_join(genoID,by=c("imputed_f"="ID_2")) %>%
                  mutate(imputedID=imputed_f,chipID=chip_f) %>%
                  select(cidB4346,imputedID,chipID) 
              }else if(traits=="bw_zscore" & cohorts=="parity1" & genome=="Child"){
                ID <- phen %>%
                  filter(iflivein1y==T & ifhealth==T & ifnottwins==T & ifnotART==T) %>%
                  filter(GA>259 & GA<300 & birth_weight>0 & PARITET==1) %>%
                  filter(
                    abs(birth_weight - mean(birth_weight, na.rm = TRUE)) <= 5 * sd(birth_weight, na.rm = TRUE)
                  ) %>%
                  inner_join(linkfile,by="cidB4346") %>%
                  inner_join(genoID,by=c("imputed_f"="ID_2")) %>%
                  mutate(imputedID=imputed_f,chipID=chip_f) %>%
                  select(cidB4346,imputedID,chipID) 
              }else if(traits=="bw_zscore" & cohorts=="parityall" & genome=="Child"){
                ID <- phen %>%
                  filter(iflivein1y==T & ifhealth==T & ifnottwins==T & ifnotART==T) %>%
                  filter(GA>259 & GA<300 & birth_weight>0) %>%
                  filter(
                    abs(birth_weight - mean(birth_weight, na.rm = TRUE)) <= 5 * sd(birth_weight, na.rm = TRUE)
                  ) %>%
                  inner_join(linkfile,by="cidB4346",relationship = "many-to-many") %>%
                  inner_join(genoID,by=c("imputed_f"="ID_2")) %>%
                  mutate(imputedID=imputed_f,chipID=chip_f) %>%
                  select(cidB4346,imputedID,chipID) 
              }else {
                print("Wrong trait name or parity name or wrong genome name!")
                print(traits)
                print(cohorts)
                print(genome)
              }
}

out <- getID(traits="gd",cohorts="parity0",genome="Mother")
write.table(out,paste0(opt$outdir,"gd_Mother_parity0/ID.tsv"),col.names = T,row.names=F,quote=F,sep="\t")
out <- getID(traits="gd",cohorts="parity1",genome="Mother")
write.table(out,paste0(opt$outdir,"gd_Mother_parity1/ID.tsv"),col.names = T,row.names=F,quote=F,sep="\t")
out <- getID(traits="gd",cohorts="parityall",genome="Mother")
write.table(out,paste0(opt$outdir,"gd_Mother_parityall/ID.tsv"),col.names = T,row.names=F,quote=F,sep="\t")

out <- getID(traits="gd",cohorts="parity0",genome="Child")
write.table(out,paste0(opt$outdir,"gd_Child_parity0/ID.tsv"),col.names = T,row.names=F,quote=F,sep="\t")
out <- getID(traits="gd",cohorts="parity1",genome="Child")
write.table(out,paste0(opt$outdir,"gd_Child_parity1/ID.tsv"),col.names = T,row.names=F,quote=F,sep="\t")
out <- getID(traits="gd",cohorts="parityall",genome="Child")
write.table(out,paste0(opt$outdir,"gd_Child_parityall/ID.tsv"),col.names = T,row.names=F,quote=F,sep="\t")


out <- getID(traits="bw_zscore",cohorts="parity0",genome="Mother")
write.table(out,paste0(opt$outdir,"bw_zscore_Mother_parity0/ID.tsv"),col.names = T,row.names=F,quote=F,sep="\t")
out <- getID(traits="bw_zscore",cohorts="parity1",genome="Mother")
write.table(out,paste0(opt$outdir,"bw_zscore_Mother_parity1/ID.tsv"),col.names = T,row.names=F,quote=F,sep="\t")
out <- getID(traits="bw_zscore",cohorts="parityall",genome="Mother")
write.table(out,paste0(opt$outdir,"bw_zscore_Mother_parityall/ID.tsv"),col.names = T,row.names=F,quote=F,sep="\t")

out <- getID(traits="bw_zscore",cohorts="parity0",genome="Child")
write.table(out,paste0(opt$outdir,"bw_zscore_Child_parity0/ID.tsv"),col.names = T,row.names=F,quote=F,sep="\t")
out <- getID(traits="bw_zscore",cohorts="parity1",genome="Child")
write.table(out,paste0(opt$outdir,"bw_zscore_Child_parity1/ID.tsv"),col.names = T,row.names=F,quote=F,sep="\t")
out <- getID(traits="bw_zscore",cohorts="parityall",genome="Child")
write.table(out,paste0(opt$outdir,"bw_zscore_Child_parityall/ID.tsv"),col.names = T,row.names=F,quote=F,sep="\t")
