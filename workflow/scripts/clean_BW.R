########### r-packages needed ###########
library(data.table)
library(dplyr)





########### Loading data needed ###########
dat = fread(snakemake@input[[1]])
gc_geno = fread(snakemake@input[[2]])





########### CLEANING PHENOFILE ###########
cat("Initial number of pregnanices in phenofile: ", nrow(dat),"\n")


#### No birth resulting in twins, triplet etc. ####
dat = dat %>% filter(is.na(plural_birth))
cat("Number of pregnanices left after removing multiple deliveries (twins, etc): ", nrow(dat),"\n")


#### Unknown sex ####
dat = dat %>% filter(sex != 0)
cat("Number of pregnanices left after removing unkown sex of child: ", nrow(dat),"\n")


#### Remove pregnancies treated with Asisted reproductive technology (ART) ####
dat = dat %>% filter(is.na(art))
cat("Number of pregnanices left after removing woman treated with ART: ", nrow(dat),"\n")


#### Gestational age ####
dat = dat %>% filter(pregnancy_duration >= 259 & pregnancy_duration < 301 & !is.na(pregnancy_duration))
cat("Number of pregnanices left after excluding pregnancies shorter than 259 days or  longer than 300 days: ", nrow(dat),"\n")


#### Remove children born with malformations ####
dat = dat %>% filter(is.na(congenital_malformations))
cat("Number of pregnanices left after removing children born with malformations: ", nrow(dat),"\n")


#### Remove stilbirths and perinatal death ####
dat = dat %>% filter(is.na(perinatal_stillborn) & is.na(early_perinatal_death) & is.na(perinatal_death))
cat("Number of pregnanices left after removing stillbirths and perinatal deaths: ", nrow(dat),"\n")


#### Create parity variable ####
dat = dat %>% filter(!is.na(n_previous_deliveries)) %>% mutate(parity = ifelse(grepl("0",n_previous_deliveries),0,1))
cat("Number of pregnanices left after removing uncertain parity: ", nrow(dat),"\n")


#### Birth weigth z-score #####
datbw = dat %>% filter(!is.na(weight_birth),weight_birth < mean(dat$weight_birth,na.rm=TRUE)+5*sd(dat$weight_birth,na.rm=TRUE), weight_birth > mean(dat$weight_birth,na.rm=TRUE)-5*sd(dat$weight_birth,na.rm=TRUE))

datbw$bw_zscore = (datbw$weight_birth - mean(datbw$weight_birth))/sd(datbw$weight_birth)
cat("Number of pregnanices left after birth weight cleaning: ", nrow(datbw),"\n")
datbw %>% group_by(parity) %>% summarize(var(weight_birth), sd(weight_birth))


#### Add genotyping chip, FID, IID ####
gc_geno = gc_geno %>% select(SENTRIXID, genotyping_chip, FID, IID)

dat_mother_bw = inner_join(datbw, gc_geno, by=c("mother_sentrix_id" = "SENTRIXID"))
dat_child_bw = inner_join(datbw, gc_geno, by=c("child_sentrix_id" = "SENTRIXID"))
dat_father_bw = inner_join(datbw, gc_geno, by=c("father_sentrix_id" = "SENTRIXID"))


cat("Only include the individuals passing genotype QC", "\n",
    "BW:" ,"\n" ,"mothers = ", nrow(dat_mother_bw), "\n", "fathers = ", nrow(dat_father_bw), "\n", "children = ", nrow(dat_child_bw),  "\n" )


#### One randrom preg per mother and father ####
set.seed(1234)
dat_mother_bw = dat_mother_bw %>% group_by(mother_id) %>% sample_n(1)

dat_father_bw = dat_father_bw %>% group_by(father_id) %>% sample_n(1)

cat("One random pregnancy per mother and father", "\n",
    "BW","\n","mothers = ", nrow(dat_mother_bw), "\n", "fathers = ", nrow(dat_father_bw), "\n")




########### Saving cleaned files ###########
fwrite(dat_mother_bw,snakemake@output[[1]])
fwrite(dat_child_bw,snakemake@output[[2]])
fwrite(dat_father_bw,snakemake@output[[3]])
