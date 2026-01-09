########### r-packages needed ###########
library(data.table)
library(dplyr)




########### Loading data needed ###########
mfr = fread(snakemake@input[[1]]) # The medical birth registry file
mfr_time_of_death = fread(snakemake@input[[2]]) #Information regarding if any pregnancy resulted in stillbirth (not existing in mfr)

parental_ids = fread(snakemake@input[[3]]) # File containing pregnancy ID, maternal ID and paternal ID

sentrix_mother = fread(snakemake@input[[4]]) # Sentrix file from genotyping - sample identifyer/positional information from array
sentrix_child = fread(snakemake@input[[5]]) # Sentrix file from genotyping - sample identifyer/positional information from array
sentrix_father = fread(snakemake@input[[6]]) # Sentrix file from genotyping - sample identifyer/positional information from array

qc_geno = fread(snakemake@input[[7]]) # Genotypes passing qc




########### CLEANING MFR ###########
cat("Initial number of pregnanices in medical birth registry (mfr): ", nrow(mfr),"\n")

#### Only spontaneous deliveries  ####
# FSTART = Initiation of delivery, if FSTART = 1 the delivery started spontaneous
mfr = mfr %>% filter(FSTART == 1 & !is.na(FSTART))
cat("Number of pregnanices left after only including spontaneous deliveries: ", nrow(mfr),"\n")


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
mfr = mfr %>% filter(SVLEN_DG >= 154 & SVLEN_DG <309 & !is.na(SVLEN_DG))
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

cat("Filtering mfr finished", "\n")


#### Calculating the maternal age when delivering child ####
mfr = mfr %>% mutate(maternalage = FAAR - MOR_FAAR)




########### Join sentrix IDs to mfr ###########

sentrix_mother$M_ID_1724 = trimws(sentrix_mother$M_ID_1724)
sentrix_child$M_ID_1724 = trimws(sentrix_child$M_ID_1724)
sentrix_father$M_ID_1724 = trimws(sentrix_father$M_ID_1724)


sentrix_mother = sentrix_mother %>% select("M_ID_1724","SENTRIX_ID") ; sentrix_father = sentrix_father %>% select("F_ID_1724","SENTRIX_ID") ; sentrix_child = sentrix_child %>% select("PREG_ID_1724","SENTRIX_ID")


mfr_mother = inner_join(mfr, sentrix_mother, by =c("M_ID_1724"))
mfr_child = inner_join(mfr, sentrix_child, by= c("PREG_ID_1724"))
mfr_father = inner_join(mfr, sentrix_father, by= c("F_ID_1724"))

cat("Adding sentrix ID to mfr", "\n", "Spliting data based on genotypes of mother, father, child", "\n", "mothers = ", nrow(mfr_mother), "\n", "fathers = ", nrow(mfr_father), "\n", "children = ", nrow(mfr_child),  "\n" )




########### Only include the individuals passing genotype QC ###########
mfr_mother = inner_join(mfr_mother,qc_geno, by =c("SENTRIX_ID" = "SENTRIXID"))
mfr_child = inner_join(mfr_child, qc_geno, by =c("SENTRIX_ID" = "SENTRIXID"))
mfr_father = inner_join(mfr_father, qc_geno, by =c("SENTRIX_ID" = "SENTRIXID"))

cat("Only include the individuals passing genotype QC", "\n","mothers = ", nrow(mfr_mother), "\n", "fathers = ", nrow(mfr_father), "\n", "children = ", nrow(mfr_child),  "\n" )




########### One randrom preg per mother and father ###########
set.seed(42)
rows <- sample(nrow(mfr_mother))
mfr_mother <- mfr_mother[rows, ]
mfr_mother = mfr_mother %>% group_by(M_ID_1724) %>% filter(row_number()==1)

set.seed(42)
rows <- sample(nrow(mfr_father))
mfr_father <- mfr_father[rows, ]
mfr_father = mfr_father %>% group_by(F_ID_1724) %>% filter(row_number()==1)

cat("One random pregnancy per mother and father", "\n","mothers = ", nrow(mfr_mother), "\n", "fathers = ", nrow(mfr_father), "\n", "children = ", nrow(mfr_child),  "\n" )




########### Saving cleaned files ###########
fwrite(mfr_mother,snakemake@output[[1]], sep=",")
fwrite(mfr_child,snakemake@output[[2]], sep=",")
fwrite(mfr_father,snakemake@output[[3]], sep=",")
