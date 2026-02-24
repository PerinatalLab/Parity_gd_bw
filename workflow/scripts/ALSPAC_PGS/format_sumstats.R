library(optparse)
library(data.table)
library(dplyr)

parsed_opts <- list(
  make_option("--gwas",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--gwas_format",
              action = "store",
              default = NA,
              type = 'character')
  
)
if (interactive()) {
  opt <- list(gwas="/mnt/scratch/karin/Parity_gd_bw_pw/results/work/GWAS/gd/regenie/step2/QC/mother/parity0.txt",
              gwas_format="results/ALSPAC_PGS_mother_GD_{parity}/formated_sumstats.txt"
  )
} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
  
}

d= fread(opt$gwas,sep=",",select = c("ID", "BETA", "SE", "A1","A2")) %>%
  mutate(SNP=ID) %>%
  select(-ID) %>%
  select(SNP,A1,A2,BETA,SE)

d$A1= toupper(d$A1)
d$A2= toupper(d$A2)

print("any missing value")
print(colSums(is.na(d)))

print("SNP without rsid")
print(sum(grepl(':', d$SNP)))

print("BETA")
print(summary(d$BETA))

print("SE")
print(summary(d$SE))


d= na.omit(d)
gwas <- d %>%
  filter(!grepl(":", SNP), !is.na(SNP)) %>%
  filter(SE>0)




fwrite(gwas, opt$gwas_format, sep= '\t')

