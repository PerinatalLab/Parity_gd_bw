library(optparse)
library(data.table)
library(dplyr)

parsed_opts <- list(
  make_option("--rdir",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--pdir",
              action = "store",
              default = NA,
              type = 'character')
  
)
if (interactive()) {
  opt <- list(rdir="results/replication/",
              pdir="/mnt/scratch/karin/topSNPs_Xiaoping/"
  )
} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
  
}


merge <- function(rdir,pdir,traits,genome,cohorts){
  
      gwas1 <- fread(paste0(rdir,traits,"_",genome,"_",cohorts,"/gwas.tsv",sep=""), header = TRUE)
      
      genome2 <- tolower(genome)
      gwas2 <- fread(paste0(pdir,traits,"/",cohorts,"_",genome2,".txt",sep=""), header = TRUE)
      
      merged <- left_join(gwas2, gwas1, by = c("CHROM", "GENPOS"), suffix = c(".d", ".r"))
      merged <- merged %>%
        mutate(
          allele_match = case_when(
            A1 == ALLELE1 & A2 == ALLELE0 ~ "same",
            A1 == ALLELE0 & A2 == ALLELE1 ~ "flip",
            TRUE ~ "mismatch"
          )) %>%
        mutate(BETA.r.harmonized = case_when(
            allele_match == "same" ~ BETA.r,
            allele_match == "flip" ~ -BETA.r,
            TRUE ~ NA_real_
          )) %>%
        select(-EXTRA) %>%
        mutate(same_sign= sign(BETA.d) == sign(BETA.r.harmonized))
      
      output <- paste0(rdir,traits,"_",genome,"_",cohorts,"/gwas_replication.tsv",sep="")
      write.table(merged,output,col.names=T,row.names=F,quote=F,sep="\t")
}


merge(rdir=opt$rdir,pdir=opt$pdir,traits="gd",genome="Mother",cohorts="parity0")
merge(opt$rdir,opt$pdir,traits="gd",genome="Mother",cohorts="parity1")
merge(opt$rdir,opt$pdir,traits="gd",genome="Mother",cohorts="parityall")
merge(opt$rdir,opt$pdir,traits="gd",genome="Child",cohorts="parity0")
merge(opt$rdir,opt$pdir,traits="gd",genome="Child",cohorts="parityall")


merge(opt$rdir,opt$pdir,traits="bw_zscore",genome="Mother",cohorts="parity0")
merge(opt$rdir,opt$pdir,traits="bw_zscore",genome="Mother",cohorts="parity1")
merge(opt$rdir,opt$pdir,traits="bw_zscore",genome="Mother",cohorts="parityall")
merge(opt$rdir,opt$pdir,traits="bw_zscore",genome="Child",cohorts="parity0")
merge(opt$rdir,opt$pdir,traits="bw_zscore",genome="Child",cohorts="parity1")
merge(opt$rdir,opt$pdir,traits="bw_zscore",genome="Child",cohorts="parityall")
