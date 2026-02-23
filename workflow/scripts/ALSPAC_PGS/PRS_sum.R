library(optparse)
library(data.table)
library(dplyr)

parsed_opts <- list(
  make_option("--PRS_chr",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--PRS",
              action = "store",
              default = NA,
              type = 'character')
  
)
if (interactive()) {
  opt <- list(PRS_chr="results/PGS_EUR_Mother_GA",
              PRS="results/PGS_EUR_Mother_GA/PGS.txt"
  )
} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
  
}

# List of PRS-CS output files for each chromosome
chr_files <- paste0(opt$PRS_chr,"/chr", 1:22, ".sscore")  # adjust names

# Initialize empty list
prs_list <- list()

for(f in chr_files){
  # Read file (adjust delimiter as needed)
  df <- fread(f)
  
  # Fix column name if needed
  if("#FID" %in% colnames(df)) colnames(df)[colnames(df)=="#FID"] <- "FID"
  
  # Detect PRS column
  prs_col <- intersect(c("SCORE1_SUM","SCORE1_AVG","NAMED_ALLELE_DOSAGE_SUM"), colnames(df))
  if(length(prs_col)==0) stop(paste("No PRS column found in", f))
  
  # If SCORE1_AVG, multiply by ALLELE_CT to get total PRS
  if("SCORE1_AVG" %in% prs_col & "ALLELE_CT" %in% colnames(df)){
    df$PRS_chr <- df$SCORE1_AVG * df$ALLELE_CT
  } else {
    df$PRS_chr <- df[[prs_col[1]]]
  }
  
  # Keep only FID, IID, PRS_chr
  prs_list[[f]] <- df %>% select(FID, IID, PRS_chr)
}

# Combine all chromosomes
all_prs <- bind_rows(prs_list)

# Sum across chromosomes for each individual
final_prs <- all_prs %>%
  group_by(FID, IID) %>%
  summarise(PRS = sum(PRS_chr), .groups = "drop")

# Save merged PRS
fwrite(final_prs, opt$PRS, sep= '\t')

