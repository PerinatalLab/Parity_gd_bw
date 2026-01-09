########### r-packages needed ###########
library(stringr)
library(data.table)




########### Function to extract h2 and genetic correlation from ldsc log files ###########
extract_ldsc_values = function(log_correlation, log_parityall) {
	# Read log files
	log_corr = readLines(log_correlation)
 	log_all = readLines(log_parityall)

  	# Variables to store the results
  	h2_1 = NA
  	se_1 = NA
  	h2_2 =  NA
  	se_2 = NA
  	gc = NA
  	se_gc = NA
  	h2_all = NA
  	se_all = NA

	# Find and extract the values
  	for (line in log_corr) {
		if (str_detect(line, "Total Observed scale h2:")) {
			if (is.na(h2_1)) {
				# h2 phenotype 1 (parity = 0)
        			h2_1 = as.numeric(str_extract(line, "\\d+\\.\\d+"))
        			se_1 = as.numeric(str_extract(line, "(?<=\\().+?(?=\\))"))
			} else {
        			# h2 phenotype 2 (parity > 0)
        			h2_2 = as.numeric(str_extract(line, "\\d+\\.\\d+"))
        			se_2 = as.numeric(str_extract(line, "(?<=\\().+?(?=\\))"))
			}
		}
		if (str_detect(line, "Genetic Correlation:")) {
			# Genetic correlation
      			gc = as.numeric(str_extract(line, "\\d+\\.\\d+"))
      			se_gc = as.numeric(str_extract(line, "(?<=\\().+?(?=\\))"))
		}
	}

	for (line in log_all) {
		if (str_detect(line, "Total Observed scale h2:")) {
			# h2 whole study population
			h2_all = as.numeric(str_extract(line, "\\d+\\.\\d+"))
	    		se_all = as.numeric(str_extract(line, "(?<=\\().+?(?=\\))"))
		}
	}

	# Dataframe for the results
	result = data.frame(c("0",">0","all"), c(h2_1,h2_2,h2_all), c(se_1, se_2, se_all), c(gc,gc,"-"), c(se_gc,se_gc,"-"))
	colnames(result) = c("Parity","h2","h2_se","genetic_correlation","genetic_correlation_se")
  	return(result)
}



########### Function for number of sig. digits for p-value ###########
sign_digits = function(x, d) {
  s = format(x, digits = d)
  if (grepl("\\.", s) && !grepl("e", s)) {
    n_sign_digits = nchar(s) -
      max(grepl("\\.", s), attr(regexpr("(^[-0.]*)", s), "match.length"))
    n_zeros = max(0, d - n_sign_digits)
    s = format(round(as.numeric(s), digits = d + n_zeros), nsmall = n_zeros)
  }
  s
}



########### Extract wanted values (h2 and genetic correlation) from log files from ldsc + CI calculated + sig. digits ###########
log_file_corr = snakemake@input[[1]]
log_file_all = snakemake@input[[2]]

ldsc_values = extract_ldsc_values(log_file_corr, log_file_all)

ldsc_values$CImin = ldsc_values$h2 - 1.96*ldsc_values$h2_se
ldsc_values$CImax = ldsc_values$h2 + 1.96*ldsc_values$h2_se

ldsc_values$genetic_correlation = sapply(ldsc_values$genetic_correlation, sign_digits, d=2)
ldsc_values$genetic_correlation_se = sapply(ldsc_values$genetic_correlation_se, sign_digits, d=2)



########### Save ldsc values  ###########
fwrite(ldsc_values, snakemake@output[[1]], sep=",")
