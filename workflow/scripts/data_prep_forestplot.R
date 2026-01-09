########### r-packages needed ###########
library(dplyr)
library(data.table)


########### Loading data ###########
p0 = fread(snakemake@input[[1]], select =c("ID","BETA","SE"))
p1 = fread(snakemake@input[[2]], select =c("ID","BETA","SE"))
interaction_pvalue = fread(snakemake@input[[3]], select =c("ID","PVALUE"))
topsnp = fread(snakemake@input[[4]], select = c("RSID", "Gene"))
LDbuddies = fread(snakemake@input[[5]])



########### Combind topsnps from mother and child, LDbuddies r2>0.3 for variants missing in our summary statistics ###########
LDbuddies = LDbuddies %>% filter(LD_r2 > 0.3)

match_indices = match(topsnp$RSID, LDbuddies$Original_snp)
topsnp$RSID[!is.na(match_indices)] = LDbuddies$Highest_LD_buddie[match_indices[!is.na(match_indices)]]

topsnp = topsnp %>% distinct()



##### Extrac topsnps from summary statistics #####
p0 = inner_join(p0, topsnp, by =c("ID" = "RSID"))
p1 = inner_join(p1, topsnp, by =c("ID" = "RSID"))

p0$Parity = "Parity 0"; p1$Parity = "Parity > 0"

dat = rbind(p0, p1)



##### Add p-value from interaction #####
dat = inner_join(interaction_pvalue,dat, by="ID")



##### CI #####
dat$CIlow = dat$BETA -1.96*dat$SE; dat$CIHigh = dat$BETA + 1.96*dat$SE



##### One sig. digits for p-value #####
sign_digits = function(x,d){
  s = format(x,digits=d)
  if(grepl("\\.", s) && ! grepl("e", s)) {
    n_sign_digits = nchar(s) -
      max( grepl("\\.", s), attr(regexpr("(^[-0.]*)", s), "match.length") )
    n_zeros = max(0, d - n_sign_digits)
    s = paste(s, paste(rep("0", n_zeros), collapse=""), sep="")
  }
  s
}

dat$PVALUE2 = sapply(dat$PVALUE,sign_digits,d=1)



#### To main- or suppfig ####
top_rows = 10

if (grepl("GWAS/gd_noPTDs", snakemake@input[[3]])) {
	tomainfig = dat %>%
                mutate(rank = rank(PVALUE)) %>%
                filter(rank <= top_rows) %>%
                ungroup() %>%
                select(-rank)
        fwrite(tomainfig, snakemake@output[[1]], sep=",")
	fwrite(dat, snakemake@output[[2]], sep=",")

} else {
	tomainfig = dat %>%
	       	mutate(rank = rank(PVALUE)) %>%
  		filter(rank <= top_rows) %>%
  		ungroup() %>%
  		select(-rank)


	tosupfig = dat %>%
  		mutate(rank = rank(PVALUE)) %>%
  		filter(rank > top_rows) %>%
  		ungroup() %>%
  		select(-rank)

	fwrite(tomainfig, snakemake@output[[1]], sep=",")
	fwrite(tosupfig, snakemake@output[[2]], sep=",") }
