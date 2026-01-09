########### r-packages needed ###########
library(dplyr)
library(data.table)




########### Loading data needed ###########
dat = fread(snakemake@input[[1]])
nearestGene = fread(snakemake@input[[2]], select = c("CHROM","GENPOS","nearestGene"))




########### Finding top SNPs ###########

### Only keep SNPs reaching sig. level ###
dat_top = dat %>% filter(PVALUE<5*10**-8) %>% arrange(PVALUE) 
dat_top_alt = dat %>% filter(PVALUE<1*10**-6) %>% arrange(PVALUE)


### Extract most sig snp in independent regions for top loci with a radius of 1.5 Mb ###
findtopsnp = function(final_pos_list,chrom_list,dataframe,cutoff) {     # function to extract
for (chrom in chrom_list) {
	chromosome = filter(dataframe, CHROM==chrom)
	positions = chromosome$GENPOS
	min_P = positions[which.min(dataframe$PVALUE)]
	final_pos_list = bind_rows(final_pos_list, data.frame("GENPOS"=min_P,"CHROM"=chrom))
	for (pos in positions) {
		if ( !(pos %in% final_pos_list[,1] & chrom %in% final_pos_list[,2])){
			same_chrom = final_pos_list$GENPOS[final_pos_list$CHR==chrom]
			if(min(abs(pos - same_chrom)) > 1.5*10**6){
				final_pos_list = bind_rows(final_pos_list, data.frame("GENPOS"=pos,"CHROM"=chrom))  
			}
		}
	}
final_pos_list$cut_off = cutoff
}
return(final_pos_list)
}


final_pos_list1 = data.frame()
chrom_list1 = unique(dat_top$CHROM)

topsnp1 = findtopsnp(final_pos_list1, chrom_list1, dat_top,"5e-8") #PVALUE<5*10**-8, apply function


final_pos_list2 = data.frame()
chrom_list2 = unique(dat_top_alt$CHROM)

topsnp2 = findtopsnp(final_pos_list2, chrom_list2, dat_top_alt,"1e-6") #PVALUE<1*10**-6, apply function


aa = rbind(topsnp2,topsnp1) # combinde 

if (nrow(aa) == 0) {
	sig_snp = data.frame()
} else {
sig_snp = inner_join(dat_top_alt, aa, by = c("GENPOS","CHROM"))

nearestGene = unique( nearestGene[ , c("GENPOS","CHROM","nearestGene")])
sig_snp = left_join(sig_snp, nearestGene, by = c("GENPOS","CHROM")) #add nearest gene
}




########### Saving top SNPs ###########
fwrite(sig_snp, snakemake@output[[1]],sep =",")
