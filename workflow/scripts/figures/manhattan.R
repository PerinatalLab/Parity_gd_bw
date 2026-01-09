########### r-packages needed ###########
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)




########### Loading data needed ###########
dat = fread(snakemake@input[[1]])
topsnp = fread(snakemake@input[[2]], select = c("ID", "cut_off","nearestGene"))




########### If top SNPs exist add gene name of the nearest gene from topsnp file ###########
if (nrow(topsnp) == 0){
	d=dat 
} else {
        topsnp = topsnp[cut_off =="5e-08",]
	d=right_join(topsnp, dat, by = "ID")
}




########### Format data for manhattan plot ###########
don <- d %>%
	    group_by(CHROM) %>% 
	    summarise(chr_len= max(GENPOS)) %>%
	    mutate(tot= cumsum(as.numeric(chr_len))-chr_len) %>% # Calculate cumulative position of each chromosome
	    select(-chr_len) %>%
	    left_join(d, ., by= 'CHROM') %>% 
	    arrange(CHROM,GENPOS) %>% # Add a cumulative position of each SNP
	    mutate(BPcum=GENPOS+tot) %>%
            ungroup()


axisdf = don %>% 
	group_by(CHROM) %>% 
	summarize(center=( max(BPcum) + min(BPcum) ) / 2 ) # Position of chrom number in x-axis
names(axisdf)= c('CHR', 'center')

HC= -log10(5*10**-8) # Significant threshold in GWAS




if (nrow(topsnp) == 0){
	don = don %>%
        mutate(plotcolor = ifelse(CHROM %% 2 != 0 , "a","b")) # If no top SNP exists, repeating 2 colors
} else {
don = don %>% 
	mutate(plotcolor = ifelse(CHROM %% 2 != 0 , "a","b")) %>% 
        mutate(plotcolor = ifelse(GENPOS %in% topsnp$GENPOS & CHROM %in% topsnp$CHROM & PVALUE< 5*10**-8, "c",plotcolor)) %>%
        mutate(Gene = ifelse(is.na(nearestGene), "",nearestGene)) # If top SNP exists, repeating 2 colours and colour top SNPs in an additional colour. 
}
 



########### Manhattan plot ###########
# Generate breaks in steps of 5
log_pvalues = -log10(don$PVALUE)
log_pvalues[is.infinite(log_pvalues)] = NA  # Replace Inf with NA
max_pvalue = max(log_pvalues, na.rm = TRUE) + 2
max_pvalue <- min(max_pvalue, 60)  # Cap the maximum value to avoid extremely large or infinite values
breaks = seq(0, ceiling(max_pvalue / 2) * 2, by = 2)

# Plot
p1= ggplot(don  %>% arrange(plotcolor), aes(x= BPcum, y= -log10(PVALUE), color = plotcolor)) +
	  geom_point(size= 0.07) +   # Show all points
	  theme_bw( ) +
          scale_colour_manual(values= c("#0571b0","#92c5de","#ca0020"), guide= F) +
	  scale_x_continuous(expand= c(0, 0.3),label = c(1:19, '', 21,''), breaks= axisdf$center) + 
	  scale_y_continuous(expand= c(0, 0.05), limits= c(0, max(-log10(don$PVALUE)) + 2), breaks = breaks, labels = abs(breaks)) + 
		  ylab('-log10(pvalue)') +
		  xlab('Chromosome') +
		  geom_hline(yintercept= 0, size= 0.25, colour= 'black') +
		  geom_hline(yintercept= c(HC), size= 0.2, linetype= 2, colour= '#878787') +
		  coord_cartesian(clip = "off") +
		  theme(legend.position= 'none',
			plot.margin = unit(c(t= 0, r=0, b= 0, l=0), 'cm'),
		       	text= element_text(family="arial", size= 9),
			axis.line= element_line(linewidth= 0.1),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
     	                panel.background = element_rect(fill = "transparent", colour = NA),
                        plot.background = element_rect(fill = "transparent", colour = NA))




########### If top SNP exists, add gene name of the nearest gene in manhattan plot ###########
if (nrow(topsnp) == 0){
	p2 = p1
} else {
p2 = p1+ geom_text_repel(dat = filter(don, Gene!= ""), aes(x= BPcum, y= -log10(PVALUE), label= Gene,fontface = "italic"),
                  size = 2.25,
                  force_pull = 0,
                  force = 0.1,
                  direction = "both",
                  hjust = 0,
                  vjust =  0.5,
		  box.padding = 0.1,
		  angle = 0,
                  segment.size = 0.3,
                  segment.square = TRUE,
                  segment.inflect = FALSE,
                  segment.colour = "black",
                  colour= "black",
                  segment.linetype = 4,
                  ylim = c(-Inf, 50),
                  xlim = c(0, Inf),
		  arrow=arrow(angle=20, length=unit(0.02, "npc")) ) 
}




########### Save figure ###########
ggsave(snakemake@output[[1]], p2, width = 15, height = 7, units="cm")


