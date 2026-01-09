########### r-packages needed ###########
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)




########### Loading data needed ###########
d= fread(snakemake@input[[1]], h= T, select= c('ID', 'CHROM', 'GENPOS', 'PVALUE'))
d$pheno= "Parity 0"
x= fread(snakemake@input[[2]], h= T, select= c('ID', 'CHROM', 'GENPOS', 'PVALUE'))
x$pheno= "Parity > 0"

topsnp0 = fread(snakemake@input[[3]], select = c("ID", "cut_off","nearestGene"))[cut_off =="5e-08",]
topsnp1 = fread(snakemake@input[[4]], select = c("ID", "cut_off","nearestGene"))[cut_off =="5e-08",]

d=right_join(topsnp0, d, by = "ID")
x=right_join(topsnp1, x, by = "ID")

d= rbind(d, x)
rm(x)




########### Format data for manhattan plot ###########
dat= arrange(d, PVALUE)
rm(d)


don <- dat %>%
            group_by(CHROM) %>%
            summarise(chr_len= max(GENPOS)) %>%
            mutate(tot= cumsum(as.numeric(chr_len))-chr_len) %>% # Calculate cumulative position of each chromosome
            select(-chr_len) %>%
            left_join(dat, ., by= 'CHROM') %>%
            arrange(CHROM,GENPOS) %>% # Add a cumulative position of each SNP
            mutate(BPcum=GENPOS+tot) %>%
            ungroup()


axisdf = don %>%
        group_by(CHROM) %>%
        summarize(center=( max(BPcum) + min(BPcum) ) / 2 ) # Position of chrom number in x-axis
names(axisdf)= c('CHR', 'center')


HC= -log10(5*10**-8) # Significant threshold in GWAS


don = don %>%
        mutate(plotcolor = ifelse(CHROM %% 2 != 0 , "a","b")) %>%
	mutate(plotcolor = paste(plotcolor, pheno, sep="_")) %>% # Repeating 2 colours, different for each parity group 
        mutate(plotcolor = ifelse(GENPOS %in% topsnp0$GENPOS & CHROM %in% topsnp0$CHROM & pheno == "Parity 0", "d",plotcolor)) %>%
        mutate(plotcolor = ifelse(GENPOS %in% topsnp1$GENPOS & CHROM %in% topsnp1$CHROM & pheno == "Parity > 0", "d",plotcolor)) %>% # Colour top SNPs in an additional colour
	mutate(Gene = ifelse(is.na(nearestGene), "",nearestGene))




########### Manhattan plot ###########
# Generate breaks in steps of 2
max_pvalue = ceiling(max(-log10(don$PVALUE))) + 2
if (max_pvalue < 20) {
	breaks = seq(ceiling(max_pvalue / 2) * -2 , ceiling(max_pvalue / 2) * 2, by = 2)
} else {
	breaks = seq(ceiling(max_pvalue / 5) * -5 , ceiling(max_pvalue / 5) * 5, by = 5)
}


#Plot
don$logpval= with(don, ifelse(pheno == "Parity 0", -log10(PVALUE), log10(PVALUE)))
p1= ggplot(don  %>% arrange(plotcolor), aes(x= BPcum, y= logpval, color = plotcolor)) +
                  geom_point(size= 0.07) +   # Show all points
                  theme_bw( ) +
                  scale_colour_manual(values= c("#7b3294","#a6dba0","#c2a5cf","#008837","#ca0020"), guide= F) +
                  scale_x_continuous(expand= c(0, 0.3),label = c(1:19, '', 21,''), breaks= axisdf$center) + # label = ifelse(axisdf$CHR== 23, 'X', axisdf$CHR)
                  scale_y_continuous(expand= c(0, 0), limits= c(min(don$logpval) - 2, max(don$logpval) + 3.2), breaks= breaks, labels= abs(breaks)) +
	          ggtitle("a)")+
                  ylab('-log10(pvalue)') +
                  xlab('Chromosome') +
                  geom_hline(yintercept= -0.1, size= 0.6, colour= '#878787') +
                  geom_hline(yintercept= c(HC, -HC), size= 0.2, linetype= 2, colour= 'grey40') +
                  coord_cartesian(clip = "off") +
                  theme(legend.position= 'none',
                        plot.margin = unit(c(t= 0, r=0, b= 0, l=0), 'cm'),
                        text= element_text(family="arial", size= 11, colour = "black"),
                        axis.line= element_line(linewidth= 0.1, color ="black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
			axis.ticks = element_line(color = "black", linewidth = 0.5),
		        axis.text = element_text(color = "black", size = 9),
		        plot.title = element_text(hjust = 0))




########### For top SNP exists, add gene name of the nearest gene in manhattan plot ###########
text_position_top = max(don %>% filter(pheno == "Parity 0") %>% pull(logpval)) + 1.5
text_position_bottom = min(don %>% filter(pheno == "Parity > 0") %>% pull(logpval)) - 0.5

p1.final = p1 + geom_text_repel(dat = filter(don, Gene!= ""), aes(x = BPcum, y = logpval, label = Gene, fontface = "italic"),
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
                  ylim = c(-20, 20),
                  xlim = c(0.05, Inf),
                  arrow=arrow(angle=20, length=unit(0.02, "npc")) 
		  ) +
             annotate("text", x =c(150000000,150000000), y = c(text_position_top,text_position_bottom), label= c("Parity = 0", "Parity > 0"), size=8/.pt, color = "grey34" ) + 
	     theme(plot.margin = unit(c(t=0, r=0, b=0.5, l=0), 'cm'), 
		   panel.background = element_rect(fill = "transparent", colour = NA),
		   plot.background = element_rect(fill = "transparent", colour = NA)
	     )




########### Save figure ###########
ggsave(snakemake@output[[1]], p1.final, width = 15, height = 7, units="cm")
