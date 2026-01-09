########### R-packages needed ###########
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggstance)



########### Loading data needed ###########
if (grepl("resources/gd", snakemake@input[[1]])) {
	aa = fread(snakemake@input[[1]]); aa$outcome = "GD"

} else if (grepl("bw_zscore", snakemake@input[[1]])) {
	aa = fread(snakemake@input[[1]]); aa$outcome = "BW"


}

if (grepl("mother", snakemake@input[[1]])) {
	aa$genome = "Maternal"
} else if (grepl("child", snakemake@input[[1]])) {
	aa$genome = "Foetal"
}



########### Order based on p-value and plotting (forest plot) ###########
if (length(unique(aa$ID)) == 0) {
	## No remaning SNP to plot for supplementary
	pp =  ggplot() + theme_void() 

	########### Save figure ###########
	ggsave(snakemake@output[[1]],pp, width = 174, height = 260, dpi = 1200, units = "mm")

} else if (length(unique(aa$ID)) < 40) {
	## Remaning SNPs fits in one plot
	aa$Parity = factor(aa$Parity, levels = c("Parity > 0","Parity 0"))
	
	new_order = aa %>%
		do(data_frame(order=levels(reorder(interaction(.$Gene, .$PVALUE, drop=TRUE), .$Gene)) )) %>%
		pull(order) # order based on p-value

	pp = aa %>%
  		mutate(order=factor(interaction(Gene, PVALUE), levels=rev(new_order))) %>%
  		ggplot(aes(x = BETA, y = order, xmin = CIlow, xmax = CIHigh, group = Parity, fill = Parity)) +
		geom_vline(aes(xintercept = 0), linetype = 2, color="darkgrey", linewidth=0.6) +
		geom_errorbarh(aes(xmin = CIlow, xmax = CIHigh),size = .6, height = .2, color = "black", position = position_dodgev(height = 0.5)) +
    	        geom_point(size = 3.3, aes(shape = Parity, fill = Parity), position = position_dodgev(height = 0.5)) +
  		ylab("") +
  		xlab("Beta (95% CI)") +
		ggtitle(paste(unique(aa$outcome), unique(aa$genome), sep=" ")) +
  		theme_bw() +
  		theme(panel.spacing.y = unit(10, "points"),
        		#axis.text.y = element_blank(),
        		axis.ticks.length.y = unit(0, "points"),
        		strip.background = element_blank(),
        		strip.placement = "outside",
        		axis.line = element_line(),
        		axis.title.x = element_text(size = 11, color="black"),
        		axis.text.x = element_text(size = 11, color ="black"),
		        axis.text.y = element_text(size = 11.5, color ="black",  face = "italic"),
			axis.title.y = element_blank(),
        		plot.title = element_text(color="black", size=14, face="bold.italic"),
        		legend.text = element_text(size=10.5),
        		legend.title = element_text(size=11),
        		legend.position = "bottom",
        		text=element_text(family="Helvetica"),
        		plot.margin = margin(r = 1, l = 1, t=0.05,unit = "cm"),
        		panel.grid.major = element_blank(),
        		panel.grid.minor = element_blank(),
        		panel.border = element_rect(color ="black", fill=NA, size=0.5),
			panel.background = element_rect(fill = "transparent", colour = NA),
			plot.background = element_rect(fill = "transparent", colour = NA),
			legend.key = element_rect(colour = "transparent", fill = "transparent"),
			legend.background = element_rect(colour = "transparent", fill = "transparent")

			) +
		scale_y_discrete(breaks= new_order, labels=gsub("\\..*", "", new_order),element_text(size = 8, color = "black", face = "italic")) +
  		geom_text(aes(x=Inf,label=as.numeric(PVALUE2)),size=3.5, hjust = 1, color="black") +
  		scale_shape_manual(labels=c("Parity > 0" ="> zero","Parity 0" = "zero"), values = c(22, 21),guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5))+
  		scale_fill_manual(labels=c("Parity > 0" ="> zero","Parity 0" = "zero"),values = c("#7b3294", "#008837"), guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5))

	########### Save figure ###########
	ggsave(snakemake@output[[1]],pp, width = 174, height = 260, dpi = 1200, units = "mm")

} else {
	## Divied into sub-plots 
	ng = 5 # nr. of groups (groups to plot)
	nr = nrow(aa)
	rows_per_group = ceiling(nr / ng) # rows per group on average
	breaks = c(quantile(aa$PVALUE, probs = seq(0, 1, by  = 1/ng)), Inf) # breaks for groups
	
	aa = aa %>%
  		mutate(subgroup = as.character(LETTERS[cut(PVALUE, breaks = breaks, labels = FALSE, include.lowest = TRUE)])) # Assign each row to a group based on P-value
        group = unique(aa$subgroup)
        
	aa$Parity = factor(aa$Parity, levels = c("Parity > 0","Parity 0"))
        	
	new_order = aa %>%
        group_by(genome) %>%
        do(data_frame(order=levels(reorder(interaction(.$Gene, .$PVALUE, drop=TRUE), .$Gene)) )) %>%
        pull(order) # order based on p-value

        pp = aa %>% 
               	mutate(order=factor(interaction(Gene, PVALUE), levels=rev(new_order))) %>%
  		ggplot(aes(x = BETA, y = order, xmin = CIlow, xmax = CIHigh, group = Parity, fill = Parity)) +
		geom_vline(aes(xintercept = 0), linetype = 2, color="darkgrey", linewidth=0.6) +
  		geom_errorbarh(aes(xmin = CIlow, xmax = CIHigh),size = .6, height = .2, color = "black", position = position_dodgev(height = 0.5)) +
		geom_point(size = 2.7, aes(shape = Parity, fill = Parity), position = position_dodgev(height = 0.5)) +
		facet_wrap(subgroup~., scales = "free",ncol = ng) +
               	ylab("") +
               	xlab("Beta (95% CI)") +
		ggtitle(paste(unique(aa$outcome), unique(aa$genome), sep=" ")) +
               	theme_bw() +
               	theme(panel.spacing.y = unit(12, "points"),
                       	axis.ticks.length.y = unit(0, "points"),
                       	strip.background = element_blank(),
                       	strip.placement = "outside",
			strip.text = element_text(size = 12, color="black"),
                       	axis.line = element_line(),
                       	axis.title.x = element_text(size = 12, color="black"),
                       	axis.text.x = element_text(size = 12, color ="black"),
			axis.title.y = element_blank(),
	                axis.text.y = element_text(size = 12.5, color ="black",  face = "italic"),
                       	plot.title = element_text(color="black", size=14, face="bold.italic"),
                       	legend.text = element_text(size=11.5),
                      	legend.title = element_text(size=12),
                       	legend.position = "bottom",
                       	text=element_text(family="Helvetica"),
                       	plot.margin = margin(r = 0.5, l = 0.5, t=0.05,unit = "cm"),
                       	panel.grid.major = element_blank(),
                       	panel.grid.minor = element_blank(),
                       	panel.border = element_rect(color ="black", fill=NA, size=0.5),
		        panel.background = element_rect(fill = "transparent", colour = NA),
		        plot.background = element_rect(fill = "transparent", colour = NA),
		        legend.key = element_rect(colour = "transparent", fill = "transparent"),
			legend.background = element_rect(colour = "transparent", fill = "transparent")
			) +
		scale_y_discrete(breaks= new_order, labels=gsub("\\..*", "", new_order),element_text(size = 9, color = "black", face = "italic")) +
  		geom_text(aes(x=Inf,label=as.numeric(PVALUE2)),size=5, hjust = 1, color="black") +
  		scale_shape_manual(labels=c("Parity > 0" ="> zero","Parity 0" = "zero"), values = c(22, 21),guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5))+
  		scale_fill_manual(labels=c("Parity > 0" ="> zero","Parity 0" = "zero"),values = c("#7b3294", "#008837"), guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5))

	########### Save figure ###########
	ggsave(snakemake@output[[1]],pp, width = 350, height = 260, dpi = 800, units = "mm")
}
