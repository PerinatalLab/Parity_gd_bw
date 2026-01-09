########### r-packages needed ###########
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggstance)
library(grid)



########### Loading data needed ###########
for (i in 1:length(snakemake@input)) {
	if (grepl("mother", snakemake@input[[i]])) {
		if (grepl("resources/gd", snakemake@input[[i]])) {
			gd_mother = fread(snakemake@input[[i]]); gd_mother$outcome = "GD"; gd_mother$genome = "Maternal"

		} else if (grepl("bw_zscore", snakemake@input[[i]])) {
			bw_mother = fread(snakemake@input[[i]]); bw_mother$outcome = "BW"; bw_mother$genome = "Maternal"


		}
	} else if (grepl("child", snakemake@input[[i]])) {
		                if (grepl("resources/gd", snakemake@input[[i]])) {
                        gd_child = fread(snakemake@input[[i]]); gd_child$outcome = "GD"; gd_child$genome = "Foetal"

                } else if (grepl("bw_zscore", snakemake@input[[i]])) {
                        bw_child = fread(snakemake@input[[i]]); bw_child$outcome = "BW"; bw_child$genome = "Foetal"


                }
	}}



########### Merge to one dataset and order each genome and phenotype based on p-value ###########
aa = rbind(bw_mother,gd_mother,bw_child,gd_child)

aa$Parity = factor(aa$Parity, levels = c("Parity > 0","Parity 0"))
aa$outcome <- factor(aa$outcome, levels = c("GD", "BW"))
aa$genome <- factor(aa$genome, levels = c("Maternal","Foetal"))

new_order = aa %>%
  group_by(genome) %>%
  arrange(genome, outcome, PVALUE) %>%
  do(data_frame(order=levels(reorder(interaction(.$Gene,.$genome, .$PVALUE, drop=TRUE), .$Gene)) )) %>%
  pull(order)



########### Plotting, forest plot ###########
pp = aa %>%
  mutate(order=factor(interaction(Gene,genome, PVALUE), levels=rev(new_order))) %>%
  ggplot(aes(x = BETA, y = order, xmin = CIlow, xmax = CIHigh, group = Parity, fill = Parity)) +
  geom_vline(aes(xintercept = 0), linetype = 2, color="darkgrey", linewidth=0.6) +
  geom_errorbarh(aes(xmin = CIlow, xmax = CIHigh),size = .6, height = .2, color = "black", position = position_dodgev(height = 0.5)) +
  geom_point(size = 3.3, aes(shape = Parity, fill = Parity), position = position_dodgev(height = 0.5)) +
  facet_wrap(genome+outcome~.,scales ="free",ncol=2) +
  ylab("") +
  xlab("   Days (95% CI)                             Z-score (95% CI)") +
  theme_bw() +
  theme(panel.spacing.y = unit(20, "points"),
        axis.ticks.length.y = unit(0, "points"),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.line = element_line(),
        strip.text.x = element_text(size = 11, color="black", margin = margin(1, 1, 1, 1, "pt")),
        axis.title.x = element_text(size = 11, color="black"),
        axis.text.x = element_text(size = 10, color ="black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 11.5, color ="black",  face = "italic"),
        plot.title = element_text(color="black", size=12, face="bold.italic"),
        legend.text = element_text(size=11),
        legend.title = element_text(size=11),
        legend.position = "bottom",
        text=element_text(family="Helvetica"),
        plot.margin = margin(r = 0.1, l = 0.1, t=0.05,unit = "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color ="black", fill=NA, size=0.5),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(colour = "transparent", fill = "transparent")

  ) +
  scale_y_discrete(breaks= new_order, labels=gsub("\\..*", "", new_order),element_text(size = 7.5, color = "black", face = "italic")) +
  geom_text(aes(x=Inf,
	    label = as.numeric(PVALUE2), 
            fontface = ifelse(PVALUE2 <= 5e-04, "bold", "plain")),
	    size=4, hjust = 1, color="black") +
  scale_shape_manual(labels=c("Parity > 0" ="> zero","Parity 0" = "zero"), values = c(22, 21),guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5))+
  scale_fill_manual(labels=c("Parity > 0" ="> zero","Parity 0" = "zero"),values = c("#7b3294", "#008837"), guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5)) + coord_cartesian(clip = "off") 



########### Save figure ###########
ggsave(snakemake@output[[1]],pp,width = 190, height = 250, dpi = 1200, units = "mm")
