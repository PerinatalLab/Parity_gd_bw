########### r-packages needed ###########
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggstance)



########### Loading data needed ###########
for (i in 1:length(snakemake@input)) {
        if (grepl("mother", snakemake@input[[i]])) {
                if (grepl("resources/gd", snakemake@input[[i]])) {
                        gd_mother = fread(snakemake@input[[i]]); gd_mother$outcome = "GD"; gd_mother$Genome = "Maternal"

                } else if (grepl("bw_zscore", snakemake@input[[i]])) {
                        bw_mother = fread(snakemake@input[[i]]); bw_mother$outcome = "BW"; bw_mother$Genome = "Maternal"


                }
        } else if (grepl("child", snakemake@input[[i]])) {
                                if (grepl("resources/gd", snakemake@input[[i]])) {
                        gd_child = fread(snakemake@input[[i]]); gd_child$outcome = "GD"; gd_child$Genome = "Foetal"

                } else if (grepl("bw_zscore", snakemake@input[[i]])) {
                        bw_child = fread(snakemake@input[[i]]); bw_child$outcome = "BW"; bw_child$Genome = "Foetal"


                }
        }}



########### Merge to one dataset, prep for plotting ###########
aa = rbind(bw_mother,gd_mother,bw_child,gd_child)

aa$gc = paste(aa$genetic_correlation, " (",aa$genetic_correlation_se,")", sep="") #genetic correlation (95% CI)

aa = aa %>% filter(Parity != "all") 

aa$Parity = factor(aa$Parity, levels = c(">0","0"))

aa$ordernames = paste(aa$Genome, aa$outcome, sep = " ")

aa$ordernames = factor(aa$ordernames, levels = c("Foetal BW","Maternal BW","Foetal GD","Maternal GD"))



########### Plot ###########
p = aa %>%
  ggplot(aes(x = h2, y = ordernames,group = Parity,fill=Parity,shape = Parity)) +
  geom_errorbarh(aes(xmin = CImin, xmax = CImax), 
                 size = 0.35, height = 0.1, color = "black", 
                 position = position_dodgev(height = 0.55))+
  geom_vline(aes(xintercept = 0), linetype = 2, colour="darkgrey", size=0.1) +
  geom_point(aes(shape = Parity, fill=Parity), size = 1.5,stroke = 0.4, position = position_dodgev(height = 0.55)) +
  ylab("") +
  xlab(bquote(h^2 ~ "(95% CI)")) +
  theme_bw() +
  theme(panel.spacing.y = unit(5, "points"),
        axis.ticks.length.y = unit(0, "points"),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.line = element_line(),
        strip.text.x = element_text(size = 8, color="black"),
        axis.title.x = element_text(size = 8, color="black"),
        axis.text.x = element_text(size = 8, color ="black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8, color ="black",  face = "italic"),
        plot.title = element_text(color="black", size=6, face="bold.italic"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "bottom",
        text=element_text(family="Helvetica"),
        plot.margin = margin(r = 1, l = 1, t=0.05,unit = "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color ="black", fill=NA, size=0.5),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(colour = "transparent", fill = "transparent"),
	plot.margin = margin(20, 50, 5.5, 5.5)
        ) +
  scale_shape_manual(labels=c("Multiparous","Nulliparous"),values=c(22, 21), guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5)) +
  scale_fill_manual(values=c("#7b3294", "#008837"),labels=c("Multiparous","Nulliparous"),guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5))+
  geom_text(aes(x=Inf,label=gc),size=2.5, hjust = -0.2, color="black") +
  coord_cartesian(clip = "off") +
  geom_text(
    x = Inf,
    y = Inf,
    label = "Genetic\ncorrelation",
    hjust = -0.1,
    vjust = 0,
    size = 2.5,
    fontface = "bold"
  ) + geom_text(
    x = -Inf,
    y = Inf,
    label = expression(h^2),
    hjust = -1, 
    vjust = -0.5,  
    size = 3,
    fontface = "bold"
  )



########### Save figure ###########
ggsave(snakemake@output[[1]],p, width = 87, height = 87, dpi = 1200, units = "mm")
