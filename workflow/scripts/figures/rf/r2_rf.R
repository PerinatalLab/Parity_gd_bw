########### r-packages needed ###########
library(dplyr)
library(data.table)
library(ggplot2)
library(ggstance)



########### Loading data needed ###########
for (i in 1:length(snakemake@input)) {
        if (grepl("mother", snakemake@input[[i]])) {
		if (grepl("parity0", snakemake@input[[i]])) {
			if (grepl("gd", snakemake@input[[i]])) {
				gd_mother_p0 = fread(snakemake@input[[i]])
			       	gd_mother_p0$outcome = "GD"; gd_mother_p0$Genome = "Maternal"; gd_mother_p0$parity = "0" 
			} else if (grepl("bw_zscore", snakemake@input[[i]])) {
				bw_mother_p0 = fread(snakemake@input[[i]])
				bw_mother_p0$outcome = "BW"; bw_mother_p0$Genome = "Maternal"; bw_mother_p0$parity = "0"
			}
		}
		if (grepl("parity1", snakemake@input[[i]])) {
                        if (grepl("gd", snakemake@input[[i]])) {
                                gd_mother_p1 = fread(snakemake@input[[i]])
                                gd_mother_p1$outcome = "GD"; gd_mother_p1$Genome = "Maternal"; gd_mother_p1$parity = "> 0"
                        } else if (grepl("bw_zscore", snakemake@input[[i]])) {
                                bw_mother_p1 = fread(snakemake@input[[i]])
                                bw_mother_p1$outcome = "BW"; bw_mother_p1$Genome = "Maternal"; bw_mother_p1$parity = "> 0"
                        }
                }

        } else if (grepl("child", snakemake@input[[i]])) {
		if (grepl("parity0", snakemake@input[[i]])) {
			if (grepl("gd", snakemake@input[[i]])) {
                           	gd_child_p0 = fread(snakemake@input[[i]])
                               	gd_child_p0$outcome = "GD"; gd_child_p0$Genome = "Foetal"; gd_child_p0$parity = "0"
                       	} else if (grepl("bw_zscore", snakemake@input[[i]])) {
                             	bw_child_p0 = fread(snakemake@input[[i]])
                                bw_child_p0$outcome = "BW"; bw_child_p0$Genome = "Foetal"; bw_child_p0$parity = "0"
			}
		}
               	if (grepl("parity1", snakemake@input[[i]])) {
                      	if (grepl("gd", snakemake@input[[i]])) {
                              	gd_child_p1 = fread(snakemake@input[[i]])
                                gd_child_p1$outcome = "GD"; gd_child_p1$Genome = "Foetal"; gd_child_p1$parity = "> 0"
                        } else if (grepl("bw_zscore", snakemake@input[[i]])) {
                                bw_child_p1 = fread(snakemake@input[[i]])
                                bw_child_p1$outcome = "BW"; bw_child_p1$Genome = "Foetal"; bw_child_p1$parity = "> 0"
			}
		}
	}
}





########### Merge to one dataset, prep for plotting ###########
aa = rbind(bw_mother_p0[2,], gd_mother_p0[2,], bw_child_p0[2,], gd_child_p0[2,], bw_mother_p1[2,], gd_mother_p1[2,], bw_child_p1[2,], gd_child_p1[2,])

colnames(aa) = c("r2","min","max","Outcome","Genome","Parity")

aa$r2 = as.numeric(aa$r2)
aa$max = as.numeric(aa$max)
aa$min = as.numeric(aa$min)
aa$Parity = factor(aa$Parity, levels = c("> 0","0"))
aa$ordernames = paste(aa$Genome, aa$Outcome, sep = " ")
aa$ordernames = factor(aa$ordernames, levels = c("Foetal BW","Maternal BW","Foetal GD","Maternal GD"))




########### Plot ###########
p = aa %>%
  ggplot(aes(x = r2, y = ordernames,group = Parity,fill=Parity,shape = Parity)) +
  geom_errorbarh(aes(xmin = min, xmax = max),
                 size = 0.35, height = 0.1, color = "black",
                 position = position_dodgev(height = 0.55))+
 # geom_vline(aes(xintercept = 0), linetype = 2, colour="darkgrey", size=0.1) +
  geom_point(aes(shape = Parity, fill=Parity), size = 1.5,stroke = 0.4, position = position_dodgev(height = 0.55)) +
  ylab("") +
  xlab(bquote(r^2 ~ "(95% CI)")) +
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
        panel.grid.major.x = element_line(size=.1, color="darkgrey", linetype = 2), #element_blank(),
	panel.grid.major.y = element_blank(), #element_line(size=.1, color="darkgrey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color ="black", fill=NA, size=0.5),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(colour = "transparent", fill = "transparent")
        ) +
  #scale_x_continuous(breaks = seq(floor(min(aa$r2)), ceiling(max(aa$r2)),by = 0.005)) +
  scale_shape_manual(labels=c("> zero","zero"),values=c(22, 21), guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5)) +
  scale_fill_manual(values=c("#7b3294", "#008837"),labels=c("> zero","zero"),guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5))

if (max(aa$r2)>0.1 ) {
	p = p + scale_x_continuous(breaks = seq(floor(min(aa$r2)), ceiling(max(aa$r2)),by = 0.05)) 
}


########### Save figure ###########
ggsave(snakemake@output[[1]],p, width = 87, height = 87, dpi = 1200, units = "mm")

fwrite(aa, snakemake@output[[2]], sep= '\t')
