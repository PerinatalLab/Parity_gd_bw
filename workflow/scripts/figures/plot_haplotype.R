library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggstance)


if (grepl("gd/lh/plot-child-parity1",snakemake@input[[1]])) {
      ## No  SNP to analyse
      pp =  ggplot() + theme_void() 
      ggsave(snakemake@output[[1]],pp, width = 174, height = 260, dpi = 100, units = "mm")
} else {



#plotdata = fread("/mnt/scratch/karin/for_review/results/effect_origin/bw_zscore/lh/plot-mother-parity0-results.csv")
#snpnames = fread("/mnt/scratch/karin/for_review/results/work/GWAS/bw_zscore/regenie/topSNP/mother/topSNP_parity0.txt")

plotdata = fread(snakemake@input[[1]])
snpnames = fread(snakemake@input[[2]])

snpnames = snpnames[snpnames$cut_off == 5e-08]
snpnames$snp= paste(snpnames$CHROM, snpnames$GENPOS, snpnames$A2, snpnames$A1, sep =':')
snpnames = snpnames %>% select(c("snp","nearestGene"))
plotdata = left_join(plotdata, snpnames, by = c("snp"))


if (grepl("mother", snakemake@input[[1]])) {
  if (grepl("gd/lh/", snakemake@input[[1]])) {
    if (grepl("parity0", snakemake@input[[1]])) {
      plotname = "Maternal genome, gestational duraiton, pregnancies of nulliparous woman"
      xname = "Days (95% CI)"
    } else if (grepl("parity1", snakemake@input[[1]])) {
      plotname = "Maternal genome, gestational duraiton, pregnancies of multiparous woman"
      xname = "Days (95% CI)"
    } else if (grepl("parityall", snakemake@input[[1]])) {
      plotname = "Maternal genome, gestational duraiton, whole study population"
      xname = "Days (95% CI)"
    }
  } else if (grepl("bw_zscore/lh/", snakemake@input[[1]])) {
    if (grepl("parity0", snakemake@input[[1]])) {
      plotname = "Maternal genome, birth weight, pregnancies of nulliparous woman"
      xname = "Z-score (95% CI)"
    } else if (grepl("parity1", snakemake@input[[1]])) {
      plotname = "Maternal genome, birth weight, pregnancies of multiparous woman"
      xname = "Z-score (95% CI)"
    } else if (grepl("parityall", snakemake@input[[1]])) {
      plotname = "Maternal genome, birth weight, whole study population"
      xname = "Z-score (95% CI)"
    }
  }
}
  
  
if (grepl("child", snakemake@input[[1]])) {
  if (grepl("gd/lh/", snakemake@input[[1]])) {
    if (grepl("parity0", snakemake@input[[1]])) {
      plotname = "Foetal genome, gestational duraiton, pregnancies of nulliparous woman"
      xname = "Days (95% CI)"
    } else if (grepl("parity1", snakemake@input[[1]])) {
      plotname = "Foetal genome, gestational duraiton, pregnancies of multiparous woman"
      xname = "Days (95% CI)"
    } else if (grepl("parityall", snakemake@input[[1]])) {
      plotname = "Foetal genome, gestational duraiton, whole study population"
      xname = "Days (95% CI)"
    }
  } else if (grepl("bw_zscore/lh/", snakemake@input[[1]])) {
    if (grepl("parity0", snakemake@input[[1]])) {
      plotname = "Foetal genome, birth weight, pregnancies of nulliparous woman"
      xname = "Z-score (95% CI)"
    } else if (grepl("parity1", snakemake@input[[1]])) {
      plotname = "Foetal genome, birth weight, pregnancies of multiparous woman"
      xname = "Z-score (95% CI)"
    } else if (grepl("parityall", snakemake@input[[1]])) {
      plotname = "Foetal genome, birth weight, whole study population"
      xname = "Z-score (95% CI)"
    }
  }
}

if (is.null(plotname) || plotname == "") {
  stop("Plotname is not set")
}

plotdata$term = factor(plotdata$term, levels= rev(c("h1", "h2", "h3")), labels= rev(c('Maternal\ntransmitted', 'Maternal\nnon-transmitted',
                                                                                  'Paternal\ntransmitted')))
 
pp = plotdata %>% ggplot(aes(x = Estimate, y = nearestGene, xmin = lo95, xmax = up95, group = term, fill = term)) +
  geom_vline(aes(xintercept = 0), linetype = 2, color="black", linewidth=0.6) +
  geom_errorbarh(aes(xmin = lo95, xmax = up95),size = .6, height = .2, color = "black", position = position_dodgev(height = 0.5)) +
  geom_point(size = 3, aes(shape = term, fill = term), position = position_dodgev(height = 0.5)) +
  ylab("") +
  xlab(paste(xname)) +
  ggtitle(paste(plotname)) +
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
        panel.grid.major = element_line(colour = "grey", linetype = "dotted", size = 0.7),
        panel.grid.minor = element_line(colour = "grey", linetype = "dotted", size = 0.7),
        panel.border = element_rect(color ="black", fill=NA, size=0.5),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(colour = "transparent", fill = "transparent")
        
  ) +
  scale_shape_manual(values = c(22, 21,23),guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5))+
  scale_fill_manual(values = c("#e66101", "#fdb863","#d7301f"), guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5))

if (grepl("child", snakemake@input[[1]])) {
ggsave(snakemake@output[[1]],pp,width = 174, height = 400, dpi = 600, units = "mm")}
 
ggsave(snakemake@output[[1]],pp,width = 174, height = 250, dpi = 600, units = "mm")

}
