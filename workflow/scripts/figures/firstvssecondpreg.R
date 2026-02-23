########### r-packages needed ###########
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)

df.2 = fread(snakemake@input[[1]])

df.2 = df.2 %>% mutate(color = ifelse(abs(Estimate1)>abs(Estimate2) ,0,1)) # Colour based on effect direction 

if (grepl("resources/forfigs/gd", snakemake@input[[1]])) {
  xname= "Days (95% CI)"
  xlimit = c(-3.5,0.6)
  ylimit = c(-3.5,0.6)

  } else {
  xname = "Z-score (95% CI)"
  xlimit = c(-0.1,0.1)
  ylimit = c(-0.1,0.1)
}


correlation = round(cor(df.2$Estimate1, df.2$Estimate2, method = 'pearson', use="complete.obs"),digits =2)
P2 = ggplot(df.2, aes(x=Estimate1, y=Estimate2,color=as.factor(color))) +
  xlab(paste("Effect size in first pregnancy", xname,sep="\n")) +
  ylab(paste("Effect size in second pregancy", xname,sep="\n")) +
  geom_abline(slope = 1,intercept =0, color = "red", linetype = "dashed")  +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = seq(-3.5,3.5, by = 0.5)) +
  scale_y_continuous(breaks = seq(-3.5,3.5, by = 0.5)) +
  coord_cartesian(xlim = xlimit, ylim = ylimit)+
  geom_errorbar(aes(ymin = clmin2, ymax = clmax2), color = "grey40",linewidth = 0.3,width =0) +
  geom_errorbarh(aes(xmin = clmin1, xmax = clmax1), color = "grey40",linewidth = 0.3,height = 0) +
  geom_point(pch = 19, size = 2) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=12),
    axis.title = element_text(size=12, color="black"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(size = 12, color="black"),
    panel.border = element_rect(color ="black", fill=NA, size=0.5),
    panel.background = element_rect(fill = "transparent", colour = NA)) +
  annotation_custom(
    grob = textGrob(paste0("Pearson correlation = ", correlation), x = 1, y = -0.5,just = c("right", "top"), , gp = gpar(fontsize = 10)), xmin = Inf, xmax = Inf, ymin = Inf, ymax = Inf
  ) + 
  scale_color_manual(
    values=c("#008837", "#7b3294"))
cat("Largets absolute effect in:","\n",table(df.2$color),"\n","\n")
cat("Number of SNPs significant in parity zero (0.05/18)","\n",table(df.2$`Pr(>|t|)1` <0.05/18),"\n","\n")
cat("Number of SNPs significant in parity one (0.05/18)","\n",table(df.2$`Pr(>|t|)2` <0.05/18),"\n","\n")

ggsave(snakemake@output[[1]],P2, dpi = 1200)
