library(dplyr)
library(data.table)
library(ggplot2)

dat0 = fread(snakemake@input[[1]])
dat1 = fread(snakemake@input[[2]])
datall = fread(snakemake@input[[3]])


if (grepl("mother", snakemake@input[[1]])) {
  xname = "Gestational duration, days" 
} else {
  xname = "Gestational duration, days"
}


datall$group = ifelse(datall$IID %in% dat0$IID, "parity = zero", "parity > zero")

p = ggplot(datall, aes(x = gd)) +
  geom_histogram(
    aes(fill = group),
    bins = 150,
    alpha = 0.5,
    position = "identity"
  ) +
  scale_fill_manual(values=c("#7b3294", "#008837")) +
  geom_histogram(
    bins = 150,
    fill = "black",
    alpha = 0.2
  ) +
  xlab(paste(xname)) +
  ylab("Count") +
  theme_bw() +
  theme(panel.spacing.y = unit(20, "points"),
        axis.ticks.length.y = unit(0, "points"),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.line = element_line(),
        strip.text.x = element_text(size = 11, color="black", margin = margin(1, 1, 1, 1, "pt")),
        axis.title.x = element_text(size = 11, color="black"),
        axis.text.x = element_text(size = 10, color ="black"),
        axis.title.y = element_text(size = 11, color="black"),
        axis.text.y = element_text(size = 11.5, color ="black",  face = "italic"),
        plot.title = element_text(color="black", size=12, face="bold.italic"),
        legend.text = element_text(size=11),
        legend.title = element_text(size=11),
        legend.position = "bottom",
        text=element_text(family="Helvetica"),
        plot.margin = margin(r = 0.1, l = 0.1, t=0.05,unit = "cm"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_rect(color ="black", fill=NA, size=0.5),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(colour = "transparent", fill = "transparent")
        
  )

ggsave(snakemake@output[[1]],p,width = 190, height = 140, dpi = 1200, units = "mm")

