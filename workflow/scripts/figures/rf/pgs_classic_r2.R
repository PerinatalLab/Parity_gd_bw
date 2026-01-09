########### r-packages needed ###########
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(patchwork)



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
                        gd_child = fread(snakemake@input[[i]]); gd_child$outcome = "GD"; gd_child$Genome = "Fetal"

                } else if (grepl("bw_zscore", snakemake@input[[i]])) {
                        bw_child = fread(snakemake@input[[i]]); bw_child$outcome = "BW"; bw_child$Genome = "Fetal"


                }
        }}



########### Merge to one dataset, prep for plotting ###########
aa = rbind(bw_mother,gd_mother,bw_child,gd_child)
colnames(aa) = c("estimate","se","tvalue","Pvalue","model","r2","n","variablemodel","outcome","Genome")

## Only keeping relevant estimates from model and extracting them for parity = 0 and parity > 0 ##
aa = aa %>% filter(variablemodel == "scale(as.numeric(SCORESUM), scale = F):as.factor(parity)1" | grepl("\\(Intercept\\)",variablemodel)& (model == "p0" | model == "p1") | grepl("scale\\(as\\.numeric\\(SCORESUM\\)\\, scale = F\\)",variablemodel)& (model == "p0" | model == "p1") 
)

plotdata = aa %>%
  mutate(
    interceptp0 = ifelse(grepl("\\(Intercept\\)",variablemodel) & model == "p0", estimate, NA),
    interceptp1 = ifelse(grepl("\\(Intercept\\)",variablemodel) & model == "p1", estimate, NA),
    slopep0 = ifelse(grepl("scale\\(as\\.numeric\\(SCORESUM\\)\\, scale = F\\)",variablemodel) & model == "p0", estimate, NA),
    slopep1 = ifelse(grepl("scale\\(as\\.numeric\\(SCORESUM\\)\\, scale = F\\)", variablemodel) & model == "p1", estimate, NA),
    r2_0 = ifelse(grepl("\\(Intercept\\)",variablemodel) & model == "p0", r2, NA),
    r2_1 = ifelse(grepl("\\(Intercept\\)",variablemodel) & model == "p1", r2, NA),
    interactionpvalue = ifelse(variablemodel == "scale(as.numeric(SCORESUM), scale = F):as.factor(parity)1", Pvalue, NA)
  ) %>%  group_by(outcome, Genome) %>% # Group by 'outcome' and 'Genome' to fill in values across rows that belong to the same outcome/Genome
  summarise(
    interceptp0 = first(na.omit(interceptp0)),
    interceptp1 = first(na.omit(interceptp1)),
    slopep0 = first(na.omit(slopep0)),
    slopep1 = first(na.omit(slopep1)),
    r2_0 = first(na.omit(r2_0)),
    r2_1 = first(na.omit(r2_1)),
    interactionpvalue = first(na.omit(interactionpvalue))
  ) %>%
  ungroup()

print(as.data.frame(plotdata))

## One sig. digits for p-value ##
sign_digits = function(x,d){
  s = format(x,digits=d)
  if(grepl("\\.", s) && ! grepl("e", s)) {
    n_sign_digits = nchar(s) -
      max( grepl("\\.", s), attr(regexpr("(^[-0.]*)", s), "match.length") )
    n_zeros = max(0, d - n_sign_digits)
    s = paste(s, paste(rep("0", n_zeros), collapse=""), sep="")
  }
  s
}

plotdata$interactionpvalue2 = sapply(plotdata$interactionpvalue,sign_digits,d=1)



########### Plotting regression lines ###########
## x-axis costume ##
custom_x_labels = data.frame(
  genome = c("Maternal", "Maternal", "Fetal", "Fetal"),
  outcome = c("GD","BW","GD","BW"),
  ylabel = c("Gestational duration","Z-score birth weight","Gestational duration","z-score birth weight")
)


## Function ##
create_facet_plot = function(data, genome_val, outcome_val,x_label) {
  # Subset data for the specific genome and outcome combination
  subset_data = data %>% filter(Genome == genome_val & outcome == outcome_val)
  
  # x and y limits
  if (outcome_val == "GD") {
	  final_range_x = c(-4,+4)
	  final_range_y = c(250,308)
  }else if(outcome_val == "BW") {
	  final_range_x = c(-4,+4)
	  final_range_y = c(-5,5)
  }
  
  # Plot
  p = subset_data %>% ggplot() +
	  coord_cartesian(xlim = final_range_x, ylim = final_range_y) +
          geom_hline(yintercept = 0, alpha = .5, color = "black") +
  	  geom_vline(xintercept = 0, alpha = .5, color = "black") +
	  xlab(paste("Polygenic score", paste(genome_val, "genome",sep =" "), sep="\n")) +
	  ylab(paste(outcome_val)) + 
	  geom_abline(intercept = subset_data$interceptp0, slope = subset_data$slopep0, color="#008837", linewidth=1.5) +
	  geom_abline(intercept = subset_data$interceptp1, slope = subset_data$slopep1, color="#7b3294", linewidth=1.5) +
	  scale_colour_manual(name = "Parity", values=c("= 0"="#008837","> 0"="#7b3294"))+
	  theme_bw() +
	  theme(
		legend.position = "none",
        	panel.grid.major = element_blank(),
        	panel.grid.minor = element_blank(),
	        axis.text = element_text(size=11),
        	axis.title = element_text(size=11, color="black"),
        	strip.background = element_blank(),
	        strip.placement = "outside",
        	strip.text.x = element_text(size = 11, color="black"),
	        panel.border = element_rect(color ="black", fill=NA, size=0.5),
        	panel.background = element_rect(fill = "transparent", colour = NA)) +
	  annotation_custom(grob = textGrob(paste0("PGS x Parity, pvalue = ",subset_data$interactionpvalue2 ), x = 1, y = -0.5,just = c("right", "top"), , gp = gpar(fontsize = 10)), xmin = Inf, xmax = Inf, ymin = Inf, ymax = Inf)

  return(p)
}

## List to store plots ##
facet_plots = list()

## Create plots for each combination of genome and outcome ##
for (g in unique(plotdata$Genome)) {
  for (o in unique(plotdata$outcome)) {
    # Get the custom x-axis label for the current combination of genome and outcome
    x_label = custom_x_labels$label[custom_x_labels$genome == g & custom_x_labels$outcome == o]
    facet_plots[[paste(g, o, sep = "_")]] = create_facet_plot(plotdata, g, o, x_label)
  }
}

## Combine  plots in specific order ##
order = c("Maternal_GD", "Maternal_BW",
  "Fetal_GD", "Fetal_BW")
ordered_plots = facet_plots[order]
combined_plot = wrap_plots(ordered_plots, ncol = 2) 



########### Save figure ###########
ggsave(snakemake@output[[1]],combined_plot,width = 190, height = 250, dpi = 1200, units = "mm")



