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



########### Merge to one dataset  ###########
aa = rbind(bw_mother,gd_mother,bw_child,gd_child)



########### x-axis costume ###########
custom_x_labels = data.frame(
  genome = c("Maternal", "Maternal", "Foetal", "Foetal"),
  outcome = c("GD","BW","GD","BW"),
  label = c("Days (95% CI)", "Z-score (95% CI)", "Days (95% CI)", "Z-score (95% CI)")
)



###########  Plot figure ###########

## Function ##
create_facet_plot = function(data, genome_val, outcome_val,x_label) {
  # Subset data for the specific genome and outcome combination
  subset_data = subset(data, genome == genome_val & outcome == outcome_val)

  # Equal x and y limits
  overall_range = range(
  c(subset_data$CImin, subset_data$CImax, subset_data$BETA, subset_data$CImin0, subset_data$CImax0, subset_data$BETA0), 
  na.rm = TRUE
)
  max_abs_range = max(abs(overall_range))
  final_range = c(-max_abs_range, max_abs_range)

  correlation = round(cor(subset_data$BETA0, subset_data$BETA, method = 'pearson', use="complete.obs"),digits =2)

  p1 = ggplot(subset_data, aes(x=BETA0, y=BETA)) +
        geom_abline(slope = 1,intercept =0, color = "grey40", linetype = "dashed") +
        geom_errorbar(aes(ymin = CImin, ymax = CImax),width =0) +
        geom_errorbarh(aes(xmin = CImin0, xmax = CImax0),height = 0) +
	geom_point(pch = 21, size = 3, fill = "#f1a340", color = "black") +
	ggtitle(paste(genome_val,outcome_val, sep ="\n") ) +
	xlab(paste("Effect size (Parity = zero)", x_label, sep="\n")) +
        ylab("Effect size (Parity > zero)") +
        geom_hline(yintercept = 0, alpha = .5, color = "black") +
        geom_vline(xintercept = 0, alpha = .5, color = "black") +
        coord_cartesian(xlim = final_range, ylim = final_range) +
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
	annotation_custom(
   	grob = textGrob(paste0("Pearson correlation = ", correlation), x = 1, y = -0.5,just = c("right", "top"), , gp = gpar(fontsize = 10)), xmin = Inf, xmax = Inf, ymin = Inf, ymax = Inf
	)

  return(p1)
}



## List to store plots #3
facet_plots = list()

## Create plots for each combination of genome and outcome ##
for (g in unique(aa$genome)) {
  for (o in unique(aa$outcome)) {
    # Get the custom x-axis label for the current combination of genome and outcome
    x_label = custom_x_labels$label[custom_x_labels$genome == g & custom_x_labels$outcome == o]
    facet_plots[[paste(g, o, sep = "_")]] = create_facet_plot(aa, g, o, x_label)
  }
}

## Combine  plots in specific order ##
order = c("Maternal_GD", "Maternal_BW",
  "Foetal_GD", "Foetal_BW")
ordered_plots = facet_plots[order]
combined_plot = wrap_plots(ordered_plots, ncol = 2) 



########### Save figure ###########
ggsave(snakemake@output[[1]],combined_plot,width = 190, height = 250, dpi = 1200, units = "mm")
