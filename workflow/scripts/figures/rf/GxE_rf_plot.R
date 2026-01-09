########## r-packages needed ###########
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
lmdata = list(bw_mother,gd_mother,bw_child,gd_child)

plotdata = c()

for (i in 1:length(lmdata)) {

interceptp0 = lmdata[[i]]$Estimate[1]
interceptp1 = lmdata[[i]]$Estimate[1] + lmdata[[i]]$Estimate[3]

slopep0 = lmdata[[i]]$Estimate[2]
slopep1 = lmdata[[i]]$Estimate[2] + lmdata[[i]]$Estimate[4]

interactionpvalue = lmdata[[i]]$`Pr(>|t|)`[4]

plotdata = rbind(plotdata, c(interceptp0,interceptp1,slopep0,slopep1,interactionpvalue,unique(lmdata[[i]]$outcome),unique(lmdata[[i]]$Genome)))

}

colnames(plotdata) = c("interceptp0","interceptp1","slopep0","slopep1","interactionpvalue","outcome","Genome")
plotdata=as.data.frame(plotdata)

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

plotdata$interactionpvalue2 = sapply(as.numeric(plotdata$interactionpvalue),sign_digits,d=1)



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
          final_range_x = c(-30,30)
          final_range_y = c(260,300)
	  yyy = 0.1 
  }else if(outcome_val == "BW") {
          final_range_x = c(-4,+4)
          final_range_y = c(-5,5)
	  yyy= 0.01
  }

  # Plot
  p = subset_data %>% ggplot() +
          coord_cartesian(xlim = final_range_x, ylim = final_range_y) +
          geom_hline(yintercept = 0, alpha = .5, color = "black") +
          geom_vline(xintercept = 0, alpha = .5, color = "black") +
          xlab(paste("Random forest genetic score", paste(genome_val, "genome",sep =" "), sep="\n")) +
          ylab(paste(outcome_val)) +
          geom_abline(intercept = as.numeric(subset_data$interceptp0), slope = as.numeric(subset_data$slopep0), color="#008837", linewidth=1.5) +
          geom_abline(intercept = as.numeric(subset_data$interceptp1), slope = as.numeric(subset_data$slopep1), color="#7b3294", linewidth=1.5) +
          scale_colour_manual(name = "Parity", values=c("= 0"="#008837","> 0"="#7b3294"))+
          theme_bw() +
          theme(
                legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text = element_text(size=11, color="black"),
                axis.title = element_text(size=11, color="black"),
                strip.background = element_blank(),
                strip.placement = "outside",
                strip.text.x = element_text(size = 11, color="black"),
		panel.border = element_rect(color ="black", fill=NA, size=0.5),
		panel.background = element_rect(fill = "transparent", colour = NA),
		plot.background = element_rect(fill = "transparent", colour = NA)) +
	geom_text( aes(x = max(final_range_x), y = max(final_range_y) - yyy), label = paste0("PGS x Parity p-value = ", subset_data$interactionpvalue2), hjust = 1, vjust = 1, size = 3)

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
ggsave(snakemake@output[[1]],combined_plot,width = 190, height = 250, dpi = 1200, units = "mm", bg = "transparent")
ggsave(snakemake@output[[2]],facet_plots[["Maternal_GD"]],width = 87, height = 87, dpi = 1200, units = "mm")
