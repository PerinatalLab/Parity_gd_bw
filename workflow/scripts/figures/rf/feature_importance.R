########### r-packages needed ###########
library(data.table)
library(dplyr)
library(ranger)
library(ggplot2)




########### Loading data ###########
filename = snakemake@input[[1]]
dat = fread(snakemake@input[[1]])
best_grid = as.data.frame(fread(snakemake@input[[2]]))
axisnames = fread(snakemake@input[[3]], header=F)
colnames(dat)[colnames(dat) %in% c("gd","bw_zscore")] = "outcome" # change column name to outcome



########### Loading "optimal model" ###########
set.seed(1234)


## Extract feature importance
optimal_rf = ranger(
		    formula = outcome ~ .,
                    data = select(dat, -IID),
                    num.trees = best_grid$num.trees,
                    mtry = best_grid$mtry.ratio,
                    min.node.size = best_grid$min.node.size,
                    splitrule = best_grid$splitrule,
                    importance = "permutation",
                    write.forest = TRUE,
                    replace = best_grid$replace,
                    sample.fraction = best_grid$sample.fraction,
                    respect.unordered.factors = "order",
                    oob.error = TRUE,       # OOB prediction errors
                    num.threads = 16,       # CPUs
                    save.memory = FALSE,    # Use only if you encounter memory problems
                    verbose = TRUE
)



########### Plot feature importance ###########
## Format data
toplot = stack(optimal_rf$variable.importance)
toplot = toplot %>% arrange(ind)
colnames(axisnames) = c("ind","Yaxis")
toplot = inner_join(toplot, axisnames, by = "ind")
toplot$Yaxis = factor(toplot$Yaxis, levels = rev(toplot$Yaxis))


## Colour when plotting
if (grepl("parityall", filename, fixed = TRUE)) {
        col = "#0571b0"
} else if (grepl("parity0", filename, fixed = TRUE)) {
        col = "#008837"
} else if (grepl("parity1", filename, fixed = TRUE)) {
        col = "#7b3294" }


## Ploting
if (nrow(toplot) < 50) {
p = ggplot(toplot, aes(Yaxis, values)) + geom_col(fill=col) + coord_flip() +
        ylab("Feature importance (permutation)") +
        xlab("")+
        theme(text= element_text(family="arial", size= 12, colour = "black"),
                        axis.line= element_line(linewidth= 0.1, color ="black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.ticks = element_line(color = "black"),
                        axis.text = element_text(color = "black",size=12),
                        panel.background = element_rect(fill = "transparent", colour = NA),
                        plot.background = element_rect(fill = "transparent", colour = NA))

	ggsave(snakemake@output[[1]],p) # Saving figure
} else { 
	## Divied into sub-plots
        ng = 4 # nr. of groups (groups to plot)
        nr = nrow(toplot)
        rows_per_group = ceiling(nr / ng) # rows per group on average
        breaks = c(quantile(toplot$values, probs = seq(0, 1, by  = 1/ng)), Inf) # breaks for groups
        toplot = toplot %>%
		mutate(subgroup = as.character(LETTERS[cut(values, breaks = breaks, labels = FALSE, include.lowest = TRUE)])) # Assign each row to a group based on feature value
	group = unique(toplot$subgroup)
	
	p = ggplot(toplot, aes(Yaxis, values)) + geom_col(fill=col) + coord_flip() +
		ylab("Feature importance (permutation)") +
        	xlab("")+
                facet_wrap(subgroup~., scales = "free",ncol = ng) +
        	theme(text= element_text(family="arial", size= 12, colour = "black"),
                        axis.line= element_line(linewidth= 0.1, color ="black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.ticks = element_line(color = "black"),
                        axis.text.y = element_text(color = "black",size=12),
                        axis.text.x = element_text(color = "black",size=9),
                        panel.background = element_rect(fill = "transparent", colour = NA),
                        plot.background = element_rect(fill = "transparent", colour = NA))

	ggsave(snakemake@output[[1]],p,width = 400, height = 260, dpi = 800, units = "mm") # Saving figure
}
