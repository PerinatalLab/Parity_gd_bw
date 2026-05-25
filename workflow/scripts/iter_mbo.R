########## r-packages needed ###########
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(patchwork)
library(ranger)

best_grid = as.data.frame(fread(snakemake@input[[1]]))
dat = fread(snakemake@input[[2]], drop = c("IID"))
names(dat) <- make.names(names(dat), unique = TRUE)                       


#### Optimal model ####
if (phenotype == "bw_zscore") {
        target = "bw_zscore"
} else {
        target = "gd"
}
r.squared2 <- c()
tot_iter = 10000                                                      

## Create progress bar
pb = txtProgressBar(min = 0, max = tot_iter, style = 3)

# Create a formula from the target string
formula_str = paste(target, "~ .")
formula_obj = as.formula(formula_str)


for(i in 1:tot_iter){
  optimal_rf = ranger(
    formula = formula_obj,
    data = dat,
    num.trees = best_grid$num.trees,
 #   mtry = max(ceiling(best_grid$mtry.ratio * task$n_features),1),              
    mtry = best_grid$mtry,
    min.node.size = best_grid$min.node.size,
    splitrule = best_grid$splitrule,
    importance = "none",
    write.forest = TRUE,
    replace = best_grid$replace,
    sample.fraction = best_grid$sample.fraction,
    respect.unordered.factors = "order",
    oob.error = TRUE,       # OOB prediction errors
    num.threads = 29,     # CPUs
    save.memory = FALSE,    # use only if you encounter memory problems.
    verbose = TRUE
  )
  # ##  save each r2
  r.squared2[i] <- optimal_rf$r.squared
  # ##   Update the progress bar
  setTxtProgressBar(pb, i)
}

# ##Close the progress bar
close(pb)

r2 = rbind(c("R2 median","5q", "95q"),
           c(quantile(r.squared2, p = 0.50), quantile(r.squared2, p = 0.025), quantile(r.squared2,, p = 0.975)))


### Variance explained, R2 ###
cat("Optimal r2: ", r2)

fwrite(r2, snakemake@output[[1]])

