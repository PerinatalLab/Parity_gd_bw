########## r-packages needed ###########
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(patchwork)
library(ranger)

best_grid = as.data.frame(fread(snakemake@input[[1]]))
dat = fread(snakemake@input[[2]])
covariates = fread(snakemake@input[[3]])
names(dat) <- make.names(names(dat), unique = TRUE)                       


#### Optimal model ####
target = "gd"
r.squared2 <- c()
predicted.outcome = data.frame(matrix(, nrow=nrow(dat), ncol=0)) # Empty data frame, one row per individual
tot_iter = 10000                                                     

## Create progress bar
pb = txtProgressBar(min = 0, max = tot_iter, style = 3)

# Create a formula from the target string
formula_str = paste(target, "~ .")
formula_obj = as.formula(formula_str)

for(i in 1:tot_iter){
  optimal_rf = ranger(
    formula = formula_obj,
    data = select(dat, -IID),
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
    num.threads = 14,     # CPUs
    save.memory = FALSE,    # use only if you encounter memory problems.
    verbose = TRUE
  )
  # ##  save each r2
  r.squared2[i] <- optimal_rf$r.squared
  predicted.outcome = cbind(predicted.outcome, optimal_rf$predictions)

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






print("Testing of PGS x Parity")

dat2 = cbind(dat$IID, dat$gd, rowMeans(predicted.outcome))
print("A")
print(head(dat2))
colnames(dat2) = c("IID","outcome","predictedOutcome")
print("B")

interactiondat = inner_join(as.data.frame(dat2), covariates, by="IID")
print("C")



## linear model, PGSxparity
interactiondat$gd <- as.numeric(
  scale(as.numeric(interactiondat$predictedOutcome), scale = FALSE)
)
interactiondat$parity <- factor(interactiondat$parity)
interactiondat$outcome <- as.numeric(interactiondat$outcome)
print("D")

mod_gxe <- lm(outcome ~ gd * parity, data = interactiondat)

saveRDS(mod_gxe, file = snakemake@output[[2]])
print("E")


## Plot interaction
newdata <- expand.grid(
  gd = seq(
    min(interactiondat$gd),
    max(interactiondat$gd),length.out = 100
  ),
  parity = levels(interactiondat$parity)
)
print("F")

pred <- predict(
  mod_gxe,
  newdata = newdata,
  interval = "confidence"
)
print("G")

pred_df <- cbind(newdata, as.data.frame(pred))
print("H")

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

pred_df$interactionpvalue2 = sapply(as.numeric(summary(mod_gxe)$coefficients[4,4]),sign_digits,d=1)
print("I")

pred_df$parity <- factor(
  pred_df$parity,
  levels = c("0", "1"),
  labels = c("= 0", "> 0")
)
print("J")

pp = ggplot(pred_df, aes(x = gd, y = fit, color = parity, fill = parity)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3, color = NA) +
  geom_vline(xintercept = 0, alpha = .5, color = "black") +
  xlab(paste("Random forest genetic score", paste("Maternal genome"), sep="\n")) +
  ylab(paste("Gestational duration")) +
  scale_colour_manual(name = "Parity", values=c("0"="#008837","1"="#7b3294"))+
  scale_fill_manual(values = c("= 0" = "#008837", "> 0" = "#7b3294"))+
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
    plot.background = element_rect(fill = "transparent", colour = NA))+
  annotation_custom(
    grob = textGrob(paste0("PGS x Parity p-value = ", pred_df$interactionpvalue2),just = c("right", "top"), , gp = gpar(fontsize = 10)), xmin = Inf, xmax = Inf, ymin = Inf, ymax = Inf
  )
print("K")

ggsave(snakemake@output[[3]], pp, width = 87, height = 87, dpi = 1200, units = "mm", bg = "transparent")
