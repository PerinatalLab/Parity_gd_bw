########## r-packages needed ###########
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(patchwork)
library(ranger)


########### Loading data needed ###########
dat = fread(snakemake@input[[1]])
best_grid = as.data.frame(fread(snakemake@input[[2]]))
covariates = fread(snakemake@input[[3]])

colnames(dat)[colnames(dat) %in% c("gd","bw_zscore")] = "outcome" # change column name to outcome




########### Loading "optimal model" ###########
set.seed(1234)

tot_iter = 10000
predicted.outcome = data.frame(matrix(, nrow=nrow(dat), ncol=0)) # Empty data frame, one row per individual

## Create progress bar
pb = txtProgressBar(min = 0, max = tot_iter, style = 3)

## Loop optimal_rf and extract predicted.outcome and feature importance
for(i in 1:tot_iter){
  optimal_rf = ranger(
    formula = outcome ~ .,
    data = select(dat, -IID),
    num.trees = best_grid$num.trees,
    mtry = best_grid$mtry.ratio,
    min.node.size = best_grid$min.node.size,
    splitrule = best_grid$splitrule,
    importance = "none",
    write.forest = TRUE,
    replace = best_grid$replace,
    sample.fraction = best_grid$sample.fraction,
    respect.unordered.factors = "order",
    oob.error = TRUE,       # OOB prediction errors
    num.threads = 16,       # CPUs
    save.memory = FALSE,    # Use only if you encounter memory problems
    verbose = TRUE
  )
  ## Save predicted outcome
  predicted.outcome = cbind(predicted.outcome, optimal_rf$predictions)
  
  ## Update the progress bar
  setTxtProgressBar(pb, i)
}

## Close the progress bar
close(pb)

## Predictions
cat("Predictions from optimal rf:",summary(rowMeans(predicted.outcome)))




########### GxE lm, for PGSxparity plot ###########
## Formating
print("Testing of PGS x Parity")
dat2 = cbind(dat$IID, dat$outcome, rowMeans(predicted.outcome))
colnames(dat2) = c("IID","outcome","predictedOutcome")

interactiondat = inner_join(as.data.frame(dat2), covariates, by="IID")



## linear model, PGSxparity
interactiondat$gd <- as.numeric(
  scale(as.numeric(interactiondat$predictedOutcome), scale = FALSE)
)
interactiondat$parity <- factor(interactiondat$parity)
interactiondat$outcome <- as.numeric(interactiondat$outcome)

mod_gxe <- lm(outcome ~ gd * parity, data = interactiondat)

saveRDS(mod_gxe, file = snakemake@output[[1]])


## Plot interaction
newdata <- expand.grid(
  gd = seq(
    min(interactiondat$gd),
    max(interactiondat$gd),length.out = 100
  ),
  parity = levels(interactiondat$parity)
)

pred <- predict(
  mod_gxe,
  newdata = newdata,
  interval = "confidence"
)

pred_df <- cbind(newdata, as.data.frame(pred))


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

pred_df$parity <- factor(
  pred_df$parity,
  levels = c("0", "1"),
  labels = c("= 0", "> 0")
)

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

ggsave(snakemake@output[[2]], pp, width = 87, height = 87, dpi = 1200, units = "mm", bg = "transparent")



