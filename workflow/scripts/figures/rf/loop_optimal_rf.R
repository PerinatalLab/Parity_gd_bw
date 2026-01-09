########### r-packages needed ###########
library(data.table)
library(dplyr)
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
mod_gxe = lm(as.numeric(outcome) ~ scale(as.numeric(predictedOutcome),scale = F) * as.factor(parity), data = interactiondat)
print(summary(mod_gxe))
coeff = summary(mod_gxe)$coefficients




########### Save ###########
fwrite(as.data.frame(dat2), snakemake@output[[1]], col.names = TRUE)# Predicted outcome
fwrite(as.data.frame(coeff), snakemake@output[[2]], col.names = TRUE) # for GxE plot
