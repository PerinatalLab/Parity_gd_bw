########### r-packages needed ###########
library(data.table)
library(dplyr)
library(ranger)
library(ggplot2)
library(mlr3verse)
library(progressr)


########### Loading data needed ###########
input_file = snakemake@input[[1]]
dat = fread(input_file, drop = c("IID"))
phenotype = snakemake@wildcards[['phenotype']]

if (phenotype == "bw_zscore") {
	target = "bw_zscore"
} else{
	target = "gd"
}

### set up progress bar ###########
handlers("progress")
handlers(global=TRUE)

########### Running random forest models ###########
set.seed(1234)

# Make regression task from data

task = as_task_regr(dat, target=target) 
# stratify on target
task$add_strata(cols=target, bins = 2L)

 # Define search space

search_space = ps(
	regr.ranger.num.trees = p_int(20, 1000),
	regr.ranger.mtry.ratio = p_dbl(0, 1),
	regr.ranger.min.node.size = p_int(1, 40),
	regr.ranger.splitrule = p_fct(c("variance", "extratrees", "maxstat")),
	regr.ranger.replace = p_fct(c(TRUE, FALSE)),
	regr.ranger.sample.fraction = p_dbl(0.1, 1)
)


learner = lrn("regr.ranger")
learner$param_set$values = list(importance="impurity", respect.unordered.factors = "order" )

learner = as_learner(
			 po("learner", learner = learner)
			 )


# Set tuner
tuner = tnr("mbo")


########### set log-level to warn
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")


at =	auto_tuner(
			tuner = tuner,
			learner = learner,
			resampling = rsmp("cv", folds=5),
			search_space = search_space,
			measure = msr("regr.rsq"),
			term_evals = 60,
            store_models = FALSE,
			)

rsmp_cv = rsmp("cv", folds=10)

rr = #with_progress(
		resample(task, at, rsmp_cv, store_models=TRUE)#,
#		enable = TRUE,
#		)

print(rr$result)

# best performing hyperparameter configuration
# Initialize an empty data frame to store the results
results <- data.frame()

# Loop over each outer resampling iteration
for (i in seq_along(rr$learners)) { #resampling$iters)) {

	# Extract the tuning result for this iteration
	tuning_result <- rr$learners[[i]]$tuning_result
	# remove the 'regr.ranger.' prefix
  params <- gsub("regr.ranger.", "", names(tuning_result$learner_param_vals[[1]]))
	values <- unlist(tuning_result$learner_param_vals)
  hyperparameters <- setNames(values, params)

  # Extract the R-squared score
	r2_score <- tuning_result$regr.rsq[1]

	# Combine the fold number, R-squared score, and hyperparameters into a data frame
	result <- data.frame(fold = i, r2_score = r2_score)
    for (j in seq_along(hyperparameters)) {
    	result[[names(hyperparameters)[j]]] <- hyperparameters[[j]]
    	}

	# Add this result to the results data frame
	results <- rbind(results, result)
}

# Print the results as a table
print(results)

fwrite(results,snakemake@output[[1]])
r.squared2 = results$r2_score

print(r.squared2)

r2 = rbind(c("R2 median","5q", "95q"),
  c(quantile(r.squared2, p = 0.50), quantile(r.squared2, p = 0.025), quantile(r.squared2,, p = 0.975)))


### Variance explained, R2 ###
cat("Optimal r2: ", r2)

fwrite(r2, snakemake@output[[2]])
