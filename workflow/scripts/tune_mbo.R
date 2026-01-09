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

handlers("progress")
handlers(global=TRUE)


########### Running random forest models ###########
set.seed(1234)

if (phenotype == "bw_zscore") {
	target = "bw_zscore"
} else {
	target = "gd"
}
# Make regression task from data
task = as_task_regr(dat, target=target) 

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

########### set log-level to warn
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")


# Set tuner
tuner = tnr("mbo")

# Use ranger oob-score as measure for hypertuning
MyMeasure = R6::R6Class(
	"MyMeasure",
	inherit = MeasureRegr,
	public = list(
		initialize = function() {
			super$initialize(
				id = "MyMeasure",
				range = c(-Inf, Inf),
				minimize = FALSE,
				predict_type = "response",
				properties = "requires_learner"
			)
		}
	),
	private = list(
		.score = function(prediction, learner, ...) {
			model = learner$state$model$regr.ranger
			if (is.null(model)) stop("Set store_models = TRUE.")
			model$model$r.squared
		}
	)
)


# hyperparameter tuning
instance = #with_progress(
			tune(
			  tuner = tuner,
			  task = task,
			  learner = learner,
			  resampling = rsmp("insample"), 
			  search_space = search_space,
			  measure = MyMeasure$new(),
			  term_evals = 60,
			  store_models = TRUE
			)#,
			#enable=TRUE,
#		)


# best performing hyperparameter configuration
print(instance$result)

# clean param names
params <- gsub("regr.ranger.", "", names(instance$result$learner_param_vals[[1]]))
values = unname(instance$result$learner_param_vals[[1]])
best_grid <- setNames(values, params)

# save best grid
fwrite(as.data.frame(best_grid), snakemake@output[[1]])



#### Optimal model ####

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
	  mtry.ratio = best_grid$mtry.ratio,            
	  min.node.size = best_grid$min.node.size,
	  splitrule = best_grid$splitrule,
	  importance = "none",    
	  write.forest = TRUE,
	  replace = best_grid$replace,
	  sample.fraction = best_grid$sample.fraction,
	  respect.unordered.factors = "order",
	  oob.error = TRUE,       # OOB prediction errors
	  num.threads = 16,     # CPUs
	  save.memory = FALSE,    # use only if you encounter memory problems.
	  verbose = TRUE         
	 	)
	# ##	save each r2
	r.squared2[i] <- optimal_rf$r.squared
	# ##	 Update the progress bar
	setTxtProgressBar(pb, i)
}

# ##Close the progress bar
close(pb)

r2 = rbind(c("R2 median","5q", "95q"),
  c(quantile(r.squared2, p = 0.50), quantile(r.squared2, p = 0.025), quantile(r.squared2,, p = 0.975)))


### Variance explained, R2 ###
cat("Optimal r2: ", r2)

fwrite(r2, snakemake@output[[2]])
