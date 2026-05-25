######## r-packages needed ###########
library(data.table)
library(dplyr)
library(ranger)
library(ggplot2)
library(mlr3verse)
library(progressr)
library(janitor)

########### Loading data needed ###########
input_file = snakemake@input[[1]]


dat = fread(input_file, drop = c("IID"))

# convert non-supported column names
names(dat) <- make.names(names(dat), unique = TRUE)                       


handlers("progress")
handlers(global=TRUE)


########### Running random forest models ###########
set.seed(1234)

target = "gd"

# Make regression task from data
task = as_task_regr(dat, target=target)

task$add_strata(cols=target, bins = 2L)

# Define search space
search_space = ps(
  regr.ranger.num.trees = p_int(200, 1000),                               
#  regr.ranger.mtry.ratio = p_dbl(0, 1),
  regr.ranger.mtry = p_int(1, task$n_features),                           
  regr.ranger.min.node.size = p_int(1, 40),
  regr.ranger.splitrule = p_fct(c("variance", "extratrees", "maxstat")),
  regr.ranger.replace = p_lgl(),                                        
  regr.ranger.sample.fraction = p_dbl(0.1, 1)
)


future::plan("multisession", workers = 6)
learner = lrn("regr.ranger")
learner$param_set$values = list(importance="impurity", respect.unordered.factors = "order" , num.threads=10) #### CHANGED
learner = as_learner(po("learner", learner=learner))                       


########### set log-level to warn
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")


# Set tuner
tuner = tnr("mbo")

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
      #model = learner$state$model$model
      model = learner$state$model$regr.ranger$model$model         
      if (is.null(model)) {
        print("Empty")
        return(NA)
      }
#      print(model$r.squared)
      model$r.squared
    }
  )
)

#print("starting tuning")
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
#  ),
#enable=TRUE,
                )


# best performing hyperparameter configuration
print(instance$result)

# clean param names
params <- gsub("regr.ranger.", "", names(instance$result$learner_param_vals[[1]]))
values = unname(instance$result$learner_param_vals[[1]])
best_grid <- setNames(values, params)

# save best grid
fwrite(as.data.frame(best_grid), snakemake@output[[1]])
print(best_grid)
rm(instance,MyMeasure,learner)
gc()

