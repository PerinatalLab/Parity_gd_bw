## random forest ##

rule check_results_all:
    'Check that all files in this snakefile are created'
    input:
        "results/checks/rf_mbo_performed.txt",
#        "results/checks/rf_cross_performed.txt"

    output:
        "results/checks/rf_allpgs_performed.txt"

    shell:
        "touch {output[0]}"



rule check_mbo_results:
    'Check that all files in tuning_mbo are created'
    input:
        expand("results/work/rf/tuning_hyperparameters/{feature}/grid_mbo_{phenotype}_{childormother}_{parity}.csv", feature = GE, phenotype = pheno_main, childormother = motherorchild,parity = parity),
        expand("results/work/rf/r2/{feature}/r2_mbo_{phenotype}_{childormother}_{parity}.csv", feature = GE, phenotype = pheno_main, childormother = motherorchild, parity = parity)

    output:
        "results/checks/rf_mbo_performed.txt"

    shell:
        "touch {output[0]}"



rule check_cross_results:
    'Check that all files in tuning_cross are created'
    input:
        expand("results/work/rf/tuning_hyperparameters/{feature}/grid_cross_{phenotype}_{childormother}_{parity}.csv", feature = GE, phenotype = pheno_main, childormother = motherorchild,parity = parity),
        expand("results/work/rf/r2/{feature}/r2_cross_{phenotype}_{childormother}_{parity}.csv", feature = GE, phenotype = pheno_main, childormother = motherorchild, parity = parity)

    output:
        "results/checks/rf_cross_performed.txt"

    shell:
        "touch {output[0]}"		



rule tuning_mbo:
    'Perform hyperparameter tuning using bayesian optimization'
    input:
        "results/work/rf/input_{feature}/{phenotype}/{childormother}/{parity}_input_rf.csv"

    output:
        "results/work/rf/tuning_hyperparameters/{feature}/grid_mbo_{phenotype}_{childormother}_{parity}.csv"

    conda:
        "../envs/rf.yml"
	
    script:
        "../scripts/tune_mbo.R"


rule tuning_mbo:
    'Run the random forest model with the tuned hyperparameters (for 10 000 iterations) to evaluate the model performance by calculating the r2 on oob predictions'
    input:
        "results/work/rf/tuning_hyperparameters/{feature}/grid_mbo_{phenotype}_{childormother}_{parity}.csv",
        "results/work/rf/input_{feature}/{phenotype}/{childormother}/{parity}_input_rf.csv"

    output:
        "results/work/rf/r2/{feature}/r2_mbo_{phenotype}_{childormother}_{parity}.csv"

    conda:
        "../envs/rf.yml"

    script:
        "../scripts/iter_mbo.R"

rule tuning_cross:
    'Use nested cross-validation to tune hyperparameters (inner fold) and evaluate model performance (outer fold)'
    input:
        "results/work/rf/input_{feature}/{phenotype}/{childormother}/{parity}_input_rf.csv"

    output:
        "results/work/rf/tuning_hyperparameters/{feature}/grid_cross_{phenotype}_{childormother}_{parity}.csv",
        "results/work/rf/r2/{feature}/r2_cross_{phenotype}_{childormother}_{parity}.csv"

    conda:
        "../envs/rf.yml"

    script:
        "../scripts/tune_cross.R"

