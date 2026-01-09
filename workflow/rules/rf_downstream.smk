##### Evaluating the PGS created with random forest #####

rule plot_r2:
    'Plotting the r2 for the random forest model with the tuned hyperparameters'
    input:
        expand("results/work/rf/r2/{{GorGE}}/r2_mbo_{pheno}_{motherorchild}_{parity}.csv", pheno = pheno_main, motherorchild = motherorchild, parity = parity_0_or_1)

    output:
        "results/figures/rf/r2_{GorGE}_plot.png",
	"results/figures/rf/r2_{GorGE}_plot.txt"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/rf/r2_rf.R"




rule optimal_rf:
    'Extracting predicted outcome (for 10 000 iterations, oob predictions) with the tuned hyperparameters'
    input:
        "results/work/rf/input_{feature}/{phenotype}/{motherorchild}/parityall_input_rf.csv",
        "results/work/rf/tuning_hyperparameters/{feature}/grid_mbo_{phenotype}_{motherorchild}_parityall.csv",
	"results/work/clean_phenotypes/{phenotype}/{motherorchild}/parityall_covar.txt"

    output:
        "results/work/rf/resources/{phenotype}/{motherorchild}/input_{feature}/mean_predicted_{phenotype}_parityall.csv",
        "results/work/rf/resources/{phenotype}/{motherorchild}/input_{feature}/PGSxparity_coeff_parityall.csv"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/rf/loop_optimal_rf.R"




rule plot_pgsxparity:
    'Plotting the potential PGSxparity interactions for each phenotype and genome'
    input:
        expand("results/work/rf/resources/{phenotype}/{motherorchild}/input_{{feature}}/PGSxparity_coeff_parityall.csv", phenotype = pheno_main, motherorchild = motherorchild)

    output:
        "results/figures/rf/GxE_input{feature}.png",
	"results/figures/rf/GxE_input{feature}_only_gd.png"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/rf/GxE_rf_plot.R"




rule plot_feature_importance:
    'Plotting feature importance for each phenotype and genome'
    input:
        "results/work/rf/input_{feature}/{phenotype}/{motherorchild}/{parity}_input_rf.csv",
        "results/work/rf/tuning_hyperparameters/{feature}/grid_mbo_{phenotype}_{motherorchild}_{parity}.csv",
        "resources/{phenotype}/{phenotype}_{motherorchild}_feature_importance_axisnames_{feature}.txt"

    output:
        "results/figures/rf/feature_importance_{phenotype}_{motherorchild}_input{feature}_{parity}.png"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/rf/feature_importance.R"




rule check_results_rf_plots:
    'Check that all files in this snakefile is created'
    input:
        expand("results/figures/rf/r2_{GE}_plot.png", GE=GE),
        expand("results/figures/rf/feature_importance_{phenotype}_{motherorchild}_input{feature}_{parity}.png", phenotype = pheno_main, motherorchild = motherorchild, feature=GE, parity=parity),
        expand("results/figures/rf/GxE_input{feature}.png", feature = GE)

    output:
        "results/checks/rf_downstream_performed.txt"

    shell:
        "touch {output[0]}"
