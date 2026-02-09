rule run_regression:
    'outcome = α+ β_1 SNP * maternal age + β_2 SNP * parity + β_3 covariate_3 + ... + ε'
    input:
        "results/work/dosage/{phenotype}/{childormother}/dosage_topSNPs.txt",
        "resources/{phenotype}/{phenotype}_{childormother}_feature_importance_axisnames_G.txt",
        "results/work/clean_phenotypes/{phenotype}/{motherorchild}/parityall_covar.txt",
        "results/work/clean_phenotypes/{phenotype}/{motherorchild}/parityall_phenofile.txt"

    output:
        "results/snpxage/{phenotype}/{childormother}_Xmaternalage.csv"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/maternalagexSNP.R"



rule check_results_maternalage:
    'Check that all files in this snakefile is created'
    input:
        expand("results/snpxage/{phenotype}/{motherorchild}_Xmaternalage.csv", phenotype = pheno_main, motherorchild = motherorchild),

    output:
        "results/checks/maternalageXsnp.txt"

    shell:
        "touch {output[0]}" 
