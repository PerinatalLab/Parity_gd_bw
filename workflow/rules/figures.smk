##### Figures from genome-wide association study #####

rule histogram:
    'histogram'
    input:
        "results/work/clean_phenotypes/{pheno}/{motherorchild}/parity0_phenofile.txt",
        "results/work/clean_phenotypes/{pheno}/{motherorchild}/parity1_phenofile.txt",
        "results/work/clean_phenotypes/{pheno}/{motherorchild}/parityall_phenofile.txt"

    output:
        "results/figures/{pheno}/{motherorchild}/histogram.png"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/histogram.R" 



rule qqplot:
    'qq-plot'
    input:
        "results/work/GWAS/{pheno}/regenie/step2/QC/{motherorchild}/{sample}.txt",

    output:
        "results/figures/{pheno}/{motherorchild}/qq_{sample}.png"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/qq.R" 



rule manhattan:
    'Manhattan plot'
    input:
        "results/work/GWAS/{pheno}/regenie/step2/QC/{motherorchild}/{sample}.txt",
        "results/work/GWAS/{pheno}/regenie/topSNP/{motherorchild}/topSNP_{sample}.txt"

    output:
        "results/figures/{pheno}/{motherorchild}/manhattan_{sample}.png"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/manhattan.R"



rule manhattan_parity:
    'Manhattan plot by parity'
    input:
        "results/work/GWAS/{pheno}/regenie/step2/QC/{motherorchild}/parity0.txt",
        "results/work/GWAS/{pheno}/regenie/step2/QC/{motherorchild}/parity1.txt",
        "results/work/GWAS/{pheno}/regenie/topSNP/{motherorchild}/topSNP_parity0.txt",
        "results/work/GWAS/{pheno}/regenie/topSNP/{motherorchild}/topSNP_parity1.txt"

    output:
        "results/figures/{pheno}/manhattan_by_parity_topsnp_names_{motherorchild}.png"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/manhattan_by_parity_with_topsnp_names.R"



rule format_data_forestplot:
    'Formatting data for forestplot of previous associated variants with the phenotypes'
    input:
        "results/work/GWAS/{phenotype}/regenie/step2/QC/{motherorchild}/parity0.txt",
        "results/work/GWAS/{phenotype}/regenie/step2/QC/{motherorchild}/parity1.txt",
        "results/work/GWAS/{phenotype}/regenie/step2/QC/{motherorchild}/interaction_parityall.txt",
        "resources/{phenotype}/{phenotype}_{motherorchild}_metasnps.txt", 
        "resources/{phenotype}/LD_proxies_{phenotype}.txt" 

    output:
        "resources/{phenotype}/{motherorchild}_tomainfig.csv",
        "resources/{phenotype}/{motherorchild}_tosupfig.csv"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/data_prep_forestplot.R"



rule forestplot_mainfig:
    'Forestplot exploring SNPxParity interactions (five lowest interaction p-value for each phenotype and genome)'
    input:
        expand("resources/{phenotype}/{motherorchild}_tomainfig.csv", phenotype = pheno_main, motherorchild = motherorchild)

    output:
        "results/figures/forest_plot.png"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/forest_topsnps_mainfig.R" 



rule forestplot_supfig:
    'Forestplot exploring SNPxParity interactions (remaning SNPs that do not make it to the main fig)' 
    input:
         "resources/{phenotype}/{motherorchild}_tosupfig.csv"

    output:
        "results/figures/supp/forestplot_{phenotype}_{motherorchild}.png"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/forest_topsnps_supfig.R"



rule prep_betabeta:
    'Extract Parity 0 and Parity > 0 betas from the parity stratified GWASs for its top SNPs' 
    input:
        "results/work/GWAS/{phenotype}/regenie/step2/QC/{motherorchild}/parity0.txt",
        "results/work/GWAS/{phenotype}/regenie/step2/QC/{motherorchild}/parity1.txt", 
        "results/work/GWAS/{phenotype}/regenie/topSNP/{motherorchild}/topSNP_parity0.txt",
        "results/work/GWAS/{phenotype}/regenie/topSNP/{motherorchild}/topSNP_parity1.txt"

    output:
        "resources/{phenotype}/{motherorchild}_to_betabeta.txt"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/beta_beta_dataprep.R"



rule betabeta:
    'Plotting Parity 0 vs Parity > 0 betas (from the parity stratified GWASs, top SNPs)'
    input:
        expand("resources/{phenotype}/{motherorchild}_to_betabeta.txt", phenotype = pheno_main, motherorchild = motherorchild)

    output:
        "results/figures/betabeta_plot.png"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/beta_beta_plot.R"



rule check_results_figures:
    'Check that all figures in this snakefile is created'
    input:
        expand("results/figures/{pheno}/{motherorchild}/qq_{sample}.png", motherorchild = motherorchild, sample = allgwasoutputs,pheno = pheno_main),
        expand("results/figures/{pheno}/{motherorchild}/manhattan_{sample}.png", motherorchild = motherorchild, sample = allgwasoutputs, pheno = pheno_main),
        expand("results/figures/{pheno}/manhattan_by_parity_topsnp_names_{motherorchild}.png", pheno = pheno_main, motherorchild = motherorchild),
        "results/figures/forest_plot.png",
        expand("results/figures/supp/forestplot_{pheno}_{motherorchild}.png", pheno = pheno_main, motherorchild = motherorchild),
        "results/figures/betabeta_plot.png"

    output:
        "results/checks/figures_created.txt"

    shell:
        "touch {output[0]}"
