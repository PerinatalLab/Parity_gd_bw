##### Prepare data for random forest #####

rule extractVariants:
    'Extract selected variants in REGIONFILE and add dosage'
    input:
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/phased_pol/phased-MoBaPsychGen-chr{chrom}.vcf.gz",
        "resources/rf/{phenotype}/REGIONFILE_{childormother}.txt"  # two-column format: CHROM     POS   (POS for LD-buddie r2>0.3 if variant do not exist in our data)

    output:
        temp("results/work/dosage/temp/{phenotype}/{childormother}/chrom/{chrom}_dosage_topSNPs.txt")

    conda:
        "../envs/bfc.yml"

    shell:
        '''
        bcftools +dosage {input[0]} -R {input[1]} > {output[0]}
        '''





rule format_header:
    'Format header for "{chrom}_dosage_topSNPs.txt" files'
    input:
        "results/work/dosage/temp/{phenotype}/{childormother}/chrom/{chrom}_dosage_topSNPs.txt"

    output:
        temp("results/work/dosage/temp/{phenotype}/{childormother}/{chrom}_dosage_topSNPs_correctcolnames.txt")

    shell:
        '''
        sed "s/[[][^]]*[]]//g" {input[0]} > {output[0]}
        '''





rule combinde_dosagefies:
    'Combind all "{chrom}_dosage_topSNPs_correctcolnames.txt" files'
    input:
        expand("results/work/dosage/temp/{{phenotype}}/{{childormother}}/{chr}_dosage_topSNPs_correctcolnames.txt", chr= chrom_list)

    output:
        "results/work/dosage/{phenotype}/{childormother}/dosage_topSNPs.txt"
    params:
        "results/work/dosage/temp/{{phenotype}}/{{childormother}}/"

    shell:
        '''
        awk 'NR == 1 || FNR > 1' {input} > {output[0]}
        '''





rule format_data:
    'Format data for the random forest and selecting features sets (G and GE)'
    input:
        "results/work/dosage/{phenotype}/{childormother}/dosage_topSNPs.txt",
        "results/work/clean_phenotypes/{phenotype}/{childormother}/{parity}_phenofile.txt",
        "results/work/clean_phenotypes/{phenotype}/{childormother}/{parity}_covar.txt"

    output:
        "results/work/rf/input_GE/{phenotype}/{childormother}/{parity}_input_rf.csv",
        "results/work/rf/input_G/{phenotype}/{childormother}/{parity}_input_rf.csv"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/rf_data_format.R"



rule check_results_rf_prep:
    'Check that all files in this snakefile are created'
    input: 
        expand("results/work/rf/input_GE/{phenotype}/{childormother}/{parity}_input_rf.csv", phenotype = pheno_main,childormother=motherorchild, parity = parity),
        expand("results/work/rf/input_G/{phenotype}/{childormother}/{parity}_input_rf.csv", phenotype = pheno_main, childormother= motherorchild, parity = parity) 
        
    output:
        "results/checks/rf_format_performed.txt"

    shell:
        "touch {output[0]}"
