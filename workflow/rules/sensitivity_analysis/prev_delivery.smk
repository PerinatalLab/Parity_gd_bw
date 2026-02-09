##### Sensitivity analysis, prev delivery #####


rule extractVariants_prevdelivery:
    'Extract selected variants in REGIONFILE and add dosage'
    input:
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/phased_pol/phased-MoBaPsychGen-chr{chrom}.vcf.gz",
        "resources/{prev}/{parity}/REGIONFILE_{childormother}.txt"  # two-column format: CHROM     POS   (POS for LD-buddie r2>0.3 if variant do not exist in our data)

    output:
        temp("results/work/dosage/temp/{prev}/{childormother}/{parity}/chrom/{chrom}_dosage.txt")

    conda:
        "../envs/bfc.yml"

    shell:
        '''
        bcftools +dosage {input[0]} -R {input[1]} > {output[0]}
        '''





rule format_header_prevdelivery:
    'Format header for "{chrom}_dosage_topSNPs.txt" files'
    input:
        "results/work/dosage/temp/{prev}/{childormother}/{parity}/chrom/{chrom}_dosage.txt"

    output:
        temp("results/work/dosage/temp/{prev}/{childormother}/{parity}/{chrom}_dosage_correctcolnames.txt")

    shell:
        '''
        sed "s/[[][^]]*[]]//g" {input[0]} > {output[0]}
        '''





rule combinde_dosagefies_prevdelivery:
    'Combind all "{chrom}_dosage_topSNPs_correctcolnames.txt" files'
    input:
        expand("results/work/dosage/temp/{{prev}}/{{childormother}}/{{parity}}/{chr}_dosage_correctcolnames.txt", chr= chrom_list)

    output:
        "results/work/dosage/{prev}/{childormother}/{parity}/dosage_all.txt"
    params:
        "results/work/dosage/temp/{prev}/{childormother}/"

    shell:
        '''
        awk 'NR == 1 || FNR > 1' {input} > {output[0]}
        '''



rule prep_mfr_gd_child:
    'Running regression for mothers with two preg in MoBa. On gd, fetal genome'
    input:
        "/mnt/work/p1724/v12/PDB1724_MFR_541_v12.csv",#mfr
        "/mnt/work/p1724/v12/parental_ID_to_PREG_ID.csv",#parentalids
        "/mnt/work/p1724/v12/20221205_DODKAT_breastmilk/20221205_PDB1724_MFR_541_v12.csv",#time of death
        "/mnt/work/p1724/v12/linkage_Mother_PDB1724.csv",#sentrix m
        "/mnt/work/p1724/v12/linkage_Child_PDB1724.csv",#sentrix c
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov-noMoBaIDs.txt",# qc geno
        expand("results/work/dosage/{prev}/{childormother}/{parity}/dosage_all.txt", childormother = ["child"], parity = parity, prev = ["prev_ptb"]),
        expand("resources/prev_ptb/{parity}/rsid_child.txt", parity = parity),        
        "results/work/clean_phenotypes/gd/child/parityall_covar.txt"

    output:
        "resources/forfigs/gd/child_prev.csv"


    conda: 
        "../envs/basicR.yml"

    script:
        "../scripts/prev_gd.R"




rule prep_mfr_gd_mother:
    'Running regression for mothers with two preg in MoBa. On gd, maternal genome'
    input:
        "/mnt/work/p1724/v12/PDB1724_MFR_541_v12.csv",#mfr
        "/mnt/work/p1724/v12/parental_ID_to_PREG_ID.csv",#parentalids
        "/mnt/work/p1724/v12/20221205_DODKAT_breastmilk/20221205_PDB1724_MFR_541_v12.csv",#time of death
        "/mnt/work/p1724/v12/linkage_Mother_PDB1724.csv",#sentrix m
        "/mnt/work/p1724/v12/linkage_Child_PDB1724.csv",#sentrix c
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov-noMoBaIDs.txt",# qc geno
        expand("results/work/dosage/{prev}/{childormother}/{parity}/dosage_all.txt", childormother = ["mother"], parity = parity, prev = ["prev_ptb"]),
        expand("resources/prev_ptb/{parity}/rsid_mother.txt", parity = parity),
        "results/work/clean_phenotypes/gd/mother/parityall_covar.txt"

    output:
        "resources/forfigs/gd/mother_prev.csv",
        "resources/forfigs/gd/mother_firstvssecond.csv"


    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/prev_gd.R"


rule prep_mfr_bw:
    'Running regression for mothers with two preg in MoBa. On bw.'
    input:
        "/mnt/work/p1724/v12/PDB1724_MFR_541_v12.csv",#mfr
        "/mnt/work/p1724/v12/parental_ID_to_PREG_ID.csv",#parentalids
        "/mnt/work/p1724/v12/20221205_DODKAT_breastmilk/20221205_PDB1724_MFR_541_v12.csv",#time of death
        "/mnt/work/p1724/v12/linkage_Mother_PDB1724.csv",#sentrix m
        "/mnt/work/p1724/v12/linkage_Child_PDB1724.csv",#sentrix c
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov-noMoBaIDs.txt",# qc geno
        expand("results/work/dosage/{prev}/{{childormother}}/{parity}/dosage_all.txt", parity = parity, prev = ["prev_bw"]),        
        expand("resources/prev_bw/{parity}/rsid_{{childormother}}.txt", parity = parity),
        "results/work/clean_phenotypes/bw_zscore/{childormother}/parityall_covar.txt"

    output:
        "resources/forfigs/bw_zscore/{childormother}_prev.csv"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/prev_bw.R"




rule plot_previous_delivery:
    'Plot betas with and without information of previous delivery'
    input:
        expand("resources/forfigs/{phenotype}/{childormother}_prev.csv", phenotype = ,pheno_main, childormother = motherorchild)

    output:
        "results/figures/prev_delivery.png"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/prev_plot.R"




rule plot_firstvssecond:
    'Plot betas for first and second preg (same mothers in both pregnancies). On gd.'
    input:
        "resources/forfigs/gd/mother_firstvssecond.csv"

    output:
        "results/figures/firstvssecond.png"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/firstvssecondpreg.R"



rule check_results_prev_ptb:
    'Check that all files in this snakefile are created'
    input:
        "results/figures/firstvssecond.png",
        "results/figures/prev_delivery.png",

    output:
        "results/checks/prev_ptb_performed.txt"

    shell:
        "touch {output[0]}"
