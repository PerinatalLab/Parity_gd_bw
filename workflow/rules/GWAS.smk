##### Genome-wide association study #####

rule reginie_step2:
    'GWAS using regenie for: whole sample, interaction(snp*parity), parity = 0 and parity > 0. Covariates are regressed out of the phenotypes and genetic markers. Linear regression is then used to test association of the residualized phenotype and the genetic marker.'
    input:
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.bed",
        multiext("results/work/clean_phenotypes/{pheno}/{motherorchild}/{parity}_", "phenofile.txt","covar.txt", "ID.txt") #Genotyped samples that are not in covar file are removed from the analysis as well as samples with missing values at any of the covariates included

    output:
        temp("results/work/GWAS/{pheno}/regenie/step2/temp/{parity}_{motherorchild}/{chrom}_{pheno}.regenie")

    params:
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc",
        "results/work/GWAS/{pheno}/regenie/step2/temp/{parity}_{motherorchild}/{chrom}"

    threads: 6

    run:
        if wildcards.parity == "parityall":
            shell('''
                /home/karin.ytterberg/soft/regenie/regenie_v3.2.5.gz_x86_64_Linux \
                --step 2 \
                --bed {params[0]} \
                --chr {wildcards.chrom} \
                --covarFile {input[2]} \
                --catCovarList genotyping_chip \
                --phenoFile {input[1]} \
                --keep {input[3]} \
                --interaction parity[=0] \
                --bsize 300 \
                --threads {threads} \
                --ignore-pred \
                --out {params[1]} \
                --verbose \
                --no-condtl
                ''')
        else:
            shell('''
                /home/karin.ytterberg/soft/regenie/regenie_v3.2.5.gz_x86_64_Linux \
                --step 2 \
                --bed {params[0]} \
                --chr {wildcards.chrom} \
                --keep {input[3]} \
                --covarFile {input[2]} \
                --catCovarList genotyping_chip \
                --phenoFile {input[1]} \
                --bsize 300 \
                --threads {threads} \
                --ignore-pred \
                --out {params[1]} \
                --verbose \
                --no-condtl
                ''')



rule concat_GWAS_results:
    'Concatenate results from regenie step 2'
    input:
        expand("results/work/GWAS/{{pheno}}/regenie/step2/temp/{{parity}}_{{motherorchild}}/{chrom}_{{pheno}}.regenie", chrom = chrom_list) 

    output:
        temp("results/work/GWAS/{pheno}/regenie/step2/temp/allchr/{motherorchild}/{parity}.txt")

    shell:
        '''
        head -1 {input[0]} > {output[0]}
        tail -n +2 -q {input} >> {output[0]}
        '''



rule gzip_results:
   'gzip results from regine step 2'
    input:
        "results/work/GWAS/{pheno}/regenie/step2/temp/allchr/{motherorchild}/parity{number}.txt"

    output:
        "results/work/GWAS/{pheno}/regenie/step2/{motherorchild}/GWAS_parity{number}.txt.gz"

    shell:
        '''
        gzip -c {input[0]} > {output[0]}
        '''



rule check_results_gwas:
    'Check that all files in this snakefile were created'
    input:
        expand("results/work/GWAS/{pheno}/regenie/step2/{motherorchild}/GWAS_{number}.txt.gz",pheno = pheno_main, number = parity, motherorchild = motherorchild)

    output:
        "results/checks/GWAS_performed.txt"

    shell:
        "touch {output[0]}"
