##### This file contains the workflow of the cleaning of phenofiles, covariates and individuals #####

rule clean_pheo_bw:
    'Cleaning mfr and only include individuals in the genotype data passing QC - birth weight'
    input:
        "/mnt/archive/moba/pheno/v12/pheno_anthropometrics_24-01-16/mfr.gz",
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov-noMoBaIDs.txt"

    output:
        "results/work/clean_phenotypes/bw_zscore/mfr_mother_cleand.csv",
        "results/work/clean_phenotypes/bw_zscore/mfr_child_cleand.csv",
        "results/work/clean_phenotypes/bw_zscore/mfr_father_cleand.csv",

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/clean_BW.R" 



rule clean_pheo_gd:
    'Cleaning mfr and only include individuals in the genotype data passing QC - gestational duration'
    input:
        "/mnt/work/p1724/v12/PDB1724_MFR_541_v12.csv",
        "/mnt/work/p1724/v12/20221205_DODKAT_breastmilk/20221205_PDB1724_MFR_541_v12.csv",
        "/mnt/work/p1724/v12/parental_ID_to_PREG_ID.csv",
        "/mnt/work/p1724/v12/linkage_Mother_PDB1724.csv",
        "/mnt/work/p1724/v12/linkage_Child_PDB1724.csv",
        "/mnt/work/p1724/v12/linkage_Father_PDB1724.csv",
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov-noMoBaIDs.txt"

    output:
        "results/work/clean_phenotypes/gd/mfr_mother_cleand.csv",
        "results/work/clean_phenotypes/gd/mfr_child_cleand.csv",
        "results/work/clean_phenotypes/gd/mfr_father_cleand.csv" 

    conda:
        "../envs/basicR.yml" 

    script:
        "../scripts/clean_GD.R" 



rule keep_relatedness:
    'Extracting family IDs and within-family IDs for the --keep comand in plink. Will be used when removing related individuals in the sample'
    input:
        "results/work/clean_phenotypes/{pheno}/mfr_{sample}_cleand.csv" 

    output:
        temp("results/work/clean_phenotypes/{pheno}/keep_plink/{sample}_keep.txt")

    run:
        dat = pd.read_csv(input[0], header = 0, sep = ",")
        dat.to_csv(output[0], sep="\t", header = False, index = False, columns = ["FID","IID"])



rule keep_pruning:
    'Extracting family IDs and within-family IDs for the --keep comand in plink. Will be used when pruning. Keeping 10 000 mohters and 10 000 fathers'
    input:
        "results/work/clean_phenotypes/{pheno}/keep_plink/mother_keep.txt",
        "results/work/clean_phenotypes/{pheno}/keep_plink/father_keep.txt"

    output:
        temp("results/work/clean_phenotypes/{pheno}/keep_plink/keep_prune.txt")

    run:
        datm = pd.read_csv(input[0], header=None, sep = "\t")
        datf = pd.read_csv(input[1], header=None,  sep = "\t")

        datsubsamplem = datm.head(10000)
        datsubsamplef = datf.head(10000)

        dat = pd.concat([datsubsamplem,datsubsamplef])

        dat.to_csv(output[0], sep="\t", header = False, index = False, columns = [0,1])



rule pruning:
    'Generate a pruned subset of SNPs that are in approximate linkage equilibrium with each other. Use in relatedness and pca'
    input:
        multiext("/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc", ".bed", ".fam", ".bim"),
        "results/work/clean_phenotypes/{pheno}/keep_plink/keep_prune.txt"

    output:
        "results/work/snps/pruned/{pheno}/pruning_snps.prune.in"

    params:
        "results/work/snps/pruned/{pheno}/pruning_snps",
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc"

    conda:
        "../envs/plink.yml"

    shell:
        "plink2 --bfile {params[1]} --keep {input[3]} --threads 6 --indep-pairwise 100 10 0.2 --memory 20000 --out {params[0]}"



rule relatedness:
    'Only include unrelated samples'
    input:
        multiext("/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc", ".bed", ".fam", ".bim"),
        "results/work/clean_phenotypes/{pheno}/keep_plink/{sample}_keep.txt",
        "results/work/snps/pruned/{pheno}/pruning_snps.prune.in"

    output:
        multiext("results/work/clean_phenotypes/{pheno}/unrelated/unrelated_{sample}", ".king.cutoff.in.id", ".king.cutoff.out.id", ".afreq")

    params:
        "results/work/clean_phenotypes/{pheno}/unrelated/unrelated_{sample}",
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc"

    conda:
        "../envs/plink.yml"

    shell:
        "plink2 --bfile {params[1]} --keep {input[3]} --extract {input[4]} --threads 6 --memory 20000 --king-cutoff 0.1327 --freq --out {params[0]}"



rule keep_pca:
    'Extracting individuals for the --keep comand in plink. Will be used in pca'
    input:
        "results/work/clean_phenotypes/{pheno}/unrelated/unrelated_{sample}.king.cutoff.in.id"

    output:
        temp("results/work/clean_phenotypes/{pheno}/keep_plink/keep_pca_{sample}.txt")


    run:
        dat = pd.read_csv(input[0], header=None, sep = "\t")

        datsubsample = dat.head(10000)

        datsubsample.to_csv(output[0], sep="\t", header = True, index = False)



rule pca:
    'Generate PC from PCA (for unrelated individuals and in uncorrelated snps)'
    input:
        multiext("/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc", ".bed", ".fam", ".bim"),#a pruned file with all autosomal chromosomes
        "results/work/clean_phenotypes/{pheno}/keep_plink/keep_pca_{sample}.txt",
        "results/work/snps/pruned/{pheno}/pruning_snps.prune.in",
        "results/work/clean_phenotypes/{pheno}/unrelated/unrelated_{sample}.afreq"

    output:
        multiext("results/work/pca/{pheno}/pca_{sample}", ".eigenval", ".eigenvec",".acount", ".eigenvec.allele")

    params:
        "results/work/pca/{pheno}/pca_{sample}",
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc"

    conda:
        "../envs/plink.yml"

    shell:
        "plink2 --bfile {params[1]} --keep {input[3]} --extract {input[4]} --memory 20000 --pca allele-wts 10 --nonfounders --maf 0.01 --freq counts --out {params[0]}"
# allele-wts is the weight/variance for each allele of the SNP



rule projection:
    'Generate PC for whole individual in sample via projection (unrealted samples)'
    input:
        multiext("/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc", ".bed", ".fam", ".bim"),
        "results/work/pca/{pheno}/pca_{sample}.acount",
        "results/work/pca/{pheno}/pca_{sample}.eigenvec.allele",
        "results/work/clean_phenotypes/{pheno}/unrelated/unrelated_{sample}.king.cutoff.in.id"

    output:
        "results/work/pca/{pheno}/projection_{sample}.sscore"

    params:
        "results/work/pca/{pheno}/projection_{sample}",
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc"

    conda:
        "../envs/plink.yml"

    shell:
        "plink2 --bfile {params[1]} --keep {input[5]} --memory 20000 --threads 6 --read-freq {input[3]} --score {input[4]} 2 5 header-read no-mean-imputation variance-standardize --score-col-nums 6-15 --out {params[0]}"
# plink_results.eigenvec.allele project onto all samples along with an allele count `plink_results.acount` file
# --score {input[2]} 2 5 sets the `ID` (2nd column) and `A1` (5th column), `score-col-nums 6-15` sets the first 10 PCs to be projected



rule prep_covar_pheno:
    'Format files to GWAS (reginie) - phenotype, covariates, ID'
    input:
        "results/work/clean_phenotypes/{pheno}/mfr_{sample}_cleand.csv",
        "results/work/pca/{pheno}/projection_{sample}.sscore",
        "results/work/clean_phenotypes/{pheno}/unrelated/unrelated_{sample}.king.cutoff.in.id"

    output:
        multiext("results/work/clean_phenotypes/{pheno}/{sample}/", "parity0_covar.txt", "parity1_covar.txt", "parityall_covar.txt", "parity0_ID.txt", "parity1_ID.txt","parityall_ID.txt","parity0_phenofile.txt", "parity1_phenofile.txt", "parityall_phenofile.txt")

    run:
        dat = pd.read_csv(input[0], header = 0, sep = ",")
        pca = pd.read_csv(input[1], header = 0, sep = "\t", usecols=["IID","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG"])
        unrelated = pd.read_csv(input[2], header = 0, sep = "\t", usecols=["IID"])

        dat = pd.merge(dat, pca, how ="inner", on = "IID" )
        dat = pd.merge(dat, unrelated, how ="inner", on = "IID" )
        dat.fillna("NA", inplace = True)

        if str(wildcards.pheno) == "gd":
        
            dat.rename(columns = {"SVLEN_DG":"gd"}, inplace= True)

            dat["parity"] = np.where(dat.PARITET_5 == 0,0,1)

            datp0 = dat.loc[dat.PARITET_5 == 0,:]
            datp1 = dat.loc[dat.PARITET_5 >= 1,:]

            dat.to_csv(output[2], sep="\t", header = True, index = False, columns = ["FID","IID", "parity","KJONN", "genotyping_chip","maternalage","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG"])
            dat.to_csv(output[5], sep="\t", header = True, index = False, columns = ["FID","IID"])
            datp0.to_csv(output[0], sep="\t", header = True, index = False, columns = ["FID","IID", "KJONN", "genotyping_chip","maternalage","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG"])
            datp0.to_csv(output[3], sep="\t", header = True, index = False, columns = ["FID","IID"])
            datp1.to_csv(output[1], sep="\t", header = True, index = False, columns = ["FID","IID","KJONN", "genotyping_chip","maternalage","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG"])
            datp1.to_csv(output[4], sep="\t", header = True, index = False, columns = ["FID","IID"])
            datp0.to_csv(output[6], sep="\t", header = True, index = False, columns = ["FID","IID","gd"])
            datp1.to_csv(output[7], sep="\t", header = True, index = False, columns = ["FID","IID","gd"])
            dat.to_csv(output[8], sep="\t", header = True, index = False, columns = ["FID","IID","gd"])
       
        else:
        
            datp0 = dat.loc[dat.parity == 0,:]
            datp1 = dat.loc[dat.parity == 1,:]

            dat.to_csv(output[2], sep="\t", header = True, index = False, columns = ["FID","IID", "parity","sex", "genotyping_chip","pregnancy_duration","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG"])
            dat.to_csv(output[5], sep="\t", header = True, index = False, columns = ["FID","IID"])
            datp0.to_csv(output[0], sep="\t", header = True, index = False, columns = ["FID","IID", "sex", "genotyping_chip","pregnancy_duration","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG"])
            datp0.to_csv(output[3], sep="\t", header = True, index = False, columns = ["FID","IID"])
            datp1.to_csv(output[1], sep="\t", header = True, index = False, columns = ["FID","IID","sex", "genotyping_chip","pregnancy_duration","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG"])
            datp1.to_csv(output[4], sep="\t", header = True, index = False, columns = ["FID","IID"])
            datp0.to_csv(output[6], sep="\t", header = True, index = False, columns = ["FID","IID",str(wildcards.pheno)])
            datp1.to_csv(output[7], sep="\t", header = True, index = False, columns = ["FID","IID",str(wildcards.pheno)])
            dat.to_csv(output[8], sep="\t", header = True, index = False, columns = ["FID","IID",str(wildcards.pheno)])



rule fertility:
    'Fertility phenotype as a sensitivity check. (ie would we see differences in genetics between first and later borns)'
    input:
        "results/work/clean_phenotypes/gd/{sample}/parityall_covar.txt"

    output:
        multiext("results/work/clean_phenotypes/fertility/{sample}/parityall_", "covar.txt", "phenofile.txt","ID.txt")

    run:
        dat_covar = pd.read_csv(input[0], header = 0, sep = "\t")

        dat_pheno = dat_covar[["FID","IID","parity"]]
        ID = dat_covar[["FID","IID"]]

        dat_covar.to_csv(output[0], sep="\t", header = True, index = False, columns = ["FID","IID", "genotyping_chip","maternalage","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG"])
        dat_pheno.to_csv(output[1], sep="\t", header = True, index = False, columns = ["FID","IID","parity"])
        ID.to_csv(output[2], sep="\t", header = True, index = False, columns = ["FID","IID"])



rule table:
    'Create descriptive table'
    input:
        "results/work/clean_phenotypes/{pheno}/{sample}/parityall_covar.txt",
        "results/work/clean_phenotypes/{pheno}/{sample}/parityall_phenofile.txt"

    output:
        "results/work/clean_phenotypes/{pheno}/{sample}/table.csv"

    conda:
       "../envs/basicR.yml"

    script:
        "../scripts/table.R"



rule check_results_phenotype:
    'Check that all files in this snakefile are created'
    input:
        expand("results/work/clean_phenotypes/{pheno}/{sample}/table.csv", pheno = pheno_all, sample = samples)

    output:
        "results/checks/phenotype_performed.txt"

    shell:
        "touch {output[0]}"
