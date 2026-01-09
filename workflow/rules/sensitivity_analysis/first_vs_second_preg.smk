##### Fertility - sensitivity check. (ie would we see differences in genetics between first and later borns) #####

rule reginie_step2_fertility:
    'GWAS'
    input:
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.bed",
        multiext("results/work/clean_phenotypes/fertility/{motherorchild}/parityall_", "phenofile.txt","covar.txt", "ID.txt") #Genotyped samples that are not in covar file are removed from the analysis as well as samples with missing values at any of the covariates included

    output:
        temp("results/work/GWAS/supp/fertility/regenie/step2/temp/{motherorchild}/{chrom}_parity.regenie"),
	temp("results/work/GWAS/supp/fertility/regenie/step2/temp/{motherorchild}/{chrom}.log")

    params:
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc",
        "results/work/GWAS/supp/fertility/regenie/step2/temp/{motherorchild}/{chrom}"

    threads: 6

    run:
        shell('''
            /home/karin.ytterberg/soft/regenie/regenie_v3.2.5.gz_x86_64_Linux \
            --step 2 \
            --bed {params[0]} \
            --chr {wildcards.chrom} \
            --covarFile {input[2]} \
            --catCovarList genotyping_chip \
            --phenoFile {input[1]} \
            --keep {input[3]} \
            --bsize 300 \
            --bt \
            --firth --approx \
            --threads {threads} \
            --ignore-pred \
            --out {params[1]} \
            --verbose
            ''')



rule concat_GWAS_results_fertility:
    'Concatenate results from regenie step 2'
    input:
        expand("results/work/GWAS/supp/fertility/regenie/step2/temp/{{motherorchild}}/{chrom}_parity.regenie", chrom = chrom_list)

    output:
        temp("results/work/GWAS/supp/fertility/regenie/step2/temp/allchr/{motherorchild}/gwas.txt")

    shell:
        '''
        head -1 {input[0]} > {output[0]}
        tail -n +2 -q {input} >> {output[0]}
        '''



rule gzip_results_fertility:
   'gzip results from regine step 2'
    input:
        "results/work/GWAS/supp/fertility/regenie/step2/temp/allchr/{motherorchild}/gwas.txt"

    output:
        "results/work/GWAS/supp/fertility/regenie/step2/{motherorchild}/GWAS.txt.gz"

    shell:
        '''
        gzip -c {input[0]} > {output[0]}
        '''



rule qc_gwas_fertility:
    'Quality control of GWAS'
    input:
        "results/work/GWAS/supp/fertility/regenie/step2/{motherorchild}/GWAS.txt.gz"

    output:
        "results/work/GWAS/supp/fertility/regenie/step2/QC/{motherorchild}/gwas_qc.txt"

    run:

        dat = pd.read_csv(input[0], header = 0, sep= " ")

        dat["MAF"] = np.where(dat.A1FREQ < 0.5, dat.A1FREQ, 1-dat.A1FREQ)
        dat["BETA"] = np.where(dat.A1FREQ < 0.5, dat.BETA, -1*dat.BETA)
        dat["A2"] = np.where(dat.A1FREQ < 0.5, dat.ALLELE0, dat.ALLELE1) #refernece
        dat["A1"] = np.where(dat.A1FREQ < 0.5, dat.ALLELE1, dat.ALLELE0) #alternative
        dat["PVALUE"] = 10 ** (-dat.LOG10P)
        dat = dat.loc[dat.MAF>0.01,:]

        dat.to_csv(output[0], sep = ",", header = True, index = False, columns = ["CHROM", "GENPOS", "ID", "A2", "A1",  "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA", "MAF", "PVALUE"])        

        

rule format_GWAS_fertility:
    'Format GWAS summary statistics to find nearest protein coding gene'
    input:
        "results/work/GWAS/supp/fertility/regenie/step2/QC/{motherorchild}/gwas_qc.txt"

    output:
        temp("results/work/UCSC_genes/supp/fertility/regenie/step2/QC/{motherorchild}/temp/gwas_qc_temp.txt"),
        temp("results/work/UCSC_genes/supp/fertility/regenie/step2/QC/{motherorchild}/temp/gwas_qc_bed.txt")

    run:
        df_list= list()
        for d in pd.read_csv(input[0], sep= ',', header= 0, chunksize= 200000):
            d= d[['CHROM', 'GENPOS', 'ID','A1', 'A2', 'MAF', 'N' ,'BETA', 'SE', 'LOG10P']]
            d.columns= ['CHR', 'POS', 'RSID','REF', 'EFF', 'EAF', 'N', 'BETA', 'SE', 'LOG10P']
            d['CHR']= np.where(d.CHR== 'X', '23', d.CHR)
            d['ID']= d.CHR.apply(str) + ':' + d.POS.apply(str) + ':' + d.REF + ':' + d.EFF
            df_list.append(d)
        d= pd.concat(df_list)
        d.to_csv(output[0], sep= '\t', header= True, index= False)
        d['start']= d.POS - 1
        d.to_csv(output[1], sep= '\t', header= False, index= False, columns= ['CHR', 'start', 'POS', 'ID'])



rule bedtools_nearest_gene_fertility:
    'Bedtools to add nearest protein coding gene'
    input:
        "results/work/UCSC_genes/supp/fertility/regenie/step2/QC/{motherorchild}/temp/gwas_qc_bed.txt",
        "results/work/UCSC_genes/bedfiles/UCSC_coding_gene_coordinates_tx_hg19.txt"

    output:
        temp("results/work/UCSC_genes/temp/supp/fertility/{motherorchild}/gwas_qc.txt")

    conda:
        "../../envs/bedtools.yml"

    shell:
        'bedtools closest -t all -a {input[0]} -b {input[1]} > {output[0]}'



rule add_rsid_nearestGene_fertility:
    'Add rsid to nearest protein coding gene'
    input:
       "results/work/UCSC_genes/supp/fertility/regenie/step2/QC/{motherorchild}/temp/gwas_qc_temp.txt",
       "results/work/UCSC_genes/temp/supp/fertility/{motherorchild}/gwas_qc.txt"

    output:
        "results/work/GWAS/supp/fertility/regenie/step2/nearest_gene/{motherorchild}/gwas_qc_nearest_gene.txt"

    run:
        nearest_gene= pd.read_csv(input[1], sep= '\t', header= None, names= ['CHR', 'X', 'POS', 'ID', 'c1', 'p1', 'p2', 'nearestGene', 'Ensembl_gene'], usecols= ['ID', 'nearestGene'])
        for d in pd.read_csv(input[0], sep= '\t', header= 0, chunksize= 200000):
            d= pd.merge(d, nearest_gene, on= 'ID', how= 'left')
            d= d[['CHR', 'POS', 'ID','nearestGene']]
            d.columns= ['CHROM', 'GENPOS', 'RSID','nearestGene']
            d.to_csv(output[0], sep= '\t', header= not os.path.isfile(output[0]), index= False, mode= 'a', columns = ["CHROM","GENPOS","RSID","nearestGene"])



rule top_SNPs_fertility:
    'Finding top SNPs for the GWASs'
    input:
        "results/work/GWAS/supp/fertility/regenie/step2/QC/{motherorchild}/gwas_qc.txt",
        "results/work/GWAS/supp/fertility/regenie/step2/nearest_gene/{motherorchild}/gwas_qc_nearest_gene.txt"

    output:
        "results/work/GWAS/supp/fertility/regenie/topSNP/fertility/{motherorchild}/topSNP_gwas_qc.txt"

    conda:
        "../../envs/basicR.yml"

    script:
        "../../scripts/topSNP.R"



rule qqplot_fertility:
    'qq-plot'
    input:
        "results/work/GWAS/supp/fertility/regenie/step2/QC/{motherorchild}/gwas_qc.txt"

    output:
        "results/figures/fertility/{motherorchild}/qq.png"

    conda:
        "../../envs/basicR.yml"

    script:
        "../../scripts/figures/qq.R"



rule manhattan_fertility:
    'Manhattan plot'
    input:
        "results/work/GWAS/supp/fertility/regenie/step2/QC/{motherorchild}/gwas_qc.txt",
        "results/work/GWAS/supp/fertility/regenie/topSNP/fertility/{motherorchild}/topSNP_gwas_qc.txt"

    output:
        "results/figures/fertility/{motherorchild}/manhattan.png"

    conda:
        "../../envs/basicR.yml"

    script:
        "../../scripts/figures/manhattan.R"

     
        
rule check_fertility:
    'Check that all figures in this snakefile is created'
    input:
        expand("results/figures/fertility/{motherorchild}/qq.png", motherorchild = motherorchild),
        expand("results/figures/fertility/{motherorchild}/manhattan.png", motherorchild = motherorchild)

    output:
        "results/checks/figures_created_fertility.txt"

    shell:
        "touch {output[0]}" 
