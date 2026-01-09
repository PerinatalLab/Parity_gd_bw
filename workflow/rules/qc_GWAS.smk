##### Quality controll for genome-wide associatin study #####

rule qc_gwas:
    'Quality control of GWAS: maf, alt and ref allele'
    input:
        "results/work/GWAS/{pheno}/regenie/step2/{motherorchild}/GWAS_parity{number}.txt.gz"

    output:
        "results/work/GWAS/{pheno}/regenie/step2/QC/{motherorchild}/parity{number}.txt"

    run:
        dat = pd.read_csv(input[0], header = 0, sep= " ")

        dat["MAF"] = np.where(dat.A1FREQ < 0.5, dat.A1FREQ, 1-dat.A1FREQ)
        dat["BETA"] = np.where(dat.A1FREQ < 0.5, dat.BETA, -1*dat.BETA)
        dat["A2"] = np.where(dat.A1FREQ < 0.5, dat.ALLELE0, dat.ALLELE1) #refernece
        dat["A1"] = np.where(dat.A1FREQ < 0.5, dat.ALLELE1, dat.ALLELE0) #alternative
        dat["PVALUE"] = 10 ** (-dat.LOG10P)

        dat = dat.loc[dat.MAF>0.01,:]

        dat.to_csv(output[0], sep = ",", header = True, index = False, columns = ["CHROM", "GENPOS", "ID", "A2", "A1", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA", "MAF", "PVALUE"])



rule one_output_per_file:
    'Whole study population GWAS results splited to seperate files'
    input:
         "results/work/GWAS/{pheno}/regenie/step2/QC/{motherorchild}/parityall.txt"

    output:
        multiext("results/work/GWAS/{pheno}/regenie/step2/QC/{motherorchild}/", "interaction_parityall.txt", "effectparity_parityall.txt", "effectSNPwithInteraction_parityall.txt", "effectSNPnoInteraction_parityall.txt")

    run:
        dat = pd.read_csv(input[0], header = 0, sep= ",")

        datInteraction = dat.loc[dat.TEST == "ADD-INT_SNPxparity"]
        datEffectParity = dat.loc[dat.TEST == "ADD-INT_parity"]
        datEffectSNP_with_interaction = dat.loc[dat.TEST == "ADD-INT_SNP"]
        datEffectSNP_no_interaction = dat.loc[dat.TEST == "ADD"]

        datInteraction.to_csv(output[0], sep = ",", header = True, index = False)
        datEffectParity.to_csv(output[1], sep = ",", header = True, index = False)
        datEffectSNP_with_interaction.to_csv(output[2], sep = ",", header = True, index = False)
        datEffectSNP_no_interaction.to_csv(output[3], sep = ",", header = True, index = False)



rule format_GWAS:
    'Format GWAS summary statistics to find nearest protein coding gene'
    input:
        "results/work/GWAS/{pheno}/regenie/step2/QC/{motherorchild}/{qc_gwas}.txt"

    output:
        temp("results/work/UCSC_genes/GWAS/{pheno}/regenie/step2/QC/{motherorchild}/temp/{qc_gwas}_temp.txt"),
        temp("results/work/UCSC_genes/GWAS/{pheno}/regenie/step2/QC/{motherorchild}/temp/{qc_gwas}_bed.txt")

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



rule format_gene_map:
    'Format gene coordinates for bedtools intersect (for finding nearest protein coding gene).'
    input:
        "resources/2023-02-03-UCSC-genes-hg19.txt",
    output:
        "results/work/UCSC_genes/bedfiles/UCSC_coding_gene_coordinates_tx_hg19.txt"
    run:
        d= pd.read_csv(input[0], sep= '\t', header= 0)
        d= d.loc[d.cdsStart!= d.cdsEnd, :]
        d['chrom']= d.chrom.str.replace(' ', '')
        d['chrom']= d.chrom.str.replace('chr', '')
        d['chrom']= d.chrom.str.replace('X', '23')
        d['chrom']= pd.to_numeric(d.chrom, errors= 'coerce')
        d.dropna(subset= ['chrom'], inplace= True)
        d = d[d['chrom'] != 23]
        d.loc[d.txStart > d.txEnd, ['txStart', 'txEnd']]= d.loc[d.txStart > d.txEnd, ['txEnd', 'txStart']].values
        x= d.groupby(['chrom', 'geneSymbol', 'name'])['txStart'].min().reset_index()
        x1= d.groupby(['chrom', 'geneSymbol', 'name'])['txEnd'].max().reset_index()
        df= pd.merge(x, x1, on= ['chrom', 'geneSymbol', 'name'])
        df.columns= ['CHR', 'geneSymbol', 'name', 'start', 'end']
        df= df.loc[df.start != df.end, :]
        df['start']= df.start - 1
        df.sort_values(by= ['CHR', 'start'], inplace= True)
        df= df[['CHR', 'start', 'end', 'geneSymbol', 'name']]
        df[['CHR', 'start', 'end']]= df[['CHR', 'start', 'end']].apply(np.int64)
        df.to_csv(output[0], sep= '\t', header= False, index= False)



rule bedtools_nearest_gene:
    'Bedtools to add nearest protein coding gene'
    input:
        "results/work/UCSC_genes/GWAS/{pheno}/regenie/step2/QC/{motherorchild}/temp/{qc_gwas}_bed.txt",
        "results/work/UCSC_genes/bedfiles/UCSC_coding_gene_coordinates_tx_hg19.txt"
    output:
        temp("results/work/UCSC_genes/temp/{motherorchild}/{pheno}/{qc_gwas}.txt")

    conda:
        "../envs/bedtools.yml"

    shell:
        'bedtools closest -t all -a {input[0]} -b {input[1]} > {output[0]}'



rule add_rsid_nearestGene:
    'Add rsid to nearest protein coding gene'
    input:
       "results/work/UCSC_genes/GWAS/{pheno}/regenie/step2/QC/{motherorchild}/temp/{qc_gwas}_temp.txt",
       "results/work/UCSC_genes/temp/{motherorchild}/{pheno}/{qc_gwas}.txt"

    output:
        "results/work/GWAS/{pheno}/regenie/step2/nearest_gene/{motherorchild}/{qc_gwas}_nearest_gene.txt"
    wildcard_constraints:
        qc_gwas = '|'.join(["interaction_parityall", "effectparity_parityall", "effectSNPwithInteraction_parityall", "effectSNPnoInteraction_parityall", "parity0", "parity1"])

    run:
        nearest_gene= pd.read_csv(input[1], sep= '\t', header= None, names= ['CHR', 'X', 'POS', 'ID', 'c1', 'p1', 'p2', 'nearestGene', 'Ensembl_gene'], usecols= ['ID', 'nearestGene'])
        for d in pd.read_csv(input[0], sep= '\t', header= 0, chunksize= 200000):
            d= pd.merge(d, nearest_gene, on= 'ID', how= 'left')
            d= d[['CHR', 'POS', 'ID','nearestGene']]
            d.columns= ['CHROM', 'GENPOS', 'RSID','nearestGene']
            d.to_csv(output[0], sep= '\t', header= not os.path.isfile(output[0]), index= False, mode= 'a', columns = ["CHROM","GENPOS","RSID","nearestGene"])



rule top_SNPs:
    'Finding top SNPs for the GWASs'
    input:
        "results/work/GWAS/{pheno}/regenie/step2/QC/{motherorchild}/{qc_gwas}.txt",
        "results/work/GWAS/{pheno}/regenie/step2/nearest_gene/{motherorchild}/{qc_gwas}_nearest_gene.txt"

    output:
        "results/work/GWAS/{pheno}/regenie/topSNP/{motherorchild}/topSNP_{qc_gwas}.txt"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/topSNP.R"



rule check_results_gwas_qc:
    'check that all files in this snakefile is created'
    input:
        expand("results/work/GWAS/{pheno}/regenie/topSNP/{motherorchild}/topSNP_{qc_gwas}.txt", motherorchild = motherorchild, qc_gwas = allgwasoutputs, pheno = pheno_main)

    output:
        "results/checks/qc_GWAS_performed.txt"

    shell:
        "touch {output[0]}"
