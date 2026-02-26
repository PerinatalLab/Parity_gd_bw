## creadit to Pol - psnavais
from functools import reduce

rule list_vcf_ids:
    "Obtain list of IID present in each chromosome."
    input:
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/phased_pol/phased-MoBaPsychGen-chr{CHR}.vcf.gz" 

    output:
        temp("results/effect_origin/aux/vcf_ids/temp/{CHR}-ids.txt")

    conda:
        "../envs/bfc.yml"

    shell:
        "bcftools query -l {input[0]} > {output[0]}"



rule merge_vcf_ids:
    "Keep only IIDs present in all chromosomes."
    input:
        expand("results/effect_origin/aux/vcf_ids/temp/{CHR}-ids.txt", CHR= chrom_list)

    output:
        "results/effect_origin/aux/vcf_ids/allchr-ids.txt"

    run:
        df_list= list()
        for infile in input:
            d= pd.read_csv(infile, header= None, names= ['IID'])
            df_list.append(d)
        d= reduce(lambda x, y: pd.merge(x, y, on = 'IID', how = 'inner'), df_list) # Joins two files in df_list at the time and return a sigle file with all ids present in all chrom (inner join)
        d.to_csv(output[0], sep= '\t', header= True, index= False)



rule create_duo_list:
    "Only keep duos (mother child) from .fam file"
    input:
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.fam"

    output:
        "results/effect_origin/aux/duo-ids.txt"

    shell:
        r"""
        awk 'BEGIN {{OFS="\t"; print "FID","Child","Mother"}} $4 != 0 {{print $1,$2,$4}}' {input[0]} > {output[0]}
        """



rule list_duo_ids:
    "Make a list of trio IDs with genotype data."
    input:
        "results/work/clean_phenotypes/{pheno}/mother/parityall_phenofile.txt",
        "results/work/clean_phenotypes/{pheno}/child/parityall_phenofile.txt",
        "results/work/pca/{pheno}/projection_mother.sscore",
        "results/work/pca/{pheno}/projection_child.sscore",
        "results/effect_origin/aux/vcf_ids/allchr-ids.txt", 
        "results/effect_origin/aux/duo-ids.txt"

    output:
        "results/effect_origin/aux/ids/{pheno}/fets_toextract.txt",
        "results/effect_origin/aux/ids/{pheno}/moms_toextract.txt",
        "results/effect_origin/aux/ids/{pheno}/duos.txt"

    run:
        dm= pd.read_csv(input[0], sep= '\t', header= 0) #phenofile mother
        dc= pd.read_csv(input[1], sep= '\t', header= 0) #phenofile children

        pcs_m = pd.read_csv(input[2], header = 0, sep = "\t", usecols=["IID","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG"])
        pcs_c = pd.read_csv(input[3], header = 0, sep = "\t", usecols=["IID","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG"])
        pcs = pd.concat([pcs_m, pcs_c], axis=0, ignore_index=True)
        print(pcs.head(10))
        x= pd.read_csv(input[4], sep= '\t', header= 0)
        
        fets= format_duos(dc, x) #function in commmon
        print(fets.shape) #dim
        moms= format_duos(dm, x) #function in common
        print(fets.shape) #dim 

        fets.to_csv(output[0], columns= ['IID'], sep= '\t', header= False, index= False)
        moms.to_csv(output[1], columns= ['IID'], sep= '\t', header= False, index= False)

        df= pd.read_csv(input[5], sep= "\t", header= 0)
        df= df.loc[df.iloc[:, 1].isin(pcs["IID"]) & df.iloc[:, 2].isin(pcs["IID"]),:]
        df.to_csv(output[2], sep= '\t', header= True, index= False)



rule format_sumstats:
    "Remove non-necessary rows from summary statistics"
    input:
        "results/work/GWAS/{pheno}/regenie/topSNP/{motherorchild}/topSNP_{parity}.txt"

    output:
        temp("results/effect_origin/aux/top_signals/{pheno}/{motherorchild}-regions_to_extract-{parity}.txt")
    run:
        df_list= list()
        for infile in input:
            d= pd.read_csv(infile, sep = ',', header= 0)
            d = d[d['cut_off'] == 5e-8]
            d['GENPOS']= d.GENPOS.apply(int).apply(str)
            d['pos2']= d.GENPOS
            d['CHROM']= d.CHROM.apply(str)
            df_list.append(d)
            d= pd.concat(df_list)
            d.sort_values(['CHROM', 'GENPOS'], inplace= True, ascending= True)
            d.to_csv(output[0], header= False, index= False, sep= '\t', columns= ['CHROM', 'GENPOS', 'pos2'])



rule get_GT_effect_origin:
    'Extract GT from VCF file for a subset of genetic variants.'
    input:
        "results/effect_origin/aux/top_signals/{pheno}/{motherorchild}-regions_to_extract-{parity}.txt",
        "results/effect_origin/aux/ids/{pheno}/{fetsormoms}_toextract.txt",
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/phased_pol/phased-MoBaPsychGen-chr{CHR}.vcf.gz"

    output:
        temp("results/effect_origin/aux/GT/{pheno}/temp/{motherorchild}-{parity}_gt{CHR}-{fetsormoms}IID")

    conda:
        "../envs/bfc.yml"

    shell:
        '''
        set -euo pipefail
        if [[ "{wildcards.pheno}" == "gd" && \
        "{wildcards.motherorchild}" == "child" && \
        "{wildcards.parity}" == "parity1" ]]; then
        touch {output[0]}

        else 

        bcftools query -S {input[1]} -R {input[0]} -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {input[2]} -o {output[0]}
        fi
        '''



rule add_header_GT_effect_origin:
    'Add header to genotype files.'
    input:
        "results/effect_origin/aux/ids/{pheno}/{fetsormoms}_toextract.txt",
        "results/effect_origin/aux/GT/{pheno}/temp/{motherorchild}-{parity}_gt{CHR}-{fetsormoms}IID"

    output:
        temp("results/effect_origin/aux/GT/{pheno}/gt-header-{motherorchild}-{parity}-{CHR}-{fetsormoms}IID")

    run:
        cols= ['chr','pos','ref','eff'] + [line.strip() for line in open(input[0], 'r')]
        d= pd.read_csv(input[1], header= None, names= cols, sep= '\t')
        d.drop_duplicates(['chr', 'pos'], keep=False, inplace= True)
        d.to_csv(output[0], sep= '\t', header= True, index= False)



rule concat_GT_effect_origin:
    'Collect GT from all CHR.'
    input:
        expand("results/effect_origin/aux/GT/{{pheno}}/gt-header-{{motherorchild}}-{{parity}}-{CHR}-{{fetsormoms}}IID", CHR= chrom_list)

    output:
        "results/effect_origin/GT/{pheno}/allchr/{motherorchild}-{parity}_GT-{fetsormoms}IID.txt"

    shell:
        '''
            set +o pipefail
            if [[ "{wildcards.pheno}" == "gd" && \
            "{wildcards.motherorchild}" == "child" && \
            "{wildcards.parity}" == "parity1" ]]; then
            touch {output[0]}

            else
            head -1 {input[0]} > {output[0]}
            cat {input} | grep -v 'chr' >> {output[0]}
            fi
        '''



rule get_allele_transmission_effect_origin:
    'Retrieve allele transmission from family trios (after phasing).'
    input:
        expand("results/effect_origin/GT/{{pheno}}/allchr/{{motherorchild}}-{{parity}}_GT-{fetsormoms}IID.txt", fetsormoms = ["fets","moms"]),
        "results/effect_origin/aux/ids/{pheno}/duos.txt"
        
    output:
        "results/effect_origin/haplotypes/{pheno}/{motherorchild}-{parity}-h1-MT", #maternal transmitted alleles
        "results/effect_origin/haplotypes/{pheno}/{motherorchild}-{parity}-h2-MnT", #the maternal nontransmitted alleles
        "results/effect_origin/haplotypes/{pheno}/{motherorchild}-{parity}-h3-PT", #the paternal transmitted alleles 

    script:
        "../scripts/phase_by_transmission.py"



rule linear_hypotheses:
	'Running regression: outcome ~ h1 + h2 + h3 + ... '
    input:
        "results/effect_origin/haplotypes/{pheno}/{motherorchild}-{parity}-h1-MT", #maternal transmitted alleles
        "results/effect_origin/haplotypes/{pheno}/{motherorchild}-{parity}-h2-MnT", #the maternal nontransmitted alleles
        "results/effect_origin/haplotypes/{pheno}/{motherorchild}-{parity}-h3-PT", #the paternal transmitted alleles 	
        "results/work/clean_phenotypes/{pheno}/{motherorchild}/parityall_phenofile.txt",
        "results/work/clean_phenotypes/{pheno}/mother/parityall_covar.txt",
        "results/work/clean_phenotypes/{pheno}/child/parityall_covar.txt",
        "results/effect_origin/aux/ids/{pheno}/duos.txt"


    output:
        "results/effect_origin/{pheno}/lh/{motherorchild}-{parity}-results.txt",
        "results/effect_origin/{pheno}/lh/plot-{motherorchild}-{parity}-results.csv"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/linear_hypotheses.R"



rule plot_linear_hypotheses:
    'Plotting the estimates from "linear_hypotheses"'
    input:
        "results/effect_origin/{pheno}/lh/plot-{motherorchild}-{parity}-results.csv",
        "results/work/GWAS/{pheno}/regenie/topSNP/{motherorchild}/topSNP_{parity}.txt"

    output:
        "results/effect_origin/{pheno}/lh/figure-{motherorchild}-{parity}.png"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/plot_haplotype.R"



rule check_results_haplotype:
    'Check that all files in this snakefile are created'
    input:
        expand("results/effect_origin/{pheno}/lh/figure-{motherorchild}-{parity}.png",pheno = pheno_main, motherorchild = motherorchild, parity=["parity1","parity0","effectSNPnoInteraction_parityall"])
    output:
        "results/checks/haplotype.txt"

    shell:
        "touch {output[0]}"
