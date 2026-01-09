##### Estimate heritability and genetic correlation (parity 0 vs parity > 0) for the phenotypes #####

rule format:
    'Format header based on ldsc requriments'
    input:
        "results/work/GWAS/{pheno}/regenie/step2/QC/{motherorchild}/{sample}.txt"

    output:
        temp("results/work/GWAS/{pheno}/regenie/step2/ldsc/{motherorchild}/{sample}.txt")

    run:
        dat = pd.read_csv(input[0], header = 0, sep= ",")

        dat.to_csv(output[0], sep = " ", header = True, index = False, columns = ["CHROM", "GENPOS", "ID", "A2", "A1", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "MAF", "PVALUE"])



rule reformat_summary_statistics:
    'Reformat sumstats according to ldsc'
    input:
        "results/work/GWAS/{pheno}/regenie/step2/ldsc/{motherorchild}/{sample}.txt"

    output:
        "results/work/GWAS/{pheno}/ldsc/{motherorchild}/ldsc_{sample}.txt.sumstats.gz"

    params:
        "results/work/GWAS/{pheno}/ldsc/{motherorchild}/ldsc_{sample}.txt"

    conda:
        "../envs/environment.yml"

    shell:
        """
        set +eu
        source /home/karin.ytterberg/miniconda3/etc/profile.d/conda.sh
        conda activate ldsc
        python2 /home/karin.ytterberg/soft/ldsc/munge_sumstats.py \
        --merge-alleles /home/karin.ytterberg/soft/ldsc/w_hm3.snplist \
        --out {params[0]} \
        --sumstats {input[0]} \
        --a1-inc \
        --chunksize 500000 \
        --snp "ID" \
        --a1 "A1" \
        --a2 "A2"
        conda deactivate
        set -eu
        """



rule heritability:
   'Estimate lambda and h2 with ldsc'
    input:
        "results/work/GWAS/{pheno}/ldsc/{motherorchild}/ldsc_{sample}.txt.sumstats.gz"

    output:
        "results/work/GWAS/{pheno}/ldsc/{motherorchild}/ldsc_{sample}_h2.log"

    params:
        "results/work/GWAS/{pheno}/ldsc/{motherorchild}/ldsc_{sample}_h2"

    conda:
        "../envs/environment.yml"

    shell:
        """
        set +eu
        source /home/karin.ytterberg/miniconda3/etc/profile.d/conda.sh
        conda activate ldsc
        python2 /home/karin.ytterberg/soft/ldsc/ldsc.py \
        --h2 {input[0]} \
        --ref-ld-chr /home/karin.ytterberg/soft/ldsc/eur_w_ld_chr/ \
        --w-ld-chr /home/karin.ytterberg/soft/ldsc/eur_w_ld_chr/ \
        --out {params[0]}
        conda deactivate
        set -eu
        """


rule genetic_correlation:
    'Genetic correlation of gestational duration between parity = zero and parity > zero using ldsc'
    input:
        "results/work/GWAS/{pheno}/ldsc/{motherorchild}/ldsc_parity1.txt.sumstats.gz",
        "results/work/GWAS/{pheno}/ldsc/{motherorchild}/ldsc_parity0.txt.sumstats.gz"

    output:
        "results/work/GWAS/{pheno}/ldsc/{motherorchild}/parity1_vs_parity0_correlation.log"

    params:
        "results/work/GWAS/{pheno}/ldsc/{motherorchild}/parity1_vs_parity0_correlation"

    conda:
        "../envs/environment.yml"

    shell:
        """
        set +eu
        source /home/karin.ytterberg/miniconda3/etc/profile.d/conda.sh
        conda activate ldsc
        python2 /home/karin.ytterberg/soft/ldsc/ldsc.py --rg {input[0]},{input[1]} --ref-ld-chr /home/karin.ytterberg/soft/ldsc/eur_w_ld_chr/ --w-ld-chr /home/karin.ytterberg/soft/ldsc/eur_w_ld_chr/ --out {params[0]}
        conda deactivate
        set -eu
        """


rule prep_figure:
    'Extract heritability and genetic correlation from logs'
    input:
        "results/work/GWAS/{pheno}/ldsc/{motherorchild}/parity1_vs_parity0_correlation.log",
        "results/work/GWAS/{pheno}/ldsc/{motherorchild}/ldsc_effectSNPnoInteraction_parityall_h2.log"

    output:
        "resources/{pheno}/{motherorchild}_to_h2.txt"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/h2_prepdata.R"



rule plot_h2:
    'Plot heritability by parity and the genetic correlation'
    input:
        expand("resources/{phenotype}/{motherorchild}_to_h2.txt", phenotype = pheno_main, motherorchild = motherorchild)

    output:
        "results/figures/h2.png"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/figures/h2.R"



rule check_results_h2:
    'Check that all files in this snakefile is created'
    input:
        expand("results/work/GWAS/{pheno}/ldsc/{motherorchild}/ldsc_{sample}_h2.log", motherorchild = motherorchild, sample = h2samples, pheno = pheno_main),
        expand("results/work/GWAS/{pheno}/ldsc/{motherorchild}/parity1_vs_parity0_correlation.log", motherorchild = motherorchild, pheno = pheno_main),
        "results/figures/h2.png"

    output:
        "results/checks/herritability_and_correlation_estimated.txt"

    shell:
        "touch {output[0]}"


rule genetic_correlations_all_to_all:
	''
	input:
		"results/work/GWAS/{pheno}/ldsc/{motherorchild}/ldsc_{sample}.txt.sumstats.gz",
		expand("results/work/GWAS/{pheno}/ldsc/{motherorchild}/ldsc_{sample}.txt.sumstats.gz", pheno= pheno_main, motherorchild= motherorchild, sample= h2samples)
	output:
		'results/work/GWAS/{pheno}/ldsc/{motherorchild}/{sample}_all_to_all_correlation.log'
	params:
		"results/work/GWAS/{pheno}/ldsc/{motherorchild}/{sample}_all_to_all_correlation"

	run:
		allfiles= [infile for infile in input if wildcards.pheno not in infile]
                allfiles= ','.join(allfiles)
                outfile= params[0] + wildcards.pheno + '_rg'
                infile= input[0]
		shell("""
		set +eu
		source /home/karin.ytterberg/miniconda3/etc/profile.d/conda.sh
		conda activate ldsc
		python2 /home/karin.ytterberg/soft/ldsc/ldsc.py --rg {infile},{allfiles} --ref-ld-chr /home/karin.ytterberg/soft/ldsc/eur_w_ld_chr/ --w-ld-chr /home/karin.ytterberg/soft/ldsc/eur_w_ld_chr/ --out {params[0]}
		conda deactivate
		set -eu
		""")

rule format_RG_all_to_all:
        ''
        input:
                'results/work/GWAS/{pheno}/ldsc/{motherorchild}/{sample}_all_to_all_correlation.log'
        output:
                temp('results/work/GWAS/{pheno}/ldsc/{motherorchild}/{sample}_all_to_all_correlation_temp')
        run:
                with open(input[0], 'r') as f:
                        x= f.readlines()
                x= x[x.index('Summary of Genetic Correlation Results\n')+1:-3]
                with open(output[0], 'w') as f:
                        f.write(''.join(x))

rule concat_RG_all_to_all:
	''
	input:
		expand('results/work/GWAS/{pheno}/ldsc/{motherorchild}/{sample}_all_to_all_correlation_temp', pheno= pheno_main, motherorchild= motherorchild, sample= h2samples)
	output:
		'results/work/GWAS/ldsc_all_to_all_correlation.txt'
	run:
		df_list= list()
		for i in input:
			d= pd.read_csv(i, delim_whitespace= True, header= 0)
			df_list.append(d)
		d= pd.concat(df_list)
		d.to_csv(output[0], sep= '\t', header= True, index= False)
