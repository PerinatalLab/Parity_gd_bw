##  credit to Pol - psnavais
from functools import reduce

rule check_results_haplotype_rf:
    'Check that all files in this snakefile are created'
    input:
        expand("results/rf/tuning_hyperparameters/G/grid_mbo_gd_mother_{sample}.csv", sample = ["h1-MT","h2-MnT","h3-PT"]),
        expand("results/rf/r2/G/r2_mbo_gd_mother_{sample}.csv",sample = ["h1-MT","h2-MnT","h3-PT"]),
	expand("results/work/rf/MT-MnT-PT/r2/G/pol_r2_mbo_G-gd_{sample}_test.csv",sample = ["h1-MT","h2-MnT","h3-PT"]),
        expand("results/work/rf/MT-MnT-PT/lm/GxE-input_G-gd-{sample}-mother.rsd",sample = ["h1-MT","h2-MnT","h3-PT"]),
        expand("results/figures/rf/MT-MnT-PT/GxE-input_G-gd-{sample}-mother.png",sample = ["h1-MT","h2-MnT","h3-PT"])

    output:
        "results/checks/haplotype_for_rf.txt"

    shell:
        "touch {output[0]}"



rule get_GT_effect_origin_rf:
    'Extract GT from VCF file for a subset of genetic variants.'
    input:
        "/mnt/scratch/karin/Parity_gd_bw_pw/resources/rf/gd/REGIONFILE_mother.txt",
        "results/effect_origin/aux/ids/gd/{fetsormoms}_toextract.txt",
        "/mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/phased_pol/phased-MoBaPsychGen-chr{CHR}.vcf.gz"

    output:
        temp("results/effect_origin/aux/GT/rf/temp/gt{CHR}-{fetsormoms}IID")

    conda:
        "../../envs/bfc.yml"

    shell:
        '''
        set -euo pipefail
        bcftools query -S {input[1]} -R {input[0]} -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {input[2]} -o {output[0]}
        '''


rule add_header_GT_effect_origin_rf:
    'Add header to genotype files.'
    input:
        "results/effect_origin/aux/ids/gd/{fetsormoms}_toextract.txt",
        "results/effect_origin/aux/GT/rf/temp/gt{CHR}-{fetsormoms}IID"

    output:
        temp("results/effect_origin/aux/GT/rf/header-gt-{CHR}-{fetsormoms}IID")

    run:
        cols= ['chr','pos','ref','eff'] + [line.strip() for line in open(input[0], 'r')]
        d= pd.read_csv(input[1], header= None, names= cols, sep= '\t')
        d.drop_duplicates(['chr', 'pos'], keep=False, inplace= True)
        d.to_csv(output[0], sep= '\t', header= True, index= False)



rule concat_GT_effect_origin_rf:
    'Collect GT from all CHR.'
    input:
        expand("results/effect_origin/aux/GT/rf/header-gt-{CHR}-{{fetsormoms}}IID", CHR= chrom_list)

    output:
        "results/effect_origin/GT/rf/allchr/GT-{fetsormoms}IID.txt"

    shell:
        ''' 
            set +o pipefail
            head -1 {input[0]} > {output[0]}
            cat {input} | grep -v 'chr' >> {output[0]}
        '''



rule get_allele_transmission_effect_origin_rf:
    'Retrieve allele transmission from family trios (after phasing).'
    input:
        expand("results/effect_origin/GT/rf/allchr/GT-{fetsormoms}IID.txt", fetsormoms = ["fets","moms"]),
        "results/effect_origin/aux/ids/gd/duos.txt"

    output:
        "results/effect_origin/haplotypes/rf/h1-MT", #maternal transmitted alleles
        "results/effect_origin/haplotypes/rf/h2-MnT", #the maternal nontransmitted alleles
        "results/effect_origin/haplotypes/rf/h3-PT", #the paternal transmitted alleles 

    script:
        "../scripts/phase_by_transmission_rf.py"



rule format_data_rf:
    'Format data for the random forest and selecting features sets (G and GE)'
    input:
        "results/effect_origin/haplotypes/rf/{sample}",
        "/mnt/work/karin/Parity_gd_bw_pw/clean_phenotypes/gd/mother/parityall_phenofile.txt",
        "/mnt/work/karin/Parity_gd_bw_pw/clean_phenotypes/gd/mother/parityall_covar.txt",

    output:
        "results/rf/MT-MnT-PT/{sample}-input_rf.csv",
        "results/rf/MT-MnT-PT/{sample}-input_rf.csv"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/rf_format_data_haplotype.R"



rule tuning_mbo:
    'Perform hyperparameter tuning using bayesian optimization'
    input:
        "results/rf/MT-MnT-PT/{sample}-input_rf.csv"

    output:
        "results/rf/MT-MnT-PT/tuning_hyperparameters/G/grid_mbo_mother_parityall_{sample}.csv",

    conda:
        "../envs/rf.yml"

    script:
        "../scripts/rf_mbo_tuning_gd.R"


rule iter_gxe_mbo:
    'Run the random forest model with the tuned hyperparameters (for 10 000 iterations) to evaluate the model performance by calculating the r2 on oob predictions + GxE'
    input:
        "results/work/rf/MT-MnT-PT/tuning_hyperparameters/G/grid_mbo_mother_parityall_{sample}.csv",
        "results/rf/MT-MnT-PT/{sample}-input_rf.csv",
        "results/work/clean_phenotypes/gd/mother/parityall_covar.txt"

    output:
        "results/work/rf/MT-MnT-PT/r2/G/pol_r2_mbo_G-gd_{sample}_test.csv",
        "results/work/rf/MT-MnT-PT/lm/GxE-input_G-gd-{sample}-mother.rsd",
        "results/figures/rf/MT-MnT-PT/GxE-input_G-gd-{sample}-mother.png"

    conda:
        "../envs/basicR.yml"

    script:
        "../scripts/rf_iter_GxE.R"                                 
