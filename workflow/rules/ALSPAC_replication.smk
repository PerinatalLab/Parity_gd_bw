## replication in ALSPAC--------------------------
## @Xiaoping Wu
## 2026-01-13


# --------------------------- Dictionaries ---------------------------
GENOME=["Mother","Child"] #maternal/fetal genome
TRAITS=["gd","bw_zscore"] ## phenotype: gestational duration/birthweight
#COHORTS=["parity0","parity1","parityall","parityallGE"]
COHORTS=["parity0","parity1","parityall"]

wildcard_constraints:
    traits="(gd|bw_zscore)",   # adjust to your allowed set
    genome="(Mother|Child)",
    cohorts="(parity0|parity1|parityall)"
    
## replication for list of SNPs
import pandas as pd
toplist="data/Moba/top_loci.txt"
SNP=pd.read_csv(toplist, header=None,sep=r"\s+",names=["CHROM", "GENPOS", "ID"])
CHR=SNP["CHROM"].unique().tolist()
print(CHR)


##--------------------------- Target rules ---------------------------
rule all:
  input:
    "results/replication/gd_Mother_parity0/gwas_replication.tsv",
    "results/replication/gd_Mother_parity1/gwas_replication.tsv",
    "results/replication/gd_Mother_parityall/gwas_replication.tsv",
    "results/replication/gd_Child_parity0/gwas_replication.tsv",
    "results/replication/gd_Child_parityall/gwas_replication.tsv",
    "results/replication/bw_zscore_Mother_parity0/gwas_replication.tsv",
    "results/replication/bw_zscore_Mother_parity1/gwas_replication.tsv",
    "results/replication/bw_zscore_Mother_parityall/gwas_replication.tsv",
    "results/replication/bw_zscore_Child_parity0/gwas_replication.tsv",
    "results/replication/bw_zscore_Child_parity1/gwas_replication.tsv",
    "results/replication/bw_zscore_Child_parityall/gwas_replication.tsv"

rule run_regenie:
  input:
    expand("results/replication/{traits}_{genome}_{cohorts}/gwas.tsv",traits=TRAITS,genome=GENOME,cohorts=COHORTS),
    expand("results/replication/{traits}_{genome}_parityallGE/gwas.tsv",traits=TRAITS,genome=GENOME),
    expand("results/replication/{traits}_{genome}_{cohorts}/rs12087657_phe.regenie",traits=TRAITS,genome=GENOME,cohorts=COHORTS),
    expand("results/replication/{traits}_{genome}_parityallGE/rs12087657_phe.regenie",traits=TRAITS,genome=GENOME,cohorts=COHORTS)
    

rule run_phen:
  input:
    expand("results/replication/{traits}_{genome}_{cohorts}/ID_unrelated.eigenvec",traits=TRAITS,genome=GENOME,cohorts=COHORTS),
    expand("results/replication/{traits}_{genome}_{cohorts}/phe.tsv",traits=TRAITS,genome=GENOME,cohorts=COHORTS),
    expand("results/replication/{traits}_{genome}_{cohorts}/cov.tsv",traits=TRAITS,genome=GENOME,cohorts=COHORTS)
##--------------------------- analysis rules -------------------------
## select ALSPAC ID for each
rule ID:
  'prepare ID with both pheno and geno'
  input:
    phen="/mnt/archive/alspac/pheno/B5316/alspac-B5316-phenotype.tsv",
    linkfile="/mnt/archive/alspac/pheno/B5316/alspac-B5316-linkage.tsv",
    genoID="/mnt/archive/alspac/alspac/B5316/imputed/swapped.sample"
  output:
    expand("results/replication/{traits}_{genome}_{cohorts}/ID.tsv",traits=TRAITS,genome=GENOME,cohorts=COHORTS)
  params:
    outdir="results/replication/"
  log:
    "log/gwas_ID.log"
  conda:
    "../envs/Rbase.yml"
  shell:
    """
    Rscript --vanilla workflow/scripts/replication/ID.R  \
      --phen {input.phen} \
      --linkfile {input.linkfile} \
      --genoID {input.genoID} \
      --outdir {params.outdir} \
      >& {log}
    """


rule IDunrelated_pca_mor:
  input:
    bed="/mnt/archive/alspac/alspac/B5316/genotyping_arrays/moms/legacy2/freeze_id.bed",
    bim="/mnt/archive/alspac/alspac/B5316/genotyping_arrays/moms/legacy2/freeze_id.bim",
    fam="/mnt/archive/alspac/alspac/B5316/genotyping_arrays/moms/legacy2/freeze_id.fam",
    ID="results/replication/{traits}_Mother_{cohorts}/ID.tsv"
  output:
    "results/replication/{traits}_Mother_{cohorts}/ID_unrelated.eigenvec"
  params:
    bed="/mnt/archive/alspac/alspac/B5316/genotyping_arrays/moms/legacy2/freeze_id",
    outdir="results/replication/{traits}_Mother_{cohorts}"
  conda:
    "../envs/plink2.yml"
  threads: 5
  resources: mem_mb=10000
  shell:
    """
    jobid=$$
    mkdir -p "tmp/job{jobid}"
    awk -F "\t"  'BEGIN {{OFS="\t"}} NR>1 {{print $3,$3}}' {input.ID} >tmp/job{jobid}/id.txt
    plink2 --bfile {params.bed} --keep tmp/job{jobid}/id.txt --indep-pairwise 50000 500 0.2  --threads {threads} --memory {resources.mem_mb} --out {params.outdir}/pruning
    plink2 --bfile {params.bed} --keep tmp/job{jobid}/id.txt --exclude {params.outdir}/pruning.prune.out --king-cutoff 0.1327  --threads {threads} --memory {resources.mem_mb} --out {params.outdir}/ID_unrelated
    plink2 --bfile {params.bed} --keep {params.outdir}/ID_unrelated.king.cutoff.in.id --exclude {params.outdir}/pruning.prune.out --pca 10  --threads {threads} --memory {resources.mem_mb} --out {params.outdir}/ID_unrelated

    """

rule IDunrelated_pca_child:
  input:
    bed="/mnt/archive/alspac/alspac/B5316/genotyping_arrays/fets/freeze_id.bed",
    bim="/mnt/archive/alspac/alspac/B5316/genotyping_arrays/fets/freeze_id.bim",
    fam="/mnt/archive/alspac/alspac/B5316/genotyping_arrays/fets/freeze_id.fam",
    ID="results/replication/{traits}_Child_{cohorts}/ID.tsv"
  output:
    "results/replication/{traits}_Child_{cohorts}/ID_unrelated.eigenvec"
  params:
    bed="/mnt/archive/alspac/alspac/B5316/genotyping_arrays/fets/freeze_id",
    outdir="results/replication/{traits}_Child_{cohorts}"
  conda:
    "../envs/plink2.yml"
  threads: 5
  resources: mem_mb=10000
  shell:
    """
    jobid=$$
    mkdir -p "tmp/job{jobid}"
    awk -F "\t" 'BEGIN {{OFS="\t"}} NR>1 {{print $3,"A"}}' {input.ID} >tmp/job{jobid}/id.txt
    
    plink2 --bfile {params.bed} --keep tmp/job{jobid}/id.txt --indep-pairwise 50000 500 0.2   --threads {threads} --memory {resources.mem_mb} --out {params.outdir}/pruning
    plink2 --bfile {params.bed} --keep tmp/job{jobid}/id.txt --exclude {params.outdir}/pruning.prune.out --king-cutoff 0.1327  --threads {threads} --memory {resources.mem_mb} --out {params.outdir}/ID_unrelated
    plink2 --bfile {params.bed} --keep {params.outdir}/ID_unrelated.king.cutoff.in.id --exclude {params.outdir}/pruning.prune.out --pca 10  --threads {threads} --memory {resources.mem_mb} --out {params.outdir}/ID_unrelated
    """

rule regenie_phen:
  input:
    phen="/mnt/archive/alspac/pheno/B5316/alspac-B5316-phenotype.tsv",
    ID="results/replication/{traits}_{genome}_{cohorts}/ID.tsv",
    pca="results/replication/{traits}_{genome}_{cohorts}/ID_unrelated.eigenvec"
  output:
    phe="results/replication/{traits}_{genome}_{cohorts}/phe.tsv",
    cov="results/replication/{traits}_{genome}_{cohorts}/cov.tsv"
  params:
    traits="{traits}",
    cohorts="{cohorts}"
  log:
    "log/regnie_phen_{traits}_{genome}_{cohorts}.log"
  conda:
    "../envs/Rbase.yml"
  shell:
    """
    Rscript --vanilla workflow/scripts/replication/regenie_phen.R  \
      --phen {input.phen} \
      --pca {input.pca} \
      --ID {input.ID} \
      --traits {params.traits} \
      --cohorts {params.cohorts} \
      --phe {output.phe} \
      --cov {output.cov} \
      >& {log}
    """

rule regenie_step2_add:
  'GWAS using regenie for: whole sample, interaction(snp*parity), parity = 0 and parity > 0. Covariates are regressed out of the phenotypes and genetic markers. Linear regression is then used to test association of the residualized phenotype and the genetic marker.'
  input:
    bgen=lambda wc: f"data/ALSPAC/filtered_{int(wc.ichr):02d}.bgen" if wc.ichr.isdigit() else f"data/ALSPAC/filtered_{wc.ichr}.bgen",
    sample="data/ALSPAC/swapped.sample",
    bgi=lambda wc: f"data/ALSPAC/filtered_{int(wc.ichr):02d}.bgen.bgi" if wc.ichr.isdigit() else f"data/ALSPAC/filtered_{wc.ichr}.bgen.bgi",
    phe="results/replication/{traits}_{genome}_{cohorts}/phe.tsv",
    cov="results/replication/{traits}_{genome}_{cohorts}/cov.tsv",
    replication_list="data/Moba/top_loci.txt"
  output:
    temp("results/replication/{traits}_{genome}_{cohorts}/chr{ichr}_phe.regenie")
  params:
    trait="{traits}",
    cohorts="{cohorts}",
    ichr="{ichr}",
    outfile="results/replication/{traits}_{genome}_{cohorts}/chr{ichr}"
  threads: 5
  resources: mem_mb=5000
  conda:
    "../envs/regenie.yml"
  shell:
    """
    jobid=$$
    mkdir -p "tmp/job{jobid}"
    trait="{params.trait}"
    cohorts="{params.cohorts}"
    chr="{params.ichr}"
    awk -v n1=$chr '$1==n1 {{print $1":"$2}}' {input.replication_list} >tmp/job{jobid}/snp.txt
    
      # base covariates
      covar_cont="maternal_age_at_delivery,PC{{1:10}}"
      covar_cat="sex_assigned_at_birth"
      
      # trait-specific covariates
      if [[ "$trait" == "bw_zscore" ]]; then
          covar_cont="maternal_age_at_delivery,GA,PC{{1:10}}"
      fi
      
      # cohort-specific covariates
      if [[ "$cohorts" == "parityall" ]]; then
          covar_cat="sex_assigned_at_birth,PARITET"
      fi
      
      # run
           regenie \
              --step 2 \
              --bgen {input.bgen} \
              --sample {input.sample} \
              --bgi {input.bgi} \
              --extract tmp/job{jobid}/snp.txt \
              --covarFile {input.cov} \
              --phenoFile {input.phe} \
              --covarColList "${{covar_cont}}" \
              --catCovarList "${{covar_cat}}" \
              --bsize 1000 \
              --qt \
              --ignore-pred \
              --lowmem \
              --lowmem-prefix {params.outfile} \
              --threads {threads} \
              --no-condtl \
              --out {params.outfile}

    """




rule regenie_step2_interactive:
  'GWAS using regenie for: whole sample, interaction(snp*parity), parity = 0 and parity > 0. Covariates are regressed out of the phenotypes and genetic markers. Linear regression is then used to test association of the residualized phenotype and the genetic marker.'
  input:
    bgen=lambda wc: f"data/ALSPAC/filtered_{int(wc.ichr):02d}.bgen" if wc.ichr.isdigit() else f"data/ALSPAC/filtered_{wc.ichr}.bgen",
    sample="data/ALSPAC/swapped.sample",
    bgi=lambda wc: f"data/ALSPAC/filtered_{int(wc.ichr):02d}.bgen.bgi" if wc.ichr.isdigit() else f"data/ALSPAC/filtered_{wc.ichr}.bgen.bgi",
    phe="results/replication/{traits}_{genome}_parityall/phe.tsv",
    cov="results/replication/{traits}_{genome}_parityall/cov.tsv",
    replication_list="data/Moba/top_loci.txt"
  output:
    temp("results/replication/{traits}_{genome}_parityallGE/chr{ichr}_phe.regenie")
  params:
    trait="{traits}",
    ichr="{ichr}",
    outfile="results/replication/{traits}_{genome}_parityallGE/chr{ichr}"
  threads: 5
  resources: mem_mb=5000
  conda:
    "../envs/regenie.yml"
  shell:
    """
    jobid=$$
    mkdir -p "tmp/job{jobid}"
    trait="{params.trait}"
    chr="{params.ichr}"
    awk -v n1=$chr '$1==n1 {{print $1":"$2}}' {input.replication_list} >tmp/job{jobid}/snp.txt
    
      # base covariates
      covar_cont="maternal_age_at_delivery,PC{{1:10}}"
      covar_cat="sex_assigned_at_birth,PARITET"
      
      # trait-specific covariates
      if [[ "$trait" == "bw_zscore" ]]; then
          covar_cont="maternal_age_at_delivery,GA,PC{{1:10}}"
      fi
      
  
      
      # run
           regenie \
              --step 2 \
              --bgen {input.bgen} \
              --sample {input.sample} \
              --bgi {input.bgi} \
              --extract tmp/job{jobid}/snp.txt \
              --covarFile {input.cov} \
              --phenoFile {input.phe} \
              --covarColList "${{covar_cont}}" \
              --catCovarList "${{covar_cat}}" \
              --interaction PARITET \
              --bsize 1000 \
              --qt \
              --ignore-pred \
              --lowmem \
              --lowmem-prefix {params.outfile} \
              --threads {threads} \
              --out {params.outfile}

    """




rule concat_gwas_add:
  'Concatenate results from regenie'
  input:
    gwas=expand("results/replication/{{traits}}_{{genome}}_{{cohorts}}/chr{ichr}_phe.regenie",ichr=CHR),
    replication_list="data/Moba/top_loci.txt"
  output:
    "results/replication/{traits}_{genome}_{cohorts}/gwas.tsv"
  params:
    "results/replication/{traits}_{genome}_{cohorts}"
  log:
    "log/replication_concat_{traits}_{genome}_{cohorts}.log"
  shell:
    """
    head -1 {params[0]}/chr1_phe.regenie  > {output[0]}
    
    sigchr=$(awk '{{print $1}}' {input.replication_list}|sort|uniq)
    for cc in ${{sigchr}}
    do
          awk 'NR>1' {params[0]}/chr${{cc}}_phe.regenie  >> {output[0]}
    done

    """
rule concat_gwas_interactive:
  'Concatenate results from regenie'
  input:
    gwas=expand("results/replication/{{traits}}_{{genome}}_parityallGE/chr{ichr}_phe.regenie",ichr=CHR),
    replication_list="data/Moba/top_loci.txt"
  output:
    "results/replication/{traits}_{genome}_parityallGE/gwas.tsv"
  log:
    "log/replication_concat_{traits}_{genome}_parityallGE.log"
  params:
    "results/replication/{traits}_{genome}_parityallGE"
  shell:
    """
    head -1 {params[0]}/chr1_phe.regenie  > {output[0]}
    
    sigchr=$(awk '{{print $1}}' {input.replication_list}|sort|uniq)
    for cc in ${{sigchr}}
    do
          awk 'NR>1' {params[0]}/chr${{cc}}_phe.regenie  >> {output[0]}
    done

    """

rule merge:
  input:
    expand("results/replication/{traits}_{genome}_{cohorts}/gwas.tsv",traits=TRAITS,genome=GENOME,cohorts=COHORTS)
  output:
    "results/replication/gd_Mother_parity0/gwas_replication.tsv",
    "results/replication/gd_Mother_parity1/gwas_replication.tsv",
    "results/replication/gd_Mother_parityall/gwas_replication.tsv",
    "results/replication/gd_Child_parity0/gwas_replication.tsv",
    "results/replication/gd_Child_parityall/gwas_replication.tsv",
    "results/replication/bw_zscore_Mother_parity0/gwas_replication.tsv",
    "results/replication/bw_zscore_Mother_parity1/gwas_replication.tsv",
    "results/replication/bw_zscore_Mother_parityall/gwas_replication.tsv",
    "results/replication/bw_zscore_Child_parity0/gwas_replication.tsv",
    "results/replication/bw_zscore_Child_parity1/gwas_replication.tsv",
    "results/replication/bw_zscore_Child_parityall/gwas_replication.tsv"
  params:
    rdir="results/replication/",
    pdir="/mnt/scratch/karin/topSNPs_Xiaoping/"
  log:
    "log/replication_merge.Rlog"
  conda:
    "../envs/Rbase.yml"
  shell:
    """
    Rscript --vanilla workflow/scripts/replication/merge.R  \
      --rdir {params.rdir} \
      --pdir {params.pdir} \
      >& {log}
    """
