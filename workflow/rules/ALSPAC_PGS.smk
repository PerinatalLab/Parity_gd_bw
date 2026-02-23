## Estimate polygenic score for gestational duration---------------------------
## @Xiaoping Wu
## 2026-01-15


# --------------------------- Dictionaries ---------------------------
CHR=list(range(1, 23))

PARITY=["parity0","parity1","effectSNPnoInteraction_parityall"]




##--------------------------- Target rules ---------------------------
 
rule run_PGS:
	input:
		expand('results/PGS/gd_Mother_{parity}/PGS.txt',parity=PARITY),
		'results/PGS/gd_Mother_parity0/PGS_parityall.txt',
		'results/PGS/gd_Mother_parity1/PGS_parityall.txt'
	   
rule run_PRScs:
  input:
    expand('results/PGS/gd_Mother_{parity}/beta/beta_pst_eff_a1_b0.5_phiauto_chr{ichr}.txt',
    ichr=CHR,parity=PARITY)
      
     
 
##--------------------------- analysis rules -------------------------

rule format_sumstats:
	'Format summary statistics according to the PRS-CS.'
	input:
		'/mnt/scratch/karin/Parity_gd_bw_pw/results/work/GWAS/gd/regenie/step2/QC/mother/{parity}.txt'
	output:
		'results/PGS/gd_Mother_{parity}/formated_sumstats.txt'
	log:
	  'log/format_sumstats_mother_GD_{parity}.log'
	conda:
		'../envs/Rbase.yml'
	shell:
		'''
		chmod +x ./workflow/scripts/ALSPAC_PGS/format_sumstats.R
		Rscript --vanilla ./workflow/scripts/ALSPAC_PGS/format_sumstats.R \
      --gwas {input[0]}  \
      --gwas_format {output[0]} \
      >& {log[0]}
		'''


rule PRScs:
  'Run PRS CS'
  input:
    sumstats='results/PGS/gd_Mother_{parity}/formated_sumstats.txt',
    ls_ref=expand('/mnt/scratch/xiaoping/resource/ldblk_1kg_eur/ldblk_1kg_chr{ichr}.hdf5',ichr=CHR),
    bim='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/rsid.bim',
    code='workflow/envs/PRScs/PRScs.py',
    gwas='/mnt/scratch/karin/Parity_gd_bw_pw/results/work/GWAS/gd/regenie/step2/QC/mother/{parity}.txt'
  output:
    'results/PGS/gd_Mother_{parity}/beta/beta_pst_eff_a1_b0.5_phiauto_chr{ichr}.txt'
  params:
    ld_ref='/mnt/scratch/xiaoping/resource/ldblk_1kg_eur',
    out_prefix='results/PGS/gd_Mother_{parity}/beta/beta',
    ichr='{ichr}',
    bim='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/rsid'
  log:
    'log/PRScs_chr{ichr}_mother_GD_{parity}.log'
  conda:
    '../envs/PRS-CS.yml'
  shell:
    """
    chmod +x ./workflow/scripts/ALSPAC_PGS/PRS-CS.sh
    ./workflow/scripts/ALSPAC_PGS/PRS-CS.sh {params.ld_ref}  {input.sumstats}  {input.gwas}  {params.out_prefix} {params.bim} {params.ichr} {input.code} &> {log} 
    """		

rule PGS_parity0:
  'Create PGS for all individuals'
  input:
    bed='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.bed',
    bim='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.bim',
    fam='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.fam',
    update_id="/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/update_id_chr{ichr}.tsv",
    dupid="/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/duplicated_id_chr{ichr}",
    beta='results/PGS/gd_Mother_parity0/beta/beta_pst_eff_a1_b0.5_phiauto_chr{ichr}.txt',
    phe='results/replication/gd_Mother_parity0/phe.tsv'
  output:
    'results/PGS/gd_Mother_parity0/score/chr{ichr}.sscore'
  params:
    bed='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}',
    score='results/PGS/gd_Mother_parity0/score/chr{ichr}'
  log:
    'log/score_chr{ichr}_gd_Mother_parity0.log'
  threads:1
  resources: mem_mb=4000
  conda:
    '../envs/plink2.yml'
  shell:
    '''
    jobid=$$
    mkdir -p "tmp/job{jobid}"
    awk 'NR>1 {{print $1,$2}}' {input.phe} >tmp/job{jobid}/id.txt
    plink2 --bfile {params.bed}           \
          --keep tmp/job{jobid}/id.txt    \
          --exclude {input.dupid}         \
          --update-name {input.update_id} \
          --score {input.beta} 2 4 6      \
          --threads {threads}             \
          --memory {resources.mem_mb}     \
          --out {params.score}            \
          &> {log}
    '''


rule PGS_parity1:
  'Create PGS for all individuals'
  input:
    bed='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.bed',
    bim='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.bim',
    fam='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.fam',
    update_id="/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/update_id_chr{ichr}.tsv",
    dupid="/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/duplicated_id_chr{ichr}",
    beta='results/PGS/gd_Mother_parity1/beta/beta_pst_eff_a1_b0.5_phiauto_chr{ichr}.txt',
    phe='results/replication/gd_Mother_parity1/phe.tsv'
  output:
    'results/PGS/gd_Mother_parity1/score/chr{ichr}.sscore'
  params:
    bed='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}',
    score='results/PGS/gd_Mother_parity1/score/chr{ichr}'
  log:
    'log/score_chr{ichr}_gd_Mother_parity1.log'
  threads:1
  resources: mem_mb=4000
  conda:
    '../envs/plink2.yml'
  shell:
    '''
    jobid=$$
    mkdir -p "tmp/job{jobid}"
    awk 'NR>1 {{print $1,$2}}' {input.phe} >tmp/job{jobid}/id.txt
    plink2 --bfile {params.bed}           \
          --keep tmp/job{jobid}/id.txt    \
          --exclude {input.dupid}         \
          --update-name {input.update_id} \
          --score {input.beta} 2 4 6      \
          --threads {threads}             \
          --memory {resources.mem_mb}     \
          --out {params.score}            \
          &> {log}
    '''
    
rule PGS_parityall:
  'Create PGS for all individuals'
  input:
    bed='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.bed',
    bim='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.bim',
    fam='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.fam',
    update_id="/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/update_id_chr{ichr}.tsv",
    dupid="/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/duplicated_id_chr{ichr}",
    beta='results/PGS/gd_Mother_effectSNPnoInteraction_parityall/beta/beta_pst_eff_a1_b0.5_phiauto_chr{ichr}.txt',
    phe='results/replication/gd_Mother_parityall/phe.tsv'
  output:
    'results/PGS/gd_Mother_effectSNPnoInteraction_parityall/score/chr{ichr}.sscore'
  params:
    bed='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}',
    score='results/PGS/gd_Mother_effectSNPnoInteraction_parityall/score/chr{ichr}'
  log:
    'log/score_chr{ichr}_gd_Mother_effectSNPnoInteraction_parityall.log'
  threads:1
  resources: mem_mb=4000
  conda:
    '../envs/plink2.yml'
  shell:
    '''
    jobid=$$
    mkdir -p "tmp/job{jobid}"
    awk 'NR>1 {{print $1,$2}}' {input.phe} >tmp/job{jobid}/id.txt
    plink2 --bfile {params.bed}           \
          --keep tmp/job{jobid}/id.txt    \
          --exclude {input.dupid}         \
          --update-name {input.update_id} \
          --score {input.beta} 2 4 6      \
          --threads {threads}             \
          --memory {resources.mem_mb}     \
          --out {params.score}            \
          &> {log}
    '''



rule concat_PGS:
	'Concat per-CHR PGS into a single PGS.'
	input:
		expand('results/PGS/gd_Mother_{{parity}}/score/chr{ichr}.sscore', ichr=CHR)
	output:
		'results/PGS/gd_Mother_{parity}/PGS.txt'
	params:
	  'results/PGS/gd_Mother_{parity}/score'
	log:
	  'log/concat_PGS_gd_Mother_{parity}.log'
	conda:
		'../envs/Rbase.yml'
	shell:
	  '''
	  chmod +x ./workflow/scripts/PRS_sum.R
		Rscript --vanilla ./workflow/scripts/ALSPAC_PGS/PRS_sum.R \
      --PRS_chr {params[0]}  \
      --PRS {output[0]} \
      > {log} 2>&1
	  '''

rule PGS_parity0toall:
  'beta from parity0, Create PGS for all individuals'
  input:
    bed='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.bed',
    bim='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.bim',
    fam='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.fam',
    update_id="/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/update_id_chr{ichr}.tsv",
    dupid="/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/duplicated_id_chr{ichr}",
    beta='results/PGS/gd_Mother_parity0/beta/beta_pst_eff_a1_b0.5_phiauto_chr{ichr}.txt',
    phe='results/replication/gd_Mother_parityall/phe.tsv'
  output:
    'results/PGS/gd_Mother_parity0/score_parityall/chr{ichr}.sscore'
  params:
    bed='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}',
    score='results/PGS/gd_Mother_parity0/score_parityall/chr{ichr}'
  log:
    'log/score_chr{ichr}_gd_Mother_parity0toparityall.log'
  threads:1
  resources: mem_mb=4000
  conda:
    '../envs/plink2.yml'
  shell:
    '''
    jobid=$$
    mkdir -p "tmp/job{jobid}"
    awk 'NR>1 {{print $1,$2}}' {input.phe} >tmp/job{jobid}/id.txt
    plink2 --bfile {params.bed}           \
          --keep tmp/job{jobid}/id.txt    \
          --exclude {input.dupid}         \
          --update-name {input.update_id} \
          --score {input.beta} 2 4 6      \
          --threads {threads}             \
          --memory {resources.mem_mb}     \
          --out {params.score}            \
          &> {log}
    '''


rule PGS_parity1toall:
  'beta from parity1, Create PGS for all individuals'
  input:
    bed='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.bed',
    bim='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.bim',
    fam='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}.fam',
    update_id="/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/update_id_chr{ichr}.tsv",
    dupid="/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/duplicated_id_chr{ichr}",
    beta='results/PGS/gd_Mother_parity1/beta/beta_pst_eff_a1_b0.5_phiauto_chr{ichr}.txt',
    phe='results/replication/gd_Mother_parityall/phe.tsv'
  output:
    'results/PGS/gd_Mother_parity1/score_parityall/chr{ichr}.sscore'
  params:
    bed='/mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr{ichr}',
    score='results/PGS/gd_Mother_parity1/score_parityall/chr{ichr}'
  log:
    'log/score_chr{ichr}_gd_Mother_parity1toparityall.log'
  threads:1
  resources: mem_mb=4000
  conda:
    '../envs/plink2.yml'
  shell:
    '''
    jobid=$$
    mkdir -p "tmp/job{jobid}"
    awk 'NR>1 {{print $1,$2}}' {input.phe} >tmp/job{jobid}/id.txt
    plink2 --bfile {params.bed}           \
          --keep tmp/job{jobid}/id.txt    \
          --exclude {input.dupid}         \
          --update-name {input.update_id} \
          --score {input.beta} 2 4 6      \
          --threads {threads}             \
          --memory {resources.mem_mb}     \
          --out {params.score}            \
          &> {log}
    '''

rule concat_PGS_parity0toall:
	'Concat per-CHR PGS into a single PGS.'
	input:
		expand('results/PGS/gd_Mother_parity0/score_parityall/chr{ichr}.sscore', ichr=CHR)
	output:
		'results/PGS/gd_Mother_parity0/PGS_parityall.txt'
	params:
	  'results/PGS/gd_Mother_parity0/score_parityall'
	log:
	  'log/concat_PGS_gd_Mother_parity0toparityall.log'
	conda:
		'../envs/Rbase.yml'
	shell:
	  '''
	  chmod +x ./workflow/scripts/ALSPAC_PGS/PRS_sum.R
		Rscript --vanilla ./workflow/scripts/ALSPAC_PGS/PRS_sum.R \
      --PRS_chr {params[0]}  \
      --PRS {output[0]} \
      > {log} 2>&1
	  '''
	  
rule concat_PGS_parity1toall:
	'Concat per-CHR PGS into a single PGS.'
	input:
		expand('results/PGS/gd_Mother_parity1/score_parityall/chr{ichr}.sscore', ichr=CHR)
	output:
		'results/PGS/gd_Mother_parity1/PGS_parityall.txt'
	params:
	  'results/PGS/gd_Mother_parity1/score_parityall'
	log:
	  'log/concat_PGS_gd_Mother_parity1toparityall.log'
	conda:
		'../envs/Rbase.yml'
	shell:
	  '''
	  chmod +x ./workflow/scripts/ALSPAC_PGS/PRS_sum.R
		Rscript --vanilla ./workflow/scripts/ALSPAC_PGS/PRS_sum.R \
      --PRS_chr {params[0]}  \
      --PRS {output[0]} \
      > {log} 2>&1
	  '''
