# strelka_via_snakemake_for_kronos

##dependencies 

Snakemake v5.10.0
Python 2.7
Python 3.7
Pandas 1.1.0

##Inputs:

###Sample Sheet:
SAMPLE_ID,ABSOLUTE_FILE_PATH

Alias or path for python 2.7 (python27_alias)
Alias or path for strelka2 (stelka2_path)


##Example execution:

snakemake -s scripts/strelkaV2.snk --cluster "qsub -q shahlab.q -l h_vmem=10G -pe ncpus 10 -cwd -V" --jobs 40


