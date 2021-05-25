SHELL = /bin/bash

CURRENT_CONDA_ENV_NAME = snakemake
# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate $(CURRENT_CONDA_ENV_NAME)

run:
	@($(CONDA_ACTIVATE); \
	snakemake --printshellcmds --cores 1 -s src/Snakemake/rules/Twist_RNA_yaml/Twist_RNA_yaml.smk; \
	snakemake --printshellcmds --cores 10 -s ./Twist_RNA.smk --use-singularity --cluster-config Config/Slurm/cluster.json)
