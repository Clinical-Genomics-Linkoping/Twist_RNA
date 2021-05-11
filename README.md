# Twist_RNA Region Östergötland implementation

---

NOTE: THIS IS A WORK IN PROGRESS SO NOT AT ALL READY IN ANY WAY YET.

---

<!-- TOC -->

- [Twist_RNA Region Östergötland implementation](#twist_rna-region-Östergötland-implementation)
- [Installation](#installation)
    - [Clone Github repo](#clone-github-repo)
    - [Be on the correct git branch](#be-on-the-correct-git-branch)
    - [Install Singularity](#install-singularity)
    - [Create Arriba Star reference](#create-arriba-star-reference)
    - [Download Star-fusion reference files](#download-star-fusion-reference-files)
    - [Download the Human genome](#download-the-human-genome)
    - [FusionCatcher](#fusioncatcher)
    - [Optional: Install Docker Engine and build Singularity containers from docker from local docker images](#optional-install-docker-engine-and-build-singularity-containers-from-docker-from-local-docker-images)
    - [Build Singularity images from gmsuppsala Dockerhub images](#build-singularity-images-from-gmsuppsala-dockerhub-images)
    - [Create snakemake conda environment](#create-snakemake-conda-environment)
    - [Create a directory for the input fastq-files](#create-a-directory-for-the-input-fastq-files)
    - [Move input fastq-files to the directory](#move-input-fastq-files-to-the-directory)
    - [Comment out demultiplexing](#comment-out-demultiplexing)
    - [Move samplesheet file and xml-file to analysis-dir](#move-samplesheet-file-and-xml-file-to-analysis-dir)
    - [Change Runfolder path for config file creation](#change-runfolder-path-for-config-file-creation)
    - [Update Config-file](#update-config-file)
- [Usage](#usage)
    - [Run `Twist_RNA.yaml` creation pipeline](#run-twist_rnayaml-creation-pipeline)
    - [Run the main pipeline (`Twist_RNA.smk`)](#run-the-main-pipeline-twist_rnasmk)

<!-- /TOC -->

# Installation

## Clone Github repo

```bash
git clone https://github.com/clinical-genomics-linkoping/Twist_RNA.git
```

## Move to `ro-implementation` branch

```bash
cd Twist_RNA
git checkout ro-implementation
```

The branch `ro-implementation` is synced regularly with `develop` branch in: https://github.com/clinical-genomics-uppsala/Twist_RNA.git

It contains changes made to get Twist_RNA running on a server at Region Östergötland in Linköping, Sweden.

## Install Singularity

Follow the installation instructions: [here](https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps 'Quick installation steps').

If you have RHEL derivative system follow [these instructions](https://sylabs.io/guides/3.0/user-guide/installation.html#install-dependencies 'Installing dependencies with yum/rpm') for installing dependencies.

These instructions were run with globally installed Singularity version 3.7.3.

## Create Arriba Star reference

<!-- TODO: Change paths in bash code blocks to relative -->

```bash
Arriba_REF="/home/Hanna/Documents/CG-Linkoping/Twist_RNA/Arriba_STAR_ref"
mkdir $Arriba_REF
cd $Arriba_REF
singularity exec \
-B $REF:/references docker://uhrigs/arriba:2.1.0 \
download_references.sh \
GRCh37+RefSeq
```

## Download Star-fusion reference files

Download from: https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.9/

```bash
SF_REF="/home/Hanna/Documents/CG-Linkoping/Twist_RNA/Star-fusion_ref"
mkdir $SF_REF
cd $SF_REF
wget -cv https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.9/GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play.tar.gz
wget -cv https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.9/GRCh37_gencode_v19_CTAT_lib_Apr032020.source.tar.gz
echo "" >> "$SF_REF"/"README.md"
md5sum "$SF_REF"/"GRCh37_gencode_v19_CTAT_lib_Apr032020.source.tar.gz" >> "$SF_REF"/"README.md"
echo "" >> "$SF_REF"/"README.md"
md5sum "$SF_REF"/"GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play.tar.gz" >> "$SF_REF"/"README.md"
```

## Download the Human genome

```bash
HG_DIR="/home/Hanna/Documents/CG-Linkoping/Twist_RNA/DATA/Test/"
HG_URL="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz"
wget -cv -O "$HG_DIR""hg19.with.mt.fa.gz" "$HG_URL"
gunzip -c "$HG_DIR""hg19.with.mt.fa.gz" > "$HG_DIR""hg19.with.mt.fa"
wget -cv -O "$HG_DIR""hg19.with.mt.md5sum.txt" http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/analysisSet/md5sum.txt
echo "" >> "$HG_DIR"/"README.md"
echo "Human reference genome downloaded from: $HG_URL" >> "$HG_DIR"/"README.md"
echo "" >> "$HG_DIR"/"README.md"
echo "md5sum:" >> "$HG_DIR"/"README.md"
echo "" >> "$HG_DIR"/"README.md"
md5sum "$HG_DIR"/"hg19.with.mt.fa.gz" >> "$HG_DIR"/"README.md"
```

## FusionCatcher

See instructions on FusionCatcher Github or follow commands in:
https://github.com/ndaniel/fusioncatcher/blob/master/data/download-human-db.sh

```bash
FUSC="/home/Hanna/Documents/CG-Linkoping/Twist_RNA/FusionCatcher"
mkdir $FUSC
cd $FUSC
wget --continue --verbose https://raw.githubusercontent.com/ndaniel/fusioncatcher/master/data/download-human-db.sh
chmod u+x download-human-db.sh
./download-human-db.sh
```

## Optional: Install Docker Engine and build Singularity containers from docker from local docker images

https://docs.docker.com/engine/install/centos/

https://stackoverflow.com/a/60316979

## Build Singularity images from gmsuppsala Dockerhub images

https://sylabs.io/guides/3.7/user-guide/build_a_container.html#overview

```bash
SIMGS_PATH="/home/Hanna/Documents/CG-Linkoping/Twist_RNA/Simgs"
mkdir $SIMGS_PATH
cd $SIMGS_PATH
sudo singularity build twist_rna_arriba.sif docker://gmsuppsala/twist_rna_arriba:develop
sudo singularity build twist_rna.sif docker://gmsuppsala/twist_rna:develop
sudo singularity build twist_rna_fusioncatcher.sif docker://olopadelab/fusioncatcher:latest
sudo singularity build twist_rna_starfusion.sif docker://gmsuppsala/twist_rna_starfusion:develop
```

## Create snakemake conda environment

Create a file named `env.yml` in the project root with this content:

```yaml
name: snakemake
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - snakemake==6.2.1
  - mamba==0.12.0
```

```bash
PROJECT_ROOT="/home/Hanna/Documents/CG-Linkoping/Twist_RNA/"
conda env create -f "$PROJECT_ROOT""env.yml"
```

## Create a directory for the input fastq-files

```bash
FASTQS="fastq_temp/RNA/"
mkdir -p $FASTQS
```

## Move input fastq-files to the directory

```bash
mv path/to/fastqs $FASTQS
```

## Comment out demultiplexing

In `src/Snakemake/workflow/Twist_RNA_workflow.smk`

```python
# This should be commented out
#include: "../rules/Fastq/demultiplex.smk"
include: "../rules/Fastq/fix_fastq_RNA.smk"
include: "../rules/Fusion/Arriba.smk"
include: "../rules/Fusion/Imbalance.smk"
#include: "../rules/Fusion/exon_splicing.smk"
include: "../rules/Fusion/Star-Fusion.smk"
include: "../rules/Fusion/FusionCatcher.smk"
include: "../rules/Fusion/Collect_fusions.smk"
include: "../rules/Fusion/Exon_skipping.smk"
include: "../rules/QC/RSeQC.smk"
include: "../rules/QC/samtools-picard-stats.smk"
include: "../rules/QC/multiqc.smk"
include: "../rules/QC/cartool.smk"
include: "../rules/QC/fastqc.smk"
include: "../rules/QC/Exon_coverage.smk"
include: "../rules/QC/Housekeeping.smk"
```

## Move samplesheet file and xml-file to analysis-dir

`RunParameters.xml`-file is obtained from the sequencer.

```bash
ANALYSIS_DIR="/home/Hanna/Documents/CG-Linkoping/Twist_RNA/"
mv /path/to/samplesheet.csv "$ANALYSIS_DIR""Test3_Pool1_Samplesheet.csv"
mv /path/to/RunParameters.xml "$ANALYSIS_DIR""RunParameters.xml"
```

## Change Runfolder path for config file creation

Around line 65 in `src/Snakemake/rules/Twist_RNA_yaml/Twist_RNA_yaml.smk` change the following:
```python
outfile.write("Runfolder: /projects/wp1/nobackup/ngs/klinik/INBOX/" + run_folder_name + "/\n\n")
```

to:

```python
outfile.write("\nRunfolder: /home/Hanna/Documents/CG-Linkoping/Twist_RNA/klinik/INBOX/" + run_folder_name + "/\n\n")
```

## Update Config-file

Update the `Config/Pipeline/configdefaults201012.yaml` file to the following:

<!-- TODO: Update config yaml when ran succesfully -->

```yaml
reference:
    # ref: "/data/ref_genomes/hg19/bwa/BWA_0.7.10_refseq/hg19.with.mt.fasta"
    ref: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/DATA/Test/hg19.fa"
    # STAR: "/projects/wp4/nobackup/workspace/jonas_test/STAR-Fusion/references/GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/"
    STAR: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/DATA/Test/STAR-Fusion/references/GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/"
    # picard_ref: "/projects/wp4/nobackup/workspace/jonas_test/STAR-Fusion/references/GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa"
    picard_ref: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/DATA/Test/STAR-Fusion/references/GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa"
    # STAR_fusion: "/projects/wp4/nobackup/workspace/jonas_test/STAR-Fusion/references/GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play/ctat_genome_lib_build_dir/"
    STAR_fusion: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/DATA/Test/STAR-Fusion/references/GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play/ctat_genome_lib_build_dir/"
    # Fusion_catcher: "/data/ref_data/fusioncatcher/human_v98/"
    Fusion_catcher: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/DATA/Test/human_v98/"
    # Arriba_index: "/projects/wp4/nobackup/workspace/jonas_test/Arriba/references2/STAR_index_GRCh37_RefSeq_hg19/"
    Arriba_index: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/DATA/Test/references2/STAR_index_GRCh37_RefSeq_hg19/"
    # Arriba_ref: "/projects/wp4/nobackup/workspace/jonas_test/Arriba/references2/GRCh37.fa"
    Arriba_ref: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/DATA/Test/references2/GRCh37.fa"
    # Arriba_gtf: "/projects/wp4/nobackup/workspace/jonas_test/Arriba/references2/RefSeq_hg19.gtf"
    Arriba_gtf: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/DATA/Test/references2/RefSeq_hg19.gtf"
    # Arriba_blacklist: "/projects/wp4/nobackup/workspace/jonas_test/Arriba/references2/blacklist_hg19_hs37d5_GRCh37_v2.1.0.tsv"
    Arriba_blacklist: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/DATA/Test/references2/blacklist_hg19_hs37d5_GRCh37_v2.1.0.tsv"
    # Arriba_refdir: "/projects/wp4/nobackup/workspace/jonas_test/Arriba/references2"
    Arriba_refdir: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/DATA/Test/references2"


configfiles:
    multiqc: "DATA/multiqc_config.yaml"

bed:
    bedfile: "DATA/Twist_RNA_pool1_2_V4.bed"
    intervals: "DATA/Twist_RNA_pool1_2_V4.interval_list"
    exonbed: "DATA/Twist_RNA_pool1_2_V4.bed"
    fpkm: "DATA/hg19_RefSeq.bed"

singularity: #Can be a docker container as well. This is converted to a singularity on the fly
    # default: "docker://gmsuppsala/twist_rna:develop"
    default: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/Simgs/twist_rna.sif"
    # default_arriba: "docker://gmsuppsala/twist_rna_arriba:develop"
    default_arriba: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/Simgs/twist_rna_arriba.sif"
    # default_starfusion: "docker://gmsuppsala/twist_rna_starfusion:develop"
    default_starfusion: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/Simgs/twist_rna_starfusion.sif"
    # Fusion_catcher: "docker://olopadelab/fusioncatcher:latest"
    Fusion_catcher: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/Simgs/twist_rna_fusioncatcher.sif"
    #default_fusioncatcher: "docker://gmsuppsala/twist_rna_fusioncatcher:develop"
    # cartool: "/projects/wp2/nobackup/Twist_Myeloid/Containers/CARTools-200206.simg"
    cartool: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/Simgs/CARTools-200206.simg"
    #execute: "singularity exec -B /data -B /projects -B /scratch "
    #multiqc: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/multiqc-1.9.simg"
    #samtools: "/projects/wp2/nobackup/Twist_Myeloid/Containers/bwa0.7.17-samtools-1.9.simg"
    #picard: "/projects/wp2/nobackup/Twist_Myeloid/Containers/bwa0.7.17-samtools-1.9.simg"
    #bwa: "/projects/wp2/nobackup/Twist_Myeloid/Containers/bwa0.7.17-samtools-1.9.simg"
    # fastqc: "/projects/wp2/nobackup/Twist_Myeloid/Containers/fastqc0.11.8.simg"
    fastqc: "/home/Hanna/Documents/CG-Linkoping/Twist_RNA/Simgs/fastqc0.11.8.simg"
    #python: "/projects/wp2/nobackup/Twist_Myeloid/Containers/python3.6.0-pysam-xlsxwriter.simg"
    #bcftools: "/projects/wp2/nobackup/Twist_Myeloid/Containers/bcftools-1.9--8.simg"
    #STAR_fusion: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/star-fusion.v1.7.0.simg"
    #Arriba: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/Arriba.simg"
    #Fusion_catcher: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/fusioncatcher_1.2.0.simg"
    #rseqc: "/projects/wp4/nobackup/workspace/somatic_dev/singularity/RSeQC_3.0.1.simg"


cartool:
    cov: "100 200 1000"
```

# Usage

## Run `Twist_RNA.yaml` creation pipeline

```bash
snakemake --printshellcmds --cores 1 -s src/Snakemake/rules/Twist_RNA_yaml/Twist_RNA_yaml.smk
```

## Run the main pipeline (`Twist_RNA.smk`)

```bash
snakemake --printshellcmds --cores 10 -s ./Twist_RNA.smk --use-singularity --cluster-config Config/Slurm/cluster.json
```
