
rule fastqc_bam:
    input:
        bam = "STAR2/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        html="qc/{sample}/{sample}_Aligned.sortedByCoord.out_fastqc.html",
        zip="qc/{sample}/{sample}_Aligned.sortedByCoord.out_fastqc.zip",
    params:
        outdir = "qc/{sample}/",
        tmpdir = "qc/",
    log:
        "logs/qc/fastqc/{sample}_Aligned.sortedByCoord.log"
    threads: 10
    container:
        config["singularity"].get("fastqc", config["singularity"].get("default", ""))
    shell:
        "(fastqc --quiet -t {threads} --outdir {params.outdir} -d {params.tmpdir} {input}) &> {log}"
    # wrapper:
    #     "0.38.0/bio/fastqc"

# rule fastqcR2:
#     input:
#         "data_processing/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed.fastq" ##one for each R1 and one for R2 should be from a samples.yaml file
#     output:
#         html="qc/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed_fastqc.html",
#         zip="qc/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#     params: ""
#     log:
#         "logs/qc/fastqc/{sample}_{seqID}_R2_trimmed.log"
#     singularity:
#         config["singularitys"]["fastqc"]
#     wrapper:
#         "0.38.0/bio/fastqc"
