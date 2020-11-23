
import sys
import subprocess

input_bed = open(sys.argv[1])
input_arriba = open(sys.argv[2])
input_starfusion = open(sys.argv[3])
input_fusioncatcher = open(sys.argv[4])
input_bam = sys.argv[5]
output_fusions = open(sys.argv[6], "w")
output_coverage_file_name = sys.argv[7]

#input_bed = open("DATA/Twist_RNA_pool1_2.bed")
#input_arriba = open("Arriba_results/R20-258.fusions.tsv")
#input_starfusion = open("Results/RNA/R20-258/Fusions/star-fusion.fusion_predictions.abridged.tsv")
#input_fusioncatcher = open("fusioncatcher/R20-258/final-list_candidate-fusion-genes.hg19.txt")
#output_fusions = open("Results/RNA/R20-258/Fusions/Fusions.tsv", "w")

output_fusions.write("Caller\tgene1\tgene2\tconfidence\tFP\tbreakpoint1\tbreakpoint2\tsplit_reads1\tsplit_reads2\tsplit_reads\tsplit_reads_unique\tSpanning_pairs\tcoverage1\tcoverage2\tFFPM\tFusion_quotient1\tFusion_quotient2\tAnnotation\n")

#Only keep fusions with one gene that are in the design
design_genes = {}
for line in input_bed :
    gene = line.strip().split("\t")[3].split("_")[0]
    if gene not in design_genes :
        design_genes[gene] = ""

#Arriba fusions
header = True
for line in input_arriba :
    if header :
        header = False
        continue
    lline = line.strip().split("\t")
    gene1 = lline[0]
    gene2 = lline[1]
    #Only keep fusions with one gene that are in the design
    if (gene1 not in design_genes and gene2 not in design_genes) :
        continue
    confidence = lline[16]
    #Only keep fusions with high or medium confidence
    #if confidence == "low" :
    #    continue
    breakpoint1 = lline[4]
    breakpoint2 = lline[5]
    split_reads1 = lline[11]
    split_reads2 = lline[12]
    discordant_mates = lline[13]
    coverage1 = lline[14]
    coverage2 = lline[15]
    #Compare fusion coverage with coverage in breakpoints
    chrom1 = "chr" + breakpoint1.split(":")[0]
    pos1 = breakpoint1.split(":")[1]
    chrom2 = "chr" + breakpoint2.split(":")[0]
    pos2 = breakpoint2.split(":")[1]
    cov1 = 0
    cov2 = 0
    subprocess.call("samtools depth -d 50000 -a -r " + chrom1 + ":" + pos1 + "-" + pos1 + " " + input_bam + " > " + output_coverage_file_name, shell=True)
    output_coverage = open(output_coverage_file_name)
    for line in output_coverage :
        cov1 = int(line.strip())
    output_coverage.close()
    subprocess.call("samtools depth -d 50000 -a -r " + chrom2 + ":" + pos2 + "-" + pos2 + " " + input_bam + " > " + output_coverage_file_name, shell=True)
    output_coverage = open(output_coverage_file_name)
    for line in output_coverage :
        cov2 = int(line.strip())
    output_coverage.close()
    q1 = (cov1 / (float(split_reads1) + float(split_reads2)))
    q2 = (cov2 / (float(split_reads1) + float(split_reads2)))
    output_fusions.write("Arriba\t" + gene1 + "\t" + gene2 + "\t" + confidence + "\t\t" + breakpoint1 + "\t" + breakpoint2 + "\t" + split_reads1 + "\t" + split_reads2 + "\t\t\t" + discordant_mates + "\t" + coverage1 + "\t" + coverage2 + "\t" + str(q1) + "\t" + str(q2) + "\t" + "\n")


#Star-fusions
header = True
for line in input_starfusion :
    if header :
        header = False
        continue
    lline = line.strip().split("\t")
    gene1 = lline[0].split("--")[0]
    gene2 = lline[0].split("--")[1]
    #Only keep fusions with one gene that are in the design
    if (gene1 not in design_genes and gene2 not in design_genes) :
        continue
    Junction_read_count = lline[1]
    Spanning_Frag_count = lline[2]
    #Flag fusions with junction_read_count < 10 and Spanning_Frag_count < 2
    confidence = ""
    if int(Junction_read_count) < 10 and int(Spanning_Frag_count) < 2 :
        confidence = "Low support"
    #Remove Fusions with very weak split read support
    if int(Junction_read_count) <= 1 :
        continue
    breakpoint1 = lline[5]
    breakpoint2 = lline[7]
    FFPM = lline[9]
    DBs = lline[14]
    #Compare fusion coverage with coverage in breakpoints
    chrom1 = breakpoint1.split(":")[0]
    pos1 = breakpoint1.split(":")[1]
    chrom2 = breakpoint2.split(":")[0]
    pos2 = breakpoint2.split(":")[1]
    cov1 = 0
    cov2 = 0
    subprocess.call("samtools depth -d 50000 -a -r " + chrom1 + ":" + pos1 + "-" + pos1 + " " + input_bam + " > " + output_coverage_file_name, shell=True)
    output_coverage = open(output_coverage_file_name)
    for line in output_coverage :
        cov1 = int(line.strip())
    output_coverage.close()
    subprocess.call("samtools depth -d 50000 -a -r " + chrom2 + ":" + pos2 + "-" + pos2 + " " + input_bam + " > " + output_coverage_file_name, shell=True)
    output_coverage = open(output_coverage_file_name)
    for line in output_coverage :
        cov2 = int(line.strip())
    output_coverage.close()
    q1 = (cov1 / (float(split_reads1) + float(split_reads2)))
    q2 = (cov2 / (float(split_reads1) + float(split_reads2)))
    output_fusions.write("StarFusion\t" + gene1 + "\t" + gene2 + "\t" + confidence + "\t\t" + breakpoint1 + "\t" + breakpoint2 + "\t\t\t" + Spanning_Frag_count + "\t\t" + Junction_read_count + "\t\t\t" + FFPM + "\t" + str(q1) + "\t" + str(q2) + "\t" + DBs + "\n")



#FusionCatcher
header = True
for line in input_fusioncatcher :
    if header :
        header = False
        continue
    lline = line.strip().split("\t")
    gene1 = lline[0]
    gene2 = lline[1]
    #Only keep fusions with one gene that are in the design
    if (gene1 not in design_genes and gene2 not in design_genes) :
        continue
    fp_filters = lline[2].split(",")
    DBs = lline[2]
    common_mapping = lline[3]
    Spanning_pairs = lline[4]
    Spanning_reads_unique = lline[5]
    Fusion_finding_method = lline[7]
    breakpoint1 = lline[8]
    breakpoint2 = lline[9]
    predicted_effect = lline[15]
    #Flag fusions with Spanning_reads_unique < 5
    confidence = ""
    if int(Spanning_reads_unique) < 5 and int(Spanning_pairs) < 5:
        confidence = "Low support"
    #MAML2 have large number of FP, remove if low evidence
    if (gene1 == "MAML2" or gene2 == "MAML2") :
            if Spanning_reads_unique < 10 :
                continue
    #Flag fusions annotated that are fusions with very high probability
    fp_db = ["banned", "bodymap2", "cacg", "1000genomes", "conjoing", "cortex", "distance1000bp", "ensembl_fully_overlapping", "ensembl_same_strand_overlapping", "gtex", "hpa", "mt", "paralogs", "refseq_fully_overlapping", "refseq_same_strand_overlapping", "rrna", "similar_reads", "similar_symbols", "ucsc_fully_overlapping", "ucsc_same_strand_overlapping"]
    fp_found = ""
    for fp in fp_db :
        if fp in fp_filters :
            fp_found = "FP"
    #Compare fusion coverage with coverage in breakpoints
    chrom1 = "chr" + breakpoint1.split(":")[0]
    pos1 = breakpoint1.split(":")[1]
    chrom2 = "chr" + breakpoint2.split(":")[0]
    pos2 = breakpoint2.split(":")[1]
    cov1 = 0
    cov2 = 0
    subprocess.call("samtools depth -d 50000 -a -r " + chrom1 + ":" + pos1 + "-" + pos1 + " " + input_bam + " > " + output_coverage_file_name, shell=True)
    output_coverage = open(output_coverage_file_name)
    for line in output_coverage :
        cov1 = int(line.strip())
    output_coverage.close()
    subprocess.call("samtools depth -d 50000 -a -r " + chrom2 + ":" + pos2 + "-" + pos2 + " " + input_bam + " > " + output_coverage_file_name, shell=True)
    output_coverage = open(output_coverage_file_name)
    for line in output_coverage :
        cov2 = int(line.strip())
    output_coverage.close()
    q1 = (cov1 / (float(split_reads1) + float(split_reads2)))
    q2 = (cov2 / (float(split_reads1) + float(split_reads2)))
    output_fusions.write("FusionCatcher\t" + gene1 + "\t" + gene2 + "\t" + confidence + "\t" + fp_found + "\t" + breakpoint1 + "\t" + breakpoint2 + "\t\t\t\t" + Spanning_reads_unique + "\t" + Spanning_pairs + "\t\t\t\t" + str(q1) + "\t" + str(q2) + "\t" + DBs + "\n")