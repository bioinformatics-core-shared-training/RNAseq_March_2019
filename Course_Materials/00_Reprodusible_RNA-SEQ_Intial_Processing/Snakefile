configfile : "config.yaml"
import glob, os






WC = glob_wildcards(os.path.join('FASTQ/', "{sample}.fastq.gz"))
SAMPLES = set(WC.sample)

print( SAMPLES)



rule all:
    input:
        expand("featureCounts/{sample}.gene_count.txt", sample=SAMPLES),
        expand("Whippet/Quant/{sample}.psi.gz", sample=SAMPLES)






rule hisat2_Genome_index:
    input:
        "Genome/dm6.fa"
    output:
        "Genome/Index/dm6.1.ht2"
    threads: 5
    shell:
        "hisat2-build -p 5 {input} Genome/Index/dm6"



rule hisat2_to_Genome:
    input:
        fastq = "FASTQ/{sample}.fastq.gz",
        genome = "Genome/Index/dm6.1.ht2"
    output:
        "hisat2/{sample}.sam"
    threads: 5
    shell:
        "hisat2 -p 5 -U {input.fastq} -x  Genome/Index/dm6  > {output}"



rule samTobam:
    input:
        temp("hisat2/{sample}.sam")
    output:
        "hisat2/{sample}.sorted.bam"
    shell:
        "samtools view -b  {input}  | samtools sort - -o {output} && samtools index {output}"


rule gene_count:
    input:
        gtf = "Gene_annotation/dm6.Ensembl.genes.gtf",
        bam = "hisat2/{sample}.sorted.bam"
    output:
        "featureCounts/{sample}.gene_count.txt"
    threads: 1
    shell:
        "featureCounts -a {input.gtf} -o {output} {input.bam}"




