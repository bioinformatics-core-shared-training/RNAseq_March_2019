include: "rules/00_download_data.skm"

#################################### Mapping and Quantification ################################
#
# In this module, we are declaring four rules that are designed to map all the reads to the  
# genome (hisat2) and count the reads that map to each gene (featureCounts). 
#
#########################################################################################    
    
rule hisat2_Genome_index:  #This is a rule and represent the first step of mapping the reads with hisat (indexing the genome)
    input:
        "Genome/dm6.fa"
    output:
        "Genome/Index/dm6.1.ht2"
    threads: 7
    conda:
        "envs/core.yaml"
    log:
        "logs/hisat2_Genome_index.log"
    shell:
        "hisat2-build -p {threads} {input} Genome/Index/dm6 2> {log}"

rule hisat2_to_Genome:
    input:
        fastq = "FASTQ/{sample}.fastq.gz",
        genome = "Genome/Index/dm6.1.ht2"
    output:
        temp("hisat2/{sample}.sam")
    threads: 3
    conda:
        "envs/core.yaml"
    shell:
        "hisat2 -p 3 -U {input.fastq} -x  Genome/Index/dm6  > {output} "

rule samTobam:
    input:
        "hisat2/{sample}.sam"
    output:
        "hisat2/{sample}.sorted.bam"
    conda:
        "envs/core.yaml"
    shell:
        "samtools view -b  {input}  | samtools sort - -o {output} && samtools index {output} "
        
rule bamstats:
    input:
        "hisat2/{sample}.sorted.bam"
    output:
        stats_txt = "QC/{sample}/{sample}.stats",
        stats_html = "QC/{sample}/{sample}.plots.html"
    params:
        "QC/{sample}/{sample}.plots"
    conda:
        "envs/core.yaml"
    shell:
        "samtools stats {input} > {output.stats_txt} && plot-bamstats -p {params} {output.stats_txt}"
    

rule featureCounts:
    input:
        gtf = "Gene_annotation/dm6.Ensembl.genes.gtf",
        bam = expand("hisat2/{sample}.sorted.bam", sample=SAMPLES)
    output:
        "featureCounts/total_samples.gene_count.txt"
    threads: 1
    conda:
        "envs/core.yaml"
    log:
        "logs/featureCounts.total.log"
    shell:
        "featureCounts -a {input.gtf} -o {output} {input.bam} 2> {log}"

        
############# Downstream analysis #############
#
# Everything below corresponds to workflows to perform different anlyses to get meaningful 
# quantitative data. On rules/ folder you can see the different snakemake modules (.skm files)
# which are `included` to be connected with the previous rules that are explicit on this
# current script. The `include` statement allows the integration of the .skm files. Notice 
# that all these snakemake scripts work under python, thus any python syntax can be used.
# 
###############################################    
    
    
include: "rules/01_stringtie.skm"    
include: "rules/02_bridge.skm"  
include: "rules/03_whippet_quant.skm"

rule get_whippet_quant:    #This is a calling point to run all whippet analysis
    input:
        expand("Whippet/Quant/{sample}.psi.gz", sample=SAMPLES)
    
include: "rules/04_whippet_delta.skm"
    

    

