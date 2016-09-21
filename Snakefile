configfile: "config.json"

rule getbundle:
    ''' Soft-link the decoy file from previusly downloaded directory
    indicated in config file'''
    output: "indexfiles/human_g1k_v37_decoy.fasta"
    input: "human_g1k_v37_decoy.fasta"
    shell: "ln -s config[gatk_bundle]/{input} {output}"

rule bowtie2:
    ''' Mapping alternative 1: using bowtie '''
    input:
        r1="trimmed_files/{dir}/{sample}.r1.allTrimmed.fq.gz",
        r2="trimmed_files/{dir}/{sample}.r2.allTrimmed.fq.gz",
        decoy="indexfiles/"human_g1k_v37_decoy.fasta"
    output:
        temp("/analysis_path/{sample}/data/{sample}.bowtie2.bam")
    params:
        bowtie="--maxins 2000 -p16",
        java="-Xmx5g",
        picard="MAX_RECORDS_IN_RAM=2500000 INPUT=/dev/stdin"
    log:
        bowtie2="logs/bowtie2.bowtie2.{sample}.log",
        picard="logs/bowtie2.sam2bam.{sample}.log",
    shell:
        "bowtie {params.bowtie} -1 {input.r1} -2 {input.r2} -x {input.decoy} 2> {log.bowtie2} |"
        "java {params.java} -jar /picard/SamFormatConverter.jar {input.picard} OUPUT={output} > {log.picard} 2>&1"


rule bwa:
    ''' Mapping alternative 1: using bwa '''
    input:
        r1="trimmed_files/{dir}/{sample}.r1.allTrimmed.fq.gz",
        r2="trimmed_files/{dir}/{sample}.r2.allTrimmed.fq.gz",
        decoy="indexfiles/"human_g1k_v37_decoy.fasta"
    output:
        temp("/analysis_path/{sample}/data/{sample}.bwa.bam")
    params:
        bwa="mem -M -t 16",
        java="-Xmx5g",
        picard="MAX_RECORDS_IN_RAM=2500000 INPUT=/dev/stdin"
    log:
        bwa="logs/bwa.bwa.{sample}.log",
        picard="logs/bwa.sam2bam.{sample}.log",
    shell:
        "bwa {params.bwa} {input.decoy} {input.r1} {input.r2}  2> {log.bwa} |"
        "java {params.java} -jar /picard/SamFormatConverter.jar {input.picard} OUPUT={output} > {log.picard} 2>&1"


rule filter:
    input:
        "/analysis_path/{sample}/data/{sample}.{mapper}.bam"
    output:
        "config[tmp]/{sample}.multimappersRemoved.{mapper}.bam
    log:
    shell:

        

