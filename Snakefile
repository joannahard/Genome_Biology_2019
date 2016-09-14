configfile: "config.json"

rule all:
    input:
       "analysis_path/{sample}/{sample}.r1.allTrimmed.fastqc", #don't know the exact names
       "analysis_path/{sample}/{sample}.r2.allTrimmed.fastqc",
       "analysis_path/{sample}/{sample}.singletts.fastqc"

rule wgaAdaptorTrimmer:
    input:
        "{sample}.{read}.fastq.gz"
    output:
        temp("{sample}.{read}.wgaTrimmed2.fq")
    log:
        "analysis_path/log/RubiconWgaTrimming.{sample}.{read}.log.txt"
    shell:
        "python2.7 scripts/wgaAdapterTrimmer.py -i {input} > {output} 2> {log};"

rule cutadapt:
    input:
        "{sample}.{read}.wgaTrimmed2.fq"
    output:
        temp("{sample}.{read}.wgaAndilluminaTrimmed.fq")
    params:
        " ".join(expand("-a {seq}", seq=config["adapterseqs"]))
    log:
        "analysis_path/logs/illuminaAndNexteraTrimming.{sample}.{read}.log.txt"
    shell:
        "cutadapt -n 3 {params} {input} > {output} 2> {log};"

rule TrimBWAstyle:
    input:
       "{sample}.{read}.wgaAndilluminaTrimmed.fq"
    output:
        temp("{sample}.{read}.wgaIlluminaAndQualityTrimmed.fq")
    log:
        "analysis_path/logs/qualityTrimming.{sample}.{read}.log.txt"
    shell:
        "perl scripts/TrimBWAstyle.pl {input} > {output} 2> {log};"

rule removeEmptyReads:
    input:
        "{sample}.r1.wgaIlluminaAndQualityTrimmed.fq",
        "{sample}.r2.wgaIlluminaAndQualityTrimmed.fq"
    output:
        "analysis_path/{sample}.r1.allTrimmed.fq.gz", 
        "analysis_path/{sample}.r2.allTrimmed.fq.gz", 
        "analysis_path/{sample}.singletts.fq.gz"
    params:
        "analysis_path/{sample}.r1.allTrimmed.fq",
        "analysis_path/{sample}.r2.allTrimmed.fq",
        "analysis_path/{sample}.singletts.fq",
    log:
        "analysis_path/logs/removeEmptyReads.{sample}.log.txt"
    shell:
        "python scripts/removeEmptyReads.py {input} {params} 2> {log};"
        "gzip {params};"

rule fastqc:
    input:
        "analysis_path/{sample}.{data}.fq.gz"
    output:
        "analysis_path/{sample}/{sample}.{data}_fastqc.zip" #don't know the exact names
    params:
        "analysis_path/{sample}.{data}_fastqc.zip" #don't know the exact names
    shell:
        "fastqc {input};"
        "mv -v {params} {output};"

