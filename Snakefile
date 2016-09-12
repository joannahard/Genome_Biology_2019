configfile = config.json

rule all:
    input:
       "analysis_path/{sample}/{sample}.r1.allTrimmed.fastqc", #don't know the exact names
       "analysis_path/{sample}/{sample}.r2.allTrimmed.fastqc",
       "analysis_path/{sample}/{sample}.singletts.fastqc"

rule wgaAdaptorTrimmer:
    input:
        "{read}.fastq.gz"
    output:
        "temp(config[tmpdir]/{sample}.{read}.wgaTrimmed2.fq)"
    log:
        "analysis_path/log/RubiconWgaTrimming.{sample}.{read}.log.txt"
    shell:
        "wgaAdapterTrimmer.py -i {input} > {output} 2> {log};"

rule cutadapt:
    input:
        "config[tmpdir]/{sample}.{read}.wgaTrimmed2.fq"
    output:
        "temp(config[tmpdir]/{sample}.{read}.wgaAndilluminaTrimmed.fq)"
    params:
        "expand("-a", config[adapterSeq])"
    log:
        "analysis_path/logs/illuminaAndNexteraTrimming.{sample}.{read}.log.txt"
    shell:
        "cutadapt -n 3 {params} {input} > {output} 2> {log};"

rule TrimBWAstyle:
    input:
        "config[tmpdir]/{sample}.{read}.wgaAndilluminaTrimmed.fq"
    output:
        "temp(config[tmpdir]/{sample}.{read}.wgaIlluminaAndQualityTrimmed.fq)"
    log:
        "analysis_path/loogs/qualityTrimming.{sample}.{read}.log.txt"
    shell:
        "TrimBWAstyle.pl {input} > {output} 2> {log};"

rule removeEmptyReads:
    input:
        "config[tmpdir]/{sample}.r1.wgaIlluminaAndQualityTrimmed.fq",
        "config[tmpdir]/{sample}.r1.wgaIlluminaAndQualityTrimmed.fq"        
    output:
        "analysis_path/{sample}.r1.allTrimmed.fq.gz",
        "analysis_path/{sample}.r2.allTrimmed.fq.gz",
        "analysis_path/{sample}.singletts.fq.gz"
    params:
        "analysis_path/{sample}.r1.allTrimmed.fq",
        "analysis_path/{sample}.r2.allTrimmed.fq",
        "analysis_path/{sample}.singletts.fq"
    log:
        "/analysis_path/removeEmptyReads.{sample}.log.txt"
    shell:
        "python removeEmptyReads.py {input} {params} 2> {log};"
        "gzip analysis_path/{sample}.*.fq"

rule fastqc:
    input:
        "analysis_path/{sample}.{data}.fq.gz"
    output:
        "analysis_path/{sample}/{sample}.{data}.fastqc" #don't know the exact names
    params:
        "analysis_path/{sample}.{data}.fastqc" #don't know the exact names
    shell:
        "fastqc {input};"
        "mv -v }params} {output};"

