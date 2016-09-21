configfile: "config.json"

rule getbundle:
    ''' Soft-link the decoy  and bundle files from previusly downloaded directory
    indicated in config file'''
    output:
        decoy = "indexfiles/human_g1k_v37_decoy.fasta", 
        bundle1 = "indexfiles/Mills_and_1000G_gold_standard.indels.b37.vcf",
        bundle2 = "indexfiles/1000G_phase1.indels.b37.vcf"
    params:
        decoy = "human_g1k_v37_decoy.fasta",
        bundle1 = "bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf",
        bundle2 = "bundle/2.8/b37/1000G_phase1.indels.b37.vcf"
    shell:
        "ln -s config[gatk_bundle]/{params.decoy} {output.decoy}",
        "ln -s config[gatk_bundle]/{params.bundle1} {output.bundle1}",
        "ln -s config[gatk_bundle]/{params.bundle2} {output.bundle2}"

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
        "/analysis_path/{sample}/data/{sample}.fixed.{mapper}.bam"
    params:
        samtools="b -q 2 -F 1028",
        java = "-Xmx5g",
        sortsam = "MAX_RECORDS_IN_RAM=2500000 SORT_ORDER=coordinate" INPUT=/dev/stdin ,
        addorreplacereadgroup = "MAX_RECORDS_IN_RAM=2500000 INPUT=/dev/stdin CREATE_INDEX=true RGID=SAMPLE RGLB=SAMPLE RGPL=ILLUMINA RGSM=SAMPLE RGCN=\¶"NA\" RGPU=\"NA\""
    log:
        samtools = "logs/filter.samtools.{sample}.{mapper}.log",
        sortsam = "logs/filter.sortsam.{sample}.{mapper}.log",
        sortsam = "logs/filter.addorreplacereadgroup.{sample}.{mapper}.log"
    shell:
        "samtools view {params.samtools} {input} 2> {log.samtools} |"
        "java {params.java} /picard/SortSam.jar {params.SortSam} 2> {log.sortsam} |"
        "java {params.java} /picard/AddOrReplaceReadGroups.jar {params.addorreplacereadgroup} OUTPUT={output} 2> {log.addorreplacereadgroup}"
        

        
rule realignertargetcreator:
    input:
        decoy = "human_g1k_v37_decoy.fasta",
        data = "/analysis_path/{sample}/data/{sample}.fixed.{mapper}.bam",
        bundle1 = "bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf",
        bundle2 = "bundle/2.8/b37/1000G_phase1.indels.b37.vcf"
    output:
        "/analysis_path/{sample}/data/{sample}.reAlignemntTargetIntervals.{mapper}.bed"
    params:
        java = "-Xmx5g",
        prog = "-T RealignerTargetCreator -nt 16"
    log:
        "logs/realigntargetcreator.{sample}.{mapper}.log"
    shell:
        "java {params.java} -jar GenomeAnalysisTK.jar {params.prog} -R {input.decoy} -I {input.data} -known {input.bundle1} -known {input.bundle1} -o {output} > {log} 2> &1;"


rule realignindels:
    input:
        data = "/analysis_path/{sample}/data/{sample}.fixed.{mapper}.bam",
        decoy = "human_g1k_v37_decoy.fasta",
        bundle1 = "bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf",
        bundle2 = "bundle/2.8/b37/1000G_phase1.indels.b37.vcf",
        target = "/analysis_path/{sample}/data/{sample}.reAlignemntTargetIntervals.{mapper}.bed"
    output:
        genomeanalysis = "/analysis_path/{sample}/data/{sample}.reAligned.bam",
#         buildbamindex = "", #no output for this????
        flagstat = "/analysis_path/{sample}/logs/Flagstat.{sample}.txt",
        bamtobed = "/analysis_path/{sample}/data/{sample}.bed"
    params:
        java = "-Xmx5g",
        genomeanalysis = "-T IndelRealigner}"
    log:
        genomeanalysis = "logs/realignindels.genomeanalysis.{sample}.{mapper}.log",
        buildbamindex = "logs/realignindels.buildbamindex.{sample}.{mapper}.log",
        flagstat = "logs/realignindels.flagstat.{sample}.{mapper}.log",
        bamtobed = "logs/realignindels.bamtobed.{sample}.{mapper}.log",
    shell:
        "java {params.java} -jar GenomeAnalysisTK.jar {params.genmeanalysis} -I {input.data} -R {input.decoy}  -known {input.bundle1} -known {input.bundle2} -o {output.genomeanalysis} > {log.genomeanalysis} 2>&1;"
        "java {params.java} -jar /picard/BuildBamIndex.jar INPUT={output.genomeanalysis} > {logs.buildbamindex} 2>&1;"
        "samtools flagstat {output.genomeanalysis} > output.flagstat 2> {log.flagstat};"
        "bamToBed -i {output.genomeanalysis} > /analysis_path/SAMPLE{sample}/data/SAMPLE{sample}·{mapper}.bed 2> {log.bamToBed}"
