configfile: "config.json"

rule getbundle:
    ''' Soft-link the reference  and bundle files from previusly downloaded directory
    indicated in config file'''
    output:
        REFERENCE = "indexfiles/human_g1k_v37_reference.fasta", 
        MILLS = "indexfiles/Mills_and_1000G_gold_standard.indels.b37.vcf",
        KGINDELS = "indexfiles/1000G_phase1.indels.b37.vcf"
    params:
        REFERENCE = "human_g1k_v37_reference.fasta",
        MILLS = "bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf",
        KGINDELS = "bundle/2.8/b37/1000G_phase1.indels.b37.vcf"
    shell:
        "ln -s config[gatk_bundle]/{params.REFERENCE} {output.REFERENCE}",
        "ln -s config[gatk_bundle]/{params.MILLS} {output.MILLS}",
        "ln -s config[gatk_bundle]/{params.KGINDELS} {output.KGINDELS}"



rule bowtie2:
    ''' Mapping alternative 1: using bowtie '''
    input:
        r1 ="/media/box2/Data/trimmed_fastq/{celltype}/{experiment}/{sample}.r1.allTrimmed.fq.gz",
        r2 ="/media/box2/Data/trimmed_fastq/{celltype}/{experiment}/{sample}.r2.allTrimmed.fq.gz",
        REFERENCE ="indexfiles/human_g1k_v37_reference.fasta"
    output:
        mapped_bam = "/analysis_path/{sample}/data/{sample}.mapped.bowtie2.bam",
    params:
        bowtie2 ="--maxins 2000 -p16",
        java ="-Xmx5g",
        sam2bam ="MAX_RECORDS_IN_RAM=2500000" 
    log:
        bowtie2 ="/analysis_path/{sample}/logs/bowtie2.{sample}.bowtie2.log",
        sam2bam ="/analysis_path/{sample}/logs/picard.sam2bam.{sample}.bowtie2.log"
    shell:
        "bowtie2 {params.bowtie2} -1 {input.r1} -2 {input.r2} -x {input.REFERENCE} 2> {log.bowtie2} |"
        "java {params.java} -jar /picard/SamFormatConverter.jar {params.sam2bam} INPUT=/dev/stdin OUPUT={output.mapped_bam} > {log.sam2bam} 2>&1"



rule bwa:
    ''' Mapping alternative 1: using bwa '''
    input:
        r1 ="/media/box2/Data/trimmed_fastq/{celltype}/{experiment}/{sample}.r1.allTrimmed.fq.gz",
        r2 ="/media/box2/Data/trimmed_fastq/{celltype}/{experiment}/{sample}.r2.allTrimmed.fq.gz",
        REFERENCE ="indexfiles/"human_g1k_v37_reference.fasta"
    output:
        mapped_bam = "/analysis_path/{sample}/data/{sample}.mapped.bwa.bam"
    params:
        bwa = "-M -t 16",
        java ="-Xmx5g",
        sam2bam ="MAX_RECORDS_IN_RAM=2500000"
    log:
        bwa = "/analysis_path/{sample}/logs/bwa.{sample}.bwa.log",
        sam2bam ="/analysis_path/{sample}/logs/picard.sam2bam.{sample}.bwa.log"
    shell:
        "bwa mem {params.bwa} {input.REFERENCE} {input.r1} {input.r2}  2> {log.bwa} |"
        "java {params.java} -jar /picard/SamFormatConverter.jar {params.sam2bam} INPUT=/dev/stdin OUPUT={output.mapped_bam} > {log.sam2bam} 2>&1"



If there are more than one pair fastqs, run:
rule merge:
    input:
        mapped_bam1 = "/analysis_path/{sample}/data/{sample}.mapped.{mapper}.bam",
        mapped_bam2 = "/analysis_path/{sample}/data/{sample}.mapped.{mapper}.bam" (fix here so it works with multiple bams)
    output:
        merged_bam = "/analysis_path/{sample}/data/{sample}.merged.{mapper}.bam",
        flagstat = "/analysis_path/{sample}/data/{sample}.flagstat.merged.{mapper}.txt"
    params:
        java = "-Xmx5g"
    log:
        merge = "/analysis_path/{sample}/logs/picard.merge.{sample}.{mapper}.log",
        buildbamindex = "/analysis_path/{sample}/logs/merge.picard.buildbamindex.{sample}.{mapper}.log"
    shell:
        "java {params.java} -jar picard/MergeSamFiles.jar INPUT={input.mapped_bam1} INPUT={input.mapped_bam2} OUTPUT={output.merged_bam} 1>&2  2> {log.merge}"
        "java {params.java} -jar /picard/BuildBamIndex.jar INPUT={output.merged_bam} > {logs.buildbamindex} 2>&1;"
        "samtools flagstat {output.merged_bam} > {output.flagstat}"


If there are one pair fastq, run:
rule mv:
    input:
        mapped_bam = "/analysis_path/{sample}/data/{sample}.mapped.{mapper}.bam"
    output:
        merged_bam = "/analysis_path/{sample}/data/{sample}.merged.{mapper}.bam",
        flagstat = "/analysis_path/{sample}/data/{sample}.flagstat.merged.{mapper}.txt"
    log:
        merge = "/analysis_path/{sample}/logs/mv.merge.{sample}.{mapper}.log",
        buildbamindex = "/analysis_path/{sample}/logs/merge.picard.buildbamindex.{sample}.{mapper}.log"
    shell:
        "mv -v {input.mapped_bam} {output.merged_bam} >&2 {log.merge}"
        "java {params.java} -jar /picard/BuildBamIndex.jar INPUT={output.merged_bam} > {logs.buildbamindex} 2>&1;"
        "samtools flagstat {output.merged_bam} > {output.flagstat}"
  
      
      

rule filter_and_fix:
    input:
        merged_bam = "/analysis_path/{sample}/data/{sample}.merged.{mapper}.bam
    output:
        fixed_bam = "/analysis_path/{sample}/data/{sample}.fixed.{mapper}.bam"
    params:
        filters ="b -q 2 -F 1028",
        java = "-Xmx5g",
        sortsam = "MAX_RECORDS_IN_RAM=2500000 SORT_ORDER=coordinate",
        addorReplaceReadgroup = "MAX_RECORDS_IN_RAM=2500000 CREATE_INDEX=true RGID=SAMPLE RGLB=SAMPLE RGPL=ILLUMINA RGSM=SAMPLE RGCN=\"NA\" RGPU=\"NA\"" 
    log:
        filters = "/analysis_path/{sample}/logs/fnf.samtools.filters.{sample}.{mapper}.log",
        sortsam = "/analysis_path/{sample}/logs/fnf.picard.sortsam.{sample}.{mapper}.log",
        addorreplacereadgroup = "/analysis_path/{sample}/logs/fnf.picard.addorreplacereadgroup.{sample}.{mapper}.log"
    shell:
        "samtools view {params.filters} {input.merged_bam} 2> {log.samtools} |"
        "java {params.java} /picard/sortsam.jar {params.sortsam} INPUT=/dev/stdin 2> {log.sortsam} |"
        "java {params.java} /picard/addOrReplaceReadGroups.jar {params.addorreplacereadgroup} INPUT=/dev/stdin OUTPUT={output.fixed_bam} 2> {log.addorreplacereadgroup}"
        

        
rule realignertargetcreator:
    input:
        fixed_bam = "/analysis_path/{sample}/data/{sample}.fixed.{mapper}.bam",
        REFERENCE = "human_g1k_v37_reference.fasta",
        MILLS = "bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf",
        KGINDELS = "bundle/2.8/b37/1000G_phase1.indels.b37.vcf"
    output:
        targets = "/analysis_path/{sample}/data/{sample}.reAlignemntTargetIntervals.{mapper}.bed"
    params:
        java = "-Xmx5g",
        target_creator = "-nt 16"
    log:
        target_creator = "/analysis_path/{sample}/logs/realigntargetcreator.{sample}.{mapper}.log"
    shell:
        "java {params.java} -jar GenomeAnalysisTK.jar -T RealignerTargetCreator {params.target_creator} -R {input.REFERENCE} -I {input.fixed_bam} -known {input.MILLS} -known {input.KGINDELS} -o {output.targets} > {log.target_creator} 2> &1;"


rule realignindels:
    input:
        fixed_bam = "/analysis_path/{sample}/data/{sample}.fixed.{mapper}.bam",
        REFERENCE = "human_g1k_v37_reference.fasta",
        MILLS = "bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf",
        KGINDELS = "bundle/2.8/b37/1000G_phase1.indels.b37.vcf",
        targets = "/analysis_path/{sample}/data/{sample}.reAlignemntTargetIntervals.{mapper}.bed"
    output:
        realigned_bam = "/analysis_path/{sample}/data/{sample}.reAligned.{mapper}.bam",
#        buildbamindex = "", #no output for this???? It is written so that an index file (.bai) is put in the same dict and with the same name as bam
        flagstat = "/analysis_path/{sample}/data/{sample}.flagstat.realigned.{mapper}.txt",
        ginkgo_bed = "/analysis_path/{sample}/data/{sample}.ginkgo.{mapper}.bed"
    params:
        java = "-Xmx5g"
        
    log:
        realign = "/analysis_path/{sample}/logs/reAlign.GATK.realignindels.{sample}.{mapper}.log",
        buildbamindex = "/analysis_path/{sample}/logs/reAlign.picard.buildbamindex.{sample}.{mapper}.log",
        bamtobed = "/analysis_path/{sample}/logs/reAlign.bam2bed.{sample}.{mapper}.log"
    shell:
        "java {params.java} -jar GenomeAnalysisTK.jar -T IndelRealigner -I {input.fixed_bam} -R {input.REFERENCE} -targetIntervals {input.targets} -o {output.realigned_bam} -known {input.MILLS} -known {input.KGINDELS} 1>&2 2> {log.realign} 2>&1;"
        "java {params.java} -jar /picard/BuildBamIndex.jar INPUT={output.realigned_bam} > {logs.buildbamindex} 2>&1;"
        "samtools flagstat {output.realigned_bam} > {output.flagstat}"
        "bamToBed -i {output.realigned_bam} > {output.ginkgo_bed} 2> {log.bamtoted}"
        
                
rule freebayes:
    input:
        realigned_bam = "/analysis_path/{sample}/data/{sample}.reAligned.{mapper}.bam",
        REFERENCE = "human_g1k_v37_reference.fasta"
    output:
        vcf = "/analysis_path/{sample}/data/{sample}.freebayes.{mapper}.vcf"
    log:
        freebayes = "/analysis_path/{sample}/logs/freebayes.{sample}.{mapper}.log"
    shell:
        "freebayes -f {input.REFERENCE} {input.realigned_bam} {output.vcf} 2> {log.freebayes}"
        

        


# how much memory
# global variable ref and bundle and different bams
# /dev/stdin??? (picard input)
# and is dev/stdin updated in rule filter
# 2>, 1>&2  2>
# var kommer {mapper} ifrån
# add replace read group might not work pga SAMPLE
# bamtobed needs samtools
# can MAX_RECORDS_IN_RAM=2500000 be removed from picard
# what is ; (in realignindels)
# rm bams after next step
#can param and log variable have same name?
