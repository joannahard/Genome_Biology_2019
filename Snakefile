configfile: "config.json"
from snakemake_helper import getMapperInput, getMergeInput


REF = config["ref"]
MILLS = config["mills"]
KGINDELS = config["kgindels"]


rule bwa:
    ''' Mapping alternative 1: using bowtie '''
    input:
        r1 = lambda wildcards: getMapperInput(config, wildcards, "r1"),
        r2 = lambda wildcards: getMapperInput(config, wildcards, "r2"),
    output:
        "results/{experiment}/{sample}/data/{sample}{extension}.mapped.bwa.sam"
    params:
        "-M -t 16"
    log:
        "results/{experiment}/{sample}/logs/bwa.{sample}{extension}.bwa.log",
    shell:
        "bwa mem {params} {REF} {input.r1} {input.r2} 2> {log} > {output}"

#runs for bwa. not bowtie because it is missing a bowtie2 index
# make tmp rules

rule bowtie2:
    ''' Mapping alternative 1: using bowtie '''
    input:
        r1 = lambda wildcards: getMapperInput(config, wildcards, "r1"),
        r2 = lambda wildcards: getMapperInput(config, wildcards, "r2"),
    output:
        "results/{experiment}/{sample}/data/{sample}{extension}.mapped.bowtie2.bam",
    params:
        bowtie2 = "--maxins 2000 -p16"
    log:
        bowtie2 = "results/{experiment}/{sample}/logs/bowtie2.{sample}{extension}.bowtie2.log",
        sam2bam = "results/{experiment}/{sample}/logs/picard.sam2bam.{sample}{extension}.bowtie2.log"
    shell:
        "bowtie2 {params.bowtie2} -1 {input.r1} -2 {input.r2} -x {REF} 2> {log.bowtie2} | "
        "picard SamFormatConverter.jar INPUT=/dev/stdin OUPUT={output} > {log.sam2bam} "

# Make rule for bowtie index




rule sam_to_bam:
    input:
        "results/{experiment}/{sample}/data/{sample}{extension}.mapped.{mapper}.sam"
    output:
        "results/{experiment}/{sample}/data/{sample}{extension}.mapped.{mapper}.bam"
    log:
        "results/{experiment}/{sample}/logs/picard.sam2bam.{sample}{extension}.{mapper}.log"
    shell:
        "picard SamFormatConverter INPUT={input} OUTPUT={output} > {log} 2>&1;"
        




rule merge:
    input:
        #mapped_bams = lambda wildcards: config["metafiles"][wildcards.sample]
#        mapped_bams = ["results/{experiment}/{sample}/data/" + s + ".mapped.bwa.bam" for s in lambda wildcards: config["metafiles"][wildcards.sample]]
       mapped_bams = lambda wildcards: getMergeInput(config, wildcards)
    output:
        merged_bam = "results/{experiment}/{sample}/data/{sample}.merged.{mapper}.bam",
        flagstat = "results/{experiment}/{sample}/data/{sample}.flagstat.merged.{mapper}.txt"
    params:
        java = "-Xmx5g",
        nbamfiles = lambda wildcards: len(config["metafiles"][wildcards.sample]),
#        picardinput = lambda wildcards: ["INPUT="+s for s in config["metafiles"][wildcards.sample]]
    log:
        merge = "results/{experiment}/{sample}/logs/merge.{sample}.{mapper}.log",
        buildbamindex = "results/{experiment}/{sample}/logs/merge.buildbamindex.{sample}.{mapper}.log"
    shell:
        '''
        if test {params.nbamfiles} -eq 1; then 
           cp {input.mapped_bams} {output.merged_bam} >&2 {log.merge};
        else
           java {params.java} -jar picard/MergeSamFiles.jar {params.picardinput} OUTPUT={output.merged_bam} 1>&2  2> {log.merge}
        fi
        java {params.java} -jar /picard/BuildBamIndex.jar INPUT={output.merged_bam} > {log.buildbamindex} 2>&1
        samtools flagstat {output.merged_bam} > {output.flagstat}
        '''

"""

rule filter_and_fix:
    input:
        "results/{experiment}/{sample}/data/{sample}.merged.{mapper}.bam"
    output:
        "results/{experiment}/{sample}/data/{sample}.fixed.{mapper}.bam"
    params:
        filters = "b -q 2 -F 1028",
        sort = "SORT_ORDER=coordinate",
        read_groups = "CREATE_INDEX=true RGID=SAMPLE RGLB=SAMPLE RGPL=ILLUMINA RGSM=SAMPLE RGCN=\"NA\" RGPU=\"NA\"" 
    log:
        filters = "results/{experiment}/{sample}/logs/fnf.samtools.filters.{sample}.{mapper}.log",
        sort = "results/{experiment}/{sample}/logs/fnf.picard.sortsam.{sample}.{mapper}.log",
        read_groups = "results/{experiment}/{sample}/logs/fnf.picard.addorreplacereadgroup.{sample}.{mapper}.log"
    shell:
        "samtools view {params.filters} {input} 2> {log.filters} |"
        "picard SortSam.jar {params.sort} INPUT=/dev/stdin 2> {log.sort} |"
        "picard addOrReplaceReadGroups.jar {params.read_groups} INPUT=/dev/stdin OUTPUT={output} 2> {log.read_groups}" 
        


rule realignertargetcreator:
    input:
        "results/{experiment}/{sample}/data/{sample}.fixed.{mapper}.bam"
    output:
        "results/{experiment}/{sample}/data/{sample}.reAlignemntTargetIntervals.{mapper}.bed"
    params:
        "-nt 16"
    log:
        target_creator = "results/{experiment}/{sample}/logs/realigntargetcreator.{sample}.{mapper}.log"
    shell:
        "GenomeAnalysisTK.jar -T RealignerTargetCreator {params} -R {REF} -I {input} -known {MILLS} -known {KGINDELS} -o {output} > {log.target_creator} 2> &1;"
        



rule realignindels:
    input:
        fixed_bam = "results/{experiment}/{sample}/data/{sample}.fixed.{mapper}.bam",
        targets = "results/{experiment/{sample}/data/{sample}.reAlignemntTargetIntervals.{mapper}.bed"
    output:
        realigned_bam = "results/{experiment}/{sample}/data/{sample}.reAligned.{mapper}.bam",
        buildbamindex = "results/{experiment}/{sample}/data/{sample}.reAligned.{mapper}.bai", 
        flagstat = "results/{experiment}/{sample}/data/{sample}.flagstat.realigned.{mapper}.txt",
        ginkgo_bed = "results/{experiment}/{sample}/data/{sample}.ginkgo.{mapper}.bed"
    log:
        realign = "results/{experiment}/{sample}/logs/reAlign.GATK.realignindels.{sample}.{mapper}.log",
        buildbamindex = "results/{experiment}/{sample}/logs/reAlign.picard.buildbamindex.{sample}.{mapper}.log",
        bamtobed = "results/{experiment}/{sample}/logs/reAlign.bam2bed.{sample}.{mapper}.log"
    shell:
        "GenomeAnalysisTK.jar -T IndelRealigner -I {input.fixed_bam} -R {REF} -targetIntervals {input.targets} -o {output.realigned_bam} -known {MILLS} -known {KGINDELS} 1>&2 2> {log.realign} 2>&1;"
        "picard BuildBamIndex.jar INPUT={output.realigned_bam} > {logs.buildbamindex} 2>&1;"
        "samtools flagstat {output.realigned_bam} > {output.flagstat}"
        "bamToBed -i {output.realigned_bam} > {output.ginkgo_bed} 2> {log.bamtoted}"
        
                
rule freebayes:
    input:
        "results/{experiment}/{sample}/data/{sample}.reAligned.{mapper}.bam"
    output:
        "results/{experiment}/{sample}/data/{sample}.freebayes.{mapper}.vcf"
    log:
        freebayes = "results/{experiment}/{sample}/logs/freebayes.{sample}.{mapper}.log"
    shell:
        "freebayes -f {REF} {input} {output} 2> {log.freebayes}"







rule getbundle:
    ''' Soft-link the reference  and bundle files from previously downloaded directory
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
        "{ref}"
        "ln -s config[gatk_bundle]/{params.REFERENCE} {output.REFERENCE};"
        "ln -s config[gatk_bundle]/{params.MILLS} {output.MILLS};"
        "ln -s config[gatk_bundle]/{params.KGINDELS} {output.KGINDELS};"

"""




#rule merge:
#    input:
#        #mapped_bams = ["results/{experiment}/{sample}/data/" + s + ".mapped.bwa.bam" for s in ["B4C6#.5","B4C6.6"]]
#        #mapped_bams = ["results/{experiment}/{sample}/data/" + s + ".mapped.bwa.bam" for s in lambda# wildcards: config["metafiles"][wildcards.sample]]
##        mapped_bams = lambda wildcards: ["results/{wildcards.experiment}/{wildcards.sample}/data/" +# s + ".mapped.bwa.bam" for s in config["metafiles"][wildcards.sample]]
#    output:
#        test = "results/{experiment}/{sample}/data/{sample}.{mergetest"
#    params: 
#        names = lambda wildcards: ["INPUT="+s for s in config["metafiles"][wildcards.sample]]
#    shell:        
#        "echo {names} > test"
#"B4C6": ["results/12MDAs/B4C6/data/B4C6.5.mapped.bwa.bam", "results/12MDAs/B4C6/data/B4C6.6.mapped.bwa.bam"],
