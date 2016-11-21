configfile: "config.yaml"
import snakemake_helper as sh

FASTQDIR = config["settings"]["seqdir"]
MAPPERS = config["settings"]["mapper"]

# read in sampleinfo file with function in snakemake_helper
sampleinfo = sh.read_sampleinfo(config)
     
TARGETS = [sampleinfo[x]["mergefmt"][0] for x in sampleinfo]
print("Will try to create all target files for:")
print(TARGETS)

# one rule to list all targets that needs to be created
rule all_vcf:
    input: expand( config["settings"]["resdir"] + "{target}.freebayes.{mapper}.vcf", target=TARGETS,mapper=MAPPERS)

# one rule to list all flagstat summary files
rule all_flagstat:
    input: expand( config["settings"]["resdir"] + "{target}.{mapper}.flagstat.summary", target=TARGETS,mapper=MAPPERS)

          
# need 2 rules for linking the fastq files in the results folder     
rule link_fastq_no_ext:
     input: FASTQDIR + "{experiment}/{sample}.{read}.allTrimmed.fq.gz"
     output: config["settings"]["resdir"] + "{experiment}/{sample}/{sample}.X.{read}.allTrimmed.fq.gz"
     shell: "ln -s {input} {output}"
         

rule link_fastq_w_ext:
     input:
        FASTQDIR + "{experiment}/{sample}.{extension}.{read}.allTrimmed.fq.gz"
     output:
        config["settings"]["resdir"] + "{experiment}/{sample}/{sample}.{extension}.{read}.allTrimmed.fq.gz"
     shell:
        "ln -s {input} {output}"

rule index_bowtie2:
    input: "{filename}"
    output: "{filename}.1.bt2"
    log: "{filename}.bowtie2-build.log"
    shell: "bowtie2-build {input} {input} > {log}"

rule index_bwa:
    input: "{filename}"
    output: "{filename}.bwt"
    log: "{filename}.bwa_index.log"
    shell: "bwa index {input} > {log}"

rule bwa:
    input:
        r1 = "{sample}.r1.allTrimmed.fq.gz",
        r2 = "{sample}.r2.allTrimmed.fq.gz",
        ref = config["ref"]["genome"],
	index = config["ref"]["genome"] + ".bwt"
    output:
        "{sample}.mapped.bwa.sam"
    params:
        "-M -t 16"
    log:
        "logs/bwa.{sample}.bwa.log",
    shell:
        "bwa mem {params} {input.ref} {input.r1} {input.r2} 2> {log} > {output}"


rule bowtie2:
    input:
        r1 = "{sample}{exension}.r1.allTrimmed.fq.gz",
        r2 = "{sample}{exension}.r2.allTrimmed.fq.gz",
        ref = config["ref"]["genome"],
        idx = config["ref"]["genome"] + ".1.bt2"
    output:
        "{sample}{exension}.mapped.bowtie2.bam",
    params:
        bowtie2 = "--maxins 2000 -p16"
    log:
        bowtie2 = "logs/bowtie2.{sample}{exension}.bowtie2.log",
        sam2bam = "logs/picard.sam2bam.{sample}{exension}.bowtie2.log"
    shell:
        "bowtie2 {params.bowtie2} -1 {input.r1} -2 {input.r2} -x {input.ref} 2> {log.bowtie2} | "
        "picard SamFormatConverter INPUT=/dev/stdin OUPUT={output} > {log.sam2bam} "

rule sam_to_bam:
    input:
        "{sample}.mapped.{mapper}.sam"
    output:
        "{sample}.mapped.{mapper}.bam"
    log:
        "logs/picard.sam2bam.{sample}.{mapper}.log"
    shell:
        "picard SamFormatConverter INPUT={input} OUTPUT={output} > {log} 2>&1;"

        

# COMMENT (merge): did we need to change the command for running picard or does java -jar work? - Removed java -jar from all!
# COMMENT (merge); Do we really want to do bam index within the merge rule, do we not have to sort as well at some step, perhaps we can index within that rule?

rule merge_bam:
    input:
        lambda wildcards: [ wildcards.dir +"/" +x+".mapped."+wildcards.mapper+".bam" for x in  sampleinfo[wildcards.sample]["outfmt"] ]
    output:
        "{dir}/{sample}.merged.{mapper}.bam"
    log:
        merge = "{dir}/logs/merge.{sample}.{mapper}.log",
        index = "{dir}/logs/merge.buildbamindex.{sample}.{mapper}.log"
    run: 
      if (len(input) > 1):
          inputstr = " ".join(["INPUT={}".format(x) for x in input])
          shell("picard MergeSamFiles {ips} OUTPUT={out} > {log}".format(param=params.java, ips=inputstr, out=output.merge, log=log.merge))
      else:
          if os.path.exists(output.merge):
              os.unlink(output.merge)
          shutil.copy(input[0], output.merge)
      shell("picard BuildBamIndex INPUT={out} > {log}".format(param=params.java, out=output.merge, log=log.index))     
     
     
"""


#old merge with shell only
rule merge:
    input:
        #mapped_bams = lambda wildcards: config["metafiles"][wildcards.sample]
        mapped_bams = ["results/{experiment}/{sample}/data/" + s + ".mapped.bwa.bam" for s in lambda wildcards: config["metafiles"][wildcards.sample]]
    output:
        merged_bam = "results/{experiment}/{sample}/data/{sample}.merged.{mapper}.bam",
        flagstat = "results/{experiment}/{sample}/data/{sample}.flagstat.merged.{mapper}.txt"
    params:
        java = "-Xmx5g",
        nbamfiles = lambda wildcards: len(config["metafiles"][wildcards.sample]),
        picardinput = lambda wildcards: ["INPUT="+s for s in config["metafiles"][wildcards.sample]]
    log:
        merge = "results/{experiment}/{sample}/logs/merge.{sample}.{mapper}.log",
        buildbamindex = "results/{experiment}/{sample}/logs/merge.buildbamindex.{sample}.{mapper}.log"
    shell:
        '''
        if test {params.nbamfiles} -eq 1; then 
           mv -v {input.mapped_bams} {output.merged_bam} >&2 {log.merge};
        else
           java {params.java} -jar picard/MergeSamFiles.jar {params.picardinput} OUTPUT={output.merged_bam} 1>&2  2> {log.merge}
        fi'
        java {params.java} -jar /picard/BuildBamIndex.jar INPUT={output.merged_bam} > {log.buildbamindex} 2>&1
        samtools flagstat {output.merged_bam} > {output.flagstat}
        '''
"""


rule filter_and_fix:
    input:
        "{sample}.merged.{mapper}.bam"
    output:
        "{sample}.fixed.{mapper}.bam"
    params:
        filters = "b -q 2 -F 1028",
        sort = "SORT_ORDER=coordinate",
        read_groups = "CREATE_INDEX=true RGID=SAMPLE RGLB=SAMPLE RGPL=ILLUMINA RGSM=SAMPLE RGCN=\"NA\" RGPU=\"NA\"" 
    log:
        filters = "logs/fnf.samtools.filters.{sample}.{mapper}.log",
        sort = "logs/fnf.picard.sortsam.{sample}.{mapper}.log",
        read_groups = "logs/fnf.picard.addorreplacereadgroup.{sample}.{mapper}.log"
    shell:
        "samtools view {params.filters} {input} 2> {log.filters} |"
        "picard SortSam.jar {params.sort} INPUT=/dev/stdin 2> {log.sort} |"
        "picard addOrReplaceReadGroups.jar {params.read_groups} INPUT=/dev/stdin OUTPUT={output} 2> {log.read_groups}" 
        
rule realignertargetcreator:
    input:
        bam = "{sample}.fixed.{mapper}.bam",
        ref = config["ref"]["genome"],
        mills =  config["ref"]["mills"],
        kgindels =  config["ref"]["kgindels"]      
    output:
        "{sample}.reAlignemntTargetIntervals.{mapper}.bed"
    params:
        "-nt 16"
    log:
        target_creator = "logs/realigntargetcreator.{sample}.{mapper}.log"
    shell:
        "GenomeAnalysisTK.jar -T RealignerTargetCreator {params} -R {input.ref} -I {input.bam} -known {input.mills} -known {input.kgindels} -o {output} > {log.target_creator} 2> &1;"
        

rule realignindels:
    input:
        fixed_bam = "{sample}.fixed.{mapper}.bam",
        targets = "{sample}.reAlignemntTargetIntervals.{mapper}.bed",
        ref = config["ref"]["genome"],
        mills =  config["ref"]["mills"],
        kgindels =  config["ref"]["kgindels"]              
    output:
        realigned_bam = "{sample}.reAligned.{mapper}.bam",
        buildbamindex = "{sample}.reAligned.{mapper}.bai", 
        ginkgo_bed = "{sample}.ginkgo.{mapper}.bed"
    log:
        realign = "logs/reAlign.GATK.realignindels.{sample}.{mapper}.log",
        buildbamindex = "logs/reAlign.picard.buildbamindex.{sample}.{mapper}.log",
        bamtobed = "logs/reAlign.bam2bed.{sample}.{mapper}.log"
    shell:
        "GenomeAnalysisTK.jar -T IndelRealigner -I {input.fixed_bam} -R {input.ref} -targetIntervals {input.targets} -o {output.realigned_bam} -known {input.mills} -known {input.kgindels} 1>&2 2> {log.realign} 2>&1;"
        "picard BuildBamIndex.jar INPUT={output.realigned_bam} > {log.buildbamindex} 2>&1;"
        "bamToBed -i {output.realigned_bam} > {output.ginkgo_bed} 2> {log.bamtobed}"

# COMMENT: Too many commands in one rule?

                
rule freebayes:
    input:
        bam = "{sample}.reAligned.{mapper}.bam",
        ref = config["ref"]["genome"]  
    output:
        "{sample}.freebayes.{mapper}.vcf"
    log:
        freebayes = "logs/freebayes.{sample}.{mapper}.log"
    shell:
        "freebayes -f {input.ref} {input.bam} {output} 2> {log.freebayes}"

rule flagstat:
     input: "{filename}.bam"
     output: "{filename}.flagstat"
     shell: "samtools flagstat {input} > {output}"


#COMMENT: this is a mock rule that only runs cat of all the flagstat runs
# for now only created as a test to see if we can force snakemake to run flagstat at multiple places.. 
     
rule sum_flagstat:
     input:
          mapped = lambda wildcards: [ wildcards.dir + "/" + x+".mapped."+wildcards.mapper+".flagstat" for x in  sampleinfo[wildcards.sample]["outfmt"] ],
          merged = "{dir}/{sample}.merged.{mapper}.flagstat",
          realigned = "{dir}/{sample}.reAligned.{mapper}.flagstat" 
     output:
          "{dir}/{sample}.{mapper}.flagstat.summary"
     shell:
          "cat {input.mapped} {input.merged} {input.realigned} > {output}"


rule getbundle:
    ''' Soft-link the reference  and bundle files from previusly downloaded directory
        indicated in config file'''
      input:
        REFERENCE = config["bundle"]["genome"],
        MILLS = config["bundle"]["mills"],
        KGINDELS = config["bundle"]["kgindels"],
      output:
        REFERENCE = config["ref"]["genome"],
        MILLS = config["ref"]["mills"],
        KGINDELS = config["ref"]["kgindels"],
      shell:
        "mkdir -p indexfiles;"
	"ln -s {input.REFERENCE} {output.REFERENCE};"
	"ln -s {input.MILLS} {output.MILLS};"
	"ln -s {input.KGINDELS} {output.KGINDELS};"

"""
OLD version OBS!
rule getbundle:
    ''' Soft-link the reference  and bundle files from previously downloaded directory
    indicated in config file'''
    output:
	 
#        REFERENCE = config["ref"]["genome"]
#        MILLS = config["ref"]["mills"]
#        KGINDELS = config["ref"]["kgindels"]
    params:
        REFERENCE = config["bundle"]["genome"]
        MILLS = config["bundle"]["mills"]
 	KGINDELS = config["bundle"]["kgindels"]
    shell:
        "mkdir -p indexfiles;"
	"ln -s {params.REFERENCE} {output.REFERENCE};"
        "ln -s {params.MILLS} {output.MILLS};"
        "ln -s {params.KGINDELS} {output.KGINDELS};"


"""
