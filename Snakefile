configfile: "config.yaml"

import shutil, os, sys, re
import snakemake_helper as sh
from snakemake.utils import R

FASTQDIR = config["settings"]["seqdir"]
MAPPERS = config["settings"]["mapper"]

# read in sampleinfo file with function in snakemake_helper
sampleinfo = sh.read_sampleinfo(config)

TARGETS = [sampleinfo[x]["mergefmt"][0] for x in sampleinfo]
print("Will try to create all target files for:")
print(TARGETS)

# one rule to list all targets that needs to be created, both vcfs and flagstats 
rule all:
    input:
        vcfs = expand( config["settings"]["resdir"] + "{target}.freebayes.{mapper}.vcf", target=TARGETS,mapper=MAPPERS),
        flagstats = expand( config["settings"]["resdir"] + "{target}.{mapper}.flagstat.summary", target=TARGETS,mapper=MAPPERS),
	qualimaps = expand( config["settings"]["resdir"] + "{target}.reAligned.{mapper}.qualimap/qualimapReport.html", target=TARGETS,mapper=MAPPERS),
	multiqc = config["settings"]["resdir"] + "multiqc_bwa/multiqc_report.html"

	
        
# need 2 rules for linking the fastq files in the results folder     
rule link_fastq_no_ext:
    input: FASTQDIR + "{experiment}/{sample}.{read}.trimmed.fq.gz"
    output: config["settings"]["resdir"] + "{experiment}/{sample}/{sample}.X.{read}.trimmed.fq.gz"
    shell: "ln -s {input} {output}"
         

rule link_fastq_w_ext:
    input:
        FASTQDIR + "{experiment}/{sample}.{extension}.{read}.trimmed.fq.gz"
    output:
        config["settings"]["resdir"] + "{experiment}/{sample}/{sample}.{extension}.{read}.trimmed.fq.gz"
    shell:
        "ln -s {input} {output}"


rule index_bowtie2:
    input: "{filename}"
    output: "{filename}.1.bt2"
    log: "{filename}.bowtie2-build.log"
    shell: "bowtie2-build {input} {input} > {log} 2>&1;"

rule index_bwa:
    input: "{filename}"
    output: "{filename}.bwt"
    log: "{filename}.bwa_index.log"
    shell: "bwa index {input} > {log}"

rule bwa:
    input:
        r1 = "{dir}/{sample}.r1.trimmed.fq.gz",
        r2 = "{dir}/{sample}.r2.trimmed.fq.gz",
        ref = config["ref"]["genome"],
	index = config["ref"]["genome"] + ".bwt"
    output:
        temp("{dir}/{sample}.mapped.bwa.sam")
    threads: 16	
    params:
        "-M"
    log:
        "{dir}/logs/bwa.{sample}.bwa.log",
    shell:
        "bwa mem {params} -t {threads} {input.ref} {input.r1} {input.r2} 2> {log} > {output}"


rule bowtie2:
    input:
        r1 = "{dir}/{sample}{exension}.r1.allTrimmed.fq.gz",
        r2 = "{dir}/{sample}{exension}.r2.allTrimmed.fq.gz",
        ref = config["ref"]["genome"],
        idx = config["ref"]["genome"] + ".1.bt2"
    output:
        temp("{dir}/{sample}{exension}.mapped.bowtie2.bam"),
    threads: 16	
    params:
        bowtie2 = "--maxins 2000",
	java = config["settings"]["javaopts"]
    log:
        bowtie2 = "{dir}/logs/bowtie2.{sample}{exension}.bowtie2.log",
        sam2bam = "{dir}/logs/picard.sam2bam.{sample}{exension}.bowtie2.log"
    shell:
        "bowtie2 {params.bowtie2} -p {threads} -1 {input.r1} -2 {input.r2} -x {input.ref} 2> {log.bowtie2} | "
        "picard {params.java} SamFormatConverter INPUT=/dev/stdin OUTPUT={output} > {log.sam2bam} 2>&1;"

rule sam_to_bam:
    input:
        "{dir}/{sample}.mapped.{mapper}.sam"
    output:
        temp("{dir}/{sample}.mapped.{mapper}.bam")
    params:
        java = config["settings"]["javaopts"]
    log:
        "{dir}/logs/picard.sam2bam.{sample}.{mapper}.log"
    shell:
        "picard {params.java} SamFormatConverter INPUT={input} OUTPUT={output} > {log} 2>&1;"
        


rule merge_bam:
    input:
        lambda wildcards: [ wildcards.dir +"/" +x+".mapped."+wildcards.mapper+".bam" for x in  sampleinfo[wildcards.sample]["outfmt"] ]
    output:
        merge = temp("{dir}/{sample}.merged.{mapper}.bam")
    log:
        merge = "{dir}/logs/merge.{sample}.{mapper}.log",
        index = "{dir}/logs/merge.buildbamindex.{sample}.{mapper}.log"
    params:
        java = config["settings"]["javaopts"]
    run:
        inputstr = " ".join(["INPUT={}".format(x) for x in input])
        shell("picard {param} MergeSamFiles {ips} OUTPUT={out} > {log} 2>&1".format(param=params.java, ips=inputstr, out=output.merge, log=log.merge))
        shell("picard {param} BuildBamIndex INPUT={out} > {log} 2>&1".format(param=params.java, out=output.merge, log=log.index))
        

rule filter_and_fix:
    input:
        "{dir}/{sample}.merged.{mapper}.bam"
    output:
        temp("{dir}/{sample}.fixed.{mapper}.bam")
    params:
        filters = "-b -q 2 -F 1028",
        sort = "SORT_ORDER=coordinate",
        read_groups = "CREATE_INDEX=true RGID={sample} RGLB={sample} RGPL=ILLUMINA RGSM={sample} RGCN=\"NA\" RGPU=\"NA\"",
	java = config["settings"]["javaopts"]	
    log:
        filters = "{dir}/logs/fnf.samtools.filters.{sample}.{mapper}.log",
        sort = "{dir}/logs/fnf.picard.sortsam.{sample}.{mapper}.log",
        read_groups = "{dir}/logs/fnf.picard.addorreplacereadgroup.{sample}.{mapper}.log"
    shell:
        "samtools view {params.filters} {input} 2> {log.filters} |"
        "picard {params.java} SortSam {params.sort} INPUT=/dev/stdin OUTPUT=/dev/stdout 2> {log.sort} |"
        "picard {params.java} AddOrReplaceReadGroups {params.read_groups} INPUT=/dev/stdin OUTPUT={output} 2> {log.read_groups}" 
        
rule realignertargetcreator:
    input:
        bam = "{dir}/{sample}.fixed.{mapper}.bam",
        ref = config["ref"]["genome"],
        mills =  config["ref"]["mills"],
        kgindels =  config["ref"]["kgindels"]      
    output:
        temp("{dir}/{sample}.reAlignemntTargetIntervals.{mapper}.bed")
    threads: 16
    log:
        target_creator = "{dir}/logs/realigntargetcreator.{sample}.{mapper}.log"
    shell:
        "GenomeAnalysisTK -T RealignerTargetCreator -nt {threads} -R {input.ref} -I {input.bam} -known {input.mills} -known {input.kgindels} -o {output} > {log.target_creator} 2>&1;"
        
rule realignindels:
    input:
        fixed_bam = "{dir}/{sample}.fixed.{mapper}.bam",
        targets = "{dir}/{sample}.reAlignemntTargetIntervals.{mapper}.bed",
        ref = config["ref"]["genome"],
        mills =  config["ref"]["mills"],
        kgindels =  config["ref"]["kgindels"]              
    output:
        realigned_bam = "{dir}/{sample}.reAligned.{mapper}.bam",
        buildbamindex = "{dir}/{sample}.reAligned.{mapper}.bai", 
        ginkgo_bed = "{dir}/{sample}.ginkgo.{mapper}.bed",
    params:
        java = config["settings"]["javaopts"],
    log:
        realign = "{dir}/logs/reAlign.GATK.realignindels.{sample}.{mapper}.log",
        buildbamindex = "{dir}/logs/reAlign.picard.buildbamindex.{sample}.{mapper}.log",
        bamtobed = "{dir}/logs/reAlign.bam2bed.{sample}.{mapper}.log"
    shell:
        "GenomeAnalysisTK {params.java} -T IndelRealigner -I {input.fixed_bam} -R {input.ref} -targetIntervals {input.targets} -o {output.realigned_bam} -known {input.mills} -known {input.kgindels} > {log.realign} 2>&1;"
        "picard {params.java} BuildBamIndex INPUT={output.realigned_bam} > {log.buildbamindex} 2>&1;"
        "bamToBed -i {output.realigned_bam} > {output.ginkgo_bed} 2> {log.bamtobed}"
        


                
rule freebayes:
    input:
        bam = "{dir}/{sample}.reAligned.{mapper}.bam",
        ref = config["ref"]["genome"]  
    output:
        "{dir}/{sample}.freebayes.{mapper}.vcf"
    log:
        freebayes = "{dir}/logs/freebayes.{sample}.{mapper}.log"
    shell:
        "freebayes -f {input.ref} {input.bam} -v {output} 2> {log.freebayes}"
        
rule flagstat:
    input: "{filename}.bam"
    output: "{filename}.flagstat"
    shell: "samtools flagstat {input} > {output}"
           
rule qualimap:
    input: "{filename}.bam"
    output:
      report = "{filename}.qualimap/qualimapReport.html",
      gr = "{filename}.qualimap/genome_results.txt",
      ish = "{filename}.qualimap/raw_data_qualimapReport/insert_size_histogram.txt",
      ch = "{filename}.qualimap/raw_data_qualimapReport/coverage_histogram.txt",      
      gc = "{filename}.qualimap/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt"
    log: "{filename}.qualimap/qualimap.log"
    params: "-sd -sdmode 0 --java-mem-size=20G -c -nw 400 -gd hg19"
    threads: 10
    shell: "qualimap bamqc -nt {threads} {params} -bam {input} -outdir $(dirname {output.report}) > {log} 2>&1;"

# since names are mixed up with bwa/bowtie before after the target files, make two rules for the different mappes. Obs! multiqc_bowtie is not implemented!
rule multiqc_bwa:
     input: expand( config["settings"]["resdir"] + "{target}{extensions}", target=TARGETS,extensions=config["multiqc"]["bwa_targets"])
     output: config["settings"]["resdir"] + "multiqc_bwa/multiqc_report.html"
     log: config["settings"]["resdir"] + "multiqc_bwa/multiqc.log"
     params:
        tempfile = config["settings"]["resdir"] + "/multiqc_bwa/multiqc_files.tmp",
	outdir = config["settings"]["resdir"] + "/multiqc_bwa/"
     threads: 1
     run:
       if not os.path.exists(params.outdir):
       	  os.mkdir(params.outdir)
 	  os.chmod(params.outdir,0o774)
       with open(params.tempfile,'w') as fh:
       	  for file in input:
	      	   fh.write(file)
     	      	   fh.write("\n")
       shell('multiqc -f -o {params.outdir} -l {params.tempfile} 2> {log} 1>&2')


     
rule sum_flagstat:
    input:
        mapped = lambda wildcards: [ wildcards.dir + "/" + x+".mapped."+wildcards.mapper+".flagstat" for x in  sampleinfo[wildcards.sample]["outfmt"] ],
        merged = "{dir}/{sample}.merged.{mapper}.flagstat",
	fixed = "{dir}/{sample}.fixed.{mapper}.flagstat",
        realigned = "{dir}/{sample}.reAligned.{mapper}.flagstat" 
    output:
        "{dir}/{sample}.{mapper}.flagstat.summary"
    shell:
        "cat {input.mapped} {input.merged} {input.fixed} {input.realigned} > {output}"
        

# Only run when the pipeline is run the first time. This rule is likely not complete... Implement when adding hg38
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
    params:
        java = config["settings"]["javaopts"]
    log:
        dict = config["ref"]["genome"]+".create_dict.log",
        fai = config["ref"]["genome"]+".create_fai.log",
    shell:
        "mkdir -p indexfiles;"
	"ln -s {input.REFERENCE} {output.REFERENCE};"
	"ln -s {input.MILLS} {output.MILLS};"
	"ln -s {input.KGINDELS} {output.KGINDELS};"
	"picard {params.java} SortSam {params.sort} INPUT=/dev/stdin OUTPUT=/dev/stdout 2> {log.dict};"
	"samtools faidx {input.REFERENCE} 2> {log.fai};"
        

"""
