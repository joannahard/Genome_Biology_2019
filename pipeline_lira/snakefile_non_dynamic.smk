import os,sys

# new version without dynamic rules. After all splits + plink are run. Read in subdirs for all chroms of all samples. 

configfile: "config_pipeline.yaml"
            
BAMDIR = config["bamdir"]
CHROMOSOMES = [str(x) for x in config["chromosomes"]]
OUTDIR = config["resdir"]


S = os.listdir(OUTDIR)
SAMPLES = [x for x in S if x in [str(y) for y in config["samples"]]]
SAMPLES_wBULK = [x for x in S if x in [str(y) for y in config["samples"] + [config["bulk"]]]]
#print(SAMPLES)
#print(SAMPLES_wBULK)

sample_dict = {}
for s in SAMPLES_wBULK:
    sample = str(s)
    sample_dict[sample] = {}
    for chr in CHROMOSOMES:
    	jobs = os.listdir(OUTDIR + sample + "/" + chr + "/jobs")
	sample_dict[sample][chr]=jobs
        
rule all:
    input: expand(OUTDIR + "{sample}/varcall_"+ config["bulk"] + "/powered.regions.bed", sample = SAMPLES)


# plink:
# Instead of their plink script, run submission of all the shell scripts in outdir/sample/job_scripts/
# has names like {chrom}_{jobid}.sh
# OBS! For next step, need to have a checkfile in folder outdir/sample/progress/
# has name like .{chrom}_{jobid}
# track final file vcf.info.rda in job folder.

rule plink:
    input:
        script = OUTDIR + "{sample}/job_scripts/{chrom}_{region_id}.sh",
        sites_file = OUTDIR + "{sample}/{chrom}/jobs/{region_id}/sites.bed"
    output:
        chkfile = OUTDIR + "{sample}/progress/.{chrom}_{region_id}",
        vcf_rda = OUTDIR + "{sample}/{chrom}/jobs/{region_id}/vcf.info.rda"
    shell:
        "bash {input.script};"
        "touch {output.chkfile}"

rule chrom_plink:
    input: lambda wc: expand(OUTDIR + "{{sample}}/{{chrom}}/jobs/{region_id}/vcf.info.rda", region_id = sample_dict[str(wc.sample)][wc.chrom])
    output: OUTDIR + "{sample}/progress/plink_{chrom}.chk"
    shell:
        "touch {output}"

        
           
"""
# compare - run for all samples vs bulk.
Creates:
folder: output/110/1/compare/
-rw-rw-r-- 1 asa asa 2127952 Nov  1 10:44 single.cell.linkage.rda
-rw-rw-r-- 1 asa asa 4508121 Nov  1 10:44 single.cell.vcf.info.rda
-rw-rw-r-- 1 asa asa     127 Nov  1 10:44 log.txt
-rw-rw-r-- 1 asa asa 5940631 Nov  1 10:45 bulk.linkage.rda
-rw-rw-r-- 1 asa asa 7149651 Nov  1 10:45 bulk.vcf.info.rda
-rw-rw-r-- 1 asa asa   36965 Nov  1 10:45 combined.somatic.bulk.rda
-rw-rw-r-- 1 asa asa    1404 Nov  1 10:45 combined.mosaic.bulk.rda
-rw-rw-r-- 1 asa asa 1065072 Nov  1 10:45 data.bulk.rda
-rw-rw-r-- 1 asa asa 8716432 Nov  1 10:46 mutations.bulk.rda
-rw-rw-r-- 1 asa asa       0 Nov  1 10:46 onek.bulk.bed
-rw-rw-r-- 1 asa asa    3886 Nov  1 10:46 all.somatic.bulk.rda
-rw-rw-r-- 1 asa asa    3892 Nov  1 10:46 candidate.somatic.bulk.rda
-rw-rw-r-- 1 asa asa    2853 Nov  1 10:46 all.mosaic.bulk.rda
-rw-rw-r-- 1 asa asa     229 Nov  1 10:46 summary.bulk.rda
Track final file: summary.bulk.rda
And scripts in:
output/110/power.bulk_job_scripts/
- contains shell scripts 1_1.sh to 1_99.sh
"""

# will have to run in 2 steps, like split, since output is both dynamic and normal.
# check for summary.bulk.rda

#print(sample_dict["110"]["1"])

# cannot expand with lambda in output - only track script1, or use dynamic output. 
# also, only track one file for vcf.info
# but, then whole pipe does not work - rule ppower does not find a way to generate all the outputs.
# and, cannot mix dynamic and non-dynamic outputs.

rule compare:
    input:
        sample_conf = OUTDIR + "{sample}/config.txt",
        bulk_conf = OUTDIR + config["bulk"] + "/config.txt",
        sample_chk = OUTDIR + "{sample}/progress/plink_{chrom}.chk",
        bulk_chk = OUTDIR + config["bulk"] + "/progress/plink_{chrom}.chk",
    output:
        compare_sum = OUTDIR + "{sample}/{chrom}/compare/summary.bulk.rda",
    shell:
        "lira compare -s {input.sample_conf} -b {input.bulk_conf} -m {wildcards.chrom} -o;"

# cannot expand on output, so cannot track all the scripts in power.bulk_job_scripts
        
# to track the files in each subdir plus the final output file (summary.bulk.rda)
# instead of dynamic input, use expand
#    
#    input: dynamic(OUTDIR + "{{sample}}/power.bulk_job_scripts/{{chrom}}_{region_id}.sh")
#    input: lambda wc: expand(OUTDIR + "{{sample}}/power.bulk_job_scripts/{{chrom}}_{region_id}.sh", region_id = sample_dict[str(wc.sample)][wc.chrom])
#        
# make one check rule that takes the output from compare and checks that all the scripts are in place.

# make on rule that checks if all the power scripts exists and prints the bash commands to a file.

rule compare_chk:
    input: OUTDIR + "{sample}/{chrom}/compare/summary.bulk.rda"
    output: script = OUTDIR + "{sample}/power.bulk_job_scripts/run_all_ppower_scripts_{chrom}.sh"
    run:
        fh = open(output.script, "w")
        region_ids = sample_dict[str(wildcards.sample)][wildcards.chrom]
        for idx in region_ids:
            path = os.path.join(OUTDIR, wildcards.sample, "power.bulk_job_scripts", "%s_%s.sh"%(wildcards.chrom,idx))
            if os.path.exists(path):
                fh.write("bash "+ path + "\n")
            else:
                raise Exception("Compare chk: Cannot find script %s"%path)
        fh.close()
        
rule compare_onesample:
    input: expand(OUTDIR + "{sample}/power.bulk_job_scripts/run_all_ppower_scripts_{chrom}.sh", chrom = CHROMOSOMES, sample = ["{sample}"])           
    output:
        OUTDIR + "{sample}/progress/all_chrom_compare.chk"
    shell:
        "touch {output}"

        
"""
# ppower - runs all the scripts in output/110/power.bulk_job_scripts/
- creates in folder output/110/1/jobs/1
-rw-rw-r-- 1 asa asa   63 Nov  1 10:54 powers.bulk.rda
-rw-rw-r-- 1 asa asa    0 Nov  1 10:54 powers.one.bulk.bed
-rw-rw-r-- 1 asa asa    0 Nov  1 10:54 powers.two.bulk.bed
# and creates a chkfile:
output/110/progress/.power_bulk_1_85
# varcall
"""

rule ppower:
    input:
        compare_done = OUTDIR + "{sample}/{chrom}/compare/summary.bulk.rda",        
        script = OUTDIR + "{sample}/power.bulk_job_scripts/run_all_ppower_scripts_{chrom}.sh"
    output:
        chkfile = OUTDIR + "{sample}/progress/ppower_{chrom}.chk"
    shell:
        "bash {input.script};"
        "touch {output.chkfile}"    

# check that all the jobs ran correctly
# check for the powers file:
# OUTDIR + "{sample}/{chrom}/jobs/{region_id}/powers.two.bulk.bed"
# also for the checkfile:
#        chkfile = OUTDIR + "{sample}/progress/.power_bulk_{chrom}_{region_id}",
# write to output file all the finished jobs.

rule ppower_chk:
    input: OUTDIR + "{sample}/progress/ppower_{chrom}.chk"
    output: chkfile = OUTDIR + "{sample}/progress/ppower_{chrom}.chk2"
    run:
        fh = open(output.chkfile, "w")
        region_ids = sample_dict[str(wildcards.sample)][wildcards.chrom]
        for idx in region_ids:
            path1 = os.path.join(OUTDIR, wildcards.sample, "progress", ".power_bulk_%s_%s"%(wildcards.chrom,idx))
            path2 = os.path.join(OUTDIR, wildcards.sample, wildcards.chrom, "jobs", idx, "powers.two.bulk.bed")
            if os.path.exists(path1) and os.path.exists(path2):
                fh.write(path1 + "\n")
                fh.write(path2 + "\n")                
            else:
                raise Exception("Ppower chk failed: Cannot find chkfile %s or ppowes file %s"%(path1,path2))
        fh.close()

# run plink with all chroms for one sample
rule ppower_onesample:
    input: expand(OUTDIR + "{sample}/progress/ppower_{chrom}.chk2", chrom = CHROMOSOMES, sample = ["{sample}"])
    output: OUTDIR + "{sample}/progress/all_chrom_ppower.chk"
    shell:
        "touch {output}"

        
        
# final varcall
# Creates subfolder:
# dir <- paste("varcall_",bulk.config$name,sep="") -> sample/varcall_bulk/
# final file in script is: powered.regions.bed
# take as input all the chrom ppower outputs
# plus the chkfiles for each chrom.

rule varcall:
    input:
        sample_conf = OUTDIR + "{sample}/config.txt",
	bulk_conf = OUTDIR + config["bulk"] + "/config.txt",        
        sample_chkfile = OUTDIR + "{sample}/progress/all_chrom_ppower.chk",
    output:
        final_call = OUTDIR + "{sample}/varcall_" + config["bulk"] + "/powered.regions.bed"
    shell:
        "lira varcall -s {input.sample_conf} -b {input.bulk_conf}"


        
