
rule gatk_dict:
    input: "{ref}.fasta"
    output: "{ref}.dict"
    shell:
        "picard CreateSequenceDictionary R={input} O={output}"

"""        
# make a bedfile with regions 1kb upstream/downstream of gSNV and sSNV        
rule make_region_list:
    input:
        sites = config["simulation"]["sites_file"]
    output:
        bed = "data/common/interval_list_for_gatk.bed"
    params:
        dist = config["settings"]["interval_dist"]
    run:
        fout = open(output.bed, "w")
        for line in open(input.sites):
            if line.startswith("chrom"): continue
            l = line.strip().split("\t")
            start = int(l[1])-int(params.dist)
            if start < 0: start = 0
            end = int(l[2])+int(params.dist)
            fout.write("%s\t%d\t%d\n"%(l[0],start,end))
        fout.close()
"""        

# GenomeAnalysisTK -T HaplotypeCaller -R /media/box2/reference_assemblies/bundle/2.8/b37/from_pall/human_g1k_v37_decoy.fasta -I /media/box2/Experiments/Joanna/LiRA/Tcell_bams/bulk/bulk.bam --emitRefConfidence GVCF -o /media/box2/Experiments/Joanna/LiRA/Tcell_bams/bulk/vcf/bulk.raw.g.vcf --variant_index_type LINEAR --variant_index_parameter 128000
rule gatk_haplotypecaller:
    input:
        bamfile = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/bams/{sample}.bam",
        baifile = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/bams/{sample}.bai",        
        reference = config["settings"]["reference"],
        ref_dict = config["settings"]["ref_dict"],
    output:
        vcf = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/vcf/{sample}.g.vcf",
    params:
        conda = config["lira"]["gatk_conda"],
        settings = config["lira"]["gatk_settings"],
        java = config["settings"]["javaopts10"],
        use_threads = config["lira"]["gatk_threads_use"]
    threads: config["lira"]["gatk_threads"]
    log:
        gatk = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/vcf/logs/{sample}_haplotypecaller.log"
    shell:
        "source activate {params.conda};"
        "GenomeAnalysisTK {params.java} -T HaplotypeCaller -nct {params.use_threads} -R {input.reference} -I {input.bamfile} -o {output.vcf} {params.settings} > {log.gatk} 2>&1;"
        
#     java -jar GenomeAnalysisTK.jar \
#   -T GenotypeGVCFs \
#   -R reference.fasta \
#   --variant 110.raw.g.vcf  \
#   --variant 111.raw.g.vcf  \
#   -o output.vcf


rule gatk_joint_calling:
    input:
        vcfs = expand("data/sim_snv{{f_SNV}}_eal{{f_EAL}}_ado{{f_ADO}}/lira/vcf/cell{nn}.g.vcf",nn=range(config["simulation"]["ncell"])),
        bulk_vcf = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/vcf/bulk.g.vcf",
        reference = config["settings"]["reference"],
        ref_dict = config["settings"]["ref_dict"],
    output:
        joint_vcf = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/vcf/all.vcf"
    params:
        conda = config["lira"]["gatk_conda"],
        java = config["settings"]["javaopts"]
    log:
        gatk = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/vcf/logs/joint_call.log"

    run:
        cmd = "source activate "+ params.conda + "; GenomeAnalysisTK "+ params.java + " -T GenotypeGVCFs -R " + input.reference + " -o " + output.joint_vcf 
        for vfile in input.vcfs:
            cmd = cmd + " --variant "+vfile
        cmd = cmd + " --variant " + input.bulk_vcf
        cmd = cmd + " >  " + log.gatk + " 2>&1;"
        shell(cmd)
        
        
"""
    shell:
        'soure activate {params.conda};'
        'command="GenomeAnalysisTK {params.java} -T GenotypeGVCFs -R {input.reference} -o {output.joint_vcf} ";'
        'for f in {input.vcfs}; do command="$command --variant $f "; done;'
        '"$command" > {log.gatk} 2>&1;'
"""

# zip and index
rule zip_vcf:
    input: "{file}.vcf"
    output: "{file}.vcf.gz"
    shell:
        "bgzip -c {input} > {output}"

rule index_vcf:
    input: "{file}.vcf.gz"
    output: "{file}.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"


