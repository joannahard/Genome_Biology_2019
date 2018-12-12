# for lira, run the other snakemake pipeline one sample at a time?


# need to have files named sample.bam - same name as in RG
# e.g. cells file need to be named cell0.bam etc.
# and bulk file bulk.bam
# OBS! cannot symlink the bai files, ambigous rule to index.bam - need to reindex them.

# check that also bai exists!
# Lira uses normalizePath when reading config!!

rule symlink_bamfiles: 
    input:
        bam = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/bams/sim_cell{nn}_sorted.bam",
        bai = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/bams/sim_cell{nn}_sorted.bai"        
    output: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/bams/cell{nn}.bam"
    shell:
        'ln -s "$(realpath {input.bam})" {output}'        
        
# ROOTDIR=os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

# OBS! Need to have the full bamfile to run lira - to be able to do phasing of gSNVs.
# but extract only from regions of interest.
# also need to change readgroup to "bulk"

rule readgroup_lira_bulk_bam:
    input:
        bam = "data/common/lira_bulk_wrong_rg.bam"
    output:
        bam = "data/common/lira_bulk.bam"
    params:
        java = config["settings"]["javaopts"],
        read_groups = "CREATE_INDEX=true RGID=bulk RGLB=bulk RGPL=ILLUMINA RGSM=bulk RGCN=\"NA\" RGPU=\"NA\""
    shell:
        "picard {params.java} AddOrReplaceReadGroups {params.read_groups} INPUT={input.bam} OUTPUT={output.bam}"

rule get_lira_bulk_regions:
    input:
        bam = config["simulation"]["homf"],
    output:
        bam = "data/common/lira_bulk_wrong_rg.bam",
    params:
        bed = "data/common/interval_list_for_gatk.bed"
    shell:
        'samtools view -b -L {params.bed} {input.bam} > {output.bam}' 

        

rule symlink_bulk_bams:
    input:
        bam = "data/common/lira_bulk.bam",
        bai = "data/common/lira_bulk.bai"        
    output: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/bams/bulk.bam",
    shell:
        'ln -s "$(realpath {input.bam})" {output}'        


rule symlink_lira_ref:
    input: config["settings"]["reference_original"]
    output: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/lira_output/ref/human_g1k_v37_decoy.fasta"
    shell:
        'ln -s {input} {output}'
            

rule lira_config:
    input:
        bamfiles = expand("data/sim_snv{{f_SNV}}_eal{{f_EAL}}_ado{{f_ADO}}/lira/bams/cell{nn}.bam",nn=range(config["simulation"]["ncell"])),
        baifiles = expand("data/sim_snv{{f_SNV}}_eal{{f_EAL}}_ado{{f_ADO}}/lira/bams/cell{nn}.bai",nn=range(config["simulation"]["ncell"])),
        bulk_bam = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/bams/bulk.bam",
        bulk_bai = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/bams/bulk.bai",
        vcf_file = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/vcf/all.vcf.gz",
        vcf_index = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/vcf/all.vcf.gz.tbi",
        ref = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/lira_output/ref/human_g1k_v37_decoy.fasta",
        config = config["lira"]["config_template"]
    output:
        config = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/config_pipeline.yaml"
    params:
        outdir = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/"
    run:
        fout = open(output.config, "w")
        cwd = os.getcwd()
        for line in open(input.config):
            # fix absolute paths to the bamfiles and output dir.
            if line.startswith("bamdir"):
                fout.write("bamdir: " + cwd + "/" + params.outdir + "bams/\n")
                continue
            if line.startswith("resdir_fullpath"):
                fout.write("resdir_fullpath: " + cwd + "/" + params.outdir + "lira_output/\n")
                continue
            if line.startswith("vcf"):
                fout.write("vcf: " + cwd + "/" + input.vcf_file + "\n")
                continue
            if line.startswith("chromosomes"):
                fout.write("chromosomes: [" + ",".join([str(x) for x in config["lira"]["chromosomes"]]) + "]\n")
                continue
            fout.write(line)            
            if line.startswith("samples"):
                break
        # write all sample names.                
        for bamfile in input.bamfiles:
            cname = os.path.split(bamfile)[1]
            cname = cname.replace(".bam","")
            fout.write(" - " + cname + "\n")

        fout.write("bulk: bulk")    
        fout.close()            


#export PATH=/media/box2/Experiments/Joanna/bin/miniconda3/envs/LiRA/bin:/media/box2/Experiments/Joanna/bin/miniconda3/envs/snake2/bin:$PATH
#export PATH=$PATH:/media/box2/Experiments/Joanna/LiRA/LiRA
#export LIRA_DIR=/media/box2/Experiments/Joanna/LiRA/LiRA
# snakemake -p -j 22 lira_output/110/split_all_chrom.chk         

rule lira_split:
    input: 
        config = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/config_pipeline.yaml",
        bam = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/bams/{sample}.bam",
        vcf = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/vcf/all.vcf.gz.tbi"        
    output:
        chkfile = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/lira_output/{sample}/split_all_chrom.chk"
    params:
        cwd = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/",
        snakefile = config["lira"]["snakefile1"],
        snake_conda_path = config["lira"]["snake_conda_path"],
        lira_conda_path = config["lira"]["lira_conda_path"],
        lira_dir = config["lira"]["lira_dir"]
    threads:
        config["lira"]["split_threads"]
    log:
        logfile = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/lira_output/logs/split_{sample}.log"
    run:
        chkfile = output.chkfile.replace(params.cwd,"")
        logfile = log.logfile.replace(params.cwd,"")
        cmd = "cd " + params.cwd + ";"
        cmd = cmd + "export PATH=" + params.lira_conda_path + ":" + params.snake_conda_path + ":$PATH:" + params.lira_dir + ";"
        cmd = cmd + "export LIRA_DIR=" + params.lira_dir + ";"
        cmd = cmd + "snakemake -p --nolock -j " + str(threads) + " --snakefile " + params.snakefile + " " +  chkfile + " > "  + logfile + " 2>&1;"
        print(cmd)
        shell(cmd)
        

# obs! have relative paths to chkfile & logfile, need to be relative to lira_output


# make an extra rule that makes sure that plink is run for the bulk sample before running varcall for each cell - otherwise conflicts when multiple varcall jobs are running plink for the bulk.
rule lira_bulk_plink:
    input:
        config = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/config_pipeline.yaml",
        chkfile_bulk = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/lira_output/bulk/split_all_chrom.chk"
    output:
        plink_file = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/lira_output/bulk/progress/plink_{chrom}.chk"
    params:
        cwd = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/",
        snakefile = config["lira"]["snakefile2"],
        snake_conda_path = config["lira"]["snake_conda_path"],
        lira_conda_path = config["lira"]["lira_conda_path"],
        lira_dir = config["lira"]["lira_dir"]
    log:
        logfile = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/lira_output/logs/plink_bulk{chrom}.log"
    run:
       	pfile = output.plink_file.replace(params.cwd,"")
	logfile = log.logfile.replace(params.cwd,"")
	cmd = "cd " + params.cwd + ";"
	cmd = cmd + "export PATH=" + params.lira_conda_path + ":" + params.snake_conda_path + ":$PATH:" + params.lira_dir + ";"
	cmd = cmd + "export LIRA_DIR=" + params.lira_dir + ";"
        cmd = cmd + "snakemake -p -j " + str(threads) + " --snakefile " + params.snakefile + " " +  pfile + " > "  + logfile + " 2>&1;"
        print(cmd)
	shell(cmd)

# and a summary rule for all the chroms, creates file lira_output/bulk/progress/all_chrom_plink.chk that is input to lira_varcall rule.
rule lira_bulk_plink_chk:
    input: expand("data/sim_snv{{f_SNV}}_eal{{f_EAL}}_ado{{f_ADO}}/lira/lira_output/bulk/progress/plink_{chrom}.chk", chrom=config['lira']["chromosomes"])
    output: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/lira_output/bulk/progress/all_chrom_plink.chk"
    shell:
        'touch {output}'
        
#export PATH=/media/box2/Experiments/Joanna/bin/miniconda3/envs/LiRA/bin:/media/box2/Experiments/Joanna/bin/miniconda3/envs/snake2/bin:$PATH
#export PATH=$PATH:/media/box2/Experiments/Joanna/LiRA/LiRA
#export LIRA_DIR=/media/box2/Experiments/Joanna/LiRA/LiRA
#snakemake -p  -j 7 lira_output/111/varcall_bulk/powered.regions.bed --snakefile snakefile_non_dynamic.smk 

rule lira_varcall:
    input:
        config = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/config_pipeline.yaml",
	chkfile = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/lira_output/cell{nn}/split_all_chrom.chk",
	chkfile_bulk = ancient("data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/lira_output/bulk/progress/all_chrom_plink.chk")
    output:
        pfile = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/lira_output/cell{nn}/varcall_bulk/powered.regions.bed"
    params:
        cwd = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/",
        snakefile = config["lira"]["snakefile2"],
        snake_conda_path = config["lira"]["snake_conda_path"],
        lira_conda_path = config["lira"]["lira_conda_path"],
        lira_dir = config["lira"]["lira_dir"]
    threads:
        config["lira"]["varcall_threads"]
    log:
        logfile = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/lira_output/logs/varcall_cell{nn}.log"
    run:
        pfile = output.pfile.replace(params.cwd,"")
        logfile = log.logfile.replace(params.cwd,"")
        cmd = "cd " + params.cwd + ";"
        cmd = cmd + "export PATH=" + params.lira_conda_path + ":" + params.snake_conda_path + ":$PATH:" + params.lira_dir + ";"
        cmd = cmd + "export LIRA_DIR=" + params.lira_dir + ";"
        cmd = cmd + "snakemake -p --nolock -j " + str(threads) + " --snakefile " + params.snakefile + " " +  pfile + " > "  + logfile + " 2>&1;"
        print(cmd)
        shell(cmd)


# make one checkfile when all cells from one simulation setting has finished.        
rule lira_one_sim:
    input:
        pfiles = expand("data/sim_snv{{f_SNV}}_eal{{f_EAL}}_ado{{f_ADO}}/lira/lira_output/cell{nn}/varcall_bulk/powered.regions.bed",nn=range(config["simulation"]["ncell"]))
    output:
        chkfile = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/lira_output/all_cells.chk"
    shell:
        "touch {output.chkfile}"


rule lira_merge:
    input: expand("data/sim_snv{{f_SNV}}_eal{{f_EAL}}_ado{{f_ADO}}/lira/lira_output/cell{nn}/varcall_bulk/out.vcf.gz",nn=range(config["simulation"]["ncell"]))
    output: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/lira/all_cell_stats.csv"
    params: script = config["lira"]["merge_script"]
    shell:
        'python {params.script} -o {output} -i {input}'
