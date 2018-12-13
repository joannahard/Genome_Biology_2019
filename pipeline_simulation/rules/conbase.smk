######################
# conbase

rule make_conbase_bamlist:
    input:
        bamfiles = expand("data/sim_snv{{f_SNV}}_eal{{f_EAL}}_ado{{f_ADO}}/bams/sim_cell{nn}_sorted.bam",nn=range(config["simulation"]["ncell"])),
        baifiles = expand("data/sim_snv{{f_SNV}}_eal{{f_EAL}}_ado{{f_ADO}}/bams/sim_cell{nn}_sorted.bai",nn=range(config["simulation"]["ncell"])),        
        bulk_bam = "data/common/homo_bulk.bam",
        bulk_bai = "data/common/homo_bulk.bai"        
    output:
        bam_list = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/conbase/bamfiles.txt"
    run:
        fh = open(output.bam_list,"w")
        fh.write("NAME\tPATH\n")
        fh.write("BULK\t" + input.bulk_bam + "\n")
        for bfile in input.bamfiles:
            cname = os.path.split(bfile)[1]
            cname = cname.replace(".bam","")
            fh.write(cname + "\t" + bfile + "\n")
        fh.close()


# OBS! need devel version of conbase to run without absolute paths and getting output in "../results"

rule conbase_stats:
    input:
        bamfiles = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/conbase/bamfiles.txt",
        reference = config["settings"]["reference"],
        sites_file = config["conbase"]["sites_file"]
    output: 
        json = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/conbase/conbase_results.json"
    params:
        path = config["conbase"]["path"],

        outprefix = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/conbase/conbase_results"
    threads: config["conbase"]["threads"]
    shell:
        "python3 {params.path} --stats {input.sites_file} {input.bamfiles} {input.reference} {threads} {params.outprefix}"

        
# python3 bin/Main.py stats <snp path> <bam path> <reference path> <number of threads> <output name>


rule conbase_analyse:
    input:
        json = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/conbase/conbase_results.json"
    output:
        pdf = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/conbase/conbase_results.pdf",
        tsv = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/conbase/conbase_results.tsv",        
        html = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/conbase/conbase_results.html"
    params:
        path = config["conbase"]["path"],
	outprefix = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/conbase/conbase_results"
    shell:
        "python3 {params.path} --analyze {input.json} {params.outprefix}"

#python3 bin/Main.py analyze <json path> <output name>

