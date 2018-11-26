######################
# conbase

rule make_monovar_bamlist:
    input:
        bamfiles = expand("data/sim_snv{{f_SNV}}_eal{{f_EAL}}_ado{{f_ADO}}/bams/sim_cell{nn}_sorted.bam",nn=range(config["simulation"]["ncell"])),
        baifiles = expand("data/sim_snv{{f_SNV}}_eal{{f_EAL}}_ado{{f_ADO}}/bams/sim_cell{nn}_sorted.bai",nn=range(config["simulation"]["ncell"]))
    output:
        bam_list = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/monovar/bamfiles.txt"
    run:
        fh = open(output.bam_list,"w")
        for bfile in input.bamfiles:
            fh.write(bfile + "\n")
        fh.close()


# monvar call:
# source activate monovar
# export PATH=$PATH:/media/box2/Experiments/Joanna/Snake_analys/j_frisen_1602/programs/monovar/src
# samtools mpileup -BQ0 -d10000 -f ref.fa -q 40 -b filenames.txt | monovar.py -p 0.002 -a 0.2 -t 0.05 -m 2 -f ref.fa -b filenames.txt -o output.vcf


rule monovar:
    input:
        bam_list = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/monovar/bamfiles.txt",
        reference = config["settings"]["reference"]
    output:
        vcf = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/monovar/monovar_output.vcf"
    params:
        mpileup_settings = config["monovar"]["mpileup_settings"],
        monovar_settings = config["monovar"]["monovar_settings"]
    threads: config["monovar"]["threads"]
    shell:
        "source activate monovar;"
        "export PATH=$PATH:/media/box2/Experiments/Joanna/Snake_analys/j_frisen_1602/programs/monovar/src;"
        "samtools mpileup {params.mpileup_settings} -f {input.reference} -b {input.bam_list} | monovar.py {params.monovar_settings} -m {threads} -f {input.reference} -b {input.bam_list} -o {output.vcf}"
