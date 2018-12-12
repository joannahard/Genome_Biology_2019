rule mpileup:
    input:
        bam = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/bams/sim_cell{nn}_sorted.bam",
        bai = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/bams/sim_cell{nn}_sorted.bai",
        ref = config["settings"]["reference"]
    output: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/mpileup/sim_cell{nn}.bam.mpileup"
    shell:
        'samtools mpileup -C50 -Osf {input.ref} {input.bam} > {output}'


rule sccaller_header:
    input: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/mpileup/sim_cell{nn}.bam.mpileup"
    output: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/output/sim_cell{nn}.hsnp.bed"
    params:
        prefix = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/output/sim_cell{nn}",
        conda = config["sccaller"]["conda"],
        dbsnp = config["sccaller"]["dbsnp"],
    shell:
        'source activate {params.conda};'
        'sccaller -a hsnp -i {input}  --snpin {params.dbsnp}  --snp dbsnp -o {params.prefix}'

rule sccaller_varcall:
    input:
        hsnp = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/output/sim_cell{nn}.hsnp.bed",
        mpile = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/mpileup/sim_cell{nn}.bam.mpileup"
    output: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/output/sim_cell{nn}.varcall.bed"
    params:
        prefix = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/output/sim_cell{nn}",
        conda = config["sccaller"]["conda"]
    shell:
        'source activate {params.conda};'
	'sccaller -a varcall -i {input.mpile} -s {input.hsnp} -o {params.prefix}'

# OBS! Input and output file is the same. Just modifies the input file with filtering cutoff.
rule sccaller_cutoff:
    input: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/output/sim_cell{nn}.varcall.bed"
    output: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/output/sim_cell{nn}.varcall.cutoff.bed"
    params:
        conda = config["sccaller"]["conda"],
    shell:
        'cp {input} {output};'
        'source activate {params.conda};'
        'sccaller -a cutoff -i {output}'

# instead of their awk command, use python to filter the file.        
# awk '$12!=0' outputheader.chr1.varcall.bed | awk '$11/$12<0.eta && $7/($6+$7)>1/8 && $8=="PASS"' > outputheader.chr1.sccaller.bed
# eta from:
##Cutoffs for L_artifact / L_H1 (if smaller than the following eta, reject H_artifact, accept H1): at alpha=0.01,eta=4.03858670759e-194; alpha=0.05,eta=3.18205160072e-44
# make filtering for both alpha values.


rule sccaller_filter:
    input: vcf = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/output/sim_cell{nn}.varcall.cutoff.bed"
    output:
        sc01 = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/output/sim_cell{nn}.sincell_filtered0.01.bed",
        sc05 = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/output/sim_cell{nn}.sincell_filtered0.05.bed",        
    run:
        # first read all lines and last line contains the eta cutoffs
        data = []
        for line in open(input.vcf):
            if not line.startswith("#"):
                data.append(line)
            else:
                l = line.replace(";",",").split(",")
                eta01 = float(l[3].split("=")[1])
                eta05 = float(l[5].split("=")[1])
        fh01 = open(output.sc01,'w')
        fh05 = open(output.sc05,'w')        
        for line in data:
            l = line.split("\t")
            if l[7]=="PASS" and float(l[11]) > 0:
                if float(l[6])/(float(l[5]) + float(l[6])) > 1/8:
                    if float(l[10])/float(l[11]) < eta01:
                        fh01.write(line)
                    if float(l[10])/float(l[11]) < eta05:
                        fh05.write(line)

        fh01.close()
        fh05.close()

        # awk '$12!=0' outputheader.chr1.varcall.bed | awk '$11/$12<0.eta && $7/($6+$7)>1/8 && $8=="PASS"' > outputheader.chr1.sccaller.bed
            
                
rule sccaller_allcells:
    input:
        eta01 = expand("data/sim_snv{{f_SNV}}_eal{{f_EAL}}_ado{{f_ADO}}/sccaller/output/sim_cell{nn}.sincell_filtered0.01.bed",nn=range(config["simulation"]["ncell"])),
        eta05 = expand("data/sim_snv{{f_SNV}}_eal{{f_EAL}}_ado{{f_ADO}}/sccaller/output/sim_cell{nn}.sincell_filtered0.05.bed",nn=range(config["simulation"]["ncell"]))        
    output: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/all_cells.chk"
    shell:
        'touch {output}'

rule sccaller_merge:
    input: expand("data/sim_snv{{f_SNV}}_eal{{f_EAL}}_ado{{f_ADO}}/sccaller/output/sim_cell{nn}.sincell_filtered{{eta}}.bed",nn=range(config["simulation"]["ncell"]))
    output: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/sccaller/all_cells_stats.{eta}.csv"
    params: script = config["sccaller"]["merge_script"]
    shell:
        'python {params.script} -o {output} -i {input}'
        
