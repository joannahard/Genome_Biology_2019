

rule merge_stats:
    input: expand( "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/simulation_stats.csv", f_SNV = config["simulation"]["f_SNV"],f_EAL = config["simulation"]["f_EAL"],f_ADO = config["simulation"]["f_ADO"])
    output:
        cb = "data/stats/stats_all_sim_conbase.csv",
        mv = "data/stats/stats_all_sim_monovar.csv"        
    params:
        path = config["parsing"]["path_sum"],
        outprefix = "data/stats/stats_all_sim"
    shell:
        "python3 {params.path} -o {params.outprefix} -i {input}"


rule parse_data:
    input:
        monovar = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/monovar/monovar_output.vcf",
        conbase = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/conbase/conbase_results.tsv",
        sim = 	"data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/bams/sim_genVals.json"
    output: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/simulation_stats.csv"
    params:
        dp = config["parsing"]["monovar_dp"],
        gq = config["parsing"]["monovar_gq"],
        nmut = config["parsing"]["monovar_nmut"],
        path = config["parsing"]["path"]
    shell:
        "python3 {params.path} -m {input.monovar} -c {input.conbase} -s {input.sim} -o {output} --monovar_gq {params.gq} --monovar_dp {params.dp} --monovar_nmut {params.nmut}"
        
