
rule make_coverage_pickle:
    input: counts_file = config["simulation"]["counts_file"]
    output: pickle = config["simulation"]["lc_pickle"]           
    params:
        limit = config["simulation"]["locus_count_limit"],
        path = config["simulation"]["path"],
        fixedCoverage = config["simulation"]["fixed"]
    run:
        import sys, pickle
        sys.path.insert(0, params.path)

        if params.fixedCoverage == True:
            from ReadsDb import createFixedLocusCounts                    
            locusCounts = createFixedLocusCounts(30)
        else:
            from ReadsDb import createLocusCounts        
            locusCounts = createLocusCounts(input.counts_file, params.limit)
        with open(output.pickle, 'wb') as f:
            pickle.dump(locusCounts, f, protocol=0)

# Skip simulation of bulk and instead extract regions from "homf"? May be needed if we also run lira
rule simulate_bulk:
    input:
        homf = config["simulation"]["homf"],
        hetf = config["simulation"]["hetf"],
        sites = config["simulation"]["sites_file"]           
    output: "data/common/homo_bulk.bam"
    params:
        path = config["simulation"]["path"],
        outprefix = "data/common/homo"
    log: "data/common/simulate_bulk.log"
    run:
        import sys, pickle
        sys.path.insert(0, params.path)
	from ReadsDb import ReadsDb
        readsDB = ReadsDb(input.hetf, input.homf, input.sites)
        print("outprefix")
        print(params.outprefix)
        readsDB.writeBulkToFile(params.outprefix)


        
# OBS! now hard-coded that ncell can only be 10 cells, otherwise, need to have pop-structure in config.
rule simulate_cells:
    input:
        homf = config["simulation"]["homf"],
        hetf = config["simulation"]["hetf"],
        pickle = config["simulation"]["lc_pickle"],
        sites = config["simulation"]["sites_file"]
    output:
        bamfiles = expand("data/sim_snv{{f_SNV}}_eal{{f_EAL}}_ado{{f_ADO}}/bams/sim_cell{nn}.bam",nn=range(config["simulation"]["ncell"])),
        json = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/bams/sim_genVals.json",
        chkfile = "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/bams/sim.chk"
    params:
        path = config["simulation"]["path"],
        tree = config["simulation"]["phyl_tree"],
        freq = config["simulation"]["use_freq"]
    log: "data/sim_snv{f_SNV}_eal{f_EAL}_ado{f_ADO}/bams/simulate_cells.log"
    run:
        import sys, pickle
	sys.path.insert(0, params.path)
        from ReadsDb import ReadsDb

        with open(input.pickle, 'rb') as f:
            locusCounts = pickle.load(f)
        readsDB = ReadsDb(input.hetf, input.homf, input.sites)

        T = params.tree
        print("tree:")
        print(T)
        #from ReadsDb import flatten, unnest
        #print("flatten")
        #print(flatten(T))
        #print("unnest")
        #print(unnest(T))
        
        outprefix = output.bamfiles[0]
        outprefix = outprefix.replace("_cell0.bam","")

        readsDB.simulateAndWriteToFile(T, float(wildcards.f_SNV),float(wildcards.f_EAL), float(wildcards.f_ADO), locusCounts, outprefix, freqs = params.freq)

        touch(output.chkfile)

        

import os        
def touch(path):        
    basedir = os.path.dirname(path)
    if not os.path.exists(basedir):
        os.makedirs(basedir)
    with open(path, 'a'):
        os.utime(path, None)


