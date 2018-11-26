rule index_bam:
    input: "{path}.bam"
    output: "{path}.bai"
    params:
        java = config["settings"]["javaopts"]
    shell:
        "picard {params.java} BuildBamIndex INPUT={input}"
            


rule sort_bam:
    input: "{path}.bam"
    output: "{path}_sorted.bam"
    params:
        java = config["settings"]["javaopts"]
    shell:
        "picard {params.java} SortSam SORT_ORDER=coordinate INPUT={input} OUTPUT={output}"
            

rule link_fasta:
    input: config["settings"]["reference_original"]
    output: config["settings"]["reference"]
    shell:
        "ln -s {input} {output}"


        
