configfile: "/scratch/proj/02-compost-microbes/src/116-hpf-transposases/smk.config"
workdir: config['workdir']
import pandas as pd

HPF_ID, = glob_wildcards(os.path.join(config['TASK'], "transposases-hpf/{HPF_ID}.faa"))
ASP_ID, = glob_wildcards(os.path.join(config['TASK'], "transposases-asp/{ASP_ID}.faa"))

# HPF_ID=HPF_ID[1]
# ASP_ID=ASP_ID[1]

rule all:
    input:
         os.path.join(config['TASK'], "needle-stats-merged/aspTrans_vs_otherHPFTrans.tsv")
    default_target: True

rule emboss_needle:
    input:
        seqA = os.path.join(config['TASK'], "genes/transposase.faa"),
        seqB = os.path.join(config['TASK'], "transposases-hpf/{HPF_ID}.faa"),
    output: 
        os.path.join(config['TASK'], "needle-stats/{HPF_ID}.txt")
    params:
        opt="-gapopen 10 -gapextend 0.5"
    threads: 1
    singularity:
        "https://depot.galaxyproject.org/singularity/emboss:6.6.0--hf7d6862_3"
    shell:
        """
        needle \
            -asequence {input.seqA} \
            -bsequence {input.seqB} \
            -outfile {output} \
            {params.opt}
        """

use rule emboss_needle as emboss_needle_asp with:
    input:
        seqA = os.path.join(config['TASK'], "genes/transposase.faa"),
        seqB = os.path.join(config['TASK'], "transposases-asp/{ASP_ID2}.faa"),
    output: 
        os.path.join(config['TASK'], "needle-stats/{ASP_ID2}.txt")

rule summarise_needle_metrics:
    input:
        expand(os.path.join(config['TASK'], "needle-stats/{HPF_ID}.txt"), HPF_ID=HPF_ID),
        expand(os.path.join(config['TASK'], "needle-stats/{ASP_ID2}.txt"), ASP_ID2=ASP_ID)
    output: 
        os.path.join(config['TASK'], "needle-stats-merged/aspTrans_vs_otherHPFTrans.tsv")
    params:
        scrip_dir="/scratch/proj/02-compost-microbes/src/116-hpf-transposases/py",
        input_dir = lambda wc, input: os.path.dirname(input[0])
    threads: 1
    conda: "python3.8"
    shell:
        """
        python {params.scrip_dir}/extract_needle_metrics.py {params.input_dir} {output}
        """