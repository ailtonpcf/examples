rule mafft_align_hmm:
    input:
        "/scratch/ref/sequences/aft2_domains_confirmed.faa"
    output:
        os.path.join(config['TASK'], "mafft/aft2_domains_confirmed.mafft")
    params:
        tmp=config['TMP']
    threads: 8
    singularity:
        "https://depot.galaxyproject.org/singularity/mafft:7.525--h031d066_0"
    shell:
        """
        export MAFFT_TMPDIR={params.tmp}
        mkdir -p {params.tmp}
        mafft --auto --thread {threads} {input} > {output}
        """

rule clipkit_hmm:
    input:
        os.path.join(config['TASK'], "mafft/aft2_domains_confirmed.mafft")
    output:
        os.path.join(config['TASK'], "clipkit/aft2_domains_confirmed.clipkit")
    singularity:
        "https://depot.galaxyproject.org/singularity/clipkit:2.4.1--pyhdfd78af_0"
    shell:
        """
        clipkit {input} -o {output}
        """

rule hmmbuild:
    input:
        os.path.join(config['TASK'], "clipkit/aft2_domains_confirmed.clipkit")
    output:
        os.path.join(config['TASK'], "transposase-hmm/aft2_domains_confirmed.hmm")
    threads: 8
    singularity:
        "https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_2"
    shell:
        """
        hmmbuild --cpu {threads} {output} {input} 
        """