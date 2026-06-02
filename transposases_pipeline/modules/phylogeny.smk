
rule fix_duplicate_headers:
    input:
        os.path.join(config['TASK'], "concat-seqs/{GENUS}.faa")
    output:
        os.path.join(config['TASK'], "results/duplicates-fixed/{GENUS}.faa")
    params:
        opt = "--by-name --force"
    singularity:
        "https://depot.galaxyproject.org/singularity/seqkit:2.7.0--h9ee0642_0"
    threads: 1
    shell:
        """
        seqkit rename {params.opt} {input} > {output}
        """

rule mafft_align_tree:
    input:
        os.path.join(config['TASK'], "results/duplicates-fixed/{GENUS}.faa")
    output:
        os.path.join(config['TASK'], "mafft/transposase-hits-from-{GENUS}hmm.mafft")
    params:
        tmp=config['TMP']
    threads: 8
    singularity:
        "https://depot.galaxyproject.org/singularity/mafft:7.525--h031d066_0"
    resources:
        mem_mb = 50000
    shell:
        """
        export MAFFT_TMPDIR={params.tmp}
        mkdir -p {params.tmp}
        mafft --auto --thread {threads} {input} > {output}
        """

rule clipkit_tree:
    input:
        os.path.join(config['TASK'], "mafft/transposase-hits-from-{GENUS}hmm.mafft")
    output:
        os.path.join(config['TASK'], "clipkit-transposase-hits/{GENUS}/hmm.clipkit")
    singularity:
        "https://depot.galaxyproject.org/singularity/clipkit:2.4.1--pyhdfd78af_0"
    shell:
        """
        clipkit {input} -o {output}
        """

rule iqtree2:
    input:
        os.path.join(config['TASK'], "clipkit-transposase-hits/{GENUS}/hmm.clipkit")
    output:
        os.path.join(config['TASK'], "iqtree2/{GENUS}/aft2_transposases.treefile")
    params:
        aln_dir = lambda wc, input: os.path.dirname(input[0]),
        out_dir = lambda wc, output: os.path.splitext(output[0])[0],
        omp     = lambda wildcards, threads: threads,
        opt     = "-redo -bb 1000 -m TEST --verbose --seed 1 --mem 100G"
    singularity:
        "https://depot.galaxyproject.org/singularity/iqtree:2.4.0--h503566f_0"
    threads: 48
    resources:
        partition = "standard,long",
        mem_mb    = 100000,
        time      = "3-00:00:00"
    shell:
        """
        export OMP_NUM_THREADS={params.omp};

        iqtree2 \
            -s {params.aln_dir} \
            -T AUTO \
            -ntmax {threads} \
            --threads-max {threads} \
            --prefix {params.out_dir} \
            {params.opt}

        touch {output}
        """