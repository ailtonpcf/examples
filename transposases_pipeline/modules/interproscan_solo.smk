include: "/scratch/proj/00-default/smk-functions/resources.py"

rule interproscan:
    input:
        os.path.join(config['TASK'], "transposases-from-hmm/{HPF}__{GENUS}.{EXT}")
    output: 
        tsv=os.path.join(config['TASK'], "interproscan/{HPF}/{HPF}-{GENUS}.{EXT}.tsv"),
        xml=os.path.join(config['TASK'], "interproscan/{HPF}/{HPF}-{GENUS}.{EXT}.xml"),
    params:
        out_dir= lambda wc, output: os.path.dirname(output.tsv),
        omp=lambda wildcards, threads: threads,
        int_dir="/home/groups/Fungal/marion/softwares/my_interproscan/interproscan-5.68-100.0",
        tmp=config['TMP_DIR'],
        opt="--formats TSV,XML"
    retries: 2
    threads: lambda wildcards, attempt: adjust_resources_medium_medium(wildcards, attempt)['threads']
    resources:
        partition = lambda wildcards, attempt: adjust_resources_medium_medium(wildcards, attempt)['partition'],
        mem_mb    = lambda wildcards, attempt: adjust_resources_medium_medium(wildcards, attempt)['mem_mb'],
        time      = lambda wildcards, attempt: adjust_resources_medium_medium(wildcards, attempt)['time'],
    shell:
        """
        export OMP_NUM_THREADS={params.omp};

        mkdir -p {params.out_dir}
        abs_in=$(readlink -f {input})
        abs_out_dir=$(readlink -f {params.out_dir})

        cd {params.int_dir}

        conda run --no-capture-output --prefix /home/groups/Fungal/conda_env/funannotate \
            ./interproscan.sh \
                -i "$abs_in" \
                -d "$abs_out_dir" \
                -goterms \
                -dp \
                -cpu {threads} \
                --tempdir {params.tmp} \
                {params.opt}
        """