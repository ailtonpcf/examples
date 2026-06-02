include: "/scratch/proj/00-default/smk-functions/resources.py"

rule funannotate_clean:
    input:
        os.path.join(config['TASK'], "genomes/{ID}.fasta")
    output: 
       os.path.join(config['TASK'], "clean/{ID}_clean.fa")
    params:
        gmes_path="/scratch/apps/gmes_linux_64_4",
    retries: 2
    threads: 1
    resources:
        partition = lambda wildcards, attempt: adjust_resources_medium(wildcards, attempt)['partition'],
        mem_mb    = lambda wildcards, attempt: adjust_resources_medium(wildcards, attempt)['mem_mb'],
        time      = lambda wildcards, attempt: adjust_resources_medium(wildcards, attempt)['time'],
    shell: 
        """
        conda run --no-capture-output --prefix /home/groups/Fungal/conda_env/funannotate \
          funannotate clean -i {input} -o {output}
        """

rule funannotate_sort:
    input:
        os.path.join(config['TASK'], "clean/{ID}_clean.fa")
    output: 
       os.path.join(config['TASK'], "sort/{ID}_sort.fa")
    params:
        gmes_path="/scratch/apps/gmes_linux_64_4",
    threads: 1
    retries: 2
    resources:
        partition = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['partition'],
        mem_mb    = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['mem_mb'],
        time      = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['time'],
    shell: 
        """
        conda run --no-capture-output --prefix /home/groups/Fungal/conda_env/funannotate \
            funannotate sort -i {input} -o {output} --minlen 1

        #the "minlen" is to avoid a weird bug.
        """

rule RepeatMasker:
    input:
        os.path.join(config['TASK'], "sort/{ID}_sort.fa")
    output: 
       os.path.join(config['TASK'], "masking/{ID}_sort.fa.masked")
    params:
        out_dir=os.path.join(config['TASK'], "masking"),
        gmes_path="/scratch/apps/gmes_linux_64_4",
        lib="/home/groups/Fungal/marion/softwares/RepBase28.03.fasta/fngrep.ref",
        omp=0
    singularity:
        "https://depot.galaxyproject.org/singularity/repeatmasker:4.1.5--pl5321hdfd78af_1"
    retries: 2
    threads: lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['threads'],
    resources:
        partition = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['partition'],
        mem_mb    = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['mem_mb'],
        time      = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['time'],
    shell: 
        """
        export OMP_NUM_THREADS={params.omp};
        RepeatMasker -s -pa {threads} -xsmall -lib {params.lib} -dir {params.out_dir} {input}
        """

rule funannotate_predict:
    input:
        os.path.join(config['TASK'], "masking/{ID}_sort.fa.masked")
    output: 
        out_dir  = directory(os.path.join(config['TASK'], "funannotate/{ID}")),
        genbank=os.path.join(config['TASK'], "funannotate/{ID}/predict_results/{ID}.gbk"),
        protein=os.path.join(config['TASK'], "funannotate/{ID}/predict_results/{ID}.proteins.fa"),
    params:
        gmes_path        = "/scratch/apps/gmes_linux_64_4",
        omp              = lambda wildcards, threads: threads,
        augustus_species = lambda wildcards: config['GENOMES'][wildcards.ID]['augustus'],
        busco_db         = lambda wildcards: config['GENOMES'][wildcards.ID]['busco_pred'],
        opt              = "",
        tmp_dir          = config['TMP_DIR']
    retries: 2
    threads: lambda wildcards, attempt: adjust_resources_medium(wildcards, attempt)['threads'],
    resources:
        partition = lambda wildcards, attempt: adjust_resources_medium(wildcards, attempt)['partition'],
        mem_mb    = lambda wildcards, attempt: adjust_resources_medium(wildcards, attempt)['mem_mb'],
        time      = lambda wildcards, attempt: adjust_resources_medium(wildcards, attempt)['time'],
    shell: 
        """
        export PATH="{params.gmes_path}:$PATH"

        conda run --no-capture-output --prefix /home/groups/Fungal/conda_env/funannotate \
            funannotate predict \
                -i {input} \
                -o {output.out_dir} \
                -s "{wildcards.ID}" \
                --augustus_species {params.augustus_species} \
                --busco_seed_species {params.augustus_species} \
                --busco_db {params.busco_db} \
                --name {wildcards.ID} \
                --cpus {threads} \
                --tmpdir {params.tmp_dir} \
                {params.opt}
        """

rule busco:
    input:
        os.path.join(config['TASK'], "funannotate/{ID}/predict_results/{ID}.proteins.fa")
    output: 
       directory(os.path.join(config['TASK'], "busco-after-funannotate/{ID}"))
    params:
        omp          = lambda wildcards, threads: threads,
        busco_db     = lambda wildcards: config['GENOMES'][wildcards.ID]['busco_pred'],
        lineages_dir = "/scratch/ref/06-busco/v5/lineages",
        opt=""
    retries: 2
    threads: lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['threads'],
    resources:
        partition = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['partition'],
        mem_mb    = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['mem_mb'],
        time      = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['time'],
    shell: 
        """
        export OMP_NUM_THREADS={params.omp};

        conda run --no-capture-output --prefix /home/groups/Fungal/conda_env/busco \
            busco -i {input} \
                -m proteins \
                -l {params.busco_db} \
                --cpu {threads} \
                -o {output} \
                {params.opt}
        """

rule antismash:
    input:
        os.path.join(config['TASK'], "funannotate/{ID}/predict_results/{ID}.gbk")
    output: 
       os.path.join(config['TASK'], "antismash/{ID}/{ID}.json"),
       os.path.join(config['TASK'], "antismash/{ID}/{ID}.gbk")
    params:
        omp     = lambda wildcards, threads: threads,
        out_dir = lambda wildcards, output: os.path.dirname(output[0])
    retries: 2
    threads: lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['threads']
    resources:
        partition = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['partition'],
        mem_mb    = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['mem_mb'],
        time      = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['time'],
    conda:
        "antismash8"
    shell: 
        """
        export OMP_NUM_THREADS={params.omp};

        antismash {input} \
            --cpus {threads} \
            --taxon fungi \
            --output-dir {params.out_dir} \
            --output-basename {wildcards.ID} \
            --genefinding-tool none
        """


rule interproscan:
    input:
        os.path.join(config['TASK'], "funannotate/{ID}/predict_results/{ID}.proteins.fa"),
    output: 
        tsv=os.path.join(config['TASK'], "interproscan/{ID}/{ID}.proteins.fa.tsv"),
        xml=os.path.join(config['TASK'], "interproscan/{ID}/{ID}.proteins.fa.xml"),
    params:
        out_dir=directory(os.path.join(config['TASK'], "interproscan/{ID}")),
        omp=lambda wildcards, threads: threads,
        int_dir="/home/groups/Fungal/marion/softwares/my_interproscan/interproscan-5.68-100.0",
        tmp_dir=os.path.join(config['TMP_DIR'], "{ID}"),
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

        mkdir -p {params.out_dir} {params.tmp_dir}
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
                --tempdir {params.tmp_dir} \
                {params.opt}
        """

rule signalP:
    input:
        os.path.join(config['TASK'], "masking/{ID}_sort.fa.masked")
    output: 
       os.path.join(config['TASK'], "signalP/{ID}/prediction_results.txt")
    params:
        omp     = lambda wildcards, threads: threads,
        out_dir = lambda wildcards, output: os.path.dirname(output[0])
    retries: 2
    threads: lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['threads']
    resources:
        partition = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['partition'],
        mem_mb    = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['mem_mb'],
        time      = lambda wildcards, attempt: adjust_resources_low(wildcards, attempt)['time'],
    shell: 
        """
        export OMP_NUM_THREADS={params.omp};

        conda run --no-capture-output --prefix /home/groups/Fungal/conda_env/signalp \
            signalp6 \
                -wp {threads} \
                --fastafile {input} \
                --organism eukarya \
                --output_dir {params.out_dir}
        """

rule funannotate_annotate:
    input:
        predictions  = os.path.join(config['TASK'], "funannotate/{ID}"),
        interproscan = os.path.join(config['TASK'], "interproscan/{ID}/{ID}.proteins.fa.xml"),
        antismash    = os.path.join(config['TASK'], "antismash/{ID}/{ID}.gbk"),
        signalP      = os.path.join(config['TASK'], "signalP/{ID}/prediction_results.txt")
    output: 
        # Explicitly declare the final, cleaned-up output files so Snakemake tracks them directly
        gff3 = os.path.join(config['TASK'], "funannotate-annotate/{ID}.gff3"),
        gbk  = os.path.join(config['TASK'], "funannotate-annotate/{ID}.gbk"),
        fa   = os.path.join(config['TASK'], "funannotate-annotate/{ID}.proteins.fa"),
        tbl  = os.path.join(config['TASK'], "funannotate-annotate/{ID}.tbl")
    wildcard_constraints:
        ID = "[A-Z0-9._]+"
    params:
        omp      = lambda wildcards, threads: threads,
        busco_db = lambda wildcards: config['GENOMES'][wildcards.ID]['busco_pred'],
        opt      = ""
    retries: 2
    threads: lambda wildcards, attempt: adjust_resources_medium(wildcards, attempt)['threads']
    resources:
        partition = lambda wildcards, attempt: adjust_resources_medium(wildcards, attempt)['partition'],
        mem_mb    = lambda wildcards, attempt: adjust_resources_medium(wildcards, attempt)['mem_mb'],
        time      = lambda wildcards, attempt: adjust_resources_medium(wildcards, attempt)['time'],
    shell: 
        """
        export OMP_NUM_THREADS={params.omp};

        # 1. Run funannotate (it will force outputs into {input.predictions}/annotate_results)
        conda run --no-capture-output --prefix /home/groups/Fungal/conda_env/funannotate \
            funannotate annotate \
                -i {input.predictions} \
                -o {input.predictions} \
                --iprscan {input.interproscan} \
                --antismash {input.antismash} \
                --signalp {input.signalP} \
                --isolate {wildcards.ID} \
                --busco_db {params.busco_db} \
                --cpus {threads} \
                --force \
                {params.opt}

        # 2. Ensure the destination directory exists
        mkdir -p $(dirname {output.gff3})

        # 3. Use bash 'mv' to pluck the nested, duplicated files out, rename them, and place them perfectly
        NESTED_DIR="{input.predictions}/annotate_results"
        
        mv "$NESTED_DIR/{wildcards.ID}_{wildcards.ID}.gff3"        "{output.gff3}"
        mv "$NESTED_DIR/{wildcards.ID}_{wildcards.ID}.gbk"         "{output.gbk}"
        mv "$NESTED_DIR/{wildcards.ID}_{wildcards.ID}.proteins.fa" "{output.fa}"
        mv "$NESTED_DIR/{wildcards.ID}_{wildcards.ID}.tbl"         "{output.tbl}"
        """