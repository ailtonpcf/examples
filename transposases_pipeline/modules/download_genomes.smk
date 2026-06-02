rule ncbi_datasets_download_genome:
    output:
        os.path.join(config['TASK'], "genomes/{ID}.fasta")
    params:
        genome_dir = lambda wildcards, output: os.path.join(config['workdir'], os.path.dirname(output[0]))
    conda: "ncbi-datasets-cli"
    retries: 2
    shell:
        """
        rm -rf {wildcards.ID}
        mkdir -p {wildcards.ID} {params.genome_dir}
        cd {wildcards.ID}
        datasets download genome accession {wildcards.ID}
        unzip ncbi_dataset.zip
        mv ncbi_dataset/data/{wildcards.ID}/*.fna {params.genome_dir}
        cd ..
        rm -rf {wildcards.ID}
        cd {params.genome_dir}
        mv {wildcards.ID}* {wildcards.ID}.fasta
        """