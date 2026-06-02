rule blast_makedb:
    input:
        "/work/qi47rin/ref/sequences/aft2_domains_confirmed.faa",
        "/work/qi47rin/ref/sequences/aft1.faa"
    output: 
       os.path.join(config['TMP'], "indexes-protein/aft2_domains_confirmed.psq")
    params:
        prefix = lambda wc, output: output[0].replace(".psq", ""),
        opt    = "-dbtype prot",
        tmp    = "/scratch/tmp"
    singularity:
        "https://depot.galaxyproject.org/singularity/blast:2.17.0--h66d330f_0"
    shell:
        """
        # Merge to a real temporary file
        cat {input}  > {params.tmp}/temp_combined.faa
        
        # Run BLAST on the real file
        makeblastdb -in {params.tmp}/temp_combined.faa -dbtype prot -out {params.prefix}
        
        # Clean up the temporary file
        rm {params.tmp}/temp_combined.faa
        """

rule blastp:
    input:
        query    = os.path.join(config['TASK'], "filtered_hits/{ID}_hits.faa"),
        database = os.path.join(config['TMP'], "indexes-protein/aft2_domains_confirmed.psq")
    output: 
       os.path.join(config['TASK'], "blastp-asp-vs-atf2/{ID}.tsv")
    params:
        database = lambda wc, input: input.database.replace(".psq", ""),
        opt      = "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs'",
    singularity:
        "https://depot.galaxyproject.org/singularity/blast:2.17.0--h66d330f_0"
    threads: 8
    shell:
        """
        blastp -query {input.query} -db {params.database} -out {output} {params.opt} -num_threads {threads}
        """