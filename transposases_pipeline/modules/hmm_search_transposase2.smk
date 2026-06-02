rule esl_translate:
    input:
        hmm = os.path.join(config['TASK'], "transposase-hmm/aft2_domains_confirmed.hmm"),
        fna = os.path.join(config['TASK'], "masking/{ID}_sort.fa.masked")
    output:
        faa = os.path.join(config['TASK'], "translated_proteins/{ID}.faa")
    params:
        opt = "-E 1e-5"
    threads: 8
    singularity:
        "https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_2"
    shell:
        """
        # 1. Translate the genome and save the entire proteome directly to a file
        esl-translate {input.fna} > {output.faa}
        """

rule add_ID_to_header:
    input:
        os.path.join(config['TASK'], "translated_proteins/{ID}.faa")
    output:
        os.path.join(config['TASK'], "translated_proteins_renamed/{ID}.faa")
    shell:
        """
        sed "s/^>\\([^ ]*\\)\\(.*\\)/>\\1|{wildcards.ID}\\2/" {input} > {output}
        """

rule hmmsearch:
    input:
        hmm = os.path.join(config['TASK'], "transposase-hmm/aft2_domains_confirmed.hmm"),
        faa = os.path.join(config['TASK'], "translated_proteins_renamed/{ID}.faa")
    output:
        tsv = os.path.join(config['TASK'], "hmmsearch/{ID}.tsv"),
    params:
        opt = "-E 1e-100"
    threads: 8
    singularity:
        "https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_2"
    shell:
        """
        # 2. Search your protein HMM against the proteins you just generated
        hmmsearch --cpu {threads} {params.opt} --tblout {output.tsv} {input.hmm} {input.faa}
        """

rule filter_matching_proteins:
    input:
        tsv = os.path.join(config['TASK'], "hmmsearch/{ID}.tsv"),
        faa = os.path.join(config['TASK'], "translated_proteins_renamed/{ID}.faa")
    output:
        filtered_faa = os.path.join(config['TASK'], "filtered_hits/{ID}_hits.faa")
    params:
        evalue = 1e-100,
    singularity:
        "https://depot.galaxyproject.org/singularity/seqkit:2.8.1--h9ee0642_0"
    shell:
        """
        # 1. Let awk handle the comment skipping and filtering safely
        awk -v ev={params.evalue} '!/^#/ && NF && $5 <= ev {{print $1}}' {input.tsv} | \
        sort -u > {output.filtered_faa}.ids

        # 2. Extract hits and append the ID safely
        if [ -s {output.filtered_faa}.ids ]; then
            seqkit grep -f {output.filtered_faa}.ids {input.faa} > {output.filtered_faa}
        else
            touch {output.filtered_faa}
        fi

        rm {output.filtered_faa}.ids
        """