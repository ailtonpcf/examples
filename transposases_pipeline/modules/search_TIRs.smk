rule parse_hits_to_bed:
    input:
        tsv = os.path.join(config['TASK'], "hmmsearch/{ID}.tsv")
    output:
        bed = os.path.join(config['TASK'], "tir_analysis/{ID}_hits.bed")
    shell:
        """
        grep -v "^#" {input.tsv} | awk 'NF {{
            chrom = ""
            coords_val = ""
            
            # Scan columns to extract the real chromosome and coordinates
            for (i=1; i<=NF; i++) {{
                if ($i ~ /^source=/) {{
                    chrom = $i;
                    gsub(/source=/, "", chrom);
                }}
                if ($i ~ /^coords=/) {{
                    coords_val = $i;
                    gsub(/coords=/, "", coords_val);
                }}
            }}
            
            # If we found both, format them properly into a BED file
            if (chrom != "" && coords_val != "") {{
                split(coords_val, pos, /\\.\\./);
                
                # Sort coordinates so start is always smaller than end
                start = (pos[1] < pos[2]) ? pos[1] : pos[2];
                end   = (pos[1] < pos[2]) ? pos[2] : pos[1];
                
                # Output standard BED: Chromosome, Start, End, Protein_ID
                print chrom "\\t" (start - 1) "\\t" end "\\t" $1;
            }}
        }}' > {output.bed}
        """

# Rule 2: Expand the boundaries by 2,000 bp to capture the flanking regions
rule expand_flanking_regions:
    input:
        bed = os.path.join(config['TASK'], "tir_analysis/{ID}_hits.bed"),
        fna = os.path.join(config['TASK'], "masking/{ID}_sort.fa.masked")
    output:
        expanded_bed = os.path.join(config['TASK'], "tir_analysis/{ID}_expanded.bed")
    singularity:
        "https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0"
    shell:
        """
        # Create a genome size file on the fly so bedtools doesn't run off chromosome edges
        awk '/^>/ {{if (seqlen) print seqlen; print substr($1,2); seqlen=0; next}} {{seqlen+=length($0)}} END {{print seqlen}}' {input.fna} | \
        paste - - > {input.fna}.sizes

        # Slop expands the BED coordinates by 2000bp symmetrically
        bedtools slop -i {input.bed} -g {input.fna}.sizes -b 2000 > {output.expanded_bed}

        rm {input.fna}.sizes
        """

# Rule 3: Extract the nucleotide sequences of these expanded blocks
rule extract_expanded_sequences:
    input:
        expanded_bed = os.path.join(config['TASK'], "tir_analysis/{ID}_expanded.bed"),
        fna = os.path.join(config['TASK'], "masking/{ID}_sort.fa.masked")
    output:
        flank_fna = os.path.join(config['TASK'], "tir_analysis/{ID}_flanks.fna")
    singularity:
        "https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0"
    shell:
        """
        # Extract the fasta sequence matching the expanded coordinates
        bedtools getfasta -fi {input.fna} -bed {input.expanded_bed} -fo {output.flank_fna} -name
        """

# Rule 4: Run a self-BLAST to find Terminal Inverted Repeats (TIRs)
rule verify_tirs:
    input:
        flank_fna = os.path.join(config['TASK'], "tir_analysis/{ID}_flanks.fna")
    output:
        blast_out = os.path.join(config['TASK'], "tir_analysis/{ID}_tir_results.txt")
    singularity:
        "https://depot.galaxyproject.org/singularity/blast:2.17.0--h66d330f_0"
    shell:
        """
        if [ -s {input.flank_fna} ]; then
            # 1. Build a blast database out of your extracted chunks
            makeblastdb -in {input.flank_fna} -dbtype nucl -out {input.flank_fna}.db

            # 2. Blast the sequences against themselves
            # We look for matches on the OPPOSITE strand (minus) to find inverted repeats
            # qcov_hsp_perc prevents full self-identity matches from drowning out small TIR flags
            blastn -query {input.flank_fna} \
                   -db {input.flank_fna}.db \
                   -strand minus \
                   -outfmt "6 qseqid sseqid length pident qstart qend sstart send evalue" \
                   -evalue 1e-3 > {output.blast_out}

            rm {input.flank_fna}.db.*
        else
            touch {output.blast_out}
        fi
        """

rule filter_active_transposons:
    input:
        blast_out = os.path.join(config['TASK'], "tir_analysis/{ID}_tir_results.txt")
    output:
        active_list = os.path.join(config['TASK'], "tir_analysis/{ID}_active_transposons.txt")
    shell:
        """
        # awk logic:
        # $1 != $2 : Exclude self-matches to look strictly at duplications across the genome
        # $3 >= 1000 : Ensure a large part of the transposon core moved (avoids background noise)
        # $4 >= 95.0 : Strict identity cutoff indicating recent, clean transposition activity
        
        awk '($1 != $2) && ($3 >= 1000) && ($4 >= 95.0) {{print $1}}' {input.blast_out} | \
        sed 's/::.*//' | sort -u > {output.active_list}
        """

rule find_disrupted_genes:
    input:
        bed = os.path.join(config['TASK'], "tir_analysis/{ID}_hits.bed"),
        gff = os.path.join(config['TASK'], "funannotate-annotate/{ID}.gff3"),
    output:
        disrupted = os.path.join(config['TASK'], "functional_analysis/{ID}_mutated_genes.txt")
    singularity:
        "https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2"
    shell:
        """
        # 1. Filter GFF3 to ONLY extract 'gene' structural blocks (Column 3)
        # 2. Pass it into bedtools intersect via a bash stream (<( ... ))
        
        bedtools intersect \
            -a <(awk '$3 == "gene"' {input.gff}) \
            -b {input.bed} \
            -wa -wb > {output.disrupted}
        """