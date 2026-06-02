# Transposase Identification and Phylogeny Pipeline

This directory contains a modular Snakemake workflow focused on downloading fungal genomes, predicting genes, profiling specific transposase domains using Hidden Markov Models (HMMs), functional validation, and phylogenetic reconstruction.

The project specifically focuses on parsing genome tables for *Aspergillus* species (and specific subgroups like `AFU_AFI_AOE`) to extract and analyze target transposases.

---

## 📋 Table of Contents
1. [Workflow Architecture](#workflow-architecture)
2. [Configuration and Input Parsing](#configuration-and-input-parsing)
3. [Pipeline Steps](#pipeline-steps)
4. [Target Outputs](#target-outputs)

---

## 🏗 Workflow Architecture

The pipeline is highly modular, leveraging Snakemake's `module` feature to incorporate rules from separate workflows located in the `modules/` subdirectory.

### Imported Modules
* **`download_genomes`**: Fetches the genome assemblies.
* **`funannotate`**: Standard eukaryotic gene prediction (optimized for fungi).
* **`hmm_build_curated`**: Builds profile HMMs from curated sequence alignments.
* **`hmm_search_transposase2`**: Searches predicted proteins against the built HMM profiles.
* **`interproscan_solo`**: Functional annotation and validation via InterProScan.
* **`phylogeny`**: Alignments and tree generation (via IQ-TREE 2).
* **`search_TIRs`**: Identifies Terminal Inverted Repeats (TIRs).
* **`blastp`**: Local sequence alignment (specifically *Aspergillus* vs ATF2 queries).

---

## ⚙️ Configuration and Input Parsing

Before running rules, the pipeline reads a standard metadata table provided in the config file. It dynamically parses target **Assembly Accessions** and sets metadata parameters for `augustus` species and `BUSCO` thresholds.

```python
# Internal Python parsing logic within the Snakefile:
for tbl in config['ASPERGILLUS']:
    df = pd.read_csv(tbl, sep="\t")
    for _, row in df.iterrows():
        ID = str(row['Assembly Accession'])
        config['GENOMES'][ID] = {
            'augustus': row.get('augustus_sp', None),
            'busco_pred': row.get('busco_prediction', None),
            'busco_comp': row.get('busco_completeness', None)
        }
```

## Verify the file depedencies:

snakemake --snakefile main.smk -nr

## Running the pipeline:

bash run_transposases.sh