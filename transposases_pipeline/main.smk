configfile: "/scratch/proj/02-compost-microbes/src/116-hpf-transposases/smk.config"
workdir: config['workdir']
import pandas as pd
include: "/scratch/proj/00-default/smk-functions/resources.py"

module download:
    snakefile: "modules/download_genomes.smk"
    config: config
use rule * from download

module gene_prediction:
    snakefile: "modules/funannotate.smk"
    config: config

module hmm_build:
    snakefile: "modules/hmm_build_curated.smk"
    config: config
use rule * from hmm_build

module hmm_search:
    snakefile: "modules/hmm_search_transposase2.smk"
    config: config
use rule * from hmm_search

module protein_annotation:
    snakefile: "modules/interproscan_solo.smk"
    config: config

module gene_prediction:
    snakefile: "modules/funannotate.smk"
    config: config
use rule * from gene_prediction

module phylogeny:
    snakefile: "modules/phylogeny.smk"
    config: config
use rule * from phylogeny

module TIRs:
    snakefile: "modules/search_TIRs.smk"
    config: config
use rule * from TIRs

module blastp:
    snakefile: "modules/blastp.smk"
    config: config
use rule * from blastp

# Initialize GENOMES dict if not present
config.setdefault('GENOMES', {})

# Read all genome tables listed in config['TABLES']
for tbl in config['ASPERGILLUS']:
    df = pd.read_csv(tbl, sep="\t")
    for _, row in df.iterrows():
        ID = str(row['Assembly Accession'])
        config['GENOMES'][ID] = {
            'augustus': row.get('augustus_sp', None),
            'busco_pred': row.get('busco_prediction', None),
            'busco_comp': row.get('busco_completeness', None)
        }

IDS = list(config['GENOMES'].keys())

rule all:
    input:
        expand(os.path.join(config['TASK'], "interproscan/transposases-from-hmm/{GENUS}.{EXT}.tsv"),GENUS="aspergillus", EXT='faa'),
        expand(os.path.join(config['TASK'], "iqtree2/{GENUS}/aft2_transposases.treefile"), GENUS='aspergillus'),
        expand(os.path.join(config['TASK'], "blastp-asp-vs-atf2/{ID}.tsv"), ID=IDS)
    default_target: True

rule concatenate_transposase_for_hmm:
    input:
        "cache/116-hpf-transposases-hpf/genes/transposase.faa",
         expand(os.path.join(config['TASK'], "filtered_hits/{ID}_hits.faa"), ID=IDS)
    output: 
       os.path.join(config['TASK'], "concat-seqs/{GENUS}.faa")
    shell:
        """
        cat {input} > {output}
        """

use rule concatenate_transposase_for_hmm as concatenate_transposase_3aspergillus with:
    input:
        "cache/116-hpf-transposases-hpf/genes/transposase.faa",
         expand(os.path.join(config['TASK'], "filtered_hits/{ID}_hits.faa"), ID=config['AFU_AFI_AOE'])
    output: 
       os.path.join(config['TASK'], "concat-seqs/afu_afis_aoerl.faa")

use rule interproscan from protein_annotation as validate_transposase with:
    input:
        os.path.join(config['TASK'], "concat-seqs/{GENUS}.faa")
    output: 
        tsv=os.path.join(config['TASK'], "interproscan/transposases-from-hmm/{GENUS}.{EXT}.tsv"),
        xml=os.path.join(config['TASK'], "interproscan/transposases-from-hmm/{GENUS}.{EXT}.xml"),

rule concatenate_found_transposases_for_phylogeny:
    input:
        expand(os.path.join(config['TASK'], "transposases-from-hmm/{HPF}__{GENUS}.faa"), HPF = IDS, GENUS="aspergillus"),
        "cache/116-hpf-transposases-hpf/genes/transposase.faa"
    output: 
       os.path.join(config['TASK'], "transposases-from-hmmsearch/merged.faa")
    shell:
        """
        cat {input} > {output}
        """

use rule interproscan from protein_annotation as annotate_transposases with:
    input:
        os.path.join(config['TASK'], "transposases-from-hmmsearch/merged.faa")
    output: 
        tsv=os.path.join(config['TASK'], "results/interproscan/merged.faa.tsv"),
        xml=os.path.join(config['TASK'], "results/interproscan/merged.faa.xml"),