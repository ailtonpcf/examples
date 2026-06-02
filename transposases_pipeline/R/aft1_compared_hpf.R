library(tidyverse)

source("src/116-hpf-transposases/R/defaults.R")
source("src/116-hpf-transposases/R/fix_phylogeny.R")

tax_file_hpf  <- "/home/ailtonpcf/draco/proj/02-compost-microbes/doc/14-phylogeny-accessions/human_pathogenic_fungi.tsv"

tax_asp_df  <- 
read_tsv(tax_file_hpf)  %>% 
rename(
    genome_id=`Assembly Accession`,
    species=`Organism Name`
)  %>% 
select(genome_id, species)


# 1. Define the directory to scan (change "." to your actual folder path)
asp_transposases_atf2 <- "blastp-asp-vs-atf2/"

# 2. List all TSV files matching the GCA_* pattern
tsv_files <- list.files(
  path = file.path(root, task, asp_transposases_atf2), 
  pattern = "^GC[A-F]_.*\\.tsv$", 
  full.names = TRUE
)

# 1. Updated to match your exact 13 -outfmt 6 columns
blast_cols <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", 
  "gapopen", "qstart", "qend", "sstart", "send", 
  "evalue", "bitscore", "qcovs"
)

combined_data <- tsv_files %>%
  map_df(~ {
    read_tsv(
      .x, 
      col_names = blast_cols, 
      show_col_types = FALSE,
      # 2. Strict type mapping for all 13 columns to prevent mixing double/character
      col_types = cols(
        qseqid   = col_character(),
        sseqid   = col_character(),
        pident   = col_double(),      
        length   = col_integer(),
        mismatch = col_integer(),
        gapopen  = col_integer(),
        qstart   = col_integer(),
        qend     = col_integer(),
        sstart   = col_integer(),
        send     = col_integer(),
        evalue   = col_double(),
        bitscore = col_double(),
        qcovs    = col_integer()      # Added matching type for your 13th column
      )
    ) %>%
      separate_wider_delim(
        cols = qseqid,
        delim = "|",
        names = c("protein_orf", "genome_id"),
        too_many = "merge"
      )
  })

# Import interproscan domains

iprsc_file  <- file.path(
    root,
    "116-hpf-transposases-curatedHMM-funannotate-v2",
    "interproscan/transposases-from-hmm",
    "aspergillus.faa.tsv"
)
library(ivs)

# Load and clean the raw data
aft2_ipr <- read_tsv(iprsc_file, col_names = FALSE) %>%
  filter(!str_detect(X4, "PANTHER|Gene3D|SUPERFAMILY")) %>%
  mutate(X7 = as.numeric(X7), X8 = as.numeric(X8))  

  # test  <- "orf130721\\|GCF_000002655.1"

# 1. Isolate and strict-filter Pfam using E-values
pfam_clean_aft2 <- aft2_ipr %>%
# filter(str_detect(X1, test))  %>% 
  filter(X4 == "Pfam") %>%
  mutate(X9 = as.numeric(X9)) %>%
  filter(X9 < 1e-3)  %>% 
  mutate(X9 = as.character(X9))



# 2. Isolate others (SMART/ProSite) where X9 is a score, not an E-value
others_clean_aft2 <- aft2_ipr %>%
# filter(str_detect(X1, test))  %>% 
  filter(X4 != "Pfam")

# 3. Combine them back and run your overlap interval engine
afu_aft1_domains <- bind_rows(pfam_clean_aft2, others_clean_aft2) %>% 
  filter(X7 != X8) %>%
  select(X1, X6, X7, X8, X9) %>%
  mutate(iv = iv(X7, X8)) %>%
  
  group_by(X1) %>%
  mutate(group = iv_identify_group(iv)) %>%
  
  # CRITICAL: Since X9 is mixed, we prioritize Pfam matches in clashes 
  # by sorting by database type or presence of 'PF' in the name
  group_by(X1, group) %>%
  arrange(desc(str_detect(X6, "PF")), .by_group = TRUE) %>%
  slice(1) %>% 
  
  ungroup() %>% 
  arrange(X1, X7) %>%
  group_by(X1) %>%
  summarise(
    X6 = paste(str_remove(X6, "\\..*"), collapse = "; ")
  )  

signif_df  <- 
combined_data  %>% 
filter(evalue < 0.05)  %>% 
filter(str_detect(sseqid, "aft1"), qcovs>90) 

hpf_aft1_metrics  <- 
tax_asp_df  %>% 
left_join(signif_df)

aft1_no_hpf  <- 
hpf_aft1_metrics  %>% 
filter(!is.na(pident))

valid_hpf  <- 
hpf_aft1_metrics  %>% 
# select(1:3,pident)  %>% 
filter(!is.na(pident))  %>% 
left_join(
    afu_aft1_domains, by=c("sseqid"="X1"))

valid_hpf  %>% 
group_by(species)  %>% 
summarise(
  mean_pident=mean(pident),
  mean_qcovs=mean(qcovs),
  mean_len = mean(length),
  n_seqs = n_distinct(protein_orf)
)  %>% 
arrange(desc(mean_pident))

genomes_oder <- valid_hpf %>% 
  group_by(species) %>% 
  summarise(median_pident = median(pident, na.rm = TRUE)) %>% 
  arrange(desc(median_pident)) %>%             # Sort from highest to lowest identity
  pull(species) %>% 
  unique()



p_aft1_group  <- 
valid_hpf  %>% 
mutate(species = factor(species, levels=genomes_oder))  %>% 
ggplot(aes(species, pident)) +
geom_boxplot() +
theme(axis.text.x = element_text(angle=90, hjust=0.9))

apereira_save_plot(
  plot      = p_aft1_group,
  plot_dir  = "compare_ddes/v1", 
  plot_name = "aft1_pident_hpf", 
  limitsize = FALSE,
  width     = 7, 
  height    =5)

