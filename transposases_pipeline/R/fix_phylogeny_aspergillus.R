library(tidyverse)
library(ape)
library(ggtree)
library(treeio)
library(httpgd)

source("src/116-hpf-transposases/R/defaults.R")

hgd()

root       <- "/home/ailtonpcf/vast/proj/02-compost-microbes/cache"
tree_dir   <- "116-hpf-transposases-curatedHMM-funannotate-v2"
tree_file  <- "iqtree2/aspergillus/aft2_transposases.treefile"

afu_tree   <- read.tree(file.path(root, tree_dir, tree_file))
afu_tree_root  <- ape::root(afu_tree, outgroup = "aft1", resolve.root = TRUE)



iprsc_file  <- file.path(
    root,
    "116-hpf-transposases-curatedHMM-funannotate-v2",
    "interproscan/transposases-from-hmm",
    "aspergillus.faa.tsv"
)

taxonomy  <- "/home/ailtonpcf/draco/proj/02-compost-microbes/doc/14-phylogeny-accessions/aspergillus.tsv"

tax_meta  <- 
read_tsv(taxonomy)  %>% 
rename(
    genome_id=`Assembly Accession`,
    species=`Organism Name`
)  %>% 
select(genome_id, species)

library(ivs)

# Load and clean the raw data
raw_ipr <- read_tsv(iprsc_file, col_names = FALSE) %>%
  filter(!str_detect(X4, "PANTHER|Gene3D|SUPERFAMILY")) %>%
  mutate(X7 = as.numeric(X7), X8 = as.numeric(X8))  

  # test  <- "orf130721\\|GCF_000002655.1"

# 1. Isolate and strict-filter Pfam using E-values
pfam_clean <- raw_ipr %>%
# filter(str_detect(X1, test))  %>% 
  filter(X4 == "Pfam") %>%
  mutate(X9 = as.numeric(X9)) %>%
  filter(X9 < 1e-3)  %>% # Keep only solid Pfam hits
  mutate(X9 = as.character(X9))



# 2. Isolate others (SMART/ProSite) where X9 is a score, not an E-value
others_clean <- raw_ipr %>%
# filter(str_detect(X1, test))  %>% 
  filter(X4 != "Pfam")

# 3. Combine them back and run your overlap interval engine
afu_ani_pfam <- bind_rows(pfam_clean, others_clean) %>% 
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
    len=max(X8),
    X6 = paste(str_remove(X6, "\\..*"), collapse = "; ")
  )


tree_tbl <- treeio::as_tibble(afu_tree_root) %>% 
mutate(label = str_remove(label, "^[0-9]+_"))  %>% 
  left_join(
    afu_ani_pfam,
    by=c("label" = "X1")
  )  %>% 
  mutate(
    X6 = replace_na(X6, "No pfam predicted"),
    genome_id = str_extract(label, "GC[A-F]_[0-9]+.[0-9]"),
    X6 = replace_na(X6, "No pfam predicted"),
)  %>% 
left_join(tax_meta)  %>% 
mutate(
  label = str_remove(label, "\\|.*"),
  species = if_else(is.na(species), label, species),
  X6 = paste(label, X6, species, sep=";")
  )


tree_tbl  %>% 
filter(label=="aft1")  %>% 
pull(X6)

tree_tbl  %>% 
filter(str_detect(label, "UniRef"))  %>% 
pull(X6)
#   tree_tbl  %>% view()

library(randomcoloR)
set.seed(42)
cols_good <- distinctColorPalette(25) 

tree_tbl  %>% 
count(X6)


tree_tbl  %>% distinct(species)

species_order  <- c(
    "Aspergillus fischeri" ,
    "Aspergillus fumigatus",
    "Aspergillus oerlinghausenensis",
    "aft1",
    "UniRef90_A0A8H6Q8Q6|afu_transposase"
)

library(randomcoloR)
set.seed(42)
cols_good <- distinctColorPalette(105) 


phylo1 <- 
  tree_tbl %>% 
    as.treedata() %>% 
    ggtree(branch.length = "branch.length") +
    geom_treescale(x=0, y=-0.5) +
    geom_tippoint(aes(size=len)) +
    geom_nodelab(aes(label = label), hjust = 1.2, vjust = -0.3, size = 3) +
    geom_tiplab(aes(label = X6)) + #align = TRUE, linetype = "dotted", offset = 0.5
    scale_color_manual(values=cols_good) +
#     scale_shape_manual(
#     name = "species",
#     values = c(
#         "Aspergillus oerlinghausenensis" = 21,  # Filled Circle
#         "Aspergillus fumigatus"          = 24,  # Filled Triangle
#         "Aspergillus fischeri"           = 22,  # Filled Square
#         "aft1"                           = 3,  # Filled Circle
#         "UniRef90_A0A8H6Q8Q6|afu_transposase" = 7    # Square with an 'X' inside
#     )
#   ) +
    xlim_tree(2) +
    labs(color = "Phylum") +
    # theme(
    #   legend.position="bottom"
    # ) +
    guides(color = guide_legend(ncol=1))

apereira_save_plot(
  plot      = phylo1,
  plot_dir  = "compare_ddes/v1", 
  plot_name = "pogo_transposases_aspergillus", 
  limitsize = FALSE,
  scale=2,
  width     = 50, 
  height    = 100)





