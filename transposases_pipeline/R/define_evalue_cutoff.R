library(tidyverse)

root       <- "/home/ailtonpcf/vast/proj/02-compost-microbes/cache"
iprsc_file  <- file.path(
    root,
    "116-hpf-transposases-curatedHMM-funannotate-v2",
    "interproscan/transposases-from-hmm",
    "aspergillus.faa.tsv"
)

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
    X6 = paste(str_remove(X6, "\\..*"), collapse = "; ")
  )

# Criando o vetor nomeado com Target Name e E-value (full sequence)
hmmer_evalues <- c(
  "orf130721|GCF_000002655.1" = 7e-276,
  "orf693855|GCF_000002655.1" = 7e-276,
  "orf411695|GCF_000002655.1" = 1.6e-275,
  "orf730477|GCF_000002655.1" = 5.2e-275,
  "orf156650|GCF_000002655.1" = 6e-275,
  "orf278067|GCF_000002655.1" = 6e-275,
  "orf654631|GCF_000002655.1" = 6e-275,
  "orf782018|GCF_000002655.1" = 7.3e-275,
  "orf383269|GCF_000002655.1" = 1.3e-274,
  "orf124324|GCF_000002655.1" = 2.8e-274,
  "orf117035|GCF_000002655.1" = 9.9e-273,
  "orf513404|GCF_000002655.1" = 9.9e-273,
  "orf417837|GCF_000002655.1" = 9.1e-222,
  "orf126868|GCF_000002655.1" = 1.9e-221,
  "orf602395|GCF_000002655.1" = 8.2e-188,
  "orf656164|GCF_000002655.1" = 2.4e-187,
  "orf249|GCF_000002655.1"    = 4.3e-161,
  "orf414106|GCF_000002655.1" = 1.3e-160,
  "orf480241|GCF_000002655.1" = 1.6e-142,
  "orf734113|GCF_000002655.1" = 1.6e-142,
  "orf828410|GCF_000002655.1" = 6.1e-141,
  "orf480231|GCF_000002655.1" = 6.9e-123,
  "orf734103|GCF_000002655.1" = 6.9e-123,
  "orf828396|GCF_000002655.1" = 1.2e-122,
  "orf671491|GCF_000002655.1" = 6e-121,
  "orf414097|GCF_000002655.1" = 8e-104,
  "orf238|GCF_000002655.1"    = 7.1e-103,
  "orf802458|GCF_000002655.1" = 1.2e-99,
  "orf656173|GCF_000002655.1" = 1.6e-79,
  "orf320297|GCF_000002655.1" = 1.2e-78,
  "orf528542|GCF_000002655.1" = 1.3e-66,
  "orf284576|GCF_000002655.1" = 1.7e-61,
  "orf411622|GCF_000002655.1" = 9.9e-54,
  "orf407816|GCF_000002655.1" = 1e-51,
  "orf282610|GCF_000002655.1" = 7.2e-46,
  "orf407810|GCF_000002655.1" = 1.1e-43,
  "orf528548|GCF_000002655.1" = 3.1e-42,
  "orf126872|GCF_000002655.1" = 5.2e-42,
  "orf417845|GCF_000002655.1" = 8.6e-42,
  "orf284570|GCF_000002655.1" = 2.4e-41,
  "orf312416|GCF_000002655.1" = 5.8e-41,
  "orf282603|GCF_000002655.1" = 5.1e-39,
  "orf315915|GCF_000002655.1" = 3.8e-32,
  "orf701032|GCF_000002655.1" = 2.3e-29,
  "orf56848|GCF_000002655.1"  = 3.1e-26,
  "orf325474|GCF_000002655.1" = 7.9e-25,
  "orf533881|GCF_000002655.1" = 1.3e-24,
  "orf667912|GCF_000002655.1" = 1.3e-24,
  "orf128738|GCF_000002655.1" = 1.4e-24,
  "orf186121|GCF_000002655.1" = 1.4e-24,
  "orf350622|GCF_000002655.1" = 1.4e-24,
  "orf437015|GCF_000002655.1" = 1.4e-24,
  "orf450279|GCF_000002655.1" = 1.4e-24,
  "orf548556|GCF_000002655.1" = 1.4e-24,
  "orf565365|GCF_000002655.1" = 1.4e-24,
  "orf730332|GCF_000002655.1" = 1.4e-24,
  "orf758918|GCF_000002655.1" = 1.4e-24,
  "orf762023|GCF_000002655.1" = 1.4e-24,
  "orf91007|GCF_000002655.1"  = 1.4e-24,
  "orf86480|GCF_000002655.1"  = 2.5e-24,
  "orf463128|GCF_000002655.1" = 2.8e-24,
  "orf734865|GCF_000002655.1" = 4.8e-24,
  "orf237159|GCF_000002655.1" = 1.4e-23,
  "orf821511|GCF_000002655.1" = 1.4e-23,
  "orf355874|GCF_000002655.1" = 1.6e-23,
  "orf299706|GCF_000002655.1" = 1.1e-18,
  "orf800827|GCF_000002655.1" = 1.4e-18,
  "orf312418|GCF_000002655.1" = 4.2e-18,
  "orf408391|GCF_000002655.1" = 3.9e-17,
  "orf602887|GCF_000002655.1" = 7e-17,
  "orf291822|GCF_000002655.1" = 7.4e-17,
  "orf523440|GCF_000002655.1" = 9.1e-17,
  "orf424020|GCF_000002655.1" = 1.8e-15,
  "orf750732|GCF_000002655.1" = 1.9e-15,
  "orf424025|GCF_000002655.1" = 9.3e-15,
  "orf408389|GCF_000002655.1" = 2.9e-14,
  "orf734861|GCF_000002655.1" = 1.1e-13,
  "orf248734|GCF_000002655.1" = 2.2e-13,
  "orf299696|GCF_000002655.1" = 2.7e-13,
  "orf320281|GCF_000002655.1" = 1.9e-12,
  "orf427619|GCF_000002655.1" = 2e-12,
  "orf427608|GCF_000002655.1" = 4.7e-12,
  "orf802444|GCF_000002655.1" = 9.7e-12,
  "orf427615|GCF_000002655.1" = 1.1e-11,
  "orf676248|GCF_000002655.1" = 1.4e-11,
  "orf427610|GCF_000002655.1" = 1.5e-11,
  "orf231822|GCF_000002655.1" = 4.2e-11,
  "orf800816|GCF_000002655.1" = 5.9e-11
)  %>% enframe()

# Exemplo de como acessar os dados:
# Para ver o E-value de um target específico:
# hmmer_evalues["orf130721|GCF_000002655.1"]

hmmer_evalues  %>% 
left_join(
    afu_ani_pfam,
    by = c("name" = "X1")
)

defined_cutoff  <- 1e-100
# It was defined 1e-100 because right below these values, proteins started to lose all 3 domains.
# Fischeri has more copies of aft2 > a fumigatus. oerligausinesis doesnt have.
# Do we find aft2 in all aspergillus species?