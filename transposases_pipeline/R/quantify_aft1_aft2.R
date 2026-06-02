library(tidyverse)

source("src/116-hpf-transposases/R/defaults.R")

abundances  <- read_tsv(
    file.path(
        root,
        "120-quantify-aft2/abundances/merged.tsv"
    )
)

tidied_abundances  <- 
abundances  %>% 
pivot_longer(
    cols=-contigname,
    names_to="sample",
    values_to="abundance"
)



p_abundance_atf1_atf2  <- 
tidied_abundances  %>% 
ggplot(aes(contigname, abundance)) +
geom_boxplot()

library(rstatix)

tidied_abundances  %>% 
wilcox_test(
    formula=abundance ~ contigname,
    ref.group="aft1",
    alternative="greater")  %>% 
add_significance()

apereira_save_plot(
  plot      = p_abundance_atf1_atf2,
  plot_dir  = "compare_ddes/v1", 
  plot_name = "aft1_aft2_abundance", 
  limitsize = FALSE,
  scale=1,
  width     = 9, 
  height    = 7)
