############### SYNOPSIS ####################
# IDEA:
# Try and find out if any of the three methods have a bias
# towards any particular cell types by looking across all GWAS
# for trends in both Tabula Muris and Mousebrain datasets. Based
# on reading Skene paper and looking at Tabula Muris data I expect to
# see MSNs enriched by Skene in Mousebrain and Pancreas enriched in

# AIMS:
#
#
#

# PLOT:
#
#
# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library


library(RColorBrewer)
library(patchwork)

library(stringr)
library(colorspace)

source(here("src/publication/lib-load_pub_lib_functions.R"))


setwd(here("src/benchmark"))

# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

### Read LDSC results
dir.data <- here("out/CELLECT-LDSC-out/relevance")
esmu.files <- list.files(dir.data,pattern=".*ESmu.continuous[.].*.txt")
esmu.files
null.files <- list.files(dir.data,pattern=".*ESmu-NULL1[.]continuous.*.txt")
null.files

files <- c(esmu.files, null.files)
files

data = tibble(File = files) %>%
  extract(File, "Trait identifier", "__([^.]*)", remove = FALSE) %>%
  mutate(Data = lapply(file.path(dir.data, files), read_tsv)) %>%
  unnest(Data) %>% 
  separate(Name, c('benchmark', 'annotation'), '__') %>% 
  group_by(`Trait identifier`, benchmark) %>% 
  mutate(rank = dense_rank(desc(Coefficient_P_value))) %>% 
  ungroup

data

### Which GWAS am I interested in?
GWAS.info <- read_csv('190713_GWAS_Timshel (tobi edit).csv')
GWAS.info

# Would be 39 but I remove 4 traits that violate LDSC assumptions too
data.only.35 <- data %>%
  left_join(GWAS.info, by="Trait identifier") %>% 
  filter(FLAG_final & GWAS.trait.representative & is.na(`GWAS violates LDSC assumptions`))

data.only.35

# ======================================================================= #
# ========================= LOAD CELL-TYPE METADATA ===================== #
# ======================================================================= #

### Read annotation metadata
df.metadata <- get_metadata(dataset_prefix="tabula_muris")
df.metadata

### Add meta data
df.ldsc_cts <- data.only.35 %>% left_join(df.metadata, by="annotation") # add meta data
df.ldsc_cts

# ======================================================================= #
# ========================= SUMMARISE AND PREP DATA ===================== #
# ======================================================================= #

# Normalise sum of ranks by the number of cells in each category
df.rank_count <- df.ldsc_cts %>% 
  filter(benchmark == 'tabula_muris.ESmu.continuous.benchmark190722',
         `Trait category` %in% c('Neurological','Neuropsychiatric') ) %>% 
  group_by(tissue) %>% 
  summarise(rank.sum = sum(-log10(Coefficient_P_value)), n.cells = n_distinct(cell_type)) %>% 
  mutate(norm.rank = rank.sum/n.cells)

df.rank_count

# Calculate the same for the null
df.null.rank_count <- df.ldsc_cts %>% 
  filter(benchmark %in% paste0('tabula_muris.ESmu-NULL1.continuous.benchmark190722.', 0:4),
         `Trait category` %in% c('Neurological','Neuropsychiatric')) %>% 
  group_by(tissue) %>% 
  summarise(rank.sum = sum(-log10(Coefficient_P_value)), n.cells = n_distinct(cell_type)) %>% 
  mutate(norm.rank = rank.sum/(n.cells * 5))

df.null.rank_count

df.plot <- df.rank_count %>% 
  mutate(norm.rank = df.rank_count$norm.rank - df.null.rank_count$norm.rank)

df.plot
# =================================================================== #
# ========================= MAKE PLOTS ============================== #
# =================================================================== #

p <- ggplot(df.plot, aes(x = factor(tissue),y=norm.rank, fill = tissue)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=7)) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(name="Paired", n = 8))(20))
p
