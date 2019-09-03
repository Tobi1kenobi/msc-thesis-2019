
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
# ============================ PARAMETERS =============================== #
# ======================================================================= #



# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

### Read LDSC results
dir.data <- here("out/CELLECT-LDSC-out/relevance")
files <- list.files(dir.data,pattern="*.txt")
files

data = tibble(File = files) %>%
  extract(File, "gwas", "__([^.]*)", remove = FALSE) %>%
  mutate(Data = lapply(file.path(dir.data, files), read_tsv)) %>%
  unnest(Data) %>% 
  separate(Name, c('benchmark', 'annotation'), '__')

# file.results <- here("results/prioritization_celltypes--tabula_muris.multi_gwas.csv.gz")
# df.ldsc_cts <- read_csv(file.results)

# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (single GWAS) ================ #
# ======================================================================= #

# gwas <- "BMI_UKBB_Loh2018"
# dataset_prefix <- "mousebrain_all"
# genomic_annotation_prefix <- get_genomic_annotation_prefix(dataset_prefix)
# 
# ### Loading BMI data
# file.ldsc_cts <- sprintf("/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/%s__%s.cell_type_results.txt", genomic_annotation_prefix, gwas)
# df.ldsc_cts <- load_ldsc_cts_results(file.ldsc_cts, dataset_prefix)
# df.ldsc_cts <- df.ldsc_cts %>% filter(sem=="sem_mean")
# 

# ======================================================================= #
# ========================= LOAD CELL-TYPE METADATA ===================== #
# ======================================================================= #

### Read annotation metadata
df.metadata <- get_metadata(dataset_prefix="tabula_muris")
df.metadata

### Add meta data
df.ldsc_cts <- data %>% left_join(df.metadata, by="annotation") # add meta data
df.ldsc_cts



# ======================================================================= #
# ========================= PREPARE DATA FOR ANALYSIS ===================== #
# ======================================================================= #

# Summing the negative logged p-values across all tissues and benchmarks
# With the aim of seeing if any particular tissues are over-represented
df.plot <- df.ldsc_cts %>%
  mutate(logPval = -log10(Coefficient_P_value)) %>% 
  filter(!grepl("[.][1-4]$",benchmark)) %>% 
  group_by(tissue, benchmark) %>% 
  summarise(sum_log_P_value = sum(logPval)) %>% 
  spread(benchmark, sum_log_P_value) %>% 
  mutate_at(vars(-1), funs(. - `tabula_muris.ESmu-NULL1.continuous.benchmark190722.0`)) %>%
  ungroup() %>% 
  gather(key = benchmark, value = sum.log.Pval.minus.null, -tissue)

df.plot


p.main <- ggplot(df.plot) +
  geom_bar(aes(y = sum.log.Pval.minus.null, x = tissue, fill=benchmark),stat="identity", position = "dodge") +
  scale_fill_manual(values=c("#b2e2e2", "#2ca25f", "#006d2c", # Null colours
                             "#fdcc8a", "#fc8d59", "#e34a33", # ESmu colours
                             "#5e3c99", # Finucane tstat
                             "#3182bd", "#08519c")) + # Skene specificity 
  theme_bw()

p.main
