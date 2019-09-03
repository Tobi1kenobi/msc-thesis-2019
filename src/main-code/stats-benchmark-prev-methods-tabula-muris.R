############### SYNOPSIS ###################
# Make table to quantify the three benchmark criteria:
#
# •	Has the most significant p-value of all the metrics AND is Bonferroni significant
# •	Cell type reaches Bonferroni significance that we deem relevant based on our knowledge of the biology of the trait
# •	No cell types reach Bonferroni significance that we deem irrelevant based on our knowledge of the biology of the trait
#
# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)
library(here)


setwd(here("src/benchmark"))

source(here("src/lib/load_functions.R")) # load sc-genetics library

`%nin%` = Negate(`%in%`)

# ==================================================================================== #
# ================== LOAD LDSC CTS RESULTS (multi GWAS, multi metric) ================ #
# ==================================================================================== #

### Read LDSC results
dir.data <- here("out/CELLECT-LDSC-out/other-methods")
files <- list.files(dir.data,pattern="*.txt")
files

data = tibble(File = files) %>%
  extract(File, "Trait identifier", "__([^.]*)", remove = FALSE) %>%
  mutate(Data = lapply(file.path(dir.data, files), read_tsv)) %>%
  unnest(Data) %>% 
  separate(Name, c('benchmark_date', 'annotation'), '__') %>% 
  extract(benchmark_date, "benchmark", "([^.]*[.][^.]*[.][^.]*)", remove = TRUE)

data

### Which GWAS am I interested in?
GWAS.info <- read_csv('190713_GWAS_Timshel (tobi edit).csv')
GWAS.info



# Would be 39 but I remove 4 traits that violate LDSC assumptions too
data.only.35 <- data %>%
  left_join(GWAS.info, by="Trait identifier") %>% 
  filter(FLAG_final & GWAS.trait.representative & is.na(`GWAS violates LDSC assumptions`))


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
# ========================= MAKE BENCHMARK TABLE ======================== #
# ======================================================================= #

# First see how many traits are in the five categories we'll use (because the others
# have difficult to define target cell types)
df.ldsc_cts %>%
  filter(`Trait category` %in% c("Neuropsychiatric", "Behavioral", "Immune",
                                 "Cardiovascular", "Lipids")) %>%
  group_by(`Trait identifier`) %>%
  summarise(n())

# Which metrics produce the most significant cell types?
most_sig.df <- df.ldsc_cts %>% 
  group_by(`Trait full name`) %>% 
  filter(Coefficient_P_value == min(Coefficient_P_value),
         Coefficient_P_value < 0.05/115,
         `Trait category` %in% c("Neuropsychiatric", "Behavioral", "Immune",
                                 "Cardiovascular", "Lipids")) %>% 
  ungroup %>% 
  select(`Trait full name`, Coefficient_P_value, benchmark)

most_sig.df
table(most_sig.df$benchmark)

# Which metrics find cell types we think are relevant to their traits?
immune.cells <- c("blood_cell","classical_monocyte", "professional_antigen_presenting_cell",
                  "classical_monocyte", "lymphocyte", "natural_killer_cell",
                  unique(filter(df.ldsc_cts, tissue %in% c('Marrow', 'Spleen', 'Thymus'))$cell_type))

brain.cells <- unique(filter(df.ldsc_cts, tissue %in% c("Brain_Non-Myeloid"))$cell_type)

cardiovascular.cells <- setdiff(unique(filter(df.ldsc_cts, tissue %in% c("Heart"))$cell_type), 
                                c("professional_antigen_presenting_cell", "leukocyte", 'unknown_cell_type'))

lipid.cells <- setdiff(unique(filter(df.ldsc_cts, tissue %in% c("Liver"))$cell_type), 
                       c("B_cell", "natural_killer_cell"))

cells.and.categories <- list(Lipids = lipid.cells,
                             Cardiovascular = cardiovascular.cells,
                             Brain = brain.cells,
                             Immune = immune.cells)

# Immune
immune.relevant.df <- df.ldsc_cts %>% 
  filter(Coefficient_P_value < 0.05/115, `Trait category` == "Immune") %>% 
  mutate(contains_relevant_cells = case_when(cell_type %in% cells.and.categories$"Immune" ~ TRUE,
                                             TRUE ~ FALSE))

# Cardiovascular
cardio.relevant.df <- df.ldsc_cts %>% 
  filter(Coefficient_P_value < 0.05/115, `Trait category` == "Cardiovascular") %>%  
  mutate(contains_relevant_cells = case_when(cell_type %in% cells.and.categories$"Cardiovascular" ~ TRUE,
                                             TRUE ~ FALSE))

# Brain
brain.relevant.df <- df.ldsc_cts %>% 
  filter(Coefficient_P_value < 0.05/115, `Trait category` %in% c("Neuropsychiatric", "Behavioral")) %>% 
  mutate(contains_relevant_cells = case_when(cell_type %in% cells.and.categories$"Brain" ~ TRUE,
                                             # TRUE ~ FALSE))

# Lipids
lipids.relevant.df <- df.ldsc_cts %>% 
  filter(Coefficient_P_value < 0.05/115, `Trait category` == "Lipids") %>% 
  mutate(contains_relevant_cells = case_when(cell_type %in% cells.and.categories$"Lipids" ~ TRUE,
                                             TRUE ~ FALSE))

all.relevant.df <- bind_rows(immune.relevant.df, cardio.relevant.df, brain.relevant.df, lipids.relevant.df)

relevant.df <- all.relevant.df %>% 
  filter(contains_relevant_cells == TRUE) %>% 
  group_by(`Trait full name`,benchmark) %>% 
  summarise(n())

not.relevant.df <- all.relevant.df %>% 
  filter(contains_relevant_cells == FALSE) %>% 
  group_by(`Trait full name`, benchmark) %>% 
  summarise(n())

table(relevant.df$benchmark)
28 - table(not.relevant.df$benchmark)
