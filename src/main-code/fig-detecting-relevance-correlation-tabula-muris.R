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
library(GGally)

source(here("src/publication/lib-load_pub_lib_functions.R"))


setwd(here("src/benchmark"))

# ===================================================================================s== #
# ================== LOAD LDSC CTS RESULTS (single GWAS, multi metric) ================ #
# ===================================================================================== #

gwas <- "EA3_Lee2018"
dataset_prefix <- "tabula_muris"
dir.data <- here("out/CELLECT-LDSC-out/ESmu-variations")


### Loading BMI data
normal_files <- list.files(dir.data, pattern = sprintf("tabula_muris.*[0-9]{6}__%s.*.txt", gwas))
null_4_files <- list.files(dir.data, pattern = sprintf("*[.]4__%s.*.txt", gwas))
files <- c(normal_files, null_4_files)

data = tibble(File = files) %>%
  mutate(Data = lapply(file.path(dir.data, files), read_tsv)) %>%
  unnest(Data) %>% 
  separate(Name, c('benchmark', 'annotation'), '__') %>%  
  mutate(bench_label = case_when(benchmark == "tabula_muris.ESmu.binary.benchmark190722" ~ 'ESmu binary',
                                 benchmark == "tabula_muris.ESmu.continuous.benchmark190722" ~ 'ESmu\ncontinuous',
                                 benchmark == "tabula_muris.ESmu.continuous-squared.benchmark190722" ~ 'ESmu\ncontinuous^2',
                                 benchmark == "tabula_muris.ESmuQ1.continuous.190730" ~ 'ESmu quantile\nnormalised v1',
                                 benchmark == "tabula_muris.ESmuQ2.continuous.190730" ~ 'ESmu quantile\nnormalised v2',
                                 benchmark == "tabula_muris.ESmu-NULL1.binary.benchmark190722.4" ~ 'ESmu binary\npermuted',
                                 benchmark == "tabula_muris.ESmu-NULL1.continuous.benchmark190722.4" ~ 'ESmu\ncontinuous\npermuted',
                                 benchmark == "tabula_muris.ESmu-NULL2.continuous.benchmark190722.4" ~ 'ESmu\ncontinuous\n+ noise')) %>% 
  select(-Coefficient, -Coefficient_std_error, -File, -benchmark) %>%
  spread(key = bench_label, value = Coefficient_P_value)

data


# ======================================================================= #
# ========================= LOAD CELL-TYPE METADATA ===================== #
# ======================================================================= #

### Read annotation metadata
df.metadata <- get_metadata(dataset_prefix="tabula_muris")
df.metadata

### Add meta data
df.ldsc_cts <- data %>% left_join(df.metadata, by="annotation") # add meta data
df.ldsc_cts 

df.ldsc_cts <- df.ldsc_cts %>% 
  mutate(relevant_tissue = case_when(tissue == "Brain_Non-Myeloid" ~ "Brain",
                                     tissue != "Brain_Non-Myeloid" ~ "Other"))

  

p <- ggpairs(df.ldsc_cts, columns = 2:9, ggplot2::aes(colour=relevant_tissue),diag = "blank",
             upper = list(continuous = wrap("cor", size = 3.2)),
             lower = list(continuous = wrap("points", alpha = 0.3,size=1))) + 
  theme_minimal() + 
  theme(strip.text = element_text(size = 12.5),
        axis.text.x = element_text(angle = 90))
p
### Save
file.out <- sprintf("figs/fig-detect-relevance-correlation.tabula-muris.190830.png")
ggsave(p, filename=file.out, width=12, height=12)
