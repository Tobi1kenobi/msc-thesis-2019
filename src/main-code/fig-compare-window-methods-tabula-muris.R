############### SYNOPSIS ###################
# Tabula Muris 100kb windows vs SNPsnap LD, R^2>0.5 windows

# IDEA:
# See if SNPsnap windows are an improvement on arbitrarily assigning
# SNPs to genes in a 100kb region either side of the gene. Intuitively
# it should improve the results but lots of SNPs are lost and this may
# impact the results.

# AIMS:
# - To observe effect on magnitude of signficance
# - To see if some of the noise is lost e.g. pancreas, non-neuronal cells in brain
# - To try and observe any improvement over 100kb window method


# PLOT:
# Multi GWAS definitely but scatter plot or something else? I only need to compare
# 2 approaches so facet grid not necessary. Maybe back-to-back ordered scatter plots
# normalised to lowest p-val with facet_wrap and two or three rows
# Can trick scatter plot into doing back-to-back by using log10(SNPsnap p-vals) and
# -log10(100kb p-vals) so no nesting facet wrap.
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

# ==================================================================================== #
# ================== LOAD LDSC CTS RESULTS (multi GWAS, SNP vs 100kb window) ========= #
# ==================================================================================== #

### Read brain and no brain LDSC results
SNP.dir.data <- here("out/CELLECT-LDSC-out/SNPsnap-windows")
SNP.files <- list.files(SNP.dir.data,pattern="*.txt") # only has ESmu continuous in the directory
SNP.files

SNP.data = tibble(File = SNP.files) %>%
  extract(File, "Trait identifier", "__([^.]*)", remove = FALSE) %>%
  mutate(Data = lapply(file.path(SNP.dir.data, SNP.files), read_tsv)) %>%
  unnest(Data) %>% 
  separate(Name, c('benchmark_date', 'annotation'), '__') %>% 
  extract(benchmark_date, "benchmark", "([^.]*[.][^.]*[.][^.]*)", remove = TRUE) %>% 
  group_by(`Trait identifier`) %>% 
  mutate(rank = dense_rank(desc(Coefficient_P_value))) %>% 
  ungroup %>% 
  mutate(window = 'Linkage disequilibrium')

SNP.data

### Read whole body results
tm.dir.data <- here("out/CELLECT-LDSC-out/tabula-muris")
tm.files <- list.files(tm.dir.data,pattern=".ESmu.continuous[.].*.txt") # only include ESmu continuous files
tm.files

tm.data = tibble(File = tm.files) %>%
  extract(File, "Trait identifier", "__([^.]*)", remove = FALSE) %>%
  mutate(Data = lapply(file.path(tm.dir.data, tm.files), read_tsv)) %>%
  unnest(Data) %>% 
  separate(Name, c('benchmark_date', 'annotation'), '__') %>% 
  extract(benchmark_date, "benchmark", "([^.]*[.][^.]*[.][^.]*)", remove = TRUE) %>% 
  group_by(`Trait identifier`) %>% 
  mutate(rank = dense_rank(desc(Coefficient_P_value))) %>% 
  ungroup %>% 
  mutate(window = '100kb')

## Join both dataframes into one


data <- bind_rows(SNP.data, tm.data)

### Which GWAS am I interested in?
GWAS.info <- read_csv('190713_GWAS_Timshel (tobi edit).csv')
GWAS.info


# Would be 39 but I remove 4 traits that violate LDSC assumptions too
data.only.35 <- data %>%
  left_join(GWAS.info, by="Trait identifier") %>% 
  filter(FLAG_final & GWAS.trait.representative & is.na(`GWAS violates LDSC assumptions`))


# gwas.interested.group <- c("HEIGHT_UKBB_Loh2018", "EA3_Lee2018", "SCZ_Pardinas2018")
# 
# gwas.bmi.paper.group <- c("HEIGHT_UKBB_Loh2018", "EA3_Lee2018", "SCZ_Pardinas2018",
#                           "BMI_UKBB_Loh2018", "INSOMNIA_Jansen2018",
#                           "INTELLIGENCE_Savage2018", "LIPIDS_LDL_Teslovich2010",
#                           "WHRadjBMI_UKBB_Loh2018")
# gwas.main = c(gwas.bmi.paper.group, "DEPRESSION_Nagel2018", "ASD_iPSYCH_PGC_Grove2018",
#               "BIP_PGC2018", "CARDIOVASCULAR_UKBB_Loh2018", "CROHNS_Jostins2012",
#               "IBD_Jostins2012", "MDD_Howard2019", "NEUROTICISM_Nagel2018", "RB_Linner_2019",
#               "WORRY_Nagel2018" )
# 
# data.only.interested.gwas <- data.only.35 %>% 
#   filter(`Trait identifier` %in% gwas.main)
# 
# data.only.interested.gwas

# ======================================================================= #
# ========================= LOAD CELL-TYPE METADATA ===================== #
# ======================================================================= #

### Read annotation metadata
df.metadata <- get_metadata(dataset_prefix="tabula_muris")
df.metadata

### Add meta data
df.ldsc_cts <- data.only.35 %>% left_join(df.metadata, by="annotation") # add meta data
df.ldsc_cts

### Order annotations
df.ldsc_cts <- df.ldsc_cts %>% mutate(annotation=factor(annotation, levels=unique(annotation[order(tissue, cell_type)]))) %>% 
  mutate(tissue = case_when(
  tissue == "Brain_Myeloid" ~ "Brain",
  tissue == "Brain_Non-Myeloid" ~ "Brain",
  TRUE ~ stringr::str_replace_all(tissue, pattern="_", replacement="\n")))

# Self-defined label formatting function
tm_annotation_formatter <- function(annotation_tissue_celltype) {
  # INPUT: annotation_tissue_celltype, string e.g. 'Marrow.Slamf1-negative_multipotent_progenitor_cell'
  # USAGE 1: df %>% mutate
  # USAGE 2: ggplot() + scale_y_discrete(label=tm_annotation_formatter)
  label <- stringr::str_split_fixed(annotation_tissue_celltype, pattern="\\.", n=Inf)
  tissues <- label[,1]
  cells <- label[,2]
  cells <- stringr::str_replace_all(cells, pattern="-", " - ") # add space to any hyphens
  cells <- stringr::str_replace_all(cells, pattern="_", " ") # convert _ to space
  tissues <- stringr::str_replace_all(tissues, pattern="_", " ") # convert _ to space
  # label <- stringr::str_to_sentence(label) # title case
  substr(cells, 1, 1) <- toupper(substr(cells, 1, 1)) # capitalize first character
  tissues <- toupper(tissues) # capitalize all characters
  out_label <- paste(tissues, cells)
  return(out_label)
}
window_comparison_plotter <- function(gwas){
  ### annotation: create clean names and then order cell types for plotting
  df.plot.tax_order <- df.ldsc_cts %>% 
    filter(`Trait identifier` == gwas) %>% 
    mutate(annotation_label_fmt = tm_annotation_formatter(annotation))
  
  
  df.plot.tax_order <- df.plot.tax_order %>%  
    mutate(annotation_label_fmt = factor(annotation_label_fmt, # IMPORTANT: Order factor by SNPsnap ranks (not sure why it does SNPsnap not 100kb here)
                                         levels=annotation_label_fmt[order(filter(df.plot.tax_order,
                                                                                  window == '100kb')$rank)])) %>% 
    filter(annotation_label_fmt %in% filter(df.plot.tax_order,window == 'Linkage disequilibrium', rank > 85)$annotation_label_fmt) # Filter all cells that aren't in the top 30 of SNPsnap significant p vals
  
  df.plot.tax_order$`Mapping method` <- factor(df.plot.tax_order$window)
  # ========================================================================== #
  # ============================= PLOT ======================================= #
  # ========================================================================== #
  
  # p <- ggplot(df.plot.tax_order) +
  #   geom_bar(aes(y = log10(Coefficient_P_value), x = rank, fill = tissue, label = annotation_label_fmt),
  #            stat = 'identity', colour='black', size = 0.2, alpha = 0.5,position = "dodge") +
  #   coord_flip() +
  #   facet_wrap(~`Trait full name`, scales = 'free')
  
  p <- ggplot(df.plot.tax_order) +
    geom_point(aes(y = -log10(Coefficient_P_value), x = annotation_label_fmt, colour =  `Mapping method`), size = 2, alpha = 0.5) +
    coord_flip() +
    facet_wrap(~`Trait full name`, scales = 'free') +
    geom_hline(yintercept = -log10(0.05/115), colour = 'black', alpha = 0.5, linetype = "dashed") +
    theme_minimal() +
    theme(axis.text=element_text(size=6),
          strip.text = element_text(size=12, face='bold'),
          legend.position = 'bottom') +
    ylab(expression(-log[10](P[CELLECT-LDSC]))) +
    xlab('Top 30 cell types\nranked by LD mapping')
  
  file.out <- sprintf("figs/fig-SNPsnap-vs-100kb.%s.tabula-muris.190831.png", gwas)
  # ggsave(file.out, p, width = 7, height = 4)
  return(p)
}

# plots <- lapply(gwas.main, window_comparison_plotter)
# plots


p_val.compare <- df.ldsc_cts %>% group_by(`Trait identifier`, window) %>%
  filter(Coefficient_P_value == min(Coefficient_P_value)) %>%
  ungroup %>% 
  select(`Trait full name`, window, Coefficient_P_value,tissue_celltype) %>%
  arrange(`Trait full name`) %>%
  spread(window, Coefficient_P_value,)

# write_csv(p_val.compare,path = 'LD-mapping-vs-physical-mapping.csv', col_names = TRUE,na = ' ')
