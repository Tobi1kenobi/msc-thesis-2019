############### SYNOPSIS ####################
# IDEA:
# Repeat the heatmap for Tabula Muris Senis. Maybe make it 
# comparative so that we can see the difference between TM and
# TMS for each cell type. However maybe not as cell type classification
# isn't identical

# AIMS:
# - See which cell types are prioritised in aged mice
# - See if diseases we know are associated with ageing have stronger signal than in TM
#

# PLOT:
# A heatmap to show 35 GWAS and as many cell types as there are. But maybe just add this to supplementary
# the really cool story would be a comparative heatmap with TM - to see if aged mice show differences
# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library


library(RColorBrewer)

library(stringr)
library(colorspace)

source(here("src/publication/lib-load_pub_lib_functions.R"))


setwd(here("src/benchmark"))

# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

### Read LDSC results
dir.data <- here("out/CELLECT-LDSC-out/tabula-muris-senis-24m")
files <- list.files(dir.data,pattern=".*ESmu.*.txt")
files

data = tibble(File = files) %>%
  extract(File, "Trait identifier", "__([^.]*)", remove = FALSE) %>%
  mutate(Data = lapply(file.path(dir.data, files), read_tsv)) %>%
  unnest(Data) %>% 
  separate(Name, c('benchmark', 'annotation'), '__')%>% 
  mutate(tms_adj_p_val =  pmin(Coefficient_P_value * 131, 1)) %>% 
  mutate(tms_sig_threshold = case_when(tms_adj_p_val < 0.001 ~ '***',
                                      tms_adj_p_val < 0.01 ~ '**',
                                      tms_adj_p_val < 0.05 ~ '*',
                                      tms_adj_p_val >= 0.05 ~ ''))

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
df.metadata.cells <- read_csv('/scratch/data-for_fast_access/pub-others/tabula_muris_senis/24mFACS_metadata-190722.csv',guess_max = 100000)
df.metadata <- df.metadata.cells %>% 
  select(-c(cells, FACS.selection, batch, cell, cell_ontology_id, cellid, method, plate, sex, well, n_genes,
            n_counts, louvain, leiden, cluster_number, cluster_names, free_annotation, mouse.id, subtissue)) %>% 
  distinct() %>% 
  rename(annotation = tissue_cell_ontology_class)

df.metadata
  

### Add meta data
df.ldsc_cts <- data.only.35 %>% left_join(df.metadata, by="annotation") # add meta data
df.ldsc_cts

### Add fdr_significant flag (within GWAS)
df.ldsc_cts.tmp.summary <- df.ldsc_cts %>% group_by(File) %>% summarise(n_obs_File=n())
df.ldsc_cts <- left_join(df.ldsc_cts, df.ldsc_cts.tmp.summary, by="File")
df.ldsc_cts <- df.ldsc_cts %>% mutate(fdr_significant = case_when(Coefficient_P_value <= 0.001/n_obs_File ~ '***',
                                                                  Coefficient_P_value <= 0.01/n_obs_File ~ '**',
                                                                  Coefficient_P_value <= 0.05/n_obs_File ~ '*',
                                                                  TRUE ~ ''))
df.ldsc_cts


# ======================================================================= #
# ======================== PRE-PROCESS ANNOTATIONS ====================== #
# ======================================================================= #

# Self-defined label formatting function
tm_annotation_formatter <- function(annotation_tissue_celltype) {
  # INPUT: annotation_tissue_celltype, string e.g. 'Marrow.Slamf1-negative_multipotent_progenitor_cell'
  # USAGE 1: df %>% mutate
  # USAGE 2: ggplot() + scale_y_discrete(label=tm_annotation_formatter)
  label <- stringr::str_split_fixed(annotation_tissue_celltype, pattern="\\.", n=Inf)[,2]
  label <- stringr::str_replace_all(label, pattern="-", " - ") # add space to any hyphens
  label <- stringr::str_replace_all(label, pattern="_", " ") # convert _ to space
  # label <- stringr::str_to_sentence(label) # title case
  substr(label, 1, 1) <- toupper(substr(label, 1, 1)) # capitalize first character
  return(label)
}

### annotation: create clean names
df.plot.tax_order <- df.ldsc_cts %>% mutate(annotation_label_fmt = tm_annotation_formatter(annotation))
df.plot.tax_order

### tissue: rename + remove some text because they contain too few cell-types
df.plot.tax_order %>% count(tissue, sort=T) # ---> potentially filter : n_annotations_in_tax >= 5
df.plot.tax_order <- df.plot.tax_order %>% mutate(tissue = case_when(
  tissue == "Brain_Myeloid" ~ "Brain",
  tissue == "Brain_Non-Myeloid" ~ "Brain",
  TRUE ~ stringr::str_replace_all(tissue, pattern="_", replacement="\n") # TRUE ~ as.character(tissue)
)
)
df.plot.tax_order %>% mutate(tissue = factor(tissue))

# ====================================================================== #
# ======================== PLOT HEATMAP ================================ #
# ====================================================================== #

multi_ES_metric_heatmap_plot <- function(df.ldsc_cts.multi_gwas) {
  
  ### Plot
  p.multi_gwas <- ggplot(df.ldsc_cts.multi_gwas, aes(y=`Trait full name`, x=annotation_label_fmt, fill=-log10(Coefficient_P_value))) + 
    facet_grid(~ tissue, scales='free_x', space="free_x") +
    geom_tile() + 
    geom_text(aes(label=fdr_significant), color="white", hjust=0.5, vjust=0.75, size = 2.2) + # add asterisk if fdr significant
    colorspace::scale_fill_continuous_sequential(palette="Greens 2", rev=TRUE, limits=c(0,-log10(0.05/131)), oob=scales::squish, na.value = "white") + # "Blues 2","Blues 3","Purples 3"
    # ^ oob: Function that handles limits outside of the scale limits (out of bounds). scales::squish "squishes" values into the range specified by limits
    # labs(fill=expression(-log[10](P[S-LDSC]))) +
    theme_minimal() + # set this first
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) + # remove grid
    theme(axis.text.x=element_text(angle=50, hjust=1, size=7)) +
    theme(legend.position="bottom") +
    facet_grid(`Trait category` ~ tissue, scales = 'free', space = 'free') +
    theme(panel.spacing = unit(0.1, "lines")) + # Decrease space between plots 
    theme(strip.text.x=element_text(size = 11,angle = 90, face = "bold"), # Controls facet label appearance
          strip.text.y=element_text(size = 8, angle = 0, face = "bold")) + 
    ylab('Complex trait') +
    xlab('Cell type') +
    labs(fill = expression(-log[10](P[CELLECT-LDSC])))
  
  ### Add spacing to heatmap
  # REF: https://stackoverflow.com/questions/40156061/white-space-between-tiles-in-heatplot-ggplot
  ### Approaches
  # 1: use facet_wrap ---> BEST
  # 2: add "NULL" values (+ order factor levels)
  # 3: patchwork / cowplot::plot_grid
  # 4: add separator line : http://www.roymfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r/ ---> only heatmap2
  # p.multi_metric <- p.multi_metric +
  #   facet_grid(~space_block, scales='free_x', space="free_x") +
  #   theme(strip.background = element_blank(), strip.text.x = element_blank()) # remov facet labels completely
  return(p.multi_gwas)
}

p <- multi_ES_metric_heatmap_plot(df.plot.tax_order)
p



# ===================================================================== #
# ======================== COMPARISON HEATMAP ========================= #
# ===================================================================== #

### Read LDSC results
dir.tm.data <- here("out/CELLECT-LDSC-out/tabula-muris")
tm.files <- list.files(dir.tm.data,pattern=".*ESmu.continuous[.].*.txt")
tm.files

tm.data = tibble(File = tm.files) %>%
  extract(File, "Trait identifier", "__([^.]*)", remove = FALSE) %>%
  mutate(Data = lapply(file.path(dir.tm.data, tm.files), read_tsv)) %>%
  unnest(Data) %>% 
  separate(Name, c('benchmark', 'annotation'), '__') %>% 
  mutate(tm_adj_p_val = pmin(Coefficient_P_value * 115, 1)) %>% 
  select(`Trait identifier`, annotation, tm_adj_p_val) %>% 
  mutate(tm_sig_threshold = case_when(tm_adj_p_val < 0.001 ~ '***',
                                      tm_adj_p_val < 0.01 ~ '**',
                                      tm_adj_p_val < 0.05 ~ '*',
                                      tm_adj_p_val >= 0.05 ~ ''))

compare.plot.df <-  df.plot.tax_order %>% left_join(tm.data) %>% 
  mutate(more_significant = case_when(str_length(tms_sig_threshold) < str_length(tm_sig_threshold) ~ 'Tabula Muris',
                                          str_length(tms_sig_threshold) > str_length(tm_sig_threshold) ~ 'Tabula Muris Senis',
                                          TRUE ~ 'Neither'))

compare.table.df <- compare.plot.df %>% 
  filter(more_significant == FALSE, tm_adj_p_val < 0.05) %>% 
  select(tms_adj_p_val, tm_adj_p_val, annotation, `Trait full name`, `Trait category`) 

multi_ES_metric_comparative_heatmap_plot <- function(df.ldsc_cts.multi_gwas) {
  
  ### Plot
  p.multi_gwas <- ggplot(df.ldsc_cts.multi_gwas, aes(y=`Trait full name`, x=annotation_label_fmt, fill=-log10(Coefficient_P_value))) + 
    facet_grid(~ tissue, scales='free_x', space="free_x") +
    geom_tile() + 
    geom_text(aes(label=fdr_significant, color = more_significant), hjust=0.5, vjust=0.75, size = 2.3) + # add asterisk if fdr significant
    colorspace::scale_fill_continuous_sequential(palette="Grays", rev=TRUE, limits=c(0,-log10(0.05/131)), oob=scales::squish, na.value = "white") + # "Blues 2","Blues 3","Purples 3"
    # ^ oob: Function that handles limits outside of the scale limits (out of bounds). scales::squish "squishes" values into the range specified by limits
    # labs(fill=expression(-log[10](P[S-LDSC]))) +
    theme_minimal() + # set this first
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) + # remove grid
    theme(axis.text.x=element_text(angle=50, hjust=1, size=7)) +
    theme(legend.position="bottom") +
    scale_color_manual(values=c("darkgrey", "red", "yellow")) +
    facet_grid(`Trait category` ~ tissue, scales = 'free', space = 'free') +
    theme(panel.spacing = unit(0.1, "lines")) + # Decrease space between plots 
    theme(strip.text.x=element_text(size = 11,angle = 90, face = "bold"), # Controls facet label appearance
          strip.text.y=element_text(size = 8, angle = 0, face = "bold")) + 
    ylab('Complex trait') +
    xlab('Cell type') +
    labs(fill = expression(-log[10](P[CELLECT-LDSC])),
         color = 'Most significant') +
    guides(colour = guide_legend(override.aes = list(size=10, label='*')))
  
  return(p.multi_gwas)
}


p.compare <- multi_ES_metric_comparative_heatmap_plot(compare.plot.df)
p.compare

heatmap.file.out <- "figs/fig-compare-all-gwas-tabula-muris-senis.190827.png"
# ggsave(p.compare, filename=heatmap.file.out, width=19, height=8)

# Older version of above plot - only show binary significance
p.compare.2 <- ggplot(compare.plot.df, aes(y=`Trait full name`, x=annotation_label_fmt, fill=more_significant)) +
  geom_tile()+
  facet_grid(`Trait category` ~ tissue, scales = 'free', space = 'free') +
  theme_minimal() + # set this first
  theme(axis.text.x=element_text(angle=50, hjust=1, size=7)) +
  theme(panel.spacing = unit(0.1, "lines")) + # Decrease space between plots
  theme(strip.text.x=element_text(size = 11,angle = 90, face = "bold"), # Controls facet label appearance
        strip.text.y=element_text(size = 8, angle = 0, face = "bold"),
        panel.grid = element_line()) +
  theme(legend.position="bottom") +
  scale_fill_manual(values=c("darkgrey", "red", "yellow")) +
  ylab('Complex trait') +
  xlab('Cell type') +
  labs(fill = 'More significant')

p.compare.2

heatmap.file.out2 <- "figs/fig-compare-all-gwas-tabula-muris-senis-only-diff.190827.png"
# ggsave(p.compare.2, filename=heatmap.file.out2, width=19, height=8)

