############### SYNOPSIS ###################
# Tabula Muris Whole vs Tabula Muris Brain + No-Brain

# IDEA:
# Compare CELLECT when the input to CELLEX is limited to only a certain region.
# In this case I have ran CELLECT with Tabula Muris, Tabula Muris[Brain_Non-Microglia],
# and Tabula_muris[-Brain_Non-Microglia] CELLEX data.
# Because lots of our analysis is performed on brain only datasets we need to
# know if the specificity method works as well when it can't "see" the whole body

# AIMS:
# - To check if order of significant cell types is preserved
# - To observe effect on magnitude of signficance
# - To confirm brain traits still have neuron emerge as significant


# PLOT:
# Will plot for maybe 2 or 3 GWAS and ESmu continuous, a back-to-back bar chart
# with Tabula Muris on left and Tabula Muris[Brain] + Tabula Muris[-Brain] 
# concatenated on the other side. In the middle maybe have a third "plot" which
# is a line connecting the ranking of the two plots
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
# ================== LOAD LDSC CTS RESULTS (multi GWAS, TM vs TM (concat)) =========== #
# ==================================================================================== #

### Read brain and no brain LDSC results
brain.dir.data <- here("out/CELLECT-LDSC-out/subset-brain")
brain.files <- list.files(brain.dir.data,pattern="*.txt") # only has ESmu continuous in the directory
brain.files

brain.data = tibble(File = brain.files) %>%
  extract(File, "Trait identifier", "__([^.]*)", remove = FALSE) %>%
  mutate(Data = lapply(file.path(brain.dir.data, brain.files), read_tsv)) %>%
  unnest(Data) %>% 
  separate(Name, c('benchmark_date', 'annotation'), '__') %>% 
  extract(benchmark_date, "benchmark", "([^.]*[.][^.]*[.][^.]*)", remove = TRUE) %>% 
  group_by(benchmark,`Trait identifier`) %>% 
  mutate(adj_p_val = Coefficient_P_value*n()) %>% 
  ungroup %>% 
  group_by(`Trait identifier`) %>% 
  mutate(Rank = dense_rank(desc(adj_p_val))) %>% 
  ungroup %>% 
  mutate(data_split = 'Brain + No brain')
  

brain.data

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
  mutate(adj_p_val = Coefficient_P_value*n()) %>% 
  mutate(Rank = dense_rank(desc(adj_p_val))) %>% 
  ungroup %>% 
  mutate(data_split = 'Whole body')

## Join both dataframes into one

data <- bind_rows(brain.data, tm.data)

### Which GWAS am I interested in?
GWAS.info <- read_csv('190713_GWAS_Timshel (tobi edit).csv')
GWAS.info


# Would be 39 but I remove 4 traits that violate LDSC assumptions too
data.only.35 <- data %>%
  left_join(GWAS.info, by="Trait identifier") %>% 
  filter(FLAG_final & GWAS.trait.representative & is.na(`GWAS violates LDSC assumptions`))

# Half the number of GWAS so plot isn't too cluttered - split into "main group" and "secondary group"
# Removed RA and MS because there are too many immune cells in different tissues that appear significant
gwas.interested.group <- c("HEIGHT_UKBB_Loh2018", "EA3_Lee2018", "SCZ_Pardinas2018")

data.only.interested.gwas <- data.only.35 %>% 
  filter(`Trait identifier` %in% gwas.interested.group)

# ======================================================================= #
# ========================= LOAD CELL-TYPE METADATA ===================== #
# ======================================================================= #

### Read annotation metadata
df.metadata <- get_metadata(dataset_prefix="tabula_muris")
df.metadata

### Add meta data
df.ldsc_cts <- data.only.interested.gwas %>% left_join(df.metadata, by="annotation") # add meta data
df.ldsc_cts

### Order annotations
df.ldsc_cts <- df.ldsc_cts %>% mutate(annotation=factor(annotation, levels=unique(annotation[order(tissue, cell_type)]))) # order annotations by their taxonomy (from MB Fig1c), and secondary order by their annotation name

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
df.plot.tax_order <- df.ldsc_cts %>% 
  mutate(annotation_label_fmt = tm_annotation_formatter(annotation)) %>% 
  mutate(significance_star = case_when(adj_p_val < 0.001 ~ "***",
                                       adj_p_val < 0.01 ~ "**",
                                       adj_p_val < 0.05 ~ "*",
                                 adj_p_val >= 0.05 ~ ""),
         `Tabula Muris subset` =  case_when(benchmark == "tabula_muris-no_neuron.ESmu.continuous" ~ "Body without brain",
                                            benchmark == "tabula_muris.ESmu.continuous" ~ "Whole body",
                                            benchmark == "tabula_muris-only_neuron.ESmu.continuous" ~ "Only brain"))


# ======================================================================= #
# ========================= PLOT RANK CHART ============================= #
# ======================================================================= #

# data.one.gwas <- filter(df.plot.tax_order, `Trait identifier` == 'HEIGHT_UKBB_Loh2018')

# Adding two empty levels either side to make space for the geom_text labels in the plot
df.plot.tax_order$data_split <- factor(df.plot.tax_order$data_split,levels=c('','Brain + No brain','Whole body',' '), ordered = TRUE)

p <- ggplot(df.plot.tax_order,
            aes(x = data_split, y = Rank, label = annotation_label_fmt, colour = `Tabula Muris subset`, group=tissue_celltype)) +
  geom_line(size=0.3, alpha = 0.5) + 
  geom_text(data=filter(df.plot.tax_order, data_split == 'Brain + No brain'),
            size = 2, hjust=1) +
  geom_text(data=filter(df.plot.tax_order, data_split == 'Whole body'),
            size = 2, hjust=0) +
  geom_text(data=filter(df.plot.tax_order, data_split == 'Whole body'),aes(label = significance_star),
            size = 3, hjust=1, fontface = 'bold') +
  geom_text(data=filter(df.plot.tax_order, data_split == 'Brain + No brain'),aes(label = significance_star),
            size = 3, hjust=0, fontface = 'bold') +
  facet_wrap(~ `Trait full name`) + 
  scale_color_brewer(palette="Dark2") +
  theme_minimal() +
  theme(axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=10, face='bold'),
        panel.grid = element_blank(),
        strip.text = element_text(size=18, face='bold')) + 
  scale_x_discrete(expand = c(0.01, 0.01),drop=FALSE) + 
  xlab('')

p

### Save
file.out <- sprintf("figs/fig-subset-data.bump-chart.tabula-muris.190901.png")
ggsave(p, filename=file.out, width=16.5, height=7.5)


# ======================================================================= #
# ========================= CALCULATE CORRELATION ======================= #
# ======================================================================= #

calc_corrs <- function(gwas_id){
  only_brain <- filter(df.plot.tax_order, `Tabula Muris subset` == 'Only brain',
                       `Trait identifier` == gwas_id) %>% 
    arrange(tissue_celltype)
  brain_cells <- unique(only_brain$tissue_celltype)
  whole_body_for_brain <- filter(df.plot.tax_order, `Tabula Muris subset` == 'Whole body',
                       `Trait identifier` == gwas_id,
                       tissue_celltype %in% brain_cells)%>% 
    arrange(tissue_celltype)
  
  no_brain <- filter(df.plot.tax_order, `Tabula Muris subset` == 'Body without brain',
                       `Trait identifier` == gwas_id)%>% 
    arrange(tissue_celltype)
    
  no_brain_cells <- unique(no_brain$tissue_celltype)
  whole_body_no_brain <- filter(df.plot.tax_order, `Tabula Muris subset` == 'Whole body',
                                 `Trait identifier` == gwas_id,
                                 tissue_celltype %in% no_brain_cells)%>% 
    arrange(tissue_celltype)

  brain_spear <- cor(only_brain$adj_p_val, whole_body_for_brain$adj_p_val, method = "spearman")
  brain_pears  <- cor(only_brain$adj_p_val, whole_body_for_brain$adj_p_val, method = "pearson")

  no_brain_spear <- cor(no_brain$adj_p_val, whole_body_no_brain$adj_p_val, method = "spearman")
  no_brain_pears <- cor(no_brain$adj_p_val, whole_body_no_brain$adj_p_val, method = "pearson")
  
  print(sprintf("Brain Spearman correlation for %s is: %s", gwas_id, brain_spear))
  print(sprintf("Brain Pearson correlation for %s is: %s", gwas_id, brain_pears))

  print(sprintf("No brain Spearman correlation for %s is: %s", gwas_id, no_brain_spear))
  print(sprintf("No brain Pearson correlation for %s is: %s", gwas_id, no_brain_pears))
  return(NA)
}

lapply(unique(df.ldsc_cts$`Trait identifier`),calc_corrs)

# ======================================================================= #
# ========================= MOST SIG CELL TYPES ========================= #
# ======================================================================= #

df.plot.tax_order %>% 
  group_by(`Trait identifier`,data_split) %>% 
  filter(adj_p_val == min(adj_p_val)) %>% 
  select(adj_p_val, `Trait identifier`,data_split,`Tabula Muris subset`, tissue_celltype, Coefficient_P_value) %>% 
  View()
