############### SYNOPSIS ###################
# Make Tabula Muris heatmap for 39 GWAS

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
files <- list.files(dir.data,pattern=".*ESmu.continuous[.].*.txt")
files

data = tibble(File = files) %>%
  extract(File, "Trait identifier", "__([^.]*)", remove = FALSE) %>%
  mutate(Data = lapply(file.path(dir.data, files), read_tsv)) %>%
  unnest(Data) %>% 
  separate(Name, c('benchmark', 'annotation'), '__')

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

### Order annotations
df.ldsc_cts <- df.ldsc_cts %>% mutate(annotation=factor(annotation, levels=unique(annotation[order(tissue, cell_type)]))) # order annotations by their taxonomy (from MB Fig1c), and secondary order by their annotation name
# ^ if df.ldsc_cts only contain 1 gwas, then you don't need the levels=unique().
df.ldsc_cts$annotation %>% levels()

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

multi_ES_metric_heatmap_plot <- function(df.ldsc_cts.multi_gwas) {
  
  ### Plot
  p.multi_gwas <- ggplot(df.ldsc_cts.multi_gwas, aes(y=`Trait full name`, x=annotation_label_fmt, fill=-log10(Coefficient_P_value))) + 
    facet_grid(~ tissue, scales='free_x', space="free_x") +
    geom_tile() + 
    geom_text(aes(label=fdr_significant), color="white", hjust=0.5, vjust=0.75, size = 2.2) + # add asterisk if fdr significant

    colorspace::scale_fill_continuous_sequential(palette="Reds 2", rev=TRUE, limits=c(0,3.361728), oob=scales::squish, na.value = "white") + # "Blues 2","Blues 3","Purples 3"
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
### Save
file.out <- sprintf("figs/fig-compare-all-gwas-heatmap.tabula-muris.190822.png")
ggsave(p, filename=file.out, width=18, height=8)
