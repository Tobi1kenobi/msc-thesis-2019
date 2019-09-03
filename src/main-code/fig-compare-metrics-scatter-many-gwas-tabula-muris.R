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

'%nin%' <- Negate('%in%')

# ==================================================================================== #
# ================== LOAD LDSC CTS RESULTS (multi GWAS, multi metric) ================ #
# ==================================================================================== #

### Read LDSC results
dir.data <- here("out/CELLECT-LDSC-out/ESmu-variations")
ESmu.variation.files <- list.files(dir.data,pattern=".*[0-9]{6}__.*.txt")
null.4.files <- list.files(dir.data,pattern="4__.*.txt")
files <- c(ESmu.variation.files, null.4.files)

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

### Order annotations
df.ldsc_cts <- df.ldsc_cts %>% mutate(annotation=factor(annotation, levels=unique(annotation[order(tissue, cell_type)]))) # order annotations by their taxonomy (from MB Fig1c), and secondary order by their annotation name
# ^ if df.ldsc_cts only contain 1 gwas, then you don't need the levels=unique().
df.ldsc_cts$annotation %>% levels()

# ======================================================================= #
# ======================== PRE-PROCESS ANNOTATIONS ====================== #
# ======================================================================= #

# Function to bbreviate annotations
word_shortener <- function(in_string){
  words = strsplit(in_string,' ')
  new_words <- lapply(words, abbreviate,named = FALSE, minlength = 6)
  return(paste(unlist(new_words), collapse = ' '))
}

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
  label <- unlist(lapply(label, word_shortener))
  return(label)
}

### annotation: create clean names
df.plot.tax_order <- df.ldsc_cts %>% 
  mutate(annotation_label_fmt = tm_annotation_formatter(annotation)) 

# Add a significant boolean for Bonf significant AND only cell type that is significant
df.plot.tax_order <- df.plot.tax_order %>% 
  group_by(`Trait identifier`, benchmark, cell_type) %>% 
    mutate(significant_show = case_when(
      Coefficient_P_value < 0.001/115 & min(Coefficient_P_value) == Coefficient_P_value ~ 6,
      Coefficient_P_value < 0.01/115 & min(Coefficient_P_value) == Coefficient_P_value ~ 5,
      Coefficient_P_value < 0.05/115 & min(Coefficient_P_value) == Coefficient_P_value ~ 4,
      Coefficient_P_value < 0.001/115 ~ 3,
      Coefficient_P_value < 0.01/115  ~ 2,
      Coefficient_P_value < 0.05/115  ~ 1,
      TRUE ~ 0 )) %>% 
  ungroup()
  

### rename + remove some text to ease visualisation
df.plot.tax_order %>% count(tissue, sort=T) # ---> potentially filter : n_annotations_in_tax >= 5
df.plot.tax_order <- df.plot.tax_order %>% mutate(tissue = case_when(
  tissue == "Brain_Myeloid" ~ "Brain",
  tissue == "Brain_Non-Myeloid" ~ "Brain",
  tissue %in% c("Tongue", "Trachea", "Bladder", "Mammary_Gland", "Lung", "Fat", "Spleen") ~ "Other",
  TRUE ~ stringr::str_replace_all(tissue, pattern="_", replacement="\n"))) %>% 
  mutate(`Trait full name` = stringr::str_replace_all(df.plot.tax_order$`Trait full name`,pattern = '(?<=[^.]) ', replace='\n')) %>% 
  mutate(bench_label = case_when(benchmark == "tabula_muris.ESmu.binary" ~ 'ESmu binary',
                                 benchmark == "tabula_muris.ESmu.continuous" ~ 'ESmu\ncontinuous',
                                 benchmark == "tabula_muris.ESmu.continuous-squared" ~ 'ESmu\ncontinuous^2',
                                 benchmark == "tabula_muris.ESmuQ1.continuous" ~ 'ESmu quantile\nnormalised v1',
                                 benchmark == "tabula_muris.ESmuQ2.continuous" ~ 'ESmu quantile\nnormalised v2',
                                 benchmark == "tabula_muris.ESmu-NULL1.binary" ~ 'ESmu binary\npermuted',
                                 benchmark == "tabula_muris.ESmu-NULL1.continuous" ~ 'ESmu\ncontinuous\npermuted',
                                 benchmark == "tabula_muris.ESmu-NULL2.continuous" ~ 'ESmu\ncontinuous\n+ noise')) %>% 
  rename(Tissue = tissue)

# Couldn't get actual mu symbol working with new lines unfortunately
df.plot.tax_order$bench_label <-  factor(df.plot.tax_order$bench_label,
                                         levels = c("ESmu binary","ESmu\ncontinuous","ESmu\ncontinuous^2","ESmu quantile\nnormalised v1",
                                                    "ESmu quantile\nnormalised v2","ESmu binary\npermuted", 'ESmu\ncontinuous\npermuted',
                                                    'ESmu\ncontinuous\n+ noise'),
                                                    ordered = TRUE)
#                                          labels=c(expression(paste("ES",mu, " binary")),
#                                                   expression(paste("ES",mu, " continuous")),
#                                                   expression(paste("ES",mu, " continuous^2")),
#                                                   expression(paste("ES",mu, " quantile normalised v1")),
#                                                   expression(paste("ES",mu, " quantile normalised v2")),
#                                                   expression(paste("ES",mu, " binary permuted")),
#                                                   expression(paste("ES",mu, " continuous permuted")),
#                                                   expression(paste("ES",mu, " continuous + noise"))))
                                                  
'%ni%' <- Negate('%in%')

# Don't plot nulls
# df.plot.tax_order <- df.plot.tax_order %>%
#   arrange(desc(Coefficient_P_value)) %>% 
#   filter(bench_label %ni% c('ESmu binary\npermuted',
#                                    'ESmu\ncontinuous\npermuted',
#                                    'ESmu\ncontinuous\n+ noise'))
# Only plots nulls
df.plot.tax_order <- df.plot.tax_order %>%
  arrange(desc(Coefficient_P_value)) %>% 
  filter(bench_label %in% c('ESmu\ncontinuous',
                            'ESmu binary\npermuted',
                            'ESmu\ncontinuous\npermuted',
                            'ESmu\ncontinuous\n+ noise'))

# Half the number of GWAS so plot isn't too cluttered - split into "main group" and "secondary group"
# Removed RA and MS because there are too many immune cells in different tissues that appear significant
gwas.bmi.paper.group <- c("HEIGHT_UKBB_Loh2018", "EA3_Lee2018", "SCZ_Pardinas2018",
                          "BMI_UKBB_Loh2018", "INSOMNIA_Jansen2018",
                          "INTELLIGENCE_Savage2018", "LIPIDS_LDL_Teslovich2010",
                          "WHRadjBMI_UKBB_Loh2018")
gwas.main = c(gwas.bmi.paper.group, "DEPRESSION_Nagel2018", "ASD_iPSYCH_PGC_Grove2018",
              "BIP_PGC2018", "CARDIOVASCULAR_UKBB_Loh2018", "CROHNS_Jostins2012",
              "IBD_Jostins2012", "MDD_Howard2019", "NEUROTICISM_Nagel2018", "RB_Linner_2019",
              "WORRY_Nagel2018" )
plot.only.maingwas <- df.plot.tax_order %>% 
  filter(`Trait identifier` %in% gwas.main)

plot.other.gwas <- df.plot.tax_order %>% 
  filter(`Trait identifier` %nin% gwas.main)


multi_ES_and_GWAS_scatter_plot <- function(df.ldsc_cts.multi_gwas_and_ES){
  p.multi <- ggplot(df.ldsc_cts.multi_gwas_and_ES, aes(x=-log10(Coefficient_P_value), y=tissue_celltype, color=Tissue)) +
    facet_grid(bench_label ~ `Trait full name`, scales = 'free', space = 'free_y') +
    geom_point(data = filter(df.ldsc_cts.multi_gwas_and_ES, significant_show > 0), size = 1, alpha = 0.4, show.legend = FALSE) +
    geom_text_repel(data = filter(df.ldsc_cts.multi_gwas_and_ES, significant_show > 3), aes(label = annotation_label_fmt),
                    size = 2.2,max.iter = 5000,segment.size = 0.2,force = 10, show.legend = FALSE) +
    geom_vline(xintercept = -log10(0.05/115), color="black", size = 0.1, linetype = "dashed") +
    geom_point(data = filter(df.ldsc_cts.multi_gwas_and_ES, significant_show == 0), size = 0.1, alpha = 0.4) +
    theme_minimal() +
    labs(x =expression(-log[10](P[CELLECT-LDSC]))) +
    theme(legend.position="bottom") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.ticks.x=element_line(color="black"),
          panel.border = element_rect(color = "black",fill= NA, size = 0.1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title=element_text(size=14,face="bold"),
          strip.text = element_text(face="bold",size=8),
          legend.text = element_text(size=12))+
    guides(colour = guide_legend(override.aes = list(size=3))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) 
    
  return(p.multi)
}

p <- multi_ES_and_GWAS_scatter_plot(plot.only.maingwas)
p
p.other  <- multi_ES_and_GWAS_scatter_plot(plot.other.gwas)
p.other
### Save
file.out <- sprintf("figs/fig-compare-null-scatter-many-gwas.tabula-muris.190830.png")
file.out.other.gwas <- sprintf("figs/fig-compare-null-scatter-many-gwas.tabula-muris-other-gwas.190830.png")
ggsave(p, filename=file.out, width=16, height=7)
ggsave(p.other, filename=file.out.other.gwas, width=16, height=7)
