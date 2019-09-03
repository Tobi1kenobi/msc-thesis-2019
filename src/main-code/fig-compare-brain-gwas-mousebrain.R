############### SYNOPSIS ###################
# Make Mouse Nervous System heatmap for 17 brain traits

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

`%nin%` = Negate(`%in%`)

# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

### Read LDSC results
dir.data <- here("out/CELLECT-LDSC-out/mousebrain")
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

# 17 brain traits that are representative and don't violate assumptions + height hackily added
brain.data.w.height <- data %>%
  left_join(GWAS.info, by="Trait identifier") %>% 
  filter(FLAG_final & GWAS.trait.representative & is.na(`GWAS violates LDSC assumptions`) & `Trait - brain`| `Trait identifier` == "HEIGHT_UKBB_Loh2018")



# ======================================================================= #
# ========================= LOAD CELL-TYPE METADATA ===================== #
# ======================================================================= #

### Read annotation metadata
df.metadata <- get_metadata(dataset_prefix="mousebrain_all")
df.metadata

### Add meta data
df.ldsc_cts <- brain.data.w.height %>% left_join(df.metadata, by="annotation") # add meta data
df.ldsc_cts

### Order annotations
df.ldsc_cts <- df.ldsc_cts %>% mutate(annotation=factor(annotation, levels=unique(annotation[order(tax_order_idx_mb_fig1c, as.character(annotation))]))) # order annotations by their taxonomy (from MB Fig1c), and se
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

### Rename
df.plot <- df.ldsc_cts %>% mutate(TaxonomyRank4_reduced1 = case_when(
  TaxonomyRank4_reduced1 == "Di- and mesencephalon excitatory neurons" ~ "Di- and mesencephalon\nexcitatory neurons",
  TaxonomyRank4_reduced1 == "Cholinergic and monoaminergic neurons" ~ "Cholinergic and\nmonoaminergic neurons",
  TaxonomyRank4_reduced1 == "Di- and mesencephalon inhibitory neurons" ~ "Di- and mesencephalon\ninhibitory neurons",
  TaxonomyRank4_reduced1 == "Telencephalon projecting excitatory neurons" ~ "Telencephalon projecting\nexcitatory neurons",
  TaxonomyRank4_reduced1 == "Telencephalon projecting inhibitory neurons" ~ "Telencephalon projecting\ninhibitory neurons",
  TaxonomyRank4_reduced1 == "Telencephalon inhibitory interneurons" ~ "Telencephalon\ninhibitory interneurons",
  TRUE ~ as.character(TaxonomyRank4_reduced1))
)


no.sig.taxons <- df.plot %>% 
  group_by(TaxonomyRank4_reduced1) %>% 
  filter(min(Coefficient_P_value) >= 0.05/265) %>% 
  ungroup %>% 
  pull(TaxonomyRank4_reduced1) %>% 
  unique()

bonf_thresh = 0.05/length(unique(df.plot$annotation))
-log10(bonf_thresh)

multi_ES_metric_heatmap_plot <- function(df.ldsc_cts.multi_gwas, x_lab_size) {
  
  ### Plot
  p.multi_gwas <- ggplot(df.ldsc_cts.multi_gwas, aes(y=`Trait full name`, x=annotation, fill=-log10(Coefficient_P_value))) +
    geom_tile() + 
    geom_text(aes(label=fdr_significant), color="white", hjust=0.5, vjust=0.75, size = 2.2, angle = 90) + # add asterisk if fdr significant
    
    colorspace::scale_fill_continuous_sequential(palette="Blues 2", rev=TRUE, limits=c(0,3.724276), oob=scales::squish, na.value = "white") + # "Blues 2","Blues 3","Purples 3"
    # ^ oob: Function that handles limits outside of the scale limits (out of bounds). scales::squish "squishes" values into the range specified by limits
    # labs(fill=expression(-log[10](P[S-LDSC]))) +
    theme_minimal() + # set this first
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) + # remove grid
    theme(axis.text.x=element_text(angle=90, hjust=1, size=x_lab_size, face = "bold")) +
    theme(legend.position="bottom") +
    facet_grid(`Trait category` ~ TaxonomyRank4_reduced1, scales = 'free', space = 'free') +
    theme(panel.spacing = unit(0.1, "lines")) + # Decrease space between plots 
    theme(strip.text.x=element_text(size = 10,angle = 90, face = "bold"), # Controls facet label appearance
          strip.text.y=element_text(size = 9, angle = 0, face = "bold"),
          axis.title=element_text(size=14,face="bold"),
          legend.text = element_text(size=12))+
    ylab('Complex trait') +
    xlab('Cell type') +
    labs(fill = expression(-log[10](P[CELLECT-LDSC])))
  return(p.multi_gwas)
}


p <- multi_ES_metric_heatmap_plot(filter(df.plot, TaxonomyRank4_reduced1 %nin% no.sig.taxons), x_lab_size = 4.5)
p.other <-  multi_ES_metric_heatmap_plot(filter(df.plot, TaxonomyRank4_reduced1 %in% no.sig.taxons), x_lab_size = 6)
p.other
p
### Save
file.out <- sprintf("figs/fig-compare-brain-gwas.mousebrain.190830.png")
file.out.other.tissues <- sprintf("figs/fig-compare-brain-gwas.mousebrain.190830-other-tissues.png")
ggsave(p, filename=file.out, width=18, height=6)
ggsave(p.other, filename=file.out.other.tissues, width=10, height=6)

