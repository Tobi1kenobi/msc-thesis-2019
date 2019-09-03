############### SYNOPSIS ###################
# Make Tabula Muris relevant cell types


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ================== IMPORTANT LEGEND INFORMATION ======================= #
# ======================================================================= #

# SEE MOUSEBRAIN

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
# ============================ PROCESS DATA =============================== #
# ======================================================================= #

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


### Create 'plotting' data frame
# my.gwas <-  "BMI_UKBB_Loh2018"

scatter_heatmap <- function(gwas.in){
  
  df.plot.tax_order <- df.ldsc_cts %>%filter(gwas == gwas.in, benchmark == "tabula_muris.ESmu.continuous.benchmark190722")  # filter for
  
  ### annotation: create clean names
  df.plot.tax_order <- df.plot.tax_order %>% mutate(annotation_label_fmt = tm_annotation_formatter(annotation))
  df.plot.tax_order
  
  ### tissue: rename + remove some text because they contain too few cell-types
  df.plot.tax_order %>% count(tissue, sort=T) # ---> potentially filter : n_annotations_in_tax >= 5
  df.plot.tax_order <- df.plot.tax_order %>% mutate(tissue = case_when(
    tissue == "Brain_Myeloid" ~ "Brain",
    tissue == "Brain_Non-Myeloid" ~ "Brain",
    TRUE ~ stringr::str_replace_all(tissue, pattern="_", replacement=" ") # TRUE ~ as.character(tissue)
    )
  )
  
  # ======================================================================= #
  # ============================ BMI point plot =============================== #
  # ======================================================================= #
  
  
  ### Create data frame with taxonomy text data
  df.tax_text_position <- df.plot.tax_order %>% 
    group_by(tissue) %>% 
    summarize(n_annotations_in_tax=n()) %>% # count how many annotations in each tax
    left_join(df.metadata %>% select(tissue) %>% distinct(), by="tissue") %>%
    # ^ OBS: we use 'distinct()' because df.taxonomy_metadata contains (intentional) duplicated rows for some combinations of {tissue, order_idx_mb_taxonomy_fig1c}.
    mutate(tissue = factor(tissue, levels=tissue[order(tissue)])) # IMPORTANT: order factor by the SAME ORDER as 'annotation' column is ordered by in df.plot.tax_order
  ### Add information of the text's position in the plot
  df.tax_text_position <- df.tax_text_position %>% 
    arrange(tissue) %>%
    mutate(pos_mid=cumsum(n_annotations_in_tax)-n_annotations_in_tax/2,
           pos_start=cumsum(n_annotations_in_tax)-n_annotations_in_tax+0.5,
           pos_end=cumsum(n_annotations_in_tax)+0.5,
           idx=1:n()) # idx is used for identifying every n'th tax
  df.tax_text_position <- df.tax_text_position %>% mutate(flag_draw_rect=if_else(idx %% 2 == 0, TRUE, FALSE)) 
  
  
  ### Set 'tissue_display_name' name
  # Some tissues have too few annotations for us to display their name
  df.plot.tax_order %>% count(tissue, sort=T) 
  df.tax_text_position <- df.tax_text_position %>% mutate(tissue_display_name = case_when(
    tissue == "Tongue" ~ "", 
    TRUE ~ as.character(tissue)
    )
  )
  
  ### Reverse tax order [only relevant for gganimate direction. x axis is plottet by annotation].
  df.plot.tax_order <- df.plot.tax_order %>% mutate(tissue = fct_rev(tissue))
  df.tax_text_position <- df.tax_text_position %>% mutate(tissue = fct_rev(tissue), tissue.text.location = log10(min(df.plot.tax_order$Coefficient_P_value))/4.2) # Second mutation makes a custom location depending on the size of plot
  
  fdr_threshold <- c(0.05, 0.01, 0.001)/nrow(df.plot.tax_order)
  p.main <- ggplot() +
    ### add ggplot 'baselayer'. This makes our 'canvas' and needs to go first (for some reason I don't fully understand...)
    geom_point(data=df.plot.tax_order, aes(x=annotation, y=-log10(Coefficient_P_value)), color="gray") +
    ### add tax info and gray rects (add gray rect before the data points, so avoid them 'covering' the points)
    geom_text(data=df.tax_text_position, aes(x=pos_mid, y=tissue.text.location, label=tissue_display_name, group=tissue), hjust="right", size=rel(3)) + # group=tissue solves gganimate problem with 'jumps' in the position. It is not needed for a static ggplot
    geom_rect(data=df.tax_text_position %>% filter(flag_draw_rect), aes(xmin=pos_start, xmax=pos_end, ymin=-15, ymax=Inf), color=NA, alpha=0.1) +
    ### cell-types
    geom_point(data=df.plot.tax_order %>% filter(fdr_significant != ''), aes(x=annotation, y=-log10(Coefficient_P_value), color=tissue)) +
    ggrepel::geom_text_repel(data=df.plot.tax_order %>% filter(fdr_significant != ''), aes(x=annotation, y=-log10(Coefficient_P_value), label=annotation_label_fmt, color=tissue, group=tissue), hjust = 0, nudge_x = 1.5, show.legend=F) + # group is to try to improve gganimate
    ### extra
    geom_hline(yintercept=-log10(fdr_threshold[1]), linetype="dashed", color="darkgray") + 
    ### axes
    labs(x="", y=expression(-log[10](P[S-LDSC]))) +
    # coord
    coord_flip(ylim = c( 0, max(-log10(df.plot.tax_order$Coefficient_P_value)) ), # This focuses the y-axis on the range of interest
               clip = 'off') +   # This keeps the labels from disappearing
    # ^ clip = 'off': it allows drawing of data points anywhere on the plot, including in the plot margins. If limits are set via xlim and ylim and some data points fall outside those limits, then those data points may show up in places such as the axes, the legend, the plot title, or the plot margins.
    # ^ clip = 'off': disable cliping so df.tax_text_position text annotation can be outside of plot | REF https://stackoverflow.com/a/51312611/6639640
    ### guides
    #...
    ### scale
    scale_x_discrete(label=tm_annotation_formatter) +
    ### Color
    #scale_fill_brewer(palette="Paired") +
    # scale_fill_brewer(palette="Dark2") +
    # scale_fill_brewer(palette="Set2") +
    ### theme
    theme_classic() + 
    theme(axis.text.y=element_text(size=rel(0.5))) + # REF: https://stackoverflow.com/a/47144823/6639640
    theme(legend.position="bottom")
  
  p.main
  
  ### Remove legend
  p.main <- p.main + guides(fill=F, color=F)
  ### Add margin to plot if displaying it directly (without pathwork)
  p.main.margin <- p.main + theme(plot.margin = unit(c(1,1,1,5), "cm")) # (t, r, b, l) widen left margin
  p.main.margin
  
  # file.out <- sprintf("figs/fig_celltypepriori.tm.pdf")
  #ggsave(p.main.margin, filename=file.out, width=9, height=8)
  
  
  # ======================================================================= #
  # =================== HEATMAP + BARPLOT =================== #
  # ======================================================================= #
  
  
  df.plot.heatmap <- df.ldsc_cts %>%filter(gwas == gwas.in, !grepl("[.][1-4]$",benchmark))
  
  set_multi_ES_metric_heatmap_plot <- function(df.ldsc_cts.one.gwas) {
    # ### Add "spacing block" column
    # df.ldsc_cts.one.gwas <- df.ldsc_cts.one.gwas %>% mutate(space_block=get_gwas_space_block(gwas))
    
    # ### Order ES metrics
    # df.plot.multi_gwas <- df.plot.multi_gwas %>% mutate(gwas=factor(gwas, levels=filter.gwas)) # Order GWAS
    # df.plot.multi_gwas
    
    ### Rename GWAS [do this after filtering, so ensure uniqueness]
    ### IMPORANTLY recode_factor() allows us to rename the factor values but KEEP the ORDER defined by filter.gwas
    
    df.ldsc_cts.one.gwas <- df.ldsc_cts.one.gwas %>%
      extract(benchmark, "ES.metric", "[.](.*)[.]b", remove = TRUE)
    
    
    ### Plot
    p.multi_metric <- ggplot(df.ldsc_cts.one.gwas, aes(x=ES.metric, y=annotation, fill=-log10(Coefficient_P_value))) + 
      geom_tile() + 
      geom_text(aes(label=fdr_significant), color="red", hjust=0.5, vjust=0.75) + # add asterisk if fdr significant
      # ^ hjust/vjust: https://stackoverflow.com/questions/7263849/what-do-hjust-and-vjust-do-when-making-a-plot-using-ggplot
      # ^ hjust="center", vjust="middle"
      # scale_fill_viridis_c(option="magma", direction=-1) + 
      # scale_fill_distiller(palette="Greys", direction=1) + # "Greys" "Blues" "Greens" "BluGn" "Reds"
      # [GOOD] scale_fill_distiller(palette="Blues", direction=1, limits=c(0,5), oob=scales::squish) + # "Greys" "Blues" "Greens" "BluGn" "Reds"
      # scale_fill_distiller(palette="Greys", direction=1, limits=c(0,1.8), na.value = "white") + # tau norm plot
      colorspace::scale_fill_continuous_sequential(palette="Blues 2", rev=TRUE, limits=c(0,3.361728), oob=scales::squish, na.value = "white") + # "Blues 2","Blues 3","Purples 3"
      # ^ oob: Function that handles limits outside of the scale limits (out of bounds). scales::squish "squishes" values into the range specified by limits
      labs(fill=expression(-log[10](P[S-LDSC]))) +
      theme_minimal() + # set this first
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank()) + # remove grid
      theme(axis.text.x=element_text(angle=90, hjust=1)) + 
      theme(legend.position="bottom") +
      theme(axis.text.y=element_blank(), 
            axis.ticks.y = element_blank(),
            axis.title = element_blank())
    
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
    # p.multi_metric
    # 
    ### Make global
    # "<<-" to set global variables | REF: https://stackoverflow.com/a/1236721
    p.multi_metric <<- p.multi_metric
    print("setting global vars: p.multi_metric")
  }
  
  
  set_multi_ES_metric_heatmap_plot(df.plot.heatmap)
  # set_h2_barplot()
  
  
  
  # ======================================================================= #
  # ============================ TOP3 cell-types ========================= #
  # ======================================================================= #
  
  # df.top <- df.ldsc_cts %>% 
  #   filter(gwas %in% filter.gwas) %>% 
  #   group_by(gwas) %>% 
  #   top_n(n=3, wt=-log10(p.value))
  # # df.top
  
  # ======================================================================= #
  # =================== COMBINE PLOTS: MAIN + HEATMAP =================== #
  # ======================================================================= #
  
  # ### Make blank plot to funciton as a 'dummy margin plot'
  p.blank <- ggplot(df.plot.tax_order, aes(x=annotation, y=-log10(Coefficient_P_value))) +
    geom_blank() +
    theme(panel.grid = element_blank(),axis.line = element_blank(),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.background = element_blank())

  # ### Patch - with barplot
  # p.patch <- 
  #   (
  #     # LEFT
  #     (p.h2.blank / 
  #        (p.blank + p.main + plot_layout(widths = c(6,10)))
  #      + plot_layout(heights = c(1,10))
  #     ) |
  #       # ) & theme(plot.margin =  unit(c(1,1,1,1), "pt")) | # dont know why "& theme()" does not recurse into nested levels
  #       # RIGHT
  #       (p.h2_no_x /
  #          p.multi_metric 
  #        + plot_layout(heights = c(1,10))
  #       )
  #     # RIGHT + LEFT
  #   ) + plot_layout(widths = c(10,3)) & theme(plot.margin =  unit(c(1,1,1,1), "pt"))
  # p.patch
  
  ### Patch (no h2 barplot)
  p.patch <- p.multi_metric + p.blank +  p.main  + plot_layout(ncol = 3, widths = c(3,1,10))
  p.patch
  
  ### Save
  file.out <- sprintf("figs/fig_relevance.tm.heatmap_scatter.190727.%s.pdf", gwas.in)
  ggsave(p.patch, filename=file.out, width=12, height=10)
  return(file.out)
}

all_gwas <- unique(df.ldsc_cts$gwas)

lapply(all_gwas, scatter_heatmap)
