#this code will generate heatmap for epistatic antigenicity of each virus
#loading the library

library(tidyverse)
library(RColorBrewer)


# Load the data
#the data should have columns: position, (ref: not to show on plot), alt, altered aa chemical group, antigenicity(either max or other desired approach)

epiant_df <- read_csv('/Users/rubayetalam/Downloads/bulk_data/cchfv_epistatic_antigenicity_v1.csv')
print(colnames(epiant_df))

#find the df columns, their types and names from the values
numeric_columns <- colnames(epiant_df[which(sapply(epiant_df, class) == 'numeric')])
antigenicity_column <- names(which(head(sapply(epiant_df[which(sapply(epiant_df, class) == 'numeric')], as.integer) == 0)[1,]))
position_column <- numeric_columns[which(numeric_columns != antigenicity_column)]
altered_aa_column <- names(which((head(sapply(epiant_df[which(sapply(epiant_df, class) == 'character')], nchar)) == 1)[1,]))
chemical_group_column <- setdiff(colnames(epiant_df), c(antigenicity_column, position_column, altered_aa_column))


#rename df columns to position, substituted_aa, chemical_group, antigenicity for ease of use in plotting
epiant_df <- epiant_df %>% rename(
  position = !!position_column,
  substituted_aa = !!altered_aa_column,
  chemical_group = !!chemical_group_column,
  antigenicity = !!antigenicity_column)


#make the plot
heatmap_epiant <- epiant_df %>% ggplot() + 
  geom_tile(aes(x=position, 
                y=substituted_aa, 
                fill= antigenicity),
            #color="white",
            #alpha = 0.2, 
            size=0.25) +
  scale_fill_distiller(palette = 'RdBu') +
  #xlim(1300, 1400) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  facet_grid(chemical_group ~ ., scales = "free_y", space = "free_y") +
  labs(
    #title = "Epistatic Antigenicity for different amino acids at each position in CCHFV GPC",
    x = "",
    y = ""
  ) +
  #coord_fixed() +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    #plot.title = element_text(size = 14, face = "bold"),
    #legend.title = element_text(size = 18),
    strip.text.y = element_text(size = 5, face = "bold"), 
    panel.spacing = unit(0.05, "mm"),
    #set thickness of axis ticks
    axis.ticks=element_line(linewidth = 0.4),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank(), 
    legend.position =  'bottom'
  )



heatmap_epiant
#save the plot
ggsave(heatmap_epiant, filename="cchfv-mod1.png", height=8.5, width=18.8, units="in", dpi=300)





make_epiant_heatmap <- function(file_path, virus_name, output_file = NULL) {
  
  # Read data
  epiant_df <- read_csv(file_path, show_col_types = FALSE)
  
  # ---- Detect columns automatically ----
  
  # Numeric columns
  numeric_columns <- colnames(epiant_df)[sapply(epiant_df, is.numeric)]
  
  # Antigenicity column (contains decimals)
  antigenicity_column <- numeric_columns[
    sapply(epiant_df[numeric_columns], function(x) any(x %% 1 != 0))
  ][1]
  
  # Position column (integer-like numeric)
  position_column <- setdiff(numeric_columns, antigenicity_column)[1]
  
  # Amino acid column (single-letter character)
  altered_aa_column <- names(epiant_df)[
    sapply(epiant_df, function(col) {
      is.character(col) && all(nchar(na.omit(col)) == 1)
    })
  ][1]
  
  # Chemical group column (remaining column)
  chemical_group_column <- setdiff(
    colnames(epiant_df),
    c(antigenicity_column, position_column, altered_aa_column)
  )[1]
  
  # ---- Rename columns ----
  epiant_df <- epiant_df %>%
    rename(
      position = !!position_column,
      substituted_aa = !!altered_aa_column,
      chemical_group = !!chemical_group_column,
      antigenicity = !!antigenicity_column
    )
  
  # ---- Plot ----
  heatmap_epiant <- epiant_df %>%
    ggplot() +
    geom_tile(aes(x = position,
                  y = substituted_aa,
                  fill = antigenicity),
              size = 0.25) +
    scale_fill_distiller(palette = "RdBu", na.value = "grey90") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(
                         floor(min(epiant_df$position, na.rm = TRUE) / 50) * 50,
                         ceiling(max(epiant_df$position, na.rm = TRUE) / 50) * 50,
                         by = 50
                       )) +
    facet_grid(chemical_group ~ ., scales = "free_y", space = "free_y") +
    labs(
      #title = paste("Epistatic Antigenicity -", virus_name),
      x = "",
      y = ""
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 8),
      strip.text.y = element_text(size = 6, face = "bold"),
      panel.spacing = unit(0.2, "mm"),
      axis.ticks = element_line(linewidth = 0.4),
      legend.position = "none"
    )
  
  # ---- Save file ----
  if (is.null(output_file)) {
    output_file <- paste0(tolower(virus_name), "-epiant-heatmap.png")
  }
  
  ggsave(
    filename = output_file,
    plot = heatmap_epiant,
    height = 3,
    width = 8.27,
    units = "in",
    dpi = 300
  )
  
  return(heatmap_epiant)
}

make_epiant_heatmap('/Users/rubayetalam/epiant_marburg.csv', 
                    'Marburg', 
                    'marburg-epiant-heatmap.png')

make_epiant_heatmap('/Users/rubayetalam/epiant_cchfv.csv', 
                    'CCHFV', 
                    'cchfv-epiant-heatmap.png')

make_epiant_heatmap('/Users/rubayetalam/epiant_borna.csv', 
                    'Bornavirus-1', 
                    'borna-epiant-heatmap.png')


#get legend from one of the plots
heatmap <- epiant_df %>% ggplot() + 
  geom_tile(aes(x=position, 
                y=substituted_aa, 
                fill= antigenicity),
            #color="white",
            #alpha = 0.2, 
            size=0.25) +
  scale_fill_distiller(palette = 'RdBu') +
  #xlim(1300, 1400) +
  #theme_void() +
  #guides(color = "none", linetype = "none", shape = "none") 
  theme(legend.position = "bottom", 
        legend.title = element_text(size = 24, face = 'bold'), 
        legend.text = element_text(size = 20), 
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(2, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        #legend.title = element_text(size=14), #change legend title font size
        #legend.text = element_text(size=10))
  )
heatmap

legend <- get_legend(heatmap)
as_ggplot(legend)
