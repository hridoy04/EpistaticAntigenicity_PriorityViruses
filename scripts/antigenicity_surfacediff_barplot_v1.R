#this code will generate the barplot comparing the surface-exposed and buried sites in a protein

# Load necessary libraries
library(tidyverse)

epi_surface_analysis <- function(
    sasa_file,
    epi_file,
    sasa_pos_col,
    sasa_value_col,
    epi_pos_col,
    epi_value_col,
    probability_col = NULL,
    rsa_threshold = 25,
    position_offset = 0,
    virus_name = "Virus Protein"
){
  
  # Load SASA
  sasa <- read_csv(sasa_file)
  
  rsa_df <- sasa %>%
    mutate(
      rsa = (.data[[sasa_value_col]] * 100 / max(.data[[sasa_value_col]], na.rm = TRUE)),
      surf = case_when(
        rsa >= rsa_threshold ~ "exposed",
        TRUE ~ "buried"
      ),
      position = .data[[sasa_pos_col]] + position_offset
    ) %>%
    select(position, sasa = all_of(sasa_value_col), rsa, surf)
  
  # Load epistatic antigenicity
  epi <- read_csv(epi_file)
  
  epi_df <- epi %>%
    rename(
      position = !!epi_pos_col,
      epi = !!epi_value_col
    )
  
  if(!is.null(probability_col)){
    epi_df <- epi_df %>%
      rename(probability = !!probability_col)
  }
  
  # Join
  merged <- epi_df %>%
    right_join(rsa_df, by = "position")
  
  merged$epi <- merged$epi %>% replace_na(0)
  
  # Plot
  p <- merged %>%
    ggplot() +
    geom_boxplot(aes(x = surf, y = epi, color = surf)) +
    labs(
      title = paste0(virus_name),
      subtitle = "Epistatic Antigenicity grouped by Surface Accessibility",
      x = "Surface Accessibility Groups",
      y = "Epistatic Antigenicity"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold")
    )
  
  # Wilcoxon test
  w <- wilcox.test(data = merged, epi ~ surf)
  
  return(list(
    data = merged,
    plot = p,
    wilcox = w
  ))
}


#example 

borna_result <- epi_surface_analysis(
  sasa_file = "bornavirus_p57_sasa.csv",
  epi_file = "epiant/borna.csv",
  sasa_pos_col = "pos",
  sasa_value_col = "sasa",
  epi_pos_col = 'position',
  epi_value_col = 'significant_antigenicity_median',
  probability_col = 'probability',
  rsa_threshold = 25,
  virus_name = "Bornavirus p57"
)

borna_result$plot
borna_result$wilcox
