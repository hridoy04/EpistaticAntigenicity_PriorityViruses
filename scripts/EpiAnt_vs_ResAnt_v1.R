

library(dplyr)
library(ggplot2)
library(ggrepel)



#1. Function to Compare Ranks
# This computes:
#   rank of epistatic antigenicity
# rank of residual antigenicity
# Δ rank
# Wilcoxon paired test

compare_antigenicity_ranks <- function(data,
                                       position_col,
                                       epi_col,
                                       res_col,
                                       positions_of_interest = NULL) {
  
  df <- data %>%
    mutate(
      rank_epi = rank({{epi_col}}, ties.method = "min"),
      rank_res = rank({{res_col}}, ties.method = "min")
    )
  
  # Filter only escape positions if provided
  if(!is.null(positions_of_interest)){
    df <- df %>% filter({{position_col}} %in% positions_of_interest)
  }
  
  df <- df %>%
    mutate(delta_rank = rank_epi - rank_res)
  
  # Statistical test
  wilcox_result <- wilcox.test(
    df$rank_epi,
    df$rank_res,
    paired = TRUE,
    exact = FALSE,
    alternative = "greater"
  )
  
  list(
    ranked_data = df,
    wilcox_test = wilcox_result
  )
}



# 2. Function for Scatter Plot
# Shows relationship between:
#   Epistatic antigenicity vs experimental escape
# Points colored by escape score and labels for escape positions.


plot_escape_scatter <- function(data,
                                x_var,
                                y_var,
                                position_var,
                                highlight_positions = NULL,
                                x_label = "Epistatic Antigenicity",
                                y_label = "Experimental Escape Score") {
  
  ggplot(data, aes({{x_var}}, {{y_var}})) +
    
    geom_point(
      aes(color = {{y_var}}, size = {{y_var}}),
      alpha = 0.9
    ) +
    
    scale_color_distiller(
      palette = "Spectral",
      name = "Experimental\nDMS score"
    ) +
    
    scale_size(range = c(2,6)) +
    
    ggrepel::geom_text_repel(
      aes(label = ifelse({{position_var}} %in% highlight_positions,
                         {{position_var}}, "")),
      size = 5,
      max.overlaps = Inf
    ) +
    
    labs(
      x = x_label,
      y = y_label
    ) +
    
    theme_bw() +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 20),
      legend.title = element_text(size = 16)
    )
}

#examples:

sc2_rank_results <- compare_antigenicity_ranks(
  data = sc2_man,
  position_col = position,
  epi_col = rank_epi,
  res_col = rank_acc,
  positions_of_interest = sc2_abescape_pos
)

plot_escape_scatter(
  data = sc2_epi_all,
  x_var = accessibility,
  y_var = Ab_escape,
  position_var = position,
  highlight_positions = sc2_abescape_pos,
  x_label = "Residual Antigenicity",
  y_label = "Experimental DMS score"
)



sc2_rank_results$wilcox_test