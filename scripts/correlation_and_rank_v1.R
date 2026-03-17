library(tidyverse)
library(patchwork)
library(corrplot)
library(RColorBrewer)

# Function to create correlation plot
create_correlation_plot <- function(df, virus_protein_name,
                                    columns, filter_df = NULL, label_col = "label",
                                    order = "alphabet", method = "color", 
                                    diag = FALSE, type = "full") {
  
  # Select specified columns
  df_selected <- df %>% 
    select(all_of(c(label_col, columns)))
  
  # Apply filtering if filter_df is provided
  if (!is.null(filter_df)) {
    filter_labels <- filter_df[[label_col]]
    df_selected <- df_selected %>% 
      filter(.data[[label_col]] %in% filter_labels)
  }
  
  # Create correlation matrix
  cor_matrix <- df_selected %>%
    select(-all_of(label_col)) %>%
    drop_na() %>%
    as.matrix() %>%
    cor()
  
  # Create correlation plot
  corrplot(cor_matrix, 
           main = virus_protein_name,
           order = order, 
           tl.col = 'black', 
           diag = diag,
           type = type,
           method = method,
           addCoef.col = 'black',
           cl.cex = 1,
           tl.cex = 2,
           pch.cex = 1.5,
           col = rev(brewer.pal(n = 10, name = "RdYlBu")))
  
  return(cor_matrix)
}

# Main function to compare ranks
compare_ranks <- function(df, columns, reference_col, filter_df = NULL, label_col = "label") {
  
  # Select and rank the specified columns
  df_ranked <- df %>% 
    select(all_of(c(label_col, columns))) %>%
    mutate(across(-all_of(label_col), ~ rank(.x, ties.method = "min"), .names = "{.col}_rank")) %>%
    select(all_of(label_col), ends_with('_rank'))
  
  # Apply filtering if filter_df is provided
  if (!is.null(filter_df)) {
    filter_labels <- filter_df[[label_col]]
    df_ranked <- df_ranked %>% 
      filter(.data[[label_col]] %in% filter_labels)
  }
  
  # Remove label column after filtering
  df_ranked <- df_ranked %>% select(-all_of(label_col))
  
  # Reorder so reference column is first
  ref_col_name <- paste0(reference_col, "_rank")
  df_ranked <- df_ranked %>% select(all_of(ref_col_name), everything())
  
  # Generate point plots
  plot_list <- list()
  other_cols <- setdiff(colnames(df_ranked), ref_col_name)
  
  for (col in other_cols) {
    p <- ggplot(df_ranked, aes(x = .data[[ref_col_name]], y = .data[[col]])) + 
      geom_point(alpha = 0.6, size = 2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(x = gsub("_rank", "", ref_col_name),
           y = gsub("_rank", "", col),
           title = paste("Rank Comparison:", gsub("_rank", "", col))) +
      theme_minimal() +
      theme(plot.title = element_text(size = 12, face = "bold"))
    
    plot_list[[col]] <- p
  }
  
  # Generate bar plot comparing reference to all others
  comparison_data <- data.frame()
  
  for (col in other_cols) {
    greater_count <- sum(df_ranked[[ref_col_name]] >= df_ranked[[col]], na.rm = TRUE)
    less_count <- sum(df_ranked[[ref_col_name]] < df_ranked[[col]], na.rm = TRUE)
    
    comparison_data <- rbind(
      comparison_data,
      data.frame(
        Method = gsub("_rank", "", col),
        Greater = greater_count,
        Less = less_count
      )
    )
  }
  
  # Create bar plot
  comparison_long <- comparison_data %>%
    pivot_longer(cols = c(Greater, Less), 
                 names_to = 'Comparison', 
                 values_to = 'Count')
  
  bar_plot <- ggplot(comparison_long, aes(x = Method, y = Count, fill = Comparison)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    labs(title = paste('Rank Comparison:', gsub("_", " ", reference_col), 
                       'vs Other Metrics'),
         x = 'Comparison Metrics',
         y = 'Count') +
    scale_fill_manual(values = c("Greater" = "#E69F00", "Less" = "#56B4E9"),
                      labels = c("Greater" = paste(reference_col, "> Other"),
                                 "Less" = paste(reference_col, "< Other"))) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 11))
   
  # Return all plots
  return(list(
    point_plots = plot_list,
    bar_plot = bar_plot,
    comparison_table = comparison_data
  ))
}


# Main function to compare ranks
compare_ranks <- function(df, columns, reference_col, filter_df = NULL, label_col = "label") {
  
  # Select and rank the specified columns
  df_ranked <- df %>% 
    select(all_of(c(label_col, columns))) %>%
    mutate(across(-all_of(label_col), ~ rank(.x, ties.method = "min"), .names = "{.col}_rank")) %>%
    select(all_of(label_col), ends_with('_rank'))
  
  # Apply filtering if filter_df is provided
  if (!is.null(filter_df)) {
    filter_labels <- filter_df[[label_col]]
    df_ranked <- df_ranked %>% 
      filter(.data[[label_col]] %in% filter_labels)
  }
  
  # Remove label column after filtering
  df_ranked <- df_ranked %>% select(-all_of(label_col))
  
  # Reorder so reference column is first
  ref_col_name <- paste0(reference_col, "_rank")
  df_ranked <- df_ranked %>% select(all_of(ref_col_name), everything())
  
  # Generate point plots
  plot_list <- list()
  other_cols <- setdiff(colnames(df_ranked), ref_col_name)
  
  for (col in other_cols) {
    p <- ggplot(df_ranked, aes(x = .data[[ref_col_name]], y = .data[[col]])) + 
      geom_point(alpha = 0.6, size = 2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(x = gsub("_rank", "", ref_col_name),
           y = gsub("_rank", "", col),
           title = paste("Rank Comparison:", gsub("_rank", "", col))) +
      theme_minimal() +
      theme(plot.title = element_text(size = 12, face = "bold"))
    
    plot_list[[col]] <- p
  }
  
  # Generate bar plot comparing reference to all others
  comparison_data <- data.frame()
  total_count <- nrow(df_ranked)
  
  for (col in other_cols) {
    greater_count <- sum(df_ranked[[ref_col_name]] >= df_ranked[[col]], na.rm = TRUE)
    less_count <- sum(df_ranked[[ref_col_name]] < df_ranked[[col]], na.rm = TRUE)
    
    comparison_data <- rbind(
      comparison_data,
      data.frame(
        Method = gsub("_rank", "", col),
        Greater = greater_count,
        Less = less_count,
        Greater_pct = (greater_count / total_count) * 100,
        Less_pct = (less_count / total_count) * 100
      )
    )
  }
  
  # Create bar plot with percentages
  comparison_long <- comparison_data %>%
    pivot_longer(cols = c(Greater_pct, Less_pct), 
                 names_to = 'Comparison', 
                 values_to = 'Percentage') %>%
    mutate(Comparison = gsub("_pct", "", Comparison))
  
  bar_plot <- ggplot(comparison_long, aes(x = Method, y = Percentage, fill = Comparison)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 3.5) +
    labs(title = paste('Rank Comparison:', gsub("_", " ", reference_col), 
                       'vs Other Metrics (%)'),
         x = 'Comparison Metrics',
         y = 'Percentage (%)') +
    scale_fill_manual(values = c("Greater" = "#E69F00", "Less" = "#56B4E9"),
                      labels = c("Greater" = paste(reference_col, ">= Other"),
                                 "Less" = paste(reference_col, "< Other"))) +
    ylim(0, 110) +  # Extra space for labels
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 11))
   
  # Return all plots
  return(list(
    point_plots = plot_list,
    bar_plot = bar_plot,
    comparison_table = comparison_data
  ))
}


#new rank plot:


# Main function to compare ranks
compare_ranks_perc <- function(df, columns, reference_col, filter_df = NULL, label_col = "label") {
  
  # Select and rank the specified columns
  df_ranked <- df %>% 
    select(all_of(c(label_col, columns))) %>%
    mutate(across(-all_of(label_col), ~ rank(.x, ties.method = "min"), .names = "{.col}_rank")) %>%
    select(all_of(label_col), ends_with('_rank'))
  
  # Apply filtering if filter_df is provided
  if (!is.null(filter_df)) {
    filter_labels <- filter_df[[label_col]]
    df_ranked <- df_ranked %>% 
      filter(.data[[label_col]] %in% filter_labels)
  }
  
  # Remove label column after filtering
  df_ranked <- df_ranked %>% select(-all_of(label_col))
  
  # Reorder so reference column is first
  ref_col_name <- paste0(reference_col, "_rank")
  df_ranked <- df_ranked %>% select(all_of(ref_col_name), everything())
  
  # Generate point plots
  plot_list <- list()
  other_cols <- setdiff(colnames(df_ranked), ref_col_name)
  
  for (col in other_cols) {
    p <- ggplot(df_ranked, aes(x = .data[[ref_col_name]], y = .data[[col]])) + 
      geom_point(alpha = 0.6, size = 2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(x = gsub("_rank", "", ref_col_name),
           y = gsub("_rank", "", col),
           title = paste("Rank Comparison:", gsub("_rank", "", col))) +
      theme_minimal() +
      theme(plot.title = element_text(size = 12, face = "bold"))
    
    plot_list[[col]] <- p
  }
  
  # Generate bar plot comparing reference to all others
  comparison_data <- data.frame()
  total_count <- nrow(df_ranked)
  
  for (col in other_cols) {
    greater_count <- sum(df_ranked[[ref_col_name]] >= df_ranked[[col]], na.rm = TRUE)
    less_count <- sum(df_ranked[[ref_col_name]] < df_ranked[[col]], na.rm = TRUE)
    
    comparison_data <- rbind(
      comparison_data,
      data.frame(
        Method = gsub("_rank", "", col),
        Greater = greater_count,
        Less = less_count,
        Greater_pct = (greater_count / total_count) * 100,
        Less_pct = (less_count / total_count) * 100
      )
    )
  }
  
  # Create bar plot with percentages
  comparison_long <- comparison_data %>%
    pivot_longer(cols = c(Greater_pct, Less_pct), 
                 names_to = 'Comparison', 
                 values_to = 'Percentage') %>%
    mutate(Comparison = gsub("_pct", "", Comparison))
  
  bar_plot <- ggplot(comparison_long, aes(x = Method, y = Percentage, fill = Comparison)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 3.5) +
    labs(title = paste('Rank Comparison:', gsub("_", " ", reference_col), 
                       'vs Other Metrics (%)'),
         x = 'Comparison Metrics',
         y = 'Percentage (%)') +
    scale_fill_manual(values = c("Greater" = "#E69F00", "Less" = "#56B4E9"),
                      labels = c("Greater" = paste(reference_col, ">= Other"),
                                 "Less" = paste(reference_col, "< Other"))) +
    ylim(0, 110) +  # Extra space for labels
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 11))
  
  # Return all plots
  return(list(
    point_plots = plot_list,
    bar_plot = bar_plot,
    comparison_table = comparison_data
  ))
}

#another bar plot type

# Main function to compare ranks
compare_ranks_onebar <- function(df, columns, reference_col, filter_df = NULL, label_col = "label") {
  
  # Select and rank the specified columns
  df_ranked <- df %>% 
    select(all_of(c(label_col, columns))) %>%
    mutate(across(-all_of(label_col), ~ rank(.x, ties.method = "min"), .names = "{.col}_rank")) %>%
    select(all_of(label_col), ends_with('_rank'))
  
  # Apply filtering if filter_df is provided
  if (!is.null(filter_df)) {
    filter_labels <- filter_df[[label_col]]
    df_ranked <- df_ranked %>% 
      filter(.data[[label_col]] %in% filter_labels)
  }
  
  # Remove label column after filtering
  df_ranked <- df_ranked %>% select(-all_of(label_col))
  
  # Reorder so reference column is first
  ref_col_name <- paste0(reference_col, "_rank")
  df_ranked <- df_ranked %>% select(all_of(ref_col_name), everything())
  
  # Generate point plots
  plot_list <- list()
  other_cols <- setdiff(colnames(df_ranked), ref_col_name)
  
  for (col in other_cols) {
    p <- ggplot(df_ranked, aes(x = .data[[ref_col_name]], y = .data[[col]])) + 
      geom_point(alpha = 0.6, size = 2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(x = gsub("_rank", "", ref_col_name),
           y = gsub("_rank", "", col),
           title = paste("Rank Comparison:", gsub("_rank", "", col))) +
      theme_minimal() +
      theme(plot.title = element_text(size = 12, face = "bold"))
    
    plot_list[[col]] <- p
  }
  
  # Generate bar plot comparing reference to all others
  comparison_data <- data.frame()
  total_count <- nrow(df_ranked)
  
  for (col in other_cols) {
    greater_count <- sum(df_ranked[[ref_col_name]] >= df_ranked[[col]], na.rm = TRUE)
    
    comparison_data <- rbind(
      comparison_data,
      data.frame(
        Method = gsub("_rank", "", col),
        Greater_pct = (greater_count / total_count) * 100
      )
    )
  }
  
  # Create bar plot with percentages
  bar_plot <- ggplot(comparison_data, aes(x = Method, y = Greater_pct)) +
    geom_bar(stat = 'identity', fill = "#E69F00") +
    geom_text(aes(label = sprintf("%.1f%%", Greater_pct)), 
              vjust = -0.5, size = 4) +
    labs(subtitle = paste('Higher Ranks in Epistatic Antigenicity \ncompared to other antigenicity metrics'),
         x = 'Comparison Antigenicity Metrics',
         y = 'Percentage (%)') +
    ylim(0, 110) +  # Extra space for labels
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold"),
          plot.subtitle = element_text(size = 18, face = "bold"))
  
  # Return all plots
  return(list(
    point_plots = plot_list,
    bar_plot = bar_plot,
    comparison_table = comparison_data
  ))
}

# Example usage:


###HIV:

# 1. Create correlation plot first:
cor_matrix <- create_correlation_plot(
  df = hivbg505_df_final,
  virus_protein_name = 'HIV(BG505)_Envelope',
  columns = c("EpistaticAntigenicity(structure)", 
              "EVEscape score", 
              "DMS score", 
              "Antigenicity(structure): individual amino acid", 
              "EpistaticSurfaceAccessibility(Sequence)",
              "EpistaticAntigenicity(sequence)"), 
  filter_df =  hiv_all_exp_mut,
  label_col = "label"
)


#try the antibody escape metrics:
hivbusy <- hivbg505_df_final %>% mutate(antibody_escape_seq = 1-`EpistaticAntigenicity(sequence)`, 
                                        antibody_escape_struct = 1-`EpistaticAntigenicity(structure)`,
                                        antibody_escape = 1/`Antigenicity(structure): individual amino acid`)
hivbg505_df_final %>% filter(position.x == 30) %>% mutate(new_antigenicity = 0.0751*`EpistaticAntigenicity(structure)`, ab_escape = 0.0751- new_antigenicity, ab_escape_scale = rescale(ab_escape, to = c(0,1))) %>% select(label, `EpistaticAntigenicity(structure)`, new_antigenicity, ab_escape, ab_escape_scale) %>% head(n=20)

hiv_busy2 <- hiv_baal3 %>% 
  mutate(new_antigenicity = `Antigenicity(structure): individual amino acid`*`EpistaticAntigenicity(structure)`, 
         ab_escape = `Antigenicity(structure): individual amino acid`^2 - new_antigenicity, 
         ab_escape_scale = rescale(ab_escape, to = c(0,1)),
         ab_escape2 = rescale(`Antigenicity(structure): individual amino acid` - `EpistaticAntigenicity(structure)`, to = c(0,1)),
         ab_escape3 = rescale(log(`Antigenicity(structure): individual amino acid` / `EpistaticAntigenicity(structure)`), to = c(0,1)), 
         ab_escape4 = rescale(log(`EpistaticAntigenicity(structure)` / `Antigenicity(structure): individual amino acid`), to = c(0,1)), 
         ab_escape5 = rescale(`Antigenicity(structure): individual amino acid` - `EpistaticAntigenicity(structure)`/`Antigenicity(structure): individual amino acid`, to = c(0,1))
         ) %>% 
  select(label, `EpistaticAntigenicity(structure)`, new_antigenicity, ab_escape, ab_escape_scale, ab_escape2,ab_escape3,ab_escape4, ab_escape5, "EpistaticAntigenicity(structure)",
         "EVEscape score", 
         "DMS score", 
         "Antigenicity(structure): individual amino acid", 
         "EpistaticSurfaceAccessibility(Sequence)",
         "EpistaticAntigenicity(sequence)"
        )

cor_matrix <- create_correlation_plot(
  df = hiv_busy2,
  virus_protein_name = 'HIV(BG505)_Envelope',
  columns = c("EpistaticAntigenicity(structure)",
              "EVEscape score", 
              "DMS score", 
              "Antigenicity(structure): individual amino acid", 
              "EpistaticSurfaceAccessibility(Sequence)",
              "EpistaticAntigenicity(sequence)", 
              "ab_escape2", 
              "ab_escape3", 
              "ab_escape4", 
              "ab_escape5"),
  filter_df =  hiv_all_exp_mut,
  label_col = "label"
)



hiv_busy3 <- hiv_baal3 %>% 
  mutate(#new_antigenicity = `Antigenicity(structure): individual amino acid`*`EpistaticAntigenicity(sequence)`, 
         #ab_escape = `Antigenicity(structure): individual amino acid`^2 - new_antigenicity, 
         #ab_escape_scale = rescale(ab_escape, to = c(0,1)),
         ab_escape2 = rescale(`Antigenicity_sequence` - `EpistaticAntigenicity(sequence)`, to = c(0,1)),
         ab_escape3 = rescale(log(`Antigenicity_sequence` / `EpistaticAntigenicity(sequence)`), to = c(0,1)), 
         ab_escape4 = rescale(log(`EpistaticAntigenicity(sequence)` / `Antigenicity_sequence`), to = c(0,1)), 
         ab_escape5 = rescale(`Antigenicity_sequence` - `EpistaticAntigenicity(sequence)`/`Antigenicity_sequence`, to = c(0,1))
  ) %>% 
  select(label, `EpistaticAntigenicity(structure)`, ab_escape2,ab_escape3,ab_escape4, ab_escape5, "EpistaticAntigenicity(structure)",
         "EVEscape score", 
         "DMS score", 
         "Antigenicity(structure): individual amino acid", 
         "EpistaticSurfaceAccessibility(Sequence)",
         "EpistaticAntigenicity(sequence)"
  )


hiv_busy4 <- hiv_baal3 %>% 
  mutate(Antigenicity_sequence = 1/(1+exp(-(Antigenicity_sequence - mean(Antigenicity_sequence, na.rm =T))/sd(Antigenicity_sequence, na.rm =T))), 
         `EpistaticAntigenicity(sequence)` = 1/(1+exp(-(`EpistaticAntigenicity(sequence)` - mean(`EpistaticAntigenicity(sequence)`, na.rm =T))/sd(`EpistaticAntigenicity(sequence)`, na.rm =T)))) %>% 
  mutate(#new_antigenicity = `Antigenicity(structure): individual amino acid`*`EpistaticAntigenicity(sequence)`, 
    #ab_escape = `Antigenicity(structure): individual amino acid`^2 - new_antigenicity, 
    #ab_escape_scale = rescale(ab_escape, to = c(0,1)),
    ab_escape2 = rescale(`Antigenicity_sequence` - `EpistaticAntigenicity(sequence)`, to = c(0,1)),
    ab_escape3 = rescale(log(`Antigenicity_sequence` / `EpistaticAntigenicity(sequence)`), to = c(0,1)), 
    ab_escape4 = rescale(log(`EpistaticAntigenicity(sequence)` / `Antigenicity_sequence`), to = c(0,1)), 
    ab_escape5 = rescale(`Antigenicity_sequence` - `EpistaticAntigenicity(sequence)`/`Antigenicity_sequence`, to = c(0,1)), 
    ab_escape6 = rescale(1 - log(`Antigenicity_sequence`/`EpistaticAntigenicity(sequence)`), to = c(0,1)), 
    ab_escape7 = rescale(1 - log(`EpistaticAntigenicity(sequence)`/`Antigenicity_sequence`), to = c(0,1)),
    ab_escape8 = rescale(1 - (`Antigenicity_sequence` - `EpistaticAntigenicity(sequence)`), to = c(0,1))
  ) %>% 
  select(label, `EpistaticAntigenicity(structure)`, 
         ab_escape2,ab_escape3,ab_escape4, ab_escape5, ab_escape6, ab_escape7, ab_escape8,
         "EpistaticAntigenicity(structure)",
         "EVEscape score", 
         "DMS score", 
         "Antigenicity(structure): individual amino acid", 
         "EpistaticSurfaceAccessibility(Sequence)",
         "EpistaticAntigenicity(sequence)", 
         "probability.x"
  )


hiv_busy5 <- hiv_baal3 %>% 
  mutate(Antigenicity_sequence = 1/(1+exp(-(Antigenicity_sequence - mean(Antigenicity_sequence, na.rm =T))/sd(Antigenicity_sequence, na.rm =T))), 
         `EpistaticAntigenicity(sequence)` = 1/(1+exp(-(`EpistaticAntigenicity(sequence)` - mean(`EpistaticAntigenicity(sequence)`, na.rm =T))/sd(`EpistaticAntigenicity(sequence)`, na.rm =T)))) %>% 
  mutate(#new_antigenicity = `Antigenicity(structure): individual amino acid`*`EpistaticAntigenicity(sequence)`, 
    #ab_escape = `Antigenicity(structure): individual amino acid`^2 - new_antigenicity, 
    #ab_escape_scale = rescale(ab_escape, to = c(0,1)),
  ) %>% 
  select(label, `EpistaticAntigenicity(structure)`, 
         starts_with('very_significant_antigenicity'),
         starts_with('significant_antigenicity'),
         "EpistaticAntigenicity(structure)",
         "EVEscape score", 
         "DMS score", 
         "Antigenicity(structure): individual amino acid", 
         "EpistaticSurfaceAccessibility(Sequence)",
         "EpistaticAntigenicity(sequence)", 
         "probability.x"
  )



hiv_busy6 <- hiv_baal3 %>% filter(position.x.x < 512) %>%  
  mutate(Antigenicity_sequence = 1/(1+exp(-(Antigenicity_sequence - mean(Antigenicity_sequence, na.rm =T))/sd(Antigenicity_sequence, na.rm =T))), 
         `EpistaticAntigenicity(sequence)` = 1/(1+exp(-(`EpistaticAntigenicity(sequence)` - mean(`EpistaticAntigenicity(sequence)`, na.rm =T))/sd(`EpistaticAntigenicity(sequence)`, na.rm =T)))) %>% 
  mutate(#new_antigenicity = `Antigenicity(structure): individual amino acid`*`EpistaticAntigenicity(sequence)`, 
    #ab_escape = `Antigenicity(structure): individual amino acid`^2 - new_antigenicity, 
    #ab_escape_scale = rescale(ab_escape, to = c(0,1)),
    ab_escape2 = rescale(`Antigenicity_sequence` - `EpistaticAntigenicity(sequence)`, to = c(0,1)),
    ab_escape3 = rescale(log(`Antigenicity_sequence` / `EpistaticAntigenicity(sequence)`), to = c(0,1)), 
    ab_escape4 = rescale(log(`EpistaticAntigenicity(sequence)` / `Antigenicity_sequence`), to = c(0,1)), 
    ab_escape5 = rescale(`Antigenicity_sequence` - `EpistaticAntigenicity(sequence)`/`Antigenicity_sequence`, to = c(0,1)), 
    ab_escape6 = rescale(1 - log(`Antigenicity_sequence`/`EpistaticAntigenicity(sequence)`), to = c(0,1)), 
    ab_escape7 = rescale(1 - log(`EpistaticAntigenicity(sequence)`/`Antigenicity_sequence`), to = c(0,1)),
    ab_escape8 = rescale(1 - (`Antigenicity_sequence` - `EpistaticAntigenicity(sequence)`), to = c(0,1))
  ) %>% 
  select(label, `EpistaticAntigenicity(structure)`, 
         ab_escape2,ab_escape3,ab_escape4, ab_escape5, ab_escape6, ab_escape7, ab_escape8,
         "EpistaticAntigenicity(structure)",
         "EVEscape score", 
         "DMS score", 
         "Antigenicity(structure): individual amino acid", 
         "EpistaticSurfaceAccessibility(Sequence)",
         "EpistaticAntigenicity(sequence)", 
         "probability.x"
  )


hiv_busy7 <- hiv_busy6 %>% 
  mutate(`EpistaticAntigenicity(sequence)` = case_when(label %in% n_gly$label ~ 1), .default = `EpistaticAntigenicity(sequence)`)


cor_matrix <- create_correlation_plot(
  df = hiv_busy3,
  virus_protein_name = 'HIV(BG505)_Envelope',
  columns = c("EpistaticAntigenicity(structure)",
              "EVEscape score", 
              "DMS score", 
              "Antigenicity(structure): individual amino acid", 
              "EpistaticSurfaceAccessibility(Sequence)",
              "EpistaticAntigenicity(sequence)", 
              "ab_escape2", 
              "ab_escape3", 
              "ab_escape4", 
              "ab_escape5"),
  filter_df =  hiv_all_exp_mut,
  label_col = "label"
)
# # Customize correlation plot appearance:
# cor_matrix <- create_correlation_plot(
#   df = flu_wsn33_final,
#   columns = c("significant_accessibility", "evescape", "max_escape_experiment"),
#   filter_df = flu_wsn33_true_antibody_escape,
#   label_col = "label",
#   order = "hclust",  # or "alphabet", "AOE", "FPC"
#   type = "upper",    # or "full", "lower"
#   diag = FALSE
# )

# 2. Then create rank comparison plots:
# Without filtering:
# results <- compare_ranks(
#   df = flu_wsn33_final,
#   columns = c("significant_accessibility", "evescape", "max_escape_experiment", "single_aa_antigenicity"),
#   reference_col = "significant_accessibility"
# )
#
# With filtering:
results <- compare_ranks(
  df = hivbg505_df_final,
  columns = c("EpistaticAntigenicity(structure)", "EVEscape score", "DMS score", 
              "Antigenicity(structure): individual amino acid", "EpistaticAntigenicity(sequence)"),
  reference_col = "EpistaticAntigenicity(structure)",
  filter_df = hiv_all_exp_mut,
  label_col = "label"
)
# 
# # View individual point plots
results$point_plots$`EVEscape score`
results$point_plots$
# 
# # View bar plot
results$bar_plot

# 
# # Combine all plots
 wrap_plots(results$point_plots, ncol = 2) / results$bar_plot + 
#   plot_layout(heights = c(2, 1))
# 
# # View comparison table
# results$comparison_table
results <- NULL
r <- c()  
r <- compare_ranks(
     df = hiv_baal3,
     columns = c(
                 'EpistaticAntigenicity(sequence)',
                 'Epistatic_Antigenicity_sequence',
                 'very_significant_antigenicity_median.y',
                 'Antigenicity_structure',
                 'Antigenicity_sequence',
                 'SurfaceAccessibility_sequence',
                 'DMS score', 
                 'EVEscape score'),
     reference_col = "Epistatic_Antigenicity_sequence",
     filter_df = hiv_all_exp_mut,
     label_col = "label"
   )


r$bar_plot
  
hw2 <- compare_ranks(
  df = hiv_busy7,
  columns = c('ab_escape2','ab_escape3', 'ab_escape4', 'ab_escape5',
              'ab_escape6', 'ab_escape7', 'ab_escape8',
              'EpistaticAntigenicity(structure)', 
              'EpistaticAntigenicity(sequence)', 
              'EpistaticSurfaceAccessibility(Sequence)', 
              'Antigenicity(structure): individual amino acid',
              'DMS score', 
              'EVEscape score', 
              'probability.x'),
  reference_col = "EpistaticAntigenicity(sequence)",
  #filter_df = hiv_all_exp_mut,
  #label_col = "label"
)


hw2$bar_plot

hw3 <- compare_ranks(
  df = hiv_busy5,
  columns = c('EpistaticAntigenicity(structure)',
              "very_significant_antigenicity_median.x",
              "significant_antigenicity_max",
              "very_significant_antigenicity_max",
              "significant_antigenicity_median",
              'EpistaticAntigenicity(sequence)', 
              'EpistaticSurfaceAccessibility(Sequence)', 
              'Antigenicity(structure): individual amino acid',
              'DMS score', 
              'EVEscape score', 
              'probability.x'),
  reference_col = "EpistaticAntigenicity(sequence)",
  filter_df = hiv_all_exp_mut,
  label_col = "label"
)

hw3$bar_plot 
hw3$comparison_table


hw4 <- compare_ranks_onebar(
  df = hiv_busy5,
  columns = c('EpistaticAntigenicity(structure)',
              "very_significant_antigenicity_median.x",
              "significant_antigenicity_max",
              "very_significant_antigenicity_max",
              "significant_antigenicity_median",
              'EpistaticAntigenicity(sequence)', 
              'EpistaticSurfaceAccessibility(Sequence)', 
              'Antigenicity(structure): individual amino acid',
              'DMS score', 
              'EVEscape score', 
              'probability.x'),
  reference_col = "EpistaticAntigenicity(sequence)",
  filter_df = hiv_all_exp_mut,
  label_col = "label"
)

hw4$bar_plot


hw5 <- compare_ranks_onebar(
  df = hiv_busy5,
  columns = c('EpistaticAntigenicity(structure)',
              'EpistaticAntigenicity(sequence)', 
              'EpistaticSurfaceAccessibility(Sequence)', 
              'DMS score', 
              'EVEscape score'),
  reference_col = "EpistaticAntigenicity(sequence)",
  filter_df = hiv_all_exp_mut,
  label_col = "label"
)

hw5$bar_plot

cor_matrix <- create_correlation_plot(
  df = hiv_busy5,
  virus_protein_name = 'HIV(BG505)_Envelope',
  columns = c('EpistaticAntigenicity(structure)',
              'EpistaticAntigenicity(sequence)', 
              
              'EpistaticSurfaceAccessibility(Sequence)', 
              'DMS score', 
              'EVEscape score', 
              'very_significant_antigenicity_sum'),
  filter_df =  hiv_all_exp_mut,
  label_col = "label"
)



#HIV N-glycosite labels issue:

hiv_baal3  %>% select(label, `EpistaticAntigenicity(sequence)`, `EVEscape score`) %>% drop_na(`EpistaticAntigenicity(sequence)`) %>% nrow()

hiv_baal3  %>% select(label, `EpistaticAntigenicity(sequence)`, `EVEscape score`) %>% drop_na(`EpistaticAntigenicity(sequence)`) %>% mutate(across(-label, ~rank(.x), .names = '{col}_rank')) %>%
  mutate(rank_diff = `EpistaticAntigenicity(sequence)_rank` - `EVEscape score_rank`) %>% filter(rank_diff >= 0) %>% summary()

hiv_baal3 %>% select(label, `EpistaticAntigenicity(sequence)`, `EVEscape score`) %>% drop_na(`EpistaticAntigenicity(sequence)`) %>% mutate(across(-label, ~rank(.x), .names = '{col}_rank')) %>%
  mutate(rank_diff = `EVEscape score_rank` - `EpistaticAntigenicity(sequence)_rank`) %>% filter(rank_diff > 0) %>% summary()


hiv_baal3 %>% select(label, `EpistaticAntigenicity(sequence)`, `EVEscape score`) %>% drop_na(`EpistaticAntigenicity(sequence)`) %>% 
  mutate(across(-label, ~rank(.x), .names = '{col}_rank')) %>% filter(label %in% unique(hiv_all_exp_mut$label)) %>% 
  mutate(rank_diff = `EpistaticAntigenicity(sequence)_rank` - `EVEscape score_rank`) %>% filter(rank_diff >= 0) %>% reframe(summary(rank_diff))

hiv_baal3 %>% select(label, `EpistaticAntigenicity(sequence)`, `EVEscape score`) %>% drop_na(`EpistaticAntigenicity(sequence)`) %>% 
  mutate(across(-label, ~rank(.x), .names = '{col}_rank')) %>%
  filter(label %in% unique(hiv_all_exp_mut$label)) %>% 
  mutate(rank_diff = `EVEscape score_rank` - `EpistaticAntigenicity(sequence)_rank`) %>% filter(rank_diff > 0) %>% reframe(summary(rank_diff))



hiv_baal3 %>% select(label, `EpistaticAntigenicity(sequence)`, `EVEscape score`, `N-glycan.x`) %>% drop_na(`EpistaticAntigenicity(sequence)`) %>% 
  mutate(across(-label, ~rank(.x), .names = '{col}_rank')) %>%
  # filter(label %in% n_gly$label) %>% 
  # filter(label %in% n_gly_exp_all) %>% 
  filter(label %in% n_nongly) %>% 
  #filter(label %in% unique(hiv_all_exp_mut$label)) %>% 
  #filter(label %in% label[which(hiv_baal3$`N-glycan.x` == TRUE)]) %>% 
  mutate(rank_diff = `EpistaticAntigenicity(sequence)_rank` - `EVEscape score_rank`) %>%
  filter(rank_diff >= 0) %>%
  nrow()
  # reframe(summary(rank_diff))


hiv_baal3 %>% select(label, `EpistaticAntigenicity(sequence)`, `EVEscape score`, `N-glycan.x`) %>% drop_na(`EpistaticAntigenicity(sequence)`) %>% 
  mutate(across(-label, ~rank(.x), .names = '{col}_rank')) %>%
  #filter(label %in% unique(hiv_all_exp_mut$label)) %>% 
  #filter(label %in% label[which(hiv_baal3$`N-glycan.x` == TRUE)]) %>% 
  # filter(label %in% n_gly$label) %>% 
  # filter(label %in% n_gly_exp_all) %>% 
  filter(label %in% n_nongly) %>% 
  mutate(rank_diff = `EVEscape score_rank` - `EpistaticAntigenicity(sequence)_rank`) %>%
  filter(rank_diff >= 0) %>%
  # nrow()
  reframe(summary(rank_diff))


  
hiv_baal3  %>% select(label, `EpistaticAntigenicity(sequence)`, `EVEscape score`) %>% drop_na(`EpistaticAntigenicity(sequence)`) %>% mutate(across(-label, ~rank(.x), .names = '{col}_rank')) %>% filter(label %in% hiv_baal3$label[which(hiv_baal3$`N-glycan.x` == TRUE)]) %>% 
  mutate(rank_diff = `EpistaticAntigenicity(sequence)_rank` - `EVEscape score_rank`) %>% filter(rank_diff >= 0) %>% summary()


## SC2

sc2_busy2 <- sc2_baal2 %>% distinct(label, .keep_all = T) %>% 
  mutate(Epistatic_Antigenicity_structure = replace_na(Epistatic_Antigenicity_structure, 0.47)) %>% 
  mutate(new_antigenicity = Antigenicity_structure*Epistatic_Antigenicity_structure, 
         ab_escape = Antigenicity_structure^2 - new_antigenicity, 
         ab_escape_scale = rescale(ab_escape, to = c(0,1)),
         ab_escape2 = rescale(Antigenicity_structure - Epistatic_Antigenicity_structure, to = c(0,1)),
         ab_escape3 = rescale(log(Antigenicity_structure / Epistatic_Antigenicity_structure), to = c(0,1)), 
         ab_escape4 = rescale(log(Epistatic_Antigenicity_structure / Antigenicity_structure), to = c(0,1)), 
         ab_escape5 = rescale(Antigenicity_structure - Epistatic_Antigenicity_structure/Antigenicity_structure, to = c(0,1))
  ) %>% 
  select(label, Epistatic_Antigenicity_structure, Antigenicity_structure, new_antigenicity, ab_escape, ab_escape_scale, ab_escape2,ab_escape3,ab_escape4, ab_escape5,
         "EVEscape score", 
         "DMS score", 
         "Epistatic_SurfaceAccessibility_sequence",
         "Epistatic_Antigenicity_sequence",
         "Antigenicity_sequence"
  )

nrow(sc2_busy2)


cor_matrix <- create_correlation_plot(
  df = sc2_busy2,
  virus_protein_name = 'sc2_spike',
  columns = c("Epistatic_Antigenicity_structure",
              "EVEscape score", 
              "DMS score", 
              "Antigenicity_structure", 
              "Epistatic_SurfaceAccessibility_sequence",
              "Epistatic_Antigenicity_sequence", 
              "ab_escape2", 
              "ab_escape3", 
              "ab_escape4", 
              "ab_escape5"),
  filter_df =  sc2_tre_high_antibody_escape,
  label_col = "label"
)

#convert Epistatic_Antigenicity_strcuture to Epistatic_Antigenicity_sequence

1/(1+exp(-((x-mean(df$a)/sd(df$a)))))
sc2_busy3 <- sc2_baal2 %>% distinct(label, .keep_all = T) %>% 
  mutate(new_antigenicity = Antigenicity_sequence*Epistatic_Antigenicity_sequence, 
         ab_escape = Antigenicity_sequence^2 - new_antigenicity, 
         ab_escape_scale = rescale(ab_escape, to = c(0,1)),
         ab_escape2 = rescale(Antigenicity_sequence - Epistatic_Antigenicity_sequence, to = c(0,1)),
         ab_escape3 = rescale(log(Antigenicity_sequence / Epistatic_Antigenicity_sequence), to = c(0,1)), 
         ab_escape4 = rescale(log(Epistatic_Antigenicity_sequence / Antigenicity_sequence), to = c(0,1)), 
         ab_escape5 = rescale(Antigenicity_sequence - Epistatic_Antigenicity_sequence/Antigenicity_sequence, to = c(0,1))
  ) %>% 
  select(label, Epistatic_Antigenicity_sequence, Antigenicity_sequence, 
         new_antigenicity, ab_escape, 
         ab_escape_scale, ab_escape2,
         ab_escape3,ab_escape4, ab_escape5,
         "EVEscape score", 
         "DMS score", 
         "Epistatic_SurfaceAccessibility_sequence",
         "Epistatic_Antigenicity_structure",
         "Antigenicity_structure"
  )

cor_matrix <- create_correlation_plot(
  df = sc2_busy3,
  virus_protein_name = 'sc2_spike',
  columns = c("Epistatic_Antigenicity_sequence",
              "EVEscape score", 
              "DMS score", 
              "Antigenicity_sequence", 
              "Epistatic_SurfaceAccessibility_sequence",
              "Epistatic_Antigenicity_structure", 
              "ab_escape2", 
              "ab_escape3", 
              "ab_escape4", 
              "ab_escape5"),
  filter_df =  sc2_tre_antibody_escape,
  label_col = "label"
)

sc2_busy4 <- sc2_baal2 %>% distinct(label, .keep_all = T) %>% 
  mutate(Antigenicity_sequence = 1/(1+exp(-(Antigenicity_sequence - mean(Antigenicity_sequence))/sd(Antigenicity_sequence))), 
         Epistatic_Antigenicity_sequence = 1/(1+exp(-(Epistatic_Antigenicity_sequence - mean(Epistatic_Antigenicity_sequence))/sd(Epistatic_Antigenicity_sequence)))) %>% 
  mutate(new_antigenicity = Antigenicity_sequence*Epistatic_Antigenicity_sequence, 
         ab_escape = Antigenicity_sequence^2 - new_antigenicity, 
         ab_escape_scale = rescale(ab_escape, to = c(0,1)),
         ab_escape2 = rescale(Antigenicity_sequence - Epistatic_Antigenicity_sequence, to = c(0,1)),
         ab_escape3 = rescale(log(Antigenicity_sequence / Epistatic_Antigenicity_sequence), to = c(0,1)), 
         ab_escape4 = rescale(log(Epistatic_Antigenicity_sequence / Antigenicity_sequence), to = c(0,1)), 
         ab_escape5 = rescale(Antigenicity_sequence - Epistatic_Antigenicity_sequence/Antigenicity_sequence, to = c(0,1))
  ) %>% 
  select(label, Epistatic_Antigenicity_sequence, Antigenicity_sequence, 
         new_antigenicity, ab_escape, 
         ab_escape_scale, ab_escape2,
         ab_escape3,ab_escape4, ab_escape5,
         "EVEscape score", 
         "DMS score", 
         "Epistatic_SurfaceAccessibility_sequence",
         "Epistatic_Antigenicity_structure",
         "Antigenicity_structure"
  )
cor_matrix <- create_correlation_plot(
  df = sc2_busy4,
  virus_protein_name = 'sc2_spike',
  columns = c("Epistatic_Antigenicity_sequence",
              "EVEscape score", 
              "DMS score", 
              "Antigenicity_sequence", 
              "Epistatic_SurfaceAccessibility_sequence",
              "Epistatic_Antigenicity_structure", 
              "ab_escape2", 
              "ab_escape3", 
              "ab_escape4", 
              "ab_escape5"),
  filter_df =  sc2_tre_high_antibody_escape,
  label_col = "label"
)



sc2_busy5 <- sc2_baal4 %>% distinct(label, .keep_all = T) %>% 
  mutate(Antigenicity_sequence = 1/(1+exp(-(Antigenicity_sequence - mean(Antigenicity_sequence))/sd(Antigenicity_sequence))), 
         Epistatic_Antigenicity_sequence = 1/(1+exp(-(Epistatic_Antigenicity_sequence - mean(Epistatic_Antigenicity_sequence))/sd(Epistatic_Antigenicity_sequence)))) %>% 
  mutate(new_antigenicity = Antigenicity_sequence*Epistatic_Antigenicity_sequence, 
         ab_escape = Antigenicity_sequence^2 - new_antigenicity, 
         ab_escape_scale = rescale(ab_escape, to = c(0,1)),
         ab_escape2 = rescale(Antigenicity_sequence - Epistatic_Antigenicity_sequence, to = c(0,1)),
         ab_escape3 = rescale(log(Antigenicity_sequence / Epistatic_Antigenicity_sequence), to = c(0,1)), 
         ab_escape4 = rescale(log(Epistatic_Antigenicity_sequence / Antigenicity_sequence), to = c(0,1)), 
         ab_escape5 = rescale(Antigenicity_sequence - Epistatic_Antigenicity_sequence/Antigenicity_sequence, to = c(0,1))
  ) %>% 
  select(position.x, label, Epistatic_Antigenicity_sequence, Antigenicity_sequence, 
         new_antigenicity, ab_escape, 
         ab_escape_scale, ab_escape2,
         ab_escape3,ab_escape4, ab_escape5,
         probability,
         semantic_score,
         "EVEscape score", 
         "DMS score", 
         "Epistatic_SurfaceAccessibility_sequence",
         "Epistatic_Antigenicity_structure",
         "Antigenicity_structure"
  )

summary(sc2_busy5)



sc2_busy6 <- sc2_rbd_baal2 %>% distinct(label, .keep_all = T) %>% 
  mutate(#Antigenicity_sequence = 1/(1+exp(-(Antigenicity_sequence - mean(Antigenicity_sequence))/sd(Antigenicity_sequence))), 
         Epistatic_Antigenicity_sequence = 1/(1+exp(-(Epistatic_Antigenicity_sequence - mean(Epistatic_Antigenicity_sequence))/sd(Epistatic_Antigenicity_sequence)))) %>% 
  mutate(new_antigenicity = Antigenicity_sequence*Epistatic_Antigenicity_sequence, 
         ab_escape = Antigenicity_sequence^2 - new_antigenicity, 
         ab_escape_scale = rescale(ab_escape, to = c(0,1)),
         ab_escape2 = rescale(Antigenicity_sequence - Epistatic_Antigenicity_sequence, to = c(0,1)),
         ab_escape3 = rescale(log(Antigenicity_sequence / Epistatic_Antigenicity_sequence), to = c(0,1)), 
         ab_escape4 = rescale(log(Epistatic_Antigenicity_sequence / Antigenicity_sequence), to = c(0,1)), 
         ab_escape5 = rescale(Antigenicity_sequence - Epistatic_Antigenicity_sequence/Antigenicity_sequence, to = c(0,1))
  ) %>% 
  select(position.x, label, Epistatic_Antigenicity_sequence, Antigenicity_sequence, 
         new_antigenicity, ab_escape, 
         ab_escape_scale, ab_escape2,
         ab_escape3,ab_escape4, ab_escape5,
         probability,
         semantic_score,
         "EVEscape score", 
         "DMS score", 
         "Epistatic_SurfaceAccessibility_sequence",
         "Epistatic_Antigenicity_structure",
         "Antigenicity_structure", 
         "antibody_escape_score1",
         "antibody_escape_score2", 
         "antibody_escape_score3"
  )

sc2_busy6 <- sc2_rbd_baal2 %>% 
  group_by(label) %>% mutate(dms_score = min(`DMS score`, na.rm = T)) %>% select(!`DMS score`) %>% distinct(label, .keep_all = T) %>% 
  ungroup() %>% 
  mutate(#Antigenicity_sequence = 1/(1+exp(-(Antigenicity_sequence - mean(Antigenicity_sequence))/sd(Antigenicity_sequence))), 
    Epistatic_Antigenicity_sequence = 1/(1+exp(-(Epistatic_Antigenicity_sequence - mean(Epistatic_Antigenicity_sequence))/sd(Epistatic_Antigenicity_sequence)))) %>% 
  mutate(new_antigenicity = Antigenicity_sequence*Epistatic_Antigenicity_sequence, 
         #ab_escape = Antigenicity_sequence^2 - new_antigenicity, 
         #ab_escape_scale = rescale(ab_escape, to = c(0,1)),
         ab_escape2 = rescale(Antigenicity_sequence - Epistatic_Antigenicity_sequence, to = c(0,1)),
         ab_escape3 = rescale(log(Antigenicity_sequence / Epistatic_Antigenicity_sequence), to = c(0,1)), 
         ab_escape4 = rescale(log(Epistatic_Antigenicity_sequence / Antigenicity_sequence), to = c(0,1)), 
         ab_escape5 = rescale(Antigenicity_sequence - Epistatic_Antigenicity_sequence/Antigenicity_sequence, to = c(0,1))
  ) %>% 
  select(position.x, label, Epistatic_Antigenicity_sequence, Antigenicity_sequence, 
         new_antigenicity,
         ab_escape2,
         ab_escape3,ab_escape4, ab_escape5,
         probability,
         semantic_score,
         "EVEscape score", 
         dms_score, 
         "Epistatic_SurfaceAccessibility_sequence",
         "Epistatic_Antigenicity_structure",
         "Antigenicity_structure", 
         "antibody_escape_score1",
         "antibody_escape_score2", 
         "antibody_escape_score3"
  )

cor_matrix <- create_correlation_plot(
  df = sc2_busy6,
  virus_protein_name = 'sc2_spike',
  columns = c("Epistatic_Antigenicity_sequence",
              "EVEscape score", 
              "DMS score", 
              "Antigenicity_sequence", 
              "Epistatic_SurfaceAccessibility_sequence",
              "Epistatic_Antigenicity_structure", 
              "ab_escape2", 
              "ab_escape3", 
              "ab_escape4", 
              "ab_escape5", 
              'probability', 
              'semantic_score', 
              "antibody_escape_score1",
              "antibody_escape_score2", 
              "antibody_escape_score3"
              ),
  filter_df =  sc2_tre_high_antibody_escape,
  label_col = "label"
)

sc2_bloommf <- sc2_dms_final %>% group_by(label) %>% 
  filter(!is.na(max_escape_experiment)) %>% mutate(dms_score= min(max_escape_experiment)) %>% 
  ungroup() %>% 
  # distinct(label, .keep_all = T) %>% 
  distinct(label, .keep_all = T) 
  #filter(position == 484, ref == 'E')
  # print(n=50)
  

# filter(!is.na(`DMS score`)) %>% 
#   group_by(label) %>% 
#   mutate(dms_score= min(`DMS score`)) %>% 
#   ungroup() %>% 
#   distinct(label, .keep_all = T) 

sc2_busy7 <- sc2_baal4 %>%
  left_join(sc2_bloommf %>% select(label, dms_score), by = 'label') %>%
  drop_na(Antigenicity_sequence, Epistatic_Antigenicity_sequence) %>% 
  mutate(Antigenicity_sequence = 1/(1+exp(-(Antigenicity_sequence - mean(Antigenicity_sequence))/sd(Antigenicity_sequence))), 
         Epistatic_Antigenicity_sequence = 1/(1+exp(-(Epistatic_Antigenicity_sequence - mean(Epistatic_Antigenicity_sequence))/sd(Epistatic_Antigenicity_sequence))),
         epiant = 1/(1+exp(-(very_significant_accessibility_median.y - mean(very_significant_accessibility_median.y))/sd(very_significant_accessibility_median.y)))
         ) %>%
  mutate(new_antigenicity = Antigenicity_sequence*Epistatic_Antigenicity_sequence, 
         #ab_escape = Antigenicity_sequence^2 - new_antigenicity, 
         #ab_escape_scale = rescale(ab_escape, to = c(0,1)),
         ab_escape2 = rescale(Antigenicity_sequence - Epistatic_Antigenicity_sequence, to = c(0,1)),
         ab_escape3 = rescale(log(Antigenicity_sequence / Epistatic_Antigenicity_sequence), to = c(0,1)), 
         ab_escape4 = rescale(log(Epistatic_Antigenicity_sequence / Antigenicity_sequence), to = c(0,1)), 
         ab_escape5 = rescale(Antigenicity_sequence - Epistatic_Antigenicity_sequence/Antigenicity_sequence, to = c(0,1))
  ) %>% 
  select(position.x, label, 
         epiant,
         Epistatic_Antigenicity_sequence, 
         Antigenicity_sequence, 
         new_antigenicity,
         ab_escape2,
         ab_escape3,ab_escape4, ab_escape5,
         probability,
         semantic_score,
         "EVEscape score", 
         dms_score, 
         "Epistatic_SurfaceAccessibility_sequence",
         "Epistatic_Antigenicity_structure",
         "Antigenicity_structure", 
         
  )


sc2_busy7 %>% filter(position.x == 484) %>% select(dms_score, label)
# # Customize correlation plot appearance:
# cor_matrix <- create_correlation_plot(
#   df = flu_wsn33_final,
#   columns = c("significant_accessibility", "evescape", "max_escape_experiment"),
#   filter_df = flu_wsn33_true_antibody_escape,
#   label_col = "label",
#   order = "hclust",  # or "alphabet", "AOE", "FPC"
#   type = "upper",    # or "full", "lower"
#   diag = FALSE
# )

# 2. Then create rank comparison plots:
# Without filtering:
# results <- compare_ranks(
#   df = flu_wsn33_final,
#   columns = c("significant_accessibility", "evescape", "max_escape_experiment", "single_aa_antigenicity"),
#   reference_col = "significant_accessibility"
# )
#


# With filtering:
results <- compare_ranks(
 df= sc2_fullspike_final_v2,
  columns = c('significant_accessibility', 'single_aa_antigenicity', 'evescape.y', 
              'max_escape_experiment', 'significant_antigenicity_median', 
              'significant_accessibility_median'), 
  reference_col = "significant_accessibility_median",
  filter_df = sc2_tre_antibody_escape,
  label_col = "label")
# 
# # View individual point plots
# results$point_plots$evescape_rank
# results$point_plots$max_escape_experiment_rank
# 
# # View bar plot
results$bar_plot

# 
# # Combine all plots
# wrap_plots(results$point_plots, ncol = 2) / results$bar_plot + 
#   plot_layout(heights = c(2, 1))
# 
# # View comparison table
# results$comparison_table


sw <- c()
sw <- compare_ranks(
  df = sc2_baal2 %>% distinct(label, .keep_all = T),
  columns = c('Epistatic_Antigenicity_structure', 
              'Epistatic_Antigenicity_sequence', 
              'Epistatic_SurfaceAccessibility_sequence', 
              'Antigenicity_structure',
              'Antigenicity_sequence',
              'SurfaceAccessibility_sequence',
              'DMS score', 
              'EVEscape score'),
  reference_col = "Epistatic_SurfaceAccessibility_sequence",
  filter_df = sc2_untrue_antibody_escape,
  label_col = "label"
)
sw$bar_plot
sw


sc2_tre_antibody_escape$label
sc2_baal2 %>% distinct(label) %>% pull(label)
sc2_untrue_antibody_escape <- data.frame(label = setdiff (sc2_baal3, sc2_tre_antibody_escape$label))


sw2 <- c()
sw2 <- compare_ranks(
  df = sc2_busy5,
  columns = c('Epistatic_Antigenicity_structure', 
              'Epistatic_Antigenicity_sequence', 
              'Epistatic_SurfaceAccessibility_sequence',
              'ab_escape2', 
              'ab_escape3',
              'ab_escape4',
              'ab_escape5',
              'Antigenicity_structure',
              'Antigenicity_sequence',
              'DMS score', 
              'EVEscape score', 
              'probability', 
              'semantic_score'),
  reference_col = "Epistatic_Antigenicity_sequence",
  filter_df = sc2_tre_high_antibody_escape,
  label_col = "label"
)

sw2$bar_plot


sw3 <- c()
sw3 <- compare_ranks(
  df = sc2_busy6,
  columns = c('Epistatic_Antigenicity_structure', 
              'Epistatic_Antigenicity_sequence', 
              'Epistatic_SurfaceAccessibility_sequence',
              'ab_escape2', 
              'ab_escape3',
              'ab_escape4',
              'ab_escape5',
              'Antigenicity_structure',
              'Antigenicity_sequence',
              'DMS score', 
              'EVEscape score', 
              'probability', 
              'semantic_score', 
              'antibody_escape_score1',
              'antibody_escape_score2',
              'antibody_escape_score3'),
  reference_col = "antibody_escape_score1",
  filter_df = sc2_tre_high_antibody_escape,
  label_col = "label"
)

sw3$bar_plot



cor_matrix <- create_correlation_plot(
  df = sc2_busy6,
  virus_protein_name = 'sc2_spike',
  columns = c("Epistatic_Antigenicity_sequence",
              "EVEscape score", 
              "dms_score",
              "Antigenicity_sequence", 
              "Epistatic_SurfaceAccessibility_sequence",
              "Epistatic_Antigenicity_structure", 
              "ab_escape2", 
              "ab_escape3", 
              "ab_escape4", 
              "ab_escape5", 
              'probability', 
              'semantic_score', 
              "antibody_escape_score1",
              "antibody_escape_score2", 
              "antibody_escape_score3"
  ),
  filter_df =  sc2_tre_high_antibody_escape,
  label_col = "label"
)
 

sw4 <- c()
sw4 <- compare_ranks(
  df = sc2_busy6,
  columns = c('Epistatic_Antigenicity_structure', 
              'Epistatic_Antigenicity_sequence', 
              'Epistatic_SurfaceAccessibility_sequence',
              'ab_escape2', 
              'ab_escape3',
              'ab_escape4',
              'ab_escape5',
              'Antigenicity_structure',
              'Antigenicity_sequence',
              'dms_score', 
              'EVEscape score', 
              'probability', 
              'semantic_score', 
              'antibody_escape_score1',
              'antibody_escape_score2',
              'antibody_escape_score3'),
  reference_col = "Epistatic_Antigenicity_sequence",
  filter_df = sc2_tre_high_antibody_escape,
  label_col = "label"
)

sw4$bar_plot

sw5 <- c()
sw5 <- compare_ranks(
  df = sc2_busy7,
  columns = c('Epistatic_Antigenicity_structure', 
              'Epistatic_Antigenicity_sequence', 
              'Epistatic_SurfaceAccessibility_sequence',
              'ab_escape2', 
              'ab_escape3',
              'ab_escape4',
              'ab_escape5',
              'Antigenicity_structure',
              'Antigenicity_sequence',
              'dms_score', 
              'EVEscape score', 
              'probability', 
              'semantic_score', 
              'epiant'
            ),
  reference_col = "Epistatic_Antigenicity_sequence",
  filter_df = sc2_tre_high_antibody_escape,
  label_col = "label"
)
sw5$bar_plot

sw5 <- compare_ranks_onebar(
  df = sc2_busy7,
  columns = c('Epistatic_Antigenicity_structure', 
              'Epistatic_Antigenicity_sequence', 
              'Epistatic_SurfaceAccessibility_sequence',
              'dms_score', 
              'EVEscape score'
  ),
  reference_col = "Epistatic_Antigenicity_sequence",
  filter_df = sc2_tre_antibody_escape,
  label_col = "label"
)
sw5$bar_plot

##FLU


# 1. Create correlation plot first:
cor_matrix <- create_correlation_plot(
  df = flu_wsn33_final,
  virus_protein_name = "FLU_HA",
  columns = c("significant_accessibility", "evescape", "max_escape_experiment", "single_aa_antigenicity"),
  filter_df = flu_wsn33_true_antibody_escape,
  label_col = "label"
)

cor_matrix <- create_correlation_plot(
  df = fl_new,
  columns = c("significant_accessibility",'significant_accessibility_median', 'significant_antigenicity_median', "evescape", "max_escape_experiment", "single_aa_antigenicity"),
  filter_df = flu_wsn33_true_antibody_escape,
  label_col = "label"
)
#
# # Customize correlation plot appearance:
# cor_matrix <- create_correlation_plot(
#   df = flu_wsn33_final,
#   columns = c("significant_accessibility", "evescape", "max_escape_experiment"),
#   filter_df = flu_wsn33_true_antibody_escape,
#   label_col = "label",
#   order = "hclust",  # or "alphabet", "AOE", "FPC"
#   type = "upper",    # or "full", "lower"
#   diag = FALSE
# )

# 2. Then create rank comparison plots:
# Without filtering:
# results <- compare_ranks(
#   df = flu_wsn33_final,
#   columns = c("significant_accessibility", "evescape", "max_escape_experiment", "single_aa_antigenicity"),
#   reference_col = "significant_accessibility"
# )
#
# With filtering:
results <- compare_ranks(
  df = flu_wsn33_final,
  columns = c("significant_accessibility", "evescape", "max_escape_experiment", "single_aa_antigenicity"),
  reference_col = "significant_accessibility",
  filter_df = flu_wsn33_true_antibody_escape,
  label_col = "label"
)
# 
# # View individual point plots
results$point_plots$evescape_rank
results$point_plots$max_escape_experiment_rank
# 
# # View bar plot
results$bar_plot
# 
# # Combine all plots
# wrap_plots(results$point_plots, ncol = 2) / results$bar_plot + 
#   plot_layout(heights = c(2, 1))
# 
# # View comparison table
# results$comparison_table

fw <- compare_ranks(
  df = flu_baal2,
  columns = c('Epistatic_Antigenicity_structure', 
              'Epistatic_Antigenicity_sequence', 
              'Epistatic_SurfaceAccessibility_sequence', 
              'Antigenicity_structure',
              'Antigenicity_sequence',
              'SurfaceAccessibility_sequence',
              'DMS score', 
              'EVEscape score'),
  reference_col = "Epistatic_Antigenicity_structure",
  filter_df = flu_wsn33_true_antibody_escape,
  label_col = "label"
)
fw$bar_plot


flu_busy2 <- flu_baal2 %>% drop_na(Antigenicity_sequence) %>% drop_na(Epistatic_Antigenicity_sequence) %>% 
  mutate(Antigenicity_sequence = 1/(1+exp(-(Antigenicity_sequence - mean(Antigenicity_sequence, na.rm = T))/sd(Antigenicity_sequence, na.rm = T))), 
         Epistatic_Antigenicity_sequence = 1/(1+exp(-(Epistatic_Antigenicity_sequence - mean(Epistatic_Antigenicity_sequence, na.rm = T))/sd(Epistatic_Antigenicity_sequence, na.rm = T)))) %>% 
  mutate(#new_antigenicity = `Antigenicity(structure): individual amino acid`*`EpistaticAntigenicity(sequence)`, 
    #ab_escape = `Antigenicity(structure): individual amino acid`^2 - new_antigenicity, 
    #ab_escape_scale = rescale(ab_escape, to = c(0,1)),
    ab_escape2 = rescale(Antigenicity_sequence - Epistatic_Antigenicity_sequence, to = c(0,1)),
    ab_escape3 = rescale(log(Antigenicity_sequence / Epistatic_Antigenicity_sequence), to = c(0,1)), 
    ab_escape4 = rescale(log(Epistatic_Antigenicity_sequence / Antigenicity_sequence), to = c(0,1)), 
    ab_escape5 = rescale(Antigenicity_sequence - Epistatic_Antigenicity_sequence/Antigenicity_sequence, to = c(0,1)),
    ab_escape6 = rescale(1- (Epistatic_Antigenicity_sequence / Antigenicity_sequence)),
    ab_escape7 = rescale(1 - (Antigenicity_sequence - Epistatic_Antigenicity_sequence))
  ) %>% 
  select(label, probability, semantic_score,
         Epistatic_Antigenicity_structure, 
         Antigenicity_structure, 
         ab_escape2,ab_escape3,ab_escape4, 
         ab_escape5,ab_escape6, ab_escape7,
         "EVEscape score", 
         "DMS score", 
         "Epistatic_SurfaceAccessibility_sequence",
         "Epistatic_Antigenicity_sequence",
         Antigenicity_sequence,
         SurfaceAccessibility_sequence
  )



cor_matrix <- create_correlation_plot(
  df = flu_busy2,
  virus_protein_name = 'flu_HA',
  columns = c("Epistatic_Antigenicity_sequence",
              "EVEscape score", 
              "DMS score",
              'probability',
              'semantic_score',
              "Antigenicity_sequence", 
              "Epistatic_SurfaceAccessibility_sequence",
              "Epistatic_Antigenicity_structure", 
              "ab_escape2", 
              "ab_escape3", 
              "ab_escape4", 
              "ab_escape5", 
              "ab_escape6",
              "ab_escape7"),
  filter_df =  flu_wsn33_true_antibody_escape,
  label_col = "label"
)

fw2 <- c()
fw2 <- compare_ranks(
  df = flu_busy2,
  columns = c('Epistatic_Antigenicity_structure', 
              'Epistatic_Antigenicity_sequence', 
              'Epistatic_SurfaceAccessibility_sequence',
              'ab_escape2', 
              'ab_escape3',
              'ab_escape4',
              'ab_escape5',
              'ab_escape6',
              'ab_escape7',
              'Antigenicity_structure',
              'Antigenicity_sequence',
              'DMS score', 
              'EVEscape score', 
              'probability', 
              'semantic_score'),
  reference_col = "Epistatic_Antigenicity_sequence",
  filter_df = flu_wsn33_true_antibody_escape,
  label_col = "label"
)
fw2$bar_plot


fw3 <- compare_ranks_onebar(
  df = flu_busy2,
  columns = c('Epistatic_Antigenicity_structure', 
              'Epistatic_Antigenicity_sequence', 
              'Epistatic_SurfaceAccessibility_sequence',
              'DMS score', 
              'EVEscape score'),
  reference_col = "Epistatic_Antigenicity_sequence",
  filter_df = flu_wsn33_true_antibody_escape,
  label_col = "label"
)
fw3$bar_plot

library(tidyverse)

# Create sample dataset
test_df <- data.frame(
  label = c("Apple", "Banana", "Cherry", "Date", "Elderberry", 
            "Fig", "Grape", "Honeydew", "Kiwi", "Lemon"),
  a = c(0.85, 0.72, 0.91, 0.65, 0.78, 
        0.82, 0.69, 0.88, 0.74, 0.80),
  b = c(0.75, 0.68, 0.88, 0.70, 0.82, 
        0.76, 0.72, 0.85, 0.79, 0.83),
  c = c(0.80, 0.71, 0.89, 0.68, 0.80, 
        0.78, 0.70, 0.86, 0.77, 0.81),
  d = c(0.77, 0.69, 0.90, 0.66, 0.79, 
        0.79, 0.71, 0.87, 0.76, 0.82)
)

print(test_df)


test_df %>% 
  mutate(across(-label, ~ rank(-.x, ties.method = "min"), .names = "{.col}_rank")) %>%
  select(label, ends_with('_rank'))

# Test the function
results <- compare_ranks(
  df = test_df,
  columns = c("a", "b", "c", "d"),
  reference_col = "a"
)

# View results
results$bar_plot
results$point_plots[[1]]
results$comparison_table

compare_ranks(
  df = sc2_epi_rank_check,
  columns = c('Antigenicity_sequence', 'max_epiant', 'med_epiant'),
  reference_col = "max_epiant",
  filter_df = scc2,
  label_col = 'position.x'
)


compare_ranks(
  df = sc2_epi_rank_check,
  columns = c('Accessibility', 'max_epi_access'),
  reference_col = "max_epi_access",
  filter_df = scc2,
  label_col = 'position'
)


#new ranking and correlation plots:

# Main function to compare ranks
compare_ranks_new <- function(df, columns, reference_col, filter_df = NULL, label_col = "label") {
  
  # Select and rank the specified columns
  df_ranked <- df %>% 
    select(all_of(c(label_col, columns))) %>%
    mutate(across(-all_of(label_col), ~ rank(.x, ties.method = "min"), .names = "{.col}_rank")) %>%
    select(all_of(label_col), ends_with('_rank'))
  
  # Apply filtering if filter_df is provided
  if (!is.null(filter_df)) {
    filter_labels <- filter_df[[label_col]]
    df_ranked <- df_ranked %>% 
      filter(.data[[label_col]] %in% filter_labels)
  }
  
  # Remove label column after filtering
  df_ranked <- df_ranked %>% select(-all_of(label_col))
  
  # Get total count for subtitle
  total_count <- nrow(df_ranked)
  
  # Reorder so reference column is first
  ref_col_name <- paste0(reference_col, "_rank")
  df_ranked <- df_ranked %>% select(all_of(ref_col_name), everything())
  
  # Generate point plots
  plot_list <- list()
  other_cols <- setdiff(colnames(df_ranked), ref_col_name)
  
  for (col in other_cols) {
    p <- ggplot(df_ranked, aes(x = .data[[ref_col_name]], y = .data[[col]])) + 
      geom_point(alpha = 0.6, size = 2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(x = gsub("_rank", "", ref_col_name),
           y = gsub("_rank", "", col),
           title = paste("Rank Comparison:", gsub("_rank", "", col))) +
      theme_minimal() +
      theme(plot.title = element_text(size = 12, face = "bold"))
    
    plot_list[[col]] <- p
  }
  
  # Generate bar plot comparing reference to all others
  comparison_data <- data.frame()
  
  for (col in other_cols) {
    higher_count <- sum(df_ranked[[ref_col_name]] >= df_ranked[[col]], na.rm = TRUE)
    
    comparison_data <- rbind(
      comparison_data,
      data.frame(
        Method = gsub("_rank", "", col),
        higher_count = higher_count,
        total = total_count,
        higher_pct = (higher_count / total_count) * 100
      )
    )
  }
  
  # Create bar plot with counts and percentages
  bar_plot <- ggplot(comparison_data, aes(x = Method, y = higher_pct)) +
    geom_bar(stat = 'identity', fill = "#E69F00") +
    geom_text(aes(label = sprintf("%.1f%% (%d/%d)", higher_pct, higher_count, total)), 
              vjust = -0.5, size = 6.5) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "black", linewidth = 1) +
    labs(
      # subtitle = paste('Higher Ranks in Epistatic Antigenicity \ncompared to other antigenicity metrics'),
         # subtitle = paste("Total no. of antibody-escape amino-acid substitutions: N =", total_count, 
         #                  "\n(supported by the literature and experimental evidence)"),
         x = 'Antigenicity Metrics',
         y = paste0('Number of Substitutions in Percentage (N = ', total_count, ')')) +
    ylim(0, 105) +  # Extra space for labels
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = 'bold'),
          axis.text.y = element_text(size = 18),
          axis.title = element_text(size = 22, face = "bold"),
          #plot.title = element_text(size = 14, face = "bold"),
          #plot.subtitle = element_text(size = 10)
          )
  
  # Return all plots
  return(list(
    point_plots = plot_list,
    bar_plot = bar_plot,
    comparison_table = comparison_data
  ))
}

fw4 <- compare_ranks_new(
  df = flu_busy2,
  columns = c('Epistatic_Antigenicity_structure', 
              'Epistatic_Antigenicity_sequence', 
              'Epistatic_SurfaceAccessibility_sequence',
              'DMS score', 
              'EVEscape score'),
  reference_col = "Epistatic_Antigenicity_sequence",
  filter_df = flu_wsn33_true_antibody_escape,
  label_col = "label"
)
fw4$bar_plot

fw5 <- compare_ranks_new(
  df = flu_busy2,
  columns = c('Antigenicity_structure', 
              'Antigenicity_sequence', 
              'Epistatic_Antigenicity_sequence',
              'DMS score', 
              'EVEscape score'),
  reference_col = "Epistatic_Antigenicity_sequence",
  filter_df = flu_wsn33_true_antibody_escape,
  label_col = "label"
)

fw5$bar_plot



sw5 <- compare_ranks_new(
  df = sc2_busy7,
  columns = c('Epistatic_Antigenicity_structure', 
              'Epistatic_Antigenicity_sequence', 
              'Epistatic_SurfaceAccessibility_sequence',
              'dms_score', 
              'EVEscape score'
  ),
  reference_col = "Epistatic_Antigenicity_sequence",
  filter_df = sc2_tre_antibody_escape,
  label_col = "label"
)
sw5$bar_plot
gc()

hw5 <- compare_ranks_new(
  df = hiv_busy5,
  columns = c('EpistaticAntigenicity(structure)',
              'EpistaticAntigenicity(sequence)', 
              'EpistaticSurfaceAccessibility(Sequence)', 
              'DMS score', 
              'EVEscape score'),
  reference_col = "EpistaticAntigenicity(sequence)",
  filter_df = hiv_all_exp_mut,
  label_col = "label"
)

hw5$bar_plot

gc()


cor_matrix <- create_correlation_plot(
  df = sc2_busy4,
  virus_protein_name = 'sc2_spike',
  columns = c("Epistatic_Antigenicity_sequence",
              "EVEscape score", 
              "DMS score",
              "Epistatic_SurfaceAccessibility_sequence",
              "Epistatic_Antigenicity_structure"),
  filter_df =  sc2_tre_antibody_escape,
  label_col = "label"
)

