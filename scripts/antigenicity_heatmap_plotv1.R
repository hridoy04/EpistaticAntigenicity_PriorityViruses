#this code will generate heatmap for epistatic antigenicity of each virus from the csv file derived from the python code that uses ESM-2

# ------------------- Load necessary libraries ------------------ #

library(tidyverse)
library(scales) # Required for the rescale() function
library(tidyverse)
library(RColorBrewer)
library(cowplot)  # only needed if extracting legend separately

#------------------- Define helper functions ------------------ #

# Amino acid classification function
classify_aa_simple <- function(aa) {
  case_when(
    aa %in% c('A', 'V', 'I', 'L', 'M', 'F', 'W', 'Y', 'P') ~ 'Hydrophobic',
    aa %in% c('S', 'T', 'C', 'N', 'Q', 'G') ~ 'Polar',
    aa %in% c('K', 'R', 'H') ~ 'Positive',
    aa %in% c('D', 'E') ~ 'Negative',
    TRUE ~ 'Unknown'
  )
}

#--------------------- Data Processing Function --------------------- #


process_epiant_data <- function(input_file, #derived from the python code that uses ESM-2; may contain many columns!! 
                                output_file, #the file that will be used for heatmap generation; should have columns: position, substituted_aa, chemical_group, epistatic_antigenicity
                                antigenicity_col = "significant_antigenicity_median") {
  
  # 1. Read input CSV
  data <- read_csv(input_file, show_col_types = FALSE)
  
  # 2. Process data
  epiant_processed <- data %>%
    distinct(label, .keep_all = TRUE) %>%
    select(position, ref, alt, label, !!sym(antigenicity_col)) %>%
    rename(
      pos = position,
      antigenicity = !!sym(antigenicity_col)
    ) %>%
    select(-c(ref, label)) %>%
    mutate(
      # Note: Ensure classify_aa_simple() is defined in your environment
      chem_alt = classify_aa_simple(alt), 
      antigenicity = rescale(antigenicity, to = c(0, 1))
    )
  
  # 3. Export to new CSV
  write_csv(epiant_processed, output_file)
  
  message(sprintf("Success! Processed data saved to: %s", output_file))
  
  # Return the dataframe invisibly in case you want to assign it to a variable
  invisible(epiant_processed)
}


# ------------------ Heatmap Generation Function ------------------ #



generate_antigenicity_heatmap <- function(file_path,
                                virus_name      = NULL,
                                output_file     = NULL,
                                position_col    = NULL,
                                substituted_col = NULL,
                                chemical_col    = NULL,
                                antigenicity_col = NULL,
                                show_legend     = FALSE,
                                height          = 3,
                                width           = 8.27,
                                dpi             = 300) {
  # ── 1. Read data ────────────────────────────────────────────────────────────
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # ── 2. Auto-detect columns (skipped if user supplies names) ─────────────────
  if (is.null(antigenicity_col)) {
    num_cols <- colnames(df)[sapply(df, is.numeric)]
    antigenicity_col <- num_cols[
      sapply(df[num_cols], function(x) any(x %% 1 != 0, na.rm = TRUE))
    ][1]
    position_col <- setdiff(num_cols, antigenicity_col)[1]
  }
  
  if (is.null(substituted_col)) {
    substituted_col <- names(df)[
      sapply(df, function(col) is.character(col) && all(nchar(na.omit(col)) == 1))
    ][1]
  }
  
  if (is.null(chemical_col)) {
    chemical_col <- setdiff(
      colnames(df),
      c(antigenicity_col, position_col, substituted_col)
    )[1]
  }
  
  message(sprintf(
    "Columns used → position: %s | substituted_aa: %s | chemical_group: %s | antigenicity: %s",
    position_col, substituted_col, chemical_col, antigenicity_col
  ))
  
  # ── 3. Rename for plotting ───────────────────────────────────────────────────
  df <- df %>%
    rename(
      position       = !!position_col,
      substituted_aa = !!substituted_col,
      chemical_group = !!chemical_col,
      antigenicity   = !!antigenicity_col
    )
  
  # ── 4. Build plot ────────────────────────────────────────────────────────────
  x_breaks <- seq(
    floor(min(df$position,   na.rm = TRUE) / 50) * 50,
    ceiling(max(df$position, na.rm = TRUE) / 50) * 50,
    by = 50
  )
  
  p <- df %>%
    ggplot(aes(x = position, y = substituted_aa, fill = antigenicity)) +
    geom_tile(linewidth = 0.25) +
    scale_fill_distiller(palette = "RdBu", na.value = "grey90") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0), breaks = x_breaks) +
    facet_grid(chemical_group ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = if (!is.null(virus_name)) paste("Epistatic Antigenicity –", virus_name) else NULL,
      x = "", y = ""
    ) +
    theme_classic() +
    theme(
      axis.text.x    = element_text(size = 10),
      axis.text.y    = element_text(size = 8),
      strip.text.y   = element_text(size = 6, face = "bold"),
      panel.spacing  = unit(0.2, "mm"),
      axis.ticks     = element_line(linewidth = 0.4),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title   = element_text(size = 12, face = "bold"),
      legend.text    = element_text(size = 10),
      legend.key.width = unit(1.5, "cm")
    )
  
  # ── 5. Save ──────────────────────────────────────────────────────────────────
  if (is.null(output_file)) {
    prefix <- if (!is.null(virus_name)) tolower(gsub(" ", "_", virus_name)) else "virus"
    output_file <- paste0(prefix, "-epiant-heatmap.png")
  }
  
  ggsave(filename = output_file, plot = p,
         height = height, width = width, units = "in", dpi = dpi)
  message("Saved → ", output_file)
  
  invisible(p)
}


#user input:

primary_epistatic_antigenicity_file <- 'CCHFV_GPC_epistatic_accessibility_antigenicity.csv'
output_antigenicity_file <- 'cchfv_epiant_processed.csv'
virus_glycoprotein_name <- 'CCHFV GPC'



#execute the functions:

process_epiant_data(primary_epistatic_antigenicity_file, output_antigenicity_file, antigenicity_col = 'significant_antigenicity_median')

generate_antigenicity_heatmap(
  file_path = output_antigenicity_file,
  virus_name = 'CCHFV GPC',
  show_legend = F,
  height          = 3,
  width           = 8.27,
  dpi             = 300
)
