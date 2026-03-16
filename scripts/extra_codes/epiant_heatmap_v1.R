epiant <- read_csv('sc2_spike_epistatic_accessibility_final_v2.csv')
epiant_sc2 <- epiant %>% distinct(label, .keep_all = T) %>% 
  select(position.x, ref.x, alt.x, label, Epistatic_Antigenicity_sequence) %>% 
  rename(ref = ref.x, pos = position.x, alt = alt.x)

classify_aa_simple <- function(aa) {
  case_when(
    aa %in% c('A', 'V', 'I', 'L', 'M', 'F', 'W', 'Y', 'P') ~ 'Hydrophobic',
    aa %in% c('S', 'T', 'C', 'N', 'Q', 'G') ~ 'Polar',
    aa %in% c('K', 'R', 'H') ~ 'Positive',
    aa %in% c('D', 'E') ~ 'Negative',
    TRUE ~ 'Unknown'
  )
}

epiant_sc2 <- epiant_sc2 %>% mutate(chem_alt = classify_aa_simple(alt)) %>% 
  select(-c(ref,label))
epiant_sc2

#make a heatmap where the x axis shows pos, y axis alt and color shows epiant, the order of alt should be based on chem_alt

epiant_sc2 %>% filter(pos >= 450 & pos < 511) %>% ggplot() + geom_tile(aes(x=pos, y=alt, fill=rescale(Epistatic_Antigenicity_sequence, to = c(0,1)))) +
  scale_fill_viridis_c(name = "Epistatic\nAntigenicity", option = "viridis") +
  scale_x_continuous(breaks = seq(450, 510, by = 5)) +
  facet_grid(chem_alt ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Epistatic Antigenicity for different amino acids at each position",
    x = "Position",
    y = "Amino Acid grouped by Chemical properties"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 11),
    strip.text.y = element_text(size = 11, face = "bold")
  )



epiant_sc2 %>%
  ggplot(aes(x = pos, y = alt, fill = Epistatic_Antigenicity_sequence)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_viridis_c(name = "Epistatic\nAntigenicity", option = "turbo") +
  facet_grid(chem_alt ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Epistatic Antigenicity by Position and Amino Acid Chemistry",
    x = "Position",
    y = "Amino Acid"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 11),
    strip.text.y = element_text(size = 11, face = "bold")
  )

print(heatmap_faceted)

library(ComplexHeatmap)


heatmap(epiant_sc2)

epiant_sc2 %>% filter(pos == 500)

#hiv


hiv_epiant <- read_csv('HIVbg505_Env_epistatic_accessibility_antigenicity_v2.csv')
epiant_hiv <- hiv_epiant %>% distinct(label, .keep_all = T) %>% 
  select(position, ref, alt, label, significant_antigenicity_median) %>% 
  rename(pos = position) %>% 
  select(-c(ref,label)) %>% 
  mutate(chem_alt = classify_aa_simple(alt), 
         significant_antigenicity_median = rescale(significant_antigenicity_median, to = c(0,1)))

epiant_hiv %>% ggplot() + geom_tile(aes(x=pos, y=alt, fill=significant_antigenicity_median)) +
  scale_fill_viridis_c(name = "Epistatic\nAntigenicity", option = "viridis") +
  scale_x_continuous(breaks = seq(0, 865, by = 50)) +
  facet_grid(chem_alt ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Epistatic Antigenicity for different amino acids at each position in HIV-1 BG505 Env",
    x = "Position",
    y = "Amino Acid grouped by Chemical properties"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 11),
    strip.text.y = element_text(size = 11, face = "bold")
  )

flu_epiant <- read_csv('flu_wsn33_epistatic_accessibility_final_v2.csv')
colnames(flu_epiant)
epiant_flu <- flu_epiant %>% distinct(label, .keep_all = T) %>% 
  select(position.x, ref.x, alt.x, label, Epistatic_Antigenicity_sequence) %>% 
  rename(pos = position.x, alt = alt.x) %>% 
  select(-c(ref.x,label)) %>% 
  mutate(chem_alt = classify_aa_simple(alt), 
         significant_antigenicity_median = rescale(Epistatic_Antigenicity_sequence, to = c(0,1)))

epiant_flu %>% ggplot() + geom_tile(aes(x=pos, y=alt, fill=significant_antigenicity_median)) +
  scale_fill_viridis_c(name = "Epistatic\nAntigenicity", option = "viridis") +
  scale_x_continuous(breaks = seq(0, 565, by = 25)) +
  facet_grid(chem_alt ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Epistatic Antigenicity for different amino acids at each position in Influenza A WSN33 HA",
    x = "Position",
    y = "Amino Acid grouped by Chemical properties"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 11),
    strip.text.y = element_text(size = 11, face = "bold")
  )


#borna

borna_epiant <- read_csv('Bornavirus_EnvelopeGlycoproteinp57_epistatic_accessibility_antigenicity.csv')

epiant_borna <- borna_epiant %>% distinct(label, .keep_all = T) %>% 
  mutate(significant_antigenicity_median = replace_na(significant_antigenicity_median, 0)) %>%
  select(position, ref, alt, label, significant_antigenicity_median) %>% 
  rename(pos = position) %>% 
  select(-c(ref,label)) %>% 
  mutate(chem_alt = classify_aa_simple(alt), 
         significant_antigenicity_median = rescale(significant_antigenicity_median, to = c(0,1)))

epiant_borna %>% ggplot() + geom_tile(aes(x=pos, y=alt, fill=significant_antigenicity_median)) +
  scale_fill_viridis_c(name = "Epistatic\nAntigenicity", option = "viridis") +
  scale_x_continuous(breaks = seq(1, 500, by = 25)) +
  facet_grid(chem_alt ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Epistatic Antigenicity for different amino acids at each position in Bornavirus-1 p57",
    x = "Position",
    y = "Amino Acid grouped by Chemical properties"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 11),
    strip.text.y = element_text(size = 11, face = "bold")
  )

#marburg virus

marburg_gp12 <- read_csv('Marburg_GP1GP2_epistatic_accessibility_antigenicity.csv')
epiant_marburg <- marburg_gp12 %>% distinct(label, .keep_all = T) %>% 
  mutate(significant_antigenicity_median = replace_na(significant_antigenicity_median, 0)) %>% 
  select(position, ref, alt, label, significant_antigenicity_median) %>% 
  rename(pos = position) %>% 
  select(-c(ref,label)) %>% 
  mutate(chem_alt = classify_aa_simple(alt), 
         significant_antigenicity_median = rescale(significant_antigenicity_median, to = c(0,1)))
epiant_marburg %>% ggplot() + geom_tile(aes(x=pos, y=alt, fill=significant_antigenicity_median)) +
  scale_fill_viridis_c(name = "Epistatic\nAntigenicity", option = "viridis") +
  scale_x_continuous(breaks = seq(1, 700, by = 50)) +
  facet_grid(chem_alt ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Epistatic Antigenicity for different amino acids at each position in Marburg GP1/GP2",
    x = "Position",
    y = "Amino Acid grouped by Chemical properties"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 11),
    strip.text.y = element_text(size = 11, face = "bold")
  )



#cchfv
cchfv_epiant <- read_csv('CCHFV_GPC_epistatic_accessibility_antigenicity.csv')
epiant_cchfv <- cchfv_epiant %>% mutate(significant_antigenicity_median = replace_na(significant_antigenicity_median, 0)) %>% 
  distinct(label, .keep_all = T) %>% 
  select(position, ref, alt, label, significant_antigenicity_median) %>% 
  rename(pos = position) %>% 
  select(-c(ref,label)) %>% 
  mutate(pos = pos +993, 
         chem_alt = classify_aa_simple(alt), 
         significant_antigenicity_median = rescale(significant_antigenicity_median, to = c(0,1))) 

epiant_cchfv %>% ggplot() + geom_tile(aes(x=pos, y=alt, fill=significant_antigenicity_median)) +
  scale_fill_viridis_c(name = "Epistatic\nAntigenicity", option = "viridis") +
  scale_x_continuous(breaks = seq(994, 1550, by = 50)) +
  facet_grid(chem_alt ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Epistatic Antigenicity for different amino acids at each position in CCHFV GPC",
    x = "Position",
    y = "Amino Acid grouped by Chemical properties"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 11),
    strip.text.y = element_text(size = 11, face = "bold")
  )

epiant_cchfv %>% ggplot() + 
  geom_tile(aes(x=pos, y=alt, fill=significant_antigenicity_median),
            color="white",
            #alpha = 0.2, 
            size=0.25) +
  scale_fill_distiller(palette = 'RdBu') +
  xlim(1300, 1400) +
  #scale_x_continuous(breaks = seq(1300, 1400, by = 25)) +
  #facet_grid(chem_alt ~ ., scales = "free_y", space = "free_y") +
  labs(
    #title = "Epistatic Antigenicity for different amino acids at each position in CCHFV GPC",
    x = "Position",
    y = "Amino acids grouped by chemical properties"
  ) +
  #coord_fixed(0.8) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    #plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18),
    strip.text.y = element_text(size = 18, face = "bold"), 
    panel.spacing = unit(0.1, "mm")
  )


epiant_cchfv %>% ggplot() + 
  geom_tile(aes(x=pos, y=alt, fill=significant_antigenicity_median),
            #color="white",
            #alpha = 0.2, 
            size=0.25) +
  scale_fill_distiller(palette = 'RdBu') +
  #xlim(1300, 1400) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(1000, 1675, by = 50)) +
  facet_grid(chem_alt ~ ., scales = "free_y", space = "free_y") +
  labs(
    #title = "Epistatic Antigenicity for different amino acids at each position in CCHFV GPC",
    x = "",
    y = ""
  ) +
  #coord_fixed() +
  theme_classic() +
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
    legend.position =  'none'
  )
ggsave(ec, filename="cchfv-mod1.png", height=8.5, width=18.8, units="in", dpi=300)

#function to do for any virus
library(tidyverse)
library(scales)

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

# Main function to process and plot epistatic antigenicity
plot_epiant_heatmap <- function(file_name, 
                                antigenicity_col = "significant_antigenicity_median",
                                title_suffix = "",
                                x_breaks = 50) {
  
  # Read CSV file
  data <- read_csv(file_name)
  
  # Process data
  epiant_processed <- data %>%
    distinct(label, .keep_all = TRUE) %>%
    select(position, ref, alt, label, !!sym(antigenicity_col)) %>%
    rename(
      pos = position,
      antigenicity = !!sym(antigenicity_col)
    ) %>%
    select(-c(ref, label)) %>%
    mutate(
      chem_alt = classify_aa_simple(alt),
      antigenicity = rescale(antigenicity, to = c(0, 1))
    )
  
  # Extract organism/protein name from file for title
  protein_name <- file_name %>%
    str_remove("\\.csv$") %>%
    str_to_upper()
  
  # Create heatmap
  heatmap <- epiant_processed %>%
    ggplot() +
    geom_tile(aes(x = pos, y = alt, fill = antigenicity)) +
    scale_fill_viridis_c(name = "Epistatic\nAntigenicity", option = "viridis") +
    scale_x_continuous(breaks = seq(0, max(epiant_processed$pos), by = x_breaks)) +
    facet_grid(chem_alt ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = paste("Epistatic Antigenicity by Amino Acid Position in", protein_name, title_suffix),
      x = "Position",
      y = "Amino Acid (grouped by Chemical properties)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 11),
      strip.text.y = element_text(size = 11, face = "bold"),
      panel.spacing = unit(1, "lines")
    )
  
  # Return both plot and processed data
  return(list(
    plot = heatmap,
    data = epiant_processed
  ))
}



# Main function to process and plot epistatic antigenicity
plot_epiant_heatmap_new <- function(file_name, 
                                antigenicity_col = "significant_antigenicity_median",
                                title_suffix = "",
                                x_breaks = 50) {
  
  # Read CSV file
  data <- read_csv(file_name)
  
  # Process data
  epiant_processed <- data %>%
    distinct(label, .keep_all = TRUE) %>%
    select(position, ref, alt, label, !!sym(antigenicity_col)) %>%
    rename(
      pos = position,
      antigenicity = !!sym(antigenicity_col)
    ) %>%
    select(-c(ref, label)) %>%
    mutate(
      chem_alt = classify_aa_simple(alt),
      antigenicity = rescale(antigenicity, to = c(0, 1))
    )
  
  # Extract organism/protein name from file for title
  protein_name <- file_name %>%
    str_remove("\\.csv$") %>%
    str_to_upper()
  
  # Create heatmap
  heatmap <- epiant_processed %>%
    ggplot() +
    geom_tile(aes(x = pos, y = alt, fill = antigenicity)) +
    scale_fill_distiller(name = "Epistatic\nAntigenicity", palette = "RdBu") +
    scale_x_continuous(breaks = seq(0, max(epiant_processed$pos), by = x_breaks)) +
    facet_grid(chem_alt ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = paste("Epistatic Antigenicity by Amino Acid Position in", protein_name, title_suffix),
      x = "Position",
      y = "Amino Acid (grouped by Chemical properties)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 11),
      strip.text.y = element_text(size = 11, face = "bold"),
      panel.spacing = unit(1, "lines")
    )
  
  # Return both plot and processed data
  return(list(
    plot = heatmap,
    data = epiant_processed
  ))
}


# Usage examples:

#cchfv

cchfv_result <- plot_epiant_heatmap('CCHFV_GPC_epistatic_accessibility_antigenicity.csv')
cchfv_result$plot
cchfv_result$data

# Example 2: Basic usage
result_hiv <- plot_epiant_heatmap_new('HIVbg505_Env_epistatic_accessibility_antigenicity_v2.csv')
result_hiv$plot

# Example 2: Custom title suffix
result_hiv <- plot_epiant_heatmap(
  'HIVbg505_Env_epistatic_accessibility_antigenicity_v2.csv',
  title_suffix = "- BG505 Strain"
)
result_hiv$plot

# Example 3: Custom x-axis breaks
result_hiv <- plot_epiant_heatmap(
  'HIVbg505_Env_epistatic_accessibility_antigenicity_v2.csv',
  x_breaks = 100
)
result_hiv$plot

# Example 4: Access processed data
epiant_hiv <- result_hiv$data
head(epiant_hiv)

# Example 5: Multiple files
files <- c(
  'HIVbg505_Env_epistatic_accessibility_antigenicity_v2.csv',
  'other_protein.csv'
)

results <- map(files, ~ plot_epiant_heatmap(.x))

# View all plots
results[[1]]$plot
results[[2]]$plot

# Save plot
result <- plot_epiant_heatmap('HIVbg505_Env_epistatic_accessibility_antigenicity_v2.csv')
ggsave('epiant_heatmap.png', result$plot, width = 14, height = 10, dpi = 300)

#flu_h3n2_k

flu_h3n2_k <- read_csv('/Users/rubayetalam/InfluenzaAH3N2K_HA_epistatic_accessibility_antigenicity.csv')
flu_h3n2_k %>% 
  # filter(ref != alt) %>% 
  mutate(position_new = position -16) %>% 
  mutate(significant_antigenicity_median = ifelse(is.na(significant_antigenicity_median), 0.2, significant_antigenicity_median)) %>% 
  mutate(significant_antigenicity_median = rescale(significant_antigenicity_median, to = c(0,1))) %>%
  filter(position >= 120, position < 251) %>% 
  ggplot() + geom_tile(aes(x=position, y=alt, fill=significant_antigenicity_median)) +
  scale_fill_viridis_c(name = "Epistatic\nAntigenicity", option = "viridis") +
  scale_x_continuous(breaks = seq(120, 250, by = 10)) +
  facet_grid(classify_aa_simple(alt) ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Epistatic Antigenicity for different amino acids at each position in Influenza A H3N2 K HA",
    x = "Position",
    y = "Amino Acid grouped by Chemical properties"
  ) +
  theme_bw()
  

flu_h3n2_j <- read_csv('/Users/rubayetalam/InfluenzaH3N2J_HA_epistatic_accessibility_antigenicity.csv')
flu_h3n2_j %>% 
  mutate(position_new = position -16) %>% 
  mutate(significant_antigenicity_median = ifelse(is.na(significant_antigenicity_median), 0.2, significant_antigenicity_median)) %>% 
  mutate(significant_antigenicity_median = rescale(significant_antigenicity_median, to = c(0,1))) %>%
  filter(position_new >= 130, position_new < 191) %>% 
  ggplot() + geom_tile(aes(x=position_new, y=alt, fill=significant_antigenicity_median)) +
  scale_fill_viridis_c(name = "Epistatic\nAntigenicity", option = "viridis") +
  scale_x_continuous(breaks = seq(130, 190, by = 5)) +
  facet_grid(classify_aa_simple(alt) ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Epistatic Antigenicity for different amino acids at each position in Influenza A H3N2 K HA",
    x = "Position",
    y = "Amino Acid grouped by Chemical properties"
  ) +
  theme_bw()


flu_h3n2_j %>% 
  mutate(position_new = position -16) %>% 
  mutate(significant_antigenicity_median = ifelse(is.na(significant_antigenicity_median), 0.2, significant_antigenicity_median)) %>% 
  mutate(significant_antigenicity_median = rescale(significant_antigenicity_median, to = c(0,1))) %>%
  filter(position_new >= 130, position_new < 191) %>%
  mutate(epiant_group = case_when(
    significant_antigenicity_median >= 0.74 ~ 'top 1%',
    significant_antigenicity_median >= 0.69 ~ 'top 5%',
    significant_antigenicity_median >= 0.60 ~ 'top 10%',
    significant_antigenicity_median >= 0.57 ~ 'top 25%',
    significant_antigenicity_median < 0.57 ~ 'below 25%'
  )) %>% 
  ggplot() + geom_tile(aes(x=position_new, y=alt, fill=epiant_group)) +
  # scale_fill_viridis_c(name = "Epistatic\nAntigenicity", option = "viridis", limits = c(0, 0.60), oob = scales::squish) +
    
  scale_x_continuous(breaks = seq(130, 190, by = 5)) +
  facet_grid(classify_aa_simple(alt) ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Epistatic Antigenicity for different amino acids at each position in Influenza A H3N2 K HA",
    x = "Position",
    y = "Amino Acid grouped by Chemical properties"
  ) +
  theme_bw()

#claude
flu_h3n2_j %>% 
  mutate(position_new = position - 16) %>% 
  mutate(significant_antigenicity_median = ifelse(is.na(significant_antigenicity_median), 0.2, significant_antigenicity_median)) %>% 
  mutate(significant_antigenicity_median = rescale(significant_antigenicity_median, to = c(0,1))) %>%
  filter(position_new >= 130, position_new < 191) %>%
  mutate(epiant_group = case_when(
    significant_antigenicity_median >= 0.69 ~ 'top 1%',
    significant_antigenicity_median >= 0.55 ~ 'top 5%',
    significant_antigenicity_median >= 0.51 ~ 'top 10%',
    significant_antigenicity_median >= 0.45 ~ 'top 25%',
    significant_antigenicity_median < 0.45 ~ 'below 25%'
  )) %>% 
  mutate(epiant_group = factor(epiant_group, levels = c('top 1%', 'top 5%', 'top 10%', 'top 25%', 'below 25%'))) %>%  # Add this line to order the factor
  ggplot() + 
  geom_tile(aes(x = position_new, y = alt, fill = epiant_group)) +
  scale_fill_manual(
    name = "Epistatic\nAntigenicity",
    values = c('below 25%' = '#440154FF', 'top 10%' = '#238A8DFF', 'top 25%' = '#3F4788FF', 'top 1%' = '#fde724', 'top 5%' = '#94D840FF')
   
  ) +
  scale_x_continuous(breaks = seq(130, 190, by = 5)) +
  facet_grid(classify_aa_simple(alt) ~ ., scales = "free_y", space = "free_y") +
  labs(
    #title = "Epistatic Antigenicity for different amino acids at each position in Influenza A H3N2 K HA",
    x = "Residue number in HA protein (H3 numbering)",
    y = "Mutated Amino Acid (grouped by Chemical properties)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20, face = 'bold'),
    axis.title.x = element_text(margin=margin(t=90), size = 20, face = 'bold'),
    #plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 16, face = 'bold'),
    strip.text.y = element_text(size = 17, face = "bold")
  )


f3j <- flu_h3n2_j %>% 
  mutate(position_new = position -16) %>% 
  mutate(significant_antigenicity_median = ifelse(is.na(significant_antigenicity_median), 0.2, significant_antigenicity_median)) %>% 
  mutate(significant_antigenicity_median = rescale(significant_antigenicity_median, to = c(0,1))) %>%
  #filter(position_new >= 130, position_new < 191) %>% 
  select(position_new , label, significant_antigenicity_median)
  # quantile(prob=90/100, na.rm =T)
  
f3j %>% quantile(prob=99/100, na.rm =T)
f3j %>% filter(position_new %in% k_mut) %>% print(n=140)
summary(f3j$significant_antigenicity_median)

flu_h3n2_ref <- read_csv('/Users/rubayetalam/InfluenzaH3N2Ref_HA_epistatic_accessibility_antigenicity.csv')

flu_h3n2_ref %>%  
  mutate(position_new = position -16) %>% 
  mutate(significant_antigenicity_median = ifelse(is.na(significant_antigenicity_median), 0.2, significant_antigenicity_median)) %>% 
  mutate(significant_antigenicity_median = rescale(significant_antigenicity_median, to = c(0,1))) %>%
  filter(position_new >= 120, position_new < 210) %>% 
  ggplot() + geom_tile(aes(x=position_new, y=alt, fill=significant_antigenicity_median)) +
  scale_fill_viridis_c(name = "Epistatic\nAntigenicity", option = "viridis") +
  scale_x_continuous(breaks = seq(120, 210, by = 10)) +
  facet_grid(classify_aa_simple(alt) ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Epistatic Antigenicity for different amino acids at each position in Influenza A H3N2 Ref HA",
    x = "Position",
    y = "Amino Acid grouped by Chemical properties"
  ) +
  theme_bw()
  
flu_h3n2_result <- plot_epiant_heatmap('flu_h3n2_k_epistatic_accessibility_final_v2.csv')


flu_j_check <- flu_h3n2_j %>% 
  mutate(position_new = position -16) %>% 
  filter(position_new %in% k_mut) %>% 
  select(position_new, ref, alt, significant_antigenicity_median, very_significant_antigenicity_median) 
hist(flu_h3n2_j$significant_antigenicity_median)

#check the plm prob:

flu_h3n2_j_plm <- read_csv('InfluenzaH3N2J_HA_DMS_scores.csv')
flu_j_prob <- flu_h3n2_j_plm %>% 
  mutate(position_new = position -16) %>% 
  # filter(position_new %in% k_mut) %>%
  select(position_new, ref, alt, probability)
