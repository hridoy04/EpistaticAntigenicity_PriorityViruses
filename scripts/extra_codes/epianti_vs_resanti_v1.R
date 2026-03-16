colnames(sc2_baal4)

sc2_baal4 %>% 
  group_by(position.x) %>%
  slice_max(Epistatic_Antigenicity_sequence, with_ties = FALSE) %>% 
  select(position.x, Antigenicity_sequence, Antigenicity_structure, 
                     Epistatic_Antigenicity_sequence, Epistatic_Antigenicity_structure, `DMS score`) %>% 
  ggplot(aes(Antigenicity_sequence, `DMS score`)) +
  geom_point() +
  geom_text(size =4.5,nudge_y =0.125,
            aes(label = ifelse(position.x %in% sc2_abescape_pos, position.x, ""))) 

#write the code to plot the relationship between Epistatic_Antigenicity_sequence and DMS score, with points colored by DMS score and labeled for positions in sc2_abescape_pos



sc2_baal4 %>% 
  group_by(position.x) %>%
  slice_max(Epistatic_Antigenicity_sequence, with_ties = FALSE) %>% 
  select(position.x, Antigenicity_sequence, Antigenicity_structure, 
                     Epistatic_Antigenicity_sequence, Epistatic_Antigenicity_structure, `DMS score`) %>% 
  ggplot(aes(Epistatic_Antigenicity_sequence, `DMS score`)) +
  geom_point() +
  geom_text(size =4.5,nudge_y =0.125,
            aes(label = ifelse(position.x %in% sc2_abescape_pos, position.x, ""))) 


sc2_baal4 %>% 
  group_by(position.x) %>%
  slice_max(Epistatic_Antigenicity_structure, with_ties = FALSE) %>% 
  drop_na() %>% 
  select(position.x, Antigenicity_sequence, Antigenicity_structure, 
         Epistatic_Antigenicity_sequence, Epistatic_Antigenicity_structure, `DMS score`) %>% 
  # mutate(`DMS score_rescaled` = scales::rescale(`DMS score`, to = c(0, 1)), 
  #        Epistatic_Antigenicity_structure_rescaled = scales::rescale(Epistatic_Antigenicity_structure, to = c(0, 1))) 
  ggplot(aes(Epistatic_Antigenicity_structure, `DMS score`)) +
  geom_point(aes(color = `DMS score`, size = `DMS score`), alpha = 0.5) +
  # scale_color_distiller(palette = 'Purples', 
  #                       limits = c(0,1),
  #                       # limits = c(min(.$`DMS score`), median(.$`DMS score`)), 
  #                       # limits = function(x) c(min(x$`DMS score`, na.rm = TRUE), median(x$`DMS score`, na.rm = TRUE)), 
  #                       direction = 1, 
  #                       oob = scales::squish) +
  geom_text(size =4.5,nudge_y =0.125,
            aes(label = ifelse(position.x %in% sc2_abescape_pos, position.x, ""))) +
  labs(title = 'SARS-CoV-2 spike protein', x = 'Epistatic Antigenicity', y = 'experimental DMS score') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 18)
  )


library(dplyr)
library(ggplot2)
library(viridis)

sc2_baal4 %>%
  group_by(position.x) %>%
  slice_max(Epistatic_Antigenicity_structure, with_ties = FALSE) %>%
  drop_na() %>%
  select(
    position.x,
    Antigenicity_sequence,
    Antigenicity_structure,
    Epistatic_Antigenicity_sequence,
    Epistatic_Antigenicity_structure,
    `DMS score`
  ) %>%
  ggplot(aes(Epistatic_Antigenicity_structure, `DMS score`)) +
  geom_point(
    aes(color = `DMS score`, size = `DMS score`),
    alpha = 0.8
  ) +
  scale_color_gradient(
    low = "white",
    high = "purple",
    limits = c(0, 6),
    oob = scales::squish
  ) +
  scale_size(range = c(2, 8), name = "DMS score") +
  geom_text(
    size = 5,
    nudge_y = 0.125,
    aes(label = ifelse(position.x %in% sc2_abescape_pos, position.x, ""))
  ) +
  labs(
    title = "SARS-CoV-2 spike protein",
    x = "Epistatic Antigenicity",
    y = "Experimental DMS score"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 24),
    plot.title = element_text(size = 22),
    
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.key.height = unit(1.2, "cm"),
    legend.key.width = unit(0.8, "cm")
  )


library(dplyr)
library(ggplot2)
library(ggrepel)

sc2_baal4 %>%
  group_by(position.x) %>%
  slice_max(Epistatic_Antigenicity_structure, with_ties = FALSE) %>%
  drop_na() %>%
  select(
    position.x,
    Antigenicity_sequence,
    Antigenicity_structure,
    Epistatic_Antigenicity_sequence,
    Epistatic_Antigenicity_structure,
    `DMS score`
  ) %>%
  ggplot(aes(Epistatic_Antigenicity_structure, `DMS score`)) +
  
  stat_density_2d(
    color = "grey60",
    linewidth = 0.5,
    bins = 6
  ) +
  
  geom_point(
    aes(color = `DMS score`, size = `DMS score`),
    alpha = 0.9
  ) +
  
  scale_color_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = median(sc2_baal4$`DMS score`, na.rm = TRUE),
    name = "DMS score",
    guide = guide_colorbar(
      barheight = unit(6, "cm"),
      barwidth = unit(1, "cm")
    )
  ) +
  
  scale_size(range = c(2, 8), name = "DMS score") +
  
  geom_text_repel(
    aes(label = ifelse(position.x %in% sc2_abescape_pos, position.x, "")),
    size = 5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  
  labs(
    title = "SARS-CoV-2 spike protein",
    x = "Epistatic Antigenicity",
    y = "Experimental DMS score"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 24),
    plot.title = element_text(size = 22),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.key.height = unit(1.2, "cm")
  )


library(dplyr)
library(ggplot2)
library(ggrepel)

median_x <- median(sc2_baal4$Antigenicity_structure, na.rm = TRUE)

sc2_baal4 %>%
  group_by(position.x) %>%
  slice_max(Epistatic_Antigenicity_sequence, with_ties = FALSE) %>%
  drop_na() %>%
  select(
    position.x,
    Antigenicity_sequence,
    Antigenicity_structure,
    Epistatic_Antigenicity_sequence,
    Epistatic_Antigenicity_structure,
    `DMS score`
  ) %>%
  ggplot(aes(Epistatic_Antigenicity_sequence, `DMS score`)) +
  
  # stat_density_2d(color = "grey60", linewidth = 0.5, bins = 6) +
  
  geom_vline(
    xintercept = median_x,
    linetype = "dashed",
    linewidth = 1
  ) +
  # geom_rect(aes(xmin = -Inf, xmax = median_x, ymin = -Inf, ymax = Inf),
  #           fill = "grey99", alpha = 0.2, inherit.aes = FALSE) +
  # 
  geom_point(
    aes(color = `DMS score`, size = `DMS score`),
    alpha = 0.9
  ) +
  
  scale_color_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    limits = c(-0.5,2),
    oob = scales::squish,
    midpoint = median(sc2_baal4$`DMS score`, na.rm = TRUE),
    name = "DMS score",
    guide = guide_colorbar(barheight = unit(6, "cm"))
  ) +
  
  scale_size(range = c(2, 8), name = "DMS score") +
  
  ggrepel::geom_text_repel(
    aes(label = ifelse(position.x %in% sc2_abescape_pos, position.x, "")),
    size = 5,
    max.overlaps = Inf
  ) +
  
  annotate(
    "text",
    x = median_x,
    y = Inf,
    label = paste0("Median = ", round(median_x, 3)),
    vjust = 1.5,
    size = 6
  ) +
  
  
  labs(
    title = "SARS-CoV-2 spike protein",
    x = "Epistatic Antigenicity",
    y = "Experimental DMS score"
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 24),
    plot.title = element_text(size = 22),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16)
  )


sc2_baal4 %>%
  group_by(position.x) %>%
  slice_max(Epistatic_Antigenicity_structure, with_ties = FALSE) %>%
  drop_na() %>%
  select(
    position.x,
    Antigenicity_sequence,
    Antigenicity_structure,
    Epistatic_Antigenicity_sequence,
    Epistatic_Antigenicity_structure,
    `DMS score`
  ) %>%
  ggplot(aes(Antigenicity_structure, `DMS score`)) +
  
  # stat_density_2d(color = "grey60", linewidth = 0.5, bins = 6) +
  
  geom_vline(
    xintercept = median_x,
    linetype = "dashed",
    linewidth = 1
  ) +
  # geom_rect(aes(xmin = -Inf, xmax = median_x, ymin = -Inf, ymax = Inf),
  #           fill = "grey99", alpha = 0.2, inherit.aes = FALSE) +
  # 
  geom_point(
    aes(color = `DMS score`, size = `DMS score`),
    alpha = 0.9
  ) +
  
  scale_color_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    limits = c(-0.5,2),
    oob = scales::squish,
    midpoint = median(sc2_baal4$`DMS score`, na.rm = TRUE),
    name = "DMS score",
    guide = guide_colorbar(barheight = unit(6, "cm"))
  ) +
  
  scale_size(range = c(2, 8), name = "DMS score") +
  
  ggrepel::geom_text_repel(
    aes(label = ifelse(position.x %in% sc2_abescape_pos, position.x, "")),
    size = 5,
    max.overlaps = Inf
  ) +
  
  annotate(
    "text",
    x = median_x,
    y = Inf,
    label = paste0("Median = ", round(median_x, 3)),
    vjust = 1.5,
    size = 6
  ) +
  
  
  labs(
    title = "SARS-CoV-2 spike protein",
    x = "Epistatic Antigenicity",
    y = "Experimental DMS score"
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 24),
    plot.title = element_text(size = 22),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16)
  )


library(ggrepel)


sc2_epi_all %>% 
  #left_join(sc2_perp, by = 'position') %>% 
  #filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
  ggplot(aes(accessibility, Ab_escape)) +
  #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
  #scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = Ab_escape), alpha = 0.9) +
  scale_color_distiller(palette = 'Spectral') +
  geom_text(size =4.5,nudge_y =0.125,
            aes(label = ifelse(position %in% sc2_abescape_pos, position, ""))) +
  labs(title = 'SARS-CoV-2 RBD', x = 'antibody accessiblity score (Epistatis based)') + 
  geom_vline(#data = ~ summarise (.x, median_accessibility = median(accessibility, na.rm = TRUE)),
             data = ~ .x %>% summarise(median_accessibility = median(accessibility, na.rm = TRUE)),
             aes(xintercept = median_accessibility),
             linetype = "dashed", color = "black", linewidth = 1) +
  geom_text(data= ~ summarise(.x, x = median(accessibility, na.rm = TRUE)),
           aes(label = 'Median', x = x, y = Inf), vjust = 1.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 18)
        )


sc2_epi_all %>% 
  ggplot(aes(accessibility, Ab_escape)) +
  
  geom_point(aes(color = Ab_escape, size = Ab_escape), alpha = 0.9) +
  
  scale_color_distiller(
    palette = "Spectral",
    name = "Experimental\nDMS score",
    guide = guide_colorbar(
      barheight = unit(6, "cm"),
      barwidth = unit(1, "cm")
    )
  ) +
  scale_size_continuous(
    name = "Experimental DMS score",
    range = c(1, 4)
  ) +
  geom_label_repel(
    aes(label = ifelse(position %in% sc2_abescape_pos, position, "")),
    size = 5,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.color = "grey40",
    max.overlaps = Inf
  ) +
  
  labs(
    #title = "SARS-CoV-2 spike",
    x = "Residual Antigenicity",
    y = "Experimental DMS score"
  ) +
  
  geom_vline(
    data = ~ summarise(.x, median_accessibility = median(accessibility, na.rm = TRUE)),
    aes(xintercept = median_accessibility),
    linetype = "dashed",
    color = "black",
    linewidth = 1
  ) +
  
  geom_text(
    data = ~ summarise(.x, x = median(accessibility, na.rm = TRUE)),
    aes(label = paste("Median: ", round(x, 3)), x = x), y = 0.40,
    vjust = 1.5,
    hjust = 0.5, 
    size = 6
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 24),
    plot.title = element_text(size = 20),
    
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.key.height = unit(1.2, "cm")
  )


plot_scatter_median <- function(data,
                                x_var,
                                y_var,
                                color_var,
                                size_var,
                                median_var,
                                median_label_y,
                                x_label,
                                label_positions = NULL) {

  
  data %>%
    ggplot(aes({{x_var}}, {{y_var}})) +
    
    geom_point(
      aes(color = {{color_var}}, size = {{size_var}}),
      alpha = 0.9
    ) +
    
    scale_color_distiller(
      palette = "Spectral",
      name = "Experimental\nDMS score",
      guide = guide_colorbar(
        barheight = unit(6, "cm"),
        barwidth = unit(1, "cm")
      )
    ) +
    
    scale_size_continuous(
      name = "Experimental DMS score",
      range = c(1, 6)
    ) +
    
    geom_label_repel(
      aes(label = ifelse(position %in% label_positions, position, "")),
      size = 5,
      box.padding = 0.5,
      point.padding = 0.4,
      segment.color = "grey40",
      max.overlaps = Inf
    ) +
    
    labs(
      #x = deparse(substitute(x_var)),
      x = x_label,
      y = "Experimental DMS score"
    ) +
    
    geom_vline(
      data = ~ summarise(.x, med = median({{median_var}}, na.rm = TRUE)),
      aes(xintercept = med),
      linetype = "dashed",
      color = "black",
      linewidth = 1
    ) +
    
    geom_text(
      data = ~ summarise(.x, med = median({{median_var}}, na.rm = TRUE)),
      aes(
        x = med,
        label = paste0("Median: ", round(med, 3))
      ),
      y = median_label_y,
      vjust = 1.5,
      hjust = 0.5,
      size = 6
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
      axis.text.y = element_text(size = 18),
      axis.title = element_text(size = 24),
      plot.title = element_text(size = 20),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 16),
      legend.key.height = unit(1.2, "cm")
    )
}

plot_scatter_median(
  sc2_epi_all,
  max_bepro_acc,
  Ab_escape,
  Ab_escape,
  Ab_escape,
  max_bepro_acc,
  median_label_y = 0.40,
  x_label = 'Epistatic Antigenicity',
  label_positions = sc2_abescape_pos
)
  

plot_scatter_median(
  sc2_epi_all,
  accessibility,
  Ab_escape,
  Ab_escape,
  Ab_escape,
  accessibility,
  median_label_y = 0.40,
  x_label = 'Epistatic Antigenicity',
  label_positions = sc2_abescape_pos
)

plot_scatter_median(hiv_epi_all_final, 
                    accessibility,
                    Ab_escape,
                    Ab_escape,
                    Ab_escape,
                    accessibility, 
                    median_label_y = 0.5, 
                    x_label = 'Residual Antigenicity', 
                    label_positions = hiv_ab_escape_sites_losalamos)


plot_scatter_median(flu_new29, 
                    Epistatic_Antigenicity_structure,
                    `DMS score`,
                    `DMS score`,
                    `DMS score`,
                    Epistatic_Antigenicity_structure, 
                    median_label_y = 0.3, 
                    x_label = 'Epistatic Antigenicity', 
                    label_positions = as.numeric(names(which(table(str_extract(flu_wsn33_true_antibody_escape$label, '\\d+')) > 1))))

plot_scatter_median(flu_new29, 
                    Antigenicity_structure,
                    `DMS score`,
                    `DMS score`,
                    `DMS score`,
                    Antigenicity_structure, 
                    median_label_y = 0.3, 
                    x_label = 'Residual Antigenicity', 
                    label_positions = as.numeric(names(which(table(str_extract(flu_wsn33_true_antibody_escape$label, '\\d+')) > 1))))

plot_delta_rank <- function(data, position_var, delta_var, tick_step = 25) {
  
  max_abs <- max(abs(dplyr::pull(data, {{delta_var}})), na.rm = TRUE)
  axis_limit <- ceiling(max_abs / tick_step) * tick_step
  axis_breaks <- seq(-axis_limit, axis_limit, by = tick_step)
  
  data %>%
    mutate(pos_factor = factor({{position_var}})) %>%
    ggplot(aes(x = reorder(position, {{delta_var}}), y = {{delta_var}}, fill = {{delta_var}} > 0)) +
    
    geom_col() +
    
    scale_fill_manual(
      values = c("blue", "red"),
      labels = c("Negative Δ", "Positive Δ"),
      name = "Δ Rank Direction"
    ) +
    
    scale_y_continuous(breaks = axis_breaks) +
    coord_flip() +  
    labs(
      x = "Position",
      y = "Δ Rank (rank of epistatic antigenicity - rank of residual antigenicity)",
      title = "Δ Rank per Position"
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16)
    )
}

plot_delta_rank(sc2_man, position, delta_rank, tick_step = 50)
plot_delta_rank(hiv_imp, position, rank_delta, tick_step = 50)
plot_delta_rank(flu_new32, position, delta_rank, tick_step = 50)

sc2_epi_all %>% pull(accessibility) %>% median(na.rm = TRUE)

binom.test(sum(delta_rank$delta_rank>0), length(delta_rank$delta_rank), p=0.5, alternative="greater")
wilcox.test(sc2_man$rank_epi,  sc2_man$rank_acc,  paired = TRUE, exact = F, alternative = "greater")


sc2_man %>%
  ggplot(aes(x = reorder(position, delta_rank), y = delta_rank, fill = delta_rank > 0)) +
  geom_col() +
  scale_fill_manual(
    values = c("blue", "red"),
    labels = c("Negative Δ", "Positive Δ"),
    name = "Δ Rank Direction"
  ) +
  scale_y_continuous(
    breaks = c(-100, -50, -25, 0, 25, 50, 75, 100, 200)  # set desired ticks
  ) +
  
  coord_flip() +   # flips X/Y so negative/positive bars appear naturally
  labs(
    x = "Position",
    y = "Δ Rank",
    title = "Divergent Δ Rank per Position"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  )


library(ggplot2)
library(dplyr)

sc2_man %>%
  mutate(position = factor(position)) %>%  # treat positions as discrete
  ggplot(aes(x = position, y = delta_rank, fill = delta_rank > 0)) +
  
  geom_col() +
  
  scale_fill_manual(
    values = c("blue", "red"),
    labels = c("Negative Δ", "Positive Δ"),
    name = "Δ Rank Direction"
  ) +
  
  scale_y_continuous(
    breaks = c(-100, -50, -25, 0, 25, 50, 100)  # set desired ticks
  ) +
  
  labs(
    x = "Position",
    y = "Δ Rank (rank_epi - rank_acc)",
    title = "Divergent Δ Rank per Position"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  )

sc2_man %>%
  ggplot(aes(x = delta_rank)) +
  geom_bar()

sc2_man %>%
  mutate(position = factor(position)) %>% 
  ggplot(aes( y = delta_rank)) +
  geom_segment(aes(x = position, xend = position, y = 0, yend = delta_rank),
               color = "grey") +
  geom_point(aes(color = delta_rank > 0), size = 3) +
  scale_color_manual(values = c("blue", "red"), labels = c("Negative Δ", "Positive Δ")) +
  labs(
    x = "Position",
    y = "Δ Rank",
    color = "Δ Rank Direction",
    title = "Delta rank per position"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )

#hiv rank test
hiv_df$rank_acc <- rank(hiv_df$disco_acc)
hiv_imp <- hiv_df %>% filter(position %in% 
                               hiv_ab_escape_sites_losalamos) %>% 
  mutate(rank_delta = rank_epi - rank_acc) %>% print(n=23)
wilcox.test(hiv_imp$rank_epi,  hiv_imp$rank_acc, paired = TRUE, exact = F, alternative = "greater")
sum(hiv_imp$rank_delta>0)

cor.test(hiv_imp$rank_epi, hiv_imp$rank_acc, method = 'kendall')


#flu rank
flu_new <- flu_wsn33_final %>% select(position, label, single_aa_antigenicity, very_significant_accessibility, evescape, max_escape_experiment) %>% 
  mutate(across(c(very_significant_accessibility, single_aa_antigenicity, evescape, max_escape_experiment), ~ rank(.x, ties.method = "min"), .names = "{.col}_rank")) %>% 
  filter(label %in% flu_wsn33_true_antibody_escape$label) %>% 
  select(ends_with('rank'), position) %>% 
  group_by(position) %>% 
  mutate(very_significant_accessibility_rank = max(very_significant_accessibility_rank)) %>% 
  select(single_aa_antigenicity_rank, very_significant_accessibility_rank, position) %>% 
  ungroup() %>% 
  distinct(position, .keep_all = T)
flu_new$rank_delta <- flu_new$very_significant_accessibility_rank - flu_new$single_aa_antigenicity_rank
wilcox.test(flu_new$very_significant_accessibility_rank, flu_new$single_aa_antigenicity_rank, paired = TRUE, exact = F, alternative = "greater")