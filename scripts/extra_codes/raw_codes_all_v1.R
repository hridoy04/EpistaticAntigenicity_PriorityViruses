
# Load the data
entropy_files <- list.files('~/prob_csv_final_viruses/', pattern = '*entropy.csv', full.names = TRUE)

entropy_virus_new <- lapply(entropy_files, function(x) {
  read_csv(x) %>%
    mutate(virus = strsplit(basename(x), '_')[[1]][1])
}) %>% bind_rows() %>% rename('position' = `...1`)

entropy_virus %>% ggplot(aes(position, entropy, color = virus)) + 
  geom_line() + 
  scale_x_continuous(breaks = seq(0, 1650, 100)) +
  #facet_wrap(~virus, ncol = 1, strip.position="right") + 
  facet_wrap(~factor(virus, 
                     levels = c('bornavirus', 'nairo', 'marburg', 'influenza',
                                'hiv', 'SARSCOV2')), 
             ncol = 1, strip.position="right") + 
  theme_bw() +
  theme (axis.text = element_text(size = 14), 
        axis.title = element_text(size = 32), 
        plot.title = element_text(size = 30),
        strip.text = element_text(size = 16),
    legend.position = 'none')


# or,
# combined_df_entropy <- entropy_files %>%
#   lapply(function(x) {
#     read_csv(x) %>%
#       mutate(virus = strsplit(basename(x), '_')[[1]][1])
#   }) %>%
#   bind_rows()

acc_files <- list.files('~/accessibility_corr', pattern = '.csv', full.names = TRUE)

acc_pred_virus <- lapply(acc_files, function(x) {
  read_csv(x) %>%
    mutate(virus = strsplit(basename(x), '_')[[1]][1])
}) %>% bind_rows() %>% rename('position' = `...1`)

ggplot(acc_pred_virus, aes(`Observed Accessibility`, `Predicted Accessibility`)) +
  geom_point(size =0.5, alpha=0.5) + geom_smooth(method = 'lm') +
  annotate('text', x=0.25, y = 0.75, label = 'r=0.83\nr2=0.7\np << 0.01',
           size = unit(10, "pt")) +
  labs(x = 'observed accessibility', y = 'predicted accessibility',
       title = 'correlation between observed vs \npredicted (semantic embeddings) accessibility') +
  theme_bw() +
  facet_wrap(~virus) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28),
        strip.text = element_text(size = 24)
  ) 

ggplot(emb_acc, aes(Y, Y_pred)) +
  geom_point() + geom_smooth(method = 'lm') +
  annotate('text', x=0.25, y = 0.75, label = 'r=0.83\nr2=0.7\np << 0.01',
           size = unit(10, "pt")) +
  labs(x = 'observed accessibility', y = 'predicted accessibility',
       title = 'correlation between observed vs \npredicted (semantic embeddings) accessibility') +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 28),
        plot.title = element_text(size = 28)
  )

#> bar
# index a b
# 1 s1 3
# 2 s2 5

df_rect <- data.frame(
  x = rep(c(2, 5, 7, 9, 12), 2),
  y = rep(c(1, 2), each = 5),
  z = factor(rep(1:5, each = 2)),
  w = rep(diff(c(0, 4, 6, 8, 10, 14)), 2)
)
df_rect

bar$a <- factor(bar$a, levels = c('s2', 's1'))
bar %>% ggplot(aes(x = 'pos', y = b, fill = a)) + 
  geom_(stat = 'identity', position = 'stack') +
  coord_flip()
bar %>% ggplot(aes('pos', b, fill = factor(a))) + 
  geom_raster()

borna <- entropy_virus %>% filter(virus == 'bornavirus')
borna <- borna %>% mutate(domain = case_when(position <= 200 ~ 's1', position > 200 ~ 's2')) %>% group_by(domain) %>% mutate(n=n())

b1 <- borna %>% ggplot(aes(position, entropy)) + 
  geom_line() + 
  theme_bw() +
  theme (axis.text = element_text(size = 14), 
         axis.title = element_text(size = 32), 
         plot.title = element_text(size = 30),
         strip.text = element_text(size = 16),
         legend.position = 'none') 

borna %>% distinct(n) %>% mutate(domain = factor(domain, levels = c('s2', 's1'))) %>% 
  ggplot(aes(y = 'pos', x = n, fill = domain)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  theme_void()

borna %>% ggplot() + 
  geom_line(aes(position, entropy)) + 
  ggplot(data = borna %>% distinct(n) %>% mutate(domain = factor(domain, levels = c('s2', 's1'))), 
         aes(y = 'pos', x = n, fill = domain)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  theme_void()

ggplot() + 
  # Bar plot (Plotted first so it's on top)
  geom_bar(
    data = borna %>% distinct(n) %>%
      mutate(domain = factor(domain, levels = c('s2', 's1'))),
    aes(x = n, y = 'pos', fill = domain),  # Shifted up
    stat = 'identity',
    position = 'stack'
  ) +
  
  # Line plot (Plotted after so it appears below)
  geom_line(data = borna, aes(position, entropy)) +
  
  theme_void()

data.frame(
  xmin = c(1, 202),   # Start positions
  xmax = c(201, max(borna$position)),  # End positions
  ymin = max(borna$entropy) * 1.05,  # Top of bars (shifted up)
  ymax = max(borna$entropy) * 1.2,   # Bottom of bars (shifted up)
  domain = factor(c('s1', 's2'), levels = c('s2', 's1'))  # Color groups
)


b2 / b1
plot(b1+b2)
library(cowplot)
library(patchwork)

borna %>% distinct(n)
  
# Combined plot with bars stacked above the line
ggplot() +
  # Bar plot (stacked, flipped)
  geom_bar(data = borna %>% distinct(n), 
           aes(x = 'pos', y = n, fill = domain), 
           stat = "identity", position = "stack", width = 0.5) +  # Adjust width
  coord_flip() +  # Flip coordinates for horizontal bars
  # Line plot (entropy vs. position)
  geom_line(data = borna, aes(x = position, y = entropy), size = 1) +
  # Styling
  #scale_y_continuous(sec.axis = dup_axis(name = "Bar Values")) +  # Creates a second y-axis
  theme_bw() +
  theme(
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 32), 
    plot.title = element_text(size = 30),
    strip.text = element_text(size = 16),
    legend.position = 'none',
    axis.title.y.right = element_text(size = 20)  # Adjust bar axis title
  ) +
  labs(x = "Position", y = "Entropy", title = "Combined Line & Bar Chart")


#Load the data for entropy correlation

# Load the data
ent_corr_files <- list.files('~', pattern = '^correlation_sw\\d{2,3}.csv$', full.names = TRUE)
ent_corr_files

ent_corr_sw_all <- lapply(ent_corr_files, function(x) {
  read_csv(x) %>%
    mutate(virus = strsplit(basename(x), '_')[[1]][1])
}) %>% bind_rows()

lapply(ent_corr_files, function(x) {
  assign(strsplit(basename(x), '\\.csv')[[1]][1], read_csv(x), envir = .GlobalEnv)}
  )


sw <- read_csv('correlation_sw10.csv')
ent_corr_sw <- sw %>% separate(position, ", ", into = c('first', 'last')) %>% 
  mutate(first = as.integer(str_replace(first, '\\(','')), 
         last = as.integer(str_replace(last, '\\)', '')), 
         position = (first + last)/2) 
ent_corr_sw 
ent_corr_sw %>% ggplot(aes(position, correlation_10)) + 
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1250, 50)) +
  geom_smooth(method = 'loess', 
              color = 'darkgoldenrod1', 
              span = 0.25) +
  geom_hline (yintercept = 0.5, 
              linetype = 'dashed',
              linewidth = 2,
              color = 'darkgreen') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 30), 
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 32), 
        legend.position = 'none')


# Load the data


vdfsm <- read_csv('virus_df.csv') %>% rename(Borna = borna, HIV= hiv, Influenza = influ, Marburg = marburg, Nairo = orthonairo, 'SARS-CoV-2' = sc2) %>% as.data.frame()
rownames(vdfsm) <- vdfsm$`...1`
vdfsm <- vdfsm %>% select(-`...1`)

corrplot(as.matrix(vdfsm), col.lim = c(0,1), order = 'alphabet', 
        method = 'color', type = 'lower', col = COL1('YlOrBr', n=200), 
        tl.col = 'black', tl.cex = (size = 3), is.corr = F, addCoef.col = 'black'
        , addgrid.col = 'white')

# Load the data for entropy categories plot

flu_ent <- read_csv('flu_entropy_final.csv')
flu_ent$virus <- 'Influenza'
flu_domain <- read_csv('influenza_domain.csv')
flu_domain_new <- flu_domain %>% mutate(virus = 'Influenza')

ggplot() + geom_point(data = flu_ent, aes(position, entropy_x, color = label)) + 
  geom_line(data = flu_ent, aes(position, entropy_x), color = 'black') +
  scale_x_continuous(breaks = seq(0, 600, 25)) +
  #geom_line(data = flu_ent, aes(position, entropy_y), color = 'grey90') +
  geom_rect(data = flu_domain, aes(xmin= pos_start, xmax = pos_end, ymin=1, ymax=1.125, fill = domain)) +
  geom_rect(data = flu_domain2, aes(xmin= pos_start, xmax = pos_end, ymin=1.125, ymax=1.5, fill = domain)) +
  geom_text(data = flu_domain2, aes(x=pos_start+(pos_end-pos_start)/2, y = 1.375, label = domain), size = 8) +
  theme_bw()

sc2_ent <- read_csv('sc2_ent_final.csv')
sc2_domain <- read_csv('sars_domain.csv')
sc2_domain_new <- sc2_domain %>% mutate(virus = 'SARS-CoV-2')
sc2_ent$virus <- 'SARS-CoV-2'

ggplot() + geom_point(data = sc2_ent, aes(position, entropy_x, color = label)) + 
  geom_line(data = sc2_ent, aes(position, entropy_x), color = 'black') +
  scale_x_continuous(breaks = seq(0, 1250, 25)) +
  #geom_line(data = flu_ent, aes(position, entropy_y), color = 'grey90') +
  geom_rect(data = sc2_domain, aes(xmin= pos_start, xmax = pos_end, ymin=1, ymax=1.125, fill = domain)) +
  #geom_rect(data = sc2_domain2, aes(xmin= pos_start, xmax = pos_end, ymin=1.125, ymax=1.5, fill = domain)) +
  #geom_text(data = sc2_domain2, aes(x=pos_start+(pos_end-pos_start)/2, y = 1.375, label = domain), size = 8) +
  theme_bw()

hiv_ent <- read_csv('hiv_entropy_final.csv')
hiv_ent$virus <- 'HIV'
hiv_domain <- read_csv('hiv_domain.csv')
hiv_domain_new <- hiv_domain %>% mutate(virus = 'HIV')

ggplot() + geom_point(data = hiv_ent, aes(position, entropy_x, color = label)) + 
  geom_line(data = hiv_ent, aes(position, entropy_x), color = 'black') +
  scale_x_continuous(breaks = seq(0, 875, 25)) +
  #geom_line(data = flu_ent, aes(position, entropy_y), color = 'grey90') +
  geom_rect(data = hiv_domain, aes(xmin= pos_start, xmax = pos_end, ymin=1, ymax=1.125, fill = domain)) +
  #geom_rect(data = sc2_domain2, aes(xmin= pos_start, xmax = pos_end, ymin=1.125, ymax=1.5, fill = domain)) +
  #geom_text(data = sc2_domain2, aes(x=pos_start+(pos_end-pos_start)/2, y = 1.375, label = domain), size = 8) +
  theme_bw()

virus_ent_all <- bind_rows(flu_ent, sc2_ent, hiv_ent)
virus_ent_domain_all <- bind_rows(flu_domain_new1, sc2_domain_new, hiv_domain_new)
write_csv(virus_ent_all, 'virus_entropy_all.csv')
write_csv(virus_ent_domain_all, 'virus_entropy_domain_all.csv')

library(hrbrthemes)

ggplot() + 
  geom_point(data = virus_ent_all, aes(position, entropy_x, color = label), size = 3) + 
  scale_color_brewer(palette="Accent") +
  geom_line(data = virus_ent_all, aes(position, entropy_x), color = 'black') +
  scale_x_continuous(breaks = seq(0, 1250, 25)) +
  scale_y_continuous(breaks = seq(0, 1.25, 0.25)) +
  labs(color = 'Entropy groups') +
  geom_rect(data = virus_ent_domain_all, 
            aes(xmin= pos_start, xmax = pos_end, ymin=1, ymax=1.25, 
                fill = domain), show.legend = FALSE) +
  geom_text(data = virus_ent_domain_all, 
            aes(x=pos_start+(pos_end-pos_start)/2, y = 1.125, 
                label = domain), size = 4) +
  facet_wrap(~virus, ncol=1, strip.position = 'right') +
  labs(x='Position', y= 'Entropy', 
       title = 'Entropy of different viruses') +
  theme_ipsum(base_size = 14, 
              plot_title_size = 20, 
              axis_title_size = 22, 
              strip_text_size = 26) +
  theme(legend.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.position = 'bottom') +
  guides(color = guide_legend(override.aes = list(size = 10)))

#acc_ent plot:
acc_ent <- read_csv('flu_acc_entropy_final.csv')
ggplot(data = acc_ent) + 
  geom_rect(xmin= -Inf, xmax = Inf, ymin = 0.23, ymax = 0.51, 
            alpha = 0.05, fill = 'slategray1') + 
  geom_point(aes(position, rsa, color = label), size = 3) + 
  scale_color_brewer(palette="Set1") +
  geom_line(aes(position, rsa), color = 'black') +
  scale_x_continuous(breaks = seq(0,575, 25)) +
  scale_y_continuous(breaks = seq(0, 0.9, 0.1)) +
  geom_rect(data = flu_domain_new1, 
            aes(xmin= pos_start, xmax = pos_end, ymin=0.85, ymax=0.9, 
                fill = domain), show.legend = FALSE) +
  geom_text(data = flu_domain_new1, 
            aes(x=pos_start+(pos_end-pos_start)/2, y = 0.875, 
                label = domain), size = 8) +
  labs(x='Position', y= 'Accessibility', 
       title = 'Entropy and Surface Accessibility Pattern') +
  theme_ipsum(base_size = 14, 
              plot_title_size = 20, 
              axis_title_size = 22, 
              strip_text_size = 26) +
  theme(legend.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.position = 'bottom') +
  guides(color = guide_legend(override.aes = list(size = 10)))

ggplot(data = hiv_acc_ent) + 
  geom_rect(xmin= -Inf, xmax = Inf, ymin = 0.23, ymax = 0.51, 
            alpha = 0.05, fill = 'slategray1') + 
  geom_point(aes(position, rsa, color = label), size = 3) + 
  scale_color_brewer(palette="Set1") +
  geom_line(aes(position, rsa), color = 'black') +
  scale_x_continuous(breaks = seq(0,875, 25)) +
  scale_y_continuous(breaks = seq(0, 0.9, 0.1)) +
  geom_rect(data = hiv_domain_new, 
            aes(xmin= pos_start, xmax = pos_end, ymin=0.85, ymax=0.9, 
                fill = domain), show.legend = FALSE) +
  geom_text(data = hiv_domain_new, 
            aes(x=pos_start+(pos_end-pos_start)/2, y = 0.875, 
                label = domain), size = 8) +
  labs(x='Position', y= 'Accessibility', 
       title = 'HIV envelope protein') +
  theme_ipsum(base_size = 14, 
              plot_title_size = 20, 
              axis_title_size = 22, 
              strip_text_size = 26) +
  theme(legend.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.position = 'bottom') +
  guides(color = guide_legend(override.aes = list(size = 10)))


ggplot(data = sc2_acc_ent) + 
  geom_rect(xmin= -Inf, xmax = Inf, ymin = 0.23, ymax = 0.57, 
            alpha = 0.05, fill = 'slategray1') + 
  geom_point(aes(position, rsa, color = label), size = 3) + 
  scale_color_brewer(palette="Set1") +
  geom_line(aes(position, rsa), color = 'black') +
  scale_x_continuous(breaks = seq(0,1250, 50)) +
  scale_y_continuous(breaks = seq(0, 0.9, 0.1)) +
  geom_rect(data = sc2_domain_new[1:6,], 
            aes(xmin= pos_start, xmax = pos_end, ymin=0.85, ymax=0.9, 
                fill = domain), show.legend = FALSE) +
  geom_text(data = sc2_domain_new[1:6,], 
            aes(x=pos_start+(pos_end-pos_start)/2, y = 0.875, 
                label = domain), size = 4) +
  labs(x='Position', y= 'Accessibility', 
       title = 'SARS-CoV-2 spike protein') +
  theme_ipsum(base_size = 14, 
              plot_title_size = 20, 
              axis_title_size = 22, 
              strip_text_size = 26) +
  theme(legend.text = element_text(size = 26), 
        legend.title = element_text(size = 28), 
        legend.position = 'bottom') +
  guides(color = guide_legend(override.aes = list(size = 10)))

acc_ent

flu_escape <- read_csv('flu_h1_evescape_sites.csv')
flu_escape <- flu_escape %>% rename(position = i) %>% select(position, max_escape_experiment)

flu_escape_new <- acc_ent %>% left_join(flu_escape, by = 'position')
ggplot(data = flu_escape) + 
  geom_point(aes(position, max_escape, color = label), size = 3) + 
  scale_color_brewer(palette="Set1") +
  geom_line(aes(position, max_escape), color = 'black') +
  scale_x_continuous(breaks = seq(0,575, 25)) +
  scale_y_continuous(breaks = seq(0, 0.1, 0.01)) +
  geom_rect(data = flu_domain_new1, 
            aes(xmin= pos_start, xmax = pos_end, ymin=0.5, ymax=0.6, 
                fill = domain), show.legend = FALSE) +
  geom_text(data = flu_domain_new1, 
            aes(x=pos_start+(pos_end-pos_start)/2, y = 0.55, 
                label = domain), size = 8) +
  labs(x='Position', y= "Bloom's escape") +
  theme_ipsum(base_size = 14, 
              plot_title_size = 20, 
              axis_title_size = 22, 
              strip_text_size = 26) +
  theme(legend.text = element_text(size = 26),
        legend.title = element_text(size = 28), 
        legend.position = 'bottom') +
  guides(color = guide_legend(override.aes = list(size = 10)))
  


bloom <- read_csv('bloom_flu_h1.csv')
bloom <- bloom %>% rename(position = i)
bloom <- bloom %>% group_by(position) %>% 
  summarise(max_escape = max(antibody_H17L19_median_mutfracsurvive))
bloom

flu_aebe <- flu_acc_ent_escape %>% left_join(bloom, by = 'position')
flu_aebe

ggplot(data = flu_aebe %>% filter(max_escape <= 0.01)) + 
  geom_point(aes(position, max_escape, color = label), size = 3) + 
  scale_color_brewer(palette="Set1") +
  geom_line(aes(position, max_escape), color = 'black') +
  scale_x_continuous(breaks = seq(0,575, 25)) +
  scale_y_continuous(breaks = seq(0, 0.01, 0.001)) +
  geom_rect(data = flu_domain_new1, 
            aes(xmin= pos_start, xmax = pos_end, ymin=0.01, ymax=0.02, 
                fill = domain), show.legend = FALSE) +
  geom_text(data = flu_domain_new1, 
            aes(x=pos_start+(pos_end-pos_start)/2, y = 0.015, 
                label = domain), size = 8) +
  labs(x='Position', y= "Bloom's escape") +
  theme_ipsum(base_size = 14, 
              plot_title_size = 20, 
              axis_title_size = 22, 
              strip_text_size = 26) +
  theme(legend.text = element_text(size = 26),
        legend.title = element_text(size = 28), 
        legend.position = 'bottom') +
  guides(color = guide_legend(override.aes = list(size = 10)))

ggplot(data = hiv_all %>% filter(position >= 30 & position < 670)) + 
  geom_point(aes(position, Ab_escape, color = label.x), size = 3) + 
  scale_color_brewer(palette="Set1") +
  geom_line(aes(position, Ab_escape), color = 'black') +
  scale_x_continuous(breaks = seq(0,675, 25)) +
  scale_y_continuous(breaks = seq(0, 0.2, 0.02)) +
  geom_rect(data = hiv_domain_new[1:11,], 
            aes(xmin= pos_start, xmax = pos_end, ymin=0.2, ymax=0.21, 
                fill = domain), show.legend = FALSE) +
  geom_text(data = hiv_domain_new[1:11,], 
            aes(x=pos_start+(pos_end-pos_start)/2, y = 0.205, 
                label = domain), size = 6) +
  labs(x='Position', y= "Bloom's escape") +
  theme_ipsum(base_size = 14, 
              plot_title_size = 20, 
              axis_title_size = 22, 
              strip_text_size = 26) +
  theme(legend.text = element_text(size = 26),
        legend.title = element_text(size = 28), 
        legend.position = 'bottom') +
  guides(color = guide_legend(override.aes = list(size = 8)))

ggplot(data = sc2_all) + 
  geom_point(aes(position, Ab_escape, color = label.x), size = 3) + 
  scale_color_brewer(palette="Set1") +
  geom_line(aes(position, Ab_escape), color = 'black') +
  #scale_x_continuous(breaks = seq(0,675, 25)) +
  #scale_y_continuous(breaks = seq(0, 0.2, 0.02)) +
  # geom_rect(data = hiv_domain_new, 
  #           aes(xmin= pos_start, xmax = pos_end, ymin=0.2, ymax=0.3, 
  #               fill = domain), show.legend = FALSE) +
  # geom_text(data = hiv_domain_new, 
  #           aes(x=pos_start+(pos_end-pos_start)/2, y = 0.25, 
  #               label = domain), size = 8) +
  labs(x='Position', y= "Bloom's escape score", title = 'SARS-CoV-2 RBD') +
  theme_ipsum(base_size = 14, 
              plot_title_size = 20, 
              axis_title_size = 22, 
              strip_text_size = 26) +
  theme(legend.text = element_text(size = 26),
        legend.title = element_text(size = 28), 
        legend.position = 'bottom') +
  guides(color = guide_legend(override.aes = list(size = 8)))

hiv_acc_ent


#make the plot with all data 

flu_acc_ent_escape %>% ggplot() + 
  geom_point(aes(entropy, rsa, color = max_escape_experiment, size = entropy)) + 
  scale_color_distiller(palette = 'Spectral') + 
  theme_bw()

flu_plot <- flu_all %>% filter(position <= 350) %>% 
  ggplot() + 
  geom_point(aes(plm_entropy, 
                 discotope_score, 
                 color = Ab_escape, 
                 size = msa_entropy), 
             position = position_jitter(width=.04)) + 
  scale_color_distiller(palette = 'Spectral') +
  theme_ipsum(base_size = 20, 
              axis_title_size = 22)
flu_plot
 
hiv_all

hiv_plot <- hiv_all %>% filter(position >= 30 & position < 700) %>% 
  ggplot() + 
  geom_point(aes(plm_entropy, accessibility, 
                 color = Ab_escape, 
                 size = msa_entropy)) + 
  scale_color_distiller(palette = 'Spectral') +
  labs(x = 'plm_entropy', y = 'accessibility', 
       title = 'Single sequence prediction for HIV Env') +
  theme_ipsum(base_size = 30, 
              axis_title_size = 32)


sc2_plot <- sc2_all %>% 
  ggplot() + 
  geom_point(aes(plm_entropy, accessibility, 
                 color = Ab_escape, 
                 size = msa_entropy)) + 
  scale_color_distiller(palette = 'Spectral') +
  labs(x = 'Entropy', y = 'Accessibility', 
       title = 'Single sequence prediction for SARS-CoV-2 RBD') +
  theme_ipsum(base_size = 30, 
              axis_title_size = 32)


library(patchwork)
library(cowplot)

sc2_plot % flu_plot % hiv_plot

sc2_flu_hiv %>% 
  ggplot() + 
  geom_point(aes(plm_entropy, accessibility, 
                 color = Ab_escape, 
                 size = msa_entropy)) + 
  scale_color_distiller(palette = 'Spectral') +
  facet_wrap(~virus,scales = 'free', ncol = 1, strip.position = 'right') +
  labs(x = 'Entropy', y = 'Accessibility', 
       title = 'Single sequence prediction for SARS-CoV-2 RBD') +
  theme_ipsum(base_size = 30, 
              axis_title_size = 32)


#matching with real antibody escape

#sc2 rbd

sc2_rbd_escape <- read_csv('spike_rbd_evescape_voc_sera.csv')
sc2_rbd_abftprnt <- read_csv('sars2_spike_antibody_footprints.csv')

sc2_rbd_final_all <- sc2_rbd_escape %>% select(i, max_sera_strain) %>% arrange(i) %>% 
  drop_na() %>% 
  filter(max_sera_strain != 'No Escape') %>% 
  left_join(sc2_rbd_abftprnt, by = 'i') %>% 
  select(i, max_sera_strain,fraction_known_epitopes) %>% 
  rename(position = i) %>% 
  left_join(sc2_all_new, by = 'position')

sc2_rbd_final_all %>% 
  ggplot(aes(position)) +  
  geom_bar(aes(y=Ab_escape), stat = 'identity', fill = 'darkorange') +
  geom_point(aes(y=plm_entropy)) +
  geom_line(aes(y = accessibility),
            color = 'darkblue', size = 1)s

sc2_rbd_final_all %>%
  ggplot(aes(factor(position))) +  
  geom_bar(aes(y = Ab_escape), stat = 'identity', fill = 'darkorange') +
  geom_point(aes(y = accessibility)) +
  #geom_line(aes(y = accessibility), color = 'darkblue', size = 1) +
  xlab("Position") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



sc2_rbd_final_all %>%
  ggplot(aes(position)) +  
  geom_bar(aes(y = Ab_escape), stat = 'identity', fill = 'darkorange') +
  geom_point(aes(y = plm_entropy)) +
  geom_line(aes(y = accessibility), color = 'darkblue', size = 1) +
  scale_x_continuous(breaks = unique(sc2_rbd_final_all$position)) +
  xlab("Position") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

sc2_rbd_abftprnt %>% select(i, fraction_known_epitopes) %>% 
  arrange(desc(fraction_known_epitopes)) %>% 
  print(n=100)

sc2_eve %>% select(i, evescape) %>% rename(position = i) %>% 
  right_join(sc2_rbd_final_all, by = 'position') +
  ggplot(aes(position,))

#for hiv
hiv_all_new %>% filter(position %in% hiv_ab_escape_sites_losalamos) %>% 
  ggplot(aes(position)) +  
  geom_bar(aes(y = Ab_escape), stat = 'identity', fill = 'darkorange') +
  geom_point(aes(y = plm_entropy)) +
  geom_line(aes(y = accessibility), color = 'darkblue', size = 1)



#for flu
# HA1 positions in H1 influenza involved in antibody escape
flu_escape_sites <- c(129, 133, 135, 137, 153, 155, 159, 165, 166, 190, 193, 225, 226, 386, 390)

flu_all %>% filter(position %in% flu_escape_sites)
flu_all %>% filter(position %in% flu_escape_sites) %>% 
  ggplot(aes(position)) +  
  geom_bar(aes(y = Ab_escape*4), stat = 'identity', fill = 'darkorange') +
  geom_point(aes(y = plm_entropy)) +
  geom_line(aes(y = accessibility), color = 'darkblue', size = 1)



#perplexity thing
sc2_perp <- read_csv('sarscov2_new_perplexity.csv')
sc2_perp <- sc2_perp %>% rename(position = `...1`) %>% mutate (position = position + 1)
sc2_all %>% 
  left_join(sc2_perp, by = 'position') %>% 
  filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
  ggplot(aes(position, perplexity)) +
  geom_point(aes(size = Ab_escape), alpha = 0.7) +
  geom_text(aes(label = ifelse(position %in% sc2_abescape_pos, position, "")))
  
sc2_all %>% 
  left_join(sc2_perp, by = 'position') %>% 
  #filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
  ggplot(aes(accessibility, perplexity)) +
  scale_x_continuous(limit = c(0,0.7), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = Ab_escape, size = msa_entropy), alpha = 0.9) +
  scale_color_distiller(palette = 'Spectral') +
  geom_text(size =4.5,nudge_y =0.125,
            aes(label = ifelse(position %in% sc2_abescape_pos, position, ""))) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 18))
sc2_all_kieran %>% 
  left_join(sc2_perp, by = 'position') %>% 
  #filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
  ggplot(aes(Accessibility, perplexity)) +
  scale_x_continuous(limit = c(0,1), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = Ab_escape, size = msa_entropy), alpha = 0.9) +
  scale_color_distiller(palette = 'Spectral') +
  geom_text(size =4.5,nudge_y =0.125,
            aes(label = ifelse(position %in% sc2_abescape_pos, position, ""))) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 18))

sc2_disc_cal %>% 
  left_join(sc2_perp, by = 'position') %>% 
  #filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
  ggplot(aes(accessibility_disco_cal, perplexity)) +
  #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
  #scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = Ab_escape, size = msa_entropy), alpha = 0.9) +
  scale_color_distiller(palette = 'Spectral') +
  geom_text(size =4.5,nudge_y =0.125,
            aes(label = ifelse(position %in% sc2_abescape_pos, position, ""))) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 18))

sc2_all %>% 
  left_join(sc2_perp, by = 'position') %>% 
  filter(perplexity >= 3 & perplexity <= 6) %>% 
  ggplot(aes(accessibility, perplexity)) +
  scale_x_continuous(limit = c(0,0.7), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limit = c(3,6), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = Ab_escape, size = msa_entropy), alpha = 0.9) +
  scale_color_distiller(palette = 'Spectral') +
  geom_text(size =4,nudge_y =0.05,
            aes(label = ifelse(position %in% sc2_abescape_pos, position, ""))) +
  theme_classic()

mutate(position = as.character(position)) %>% 
  ggplot(aes(position)) +
  geom_bar(aes(y = perplexity), stat = 'identity', fill = 'grey70')+
  geom_point(aes(y = accessibility*6), color = 'darkorange', size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8), 
        axis.text.y = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        plot.title = element_text(size = 18)
        
    )
  
  ggplot(aes(position)) +  
  geom_bar(aes(y = accessibility), stat = 'identity', fill = 'darkorange') +
  #geom_point(aes(y = plm_entropy)) +
  geom_line(aes(y = perplexity), color = 'darkblue', size = 1) +
  scale_x_continuous(breaks = seq(0,1250, 50)) +
  labs(x='Position', y= 'Perplexity', 
       title = 'SARS-CoV-2 RBD perplexity') +
  theme_ipsum(base_size = 30, 
              axis_title_size = 32)

bind_rows(hiv_all %>% filter(!is.na(Ab_escape)), sc2_all_check, flu_all) %>% 
    ggplot(aes(accessibility, plm_entropy)) + 
    geom_point(aes(color = Ab_escape, size = msa_entropy)) + 
    scale_color_distiller(palette = 'Spectral') +
    theme_bw()

sc2_disc_cal %>% 
  left_join(sc2_perp, by = 'position') %>% 
  #filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
  ggplot(aes(accessibility_disco_cal, perplexity)) +
  #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
  #scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = Ab_escape, size = msa_entropy), alpha = 0.9) +
  scale_color_distiller(palette = 'Spectral') +
  geom_text(size =4.5,nudge_y =0.125,
            aes(label = ifelse(position %in% sc2_abescape_pos, position, ""))) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 18))


#hiv

hiv_perp <- read_csv('hiv_env_glycoprotein_new_perplexity.csv')
hiv_acc_disco3 <- read_csv('hiv_env_dis3.csv')
hiv_perp <- hiv_perp %>% rename(position = `...1`) %>% mutate (position = position + 1)
hiv_acc_disco3 <- hiv_acc_disco3 %>% rename(position = res_id)
hiv_env_variable_region <- hiv_perp %>% left_join(hiv_acc_disco3, by = 'position') %>% 
                                        right_join(hiv_all_new %>% filter(position %in% 
                                                                            hiv_hypervar_region), by = 'position')


hiv_env_variable_region
hiv_env_variable_region %>%
  #filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
  ggplot(aes(`DiscoTope-3.0_score`, plm_entropy)) +
  #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
  #scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = Ab_escape), alpha = 0.9) +
  scale_color_distiller(palette = 'Spectral') +
  geom_text(size = 8, nudge_y =0.01,
            aes(label = ifelse(position %in% hiv_ab_escape_sites_losalamos, position, ""))) +
  labs(x='Antigenicity (single aa residue)', y = 'PLM entropy', legend.title = 'Antibody escape (DMS score)') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 28),
        legend.title  = element_text(size = 28),
        plot.title = element_text(size = 18))

hiv_epi_all_final %>%
  drop_na() %>% 
  #filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
  ggplot(aes(max_epi_acc, Ab_escape)) +
  #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
  #scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = Ab_escape), alpha = 0.9) +
  scale_color_distiller(palette = 'Spectral') +
  geom_text(size = 6, nudge_y =0.01,
            aes(label = ifelse(position %in% hiv_ab_escape_sites_losalamos, position, ""))) +
  labs(x='Epistatic Antigenicity(combined effect of aa residues)', y = 'Experimental Deep Mutationl Scanning', legend.title = 'Antibody escape (DMS score)') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 28),
        legend.title  = element_text(size = 28),
        plot.title = element_text(size = 18))

hiv_epi_all_final %>%
  drop_na() %>% 
  #filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
  ggplot(aes(accessibility, Ab_escape)) +
  #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
  #scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = Ab_escape, size = Ab_escape), alpha = 0.2) +
  geom_point(aes())
  #geom_jitter() +
  scale_color_distiller(palette = 'Spectral') +
  geom_text(size = 4, nudge_y =0.005,
            aes(label = ifelse(position %in% hiv_ab_escape_sites_losalamos, position, ""))) +
  labs(x='Antigenicity (single aa effect only)', y = 'Experimental Deep Mutational Scanning', 
       legend.title = 'Antibody escape (DMS score)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 28),
        legend.title  = element_text(size = 28),
        plot.title = element_text(size = 18))
  
 
    #filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
ggplot() +
    #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
    #scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
    geom_point(data = hiv_epi_all_final %>% drop_na() %>% filter(!position %in% hiv_ab_escape_sites_losalamos), 
               aes(accessibility, Ab_escape, color = Ab_escape, size = Ab_escape),  
               alpha = 0.4) +
  geom_point(data = hiv_epi_all_final %>% drop_na() %>% filter(position %in% hiv_ab_escape_sites_losalamos), 
             aes(accessibility, Ab_escape, color = Ab_escape, size = Ab_escape), 
             shape = 19) +
  geom_vline(xintercept = mean(hiv_epi_all_final$accessibility, na.rm=TRUE)) +
  #geom_jitter() +
  scale_color_distiller(palette = 'Spectral') +
    geom_text(data = hiv_epi_all_final %>% drop_na(), size = 6, fontface = 'bold', nudge_y =0.005,
              aes(accessibility, Ab_escape,label = ifelse(position %in% hiv_ab_escape_sites_losalamos, position, ""))) +
    labs(x='Antigenicity (single aa effect only)', y = 'Experimental Deep Mutational Scanning', 
         legend.title = 'Antibody escape (DMS score)') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
          axis.text.y = element_text(size = 18), 
          axis.title = element_text(size = 28),
          legend.title  = element_text(size = 28),
          plot.title = element_text(size = 18))

ggplot() +
  #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
  #scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
  geom_point(data = hiv_epi_all_final %>% drop_na() %>% filter(!position %in% hiv_ab_escape_sites_losalamos), 
             aes(accessibility, Ab_escape, fill = Ab_escape, size = Ab_escape), shape = 21, color = 'white',
             alpha = 0.3) +
  geom_point(data = hiv_epi_all_final %>% drop_na() %>% filter(position %in% hiv_ab_escape_sites_losalamos), 
             aes(accessibility, Ab_escape,  fill = Ab_escape, size = Ab_escape), color = 'black',
             shape = 21) +
  geom_vline(xintercept = median(hiv_epi_all_final$accessibility, na.rm=TRUE), color = 'black', lwd = 1.5) +
  #geom_jitter() +
  scale_fill_distiller(palette = 'Spectral') +
  geom_text(data = hiv_epi_all_final %>% drop_na(), size = 6, fontface = 'bold', nudge_y =0.005,
            aes(accessibility, Ab_escape,label = ifelse(position %in% hiv_ab_escape_sites_losalamos, position, ""))) +
  labs(x='Antigenicity (single aa effect only)', y = 'Experimental Deep Mutational Scanning', 
       legend.title = 'Antibody escape (DMS score)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 28),
        legend.title  = element_text(size = 28),
        plot.title = element_text(size = 18))


ggplot() +
  #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
  #scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
  geom_point(data = hiv_epi_all_final %>% drop_na() %>% filter(!position %in% hiv_ab_escape_sites_losalamos), 
             aes(max_epi_acc, Ab_escape, fill = Ab_escape, size = Ab_escape), color = 'white', shape = 21, 
             alpha = 0.3) +
  geom_point(data = hiv_epi_all_final %>% drop_na() %>% filter(position %in% hiv_ab_escape_sites_losalamos), 
             aes(max_epi_acc, Ab_escape, fill = Ab_escape, size = Ab_escape), color = 'black', shape = 21) +
  geom_vline(xintercept = median(hiv_epi_all_final$max_epi_acc, na.rm=TRUE), color = 'black', lwd = 1.5) +
  #geom_jitter() +
  scale_fill_distiller(palette = 'Spectral') +
  geom_text(data = hiv_epi_all_final %>% drop_na(), size = 6, fontface = 'bold', nudge_y =0.002,
            aes(max_epi_acc, Ab_escape,label = ifelse(position %in% hiv_ab_escape_sites_losalamos, position, ""))) +
  labs(x='Epistatic Antigenicity(combined effect of aa residues)', y = 'Experimental Deep Mutational Scanning', 
       legend.title = 'Antibody escape (DMS score)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 28),
        legend.title  = element_text(size = 28),
        plot.title = element_text(size = 18))
#epistatic thing: kieran

acc_epi_max <- acc_epi_kieran %>% group_by(pos) %>% select(pos, starts_with('very_significant_')) %>% 
  summarise(max_surface_acc= max(very_significant_accessibility, na.rm=T), 
            max_bepro_acc = max(very_significant_BEPro_accessibility, na.rm = T), 
            max_disco_acc = max(very_significant_discotope, na.rm=T)) %>% 
  rename(position = pos)

sc2_epi_all <- sc2_all %>% left_join(acc_epi_max, by = 'position')
sc2_epi_all %>% 
  left_join(sc2_perp, by = 'position') %>% 
  #filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
  ggplot(aes(max_bepro_acc, perplexity)) +
  #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
  #scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = Ab_escape, size = msa_entropy), alpha = 0.9) +
  scale_color_distiller(palette = 'Spectral') +
  geom_text(size =4.5,nudge_y =0.125,
            aes(label = ifelse(position %in% sc2_abescape_pos, position, ""))) +
  labs(title = 'SARS-CoV-2 RBD', x = 'antibody accessiblity score (Epistatis based)') + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 18))

sc2_epi_all %>% 
  left_join(sc2_perp, by = 'position') %>% 
  filter(perplexity >= 4 & perplexity <= 7.5 & max_bepro_acc >= 0.7 & position %in% sc2_abescape_pos) 
  
sc2_all 
sc2_all %>% select(-c(label.x, label.y, virus)) %>% 
  gather(para, value, -position) %>% 
  ggplot(aes(x=position, y= value)) + 
  geom_point() + 
  facet_wrap(~para, nrow = 4) +
  theme_bw()
sc2_epi_all %>% 
  left_join(sc2_perp, by = 'position') %>% 
  select(-c(label.x, label.y, virus)) %>% 
  gather(para, value, -position) %>% 
  ggplot(aes(x=position, y= value)) + 
  geom_point() + 
  facet_wrap(~para, nrow = 8, scale = 'free', strip.position="right") +
  theme_bw()

flu_perp <- flu_perp %>% rename(position = `...1`)
flu_true_ent %>% 
  left_join(flu_all, by = 'position') %>% 
  


#flu
flu_all %>% 
  #left_join(sc2_perp, by = 'position') %>% 
  #filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
  ggplot(aes(`DiscoTope-3.0_score`, plm_entropy)) +
  #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
  #scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = Ab_escape, size = msa_entropy), alpha = 0.9) +
  scale_color_distiller(palette = 'Spectral') +
  geom_text(size =4.5,
            aes(label = ifelse(position %in% flu_real$real, position, ""))) +
  labs(title = 'Influenza H1 HA protein', x = 'antibody accessiblity score (Discotope)', y = "PLM Entropy") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 18)
        )



#epistatic accessibility
#flu

flu_epi_acc <- read_csv('flu_ha_epistatic_accessibility.csv') %>% select(position, ref, significant_accessibility) %>% 
  group_by(position) %>% summarise(max_epi_acc = mean(significant_accessibility, na.rm = TRUE))
flu_epi_all_final <- flu_all %>% left_join(flu_epi_acc, by = 'position')

flu_epi_all_final %>% 
  left_join(flu_perp, by = 'position') %>% 
  #filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
  # ggplot(aes(max_epi_acc, perplexity)) +
  ggplot(aes(max_epi_acc, perplexity)) +
  #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
  scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = Ab_escape)) +
  scale_color_distiller(palette = 'Spectral') +
  geom_text(size =4.5,
            aes(label = ifelse(position %in% flu_real$real, position, ""))) +
  labs(title = 'Influenza H1 HA protein', x = 'Epistatic Antigenicity', y = "PLM-perplexity") +
  # labs(title = 'Influenza H1 HA protein', x = 'Antigenicity of individual amino-acid (Discotope)', y = "PLM-perplexity") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 18)
  )


#hiv
hiv_epi_acc <- read_csv('hiv_env_epistatic_accessibility.csv') %>% select(position, ref, significant_accessibility) %>% 
  group_by(position) %>% summarise(max_epi_acc = mean(significant_accessibility, na.rm = TRUE))
hiv_epi_all_final <- hiv_all %>% left_join(hiv_epi_acc, by = 'position')


hiv_ultimate %>% 
  #left_join(hiv_plm_perplexity, by = 'position') %>% 
  #filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
  ggplot(aes(max_escape_experiment, Epistatic_AntigenicitySignificant)) +
  #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
  #scale_y_continuous(limit = c(0,15), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = max_escape_experiment), alpha = 0.9) +
  scale_color_distiller(palette = 'Spectral') +
  geom_text(size =4.5,
            aes(label = ifelse(position %in% hiv_ab_escape_sites_losalamos, position, ""))) +
  labs(title = 'HIV Env protein', x = 'Epistatic antigenicity', y = "Antibody escape (Bloom)") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 18)
  )

hiv_epi_all_final %>% 
  left_join(hiv_perp, by = 'position') %>% 
  #filter(perplexity >= 3 & perplexity <= 6 & accessibility >= 0) %>% 
  ggplot(aes(accessibility, plm_entropy)) +
  #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
  #scale_y_continuous(limit = c(0,13.5), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = Ab_escape), alpha = 0.9) +
  scale_color_distiller(palette = 'Spectral') +
  geom_text(size =4.5,
            aes(label = ifelse(position %in% hiv_ab_escape_sites_losalamos, position, ""))) +
  labs(title = 'HIV Env protein', x = 'Epistatic antigenicity', y = "PLM Entropy") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 18)
  )

hiv_ultimate %>% 
  #left_join(hiv_plm_perplexity, by = 'position') %>% 
  #filter(label %in% ) %>% 
  ggplot(aes(Epistatic_AntigenicitySignificant, probability)) +
  #scale_x_continuous(limit = c(0,0.4), breaks = seq(0, 0.4, 0.05)) +
  #scale_y_continuous(limit = c(0,15), breaks = seq(0, 15, 1.5)) +
  geom_point(aes(color = max_escape_experiment), alpha = 0.9) +
  scale_color_distiller(palette = 'Spectral') +
  # geom_text(size =4.5,
  #           aes(label = ifelse(position %in% hiv_ab_escape_sites_losalamos, position, ""))) +
  labs(title = 'HIV Env protein', x = 'Epistatic antigenicity', y = "Antibody escape (Bloom)") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 18)
  )
colnames(hiv_ultimate) == is.numeric
cor.test(hiv_ultimate$max_escape_experiment, hiv_ultimate$Epistatic_AntigenicitySignificant, method = 'pearson')


#ranking epistasis
sc2_df <- sc2_epi_rank_check
sc2_df$position[sort(sc2_df$max_epi_access)]
sc2_df$rank_epi <- rank(sc2_df$max_epi_access, ties.method = "first", na.last = "keep")
sc2_df$rank_acc <- rank(sc2_df$Accessibility, ties.method = "first", na.last = "keep")
sc2_df[order(sc2_df$rank_epi), ]
sc2_df %>% filter(position %in% sc2_pos) %>% mutate(rank_delta = rank_acc - rank_epi) %>% arrange(rank_delta)
sc2_man <- sc2_df %>% filter(position %in% sc2_abescape_pos)  %>% select(rank_epi, rank_acc, position) %>% pivot_longer(cols = c(rank_epi, rank_acc))
sc2_man <- sc2_df %>% filter(position %in% sc2_abescape_pos)
sc2_man1 <- sc2_df %>% filter(position>330 & position < 530) 
delta_rank <- sc2_man1 %>% mutate(delta_rank = rank_epi - rank_acc)
sum(delta_rank$delta_rank)
cor.test(sc2_df$rank_epi, sc2_df$rank_acc, method = 'kendall')

hiv_corr <- Filter(is.numeric, hiv_ultimate) 
hiv_corr %>% select(!starts_with('summary_')) %>% 
  select(!c(7, 13, 15, 18, 24)) %>% 
  drop_na() %>% as.matrix %>% 
  cor %>% 
  corrplot(order="alphabet", 
           tl.color = 'black', 
           diag = FALSE,
           method = 'color', 
           col=rev(brewer.pal(n=10, name="RdYlBu")))

hiv_weight <- read_csv('/Users/rubayetalam/hiv_bg505_env_epistatic_accessibility_weighted .csv')
hiv_premium <- hiv_ultimate %>% 
                left_join(hiv_weight %>% select(!c(ref, alt, position)),  
                          by = 'label')
colnames(hiv_premium)
hiv_premium

hiv_corra <- Filter(is.numeric, hiv_premium) %>% filter(position %in% hiv_ab_escape_sites_losalamos)
hiv_corrraa <- hiv_premium %>% filter(label %in% hiv_all_exp_mut$label)  
Filter(is.numeric, hiv_corrraa)
hiv_corra %>% select(!starts_with('summary_')) %>% 
  select(!c(7, 13, 15, 18, 24)) %>% 
  drop_na() %>% as.matrix %>%
  cor %>% 
  corrplot(order="alphabet", 
           #tl.color = 'black', 
           diag = FALSE,
           method = 'number', 
           col=rev(brewer.pal(n=10, name="RdYlBu")))


hiv_corra <- Filter(is.numeric, hiv_premium) 
hiv_corra %>% select(!starts_with('summary_')) %>% 
  select(!c(7, 13, 15, 18, 24)) %>% 
  drop_na() %>% 
  filter(max_escape_experiment >= 0.2) %>% 
  as.matrix %>%
  cor %>% 
  corrplot(order="alphabet", 
           #tl.color = 'black', 
           diag = FALSE,
           method = 'number', 
           col=rev(brewer.pal(n=10, name="RdYlBu")))

Filter(is.numeric, hiv_corrraa) %>% 
  select(!starts_with('summary_')) %>% 
  select(!c(7, 13, 15, 18, 24)) %>% 
  drop_na() %>% 
  filter(max_escape_experiment >= 0.05) %>% 
  as.matrix %>%
  cor %>% 
  corrplot(order="alphabet", 
           #tl.color = 'black', 
           diag = FALSE,
           method = 'color', 
           col=rev(brewer.pal(n=10, name="RdYlBu")))


Filter(is.numeric, hivbg505_df_final) %>% 
  select(!starts_with('summary_')) %>% 
  #select(!c(7, 13, 15, 18, 24)) %>% 
  drop_na() %>% 
  #filter(max_escape_experiment >= 0.05) %>% 
  as.matrix %>%
  cor %>% 
  corrplot(order="alphabet", 
           #tl.color = 'black', 
           diag = FALSE,
           method = 'color', 
           col = rev(brewer.pal(n=10, name="RdYlBu"))
           )

summary(hiv_premium$Epistatic_AntigenicitySignificant)
binom.test(sum(delta_rank$delta_rank>0), length(delta_rank$delta_rank), p=0.5, alternative="greater")

Filter(is.numeric, hiv_corrraa) %>% 
  select(!starts_with('summary_')) %>% 
  select(!c(7, 13, 15, 18, 24)) %>% 
  select(c(probability, very_significant_accessibility_sum, significant_accessibility_median, SingleAA_Antigenicity, aa_prob, evescape, max_escape_experiment)) %>% 
  drop_na() %>% 
  #filter(max_escape_experiment >= 0.05) %>% 
  as.matrix %>%
  cor %>% 
  corrplot(order="alphabet", 
           #tl.color = 'black', 
           diag = FALSE,
           type = 'upper',
           method = 'number', 
           col=rev(brewer.pal(n=10, name="RdYlBu")))

Filter(is.numeric, hiv_corrraa) %>% 
  select(!starts_with('summary_')) %>% 
  select(!c(7, 13, 15, 18, 24)) %>% 
  select(c(probability, very_significant_accessibility_sum, SingleAA_Antigenicity, aa_prob, evescape, max_escape_experiment)) %>% 
  drop_na() %>% 
  rename('MSA entropy' = aa_prob,
         'Antigenicity (single aa)' = SingleAA_Antigenicity,
         'Epistatic Antigenicity' = very_significant_accessibility_sum,
         'PLM entropy' = probability,
         'DMS score' = max_escape_experiment,
         'EVEscape value' = evescape) %>% 
  select('Epistatic Antigenicity', 'Antigenicity (single aa)', 'DMS score', 'EVEscape value') %>% 
  #filter(max_escape_experiment >= 0.05) %>% 
  as.matrix %>%
  cor %>% 
  corrplot(#order="alphabet", 
           tl.col = 'black', 
           diag = FALSE,
           #type = 'upper',
           method = 'color',
           addCoef.col = 'black',
           cl.cex = 1,
           tl.cex = 2,
           pch.cex = 1.5,
           col=rev(brewer.pal(n=10, name="RdYlBu")))

hiv_corrraa

wilcox.test(value ~ name, data = sc2_man, distribution = "exact")
wilcox.test(sc2_df$Accessibility, sc2_df$max_epi_access)
wilcox.test(sc2_man$rank_epi,  sc2_man$rank_acc,  paired = TRUE, exact = F, alternative = "greater")


#hiv rank test
hiv_df$rank_acc <- rank(hiv_df$disco_acc)
hiv_imp <- hiv_df %>% filter(position %in% hiv_ab_escape_sites_losalamos) %>% mutate(rank_delta = rank_epi - rank_acc) %>% print(n=23)
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


flu_new23 <- flu_epi_all_final %>% select(position, `DiscoTope-3.0_score`, max_epi_acc) %>% 
  mutate(across(c(`DiscoTope-3.0_score`, max_epi_acc), ~ rank(.x, ties.method = "min"), .names = "{.col}_rank")) %>% 
  filter(position %in% str_extract(flu_wsn33_true_antibody_escape$label, '\\d+')) 
wilcox.test(flu_new23$max_epi_acc_rank, flu_new23$`DiscoTope-3.0_score_rank`, paired = TRUE, exact = F, alternative = "greater")
  
flu_new23$delta_rank <- flu_new23$max_epi_acc_rank - flu_new23$`DiscoTope-3.0_score_rank`
length(which(flu_new23$delta_rank < 0))


flu_new24 <- flu_busy2 %>% select(label, Antigenicity_sequence, Epistatic_Antigenicity_sequence) %>% 
  mutate(position = str_extract(label, '\\d+')) %>% 
  group_by(position) %>% 
  mutate(Epistatic_Antigenicity_sequence = max(Epistatic_Antigenicity_sequence)) %>%
  ungroup() %>% 
  mutate(across(c(Antigenicity_sequence, Epistatic_Antigenicity_sequence), ~ rank(.x, ties.method = "min"), .names = "{.col}_rank")) %>% 
  filter(label %in% flu_wsn33_true_antibody_escape$label) %>%  
  select(ends_with('rank'), position) %>% 
  distinct(position, .keep_all = T) %>%
  mutate(delta_rank = Epistatic_Antigenicity_sequence_rank - Antigenicity_sequence_rank)
  
  # group_by(label) %>% 
  # mutate(Epistatic_Antigenicity_sequence_rank = max(Epistatic_Antigenicity_sequence_rank)) %>% 
  # ungroup() %>% 
  # distinct(label, .keep_all = T) %>% 
  # mutate(rank_delta = Epistatic_Antigenicity_sequence_rank - Antigenicity_sequence_rank) %>% 
wilcox.test(flu_new24$Epistatic_Antigenicity_sequence_rank, flu_new24$Antigenicity_sequence_rank, paired = TRUE, exact = F, alternative = "greater")

flu_new25 <- flu_busy2 %>% select(label, Antigenicity_sequence, Epistatic_Antigenicity_sequence) %>% 
  mutate(position = str_extract(label, '\\d+')) %>% 
  group_by(position) %>% 
  mutate(Epistatic_Antigenicity_sequence = max(Epistatic_Antigenicity_sequence)) %>%
  ungroup() %>% 
  mutate(across(c(Antigenicity_sequence, Epistatic_Antigenicity_sequence), ~ rank(.x, ties.method = "min"), .names = "{.col}_rank")) %>% 
  filter(label %in% flu_wsn33_true_antibody_escape$label) %>%  
  select(ends_with('rank'), position) %>% 
  # distinct(position, .keep_all = T) %>%
  mutate(delta_rank = Epistatic_Antigenicity_sequence_rank - Antigenicity_sequence_rank)

sum(flu_new24$delta_rank)

wilcox.test(flu_new25$Epistatic_Antigenicity_sequence_rank, flu_new25$Antigenicity_sequence_rank, paired = TRUE, exact = F, alternative = "greater")



flu_new26 <- flu_busy2 %>% select(label, Antigenicity_sequence, Epistatic_Antigenicity_sequence) %>% 
  # mutate(position = str_extract(label, '\\d+')) %>% 
  # group_by(position) %>% 
  # mutate(Epistatic_Antigenicity_sequence = max(Epistatic_Antigenicity_sequence)) %>%
  # ungroup() %>% 
  mutate(across(c(Antigenicity_sequence, Epistatic_Antigenicity_sequence), ~ rank(.x, ties.method = "min"), .names = "{.col}_rank")) %>% 
  filter(label %in% flu_wsn33_true_antibody_escape$label) %>%  
  select(ends_with('rank')) %>% 
  # distinct(position, .keep_all = T) %>%
  mutate(delta_rank = Epistatic_Antigenicity_sequence_rank - Antigenicity_sequence_rank)

sum(flu_new26$delta_rank == 0)
wilcox.test(flu_new26$Epistatic_Antigenicity_sequence_rank, flu_new26$Antigenicity_sequence_rank, paired = TRUE, exact = F, alternative = "greater")

wilcox.test(flu_new27$Epistatic_Antigenicity_sequence_rank, flu_new27$Antigenicity_sequence_rank, paired = TRUE, exact = F, alternative = "greater")

flu_new29 <- flu_busy2 %>% select(label, Antigenicity_structure, Epistatic_Antigenicity_structure, `DMS score`) %>% 
  mutate(position = str_extract(label, '\\d+')) %>% 
  group_by(position) %>% 
  mutate(Epistatic_Antigenicity_structure = max(Epistatic_Antigenicity_structure, na.rm = T)) %>%
  ungroup() %>%
  # mutate(across(c(Antigenicity_structure, Epistatic_Antigenicity_structure), ~ rank(.x, ties.method = "min"), .names = "{.col}_rank")) %>% 
  # filter(label %in% flu_wsn33_true_antibody_escape$label) %>%  
  # select(ends_with('rank'), position) %>% 
  distinct(position, .keep_all = T) 
  # mutate(delta_rank = Epistatic_Antigenicity_structure_rank - Antigenicity_structure_rank)


flu_new30 <- flu_new29 %>% 
  mutate(across(c(Antigenicity_structure, Epistatic_Antigenicity_structure), ~ rank(.x, ties.method = "min"), .names = "{.col}_rank")) %>% 
  distinct(position, .keep_all = T) %>% 
  filter(position %in% as.numeric(names(table(str_extract(flu_wsn33_true_antibody_escape$label, '\\d+'))))) %>% 
  mutate(delta_rank = Epistatic_Antigenicity_structure_rank - Antigenicity_structure_rank)


flu_new31 <- flu_busy2 %>% select(label, Antigenicity_sequence, Epistatic_Antigenicity_sequence, `DMS score`) %>% 
  mutate(position = str_extract(label, '\\d+')) %>% 
  group_by(position) %>% 
  mutate(Epistatic_Antigenicity_sequence = max(Epistatic_Antigenicity_sequence, na.rm = T)) %>%
  ungroup() %>%
  # mutate(across(c(Antigenicity_structure, Epistatic_Antigenicity_structure), ~ rank(.x, ties.method = "min"), .names = "{.col}_rank")) %>% 
  # filter(label %in% flu_wsn33_true_antibody_escape$label) %>%  
  # select(ends_with('rank'), position) %>% 
  distinct(position, .keep_all = T) 


flu_new32 <- flu_new31 %>% 
  mutate(across(c(Antigenicity_sequence, Epistatic_Antigenicity_sequence), ~ rank(.x, ties.method = "min"), .names = "{.col}_rank")) %>% 
  distinct(position, .keep_all = T) %>% 
  filter(position %in% as.numeric(names(which(table(str_extract(flu_wsn33_true_antibody_escape$label, '\\d+')) > 1)))) %>% 
  mutate(delta_rank = Epistatic_Antigenicity_sequence_rank - Antigenicity_sequence_rank)
# mutate(delta_rank = Epistatic_Antigenicity_structure_rank - Antigenicity_structure_rank)
wilcox.test(flu_new32$Epistatic_Antigenicity_sequence_rank, flu_new32$Antigenicity_sequence_rank, paired = TRUE, exact = F, alternative = "greater")


Filter(is.numeric, flu_wsn33_final) %>% 
  #select(!starts_with('summary_')) %>% 
  #select(!c(7, 13, 15, 18, 24)) %>% 
  #drop_na() %>% 
  #filter(max_escape_experiment >= 0.05) %>% 
  as.matrix %>%
  cor %>% 
  corrplot(order="alphabet", 
           #tl.color = 'black', 
           diag = FALSE,
           method = 'color', 
           col=rev(brewer.pal(n=10, name="RdYlBu")))
colnames(flu_wsn33_final)


FluWSN <- flu_wsn33_final %>% 
  filter(label %in% flu_wsn33_true_antibody_escape$label) 
FluWSN

Filter(is.numeric, FluWSN) %>% 
  select(probability, significant_accessibility, single_aa_antigenicity, msa_aa_prob, evescape, max_escape_experiment) %>% 
  drop_na() %>% 
  as.matrix %>% 
  cor %>% 
  corrplot(order="alphabet", 
           #tl.color = 'black', 
           diag = FALSE,
           type= 'upper',
           method = 'number', 
           col=rev(brewer.pal(n=10, name="RdYlBu")))

 write_csv(flu_wsn33_final, 'flu_wsn33_final_all.csv') 
 
 
 Filter(is.numeric, SC2_Wuhan_rbd) %>% 
   select(probability, very_significant_accessibility, single_aa_antigenicity, msa_aa_prob, evescape, max_escape_experiment_bloom) %>% 
   drop_na() %>% 
   as.matrix %>% 
   cor %>% 
   corrplot(order="alphabet", 
            #tl.color = 'black', 
            diag = FALSE,
            type = 'upper',
            method = 'number', 
            col=rev(brewer.pal(n=10, name="RdYlBu")))
 
 
 Filter(is.numeric, sc2_high_abescape_final) %>% 
   #select(!starts_with('summary_')) %>% 
   #select(!c(7, 13, 15, 18, 24)) %>% 
   select(c(probability, very_significant_accessibility, single_aa_antigenicity, msa_aa_prob, evescape, max_escape_experiment)) %>% 
   drop_na() %>% 
   rename('MSA entropy' = msa_aa_prob,
          'Antigenicity (single aa)' = single_aa_antigenicity,
          'Epistatic Antigenicity' = very_significant_accessibility,
          'PLM entropy' = probability,
          'DMS score' = max_escape_experiment,
          'EVEscape value' = evescape) %>% 
   #filter(max_escape_experiment >= 0.05) %>% 
   as.matrix %>%
   cor %>% 
   corrplot(order="alphabet", 
            tl.col = 'black', 
            diag = FALSE,
            #type = 'upper',
            method = 'color',
            addCoef.col = 'black',
            cl.cex = 1,
            tl.cex = 2,
            pch.cex = 1.5,
            col=rev(brewer.pal(n=10, name="RdYlBu")))
 
 
 Filter(is.numeric, sc2_fullspike_final_v2) %>% 
   select(!starts_with('summary_')) %>% 
   #select(!c(7, 13, 15, 18, 24)) %>% 
   #select(c(probability, very_significant_accessibility_sum, significant_accessibility_median, SingleAA_Antigenicity, aa_prob, evescape, max_escape_experiment)) %>% 
   drop_na() %>% 
   #filter(max_escape_experiment >= 0.05) %>% 
   as.matrix %>%
   cor %>% 
   corrplot(order="alphabet", 
            #tl.color = 'black', 
            diag = FALSE,
            type = 'upper',
            method = 'color', 
            col=rev(brewer.pal(n=10, name="RdYlBu")))
 
 
 sc2c <- Filter(is.numeric, sc2_high_abescape_final) %>% mutate(virus = 'SARS-CoV-2 Spike') %>% 
   select(c(virus, probability, very_significant_accessibility, single_aa_antigenicity, msa_aa_prob, evescape, max_escape_experiment)) %>% 
   drop_na() %>% 
   rename('MSA probability' = msa_aa_prob,
          'Antigenicity (single aa)' = single_aa_antigenicity,
          'Epistatic Antigenicity' = very_significant_accessibility,
          'PLM probability' = probability,
          'DMS score' = max_escape_experiment,
          'EVEscape value' = evescape)
 
 hivc <- Filter(is.numeric, hiv_corrraa) %>% mutate(virus = 'HIV Env') %>% select(!starts_with('summary_')) %>% 
   select(!c(7, 13, 15, 18, 24)) %>% 
   select(c(virus, probability, very_significant_accessibility_sum, SingleAA_Antigenicity, aa_prob, evescape, max_escape_experiment)) %>% 
   drop_na() %>% 
   rename('MSA probability' = aa_prob,
          'Antigenicity (single aa)' = SingleAA_Antigenicity,
          'Epistatic Antigenicity' = very_significant_accessibility_sum,
          'PLM probability' = probability,
          'DMS score' = max_escape_experiment,
          'EVEscape value' = evescape) 
   
 
 
 fluc <- Filter(is.numeric, FluWSN) %>% mutate(virus = 'Influenza H1 HA') %>% 
   select(c(virus, probability, significant_accessibility, single_aa_antigenicity, msa_aa_prob, evescape, max_escape_experiment)) %>% 
   drop_na() %>% 
   rename('MSA probability' = msa_aa_prob,
          'Antigenicity (single aa)' = single_aa_antigenicity,
          'Epistatic Antigenicity' = significant_accessibility,
          'PLM probability' = probability,
          'DMS score' = max_escape_experiment,
          'EVEscape value' = evescape)
   
   
tvc <- bind_rows(sc2c, hivc, fluc)
   
tvcs <- split(tvc[,2:7], tvc$virus)
tvcm <- lapply(tvcs, function(x) cor(as.matrix(x)))
names(tvcm[1])

par(mfrow=c(1,3))
lapply(seq_along(tvcm), function(x){
  corrplot(x, order="alphabet", 
             tl.col = 'black', 
             diag = FALSE,
             #type = 'upper',
             method = 'color',
             addCoef.col = 'black',
             cl.cex = 2,
             tl.cex = 2,
             pch.cex = 1.5,
             col=rev(brewer.pal(n=10, name="RdYlBu")))
  mtext(paste(x),line=1,side=3)})
  
names(tvcm)
par(mfrow=c(3,1))
lapply(seq_along(tvcm), function(x){
  corrplot(tvcm[[x]], order="alphabet", 
           tl.col = 'black', 
           diag = FALSE,
           type = 'upper',
           method = 'color',
           addCoef.col = 'black',
           cl.cex = 1,
           tl.cex = 2,
           pch.cex = 1.5,
           col=rev(brewer.pal(n=10, name="RdYlBu")))
  #mtext(paste(names(tvcm[x])),line=2,side=3, adj = 0.75, cex =2)
  })

#filter(max_escape_experiment >= 0.05) %>% 
tvc %>% select(!virus) %>% 
   as.matrix %>%
   cor %>% 
  
  
hiv_ultimate %>% filter(label %in% hiv_all_exp_mut)


Filter(is.numeric, fl_new) %>% 
  select(!starts_with('summary_')) %>% 
  #select(!c(7, 13, 15, 18, 24)) %>% 
  drop_na() %>% 
  #filter(max_escape_experiment >= 0.05) %>% 
  as.matrix %>%
  cor %>% 
  corrplot(order="alphabet",
           main = 'Influenza H1N1 HA protein',
           sub = 'Correlation between metrics linked to antigenicity and site variation',
           #tl.color = 'black', 
           diag = FALSE,
           method = 'color', 
           col=rev(brewer.pal(n=10, name="RdYlBu")))

Filter(is.numeric, hivbg505_df_final) %>% 
  select(!starts_with('summary_')) %>% 
  #select(!c(7, 13, 15, 18, 24)) %>% 
  drop_na() %>% 
  #filter(max_escape_experiment >= 0.05) %>% 
  as.matrix %>%
  cor %>% 
  corrplot(order="alphabet",
           main = 'HIV (BG505) Envelope protein',
           sub = 'Correlation between metrics linked to antigenicity and site variation',
           #tl.color = 'black', 
           diag = FALSE,
           method = 'color', 
           col=rev(brewer.pal(n=10, name="RdYlBu")))
 

Filter(is.numeric, sc2_fullspike_final_v2) %>% 
  select(!starts_with('summary_')) %>% 
  #select(!c(7, 13, 15, 18, 24)) %>% 
  drop_na() %>% 
  #filter(max_escape_experiment >= 0.05) %>% 
  as.matrix %>%
  cor %>% 
  corrplot(order="alphabet",
           main = 'SaRS-CoV-2 Spike protein',
           sub = 'Correlation between metrics linked to antigenicity and site variation',
           #tl.color = 'black', 
           diag = FALSE,
           method = 'color', 
           col=rev(brewer.pal(n=10, name="RdYlBu")))


#sc2 aa prob check


sc2_baal2 %>% distinct(label, .keep_all = T)  %>%  filter(ref.x != alt.x) %>% select(probability, msa_aa_prob) %>% drop_na() %>% 
  filter(msa_aa_prob > 0.05, msa_aa_prob < 0.25) %>% cor()
sc2_baal2 %>% distinct(label, .keep_all = T)  %>%  
  filter(ref.x != alt.x) %>% 
  filter(label %in% sc2_tre_antibody_escape$label) %>% 
  select(probability, msa_aa_prob) %>% 
  drop_na() %>% 
  filter(msa_aa_prob > 0) %>%
  cor()

sc2_baal2 %>% distinct(label, .keep_all = T)  %>%  
  filter(ref.x != alt.x) %>% 
  select(probability, msa_aa_prob) %>% 
  drop_na() %>% 
  filter(msa_aa_prob > 0.005) %>%
  # ggplot() + geom_point(aes(probability, msa_aa_prob)) +
  # geom_smooth(method = 'lm', aes(probability, msa_aa_prob), color = 'red', se=F) +
  # labs(x='PLM amino-acid substitution likelihood',
  #      y = 'MSA amino-acid substitution probability',
  #      subtitle = 'SARS-CoV-2 Spike') +
  # theme_ipsum()
  cor()

hiv_baal3 %>% distinct(label, .keep_all = T)  %>% filter(ref.x.x != alt.x.x, SurfaceAccessibility_sequence > 0.5) %>% select(probability.x, aa_prob.x) %>% drop_na() %>% 
  filter(aa_prob.x > 0) %>% cor()


hiv_baal3 %>% distinct(label, .keep_all = T)  %>% 
  filter(ref.x.x != alt.x.x) %>% 
  filter(label %in% hiv_all_exp_mut$label) %>%
  select(probability.x, aa_prob.x) %>% drop_na() %>% 
  filter(aa_prob.x > 0) %>% 
  # ggplot() + geom_point(aes(probability.x, aa_prob.x)) +
  # geom_smooth(method = 'lm', aes(probability.x, aa_prob.x), color = 'red', se=F) +
  # labs(x='PLM amino-acid substitution likelihood',
  #      y = 'MSA amino-acid substitution probability',
  #      subtitle = 'HIV Env') +
  # theme_ipsum()
  cor()
  
flu_baal3 <- flu_baal2 %>% select(!msa_aa_prob) %>% left_join(flu_wsn33_msa_prob, by = 'label')

flu_baal3 %>% distinct(label, .keep_all = T)  %>% 
  filter(ref.x != alt.x) %>% 
  #filter(label %in% flu_wsn33_true_antibody_escape$label) %>% 
  select(probability, msa_aa_prob) %>% drop_na() %>% 
  filter(msa_aa_prob > 0) %>%
  # ggplot() + 
  # geom_point(aes(probability, msa_aa_prob)) +
  # geom_smooth(method = 'lm', aes(probability, msa_aa_prob), color = 'red', se=F) +
  # labs(x='PLM amino-acid substitution likelihood',
  #      y = 'MSA amino-acid substitution probability',
  #      subtitle = 'Influenza H1N1 HA') +
  # theme_ipsum()
  cor()


sc2_baal2 %>% distinct(label, .keep_all = T)  %>% filter(ref.x != alt.x, label %in% sc2_tre_antibody_escape$label) %>% 
  select(probability, msa_aa_prob) %>% drop_na() %>% 
  filter(msa_aa_prob <0.95) %>% 
  ggplot() + geom_point(aes(probability, msa_aa_prob)) + 
  labs(title = 'SARS-CoV-2 Spike protein', 
       subtitle = 'Correlation between PLM amino-acid substitutions likelihood \nand MSA amino-acid substitutions probability', 
       x = "PLM amino-acid substitutions likelihood", 
       y = 'MSA amino-acid substitutions probability') +
  theme_bw() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 24)
        )


flu_baal2 %>% distinct(label, .keep_all = T)  %>% filter(ref.x != alt.x) %>% 
  select(probability, msa_aa_prob) %>% drop_na() %>% 
  #filter(msa_aa_prob <0.95) %>% 
  ggplot() + geom_point(aes(probability, msa_aa_prob)) + 
  labs(title = 'Influenza H1N1 HA protein', 
       subtitle = 'Correlation between PLM amino-acid substitutions likelihood \nand MSA amino-acid substitutions probability', 
       x = "PLM amino-acid substitutions likelihood", 
       y = 'MSA amino-acid substitutions probability') +
  theme_bw() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 24)
  )

hiv_baal3 %>% distinct(label, .keep_all = T)  %>% filter(ref.x.x != alt.x.x) %>% 
  select(probability.x, aa_prob.x) %>% drop_na() %>% 
  #filter(msa_aa_prob <0.95) %>% 
  ggplot() + geom_point(aes(probability.x, aa_prob.x)) + 
  labs(title = 'HIV Envelope protein', 
       subtitle = 'Correlation between PLM amino-acid substitutions likelihood \nand MSA amino-acid substitutions probability', 
       x = "PLM amino-acid substitutions likelihood", 
       y = 'MSA amino-acid substitutions probability') +
  theme_bw() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 24)
  )



sc2_wuhan_msa_prob %>% filter(msa_aa_prob > 0.1, msa_aa_prob < 0.25) %>% ggplot(aes(x= seq_along(msa_aa_prob), y = msa_aa_prob)) + geom_point() + geom_label(aes(label = label))

sc2_wuhan_msa_prob %>% filter(msa_aa_prob > 0.1, msa_aa_prob < 0.25) %>% left_join(sc2_wuhan_plm_prob, by = 'label') %>% 
  ggplot(aes(x= msa_aa_prob, y = probability)) + geom_point() + geom_label(aes(label = label))
