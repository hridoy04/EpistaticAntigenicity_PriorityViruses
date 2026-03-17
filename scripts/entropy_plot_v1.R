## this code assigns domains onto the lineplot to show entropy and other metrics

library(ggplot2)
library(dplyr)

# Define the bar segments with explicit start and end positions
bar_data <- data.frame(
  xmin = c(1, 4, 7),   # Start positions
  xmax = c(4, 7, 10),   # End positions
  ymin = 0,            # Bottom of the bar (aligned)
  ymax = 1,            # Height of the bar (fixed height)
  domain = factor(c('s1', 's2', 's1'), levels = c('s2', 's1')) # Order
)

bar_data
ggplot() + 
  # Large single bar split into s1 → s2 → s1 segments
  geom_rect(
    data = bar_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = domain),
    color = "black" # Optional border for better visibility
  ) +
  theme_void()

ggplot() + 
  geom_rect(
    data = marg_dom,
    aes(xmin = min(pos_start), xmax = max(pos_end), ymin = ymin, ymax = ymax),
    fill = 'grey99',
    color = 'black' # Optional border for better visibility
  ) +
  # Large single bar split into s1 → s2 → s1 segments
  geom_rect(
    data = marg_dom,
    aes(xmin = pos_start, xmax = pos_end, ymin = ymin, ymax = ymax, fill = domain),
    alpha = 0.15,
    color = "black" # Optional border for better visibility
  ) +
  geom_line(data = marburg, aes(x = position, y = entropy), size = 1) +
  theme_void()

marg_dom
marburg <- entropy_virus %>% filter(virus == 'marburg')


t(mapply(function(x,y) marburg %>% filter(position >=x & position < y), marg_dom$pos_start, marg_dom$pos_end)) %>% as_tibble()
t(mapply(function(x,y) marburg %>% filter(position >=x & position < y), marg_dom$pos_start, marg_dom$pos_end)) %>% as_tibble() %>% unlist()

t(mapply(function(x,y,z) marburg %>% mutate(domain = case_when(position >=x & position < y ~ z, TRUE ~ NA)) %>% filter(!is.na(domain)),  
         marg_dom$pos_start, marg_dom$pos_end, marg_dom$domain)) %>% as_tibble() %>% unnest()


mapply(function(x,y,z) marburg %>% filter(position >=x & position < y) %>% mutate(domain = z), 
         marg_dom$pos_start, marg_dom$pos_end, marg_dom$domain)



# Correct approach using mapply inside mutate
marburg_with_domains <- marburg %>%
  mutate(
    domain = reduce(
      mapply(
        function(x, y, z) {
          if_else(position >= x & position <= y, z, NA_character_)
        },
        marg_dom$pos_start,
        marg_dom$pos_end,
        marg_dom$domain,
        SIMPLIFY = FALSE
      ),
      coalesce
    ) %>% 
      replace_na("uncharacterized")
  )
marburg_with_domains
#make a function in R with name 'domain_assignment for this domain assignment
domain_assignment <- function(original_df, domain_df) {
  domain_assigned_df <- 
    original_df %>%
    mutate(
      domain = reduce(
        mapply(
          function(x, y, z) {
            if_else(position >= x & position <= y, paste(z,x,y, sep = '_'), NA_character_)
          },
          domain_df$pos_start,
          domain_df$pos_end,
          domain_df$domain,
          SIMPLIFY = FALSE
        ),
        coalesce
      ) %>% 
        replace_na("uncharacterized")
    )
  return(domain_assigned_df)
}


hiv <- entropy_virus %>% filter(virus == 'hiv')
sars <- entropy_virus %>% filter(virus == 'SARSCOV2')
influenza <- entropy_virus %>% filter(virus == 'influenza')
marburg <- entropy_virus %>% filter(virus == 'marburg')
borna <- entropy_virus %>% filter(virus == 'bornavirus')
nairo <- entropy_virus %>% filter(virus == 'nairo')
#do that for new entropy_virus dataset 



entropy_virus <- entropy_virus %>% 
  mutate(virus = case_when(virus == 'Sars-CoV-2' ~ 'SARS-CoV-2', TRUE ~ virus))

hiv <- entropy_virus %>% filter(virus == 'HIV')
sars <- entropy_virus %>% filter(virus == 'SARS-CoV-2')
influenza <- entropy_virus %>% filter(virus == 'Influenza')
marburg <- entropy_virus %>% filter(virus == 'Marburg')
borna <- entropy_virus %>% filter(virus == 'Bornavirus')
nairo <- entropy_virus %>% filter(virus == 'Orthonairovirus')

#for the entropy_virus_new, change the name with mutate to the capital letter based name and do filtering for each virus

entropy_virus_new <- entropy_virus_new %>% 
  mutate(virus = case_when(virus == 'SARSCOV2' ~ 'SARS-CoV-2', 
                           virus == 'hiv' ~ 'HIV', 
                           virus == 'influenza' ~ 'Influenza', 
                           virus == 'marburg' ~ 'Marburg', 
                           virus == 'bornavirus' ~ 'Bornavirus', 
                           virus == 'nairo' ~ 'Orthonairovirus', 
                           TRUE ~ virus))

hiv <- entropy_virus_new %>% filter(virus == 'HIV')
sars <- entropy_virus_new %>% filter(virus == 'SARS-CoV-2')
influenza <- entropy_virus_new %>% filter(virus == 'Influenza')
marburg <- entropy_virus_new %>% filter(virus == 'Marburg')
borna <- entropy_virus_new %>% filter(virus == 'Bornavirus')
nairo <- entropy_virus_new %>% filter(virus == 'Orthonairovirus')




sars
sars_domain <- read_csv('sars_domain.csv')
sars_all_domain <- domain_assignment(sars, sars_domain)
sd <- sars_all_domain %>% distinct(domain, .keep_all=T) %>% 
  separate(domain, '_', into = c('domain', 'pos_start', 'pos_end'), convert = TRUE) %>% 
  filter(domain != 'uncharacterized')
#same for hiv
hiv_domain <- read_csv('hiv_domain.csv')
hiv_all_domain <- domain_assignment(hiv, hiv_domain)
hd <- hiv_all_domain %>% distinct(domain, .keep_all=T) %>% 
  separate(domain, '_', into = c('domain', 'pos_start', 'pos_end'), convert = TRUE) %>% 
  filter(domain != 'uncharacterized')

#same for influenza
influenza_domain <- read_csv('influenza_domain.csv')
influenza_all_domain <- domain_assignment(influenza, influenza_domain)
id <- influenza_all_domain %>% distinct(domain, .keep_all=T) %>% 
  separate(domain, '_', into = c('domain', 'pos_start', 'pos_end'), convert = TRUE) %>% 
  filter(domain != 'uncharacterized')

#same for marburg
marburg_domain <- read_csv('marg.csv')
marburg_all_domain <- domain_assignment(marburg, marburg_domain)
md <- marburg_all_domain %>% distinct(domain, .keep_all=T) %>% 
  separate(domain, '_', into = c('domain', 'pos_start', 'pos_end'), convert = TRUE) %>% 
  filter(domain != 'uncharacterized')

#same for bornavirus
borna_domain <- read_csv('borna_domain.csv')
bornavirus_all_domain <- domain_assignment(borna, borna_domain)
 bd <- bornavirus_all_domain %>% distinct(domain, .keep_all=T) %>% 
  separate(domain, '_', into = c('domain', 'pos_start', 'pos_end'), convert = TRUE) %>% 
  filter(domain != 'uncharacterized')

#same for nairovirus
nairo_domain <- read_csv('nairo_gpc_domain.csv')
nairovirus_all_domain <- domain_assignment(nairo, nairo_domain)
nd <- nairovirus_all_domain %>% distinct(domain, .keep_all=T) %>% 
  separate(domain, '_', into = c('domain', 'pos_start', 'pos_end'), convert = TRUE) %>% 
  filter(domain != 'uncharacterized')

marburg_with_domains %>% print(n=500)
mapply(function(x,y) marburg %>% mutate() %>% pull(entropy), marg_dom$pos_start, marg_dom$pos_end)

ggplot() + 
  geom_rect(
    data = id,
    aes(xmin = min(pos_start), xmax = max(pos_end), ymin = 1, ymax = 1.5),
    fill = 'grey99',
    alpha = 1,
    color = 'black' # Optional border for better visibility
  ) +
  # Large single bar split into s1 → s2 → s1 segments
  geom_rect(
    data = id,
    aes(xmin = pos_start, xmax = pos_end, ymin = 1, ymax = 1.5, 
        fill = domain),
    position = 'identity',
    
    alpha = 0.15,
    show.legend = TRUE,
    color = "black" # Optional border for better visibility
  ) +
  geom_line(data = influenza_all_domain, aes(x = position, y = entropy), size = 1) +
  theme_void()

high_seq_virus <- bind_rows(influenza_all_domain, hiv_all_domain, sars_all_domain)
high_dom <- high_seq_virus %>% filter(domain != 'uncharacterized') %>% 
  distinct(virus, domain, .keep_all=T) %>% 
  separate(domain, '_', into = c('domain', 'pos_start', 'pos_end'), convert = TRUE)

unique(six_virus$virus)
six_virus <- bind_rows(influenza_all_domain, hiv_all_domain, sars_all_domain, 
                       marburg_all_domain, bornavirus_all_domain, nairovirus_all_domain)
six_virus_dom <- six_virus %>% filter(domain != 'uncharacterized') %>% 
  distinct(virus, domain, .keep_all=T) %>% 
  separate(domain, '_', into = c('domain', 'pos_start', 'pos_end'), convert = TRUE)
six_virus
six_virus_dom
 
ggplot()  + 
  # Large single bar split into s1 → s2 → s1 segments
  geom_rect(
    data = high_dom,
    aes(xmin = pos_start, xmax = pos_end, ymin = 1, ymax = 1.5, 
        fill = domain),
    alpha = 0.55,
    show.legend = FALSE,
    color = "black" # Optional border for better visibility
  ) +
  geom_text(data = high_dom, aes(x=pos_start+(pos_end-pos_start)/2, y = 1.25, label = domain), size = 3) +
  geom_line(data = high_seq_virus, aes(x = position, y = entropy), size = 1) +
  facet_wrap(~ virus, ncol=1) +
  theme_void()


ggplot()  + 
  # Large single bar split into s1 → s2 → s1 segments
  geom_rect(
    data = sd,
    aes(xmin = pos_start, xmax = pos_end, ymin = 1, ymax = 1.5, 
        fill = domain),
    alpha = 0.55,
    show.legend = FALSE,
    color = "black" # Optional border for better visibility
  ) +
  geom_text(data = sd, aes(x=pos_start+(pos_end-pos_start)/2, y = 1.25, label = domain), size = 3) +
  geom_line(data = sars, aes(x = position, y = entropy), size = 1) +
  facet_wrap(~ virus, ncol=1) +
  theme_void()
six_virus_dom
six_virus

ggplot()  + 
  # Large single bar split into s1 → s2 → s1 segments
  geom_rect(
    data = six_virus_dom,
    aes(xmin = pos_start, xmax = pos_end, ymin = 1, ymax = 1.25, 
        fill = domain),
    alpha = 0.55,
    show.legend = FALSE,
    color = "black" # Optional border for better visibility
  ) +
  geom_text(data = six_virus_dom, aes(x=pos_start+(pos_end-pos_start)/2, y = 1.125, label = domain), size = 3) +
  geom_line(data = six_virus, aes(x = position, y = entropy), size = 1) +
  scale_x_continuous(breaks = seq(0, 1250, by = 50)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  facet_wrap(~ factor(virus, levels = c('Bornavirus', 'Orthonairovirus', 'Marburg', 
                                         'Influenza', 'HIV', 'SARS-CoV-2')), 
             strip.position = 'left', 
             ncol=1) +
  labs(x = 'Position', y = 'Entropy') +
  theme_bw() +
  theme(strip.text = element_text(size = 10),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12))

write_csv(six_virus, 'six_virus_entropy.csv')
write_csv(six_virus_dom, 'domains_six_virus.csv')

