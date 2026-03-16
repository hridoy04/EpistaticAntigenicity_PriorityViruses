#this code is for merging multiple dataframes of viruses into a single dataframe


# Load necessary libraries
library(tidyverse)


# Load the dataframes
#HIV


#dataframe for plm_logit probabilities and shannon entropy/perplexity 

hiv_plm_prob <- read_csv('hiv_env_glycoprotein_new.csv')
glimpse(hiv_plm_prob)
hiv_plm_prob %>% pull(probability)

hiv_plm_shannon_entropy <- read_csv('hiv_env_glycoprotein_new_entropy.csv')
glimpse(hiv_plm_shannon_entropy)
hiv_plm_shannon_entropy <- hiv_plm_shannon_entropy %>% 
                           rename(position = `...1`) %>% 
                           mutate(position = position + 1) %>% 
                           slice(rep(1:n(), each = 20))

hiv_plm_perplexity <- read_csv('hiv_env_glycoprotein_new_perplexity.csv')
glimpse(hiv_plm_perplexity)
hiv_plm_perplexity <- hiv_plm_perplexity %>% 
                       rename(position = `...1`) %>% 
                       mutate(position = position + 1) %>% 
                       slice(rep(1:n(), each = 20))

#hiv_plm_perplexity <- tibble(rep(hiv_plm_perplexity, each = 20))


#dataframe of surface accessibility from single sequence (single aa based)
hiv_surface_accessibility <- read_csv('hiv_env_acc.csv')
glimpse(hiv_surface_accessibility)
hiv_surface_accessibility <- hiv_surface_accessibility %>% 
  select(seq, n, rsa) %>% 
  rename(position = n, rel_surf_acc = rsa, ref = seq) %>% 
  slice(rep(1:n(), each = 20))

#dataframe for antigenicity from single sequence (single aa based)
hiv_antigenicity <- read_csv('hiv_env_disco3.csv')
glimpse(hiv_antigenicity)
hiv_antigenicity <- hiv_antigenicity %>% 
  select(res_id, residue, `DiscoTope-3.0_score`) %>% 
  rename(position = res_id, antigenicity = `DiscoTope-3.0_score`, ref = residue) %>% 
  slice(rep(1:n(), each = 20)) %>% 
  filter(row_number() <= n()-20)

#dataframe for epistatic antigenicity (multiple aa change effect on a single aa change)
hiv_epistatic_antigenicity <- read_csv('hiv_env_epistatic_accessibility.csv')
glimpse(hiv_epistatic_antigenicity)

hiv_singleseq_analyses_df <- bind_cols(hiv_plm_shannon_entropy, hiv_plm_perplexity, hiv_surface_accessibility, hiv_antigenicity, hiv_epistatic_antigenicity) %>% 
  select(!c(`position...3`, `position...6`, `position...8`, `ref...5`, `position...13`, `ref...14`, `...11`)) %>% 
  rename(position = `position...1`,  ref = `ref...9`) %>% 
  right_join(hiv_plm_prob, by = c('label', 'position', 'ref', 'alt'))

hiv_singleseq_analyses_df

#dataframe for sequence based msa probabilities
hiv_msa_prob <- read_csv('hiv_msa_aa_freq.csv')
glimpse(hiv_msa_prob)
hiv_msa_prob <- head(hiv_msa_prob, -20)
hiv_msa_prob <- hiv_msa_prob %>% 
  rename(label = aa_change) 

## Dataframes for experimental works ##
#dataframe for DMS experiments
hiv_dms_prob1 <- read_csv('DMS_Dingens2019a_hiv_env_antibodies_x10.csv')
glimpse(hiv_dms_prob1)
hiv_dms_score1 <- hiv_dms_prob1 %>% unite(label, c(wt,i,mut), sep = '') %>%
                  select(label, starts_with('summary_'))
glimpse(hiv_dms_score1)


hiv_dms_prob2 <- read_csv('DMS_Haddox2018_hiv_BG505_env_replication_pref.csv')
glimpse(hiv_dms_prob2)
hiv_dms_score2 <- hiv_dms_prob2 %>% select(!site) %>% rename(ref = wt, alt = mut, position = i, label = mutation)
glimpse(hiv_dms_score2)
length(hiv_dms_score1$label[5])

hiv_dms_score_all <- hiv_dms_score1 %>% left_join(hiv_dms_score2, by = 'label')
str(hiv_dms_score_all)
hiv_dms_score_all %>% arrange(position) %>% pull(label)

#evescape hiv env escape data
hiv_evescape_score <- read_csv('hiv_env_evescape.csv')
glimpse(hiv_evescape_score)
hiv_eve <- hiv_evescape_score %>% mutate(label = paste0(wt, as.character(i), mut))

#merge all dfs into one DF

HIV_FINAL_DF_v1 <- hiv_singleseq_analyses_df %>% 
  right_join(hiv_msa_prob, by = c('label')) %>% 
  left_join(hiv_dms_score_all, by = 'label')

write_csv(HIV_FINAL_DF_v1, 'HIV_FINAL_DF_v1.csv')

hiv_singleseq_analyses_df %>% 
  arrange(position) %>% 
  left_join(hiv_msa_prob %>% select(!`...1`), by = c('label')) %>% 
  left_join(hiv_dms_score_all, by = c('label'))
  # pull('summary_VRC34-medianmutfracsurvive')


merge(hiv_1, hiv_dms_score_all, by = c('position', 'label'))

hiv_2 <- read_csv('HIV_BG505_Env_new.csv') 
hiv_2 %>% left_join(hiv_dms_score_all, by = 'label') %>% filter(position.y >= 30) %>% tail() %>% select(starts_with('summary'))

#baal
baal <- read_csv('hiv_env_evescape.csv')
glimpse(baal)
baal <- baal %>% mutate(label = paste0(wt,i,mut)) %>% rename(position = i, ref = wt, alt = mut)

hiv_dm
hiv_exp <- merge(baal,hiv_dms_score_all, by = c('position', 'ref', 'alt', 'label')) %>% arrange(position) %>% tibble()
hiv_exp$label[1:20]

#FLU

#dataframe for plm_logit probabilities and shannon entropy/perplexity

flu_plm_prob <- read_csv('flu_h3n2_ha_new.csv')




#FLUU


flu_wsn33_plm_prob <- read_csv('Influenza_H1_WSN_1933_new.csv')
glimpse(flu_wsn33_plm_prob)

flu_wsn33_plm_epi_antigenicity <- read_csv('fluh1_wsn1933_ha_alpha_epistatic_accessibility.csv')
glimpse(flu_wsn33_plm_epi_antigenicity)

flu_wsn33_single_antigenicity <- read_csv('flu_wsn1933_h1_ha_alpha_discotope.csv')
flu_wsn33_single_antigenicity <- flu_wsn33_single_antigenicity %>% rename(ref = residue, position = res_id, single_aa_antigenicity = `DiscoTope-3.0_score`)
glimpse(flu_wsn33_single_antigenicity)

flu_wsn33_singleseq_score <- flu_wsn33_plm_prob %>% left_join(flu_wsn33_plm_epi_antigenicity, by = c('position', 'ref', 'alt', 'label')) %>% 
  left_join(flu_wsn33_single_antigenicity, by = c('position', 'ref'))
glimpse(flu_wsn33_singleseq_score)


flu_wsn33_msa_prob <- read_csv('flu_msa_aa_freq.csv')
flu_wsn33_msa_prob <- flu_wsn33_msa_prob %>% 
  rename(label = aa_change, position = `...1`, msa_aa_prob = aa_prob) %>% 
  select(!position)
  
glimpse(flu_wsn33_msa_prob)


flu_wsn33_exp <- read_csv('flu_h1_evescape.csv')
glimpse(flu_wsn33_exp)
flu_wsn33_exp <- flu_wsn33_exp %>% mutate(label = paste0(wt,i,mut)) %>% rename(position = i, ref = wt, alt = mut)


flu_wsn33_final <- flu_wsn33_singleseq_score %>% left_join(flu_wsn33_msa_prob, by = 'label') %>% 
  left_join(flu_wsn33_exp, by = c('label', 'position', 'ref', 'alt'))
glimpse(flu_wsn33_final)
write_csv(flu_wsn33_final, 'flu_wsn33_final.csv')

flu_wsn33_true_antibody_escape <- read_csv('h1_label_abescape_exp.csv')
glimpse(flu_wsn33_true_antibody_escape)



flu_baal <- read_csv('InfluenzaWSN33_H1HA_epistatic_accessibility_antigenicity_v2.csv')
fl_new
flu_baal <- fl_new %>% left_join(flu_baal, by = 'label')
nrow(flu_baal)


flu_acc <- read_csv('InfluenzaWSN33_H1HARSA.csv')
flu_ant <- read_csv('/Users/rubayetalam/InfluenzaWSN33_H1HAAntigenicity.csv')

flu_acc_ant <- cbind(flu_acc, flu_ant) %>% slice(rep(1:n(), each = 20))
nrow(flu_acc_ant)

flu_baal2 <- cbind(flu_baal, flu_acc_ant)
colnames(flu_baal2)


flu_baal2 <- flu_baal2 %>% select(!position) %>% rename(Epistatic_Antigenicity_structure = significant_accessibility,
                     Epistatic_Antigenicity_sequence = significant_antigenicity_median.y,
                     Epistatic_SurfaceAccessibility_sequence = significant_accessibility_median.y,
                     Antigenicity_structure = single_aa_antigenicity,
                     Antigenicity_sequence = Antigenicity,
                     SufaceAccessibility_sequence = SurfaceAccessibility,
                     'DMS score' = max_escape_experiment,
                     'EVEscape score' = evescape)

flu_baal2 <- flu_baal2 %>% rename(SurfaceAccessibility_sequence = SufaceAccessibility_sequence) 

write_csv(flu_baal2, 'flu_wsn33_epistatic_accessibility_final_v2.csv')
write_csv(flu_baal2, 'flu_wsn33_epistatic_antigenicity_final.csv')

Filter(is.numeric, flu_baal2 %>% 
         filter(label %in% flu_wsn33_true_antibody_escape$label)) %>% 
  select(Epistatic_Antigenicity_structure, 
         Epistatic_Antigenicity_sequence, 
         Epistatic_SurfaceAccessibility_sequence, 
         Antigenicity_structure,
         Antigenicity_sequence,
         SurfaceAccessibility_sequence,
         'DMS score', 'EVEscape score') %>% 
  drop_na() %>% 
  #filter(max_escape_experiment >= 0.05) %>% 
  as.matrix %>%
  cor %>% 
  corrplot(#order="alphabet", 
    #main = 'HIV Envelope protein',
    tl.col = 'black', 
    diag = FALSE,
    #type = 'upper',
    method = 'color',
    addCoef.col = 'black',
    cl.cex = 1,
    tl.cex = 1,
    pch.cex = 1.5,
    col=rev(brewer.pal(n=10, name="RdYlBu")))




#SARS-CoV-2

sc2_wuhan_plm_prob <- read_csv('sc2_spike_new.csv')
glimpse(sc2_wuhan_plm_prob)

sc2_wuhan_plm_epi_antigenicity <- read_csv('sc2_spike_epistatic_accessibility.csv')
glimpse(sc2_wuhan_plm_epi_antigenicity)

sc2_wuhan_single_antigenicity <- read_csv('sc2_bepro_antigenicity.csv')
glimpse(sc2_wuhan_single_antigenicity)
sc2_wuhan_single_antigenicity <- sc2_wuhan_single_antigenicity %>% rename(single_aa_antigenicity = Accessibility)


sc2_wuhan_singleseq_score <- sc2_wuhan_plm_prob %>% left_join(sc2_wuhan_plm_epi_antigenicity, by = c('position', 'ref', 'alt', 'label')) %>% 
  left_join(sc2_wuhan_single_antigenicity, by = c('position', 'ref'))
glimpse(sc2_wuhan_singleseq_score)

sc2_wuhan_msa_prob <- read_csv('sc2_msa_aa_freq.csv')
glimpse(sc2_wuhan_msa_prob)
sc2_wuhan_msa_prob <- sc2_wuhan_msa_prob %>% 
  rename(label = aa_change, position = `...1`, msa_aa_prob = aa_prob) %>% 
  select(!position)


sc2_wuhan_exp <- read_csv('spike_rbd_evescape.csv')
glimpse(sc2_wuhan_exp)
sc2_wuhan_exp <- sc2_wuhan_exp %>% mutate(label = paste0(wt,i,mut)) %>% rename(position = i, ref = wt, alt = mut)

sc2_wuhan_eve <- read_csv('full_spike_evescape.csv')
glimpse(sc2_wuhan_eve)
sc2_wuhan_eve <- sc2_wuhan_eve %>% mutate(label = paste0(wt,i,mut)) %>% rename(position = i, ref = wt, alt = mut)

sc2_fullspike_wuhan_final <- sc2_wuhan_singleseq_score %>% left_join(sc2_wuhan_msa_prob, by = 'label') %>% 
  left_join(sc2_wuhan_eve, by = c('label', 'position', 'ref', 'alt'))
glimpse(sc2_fullspike_wuhan_final)

sc2_rbd_wuhan_final <- sc2_wuhan_singleseq_score %>% left_join(sc2_wuhan_msa_prob, by = 'label') %>% 
  left_join(sc2_wuhan_exp, by = c('label', 'position', 'ref', 'alt'))
sc2_rbd_wuhan_final
#other bloom's dms score on other variants like xbb, ba, delta, kp:
#xbb.1.5: mean of escape score from 10 polyclonal Abs sera

sc2_xbb15_dms <- read_csv('/Users/rubayetalam/Downloads/summary.csv')
xbb_dms_final <- sc2_xbb15_dms %>% rename(ref = wildtype, alt = mutant, position = site, max_escape_experiment = `human sera escape`) %>% 
  unite(label, c(ref, position, alt), sep = '', remove = F) %>% 
  select(position, ref, alt, label, max_escape_experiment)
glimpse(xbb_dms_final)

#delta 
delta_files <- list.files('/Users/rubayetalam/Downloads/bloom_dms/delta')
delta_filename <- paste('/Users/rubayetalam/Downloads/bloom_dms/delta', delta_files, sep = '/')
delta_dms_combined_df <- bind_rows(lapply(delta_filename, read_csv))
glimpse(delta_dms_combined_df)

delta_dms_final <- delta_dms_combined_df %>% group_by(site, wildtype, mutant, mutation) %>% 
  summarise(max_escape_experiment = max(escape_mean, na.rm = TRUE)) %>% 
  rename(label = mutation, ref=wildtype, alt= mutant, position = site) %>% ungroup()
glimpse(delta_dms_final)

#BA.1
ba1_files <- list.files('/Users/rubayetalam/Downloads/bloom_dms/ba1')
ba1_filename <- paste('/Users/rubayetalam/Downloads/bloom_dms/ba1', ba1_files, sep = '/')
ba1_dms_combined_df <- bind_rows(lapply(ba1_filename, read_csv))
glimpse(ba1_dms_combined_df)
ba1_dms_final <- ba1_dms_combined_df %>% group_by(site, wildtype, mutant, mutation) %>% 
  summarise(max_escape_experiment = max(escape_mean, na.rm = TRUE)) %>% 
  rename(label = mutation, ref=wildtype, alt= mutant, position = site) %>% ungroup()
ba1_dms_final$position <- as.numeric(ba1_dms_final$position)
glimpse(ba1_dms_final)

#BA.2 not done by bloom
#KP.3 in preprint still

#wuhan rbd
wuhan_dms_rbd_final <- sc2_wuhan_exp %>% select(position, ref, alt, label, max_escape_experiment_bloom, max_escape_experiment_xie) %>% 
  mutate(max_escape_experiment = if_else(max_escape_experiment_bloom>=max_escape_experiment_xie, max_escape_experiment_bloom, max_escape_experiment_xie)) %>%
  select(!c(max_escape_experiment_bloom, max_escape_experiment_xie)) 
#wuhan_dms_rbd_final$position <- as.character(wuhan_dms_rbd_final$position)
glimpse(wuhan_dms_rbd_final)

#final merge of all dms data
SC2_DMS_all <- bind_rows(wuhan_dms_rbd_final, xbb_dms_final, delta_dms_final, ba1_dms_final)
glimpse(SC2_DMS_all %>% arrange(position) %>% filter(position == 3))

write_csv(SC2_DMS_all, 'sc2_all_dms_escape.csv')

sc2_dms_all_final <- SC2_DMS_all %>% drop_na() %>% group_by(position, alt) %>% 
  summarise(max_escape_experiment = max(max_escape_experiment, na.rm = TRUE)) %>% ungroup() %>% 
  arrange(position) 
glimpse(sc2_dms_all_final)
sc2_dms_all_final %>% filter(position == 484)

#matching with refseq label:
pos_diff <- sc2_dms_final %>% distinct(position, .keep_all = T) %>% right_join(sc2_refseq_spike, by = c('position', 'ref')) %>% 
  filter(is.na(max_escape_experiment)) %>% pull(position)

sc2_refseq_spike <- read_csv('sc2_wuhan_ref.csv')
sc2_refseq_spike
sc2_refseq_spike <- sc2_refseq_spike %>% rename(position = pos, ref = refseq_wuhan) %>% select(!`...1`) 

sc2_dms_final <- sc2_dms_all_final %>% left_join(sc2_refseq_spike, by = c('position')) %>% 
  unite(label, c(ref, position, alt), sep = '', remove = F) %>% 
  select(position, ref, alt, label, max_escape_experiment) %>% arrange(position) 
sc2_dms_final %>% filter(position == 484)

sc2_true_antibody_escape <- read_csv('~/SC2_antigenic/sc2_spike_abescape_mutations.csv')
glimpse(sc2_true_antibody_escape)


sc2_true_high_antibody_escape <- read_csv('sc2_spike_high_abescape_mutations.csv')
glimpse(sc2_true_high_antibody_escape)


write_csv(sc2_wuhan_final, 'sc2_wuhan_final.csv')

SC2_Wuhan_full <- sc2_fullspike_wuhan_final %>% filter(label %in% sc2_true_high_antibody_escape$`Mutation\n`)
SC2_Wuhan_rbd <- sc2_rbd_wuhan_final %>% filter(label %in% sc2_true_high_antibody_escape$`Mutation\n`)
SC2_Wuhan_full


#merging all into one df
sc2_fullspike_final <- sc2_wuhan_singleseq_score %>% left_join(sc2_wuhan_msa_prob, by = 'label') %>% 
  left_join(sc2_wuhan_eve, by = c('label', 'position', 'ref', 'alt')) %>% 
  left_join(sc2_dms_final, by = c('label', 'position', 'ref', 'alt'))
glimpse(sc2_fullspike_final)

sc2_high_abescape_final <- sc2_fullspike_final %>% filter(label %in% sc2_true_high_antibody_escape$`Mutation\n`)
sc2_abescape_final <- sc2_fullspike_final %>% filter(label %in% sc2_true_antibody_escape$`Mutation\n`)
sc2_fullspike_abescape_final$evescape


#adding new df for antigenicity
sc2_new_antigenicity <- read_csv('/Users/rubayetalam/three_viruses/sc2/sc2_spike_epistatic_accessibility_antigenicity.csv')
glimpse(sc2_new_antigenicity)

sc2_dms_final <- read_csv('~/SC2_antigenic/sc2_all_dms_escape.csv')


sc2_wuhan_eve <- read_csv('~/SC2_antigenic/full_spike_evescape.csv') %>% mutate(label = paste0(wt,i,mut)) %>% rename(position = i, ref = wt, alt = mut)

sc2_old_antigenicity <- read_csv('~/SC2_antigenic/sc2_wuhan_final.csv')

sc2_fullspike_final_v2 <- sc2_old_antigenicity %>% 
  left_join(sc2_wuhan_eve, by = c('label', 'position', 'ref', 'alt')) %>% 
  left_join(sc2_dms_final, by = c('label', 'position', 'ref', 'alt')) %>%
  left_join(sc2_new_antigenicity, by = c('label', 'position', 'ref', 'alt'))
glimpse(sc2_fullspike_final_v2)
head(sc2_fullspike_final_v2, 20)


sc2_full_spike_final_v3 <- sc2_fullspike_final_v2 %>% distinct(label, .keep_all = T)

sc2_fullspike_final_v4 <- sc2_old_antigenicity %>% 
  left_join(sc2_wuhan_eve, by = c('label', 'position', 'ref', 'alt')) %>% 
  inner_join(sc2_dms_final %>% drop_na(max_escape_experiment), by = c('label', 'position', 'ref', 'alt')) %>% 
  group_by(label) %>% mutate(max_escape_experiment = max(max_escape_experiment)) %>% ungroup() %>% 
  distinct(label, .keep_all = T) %>% 
  left_join(sc2_new_antigenicity, by = c('label', 'position', 'ref', 'alt')) 
head(sc2_fullspike_final_v4)

glimpse(sc2_fullspike_final_v4)
sc2_fullspike_final_v4 %>% filter(is.na(evescape.x)) %>% pull(label)


sc2_tre_high_antibody_escape <- read_csv('/Users/rubayetalam/SC2_antigenic/sc2_spike_high_abescape_mutations.csv')
sc2_tre_high_antibody_escape <- sc2_tre_high_antibody_escape %>% select(`Mutation\n`) %>% rename(label = `Mutation\n`)
sc2_tre_high_antibody_escape


sc2_tre_antibody_escape <- sc2_true_antibody_escape %>% select(`Mutation\n`) %>% rename(label = `Mutation\n`)


sc2_baal <- read_csv('sc2_spike_epistatic_accessibility_antigenicity_v2.csv')
colnames(sc2_baal)
sc2_acc <- read_csv('sc2_spikeRSA.csv')
sc2_ant <- read_csv('sc2_spikeAntigenicity.csv')

sc2_acc

sc2_acc_ant <- cbind(sc2_acc, sc2_ant) %>% slice(rep(1:n(), each = 20)) %>% select(!c(position, ref))

sc2_baal2 <- sc2_fullspike_final_v2 %>% left_join(cbind(sc2_baal, sc2_acc_ant), by = 'label')


sc2_baal2 <- sc2_baal2 %>% rename('Epistatic_Antigenicity_structure' = significant_accessibility,
                     'Epistatic_Antigenicity_sequence' = significant_antigenicity_median.y,
                     'Epistatic_SurfaceAccessibility_sequence' = significant_accessibility_median.y,
                     'Antigenicity_structure' = single_aa_antigenicity,
                     'Antigenicity_sequence' = Antigenicity,
                     'SurfaceAccessibility_sequence' = SurfaceAccessibility,
                     'DMS score' = max_escape_experiment,
                     'EVEscape score' = evescape.x)

sc2_baal

sc2_baal4 <- sc2_fullspike_final_v4 %>% left_join(cbind(sc2_baal, sc2_acc_ant), by = 'label')


sc2_baal4 <- sc2_baal4 %>% rename('Epistatic_Antigenicity_structure' = significant_accessibility,
                                  'Epistatic_Antigenicity_sequence' = significant_antigenicity_median.y,
                                  'Epistatic_SurfaceAccessibility_sequence' = significant_accessibility_median.y,
                                  'Antigenicity_structure' = single_aa_antigenicity,
                                  'Antigenicity_sequence' = Antigenicity,
                                  'SurfaceAccessibility_sequence' = SurfaceAccessibility,
                                  'DMS score' = max_escape_experiment,
                                  'EVEscape score' = evescape.x)

Filter(is.numeric, sc2_baal4 %>% filter(label %in% sc2_tre_high_antibody_escape$label)) %>% 
  select(Epistatic_Antigenicity_structure, 
         Epistatic_Antigenicity_sequence, 
         Epistatic_SurfaceAccessibility_sequence, 
         Antigenicity_structure,
         Antigenicity_sequence,
         SurfaceAccessibility_sequence,
         `DMS score`, `EVEscape score`) %>% 
  drop_na() %>% 
  #filter(max_escape_experiment >= 0.05) %>% 
  as.matrix %>%
  cor %>% 
  corrplot(#order="alphabet", 
    #main = 'SARS-CoV-2 Spike protein',
    tl.col = 'black', 
    diag = FALSE,
    #type = 'upper',
    method = 'color',
    addCoef.col = 'black',
    cl.cex = 1,
    tl.cex = 1,
    pch.cex = 1.5,
    col=rev(brewer.pal(n=10, name="RdYlBu")))

colnames(sc2_baal4)

write_csv(sc2_baal2, 'sc2_spike_epistatic_accessibility_final_v2.csv')
write_csv(sc2_baal4, 'sc2_spike_epistatic_antigenicity_final.csv')

##HIV
#comparison to real data
#HIV bg505 env

HIV_BG505Env_EpiAcc <- read_csv('hivbg505_env_alpha_epistatic_accessibility.csv')
HIV_BG505Env_EpiAcc <- HIV_BG505Env_EpiAcc %>% rename(Epistatic_AntigenicitySignificant = significant_accessibility, 
                                                      Epistatic_AntigenicityVSignificant = very_significant_accessibility)

HIV_PLMProb <- read_csv('HIV_BG505_Env_new.csv')
HIV_PLMProb


HIV_MSAProb <- read_csv('hiv_bg505_msa_aa_freq.csv')
HIV_MSAProb <- HIV_MSAProb %>% rename(label = aa_change)
HIV_MSAProb$label[1:20]
HIV_MSAProb


HIV_PLMEnt <- read_csv('HIV_BG505_Env_new_entropy.csv')
glimpse(HIV_PLMEnt)
HIV_PLMEnt <- HIV_PLMEnt %>% 
  rename(position = `...1`) %>% 
  mutate(position = position + 1)

HIV_PLMPerp <- read_csv('HIV_BG505_Env_new_perplexity.csv')
glimpse(HIV_PLMPerp)
HIV_PLMPerp <- HIV_PLMPerp %>% 
  rename(position = `...1`) %>% 
  mutate(position = position + 1) 

HIV_Antigenicity <- read_csv('hiv_bg505_env_alpha_discotope.csv')
glimpse(HIV_Antigenicity)
HIV_Antigenicity <- HIV_Antigenicity %>% 
  rename(position = res_id, SingleAA_Antigenicity = `DiscoTope-3.0_score`, ref = residue) 

hiv_bg505_pos_value <- list(HIV_PLMEnt, HIV_PLMPerp, HIV_Antigenicity) %>% 
  reduce(right_join, by = 'position') %>% 
  slice(rep(1:n(), each = 19))
hiv_pos <- hiv_bg505_pos_value %>% 
  filter(position >=30 & position <= 699) %>% 
  select(!c(position, ref))
hiv_pos
hiv_ultimate <- NULL
hiv_ultimate <- right_join(HIV_PLMProb, HIV_MSAProb, by = c('label')) %>% right_join(HIV_BG505Env_EpiAcc, by = c('label', 'position', 'ref', 'alt')) %>% 
  right_join(hiv_exp, by = c('label', 'position', 'ref', 'alt')) %>% bind_cols(hiv_ultimate, hiv_pos)
glimpse(hiv_ultimate)

hiv_weight <- read_csv('hiv_bg505_env_epistatic_accessibility_weighted.csv')


write_csv(hiv_ultimate, 'hiv_bg505_env_alpha_ultimate.csv')

hiv_ultimate %>% ggplot() + geom_point(aes(position,probability))
hiv_ultimate %>% select(position, ref, alt, label) %>% distinct(position, .keep_all = T) %>% print(n=100)

hiv_premium <- hiv_ultimate %>% left_join(hiv_weight %>% select(label, sig_like_weighted), by = 'label') 
nrow(hiv_premium)


hivbg505_df <- read_csv("HIV_FINAL_DF_v1.csv")
hivbg <- read_csv('~/three_viruses/hiv/HIVbg505_Env_epistatic_accessibility_antigenicity.csv')

hivbg505_df_final <- hiv_ultimate %>% left_join(hivbg, by = 'label')


hivbg505_df_final <-  hivbg505_df_final %>% 
  rename('EVEscape score' = 'evescape',
         'DMS score' = 'max_escape_experiment',
         'EpistaticAntigenicity(structure)' = 'Epistatic_AntigenicitySignificant',
         'EpistaticAntigenicity(sequence)' = 'significant_antigenicity_median',
         'LinkedRSA(sequence)' = 'significant_accessibility_median', 
         'Antigenicity(structure): individual amino acid' = 'SingleAA_Antigenicity') 
glimpse(hivbg505_df_final)

hivcc <- hiv_corrraa %>% left_join(hivbg, by = 'label')
colnames(hivcc) 


#adding latest epistatic antigenicity data
hivbg505_df_v2 <- read_csv('/Users/rubayetalam/HIVbg505_Env_epistatic_accessibility_antigenicity_v2.csv')
colnames(hivbg505_df_v2)

hiv_baal <- hiv_premium %>% left_join(hivbg505_df_v2, by = 'label') %>% left_join(hivbg505_df_final, by = 'label')
colnames(hiv_baal)

hiv_baal2 <- hiv_baal %>% rename(Epistatic_Antigenicity_structure = very_significant_accessibility_sum.x,
                    Epistatic_Antigenicity_sequence = significant_antigenicity_sum,
                    Epistatic_SurfaceAccessibility_sequence = very_significant_accessibility_sum.y
                    )
nrow(hiv_baal2)

hiv_seq_anti <- read_csv('/Users/rubayetalam/HIVbg505_EnvAntigenicity.csv')
hiv_seq_anti

hiv_seq_acc <- read_csv('/Users/rubayetalam/HIVbg505_EnvRSA.csv')
hiv_seq_acc

hiv_seq_acc_anti <- cbind(hiv_seq_acc, hiv_seq_anti)
hiv_seq_acc_anti

hiv_seq_acc_anti_final <- hiv_seq_acc_anti %>% slice(rep(1:n(), each = 19))
glimpse(hiv_seq_acc_anti_final)
hiv_accant_final <- hiv_seq_acc_anti_final %>% filter(position >= 30 & position <700)
glimpse(hiv_accant_final)

summary(hiv_baal2$position.x.x)

hiv_baal3 <- cbind(hiv_baal2, hiv_accant_final) 
hiv_baal3 <- hiv_baal3 %>% rename('Antigenicity_structure' = SingleAA_Antigenicity, 
                                  'Antigenicity_sequence' = Antigenicity,
                    'SurfaceAccessibility_sequence' = SurfaceAccessibility)
  
write_csv(hiv_baal3, 'HIV_BG505_Env_Antigenicity_final.csv')

  
glimpse(hiv_baal3)

Filter(is.numeric, hiv_baal2 %>% filter(label %in% hiv_corrraa$label)) %>% 
  select(!starts_with('summary_')) %>% 
  select(!c(7, 13, 15, 18, 24)) %>% 
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

Filter(is.numeric, hiv_baal3 %>% filter(label %in% hiv_corrraa$label)) %>% 
  select(Epistatic_Antigenicity_structure, 
         Epistatic_Antigenicity_sequence, 
         Epistatic_SurfaceAccessibility_sequence, 
         Antigenicity_structure,
         Antigenicity_sequence,
         SurfaceAccessibility_sequence,
         'DMS score', 'EVEscape score') %>% 
  drop_na() %>% 
  #filter(max_escape_experiment >= 0.05) %>% 
  as.matrix %>%
  cor %>% 
  corrplot(#order="alphabet", 
    main = 'HIV Envelope protein',
    tl.col = 'black', 
    diag = FALSE,
    #type = 'upper',
    method = 'color',
    addCoef.col = 'black',
    cl.cex = 1,
    tl.cex = 1.5,
    pch.cex = 1.5,
    col=rev(brewer.pal(n=10, name="RdYlBu")))


#new antigenicity trying 

sc2_epitope <- read_csv('/Users/rubayetalam/mutation_epitope_scores.csv')
sc2_epitope <- sc2_epitope %>% mutate(Position = Position + 327, 
                                      label = paste0(Reference, as.character(Position), Mutated)) %>% 
                                      select(label, Epitope_Score)

sc2_rbd_baal1 <- sc2_baal2 %>% right_join(sc2_epitope, by = 'label') 
sc2_rbd_baal2 <- sc2_rbd_baal1 %>% mutate(Antigenicity_sequence = 1/(1+exp(-(Antigenicity_sequence - mean(Antigenicity_sequence))/sd(Antigenicity_sequence))),
                         Epitope_Score = 1/(1+exp(-(Epitope_Score - mean(Epitope_Score))/sd(Epitope_Score))),
                         antibody_escape_score1 = (Antigenicity_sequence - Epitope_Score), 
                         antibody_escape_score2 = log(Epitope_Score/Antigenicity_sequence), 
                         antibody_escape_score3 = 1 - (Epitope_Score/Antigenicity_sequence) 
                         )

                        

