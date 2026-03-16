#for maximum epistatic antigenicity how surface accessibility works

#bornavirus
borna_epiant_surf %>% 
  ggplot() + 
  geom_boxplot(aes(x=surf, y=epiant, color = surf)) +
  labs(title='Borna Virus Epistatic Antigenicity grouped by surface accessibility',
       x='Surface Accessibility Groups',
       y='Epistatic Antigenicity (maxium at each site)') +
  theme_minimal() +
  theme(plot.title = element_text(size=14, face='bold'),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=12, face='bold'),
        axis.title.y = element_text(size=12, face='bold'))

wilcox.test(data = borna_epiant_surf, epiant ~ surf)


#cchfv

cchfv_sasa <- read_csv('cchfv_gpc_sasa.csv')
max(cchfv_sasa$sasa)
cchfv_rsa <- cchfv_sasa %>%
  mutate(
    rsa = (sasa * 100 / max(sasa)),
    surf = case_when(
      rsa >= 25 ~ 'exposed',
      rsa < 25 ~ 'buried'
    )
  ) %>%
  select(pos, sasa, rsa, surf) %>%
  rename(position = pos)

cchfv_rsa

cchfv_epiant <- read_csv('CCHFV_GPC_epistatic_accessibility_antigenicity.csv')
cchfv_epiant <- cchfv_epiant %>% group_by(position) %>% 
  summarise(epi = max(significant_antigenicity_median, na.rm=TRUE)) %>% 
  mutate(position = position + 993)

cchfv_epiant_max <- cchfv_epiant %>% group_by(position) %>% mutate(max_epi = max(significant_antigenicity_median, na.rm=TRUE)) %>% 
  select(position, label, ref, alt, significant_antigenicity_median) %>% 
  slice_max(significant_antigenicity_median) %>% 
  mutate(position = position + 993) %>% 
  ungroup()

tail(cchfv_epiant)

cchfv_epiant_surf <- cchfv_epiant %>%
  right_join(cchfv_rsa, by = 'position') %>%
  select(position, epi, rsa, surf)

cchfv_epiant_surf
cchfv_epiant_surf %>% ggplot() + geom_boxplot(aes(x=surf, y=epi, color = surf))  +
  labs(title='CCHFV GPC Epistatic Antigenicity grouped by surface accessibility',
       x='Surface Accessibility Groups',
       y='Epistatic Antigenicity (maxium at each site)') +
  theme_minimal() +
  theme(plot.title = element_text(size=14, face='bold'),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=12, face='bold'),
        axis.title.y = element_text(size=12, face='bold'))

wilcox.test(data = cchfv_epiant_surf, epi ~ surf)

#marburg

marburg_sasa <- read_csv('marburg_sasa.csv')
marburg_sasa
marburg_rsa <- marburg_sasa %>%
  mutate(
    rsa = (sasa * 100 / max(sasa)),
    surf = case_when(
      rsa >= 15 ~ 'exposed',
      rsa < 15 ~ 'buried'
    )
  ) %>%
  select(position, sasa, rsa, surf)
marburg_rsa

marburg_epiant <- read_csv('Marburg_GP1GP2_epistatic_accessibility_antigenicity.csv')
marburg_epiant_max <- marburg_epiant %>% 
  group_by(position) %>% 
  mutate(max_epi = max(significant_antigenicity_median, na.rm=TRUE)) %>% 
  select(position, label, ref, alt, significant_antigenicity_median) %>% 
  slice_max(significant_antigenicity_median) %>% 
  ungroup() 

  # filter(position %in% marv_mutations_list_all) %>% 
  # print(n=50)
marburg_epiant <- marburg_epiant %>% group_by(position) %>% 
  summarise(epi = max(significant_antigenicity_median, na.rm=TRUE))


marburg_epiant_surf <- marburg_epiant %>%
  right_join(marburg_rsa, by = 'position') %>%
  select(position, epi, rsa, surf)


marburg_epiant_surf %>% ggplot() + geom_boxplot(aes(x=surf, y=epi, color = surf)) +
  labs(title='Marburg Virus Epistatic Antigenicity grouped by surface accessibility',
       x='Surface Accessibility Groups',
       y='Epistatic Antigenicity (maxium at each site)') +
  theme_minimal() +
  theme(plot.title = element_text(size=14, face='bold'),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=12, face='bold'),
        axis.title.y = element_text(size=12, face='bold'))

wilcox.test(data = marburg_epiant_surf, epi ~ surf)

## for highest aa change likelihood the epistatic accessibility works

#bornavirus
borna_highlikelihood <- read_csv('epiant/borna.csv')
borna_highlikelihood <- borna_highlikelihood %>% 
  select(position, label, probability, significant_antigenicity_median) %>% 
  rename(epi = significant_antigenicity_median, pos = position)
borna_highlikelihood
borna_rsa

borna_highlikelihood_surf <- borna_highlikelihood %>%
  right_join(borna_rsa, by = 'pos') %>%
  select(pos, probability, epi, rsa, surf)
borna_highlikelihood_surf %>% ggplot() + geom_boxplot(aes(x=surf, y=epi, color = surf)) +
  labs(title='Bornavirus-1 p57', 
       subtitle = 'Epistatic Antigenicity grouped by surface accessibility for highest AA change likelihood',
       x='Surface Accessibility Groups',
       y='Epistatic Antigenicity') +
  theme_minimal() +
  theme(plot.title = element_text(size=14, face='bold'),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=12, face='bold'),
        axis.title.y = element_text(size=12, face='bold'))

wilcox.test(data = borna_highlikelihood_surf, epi ~ surf)

#cchfv
cchfv_highlikelihood <- read_csv('epiant/cchfv.csv')
cchfv_highlikelihood
cchfv_highlikelihood <- cchfv_highlikelihood %>% select(position, label, probability, significant_antigenicity_median) %>% 
  rename(epi = significant_antigenicity_median, pos = position)
cchfv_highlikelihood
cchfv_highlikelihood_surf <- cchfv_highlikelihood %>%
  right_join(cchfv_rsa, by = 'pos') %>%
  select(pos, probability, epi, rsa, surf)
cchfv_highlikelihood_surf %>% ggplot() + geom_boxplot(aes(x=surf, y=epi, color = surf)) +
  labs(title='CCHFV GPc', 
       subtitle = 'Epistatic Antigenicity grouped by surface accessibility for highest AA change likelihood at each site',
       x='Surface Accessibility Groups',
       y='Epistatic Antigenicity') +
  theme_minimal() +
  theme(plot.title = element_text(size=14, face='bold'),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=12, face='bold'),
        axis.title.y = element_text(size=12, face='bold'))

wilcox.test(data = cchfv_highlikelihood_surf, epi ~ surf)

#marburg

marburg_highlikelihood <- read_csv('epiant/mar.csv')
marburg_highlikelihood <- marburg_highlikelihood %>% select(position, label, probability, significant_antigenicity_median) %>% 
  rename(epi = significant_antigenicity_median, pos = position)
marburg_highlikelihood

marburg_highlikelihood_surf <- marburg_highlikelihood %>%
  right_join(marburg_rsa %>% rename(pos = position), by = 'pos') %>%
  select(pos, probability, epi, rsa, surf)
marburg_highlikelihood_surf$epi <- marburg_highlikelihood_surf$epi %>% replace_na(0)
glimpse(marburg_highlikelihood_surf)

marburg_highlikelihood_surf %>% ggplot() + 
  geom_boxplot(aes(x=surf, y=epi, color = surf)) +
  labs(title='Marburg Virus GP1GP2', 
       subtitle = 'Epistatic Antigenicity grouped by surface accessibility for highest AA change likelihood at each site',
       x='Surface Accessibility Groups',
       y='Epistatic Antigenicity') +
  theme_minimal() +
  theme(plot.title = element_text(size=14, face='bold'),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=12, face='bold'),
        axis.title.y = element_text(size=12, face='bold'))

wilcox.test(data = marburg_highlikelihood_surf, epi ~ surf)



marburg_highlikelihood_surf %>% ggplot(aes(epi, rsa)) + 
  geom_point(aes(color = surf))

borna_highlikelihood_surf %>% ggplot(aes(epi, rsa)) + 
  geom_point(aes(color = surf))


#CHECK the heatmap for the positions that are associated with neutralization resistance in literature

marburg_mut <- read_csv('marburg_mutations.csv')
marburg_mut$Substitution
# List of individual mutations and deletions found in the text
marv_mutations <- c(
  "S220P", "P455L", "C226Y", "D449A", 
  "N453A", "D457A", "D459A", "P446A", "F447A", 
  "I452A", "W439A", "F445A", 
  "K465E", "N297Q", "F72A", "W70A", "H124S", 
  "W598A"
)

marv_mutations_list <- c(
  "S220P", "P455L", "C226Y", "D449A", 
  "N453A", "D457A", "D459A", "P446A", "F447A", 
  "I452A", "W439A", "F445A", 
  "K465E", "N297Q", "F72A", "W70A", "H124S", 
  "W598A", 'Q128S', 'N129S', 'K465E' 
)
extra_marg <- c(99, 100, 154, 63, 67, 70, 72, 95, 125, 63,64, 67,68, 70,72, 95, 97, 100, 106, 124, 125, 126, 127, 128, 129, 130, 131, 154)

marv_mutations_list_all <- unique(c(str_extract(marv_mutations, '\\d{2,3}'), extra_marg))
marv_mutations_list_all

marburg_all_likelihood <- read_csv('Marburg_GP1GP2_epistatic_accessibility_antigenicity.csv')

marburg_all_likelihood %>% 
  filter(position %in% marv_mutations_list_all) %>%
  select(label, position, ref, alt, significant_antigenicity_median) %>%
  rename(epi = significant_antigenicity_median) %>% 
  ggplot() +
  geom_tile(aes(x=as.factor(position), y=alt, fill=epi)) +
  scale_fill_viridis_c() +
  labs(title='Marburg Virus GP1GP2 Epistatic Antigenicity for Mutations Associated with Neutralization Resistance',
       x='Position',
       y='Amino Acid Substitution',
       fill='Epistatic Antigenicity') +
  theme_minimal() 


cchfv_mutations <- # List of individual residues from the antibody interaction table
  c(
    "C1165", "T1166", "H1187", "W1191", "N1194", "W1197", "C1198", "W1199", 
    "G1200", "V1201", "G1202", "T1203", "M1143", "S1145", "P1146", "V1147", 
    "F1148", "E1149", "K1225", "E1227", "Y1228", "I1229", "K1230", "P1276", 
    "E1277", "L1307", "Q1308", "S1309", "Y1321", "H1322", "T1346", "H1405", 
    "T1408", "Q1410"
  )
cchfv_mutations

cchfv_mutations_list <- c(
  "W1090R", "V1093E", "V1124E", "V1124G", "S1128R", "S1128N", "A1163T", 
  "W1185A", "R1189T", "W1199R", "W1199S", "G1204A", "G1204D", "I1229N", 
  "K1297N", "L1307R", "S1309N", "T1311I", "T1311N", "V1314M", "D1352N", 
  "D1352H", "D1352E", "G1353D", "G1353S", "N1386T", "N1386K", "K1393A", 
  "K1393N", "L1450Q", "I1462A", "E1500D", "D1504E", "L1518A", "F1520A", 
  "I1522N"
)
              

cchfv_all_likelihood %>% 
  mutate(position = position + 993) %>%
  filter(position %in% str_extract(cchfv_mutations_list, '\\d{3,4}')) %>% 
  select(label, position, ref, alt, significant_antigenicity_median) %>%
  rename(epi = significant_antigenicity_median) %>%
  ggplot() +
  geom_tile(aes(x=as.factor(position), y=alt, fill=epi)) +
  scale_fill_viridis_c() +
  labs(title='CCHFV GPC Epistatic Antigenicity for Mutations Associated with Neutralization Resistance',
       x='Position',
       y='Amino Acid Substitution',
       fill='Epistatic Antigenicity') +
  theme_minimal()

cchfv_all_likelihood %>% 
  mutate(position = position + 993) %>%
  
  filter(label %in% cchfv_mutations_list)

marv_mutations
