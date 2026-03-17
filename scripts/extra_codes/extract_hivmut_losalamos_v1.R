env_hiv$`bNAb signature predictions`[which(!is.na(env_hiv$`bNAb signature predictions`))]

env_hiv %>% 
  select( "HXB2 position", "HXB2 AA", `bNAb signature predictions`) %>% 
  filter(!is.na(`bNAb signature predictions`)) %>% 
  filter(str_detect(`bNAb signature predictions`, 'increased resistance to neutralization')) %>%
  mutate(
    position = str_extract_all(`bNAb signature predictions`, '\\b[A-Z]\\d{2,3}\\b is associated with increased resistance to neutralization')) %>% 
  pull(position) 


my_fun <- function(x) deparse(substitute(x))

my_fun(T105E)

lapply(c(H84Y), my_fun)
l <- c(H84Y, T105E, K146Q, I175M, D229N, K230E, K231T, N233N, P237P, P239S, N300N, R348R, S362H, S362P, L367I, L367L, N384N, Q400T, V439E, G468E, Q637G)

hiv_ab_expmut <- c('H84Y', 'T105E', 'K146Q', 'I175M', 'D229N', 'K230E', 'K231T', 'N233N', 'P237P', 'P239S', 'N300N', 'R348R', 'S362H', 'S362P', 'L367I', 'L367L', 'N384N', 'Q400T', 'V439E', 'G468E', 'Q637G')


tibble::lst(H84Y, T105E, K146Q, I175M, D229N, K230E, K231T, N233N, P237P, P239S, N300N, R348R, S362H, S362P, L367I, L367L, N384N, Q400T, V439E, G468E, Q637G)


c(H84Y, T105E, K146Q)
for (i in (H84Y, T105E, K146Q, I175M, D229N, K230E, K231T, N233N, P237P, P239S, N300N, R348R, S362H, S362P, L367I, L367L, N384N, Q400T, V439E, G468E, Q637G)) {print (i)}

env_hiv %>% 
  select( "HXB2 position", "HXB2 AA", starts_with('Mutations associated with resistance to neutralization')) %>% 
  filter(!rowSums(is.na(.)) == 4) %>% 
  # mutate(across(everything(), ~ replace_na(.x, "NM"))) %>%
  # mutate(mutation = paste(., collapse = "; ")) 
  unite(mutation, starts_with('Mutations'), sep = "; ", na.rm = TRUE) %>% 
  write_csv('hiv_ref_ab_escape.csv')
  

  

hiv_all_exp_mut <- c(hiv_exp_mut, hiv_ab_expmut)


env_hiv[31, 7:8]

hiv_ultimate %>% select(position, ref) %>% distinct(position, .keep_all = T) %>% print(n=670)

hiv_mutt <- read_csv('mutation_match_hiv.csv')
str_match(hiv_mutt$escapy_mutation_label[1], '\\b[A-Z]\\d{2,3}\\b[A-Z]\\b')

hiv_exp_mut <-unlist(str_extract_all(hiv_mutt$escapy_mutation_label, "\\b[A-Z]\\d{2,3}[A-Z]\\b"))

hiv_all_exp_mut <- data.frame('label' = hiv_all_exp_mut)
hiv_all_exp_mut <- hiv_all_exp_mut %>% mutate(position = str_extract_all(label, '\\d{2,3}'))
write_csv(hiv_all_exp_mut, 'hiv_all_exp_mut.csv')                           
                           