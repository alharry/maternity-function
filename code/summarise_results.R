# Housekeeping ----------------------------------------------------------------

source(here::here("code", "setup.R"))
source(here::here("code", "functions.R"))

# Summarise simulation study --------------------------------------------------
data <- read_rds(here("data", "simulations.Rds")) %>% 
  group_by(id) %>% nest() %>% 
  # How many iterations in each scenario had zero pregnant females and were excluded? 
  mutate(zero_maternal = map_dbl(data, ~ sum(.$no_maternal))) %>% 
  unnest(cols = c(data)) %>% 
  ungroup() %>% 
  filter(no_maternal == 0) %>% select(-no_maternal) %>% 
  # How many other iterations in each scenario failed due to lack of pregnant females? 
  mutate(few_maternal = ifelse(is.na(method), 1, 0)) %>% 
  group_by(id) %>% 
  nest() %>% 
  mutate(low_maternal = map_dbl(data, ~ sum(.$few_maternal))) %>% 
  unnest(cols = c(data)) %>% 
  ungroup() %>% 
  filter(!is.na(method)) %>% 
  # Calculate total number of iterations that failed to start due to lack of data
  mutate(maternal_data_fail = low_maternal + zero_maternal) %>% 
  # Calculate total number of iterations done for each simulation
  mutate(ntotal = 250 - maternal_data_fail) %>% 
  select(-low_maternal, -few_maternal, -zero_maternal) %>% 
  # Now tally parameter boundary failures and non-convergence failures
  group_by(id, iteration, method) %>%
  nest() %>% 
  # Remove iterations where boundary was hit
  mutate(hit_boundary = map_lgl(data, par_check)) %>%
  filter(hit_boundary %in% FALSE) %>% 
  select(-hit_boundary) %>% 
  unnest(cols = c(data)) %>% 
  # Remove iterations where convergence failed
  filter(convergence == 0) %>% 
  mutate(method_long = case_when(
    method %in% "2" ~ "2PLF - maturity",
    method %in% "2a" ~ "2PLF - maternity",
    method %in% "3" ~ "3PLF - free",
    method %in% "3a" ~ "3PLF - fixed"
  )) %>%
  group_by(species, method_long, repro_freq, mesh_name, Nsamples, ntotal) %>%
  nest() %>% 
  mutate(nsims = map_dbl(data, nrow) / 4) %>%
  mutate(conv_success = 100 * (1 - (ntotal - nsims) / ntotal)) %>%
  unnest(cols = c(data)) %>%
  group_by(species, method_long, repro_freq, mesh_name, Nsamples, ntotal, conv_success, par_name) %>%
  summarise(
    sd_re = sd(rel_error), rel_error = mean(rel_error),
    sd_ae = sd(abs_error), abs_error = mean(abs_error),
    int_coverage = sum(int_coverage),
    conv_success = mean(conv_success), 
    prop_pregnant = mean(prop_pregnant),
    n_maternal = median(n_maternal),
    n_immature = median(n_immature)
  ) %>%
  mutate(int_coverage = 100 * int_coverage / ntotal) %>% 
  ungroup() %>%
  mutate(method_long = fct_relevel(
    method_long, "3PLF - free",
    "3PLF - fixed",
    "2PLF - maternity",
    "2PLF - maturity"
  )) %>%
  mutate(mesh_name_long = case_when(
    mesh_name %in% "high" ~ "High selectivity",
    mesh_name %in% "low" ~ "Low selectivity"
  ))

write_rds(data, here("data", "simulations-summarised.Rds"))

# Pivot table of simulation study --------------------------------------------------
data <- read_rds(here("data", "simulations.Rds")) %>%
  mutate(method_long = case_when(
    method %in% "2" ~ "2PLF - maturity",
    method %in% "2a" ~ "2PLF - maternity",
    method %in% "3" ~ "3PLF - free",
    method %in% "3a" ~ "3PLF - fixed"
  )) %>%
  mutate(method_long = fct_relevel(
    method_long, "3PLF - free",
    "3PLF - fixed",
    "2PLF - maternity",
    "2PLF - maturity"
  )) %>%
  filter(method %in% c("3", "3a")) %>%
  filter(convergence == 0) %>%
  select(id, species, method_long, iteration, Nsamples, mesh_name, repro_freq, par_name, par_est, rel_error) %>%
  pivot_wider(names_from = par_name, values_from = c(par_est, rel_error)) %>%
  ungroup() %>%
  mutate(repro_freq = as.factor(repro_freq)) %>%
  arrange(repro_freq) %>%
  mutate(ratio = par_est_m50 / par_est_m95)

write_rds(data, here("data", "simulations-pivot.Rds"))

