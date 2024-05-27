# Housekeeping ----------------------------------------------------------------

source(here::here("code", "setup.R"))
source(here::here("code", "functions.R"))

compile("code/logistic3.cpp")
dyn.load(dynlib("code/logistic3"))
compile("code/logistic2.cpp")
dyn.load(dynlib("code/logistic2"))
compile("code/logistic2maternal.cpp")
dyn.load(dynlib("code/logistic2maternal"))


sim_pars_galeus <- list(
  species = "School shark",
  
  # Maximum age
  AMAX = c(27, 54), 
  # Natural mortality (Grant et al 1979)
  M = 0.1075,
  # Fishing mortality (Nominally set to 25 % of M)
  F_mort = 0.25 * 0.1075,
  # von Bertalanffy growth parameters (Grant et al 1979)
  K = 0.1600,
  Linf = 1618.3,
  t0 = -1.2818,
  cv_l = 0.05,
  # Maternity parameters (Walker 2005)
  m50 = 1421,
  m95 = 1488,
  c = c(1 / 3), # Proportion maternal
  # Maturity parameters (Walker 2005)
  l50 = 1349,
  l95 = 1502,
  # Uterine fecundity as function of length (Walker 2005)
  f = -46.0,
  g = 0.0491,
  # Selectivity parameters (Punt and Walker 1998)
  # Correspond to peak sel at 25% and 75% length maternity
  mesh = c(7.53),
  theta1 = 192,
  theta2 = 67595,
  
  # Number of samples
  Nsamples = c(100, 1000, 2500)
)

data <-
  # Generate data frame of all combinations of scenarios
  cross_df(sim_pars_galeus) %>%
  # Create some additional variables
  mutate(mesh_name = case_when(
    mesh == 5.68 | mesh == 7.27 ~ "low",
    mesh == 6.57 | mesh == 7.53 ~ "high"
  )) %>%
  mutate(repro_freq = case_when(
    c == 1 ~ "1 (Annual)",
    c == 3 / 4 ~ "0.75 (Annual / biennial)",
    c == 2 / 3 ~ "0.66 (Annual / biennial)",
    c == 0.5 ~ "0.5 (Biennial)",
    c == 1 / 3 ~ "0.33 (Triennial)",
    c == 0.25 ~ "0.25 (Quadrennial)"
  )) %>%
  # Give unique row numbers to each simulation
  mutate(id = paste(row_number(), species, Nsamples, mesh_name, round(c, 2), sep = "-")) %>%
  # Make into nested list column
  group_by(id) %>%
  nest() %>%
  ungroup() %>% 
  rename(pars = data) %>%
  # Create and number n repliactes of each row
  mutate(reps = 1) %>% 
  uncount(reps) %>%
  group_by(id) %>%
  mutate(iteration = row_number()) %>%
  ungroup() %>% 
  # Generate simulated length data
  mutate(data = map2(pars, id, generate_data)) %>%   mutate(data = map2(pars, data, generate_maternal_data)) %>% # How many maternal females are there in the simulated data
  mutate(n_maternal = map_dbl(data, ~ sum(.$z))) %>%
  # Are there zero maternal females?
  mutate(no_maternal = ifelse(n_maternal == 0, 1, 0)) %>% 
  # How many immature females are in the simulated data
  mutate(n_immature = map_dbl(data, ~ sum(.$y == 0))) %>%
  # What is the range of x values (for plotting)?
  mutate(xmin = map_dbl(data, ~ min(.$x))) %>%
  mutate(xmax = map_dbl(data, ~ max(.$x))) %>% 
  unnest(cols = c(pars, data)) %>% 
  select(mesh_name, AMAX, Nsamples, x, z, y) %>% 
  mutate(maternal = ifelse(z == 0, "Non maternal", "Maternal")) %>% 
  mutate(mesh_name = factor(mesh_name, levels = c("low", "high"))) %>% 
  mutate(AMAX = as.factor(AMAX))


p <- ggplot() + 
  geom_density(data = data, aes(x = x, group = AMAX, fill = AMAX), alpha = 0.5) + 
  facet_wrap(~Nsamples) + 
  labs(x = "Length")
