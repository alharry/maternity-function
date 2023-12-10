# Housekeeping ----------------------------------------------------------------

source(here::here("code", "setup.R"))
source(here::here("code", "functions.R"))

compile("code/logistic3.cpp")
dyn.load(dynlib("code/logistic3"))
compile("code/logistic2.cpp")
dyn.load(dynlib("code/logistic2"))
compile("code/logistic2maternal.cpp")
dyn.load(dynlib("code/logistic2maternal"))

# Parameter section  ----------------------------------------------------------

sim_pars_galeus <- list(
  species = "School shark",
  
  # Maximum age
  AMAX = 54, 
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
  c = c(1, 3 / 4, 2 / 3, 1 / 2, 1 / 3, 1 / 4), # Proportion maternal
  # Maturity parameters (Walker 2005)
  l50 = 1349,
  l95 = 1502,
  # Uterine fecundity as function of length (Walker 2005)
  f = -46.0,
  g = 0.0491,
  # Selectivity parameters (Punt and Walker 1998)
  # Correspond to peak sel at 25% and 75% length maternity
  mesh = c(7.27, 7.53),
  theta1 = 192,
  theta2 = 67595,
  
  # Number of samples
  Nsamples = c(50, 100, 250, 500, 1000, 2500)
)

sim_pars_mustelus <- list(
  species = "Gummy shark",
  
  # Maximum age
  AMAX = 20, 
  # Natural mortality (Max age reported in Moulton et al 1992, Walker 1992 M)
  M = 0.1970,
  # Fishing mortality (Nominally set to 25 % of M)
  F_mort = 0.25 * 0.1970,
  # von Bertalanffy growth parameters (Moulton et al 1992, BS 1973-1976)
  K = 0.123,
  Linf = 2019,
  t0 = -1.55,
  cv_l = 0.05,
  # Maternity parameters (Walker 2007 EKI 73-76)
  m50 = 1129,
  m95 = 1344,
  c = c(1, 3 / 4, 2 / 3, 1 / 2, 1 / 3, 1 / 4), # Proportion maternal
  # Maturity parameters (Walker 2007 EKI 73-76)
  l50 = 1105,
  l95 = 1293,
  # Uterine fecundity as function of length (Walker 2007) F = f*exp(g x length)
  f = 0.2804,
  g = 0.00286,
  # Selectivity parameters (Kirkwood and Walker 1986)
  # Correspond to peak sel at 25% and 75% length maternity
  mesh = c(5.68, 6.57),
  theta1 = 184.3,
  theta2 = 29739,

  # Number of samples
  Nsamples = c(50, 100, 250, 500, 1000, 2500)
)

# Simulation study  -----------------------------------------------------------
data <-
  # Generate data frame of all combinations of scenarios
  rbind(cross_df(sim_pars_mustelus), cross_df(sim_pars_galeus)) %>%
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
  mutate(reps = 300) %>% 
  uncount(reps) %>%
  group_by(id) %>%
  mutate(iteration = row_number()) %>%
  ungroup() %>% 
  # Generate simulated length data
  mutate(data = map2(pars, id, generate_data))

write_rds(data, "data/generated_data.rds")
data <- read_rds("data/generated_data.rds")


data <- data %>%
  slice(41401:43200) %>% 
  # Simulate maternal data
  mutate(data = map2(pars, data, generate_maternal_data)) %>% # How many maternal females are there in the simulated data
  mutate(n_maternal = map_dbl(data, ~ sum(.$z))) %>%
  # Are there zero maternal females?
  mutate(no_maternal = ifelse(n_maternal == 0, 1, 0)) %>% 
  # How many immature females are in the simulated data
  mutate(n_immature = map_dbl(data, ~ sum(.$y == 0))) %>%
  # What is the range of x values (for plotting)?
  mutate(xmin = map_dbl(data, ~ min(.$x))) %>%
  mutate(xmax = map_dbl(data, ~ max(.$x))) %>%
    # Run simulations
    group_by(id) %>%
    mutate(iteration = row_number()) %>%
    ungroup() %>%
    mutate(results = pmap(., run_sims)) %>%
    mutate(xmin = map_dbl(data, ~ min(.$x))) %>%
    mutate(xmax = map_dbl(data, ~ max(.$x))) %>%
    select(-data) %>%
    unnest(cols = c(pars)) %>%
    # Not sure what happened here, but couldn't unnest...
    mutate(results = map(results, ~ as_tibble(.))) %>%
    unnest(cols = c(results)) %>%
    select(-pars) %>% 
    group_by(id) %>%
    nest() %>%
    # Plot data
    pwalk(., plot_sim_data) %>%
    mutate(r_0 = map_dbl(data, ~ integrate(r_0, 0, Inf, df = .[1, ])$value)) %>%
    unnest(cols= c(data)) %>% 
    group_by(id, iteration, method) %>%
    mutate(par_true = case_when(
      par_name %in% "r_0" ~ r_0,
      par_name %in% "m50" ~ m50,
      par_name %in% "m95" ~ m95,
      par_name %in% "c" ~ c
    )) %>%
    # Median relative error (bias)
    mutate(rel_error = 100 * ((par_est - par_true) / par_true)) %>%
    # Median absolute relative error (precision)
    mutate(abs_error = abs(rel_error)) %>%
    # 50 % interval coverage
    mutate(int_coverage = ifelse(par_true >= lower & par_true <= upper, TRUE, FALSE)) %>%
    select(-xmin, -xmax, -data)

write_rds(data, here("data", "simulations-school-7.Rds"))


data <- read_rds("data/simulations-gummy.rds")
data1 <- read_rds("data/simulations-school-1.rds")
data2 <- read_rds("data/simulations-school-2.rds")
data3 <- read_rds("data/simulations-school-3.rds")
data4 <- read_rds("data/simulations-school-4.rds")
data5 <- read_rds("data/simulations-school-5.rds")
data6 <- read_rds("data/simulations-school-6.rds")
data7 <- read_rds("data/simulations-school-7.rds")

data <- rbind(data, data1) %>% rbind(data2) %>% rbind(data3) %>% 
  rbind(data4) %>% rbind(data5) %>% rbind(data6) %>% rbind(data7)

write_rds(data, here("data", "simulations.Rds"))


