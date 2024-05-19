

sel_pars_gummy <- list(mesh = seq(4, 9, 0.01),
                 
                 theta1 = 184.3,
                 theta2 = 29739,
                 m50 = 1129,
                 m95 = 1344)

sel_pars_school <- list(mesh = seq(4, 9, 0.01),
                       
                        theta1 = 192,
                        theta2 = 67595,
                        m50 = 1421,
                        m95 = 1488)

peak_sel <- function(len, df){
  1 - (len / (df$rho * df$eta))^df$rho * exp(df$rho - len / df$eta)
}

sel_df <- cross_df(sel_pars_gummy) %>% 
  mutate(eta = -0.5 * (theta1 * mesh - sqrt(theta1^2 * mesh^2 + 4 * theta2))) %>% 
  mutate(rho = (theta1 * mesh) / eta) %>% 
  mutate(id = row_number()) %>%
  group_by(id) %>% nest(.key = "df")  %>% 
  mutate(peak_sel = map_dbl(df, ~ optimize(peak_sel, lower = 0, upper = 2500, df = .)$minimum)) %>% 
  unnest() %>% 
  mutate(peak_sel = round(peak_sel, 0)) %>% 
  mutate(P = (1 / (1 + exp(-log(19) * ((peak_sel - m50) / (m95 - m50)))))) %>% 
  mutate(P = round(P, 2))

# Selectivity parameters (Kirkwood and Walker 1986)

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
  c = c(1 / 2), # Proportion maternal
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
  Nsamples = c(1000)
)




data <-
  # Generate data frame of all combinations of scenarios
  cross_df(sim_pars_mustelus) %>%
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
  select(mesh_name, x, z, y) %>% 
  mutate(maternal = ifelse(z == 0, "Non maternal", "Maternal")) %>% 
  mutate(mesh_name = factor(mesh_name, levels = c("low", "high")))

# True underlying length at age and relative selectivity
sel <- tibble(a = seq(0, 54, 0.1)) %>% 
  expand_grid(mesh = sim_pars_mustelus$mesh) %>%
  arrange(mesh) %>% 
  mutate(len = sim_pars_mustelus$Linf *(1 - exp(-sim_pars_mustelus$K *(a - sim_pars_mustelus$t0)))) %>% 
  mutate(eta = -0.5 * (sim_pars_mustelus$theta1 * mesh - sqrt(sim_pars_mustelus$theta1^2 * mesh^2 + 4 * sim_pars_mustelus$theta2))) %>%
  mutate(rho = (sim_pars_mustelus$theta1 * mesh) / eta) %>% 
  mutate(sel  = (len / (rho * eta))^rho * exp(rho - len / eta)) %>% 
  mutate(mesh_name = case_when(
    mesh == 5.68 | mesh == 7.27 ~ "low",
    mesh == 6.57 | mesh == 7.53 ~ "high"
  )) %>% 
  mutate(mesh_name = factor(mesh_name, levels = c("low", "high")))



p <- ggplot() + 
  geom_histogram(data = data, aes(x = x, y=..count../sum(..count..), fill = maternal), col = "black", bins = 20, size = 0.1) + 
  scale_fill_manual(values = c("#fcbba1", "grey35")) +
  guides(fill = guide_legend(title = NULL)) +
  theme_classic() + 
  theme(panel.grid = element_blank(), axis.ticks = element_blank(),
        legend.position = c(0.1, 0.9),
        strip.text.x = element_blank() , 
        strip.background = element_blank()
        ) + 
  labs(x = NULL, y = NULL) +
  scale_x_continuous(labels = NULL, limits = c(500, 1800)) + scale_y_continuous(labels = NULL) + 
  facet_wrap(~mesh_name, scales = "free_y")



max_y <- max(ggplot_build(p)$data[[1]]$y)

p <- p + geom_line(data = sel, aes(x = len, y = sel * max_y))

ggsave("img/sel.pdf", p, width = 9, height = 3)

