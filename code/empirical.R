# Housekeeping ----------------------------------------------------------------

source(here::here("code", "setup.R"))
source(here::here("code", "functions.R"))

# Sandbar shark analysis ------------------------------------------------------

compile("code/logistic3.cpp")
dyn.load(dynlib("code/logistic3"))

# Load data
# Complete dataset from Ivy, includes some smaller females that were missing from previous spreadsheet
data_complete <- readxl::read_xlsx(here("data", "plumbeus-unprocessed.xlsx"))

# Incomplete dataset that has full repro details
data_incomplete <- readxl::read_xls(here("data", "plumbeus-extra-unprocessed.xls"), sheet = 6) %>%
  select(
    ID = `SBAR .`, `OVARY L (mm)`, `OVARY WIDTH`, `OVARY WEIGHT (g)`, `TOTAL OOCYTES`, `. VIT`,
    `DIAMETER LARGEST OO`, `UTERUS WIDTH`, `NIDAMENTAL WIDTH`, `PUPS?`, `. LEFT`, `. RIGHT`
  ) %>%
  arrange(ID) %>%
  # Remove #358 as it is a duplicate, doesn't seem to have been included in final analysis
  filter(!row_number() %in% 281)

# Load data for pups
pups <- read_csv(here("data", "plumbeus-pups-unprocessed.csv")) %>% 
  rename_all(tolower) %>%
  select(ID = maternal.id, puplen = `stl.cm`) %>% 
  mutate(puplen = as.numeric(puplen))

# Treat July as beginning of a new year / season

# SBAR1242 stage 4 just ovulated and pregnant, but June so count as resting for this year
# SBAR1248 stage 4 just ovulated and pregnant, but June so count as resting for this year
# SBAR1244 stage 4 just ovulated and pregnant, but June so count as resting for this year
# SBAR471 non-pregnant, perhaps pupped late?
# SBAR353 stage 6 postpartum but month is September so would have been maternal last season
# SBAR303 possible stage 6 but definitely not in maternal condition
# SBAR376 stage 6 postpartum but month is September so would have been maternal last season
# SBAR219 stage 6 postpartum but month is August so would have been maternal last season
# SBAR478 stage 6 postpartum but month is November so would have been maternal last season

data <- left_join(data_complete, data_incomplete) %>%
  mutate(month = month(Date)) %>%
  rename(MOD = `DIAMETER LARGEST OO`) %>%
  mutate(z = ifelse(Stage %in% c(4, 5, 6), 1, 0)) %>%
  mutate(z = ifelse(ID %in% c(
    "SBAR1242", "SBAR1248", "SBAR1244", "SBAR471",
    "SBAR353", "SBAR303", "SBAR376",
    "SBAR219", "SBAR478"
  ), 0, z)) %>%
  left_join(pups) %>% 
  select(ID, month, MOD, Stage, x = FL, z, puplen) %>%
  group_by(ID, month, MOD, Stage, x, z) %>% 
  summarise(puplen = mean(puplen)) %>% 
  filter(!is.na(x), !is.na(z)) 

# Look at MOD and embryo length vs month to check maternal classifications
ggplot() +
  geom_point(data = filter(data, !Stage %in% c(1, 2)), aes(x = month, y = MOD, col = as.factor(Stage))) +
  facet_wrap(~z, nrow = 2) +
  scale_x_continuous(breaks = seq(0, 12, 2))

ggplot() +
  geom_point(data = filter(data, Stage %in% c(4, 5)), aes(x = month, y = puplen, col = as.factor(Stage))) +
  scale_x_continuous(breaks = seq(0, 12, 2))



# Load data from Andrew
data_complete2 <- readxl::read_xlsx(here("data", "plumbeus-unprocessed-piercy.xlsx")) %>% 
  select(month, MOD, puplen = emb_tl, Stage = stage, x = FL, z = mat_stage) %>% 
  mutate(MOD = as.numeric(MOD))

# Look at MOD and embryo length vs month to check maternal classifications
ggplot() +
  geom_point(data = filter(data_complete2, Stage %in% c(3:6), !is.na(MOD)), aes(x = month, y = MOD, col = as.factor(Stage))) +
  facet_wrap(~z, nrow = 2) +
  scale_x_continuous(breaks = seq(0, 12, 2))

ggplot() +
  geom_point(data = filter(data_complete2, Stage %in% c(4, 5)), aes(x = month, y = puplen, col = as.factor(Stage))) +
  scale_x_continuous(breaks = seq(0, 12, 2))

data2 <- select(data_complete2, month, Stage, x, z) %>%
  filter(!is.na(x), !is.na(z))

# Joing Ivy and Andrew's data

# Tidy up data 
data <- data %>% 
  ungroup() %>% select(month, Stage, x, z) %>%
  rbind(data2) %>% 
  write_rds(here("data", "empirical-plumbeus.Rds"))

# Run analyses ----------------------------------------------------------------

# 3PLF Free method
method <- "3"

# Parameter list
pars <- list(m50 = 160, m95 = 170, c = 0.5)

# Create objective function
obj <- MakeADFun(
  data = data,
  parameters = pars,
  DLL = "logistic3"
)

# Optimize
opt <- nlminb(obj$par, obj$fn, obj$gr, lower = c(0, 0, 0), upper = c(Inf, Inf, 1))

# Calculate AIC
opt_aic <- 2 * length(opt$par) + 2 * opt$objective

# Calculate 95% bootstrap confidence intervals
output <- bootstraps(data, 10000) %>%
  mutate(n_pregnant = map_dbl(splits, ~ sum(analysis(.)$z))) %>%
  filter(n_pregnant > 0) %>%
  mutate(opt = map(splits, boot_helper,
                   obj = obj,
                   prop_pregnant = prop_pregnant, method = method
  )) %>%
  mutate(convergence = map_dbl(opt, ~ .$convergence)) %>%
  filter(convergence == 0) %>%
  mutate(pars = map(opt, ~ .$par)) %>%
  mutate(par_name = map(opt, ~ .$par %>% names()))

# Generate confidence intervals for fitted model
preds1 <- output %>%
  mutate(x = map(splits, ~ seq(min(analysis(.)$x) %>% floor(), max(analysis(.)$x) %>% ceiling(), 1))) %>%
  mutate(y = map2(pars, x, ~ (.x[3] / (1 + exp(-log(19) * ((.y - .x[1]) / (.x[2] - .x[1]))))))) %>%
  select(id, x, y) %>%
  unnest(cols = c(x, y)) %>%
  group_by(x) %>%
  summarise(lower = quantile(y, 0.05), upper = quantile(y, 0.975)) %>%
  mutate(y = (opt$par[3] / (1 + exp(-log(19) * ((x - opt$par[1]) / (opt$par[2] - opt$par[1])))))) %>%
  mutate(label = "3PLF - free") %>%
  mutate(species = "C. plumbeus")

# Get 95% confidence intervals around each estimated parameter
summary1 <- output %>%
  unnest(c(pars, par_name)) %>%
  group_by(par_name) %>%
  summarise(
    lower = quantile(pars, 0.05 / 2),
    upper = quantile(pars, 1 - 0.05 / 2)
  ) %>%
  cbind(est = c(opt$par[3], opt$par[-3])) %>%
  mutate(method = "3PLF - free") %>%
  mutate(species = "C. plumbeus") %>%
  mutate(AIC = opt_aic)

# 3PLF fixed method
method <- "3a"

# Parameter list
prop_pregnant <- 0.5
pars <- list(m50 = 160, m95 = 170, c = prop_pregnant)

# Create objective function
obj2 <- MakeADFun(
  data = data,
  parameters = pars,
  DLL = "logistic3",
  map = list(c = factor(NA))
)

# Optimize
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr, lower = c(0, 0, 0), upper = c(Inf, Inf, 1))

# Calculate AIC
opt2_aic <- 2 * length(opt2$par) + 2 * opt2$objective

# Calculate 95% bootstrap confidence intervals
output <- bootstraps(data, 10000) %>%
  mutate(n_pregnant = map_dbl(splits, ~ sum(analysis(.)$z))) %>%
  filter(n_pregnant > 0) %>%
  mutate(opt = map(splits, boot_helper,
                   obj = obj2,
                   prop_pregnant = prop_pregnant, method = method
  )) %>%
  mutate(convergence = map_dbl(opt, ~ .$convergence)) %>%
  filter(convergence == 0) %>%
  mutate(pars = map(opt, ~ .$par)) %>%
  mutate(pars = map(pars, ~ c(., prop_pregnant))) %>%
  mutate(par_name = map(opt, ~ .$par %>% names())) %>%
  mutate(par_name = map(par_name, ~ c(., "c")))

# Generate confidence intervals for fitted model
preds2 <- output %>%
  mutate(x = map(splits, ~ seq(min(analysis(.)$x) %>% floor(), max(analysis(.)$x) %>% ceiling(), 1))) %>%
  mutate(y = map2(pars, x, ~ (.x[3] / (1 + exp(-log(19) * ((.y - .x[1]) / (.x[2] - .x[1]))))))) %>%
  select(id, x, y) %>%
  unnest(cols = c(x, y)) %>%
  group_by(x) %>%
  summarise(lower = quantile(y, 0.05), upper = quantile(y, 0.975)) %>%
  mutate(y = (prop_pregnant / (1 + exp(-log(19) * ((x - opt$par[1]) / (opt$par[2] - opt$par[1])))))) %>%
  mutate(label = "3PLF - fixed (PMax = 0.5)") %>%
  mutate(species = "C. plumbeus")

# Get 95% confidence intervals around each estimated parameter
summary2 <- output %>%
  unnest(c(pars, par_name)) %>%
  group_by(par_name) %>%
  summarise(
    lower = quantile(pars, 0.05 / 2),
    upper = quantile(pars, 1 - 0.05 / 2)
  ) %>%
  cbind(est = c(prop_pregnant, opt2$par)) %>%
  mutate(method = "3PLF - fixed") %>%
  mutate(species = "C. plumbeus") %>%
  mutate(AIC = opt2_aic)

# 3PLF fixed method (PMax = 0.33)
method <- "3a"

# Parameter list
prop_pregnant <- 1/3
pars <- list(m50 = 160, m95 = 170, c = prop_pregnant)

# Create objective function
obj3 <- MakeADFun(
  data = data,
  parameters = pars,
  DLL = "logistic3",
  map = list(c = factor(NA))
)

# Optimize
opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr, lower = c(0, 0, 0), upper = c(Inf, Inf, 1))

# Calculate AIC
opt3_aic <- 2 * length(opt3$par) + 2 * opt3$objective

# Calculate 95% bootstrap confidence intervals
output <- bootstraps(data, 10000) %>%
  mutate(n_pregnant = map_dbl(splits, ~ sum(analysis(.)$z))) %>%
  filter(n_pregnant > 0) %>%
  mutate(opt = map(splits, boot_helper,
                   obj = obj3,
                   prop_pregnant = prop_pregnant, method = method
  )) %>%
  mutate(convergence = map_dbl(opt, ~ .$convergence)) %>%
  filter(convergence == 0) %>%
  mutate(pars = map(opt, ~ .$par)) %>%
  mutate(pars = map(pars, ~ c(., prop_pregnant))) %>%
  mutate(par_name = map(opt, ~ .$par %>% names())) %>%
  mutate(par_name = map(par_name, ~ c(., "c")))

# Generate confidence intervals for fitted model
preds3 <- output %>%
  mutate(x = map(splits, ~ seq(min(analysis(.)$x) %>% floor(), max(analysis(.)$x) %>% ceiling(), 1))) %>%
  mutate(y = map2(pars, x, ~ (.x[3] / (1 + exp(-log(19) * ((.y - .x[1]) / (.x[2] - .x[1]))))))) %>%
  select(id, x, y) %>%
  unnest(cols = c(x, y)) %>%
  group_by(x) %>%
  summarise(lower = quantile(y, 0.05), upper = quantile(y, 0.975)) %>%
  mutate(y = (prop_pregnant / (1 + exp(-log(19) * ((x - opt$par[1]) / (opt$par[2] - opt$par[1])))))) %>%
  mutate(label = "3PLF - fixed (PMax = 0.33)") %>%
  mutate(species = "C. plumbeus")

# Get 95% confidence intervals around each estimated parameter
summary3 <- output %>%
  unnest(c(pars, par_name)) %>%
  group_by(par_name) %>%
  summarise(
    lower = quantile(pars, 0.05 / 2),
    upper = quantile(pars, 1 - 0.05 / 2)
  ) %>%
  cbind(est = c(prop_pregnant, opt3$par)) %>%
  mutate(method = "3PLF - fixed") %>%
  mutate(species = "C. plumbeus") %>%
  mutate(AIC = opt3_aic)

# Model predictions
preds <- rbind(preds1, preds2, preds3)

# Save analysis
write_rds(preds, here("data", "empirical-predictions.Rds"))

# Model summaries
summary <- rbind(summary1, summary2, summary3) %>%
  mutate_at(c("est", "lower", "upper"), signif, 3) %>%
  mutate(lower = ifelse(par_name %in% "c" & method %in% "3PLF - fixed", NA, lower)) %>%
  mutate(upper = ifelse(par_name %in% "c" & method %in% "3PLF - fixed", NA, upper)) %>%
  mutate(par = ifelse(!is.na(upper), paste0(est, " (", lower, " - ", upper, ")"), paste0(est, "*"))) %>%
  select(method, AIC, par_name, par) %>%
  pivot_wider(names_from = par_name, values_from = par) %>%
  rename(PMAX = c, L50 = m50, L95 = m95) %>%
  relocate(AIC, .after = last_col()) %>%
  mutate(delta_AIC = AIC - min(AIC)) %>%
  mutate(w_AIC = 100 * exp(-delta_AIC / 2) / sum(exp(-delta_AIC / 2))) %>%
  arrange(delta_AIC)

write_rds(summary, here("data", "empirical-summary.Rds"))
