# Utility functions -----------------------------------------------------------

## Parameter transformation
inv_logit <- function(x) {
  1 / (1 + exp(-x))
}
logit <- function(x) {
  -log((1 - x) / x)
}

## Save data
save_sim_data <- function(df, id) {
  dir_out <- paste0(here(), "/data/simulations/", id)
  if (!dir.exists(dir_out)) {
    dir.create(dir_out)
  }
  write_rds(df, paste0(dir_out, "/population-structure.rds"))
}

## The continuous time Lotka equation (Xiao and Walker 2000)
continuous_lotka <- function(a, Lambda, df) {

  # Length at age function
  len <- function(a) {
    df$Linf * (1 - exp(-df$K * (a - df$t0)))
  }

  # Maternity at age function
  maternity <- function(a) {
    (df$c / (1 + exp(-log(19) * ((len(a) - df$m50) / (df$m95 - df$m50)))))
  }

  # Fecundity at age function
  fecundity <- function(a) {
    if (df$species == "School shark") {
      df$f + df$g * len(a)
    } else {
      df$f * exp(df$g * len(a))
    }
  }

  # Birth function
  birth <- function(a) {
    if (df$species == "School shark") {
      ifelse(a >= 11.2404, fecundity(a) * maternity(a) * (1 / 2), 0)
    } else {
      ifelse(a >= 3.9694, fecundity(a) * maternity(a) * (1 / 2), 0)
    }
  }

  # Intrinsic rate of population decrease with age
  LAMBDA <- function(a) {
    rep(Lambda, length(a))
  }

  # Sum to optimize
  birth(a) * exp(sapply(a, function(x) {
    integrate(LAMBDA, 0, x)$value
  }))
}

## Function for solving Lotka equation
solve_lotka <- function(df, Lambda = 0) {
  Integral <- integrate(continuous_lotka, 0, Inf, Lambda = Lambda, df = df)$value
  dif <- abs(1 - Integral)
  return(dif)
}

## Function for optimizing Lotka equation
optimize_lotka <- function(df, lower = -2, upper = 0) {
  lambda <- optimize(solve_lotka, c(lower, upper), df = df)$minimum
  return(list(lambda = lambda))
}

## Function for generating simulated length and age data
generate_data <- function(df, id, Return = FALSE) {
  df <- df %>%
    mutate(lambda = optimize_lotka(df)$lambda)

  # Convert selectivity parameters
  df <- df %>%
    mutate(eta = -0.5 * (theta1 * mesh - sqrt(theta1^2 * mesh^2 + 4 * theta2))) %>%
    mutate(rho = (theta1 * mesh) / eta)

  # Initialise an age structure based on the
  # stable age distribution
  N_stable <- -df$lambda * exp(df$lambda * 0:100)
  N_stable <- ((N_stable / sum(N_stable)) * df$N0) %>% round()

  ages <- 0:100

  dat <- tibble(N = c(rep(1, length(N_stable)), N_stable), age = c(ages, ages)) %>%
    group_by(age) %>%
    expand(ID = full_seq(N, 1), age) %>%
    filter(age %in% ages[which(N_stable != 0)]) %>%
    select(-ID) %>%
    ungroup() %>%
    mutate(len = df$Linf * (1 - exp(-df$K * (age - df$t0)))) %>%
    mutate(len = len + rnorm(n(), 0, df$cv_l * len)) %>%
    mutate(sel = (len / (df$rho * df$eta))^df$rho * exp(df$rho - len / df$eta)) %>%
    save_sim_data(id = id)
  
  return(dat)
}

## Randomly sample from length structure
sample_data <- function(df, data){
  data %>% 
  sample_n(df$Nsamples, weight = sel) %>%
  arrange(len) %>%
  select(x = len) %>%
  mutate(y = NA) %>%
  arrange(x)
}

## Optimize objective function
opt_fun <- function(obj) {
  nlminb(obj$par, obj$fn, obj$gr, lower = c(0, 0, 0), upper = c(10000, 10000, 1))
}

## Generate maternal data by simulating from logistic functions
generate_maternal_data <- function(pars, data){
  data_mature <- MakeADFun(data, list(l50 = pars$l50, l95 = pars$l95), DLL = "logistic2")$simulate(complete = TRUE) %>% 
    as_tibble()
  
  MakeADFun(data %>% mutate(z = NA), list(c = pars$c, m50 = pars$m50, m95 = pars$m95), DLL = "logistic3")$simulate(complete = TRUE) %>% 
    as_tibble() %>% select(-y) %>% 
    bind_cols("y" = data_mature$y)
}

# #
# id = data$id
# iteration <- data$iteration
# n_maternal <- data$n_maternal
# pars <- data$pars[[1]]
# data <- data$data[[1]]


## Run simulations
run_sims <- function(data, pars, n_maternal, id, iteration, ...) {

  if(n_maternal > 0){

  # Create TMB objective functions for 2 parameter model
  tmb_obj2 <- MakeADFun(data, list(l50 = pars$l50, l95 = pars$l95), DLL = "logistic2", silent = TRUE, hessian = T)
  tmb_opt2 <- opt_fun(tmb_obj2)
  
  # Procedure for 'guessing' the proportion pregnant
  # Calculate proportion of pregnant females at lengths > 99% maturity
  # If there are 0 pregnant females or that doesn't work go to lengths > 95%, 50%, 25% 
  # then % of all mature females
  l_50 <- tmb_opt2$par[[1]]
  l_95 <- tmb_opt2$par[[2]]
  beta1 = l_50 * -log(19) / (l_95 - l_50) 
  beta2 = log(19) / (l_95 - l_50)
  l_99 = log(1 /0.01 - 1) / beta2 - beta1/beta2
  l_25 = log(1/0.75 -1) / beta2 - beta1/beta2
  
  prop_pregnant_99 <- min(sum(data$z[which(data$x > l_99)]) /
                            sum(data$y[which(data$x > l_99)]), 1)
  prop_pregnant_95 <- min(sum(data$z[which(data$x > l_95)]) /
                            sum(data$y[which(data$x > l_95)]), 1)
  prop_pregnant_50 <- min(sum(data$z[which(data$x > l_50)]) /
                            sum(data$y[which(data$x > l_50)]), 1)
  prop_pregnant_25 <- min(sum(data$z[which(data$x > l_25)]) /
                            sum(data$y[which(data$x > l_25)]), 1)
  prop_pregnant_all = sum(data$z[which(data$y == 1)]) /
                                 sum(data$y[which(data$y == 1)])
  prop_pregnant = ifelse(is.na(prop_pregnant_99) | prop_pregnant_99 == 0,  prop_pregnant_95, prop_pregnant_99)
  prop_pregnant = ifelse(is.na(prop_pregnant_95) | prop_pregnant_95 == 0,  prop_pregnant_50, prop_pregnant)
  prop_pregnant = ifelse(is.na(prop_pregnant_50) | prop_pregnant_50 == 0, prop_pregnant_25, prop_pregnant)
  prop_pregnant = ifelse(is.na(prop_pregnant_25) | prop_pregnant_25 == 0, prop_pregnant_all, prop_pregnant)

  if(prop_pregnant > 0){
  # Create 3 parameter objective function
  tmb_obj3 = MakeADFun(data, list(c = prop_pregnant, m50 = pars$m50, m95 = pars$m95), 
                       DLL = "logistic3", silent = TRUE) 
  tmb_opt3 <- opt_fun(tmb_obj3)

  # Create 3 parameter w/ fixed asymptote objective function
  tmb_obj3a = MakeADFun(data, list(c = prop_pregnant, m50 = pars$m50, m95 = pars$m95), 
                        DLL = "logistic3", map = list(c = factor(NA)), silent = TRUE)
  tmb_opt3a <- opt_fun(tmb_obj3a)
  
  # Create 2 parameter with maternal data objective function
  tmb_obj2a = MakeADFun(data, list(c = prop_pregnant, m50 = pars$m50, m95 = pars$m95), 
                        DLL = "logistic2maternal", silent = TRUE)
  tmb_opt2a <- opt_fun(tmb_obj2a)
  
  print(paste0("Simulation ", id, ". Iteration ", iteration, " complete!"))
  
  results <- tibble(data = list(data), method = c("2", "2a", "3", "3a"), 
                    obj = list(tmb_obj2, tmb_obj2a, tmb_obj3, tmb_obj3a),
                    opt = list(tmb_opt2, tmb_opt2a, tmb_opt3, tmb_opt3a)) %>% 
    mutate(prop_pregnant = prop_pregnant) %>% 
    # Check for convergence
    mutate(convergence = map_dbl(opt, ~ .$convergence)) %>% 
    # Run confidence intervals
    mutate(boot_ci = pmap(., resampler, id = id, iteration = iteration)) %>%
    # Tidy up data
    mutate(par_est = map(opt, ~ .$par[order(names(.$par))])) %>%
    mutate(par_name = map(opt, ~ .$par %>% names() %>% sort())) %>%
    mutate(lower = map_if(boot_ci, is_tibble, ~ .$lower)) %>%
    mutate(upper = map_if(boot_ci, is_tibble, ~ .$upper)) %>%
    # Check this, assume it can be removed to save space but untested
    select(-boot_ci, -opt, -obj) %>% 
    mutate(pars = map(method, ~pars)) %>%
    mutate(par_est = pmap(., util_func_1)) %>%
    mutate(par_name = pmap(., util_func_2)) %>%
    mutate(lower = pmap(., util_func_3)) %>%
    mutate(upper = pmap(., util_func_4)) %>%
    # Calculate Reproductive Output
    mutate(r_0_est = pmap_dbl(., util_func_5)) %>%
    # Tidy up data
    mutate(par_est = map2(par_est, r_0_est, ~ c(.x, .y))) %>%
    mutate(par_name = map(par_name, ~ c(., "r_0"))) %>%
    mutate(lower = map(lower, ~ c(., NA))) %>%
    mutate(upper = map(upper, ~ c(., NA))) %>%
    unnest(c(par_name, par_est, lower, upper)) %>%
    unnest(c(par_est, lower, upper)) %>%
    select(-r_0_est)
   
  return(results)
  }else{return(NA)}} 
 
  else{return(NA)}
  
}

## Plot simulation (first 6)
plot_sim_data <- function(id, data){
  df <- data %>% filter(!is.na(par_est))
  
  dir_out <- paste0(here(), "/data/simulations/", id)
  if (!dir.exists(dir_out)) {
    dir.create(dir_out)
  }
  
  pars <- df[1,]
  
  len_hist <- data %>% select(iteration, data) %>% group_by(iteration) %>%
    slice(1) %>% ungroup() %>% slice(1:6) %>% unnest(data)
    
  
  # Plot 1. Simulated population length structure
  (ggplot(data = len_hist) + geom_histogram(aes(x = x, fill = as.factor(z)), col = "black") + 
      facet_wrap(~iteration) +
    labs(title = "Simulated datasets - iterations 1 to 6", x = "Length (mm)", y = "Count")) %>% 
    ggsave(paste(dir_out, "p1.png", sep = "/"), ., width = 8, height = 6)  
    

  # Plot 3. True underlying length at age and relative selectivity
  data3 <- tibble(a = seq(0:pars$AMAX)) %>% 
    mutate(len = pars$Linf *(1 - exp(-pars$K *(a - pars$t0)))) %>% 
    mutate(mat = (1 / (1 + exp(-log(19) * ((len  - pars$l50) / (pars$l95 - pars$l50)))))) %>% 
    mutate(matern = (pars$c / (1 + exp(-log(19) * ((len  - pars$m50) / (pars$m95 - pars$m50)))))) %>% 
    mutate(eta = -0.5 * (pars$theta1 * pars$mesh - sqrt(pars$theta1^2 * pars$mesh^2 + 4 * pars$theta2))) %>%
    mutate(rho = (pars$theta1 * pars$mesh) / eta) %>% 
    mutate(sel  = (len / (rho * eta))^rho * exp(rho - len / eta))
    
  (ggplot() + 
    geom_line(data = data3, aes(x = a, y = len)) + 
    geom_line(data = data3, aes(x = a, y = sel * pars$Linf), col = "red") + 
    labs(title = "Sampled length at age data", x = "Age (years)", y = "Length (mm)")) %>% 
    ggsave(paste(dir_out, "p2.png", sep = "/"), ., width = 6, height = 4)  
  
  # Plot 3. True underlying length at maturity, maternity & selectivity
  (ggplot() + 
    geom_line(data = data3, aes(x = len, y = mat), linetype = "dashed") +
    geom_line(data = data3, aes(x = len, y = matern), linetype = "dotted") +
    geom_line(data = data3, aes(x = len, y = sel), col = "red") +
    labs(title = "True underlying length at maturity, maternity & selectivity", 
         x = "Length (mm)", y = "Proportion")) %>% 
    ggsave(paste(dir_out, "p3.png", sep = "/"), ., width = 6, height = 4)  

  data4 <- df %>% 
    select(iteration, method, par_est, par_name, xmin, xmax) %>% 
    spread(key = par_name, value = par_est) %>% 
    mutate(x = map2(xmin, xmax, ~ seq(.x, .y, 5))) %>% 
    mutate(y = pmap(., ~ (..5) / (1 + exp(-log(19) * (((..9) - (..6)) / (..7 - ..6)))))) %>% 
    select(iteration, method, x, y) %>% 
    unnest(c(x,y))
  
  data5 <- data3 %>% 
    mutate(method = "2") %>% 
    group_by(method) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(method = map(method, ~ as.list(unique(data4$method)))) %>% 
    unnest(method) %>% 
    unnest(method) %>% 
    unnest(data)
  
  (ggplot() + 
    geom_line(data = data4, aes(x = x, y = y, group = iteration), alpha = 0.1) + 
    geom_line(data = data5, aes(x = len, y = matern), col = "red") +
    facet_wrap(~method) + ylim(0, 1) +
    labs(title = "Simulated vs true length at maternity", 
         x = "Length (mm)", y = "Proportion")) %>% 
    ggsave(paste(dir_out, "p4.png", sep = "/"), ., width = 8, height = 6)  
    
}

util_func_1 <- function(method, par_est, prop_pregnant, convergence, ...) {

  if (method %in% c("3a", "2")) {
    par_est <- as.list(par_est)
    par_est$c <- prop_pregnant
  } else {
    par_est <- as.list(par_est)
    if(method %in% "2a"){
      par_est$c = 1
    }
    
  }

  
  names(par_est) <- ifelse(names(par_est) == "l50", "m50", names(par_est))
  names(par_est) <- ifelse(names(par_est) == "l95", "m95", names(par_est))
if(convergence == 0){
  return(par_est)}else{return(list(m50 = NA, m95 = NA, c = NA))}
}

util_func_2 <- function(method, par_name, convergence, ...) {
  if (!method %in% "3") {
    par_name <- c(par_name, "c")
  }
  par_name <- ifelse(par_name %in% "l50", "m50", par_name)
  par_name <- ifelse(par_name %in% "l95", "m95", par_name)
  
  if(convergence == 1){
    par_name <- c("m50", "m95", "c")}
  return(par_name)
}

util_func_3 <- function(method, lower, convergence, ...) {
  if (!method == "3") {
    lower <- as.list(c(lower, NA))
  } else {
    lower <- as.list(lower)
  }
  if(convergence == 1){
    lower <- list(NA, NA, NA)}
    return(lower)
}

util_func_4 <- function(method, upper, convergence, ...) {
  if (!method == "3") {
      upper <- as.list(c(upper, NA))
  } else {
    upper <- as.list(upper)
  }
  if(convergence == 1){
    upper <- list(NA, NA, NA)}
  return(upper)
}

util_func_5 <- function(pars, par_est, convergence, ...) {
  if(convergence == 0){
  integrate(r_0, 0, Inf, df = pars, df2 = par_est)$value} else
  {NA_real_}
}


# Need to have a look at this, in particular whether non-convergent bootstraps
# should be filtered out... 

resampler <- function(data, convergence, obj, method, prop_pregnant, id, iteration, ...) {
  if (convergence == 0) {
    output <- bootstraps(data, 250) %>%
      mutate(n_pregnant = map_dbl(splits, ~ sum(analysis(.)$z))) %>%
      filter(n_pregnant > 0) %>%
      mutate(opt = map(splits, boot_helper,
        obj = obj,
        prop_pregnant = prop_pregnant, method = method
      )) %>%
      mutate(convergence = map_dbl(opt, ~ .$convergence)) %>%
      filter(convergence == 0) %>%
      mutate(pars = map(opt, ~ .$par)) %>%
      mutate(par_name = map(opt, ~ .$par %>% names())) %>%
      unnest(c(pars, par_name)) %>%
      group_by(par_name) %>%
      summarise(
        lower = quantile(pars, 0.5 / 2),
        upper = quantile(pars, 1 - 0.5 / 2)
      )

    print(paste0("Boostrapping complete for ", id, ". Iteration ", iteration, "- method", method))

    return(output)
  }

  else {
    print(paste0("Non convergence for ", id, ". Iteration ", iteration, "- method", method))
    return(NA)
  }
}

boot_helper <- function(splits, obj, method, prop_pregnant) {
  dll <- case_when(
    method %in% "2" ~ "logistic2",
    method %in% "2a" ~ "logistic2maternal",
    method %in% "3" ~ "logistic3",
    method %in% "3a" ~ "logistic3"
  )
  par_map <- if (method %in% "3a") {
    list(c = factor(NA))
  } else {
    NULL
  }
  par_list <- as.list(obj$par)
  if (method %in% "3a") {
    par_list$c <- prop_pregnant
  }

  obj <- MakeADFun(data = analysis(splits), parameters = par_list, DLL = dll, map = par_map, silent = TRUE)
  opt <- opt_fun(obj)
  return(opt)
}

# Reproductive output
r_0 <- function(a, df, df2 = NULL) {

  # Change parameters to estimated / fixed values
  if (!is.null(df2)) {
    names(df2) <- ifelse(names(df2) == "l50", "m50", names(df2))
    names(df2) <- ifelse(names(df2) == "l95", "m95", names(df2))
    df$m50 <- df2$m50
    df$m95 <- df2$m95
    df$c <- df2$c
  }


  # Length at age function
  len <- function(a) {
    df$Linf * (1 - exp(-df$K * (a - df$t0)))
  }

  # Maternity at age function
  maternity <- function(a) {
    (df$c / (1 + exp(-log(19) * ((len(a) - df$m50) / (df$m95 - df$m50)))))
  }

  # Fecundity at age function
  fecundity <- function(a) {
    if (df$species == "School shark") {
      df$f + df$g * len(a)
    } else {
      df$f * exp(df$g * len(a))
    }
  }

  # Birth function
  birth <- function(a) {
    if (df$species == "School shark") {
      ifelse(a >= 11.2404, fecundity(a) * maternity(a) * (1 / 2), 0)
    } else {
      ifelse(a >= 3.9694, fecundity(a) * maternity(a) * (1 / 2), 0)
    }
  }

  # Natural mortality
  mortality <- function(a) {
    rep(-df$M, length(a))
    
    
  }

  # Sum to optimize
  birth(a) * exp(sapply(a, function(x) {
    integrate(mortality, 0, x)$value
  }))
}

# Did any parameters hit boundaries?
par_check <- function(df){
  any(ifelse(filter(df, par_name %in% c("m50", "m95"))$par_est %in% c(0, 10000), TRUE, FALSE))
}


TMBAIC=function(opt, p=2, n=Inf){
  k = length(opt[["par"]])
  if( all(c("par","objective") %in% names(opt)) ) negloglike = opt[["objective"]]
  if( all(c("par","value") %in% names(opt)) ) negloglike = opt[["value"]]
  Return = p*k + 2*negloglike + 2*k*(k+1)/(n-k-1)
  return( Return )
}


# Revised functions


# Length-transition matrix
len_transition_matrix <- function(nlencl, mid, df, upper, lower){
  
  LTM <- matrix(ncol = nlencl, nrow = nlencl) # to, from
  for (ii in 1:nlencl) {  # starting length class
    
    l_mid <- mid[ii]
    end_len = l_mid + (df$Linf - l_mid) * 
      (1 - exp(-df$K)) # expected 
    
    tempsum = 0
    for (i in 1:nlencl) { # ending length class
      
      temp_sd <- end_len * df$cv_l
      
      if (i == 1) {
        # x+1 is upper bound of current 1 mm length class
        LTM[i,ii] = pnorm(upper[i], mean = end_len,
                          sd = temp_sd, lower.tail = T)
      } # if
      else {
        # upper bound - lower bound
        temp_upper = pnorm(upper[i], mean = end_len,
                           sd = temp_sd, lower.tail = T)
        temp_lower = pnorm(lower[i], mean = end_len,
                           sd = temp_sd, lower.tail = T)
        LTM[i,ii] = temp_upper - temp_lower
        
      } # else
      
      tempsum = tempsum + LTM[i,ii]
      
    } #i
    
    # ensure each row of the LTM sums to 1
    for (i in 1:nlencl) { # ending length class
      LTM[i,ii] = LTM[i,ii] / tempsum
    } #i
  } #ii
  
  return(LTM)
}


generate_data <- function(df, id) {
  
  # Convert selectivity parameters
  df <- df %>%
    mutate(eta = -0.5 * (theta1 * mesh - sqrt(theta1^2 * mesh^2 + 4 * theta2))) %>%
    mutate(rho = (theta1 * mesh) / eta)
  
  # Mean length at age
  a = 0:df$AMAX
  len <- df$Linf * (1 - exp(-df$K * (a - df$t0)))
  
  # Length and sd at recruitment (Age 0)
  l_rec <- len[1]
  sd_l_rec = l_rec * df$cv_l
  
  # initialize length structure of model and length bins
  l_max = 2500
  l_inc = 10
  lower = seq(0, l_max - l_inc, l_inc)
  upper = lower + l_inc
  mid = lower + (l_inc/2)
  nlencl = length(mid)
  
  # Length distribution of recruits
  l_rec_dist = pnorm(upper, mean = l_rec, sd = sd_l_rec, lower.tail = T) -
    pnorm(lower, mean = l_rec, sd = sd_l_rec, lower.tail = T)
  l_rec_dist = l_rec_dist / sum(l_rec_dist)
  
  # Gillnet selectivity
  sel = (mid / (df$rho * df$eta))^df$rho * exp(df$rho - mid / df$eta)
  sel = sel / max(sel)
  
  # Generate a length transition matrix
  LTM <- len_transition_matrix(nlencl, mid, df, upper, lower)
  
  # Alex's per-recruit model
  
  # Pre recruit numbers surviving after natural mortality
  f_n_per_rec <- rep(0, nlencl)
  f_n_per_rec_a <- matrix(0L, nrow = (df$AMAX + 1), ncol = nlencl)
  f_s_per_rec_a <- catch <- f_n_per_rec_a
  catch_len = rep(0, nlencl)
  
  # Get length distribution of recruits
  f_n_per_rec_a[1,] = l_rec_dist
  f_n_per_rec = f_n_per_rec + f_n_per_rec_a[1,]
  
  # Fishing mortality
  f_len= sel * df$F_mort
  
  # Total mortality experienced by the stock, due to fishing
  z_len = df$M + f_len
  
  # calculate catch at length for timestep (in numbers)
  catch[1,] = f_n_per_rec_a[1,] * (f_len / z_len) * (1 - exp(-z_len))  
  catch_len = catch_len + catch[1,]
  
  # apply mortality to calculate survival
  for (t in 2:(df$AMAX + 1)) {
    if (t < (df$AMAX + 1)) {
      f_s_per_rec_a[t,] = f_n_per_rec_a[t-1,] * exp(-z_len)
    } else {
      f_s_per_rec_a[t,] = f_n_per_rec_a[t-1,] * exp(-z_len) / 
        (1 - exp(-z_len))
    }
    
    # apply growth - matrix multiplication - faster
    M1=t(f_s_per_rec_a[t,])
    M2=t(LTM[,])
    f_n_per_rec_a[t, ] = as.numeric(M1 %*% M2)
    f_n_per_rec = f_n_per_rec + f_n_per_rec_a[t,]
    
    # calculate catch at length for timestep (in numbers)
    catch[t,] = f_n_per_rec_a[t,] * (f_len / z_len) * (1 - exp(-z_len))  
    
    catch_len = catch_len + catch[t,]
    
  }
  
  
  if (!is.numeric(sum(catch))) {
    cat("Problem - OperatingModel. sum(Catch) not numeric")
  }
  
  # Results
  pr_results <- list(catch = catch,
                  catch_len = catch_len,
                  f_n_per_rec = f_n_per_rec)
  
  # Generate required catch at length
  data <- tibble(x = mid) %>%
    slice_sample(n = df$Nsamples, weight_by = pr_results$catch_len, replace = TRUE) %>% 
    mutate(y = NA)
  
  print(id)
  return(data)
}

