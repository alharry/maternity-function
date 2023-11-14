

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
