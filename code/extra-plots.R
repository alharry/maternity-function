data <- read_rds(here("data", "simulations-summarised.Rds")) %>% 
  #filter(!(method_long %in% "2PLF - maternity" & repro_freq %in% "1 (Annual)")) %>% 
  filter(method_long %in% c("3PLF - estimated", "3PLF - fixed")) %>% 
  mutate(Nsamples_sel = fct_cross(factor(mesh_name, levels = c("low", "high")), as.factor(Nsamples))) %>% 
  mutate(n_maternal_cut = cut(n_maternal, breaks = seq(0, 1600, 20))) %>% 
  filter(par_name %in% "r_0") %>% group_by(species, Nsamples_sel, repro_freq) %>% 
  slice_min(abs_error)

data %>% 
  ggplot() + 
  geom_tile(aes(x = as.factor(Nsamples_sel), y = repro_freq, fill = method_long)) + 
  #scale_fill_brewer(palette = "Spectral") + 
  scale_fill_viridis_d() + 
  guides(fill = guide_legend(title = "Method")) + 
  facet_wrap(~species, nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Selectivity: Sample size", y = bquote(atop(italic(P)[Max],("Reproductive frequency")))) 



(p <- ggplot(data = data, aes(x = n_maternal, y = as.numeric(method_long))) + 
  geom_point() + 
  geom_smooth() + 
  scale_x_continuous(breaks = seq(0,1600,200)) + 
    facet_wrap(~species))


(p <- ggplot(data = data, aes(x = n_maternal, y = abs_error, col = repro_freq)) + 
  geom_point() + 
  scale_color_brewer(palette = "Spectral") +
  scale_y_continuous(breaks = seq(0, 90, 10)) +
  facet_grid( ~ species) + 
  guides(colour = guide_legend(bquote(atop(italic(P)[Max],("Reproductive frequency"))))) + 
  theme_bw() + 
  labs(x = "Median number of maternal females", y = "Mean % absolute relative error") +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center") +
  theme(legend.box = "vertical", legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center"))


data <- read_rds(here("data", "simulations-pivot.Rds")) %>% 
  filter(method_long %in% "3PLF - estimated")

(p <- ggplot() + 
    geom_point(data = data, aes(x = n_maternal, y = abs(rel_error_c), fill = repro_freq), pch =21, alpha = 0.5) + 
    scale_fill_brewer(palette = "Spectral") +
    #scale_y_continuous(breaks = seq(0, 90, 10)) +
    facet_grid( ~ species) + 
    guides(colour = guide_legend(bquote(atop(italic(P)[Max],("Reproductive frequency"))))) + 
    theme_bw() + 
    labs(x = "Median number of maternal females", y = "Mean % absolute relative error") +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center") +
    theme(legend.box = "vertical", legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center"))


data <- read_rds(here("data", "simulations-summarised.Rds")) %>% 
  filter(method_long %in% c("3PLF - estimated")) %>% 
  filter(abs_error <= 10) %>% 
  filter(par_name %in% "c")

(p <- ggplot(data = data, aes(x = n_maternal, y = abs_error, col = repro_freq)) + 
    geom_point() + 
    scale_color_brewer(palette = "Spectral") +
    scale_y_continuous(breaks = seq(0, 90, 10)) +
    scale_x_continuous(breaks = seq(0, 1500, 200)) +
    facet_grid( ~ species) + 
    guides(colour = guide_legend(bquote(atop(italic(P)[Max],("Reproductive frequency"))))) + 
    theme_bw() + 
    labs(x = "Median number of maternal females", y = "Mean % absolute relative error") +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center") +
    theme(legend.box = "vertical", legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center"))
