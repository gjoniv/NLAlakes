library(tidyverse)
library(brms)
library(tidybayes)

d_temp = read_csv("ccsr_colonybased_summaries_2022_10_27.csv")

d_temp_c = d_temp %>% select(uid, site_id, temp_mean) %>% 
  mutate(temp_c = scale(temp_mean))

d_cell_based = read_csv("ccsr_cellbased_summaries_2020_11_09.csv") %>% 
  select(-temp_mean) %>% 
  left_join(d_temp_c) %>% 
  mutate(zoo_total_biomass = zoo_mean_biomass*zoo_total_density)

d_cell_based_centered = d_cell_based %>% 
  mutate(log_nit = log(nitrate_mean),
         log_phos = log(phos_total_mean),
         log_zoo = log(zoo_total_biomass),
         log_biovol = log(phyto_mean_biovol),
         log_phyto = log(phyto_total_density)) %>% 
  mutate(nit_c = scale(log_nit),
         pho_c = scale(log_phos),
         zoo_c = scale(log_zoo),
         biovol_c = scale(log_biovol))

mod_a = brm(log_phyto ~ biovol_c*temp_c,
            family = gaussian(),
            data = d_cell_based_centered,
            prior = c(prior(normal(0, 1), class = "b"),
                      prior(normal(-0.875, 0.2), coef = "biovol_c"),
                      prior(normal(10, 2), class = "Intercept"),
                      prior(exponential(1), class = "sigma")))

mod_b = update(mod_a, formula = .~biovol_c*temp_c*pho_c, newdata = d_cell_based_centered)
mod_c = update(mod_a, formula = .~biovol_c*temp_c*zoo_c, newdata = d_cell_based_centered)
mod_d = update(mod_a, formula = .~biovol_c*temp_c*pho_c*zoo_c, newdata = d_cell_based_centered)

saveRDS(mod_a, file = "models/mod_a.rds")
saveRDS(mod_b, file = "models/mod_b.rds")
saveRDS(mod_c, file = "models/mod_c.rds")
saveRDS(mod_d, file = "models/mod_d.rds")

waic(mod_a)
waic(mod_b)
waic(mod_c)
waic(mod_d)


plot(conditional_effects(mod_d, effects = "biovol_c"), points = T)

