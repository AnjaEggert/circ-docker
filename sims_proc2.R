##----clear workspace----

rm(list=ls())


##----libraries/sources----

library(tidyverse)

##----further processing of simulation results on Linux server------------##

load(file = "sims_proc1.rdata")

t   <- enframe(sims_sine_rob_proc1)

tt  <- t %>% unnest_wider(value)

ttt <- tt %>%

# exclude all columns with nested lists
dplyr::select(-FPR, -TPR, -prec, -amp_est, -sigma_est) %>%
tidyr::separate(col = "name",
                  into = c("s_int", "frac_arr", "amp_sdlog",
                           "delta", "outliers", "robust"),
                  sep = "_")

sims_sine_rob_proc2 <- ttt %>% select(-robust) %>% mutate_if(is.character, as.numeric) %>%
  mutate(s_int = as.integer(round(s_int)))

save(sims_sine_rob_proc2, file="sims_proc2.rdata")

