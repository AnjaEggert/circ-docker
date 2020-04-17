##----clear workspace----

rm(list=ls())

##----libraries/sources----

library(tidyverse)
library(rain)
library(ROCR)
library(furrr) # future_map()
library(HarmonicRegression)

source("rhythm-generators.R")
source("helper-functions.R")
source("simulation-run.R")

##----my_experiment() = simulation_run() + regularize_pvals()----

# following argument fixed for this experiment:
# samp_dur     = 47
# phase_const  = 10
# phase_rhythm = 14
# comps        = 1

# following argument default for this experiment:
# nsim         = 20000
# amp_meanlog  = log(0.2)
# noise_const  = 0.2
# noise_rhythm = 0.2

reg_seq <- c(0.99, 0.98, 0.9, 0.8)

my_experiment <- function (pars,
                           samp_dur    = 47,
                           nsim        = 100,
                           amp_meanlog = log(0.2),
                           noise_level = 0.2) {

  sims <- simulation_run(samp_interv  = pars$samp_interv,
                         samp_dur     = samp_dur,
                         nsim         = nsim,
                         frac_arr     = pars$frac_arr,
                         amp_meanlog  = amp_meanlog,
                         amp_sdlog    = pars$amp_sdlog,
                         delta        = pars$delta,
                         phase_const  = 10,
                         phase_rhythm = 14,
                         noise_const  = noise_level,
                         noise_rhythm = noise_level,
                         outl         = pars$outl,
                         #comps       = 1 (is default),
                         #comps       = 2/pi*(-1)^(1:10)/(1:10),
                         robust       = pars$robust
                         )

  nonc_pvals     <- sims$hreg_nonc_pvals
  reg_nonc_pvals <- regularize_pvals(nonc_pvals, breakpoints = reg_seq)
  
  #sigma_est      <- sqrt(sims$ssr/((samp_dur + 1)/pars$samp_interv - 3))

  list(
    frac_really_zero    = sims$frac_really_zero,
    frac_really_rhy     = sims$frac_really_rhy,
    frac_really_rhy_saw = sims$frac_really_rhy_saw,
    cen_pvals           = sims$hreg_cen_pvals,
    nonc_pvals          = nonc_pvals,
    reg_nonc_pvals      = reg_nonc_pvals,
    theanswers          = sims$theanswers,
    amps                = sims$amps,
    amp_est             = sims$amp_est,
    sigma_est           = sims$sigma_est
  )
}

##----defining parameter space to apply function my_experiment()----

parcombs <- cross(list(samp_interv = c(2),
                       frac_arr    = c(0.4, 0.5),
                       amp_sdlog   = c(0.7),
                       delta       = c(0.1),
                       outl        = c(0, 0.02),
                       robust      = TRUE
                       )
                  )

names(parcombs) <- parcombs %>% map(lift(str_c, sep = "_"))

##----run simulation in parallel mode on 90 cores----

plan(multicore, workers = 90)

sims_sine_rob <- parcombs %>%
  future_map(my_experiment, amp_meanlog = log(0.1), noise_level = 0.1)

attr(sims_sine_rob, "global_pars") <- list(amp_meanlog = quote(log(0.1)),
                                              noise_level = 0.1)

##----save data----

save(sims_sine_rob, file = "sims.rdata")

##----save R and package versions----

sessionInfo()
