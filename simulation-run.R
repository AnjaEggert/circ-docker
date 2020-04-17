##----generating rhythmic data----

simulation_run <- function (samp_interv, samp_dur, nsim, frac_arr,
                            amp_meanlog, amp_sdlog, delta,
                            phase_const, phase_rhythm,
                            noise_const, noise_rhythm,
                            outl, robust,
                            comps = 1) {
                            #comps = 2/pi*(-1)^(1:10)/(1:10)) {
  
  ## the time points
  thetime <- seq(0, samp_dur, samp_interv)
  
  ## number of arrhythmic and rhythmic features
  n_arr <- round(nsim*frac_arr)
  n_rhy <- nsim - n_arr
  
  ## generating rhythmic amplitudes
  amps <- rlnorm(n=n_rhy, meanlog = amp_meanlog, sdlog = amp_sdlog)
  
  ## generating data with mult error and outliers
  ts_const  <- generate_fourier_wave_mult_outl(amps      = rep(0, n_arr),
                                               phases = rep(phase_const, n_arr),
                                               thetime   = thetime,
                                               noise     = noise_const,
                                               frac_outl = outl,
                                               comps     = comps)
  
  ts_rhythm <- generate_fourier_wave_mult_outl(amps      = amps,
                                               phases = rep(phase_rhythm,
                                                            n_rhy),
                                               thetime   = thetime,
                                               noise     = noise_rhythm,
                                               frac_outl = outl,
                                               comps     = comps)
  
  ## harmonic regression
  hreg_const  <- harmonic.regression(inputts      = t(ts_const),
                                     inputtime    = thetime,
                                     amp_delta    = delta,
                                     a_over_sigma = delta/noise_const,
                                     robust       = robust)
  
  hreg_rhythm <- harmonic.regression(inputts      = t(ts_rhythm),
                                     inputtime    = thetime,
                                     amp_delta    = delta,
                                     a_over_sigma = delta/noise_rhythm,
                                     robust       = robust)
  
  ## fraction of really rythmic features, given as vector of ordered factors
  ## 1 for constant or amp/noise <= delta/noise_rhythm
  ## 0 for rhythmic or amp/noise > delta/noise_rhythm
  ## for sawtooth simulations use adjusted amplitudes, i.e. 1st harmonic
  
  theanswers <-
    ordered(c(rep(1, n_arr),
              as.integer(amps/noise_rhythm <= delta/noise_rhythm)))
  
  theanswers_saw <-
    ordered(c(rep(1, n_arr),
              as.integer((amps*2)/(pi*noise_rhythm) <= delta/noise_rhythm)))
  
  ## save results
  list(amps                = c(rep(0, n_arr), amps),
       hreg_cen_pvals      = c(hreg_const$pvals,
                               hreg_rhythm$pvals),
       hreg_nonc_pvals     = c(hreg_const$pvals_noncentral,
                               hreg_rhythm$pvals_noncentral),
       amp_est             = c(hreg_const$pars$amp,
                               hreg_rhythm$pars$amp),
       sigma_est           = c(hreg_const$sigma_hat,
                               hreg_rhythm$sigma_hat),
       theanswers          = theanswers,
       frac_really_zero    = n_arr/nsim,
       frac_really_rhy     = unname(myfrac(theanswers == 0)),
       frac_really_rhy_saw = unname(myfrac(theanswers_saw == 0)))
}

##----save R and package versions----

sessionInfo()
