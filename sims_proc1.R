##----clear workspace----

rm(list=ls())

##----libraries/sources----

library(tidyverse)
library(qvalue)       #qvalue()
library(ROCR)         #prediction(), performance()
#library(pwrutilities) #myfrac(), here added manually
#library(BioNet)       #fitBumModel(), fits a beta-uniform mixture model
library(furrr)

options(future.globals.maxSize = 1000 * 1024^2) # to set maximum allowed size to 1GB
plan(multicore, workers = 90)

#' Print the fraction of TRUEs in a vector of logicals
#'
#' @param vec Vector of logicals
#'
#' @return Fraction of TRUE values in vec
#' @export
#'
#' @examples
#' set.seed(123)
#' myfrac(rnorm(10) > 0)
#' myfrac(TRUE)
#' myfrac(FALSE)
myfrac <- function (vec) {
  tally <- table(vec)
  if (length(tally) == 2L)
    return(tally["TRUE"]/sum(tally))
  if (names(tally) == "FALSE")
    return(stats::setNames(0, "TRUE"))
  stats::setNames(1, "TRUE")
}

##----load simulation data----

load(file = "sims.rdata")

noise <- attributes(sims_sine_rob)$global_pars$noise_level

##----functions to calculate pi0----

calculate_pi0_boot_cen  <- function (cen_pvals) {
  pi0_boot_cen <- try(qvalue::pi0est(cen_pvals, pi0.method = "bootstrap")$pi0, silent = TRUE)
  if (class(pi0_boot_cen) == "try-error") {
    return(NA)
  }
  pi0_boot_cen
}

calculate_pi0_boot_nonc <- function (nonc_pvals) {
  pi0_boot_nonc <- try(qvalue::pi0est(nonc_pvals, pi0.method = "bootstrap")$pi0, silent = TRUE)
  if (class(pi0_boot_nonc) == "try-error") {
    return(NA)
  }
  pi0_boot_nonc
}

calculate_pi0_boot_reg_nonc <- function (reg_nonc_pvals) {
  pi0_boot_reg_nonc <- try(qvalue::pi0est(reg_nonc_pvals, pi0.method = "bootstrap")$pi0, silent = TRUE)
  if (class(pi0_boot_reg_nonc) == "try-error") {
    return(NA)
  }
  pi0_boot_reg_nonc
}

#calculate_pi0_bum_reg_nonc <- function (reg_nonc_pvals) {
#  names(reg_nonc_pvals) <- 1:length(reg_nonc_pvals)
#  bum <- suppressWarnings(fitBumModel(reg_nonc_pvals, plot=FALSE))
#  if (class(bum) == "character") {
#    pi0_bum_reg_nonc <- NA
#  } else {
#    pi0_bum_reg_nonc <- piUpper(bum)
#  }
#  pi0_bum_reg_nonc
#}

calculate_pi0_hreg_est_nonc <- function (parcomb_list_el, parcomb_name) {
  params <- tibble(parcomb_name = parcomb_name) %>% 
    tidyr::separate(col = "parcomb_name",
                    into = c("s_int", "frac_arr", "amp_sdlog",
                             "delta", "outliers", "robust"),
                    sep = "_")
  delta <- as.numeric(params$delta)
  amp_over_sigma_est <- parcomb_list_el$amp_est/parcomb_list_el$sigma_est
  unname(myfrac(amp_over_sigma_est >= delta/noise))
}

##----functions to calculate ROC & PR curve & AUC----

calculate_roc <- function (nonc_pvals, theanswers) {
  pred <- prediction(1 - nonc_pvals, theanswers)
  performance(pred, measure = "tpr", x.measure = "fpr")
}

calculate_pr <- function (nonc_pvals, theanswers) {
  pred <- prediction(1 - nonc_pvals, theanswers)
  performance(pred, measure = "tpr", x.measure = "prec")
}

calculate_auc <- function (nonc_pvals, theanswers) {
  pred <- prediction(1 - nonc_pvals, theanswers)
  performance(pred, measure = "auc")
}

calculate_aFDR <- function (pi0_boot_reg_nonc, 
                            parcomb_list_el, maxFDR) {
  if (is.na(pi0_boot_reg_nonc)) {
    return(NA)
  }
  nonc_pvals <- parcomb_list_el$nonc_pvals
  theanswers <- parcomb_list_el$theanswers
  qvals          <- qvalue(nonc_pvals, pi0 = pi0_boot_reg_nonc)$qvalues
  theanswers_sel <- theanswers[qvals < maxFDR] == 0L
  if (length(theanswers_sel) == 0L || all(is.na(theanswers_sel))) {
    aFDR <- NA
  } else {
    aFDR <- unname(myfrac(theanswers_sel))
  }
  aFDR
}

calculate_power <- function (FPR, TPR, frac_really_rhy, maxFDR) {
  df           <- data.frame(FPR,TPR)
  colnames(df) <- c("FPR","TPR")
  curve1       <- function(x) x*(frac_really_rhy)/(1-frac_really_rhy)*(1 - maxFDR)/maxFDR
  curve2       <- approxfun(df$FPR, df$TPR, rule = 2)
  point_x      <- try(uniroot(function(x) curve1(x) - curve2(x), c(0,1))$root, silent = TRUE)
  if (class(point_x) == "try-error") {
    point_x <- NA
    point_y <- NA
  } else {
    point_y <- curve1(point_x)
  }
  list(x = point_x, y = point_y)
}


##----applying functions to calculate pi0----

pi0_boot_cen      <- sims_sine_rob %>%
  future_map_dbl(~ calculate_pi0_boot_cen(.x$cen_pvals))

pi0_boot_nonc     <- sims_sine_rob %>%
  future_map_dbl(~ calculate_pi0_boot_nonc(.x$nonc_pvals))

pi0_boot_reg_nonc <- sims_sine_rob %>%
  future_map_dbl(~ calculate_pi0_boot_reg_nonc(.x$reg_nonc_pvals))

#pi0_bum_reg_nonc  <- sims_sine_rob %>%
#  future_map_dbl(~ calculate_pi0_bum_reg_nonc(.x$reg_nonc_pvals))

pi0_hreg_est_nonc <- sims_sine_rob %>% 
  imap(calculate_pi0_hreg_est_nonc)


##----applying functions to calculate roc, pr and auc----

FPR <- sims_sine_rob %>%
  future_map(~ calculate_roc(.x$nonc_pvals, .x$theanswers)@x.values[[1]])

TPR <- sims_sine_rob %>%
  future_map(~ calculate_roc(.x$nonc_pvals, .x$theanswers)@y.values[[1]])

prec <- sims_sine_rob %>%
  future_map(~ calculate_pr(.x$nonc_pvals, .x$theanswers)@x.values[[1]])

AUC <- sims_sine_rob %>%
  future_map_dbl(~ calculate_auc(.x$nonc_pvals, .x$theanswers)@y.values[[1]])

frac_really_rhy   <- sims_sine_rob %>%
  future_map_dbl(~ .x$frac_really_rhy)

amp_est <- sims_sine_rob %>%
  future_map(~ .x$amp_est)

sigma_est <- sims_sine_rob %>%
  future_map(~ .x$sigma_est)

aFDR_005 <- future_map2(pi0_boot_reg_nonc, 
                    sims_sine_rob,
                    calculate_aFDR, 0.05)

aFDR_010 <- future_map2(pi0_boot_reg_nonc, 
                        sims_sine_rob,
                        calculate_aFDR, 0.10)

aFDR_015 <- future_map2(pi0_boot_reg_nonc, 
                        sims_sine_rob,
                        calculate_aFDR, 0.15)

##----FPR, TPR need to be calcluated first----
##----then applying function to calculate power = y-value----

df <- purrr::transpose(list(FPR = FPR,
                            TPR = TPR,
                            frac_really_rhy = frac_really_rhy))

maxPower_005 <- df %>%
  future_map_dbl(~ calculate_power(.x$FPR, .x$TPR, .x$frac_really_rhy, 0.05)$y)
maxPower_010 <- df %>%
  future_map_dbl(~ calculate_power(.x$FPR, .x$TPR, .x$frac_really_rhy, 0.10)$y)
maxPower_015 <- df %>%
  future_map_dbl(~ calculate_power(.x$FPR, .x$TPR, .x$frac_really_rhy, 0.15)$y)

aPower_005 <- df %>%
  future_map_dbl(~ calculate_power(.x$FPR, .x$TPR, .x$frac_really_rhy, aFDR_005)$y)
aPower_010 <- df %>%
  future_map_dbl(~ calculate_power(.x$FPR, .x$TPR, .x$frac_really_rhy, aFDR_010)$y)
aPower_015 <- df %>%
  future_map_dbl(~ calculate_power(.x$FPR, .x$TPR, .x$frac_really_rhy, aFDR_015)$y)

##----combine results----

sims_sine_rob_proc1 <- purrr::transpose(list(
                                 pi0_boot_cen      = pi0_boot_cen,
                                 pi0_boot_nonc     = pi0_boot_nonc,
                                 pi0_boot_reg_nonc = pi0_boot_reg_nonc,
#                                 pi0_bum_reg_nonc  = pi0_bum_reg_nonc,
                                 pi0_hreg_est_nonc = pi0_hreg_est_nonc,
                                 frac_really_rhy   = frac_really_rhy,
                                 FPR       = FPR,
                                 TPR       = TPR,
                                 prec      = prec,
                                 AUC       = AUC,
                                 aFDR_005  = aFDR_005,
                                 aFDR_010  = aFDR_010,
                                 aFDR_015  = aFDR_015,
                                 maxPower_005 = maxPower_005,
                                 maxPower_010 = maxPower_010,
                                 maxPower_015 = maxPower_015,
                                 aPower_005   = aPower_005,
                                 aPower_010   = aPower_010,
                                 aPower_015   = aPower_015,
                                 amp_est      = amp_est,
                                 sigma_est    = sigma_est))

##----save data----

save(sims_sine_rob_proc1,
     file = "sims_proc1.rdata")


##----save R and package versions----

sessionInfo()
