break_in_mgi <- function (df, theid="entrez_id",
                          thefunc=e2me) {
  df %>% left_join(thefunc(unique(pull(df, !! as.name(theid)))), by=theid) %>%
    #dplyr::select(-entrez_id) %>%
    dplyr::select(mgi_symbol, everything())
}


myfrac <- function (vec) {
  tally <- table(vec)
  if (length(tally) == 2L)
    return(tally["TRUE"]/sum(tally))
  if (names(tally) == "FALSE")
    return(stats::setNames(0, "TRUE"))
  stats::setNames(1, "TRUE")
}


plot_many_fit <- function (mgi_symbols, tissue = "liver",
                           df_ave = zhang_long_pval,
                           facet_ncol=NULL, thetitle=NULL, strict_scale=TRUE,
                           pval_cutoff = 0.1,
                           ylab_str="Affy expression") {

  # mutate(
  #   tcols=ifelse(!! circ_do, "red", "black"),
  #   tface=ifelse(!! de, "bold.italic", "italic"))

  df_ave %>% dplyr::filter(tissue == !!tissue) %>%
    dplyr::filter(mgi_symbol %in% mgi_symbols) %>%
    ## respect the ordering of incoming mgi symbols
    mutate(mgi_symbol=factor(mgi_symbol, levels=mgi_symbols)) %>%

    ggplot(aes(ct, affy_exp)) +


    stat_hrfit(aes(alpha = ifelse(qval_hreg < pval_cutoff &
                                    !is.na(qval_hreg),
                                  1 - qval_hreg/pval_cutoff,
                                  0)),
               geom="line", size=0.75) +

    geom_point(color="black", fill = "dodgerblue", shape=21) +


    my_x(thebreaks=my_circbreaks(12)) +

    ## if strict, begin all y axes at 0!
    scale_y_continuous(limits=c(ifelse(strict_scale, 0, NA), NA)) +
    scale_alpha_identity() +

    ## annotations
    xlab("CT (hrs)") + ylab(ylab_str) +
    ggtitle(thetitle) +

    ## faceting
    facet_wrap(~ mgi_symbol, scales="free_y", ncol=facet_ncol) +

    ## theming
    pw_theme + pw_theme_lattice +
    theme(strip.text = element_text(face = "italic"),
          panel.grid.major.x = element_line(color = "grey90"),
          strip.background = element_rect(fill = "gray93", color = NA),
          panel.spacing.y = unit(1, "mm"))

}


plot_regexp_many_fit <- function (regexp, df_ave = zhang_long_pval, ...) {
  df_ave$mgi_symbol %>%
    str_subset(regexp) %>% unique %>%
    plot_many_fit(df_ave = df_ave, ...)
}


clean_tissue_str_helper <- function (df) {
  df %>%
    mutate(tissue = str_replace(tissue, "_", " ")) %>%
    mutate(tissue = str_to_sentence(tissue)) %>%
    ## respect the ordering of incoming
    mutate(tissue = factor(tissue, levels = unique(tissue)))
}


plot_tissue_many_fit <- function (mgi_symbol,
                                  tissues = "all",
                                  df_ave = zhang_long_pval,
                                  facet_ncol = NULL, thetitle = NULL,
                                  strict_scale = TRUE,
                                  pval_cutoff = 0.1,
                                  ylab_str = "Affy expression") {

  if (tissues == "all") {
    tissues = unique(df_ave$tissue)
  }

  # mutate(
  #   tcols = ifelse(!! circ_do, "red", "black"),
  #   tface = ifelse(!! de, "bold.italic", "italic"))

  df_ave %>% dplyr::filter(mgi_symbol == !!mgi_symbol) %>%
    dplyr::filter(tissue %in% !!tissues) %>%
    clean_tissue_str_helper %>%

    ggplot(aes(ct, affy_exp)) +


    stat_hrfit(aes(alpha = ifelse(qval_hreg < pval_cutoff &
                                    !is.na(qval_hreg),
                                  1 - qval_hreg/pval_cutoff,
                                  0)),
               geom = "line", size = 0.75) +

    geom_point(color = "black", fill = "dodgerblue", shape = 21) +


    my_x(thebreaks = my_circbreaks(12)) +

    ## if strict, begin all y axes at 0!
    scale_y_continuous(limits = c(ifelse(strict_scale, 0, NA), NA)) +
    scale_alpha_identity() +

    ## annotations
    xlab("CT (hrs)") + ylab(ylab_str) +
    ggtitle(thetitle) +

    ## faceting
    facet_wrap(~ tissue, scales = "free_y", ncol = facet_ncol) +

    ## theming
    pw_theme + pw_theme_lattice +
    theme(strip.text = element_text(face = "italic"),
          panel.grid.major.x = element_line(color = "grey90"),
          strip.background = element_rect(fill = "gray93", color = NA),
          panel.spacing.y = unit(1, "mm"))
}

# redistribute_pvals <- function (pvals, ref_low = 0.5, ref_high = 0.6,
#                                 breakpoint = 0.9) {
#   candidate_ind <- which(pvals > breakpoint)
#   n_unif_reference <-
#     round(length(which(pvals > ref_low & pvals < ref_high)) *
#             (1 - breakpoint)/(ref_high - ref_low))
#   if (length(candidate_ind) <= n_unif_reference)
#     stop("redistribute_pvals(): No right tail peak")
#   to_replace <- sample(candidate_ind,
#                        length(candidate_ind) - n_unif_reference,
#                        replace = FALSE)
#   replace(pvals, to_replace, runif(length(to_replace)))
# }
#
# regularize_pvals <- function (pvals, breakpoints = c(0.8, 0.95)) {
#   while (length(breakpoints) >= 1L) {
#     breakpoint <- head(breakpoints, 1)
#     breakpoints <- tail(breakpoints, -1)
#     pvals_tent <- try(redistribute_pvals(pvals, breakpoint = breakpoint),
#                       silent = TRUE)
#     if (class(pvals_tent) == "try-error")
#       return(pvals)
#     else
#       pvals <- pvals_tent
#   }
#   pvals
# }

redistribute_pvals <- function (pvals, breakpoint = 0.9, next_breakpoint = 0.8,
                                ref_low = 0.6, ref_high = 0.7) {
  candidate_ind <- which(pvals > breakpoint)
  n_unif_reference <-
    round(length(which(pvals > ref_low & pvals < ref_high)) *
            (1 - breakpoint)/(ref_high - ref_low))
  if (length(candidate_ind) <= n_unif_reference)
    stop("redistribute_pvals(): No right tail peak")
  to_replace <- sample(candidate_ind,
                       length(candidate_ind) - n_unif_reference,
                       replace = FALSE)
  replace(pvals, to_replace, runif(length(to_replace),
                                   min = next_breakpoint, max = 1))
}

regularize_pvals <- function (pvals, breakpoints = c(0.99, 0.98, 0.9, 0.8),
                              ...) {
  while (length(breakpoints) >= 1L) {
    breakpoint <- head(breakpoints, 1)
    breakpoints <- tail(breakpoints, -1)
    if (length(breakpoints) == 0L)
      next_breakpoint = 0
    else
      next_breakpoint = head(breakpoints, 1)
    pvals_tent <- try(redistribute_pvals(pvals, breakpoint = breakpoint,
                                         next_breakpoint = next_breakpoint,
                                         ...),
                      silent = TRUE)
    if (class(pvals_tent) == "try-error")
      return(pvals)
    else
      pvals <- pvals_tent
  }
  pvals
}

simple_pi0 <- function (pvals, upper = 0.6, lower = 0.5, tail = 0.9) {
  (round(length(pvals[pvals > lower & pvals < upper])*
           tail/(upper - lower)) +
     length(pvals[pvals >= tail]))/length(pvals)
}


meta_fun <- function (pvals, fun, strict = TRUE, ...) {
  extra_args <- list(...)
  if ("r" %in% names(extra_args) && extra_args[["r"]] > length(pvals)) {
    if (strict) return(NA)
    extra_args[["r"]] <- length(pvals)
  }
  ifelse(length(pvals) == 1L, pvals,
         do.call(fun, c(list(pvals), extra_args))$p)
}


get_pi0_reg <- function (pvals, reg_seq = c(0.99, 0.98, 0.9, 0.8)) {
  pi0est(regularize_pvals(pvals[!is.na(pvals)],
                          breakpoints = reg_seq),
         pi0.method = "bootstrap")$pi0
}


