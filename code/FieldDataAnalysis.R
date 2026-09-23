setwd('D:/R') # set current working directory
rm(list = ls()) # clear the environment
cat("\014") # clear the console
if (!require(pacman)) install.packages("pacman")


# custom function -------------------------------------------------------------------
# Define a piecewise normalization function (optimal at the median, poorer at both ends)
# normalize_seg <- function(x, n) { # n is the optimal median
#   Vmin <- min(x, na.rm = TRUE); Vmax <- max(x, na.rm = TRUE)
#   d1 <- abs(n - Vmin); d2 <- abs(Vmax - n)
#   if (Vmin > n) (Vmax - x) / (Vmax - Vmin) else if (Vmax < n) (x - Vmin) / (Vmax - Vmin)
#   else if (d1 > d2) ifelse(x <= n, (x - Vmin) / (n - Vmin), (n * 2 - x - Vmin) / (n - Vmin))
#   else ifelse(x <= n, (x + Vmax - n * 2) / (Vmax - n), (Vmax - x) / (Vmax - n))}
# All pH values were > 7 in the field experiment, so a simple range normalization method was used directly

# statistical analysis function that performs ANOVA and subsequent multiple comparisons and outputs the results
analyze_variable_tworeg <- function(data, response_var, formula_effects, random_factors) {
  create_base_df <- function() {
    df <- data.frame(Variable = response_var)
    return(df)}
  safe_analyze <- safely(function() {
    data_clean <- data[!is.na(data[[response_var]]), ]
    valid_random_effects <- c()
    for (factor in random_factors) {
      levels <- length(unique(na.omit(data_clean[[factor]])))
      if (levels > 1) {valid_random_effects <- c(valid_random_effects, paste("(1 |", factor, ")"))}}
    if (length(valid_random_effects) == 0) return(NULL)
    random_effects_used <- paste(valid_random_effects, collapse = " + ")
    used_formula <- as.formula(paste(response_var, "~", formula_effects, "+", random_effects_used))
    formula_terms <- attr(terms(used_formula), "term.labels")
    fixed_effects_terms <- formula_terms[!grepl("\\|", formula_terms)]
    random_effects_terms <- formula_terms[grepl("\\|", formula_terms)]
    model_lmer <- suppressMessages(lmer(used_formula, data = data_clean, REML = TRUE, control = lmerControl(optimizer = "Nelder_Mead")))
    residuals_lmer <- residuals(model_lmer)
    base_info <- create_base_df()

    loglik <- logLik(model_lmer)
    mse <- mean(residuals_lmer^2)
    table0 <- cbind(base_info, data.frame(
      UsedFormula = gsub("\\s+", " ", paste(deparse(used_formula), collapse = "")),
      AIC = AIC(model_lmer), BIC = BIC(model_lmer),
      logLik = as.numeric(loglik), npar = attr(loglik, "df"),
      Deviance = -2 * as.numeric(loglik), DfResid = df.residual(model_lmer),
      R2marginal <- r2(model_lmer, tolerance = 1e-1000)$R2_marginal,
      R2conditional <- r2(model_lmer, tolerance = 1e-1000)$R2_conditional,
      MSE = mse, RMSE = sqrt(mse)))
    colnames(table0) <- c(
      names(base_info), "UsedFormula", "AIC", "BIC", "logLik", "npar", "Deviance",
      "DfResid", "R2marginal", "R2conditional", "MSE", "RMSE")

    anova_table <- car::Anova(model_lmer, type = 3)
    anova_sig <- symnum(
      anova_table$`Pr(>Chisq)`, corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.01, 0.05, 1), symbols = c("**", "*", ""))
    all_terms <- rownames(anova_table)
    anova_table_F <- anova(model_lmer)
    anova_table_F_aligned <- anova_table_F[match(all_terms, rownames(anova_table_F)), , drop = FALSE]
    anova_F_significance <- symnum(
      anova_table_F_aligned$`Pr(>F)`, corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.01, 0.05, 1), symbols = c("**", "*", ""))
    table1 <- cbind(base_info, data.frame(
      Term = all_terms,
      DfChi = anova_table$Df, Chisq = anova_table$Chisq,
      PrChi = anova_table$`Pr(>Chisq)`, SigChi = as.character(anova_sig),
      SumSq = anova_table_F_aligned$`Sum Sq`, MeanSq = anova_table_F_aligned$`Mean Sq`,
      NumDF = anova_table_F_aligned$NumDF, DenDF = anova_table_F_aligned$DenDF,
      FValue = anova_table_F_aligned$`F value`, PrF = anova_table_F_aligned$`Pr(>F)`,
      SigF = anova_F_significance))

    desired_order <- gsub("1 \\| ", "", random_effects_terms)
    desired_order <- c("Residual", desired_order)
    random_effects <- VarCorr(model_lmer)
    random_effects_summary <- as.data.frame(random_effects)
    random_effects_summary <- random_effects_summary[match(desired_order, random_effects_summary$grp), ]
    colnames(random_effects_summary)[colnames(random_effects_summary) == "vcov"] <- "var"
    colnames(random_effects_summary)[colnames(random_effects_summary) == "sdcor"] <- "sd"
    random_effect_n <- ngrps(model_lmer)
    random_effect_n <- random_effect_n[desired_order[-1]]
    random_effects_comparison <- suppressMessages(ranova(model_lmer))
    table2 <- cbind(base_info, data.frame(
      Term = as.character(random_effects_summary$grp),
      SampleSize = c(NA, random_effect_n), Variance = random_effects_summary$var,
      StdDev = random_effects_summary$sd, npar = random_effects_comparison$npar,
      logLik = random_effects_comparison$logLik, AIC = random_effects_comparison$AIC,
      LRT = random_effects_comparison$LRT, Df = random_effects_comparison$Df,
      PValue = random_effects_comparison$`Pr(>Chisq)`,
      Sig = symnum(as.numeric(
        random_effects_comparison$`Pr(>Chisq)`), corr = FALSE, na = FALSE,
        cutpoints = c(0, 0.01, 0.05, 1), symbols = c("**", "*", ""))))

    overall_slope_emm <- emtrends(model_lmer, specs = ~ 1, var = fixed_effects_terms[2])
    overall_slope_df <- as.data.frame(summary(overall_slope_emm, infer = TRUE))
    overall_slope_df$Sig <- cut(overall_slope_df$p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("**", "*", ""), right = FALSE)
    overall_intercept_emm <- emmeans(model_lmer, specs = ~ 1, at = list(Period2 = 0))
    overall_intercept_df <- as.data.frame(summary(overall_intercept_emm, infer = TRUE))
    table3 <- data.frame(base_info, overall_slope_df[, 2:9], overall_intercept_df[, 2:8])
    colnames(table3) <- c(
      names(base_info),
      "Slope", "SloSE", "SloDf", "SloCLlower", "SloCLupper", "SloTratio", "SloP", "SloPSig",
      "Intercept", "IntSE", "IntDf", "IntCLlower", "IntCLupper", "IntTratio", "IntP")

    slopes_emm <- emtrends(model_lmer, specs = fixed_effects_terms[1], var = fixed_effects_terms[2])
    slope_groups <- multcomp::cld(slopes_emm, Letters = letters, reversed = TRUE)
    slope_groups$.group <- trimws(slope_groups$.group)
    slope_groups$Sig_noA <- if (all(slope_groups$.group == "a")) "" else slope_groups$.group
    all_levels_A <- data.frame(x = levels(data_clean[[fixed_effects_terms[1]]]))
    names(all_levels_A) <- fixed_effects_terms[1]
    slope_df <- all_levels_A %>% left_join(slope_groups, by = fixed_effects_terms[1]) %>% left_join(
      as.data.frame(summary(slopes_emm, infer = c(TRUE, TRUE))) %>%
        dplyr::select(all_of(fixed_effects_terms[1]), t.ratio, p.value), by = fixed_effects_terms[1])
    slope_df$Sig_stars <- cut(slope_df$p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("**", "*", ""), right = FALSE)
    slope_df_sorted <- slope_df[match(levels(data_clean[[fixed_effects_terms[1]]]), slope_df[[fixed_effects_terms[1]]]), ]
    intercepts_emm <- suppressMessages(emmeans(
      model_lmer, specs = fixed_effects_terms[1], at = setNames(list(0), fixed_effects_terms[2])))
    intercepts_df <- all_levels_A %>% left_join(
      as.data.frame(summary(intercepts_emm, infer = c(TRUE, TRUE))), by = fixed_effects_terms[1])
    intercepts_df_sorted <- intercepts_df[match(levels(data_clean[[fixed_effects_terms[1]]]), intercepts_df[[fixed_effects_terms[1]]]), ]
    table4 <- data.frame(base_info, slope_df_sorted[, 1:11], intercepts_df_sorted[, 2:8])
    colnames(table4) <- c(
      names(base_info), fixed_effects_terms[1],
      "Slope", "SloSE", "SloDf", "SloCLlower", "SloCLupper", "SloSig", "SloSignoA", "SloTratio", "SloP", "SloPSig",
      "Intercept", "IntSE", "IntDf", "IntCLlower", "IntCLupper", "IntTratio", "IntP")

    tables <- list(
      table0 = table0, table1 = table1, table2 = table2,
      table3 = table3, table4 = table4)
    tables <- lapply(tables, function(df) {
      rownames(df) <- NULL
      df[] <- lapply(df, function(x) {
        if (is.numeric(x)) {x[is.nan(x) | is.infinite(x)] <- NA_real_
        } else if (is.character(x) || is.factor(x)) {
          x <- as.character(x)
          x[is.na(x)] <- ""}
        x})
      df})
    return(tables)})

  result <- safe_analyze()
  if (!is.null(result$error)) {
    warning(sprintf(
      "Analysis failed for %s - Error: %s",
      response_var, result$error$message))
    return(NULL)}
  return(result$result)
}

analyze_variable_one <- function(
  data, response_var, formula_effects, random_factors) {
  safe_analyze <- safely(function() {
    data_clean <- data[!is.na(data[[response_var]]), ]
    valid_random_effects <- c()
    for (factor in random_factors) {
      levels <- length(unique(na.omit(data_clean[[factor]])))
      if (levels > 1) {valid_random_effects <- c(valid_random_effects, paste("(1 |", factor, ")"))}}
    if (length(valid_random_effects) == 0) return(NULL)
    random_effects_used <- paste(valid_random_effects, collapse = " + ")
    used_formula <- as.formula(paste(response_var, "~", formula_effects, "+", random_effects_used))
    formula_terms <- attr(terms(used_formula), "term.labels")
    fixed_effects_terms <- formula_terms[!grepl("\\|", formula_terms)]
    random_effects_terms <- formula_terms[grepl("\\|", formula_terms)]
    model_lmer <- suppressMessages(lmer(used_formula, data = data_clean, REML = TRUE, control = lmerControl(optimizer = "Nelder_Mead")))
    residuals_lmer <- residuals(model_lmer)
    base_info <- data.frame(Variable = response_var)

    loglik <- logLik(model_lmer)
    mse <- mean(residuals_lmer^2)
    table0 <- cbind(base_info, data.frame(
      UsedFormula = gsub("\\s+", " ", paste(deparse(used_formula), collapse = "")),
      AIC = AIC(model_lmer), BIC = BIC(model_lmer),
      logLik = as.numeric(loglik), npar = attr(loglik, "df"),
      Deviance = -2 * as.numeric(loglik), DfResid = df.residual(model_lmer),
      R2marginal <- r2(model_lmer, tolerance = 1e-1000)$R2_marginal,
      R2conditional <- r2(model_lmer, tolerance = 1e-1000)$R2_conditional,
      MSE = mse, RMSE = sqrt(mse)))
    colnames(table0) <- c(
      names(base_info), "UsedFormula", "AIC", "BIC", "logLik", "npar", "Deviance",
      "DfResid", "R2marginal", "R2conditional", "MSE", "RMSE")

    anova_table <- car::Anova(model_lmer, type = 3)
    anova_sig <- symnum(
      anova_table$`Pr(>Chisq)`, corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.01, 0.05, 1), symbols = c("**", "*", ""))
    all_terms <- rownames(anova_table)
    anova_table_F <- anova(model_lmer)
    anova_table_F_aligned <- anova_table_F[match(all_terms, rownames(anova_table_F)), , drop = FALSE]
    anova_F_significance <- symnum(
      anova_table_F_aligned$`Pr(>F)`, corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.01, 0.05, 1), symbols = c("**", "*", ""))
    table1 <- cbind(base_info, data.frame(
      Term = all_terms,
      DfChi = anova_table$Df, Chisq = anova_table$Chisq,
      PrChi = anova_table$`Pr(>Chisq)`, SigChi = as.character(anova_sig),
      SumSq = anova_table_F_aligned$`Sum Sq`, MeanSq = anova_table_F_aligned$`Mean Sq`,
      NumDF = anova_table_F_aligned$NumDF, DenDF = anova_table_F_aligned$DenDF,
      FValue = anova_table_F_aligned$`F value`, PrF = anova_table_F_aligned$`Pr(>F)`,
      SigF = anova_F_significance))

    desired_order <- gsub("1 \\| ", "", random_effects_terms)
    desired_order <- c("Residual", desired_order)
    random_effects <- VarCorr(model_lmer)
    random_effects_summary <- as.data.frame(random_effects)
    random_effects_summary <- random_effects_summary[match(desired_order, random_effects_summary$grp), ]
    colnames(random_effects_summary)[colnames(random_effects_summary) == "vcov"] <- "var"
    colnames(random_effects_summary)[colnames(random_effects_summary) == "sdcor"] <- "sd"
    random_effect_n <- ngrps(model_lmer)
    random_effect_n <- random_effect_n[desired_order[-1]]
    random_effects_comparison <- suppressMessages(ranova(model_lmer))
    table2 <- cbind(base_info, data.frame(
      Term = as.character(random_effects_summary$grp),
      SampleSize = c(NA, random_effect_n), Variance = random_effects_summary$var,
      StdDev = random_effects_summary$sd, npar = random_effects_comparison$npar,
      logLik = random_effects_comparison$logLik, AIC = random_effects_comparison$AIC,
      LRT = random_effects_comparison$LRT, Df = random_effects_comparison$Df,
      PValue = random_effects_comparison$`Pr(>Chisq)`,
      Sig = symnum(as.numeric(
        random_effects_comparison$`Pr(>Chisq)`), corr = FALSE, na = FALSE,
        cutpoints = c(0, 0.01, 0.05, 1), symbols = c("**", "*", ""))))

    summary_stats_A <- data_clean %>%
      group_by(!!sym(fixed_effects_terms[1]), .drop = FALSE) %>% summarise(
        Mean = mean(.data[[response_var]], na.rm = TRUE), SD = sd(.data[[response_var]], na.rm = TRUE),
        n = sum(!is.na(.data[[response_var]])), SE = ifelse(n > 0, SD / sqrt(n), NA_real_),
        lower = ifelse(n > 0, Mean - qt(0.975, n - 1) * SE, NA_real_),
        upper = ifelse(n > 0, Mean + qt(0.975, n - 1) * SE, NA_real_), .groups = "drop")
    emm_A <- suppressMessages(emmeans(model_lmer, reformulate(fixed_effects_terms[1])))
    cld_A_results <- suppressMessages(
      multcomp::cld(emm_A, adjust = "tukey", Letters = letters, sort = TRUE, reverse = TRUE))
    cld_A_results$.group <- gsub(" ", "", cld_A_results$.group)
    cld_A_results$Sig_noA <- if (all(cld_A_results$.group == "a")) "" else cld_A_results$.group
    all_levels_A <- data.frame(x = levels(data_clean[[fixed_effects_terms[1]]]))
    names(all_levels_A) <- fixed_effects_terms[1]
    cld_A_resultsNew <- merge(
      all_levels_A, cld_A_results, by = fixed_effects_terms[1], all.x = TRUE)
    cld_A_sorted <- cld_A_resultsNew[match(levels(data_clean[[fixed_effects_terms[1]]]), cld_A_resultsNew[[fixed_effects_terms[1]]]), ]
    table3 <- data.frame(
      base_info,
      summary_stats_A[[fixed_effects_terms[1]]],
      summary_stats_A[, c("Mean", "SD", "SE", "n", "lower", "upper")],
      Sig = cld_A_sorted$.group, Sig_noA = cld_A_sorted$Sig_noA,
      cld_A_sorted[, c("emmean", "SE", "df", "lower.CL", "upper.CL")])
    colnames(table3) <- c(
      names(base_info), fixed_effects_terms[1], "Mean", "SD", "SE", "n",
      "CLlower", "CLupper", "Sig", "SignoA",
      "emmean", "SEemm", "dfemm", "CLloweremm", "CLupperemm")

    tables <- list(table0 = table0, table1 = table1, table2 = table2, table3 = table3)
    tables <- lapply(tables, function(df) {
      rownames(df) <- NULL
      df[] <- lapply(df, function(x) {
        if (is.numeric(x)) {x[is.nan(x) | is.infinite(x)] <- NA_real_
        } else if (is.character(x) || is.factor(x)) {
          x <- as.character(x)
          x[is.na(x)] <- ""}
        x})
      df})
    return(tables)})

  result <- safe_analyze()
  if (!is.null(result$error)) {
    warning(sprintf(
      "Analysis failed for %s - Error: %s",
      response_var, result$error$message))
    return(NULL)}
  return(result$result)
}

create_three_line <- function(tukey_table, factor_var, factor_levels) {
  if (is.null(tukey_table)) return(NULL)
  three_line_table <- tukey_table %>% mutate(
    Mean = sprintf("%.2f", as.numeric(Mean)), SE = sprintf("%.2f", as.numeric(SE)),
    value = if_else(
      is.na(SignoA) | SignoA == "", paste0(Mean, " ± ", SE), paste0(Mean, " ± ", SE, " ", SignoA))) %>%
    dplyr::select(all_of(factor_var), Variable, value) %>%
    pivot_wider(names_from = Variable, values_from = value) %>%
    arrange(factor(!!sym(factor_var), levels = factor_levels))
  return(three_line_table)}

nice_breaks_2_3 <- function(x) {
  rng <- range(x, na.rm = TRUE, finite = TRUE)
  lo <- rng[1]
  hi <- rng[2]
  if (!is.finite(lo) || !is.finite(hi)) return(NULL)
  if (lo == hi) return(lo)
  span <- hi - lo
  nice_base <- c(1, 2, 2.5, 5, 10)
  expo <- seq(floor(log10(span)) - 2, ceiling(log10(span)) + 2)
  steps <- sort(unique(as.vector(outer(nice_base, 10^expo))))
  cand_breaks <- list()
  cand_info <- data.frame(
    id = integer(), n = integer(), coverage = double(), center_dev = double(), step = double())
  for (s in steps) {
    b0 <- ceiling(lo / s) * s
    b1 <- floor(hi / s) * s
    if (b0 > b1) next
    b <- seq(b0, b1, by = s)
    digits <- max(0, ceiling(-log10(s)) + 2)
    b <- unique(round(b, digits))
    b <- b[b >= lo & b <= hi]
    n_b <- length(b)
    if (n_b >= 2 && n_b <= 3) {
      cand_breaks[[length(cand_breaks) + 1]] <- b
      cand_info <- rbind(cand_info, data.frame(
        id = length(cand_breaks), n = n_b, coverage = (max(b) - min(b)) / span, center_dev = abs(mean(b) - mean(c(lo, hi))), step = s))}
  }
  if (nrow(cand_info) > 0) {
    cand_info <- cand_info[order(-cand_info$n, -cand_info$coverage, cand_info$center_dev, cand_info$step), ]
    return(cand_breaks[[cand_info$id[1]]])
  }
  b <- pretty(c(lo, hi), n = 5)
  b <- unique(b[b >= lo & b <= hi])
  if (length(b) >= 3) {
    idx <- unique(round(seq(1, length(b), length.out = 3)))
    return(b[idx])}
  if (length(b) == 2) return(b)
  return(c(lo, hi))
}

my_palette <- c("#D3A7AE", "#BCAAD4", "#7EBEB1", "#ADB974", "#D6AB85", "#F594D1", "#ACC7E2")

create_raincloud <- function(dataforplot, y_title = "Response", sigtext = NULL, point_size = 2) {
  if (!all(c("xTreat", "yResp") %in% names(dataforplot))) {stop("dataforplot lack xTreat yResp")}
  if (!is.factor(dataforplot$xTreat)) {dataforplot$xTreat <- as.factor(dataforplot$xTreat)}
  treat_levels <- levels(dataforplot$xTreat <- droplevels(as.factor(dataforplot$xTreat)))
  num_levels <- length(treat_levels)
  color_values <- if (num_levels == 2) {c('#D3A7AE', '#95B7DA')} else {my_palette[1:num_levels]}
  p <- ggplot(data = dataforplot, mapping = aes(x = xTreat, y = yResp, fill = xTreat))
  y_range <- ggplot_build(p + scale_y_continuous(expand = expansion(mult = 0.2, add = 0)))$layout$panel_params[[1]]$y.range
  p <- p + sm_raincloud(
    data = subset(dataforplot, xTreat == treat_levels[1]),
    which_side = 'left', position = position_nudge(x = -0.0), show.legend = FALSE,
    boxplot.params = list(outlier.shape = NA), violin.params = list(width = 0.5), point.params = list(
      show.legend = TRUE, alpha = 0.3, size = point_size, stroke = 0,
      position = sdamr::position_jitternudge(nudge.x = +0.1, seed = 10, jitter.width = 0.06))) +
    sm_raincloud(
      data = subset(dataforplot, xTreat != treat_levels[1]),
      which_side = 'right', position = position_nudge(x = +0.0), show.legend = FALSE,
      boxplot.params = list(outlier.shape = NA), violin.params = list(width = 0.5), point.params = list(
        show.legend = TRUE, alpha = 0.3, size = point_size, stroke = 0,
        position = sdamr::position_jitternudge(nudge.x = -0.1, seed = 10, jitter.width = 0.06))) +
    scale_fill_manual(values = color_values) + labs(y = y_title) + theme_classic() + theme(
      plot.background = element_blank(), panel.background = element_blank(),
      axis.text = element_text(size = 14), axis.text.y = element_text(size = 13),
      axis.title = element_text(size = 15), axis.title.x = element_blank(), legend.position = 'none',
      plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) +
    scale_y_continuous(expand = c(0, 0), breaks = nice_breaks_2_3, labels = scales::label_number(trim = TRUE, big.mark = "")) +
    coord_cartesian(ylim = y_range, clip = "off")
  if (!is.null(sigtext) && length(sigtext) > 0) {
    y_pos <- y_range[1]
    x_pos <- if (length(sigtext) == 1) 1.5 else seq_along(sigtext)
    p <- p + annotate("text", y = y_pos, x = x_pos, label = sigtext, size = 13/2.8346, vjust = -0.4)}
  return(p)
}

create_forest <- function(dataforplot, y_title = "Response", sigtext = NULL) {
  if (!all(c("xTreat", "yResp") %in% names(dataforplot))) {stop("dataforplot lack xTreat yResp")}
  if (!is.factor(dataforplot$xTreat)) {dataforplot$xTreat <- as.factor(dataforplot$xTreat)}
  treat_levels <- levels(dataforplot$xTreat <- droplevels(as.factor(dataforplot$xTreat)))
  num_levels <- length(treat_levels)
  color_values <- if (num_levels == 2) {c('#D3A7AE', '#95B7DA')} else {my_palette[1:num_levels]}
  set.seed(123)
  p <- ggplot(data = dataforplot, mapping = aes(x = xTreat, y = yResp, color = xTreat, group = xTreat)) +
    sm_forest(
      errorbar_type = "ci", points = FALSE, refLine = FALSE, avgPoint.params = list(shape = 16, size = 3),
      err.params = list(size = 1, color = NULL, linetype = "solid", show.legend = FALSE))
  set.seed(123)
  y_range <- ggplot_build(p + scale_y_continuous(expand = expansion(mult = 0.2, add = 0)))$layout$panel_params[[1]]$y.range
  set.seed(123)
  p <- p + scale_color_manual(values = color_values) + labs(y = y_title) + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(size = 14), axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 15), axis.title.x = element_blank(),
    legend.title = element_blank(), legend.text = element_text(size = 12),
    legend.position = c(1,1), legend.justification = c(1,1), legend.direction = "horizontal",
    legend.background = element_blank(), legend.key = element_blank(),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 0.5/1) +
    theme(legend.position = 'none') +
    scale_y_continuous(expand = c(0, 0), breaks = nice_breaks_2_3, labels = scales::label_number(trim = TRUE, big.mark = "")) +
    coord_cartesian(ylim = y_range, clip = "off")
  if (!is.null(sigtext) && length(sigtext) > 0) {
    y_pos <- y_range[1]
    x_pos <- if (length(sigtext) == 1) 1.5 else seq_along(sigtext)
    set.seed(123)
    p <- p + annotate("text", y = y_pos, x = x_pos, label = sigtext, size = 13/2.8346, vjust = -0.4)}
  set.seed(123)
  return(p)
}

create_raincloud_sep <- function(
    dataforplot, y_title = "Response", sigtext = NULL, point_size = 2, dodge_width = 0.9) {
  if (!all(c("xTreat", "yResp") %in% names(dataforplot))) {stop("dataforplot lack xTreat yResp")}
  if (!is.factor(dataforplot$xTreat)) {dataforplot$xTreat <- as.factor(dataforplot$xTreat)}
  if (!as.character(NameLevel) %in% names(dataforplot)) {stop(paste0("dataforplot lack ", NameLevel))}
  treat_levels <- levels(dataforplot$xTreat <- droplevels(as.factor(dataforplot$xTreat)))
  num_levels <- length(treat_levels)
  color_values <- if (num_levels == 2) {c('#D3A7AE', '#95B7DA')} else {my_palette[1:num_levels]}
  p <- ggplot(data = dataforplot, mapping = aes(x = !!sym(NameLevel), y = yResp, fill = xTreat)) +
    geom_stripped_cols(odd = 'transparent', even = 'grey98')
  y_range <- ggplot_build(p + scale_y_continuous(expand = expansion(mult = 0.2, add = 0)))$layout$panel_params[[1]]$y.range
  if (num_levels == 2) {
    p <- p + sm_raincloud(
      data = subset(dataforplot, xTreat == treat_levels[1]),
      which_side = 'left', position = position_nudge(x = -0.13), show.legend = FALSE,
      boxplot.params = list(outlier.shape = NA), violin.params = list(width = 0.5), point.params = list(
        show.legend = TRUE, alpha = 0.3, size = point_size, stroke = 0,
        position = sdamr::position_jitternudge(nudge.x = -0.06, seed = 10, jitter.width = 0.06))) +
      sm_raincloud(
        data = subset(dataforplot, xTreat == treat_levels[2]),
        which_side = 'right', position = position_nudge(x = +0.13), show.legend = FALSE,
        boxplot.params = list(outlier.shape = NA), violin.params = list(width = 0.5), point.params = list(
          show.legend = TRUE, alpha = 0.3, size = point_size, stroke = 0,
          position = sdamr::position_jitternudge(nudge.x = +0.06, seed = 10, jitter.width = 0.06)))
  } else {
    p <- p + sm_raincloud(
      which_side = 'right', position = position_dodge(width = dodge_width), show.legend = FALSE,
      boxplot.params = list(outlier.shape = NA, width = 0.4), violin.params = list(width = 0.8), point.params = list(
        show.legend = TRUE, alpha = 0.3, size = point_size, stroke = 0,
        position = position_jitterdodge(dodge.width = dodge_width, seed = 10, jitter.width = 0.03)))}
  p <- p +
    scale_fill_manual(values = color_values) + labs(y = y_title) + theme_classic() + theme(
      plot.background = element_blank(), panel.background = element_blank(),
      axis.text = element_text(size = 14), axis.text.y = element_text(size = 13),
      axis.title = element_text(size = 15), axis.title.x = element_blank(),
      legend.title = element_blank(), legend.text = element_text(size = 12),
      legend.position = c(1,1), legend.justification = c(1,1), legend.direction = "horizontal",
      legend.background = element_blank(), legend.key = element_blank(),
      plot.margin = margin(0, 0, 0, 0), aspect.ratio = 0.8/length(unique(dataforplot[[NameLevel]]))) +
    theme(legend.position = 'none') +
    scale_y_continuous(expand = c(0, 0), breaks = nice_breaks_2_3, labels = scales::label_number(trim = TRUE, big.mark = "")) +
    scale_x_discrete(expand = expansion(0)) +
    coord_cartesian(ylim = y_range, xlim = c(0.5, length(unique(dataforplot[[NameLevel]])) + 0.5), clip = "off")
  if (!is.null(sigtext) && length(sigtext) > 0) {
    y_pos <- y_range[1]
    if (num_levels == 2) {x_pos <- seq_along(sigtext)} else {
      group_centers <- seq_along(unique(dataforplot[[NameLevel]]))
      offsets <- seq((-dodge_width/2 + dodge_width/(2*num_levels)), (dodge_width/2 - dodge_width/(2*num_levels)), length.out = num_levels)
      pos_grid <- expand.grid(offset = offsets, center = group_centers)
      x_pos <- pos_grid$center + pos_grid$offset}
    p <- p + annotate("text", y = y_pos, x = x_pos, label = sigtext, size = 13/2.8346, vjust = -0.4)}
  return(p)
}

create_forest_sep <- function(dataforplot, y_title = "Response", sigtext = NULL, dodge_width = 0.9) {
  if (!all(c("xTreat", "yResp") %in% names(dataforplot))) {stop("dataforplot lack xTreat yResp")}
  if (!is.factor(dataforplot$xTreat)) {dataforplot$xTreat <- as.factor(dataforplot$xTreat)}
  if (!as.character(NameLevel) %in% names(dataforplot)) {stop(paste0("dataforplot lack ", NameLevel))}
  treat_levels <- levels(dataforplot$xTreat <- droplevels(as.factor(dataforplot$xTreat)))
  num_levels <- length(treat_levels)
  color_values <- if (num_levels == 2) {c('#D3A7AE', '#95B7DA')} else {my_palette[1:num_levels]}
  set.seed(123)
  p0 <- ggplot(data = dataforplot, mapping = aes(x = !!sym(NameLevel), y = yResp, color = xTreat, group = xTreat)) +
    sm_forest(
      position = position_dodge(width = dodge_width), errorbar_type = "ci",
      avgPoint.params = list(shape = 16, size = 3), points = FALSE, refLine = FALSE,
      err.params = list(size = 1, color = NULL, linetype = "solid", show.legend = FALSE))
  set.seed(123)
  y_range <- ggplot_build(p0 + scale_y_continuous(expand = expansion(mult = 0.2, add = 0)))$layout$panel_params[[1]]$y.range
  set.seed(123)
  p <- ggplot(data = dataforplot, mapping = aes(x = !!sym(NameLevel), y = yResp, color = xTreat, group = xTreat)) +
    geom_stripped_cols(odd = 'transparent', even = 'grey98', color = NA) +
    sm_forest(
      errorbar_type = "ci", points = FALSE, refLine = FALSE,
      position = position_dodge(width = dodge_width), avgPoint.params = list(shape = 16, size = 3),
      err.params = list(size = 1, color = NULL, linetype = "solid", show.legend = FALSE)) +
    scale_color_manual(values = color_values) + labs(y = y_title) + theme_classic() + theme(
      plot.background = element_blank(), panel.background = element_blank(),
      axis.text = element_text(size = 14), axis.text.y = element_text(size = 13),
      axis.title = element_text(size = 15), axis.title.x = element_blank(),
      legend.title = element_blank(), legend.text = element_text(size = 12),
      legend.position = c(1,1), legend.justification = c(1,1), legend.direction = "horizontal",
      legend.background = element_blank(), legend.key = element_blank(),
      plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1.2/length(unique(dataforplot[[NameLevel]]))) +
    theme(legend.position = 'none') +
    scale_y_continuous(expand = c(0, 0), breaks = nice_breaks_2_3, labels = scales::label_number(trim = TRUE, big.mark = "")) +
    scale_x_discrete(expand = expansion(0)) + coord_cartesian(
      ylim = y_range, xlim = c(0.5, length(unique(dataforplot[[NameLevel]])) + 0.5), clip = "off")
  if (!is.null(sigtext) && length(sigtext) > 0) {
    y_pos <- y_range[1]
    if (num_levels == 2) {x_pos <- seq_along(sigtext)} else {
      group_centers <- seq_along(unique(dataforplot[[NameLevel]]))
      offsets <- seq((-dodge_width/2 + dodge_width/(2*num_levels)), (dodge_width/2 - dodge_width/(2*num_levels)), length.out = num_levels)
      pos_grid <- expand.grid(offset = offsets, center = group_centers)
      x_pos <- pos_grid$center + pos_grid$offset}
    set.seed(123)
    p <- p + annotate("text", y = y_pos, x = x_pos, label = sigtext, size = 13/2.8346, vjust = -0.4)}
  set.seed(123)
  return(p)
}

create_regression <- function(dataforplot, stat_data, y_title = "Response", x_title = "Variable") {
  if (!all(c("xTreat", "yResp") %in% names(dataforplot))) {stop("dataforplot lack xTreat yResp")}
  if (!is.factor(dataforplot$xTreat)) {dataforplot$xTreat <- as.factor(dataforplot$xTreat)}
  if (!as.character(NameLevel) %in% names(dataforplot)) {stop(paste0("dataforplot lack ", NameLevel))}
  required_stat_cols <- c("Slope", "Intercept", "xTreat", "Sigused1", "Sigused2")
  if (!all(required_stat_cols %in% names(stat_data))) {stop(paste0("stat_data lack ", paste(required_stat_cols, collapse = ", ")))}
  treat_levels <- levels(dataforplot$xTreat <- droplevels(as.factor(dataforplot$xTreat)))
  num_levels <- length(treat_levels)
  color_values <- if (num_levels == 2) {c('#D3A7AE', '#95B7DA')} else {my_palette[1:num_levels]}
  p <- ggplot(data = dataforplot, aes(x = !!sym(NameLevel), y = yResp, color = xTreat))
  x_range <- ggplot_build(p)$layout$panel_params[[1]]$x.range
  x_anno <- (mean(x_range) + x_range[1])/2
  y_range <- ggplot_build(p + scale_y_continuous(expand = expansion(mult = c(0.2, 0.5), add = 0)))$layout$panel_params[[1]]$y.range
  p <- p + geom_point(size = 2, alpha = 0.2) +
    geom_abline(
      data = stat_data, size = 1, show.legend = FALSE,
      aes(slope = Slope, intercept = Intercept, color = xTreat, linetype = ifelse(is.na(Sigused1) | Sigused1 == "", "twodash", "solid"))) +
    scale_linetype_identity() + scale_color_manual(values = color_values) +
    labs(x = x_title, y = y_title) + theme_classic() + theme(
      plot.background = element_blank(), panel.background = element_blank(),
      axis.text = element_text(size = 14), axis.text.y = element_text(size = 13),
      axis.title = element_text(size = 15),
      legend.title = element_blank(), legend.text = element_text(size = 12),
      legend.position = c(1,1), legend.justification = c(1,1), legend.direction = "horizontal",
      legend.background = element_blank(), legend.key = element_blank(),
      plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1.7/3) +
    theme(legend.position = 'none') +
    scale_y_continuous(expand = c(0, 0), breaks = nice_breaks_2_3, labels = scales::label_number(trim = TRUE, big.mark = "")) +
    coord_cartesian(ylim = y_range, clip = "off") +
    annotate(
      "text", x = x_anno, hjust = 0, y = Inf, color = color_values, size = 12/2.8346, vjust = seq(1.1, by = 1.1, length.out = nrow(stat_data)),
      label = paste0(
        "Slope = ", paste0(format(stat_data$Slope, scientific = TRUE, digits = 2)),
        ifelse(is.na(stat_data$Sigused1), "", paste0(stat_data$Sigused1)), " ",
        ifelse(is.na(stat_data$Sigused2), "", paste0(stat_data$Sigused2))))
  return(p)
}


# Data analysis -----------------------------------------------------------
setwd('D:/R') # set current working directory
rm(list = ls()) # clear the environment
cat("\014") # clear the console

pacman::p_load(openxlsx, tidyverse)
dataOrigin <- read.xlsx("Input/FieldData.xlsx", "dataall(plain)") %>% mutate(
  Site = factor(Site, levels = c("Zhangbei", "Ulanqab", "Youyu", "Chifeng")),
  Year = as.numeric(Year),
  Crop = factor(Crop, levels = c("Oat", "Soybean")),
  Treat_All = factor(Treat_All, levels = c("CK", "PP-P", "PP-F", "PLA-P", 'PLA-F')),
  Rep = factor(Rep, levels = c('1', '2', '3', '4')),
  Treat_Is = factor(if_else(
    Treat_All == "CK", "Without", "With"), levels = c("Without", "With"))
) %>% rename_with(~ str_replace(., "^(\\d)", "X\\1"))
dataOrigin_202125MonoMaturity <- dataOrigin %>% filter(
  Year %in% 2021:2025, Period == "Maturity", Pattern == "monoculture")
dataOrigin_202125MonoMaturity_clean <- dataOrigin_202125MonoMaturity %>% dplyr::select(where(~!all(is.na(.))))

## SQI and EMF calculation ---------------------------------------------------------------
# Data normalization
dataOrigin_202125MonoMaturity_norm <- dataOrigin_202125MonoMaturity_clean %>%
  mutate(across(Height:PER, ~ {
    if(all(is.na(.))) return(.)
    reverse_cols <- # indicators smaller values are better
      c("BD", "pH", "POX", "PER")
    col_name <- cur_column()
    if(col_name %in% reverse_cols) {(max(., na.rm = TRUE) - .) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))
    } else {(. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))}}, .names = "{.col}"))

# Calculate SQI: normalize indicators, assign weights using PCA, and perform weighted summation
data_for_pca_202125MonoMaturity <- dataOrigin_202125MonoMaturity_norm %>%
  dplyr::select(SWC:PER) %>%
  filter(if_all(everything(), ~ !is.na(.) & !is.infinite(.)))
# weights are based on the norm values from PCA (comprehensive loadings, magnitude (length) of the vector representing the variable in the multi-dimensional space）
weights_202125MonoMaturity <- prcomp(data_for_pca_202125MonoMaturity, center = TRUE, scale. = TRUE) |>
  (\(pca) {kept <- pca$sdev^2 > 1 # Retain PCs with eigenvalue > 1 (Kaiser criterion)
  contributions <- apply(pca$rotation[, kept], 1, \(u) sqrt(sum((u^2) * pca$sdev[kept]^2)))
  contributions / sum(contributions)})()

# Calculate EMF: normalize indicators using the group mean method
dataOrigin_202125MonoMaturity_norm <- dataOrigin_202125MonoMaturity_norm %>% mutate(
  PlantProd = rowMeans(dplyr::select(., Yield, Biomass), na.rm = TRUE),
  PlantGrow = rowMeans(dplyr::select(
    ., Height, SpikeLength, SpikeNumber, KernelsPerSpike, Weight1000Grain,
    PodNumberPerPlant, SeedNumberPerPod, Weight100Seed), na.rm = TRUE),
  PlantForBelow = rowMeans(dplyr::select(
    ., BiomassBelowground, RootLength, RootSurfaceArea, RootVolume), na.rm = TRUE),
  SoilNutriTurnover = rowMeans(dplyr::select(., NH4, NO3, AvailP), na.rm = TRUE),
  SoilActi = rowMeans(dplyr::select(., BG, BX, CBH, LAP, NAG, ALP), na.rm = TRUE),
  SoilPhy = rowMeans(dplyr::select(., BD, SWC), na.rm = TRUE),
  SoilChem = rowMeans(dplyr::select(., pH, EC), na.rm = TRUE),
  SoilCSta = rowMeans(dplyr::select(., POX, PER), na.rm = TRUE)) %>% mutate(
    Support = rowMeans(dplyr::select(
      ., PlantGrow, PlantForBelow, SoilNutriTurnover, SoilActi), na.rm = TRUE),
    Regulation = rowMeans(dplyr::select(., SoilPhy, SoilChem, SoilCSta), na.rm = TRUE)) %>%
  mutate(EMF = rowMeans(dplyr::select(., PlantProd, Support, Regulation), na.rm = TRUE)) %>%
  mutate(across(c(PlantProd:EMF), ~ ifelse(is.nan(.), NA, .)))

dataTotal_202125MonoMaturity <- dataOrigin_202125MonoMaturity_clean %>% mutate(SQI = {
  norm_data <- dataOrigin_202125MonoMaturity_norm %>%
    dplyr::select(SWC:PER) %>% as.matrix()
  sapply(1:nrow(norm_data), function(i) {
    row_data <- norm_data[i, ]
    valid <- !is.na(row_data)
    if (sum(valid) > 0) { # calculate row by row using only non-missing values
      available_weights <- weights_202125MonoMaturity[valid]
      adjusted_weights <- available_weights / sum(available_weights) # rescale the weights so that they sum to 1 (preserving the weight proportions)
      sum(row_data[valid] * adjusted_weights)} else NA})
}) %>% bind_cols(
  dataOrigin_202125MonoMaturity_norm %>% dplyr::select(
    PlantProd, PlantGrow, PlantForBelow,
    SoilNutriTurnover, SoilActi, SoilPhy, SoilChem, SoilCSta,
    Support, Regulation, EMF))

dataTotal_202125MonoMaturity <- dataTotal_202125MonoMaturity %>%
  mutate(Year2 = case_when( # Create a new column for subsequent regression analysis
    Year == "2021" ~ 1, Year == "2022" ~ 2,
    Year == "2023" ~ 3, Year == "2024" ~ 4, Year == "2025" ~ 5))
# Export the results
# if (!dir.exists("output")) {dir.create("output")}
# write.csv(dataTotal_202125MonoMaturity, "output/dataTotal202125MonoMaturity.csv", row.names = FALSE)

## Overall / Crop- / Site- / Time-specific statistical analysis -------------------------------------------------------------------
pacman::p_load(openxlsx, tidyverse, lmerTest, nortest, car, emmeans, multcomp,
               multcompView, DHARMa, performance, furrr, smplot2, sdamr, GGally)
dataTotal_202125MonoMaturity <- read.csv(
  "output/dataTotal202125MonoMaturity.csv") %>% mutate(
    Site = factor(Site, levels = c("Zhangbei", "Ulanqab", "Youyu", "Chifeng")),
    Year = as.numeric(Year),
    Year2 = as.numeric(Year2),
    Crop = factor(Crop, levels = c("Oat", "Soybean")),
    Treat_All = factor(Treat_All, levels = c("CK", "PP-P", "PP-F", "PLA-P", 'PLA-F')),
    Rep = factor(Rep, levels = c('1', '2', '3', '4')),
    PlotID = interaction(Site, Crop, Treat_All, Rep, drop = TRUE, lex.order = TRUE),
    Treat_Is = factor(Treat_Is, levels = c("Without", "With"))
  ) %>% rename_with(~ str_replace(., "^(\\d)", "X\\1"))
variables_original <-
  c("Yield", "SQI", "EMF")
variables <- c(variables_original)

### Overall statistical analysis -------------------------------------------------------------------
FactorAs <- c("Treat_All", "Treat_Is")
formula_Randoms <- c("Site", "Year", "Crop", "PlotID")
for (factor in FactorAs) {
  DataForAnal <- dataTotal_202125MonoMaturity
  FactorA <- factor
  ideal_A_order <- switch(
    FactorA,
    "Treat_All" = c("CK", "PP-P", "PP-F", "PLA-P", "PLA-F"),
    "Treat_Is"  = c("Without", "With"))
  existing_A <- ideal_A_order[ideal_A_order %in% unique(DataForAnal[[FactorA]])]
  DataForAnal[[FactorA]] <- factor(DataForAnal[[FactorA]], levels = existing_A)
  factorA_levels <- existing_A
  plan(multisession, workers = availableCores() - 1)
  results_list <- future_map(variables, function(var) {
    if (all(is.na(DataForAnal[[var]]))) return(NULL)
    analyze_variable_one(
      data = DataForAnal, response_var = var,
      formula_effects = FactorA, random_factors = formula_Randoms)},
    .progress = TRUE, .options = furrr_options(seed = TRUE)) %>% set_names(variables)
  results_list <- results_list %>% compact()
  plan(sequential)
  all_tables <- map(0:3, function(i) {
    table_key <- paste0("table", i)
    map_dfr(results_list, ~ .x[[table_key]])
  }) %>% set_names(paste0("table", 0:3))
  three_line_A_original <- create_three_line(all_tables$table3 %>% filter(
    Variable %in% variables_original), FactorA, factorA_levels)
  wb <- createWorkbook()
  sheet_names_original <-
    c("ModelEval", "FixedANOVA", "RandomTest", "Tukey")
  walk2(seq_along(sheet_names_original), sheet_names_original, ~ {
    table_data <- all_tables[[.x]] %>%
      filter(Variable %in% variables_original) %>%
      arrange(match(Variable, variables_original))
    if (nrow(table_data) > 0) {
      addWorksheet(wb, .y)
      writeData(wb, .y, table_data)}})
  addWorksheet(wb, "Line3")
  if (!is.null(three_line_A_original)) {writeData(wb, "Line3", three_line_A_original)}
  output_filename <- paste0("Output/Sta_202125_Overview_", FactorA, ".xlsx")
  saveWorkbook(wb, output_filename, overwrite = TRUE)
}

#### Overall statistical visualization ---------------------------------------------------------------------
dataSta99Overview99Treat_Is <- read.xlsx(
  "output/Sta_202125_Overview_Treat_Is.xlsx",
  "FixedANOVA", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  filter(Term != "(Intercept)") %>% mutate(Sigused = SigF)
dataSta99Overview99Treat_All <- read.xlsx(
  "output/Sta_202125_Overview_Treat_All.xlsx",
  "Tukey", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  mutate(Sigused = SignoA)
variables_original <-
  c("Yield", "SQI", "EMF")
resp_df <- data.frame(
  NameyResp = variables_original,
  y_title_text =
    c("Yield ~ (t ~ ha^-1)", "SQI", "EMF"), stringsAsFactors = FALSE)

NameLevel <- "Overview"
NamexTreat <- "Treat_Is"
# NamexTreat <- "Treat_All"
# if (!dir.exists("Output")) dir.create("Output", recursive = TRUE)
for (i in 1:nrow(resp_df)) {
  NameyResp <- resp_df$NameyResp[i]
  y_title <- parse(text = resp_df$y_title_text[i])
  stat_data <- get(paste0("dataSta99", NameLevel, "99", NamexTreat)) %>% filter(Variable == NameyResp)
  if (nrow(stat_data) == 0) {message("skip: ", NameyResp, " - stat_data is empty"); next}
  sigtext <- stat_data %>% pull(Sigused) %>% ifelse(is.na(.), "", .)
  dataforplot <- dataTotal_202125MonoMaturity %>% mutate(
    xTreat = !!sym(NamexTreat), yResp = !!sym(NameyResp)) %>% filter(!is.na(yResp))
  OutputName1 <- paste0("Output/Raincloud-", NameLevel, "-", NameyResp, "-", NamexTreat, ".pdf")
  p1 <- create_raincloud(
    dataforplot = dataforplot, y_title = y_title, sigtext = sigtext,
    point_size = 2.5)
  # point_size = 2)
  cairo_pdf(OutputName1, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54)
  print(p1); dev.off()
  OutputName2 <- paste0("Output/Forest-", NameLevel, "-", NameyResp, "-", NamexTreat, ".pdf")
  p2 <- create_forest(dataforplot = dataforplot, y_title = y_title, sigtext = sigtext)
  cairo_pdf(OutputName2, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54)
  print(p2); dev.off()
}

### Crop-specific statistical analysis -------------------------------------------------------------------
FactorAs <- c("Treat_All", "Treat_Is")
FactorB <- "Crop"
formula_Randoms <- c("Site", "Year", "PlotID")

for (factor in FactorAs) {
  DataForAnal <- dataTotal_202125MonoMaturity
  FactorA <- factor
  ideal_A_order <- switch(
    FactorA,
    "Treat_All" = c("CK", "PP-P", "PP-F", "PLA-P", "PLA-F"),
    "Treat_Is"  = c("Without", "With"))
  existing_A <- ideal_A_order[ideal_A_order %in% unique(DataForAnal[[FactorA]])]
  DataForAnal[[FactorA]] <- factor(DataForAnal[[FactorA]], levels = existing_A)
  factorA_levels <- existing_A
  ideal_B_order <- levels(DataForAnal[[FactorB]])
  existing_B <- ideal_B_order[ideal_B_order %in% unique(DataForAnal[[FactorB]])]
  DataForAnal[[FactorB]] <- factor(DataForAnal[[FactorB]], levels = existing_B)
  factorB_levels <- existing_B
  workbook_list <- list()
  for (factorb in factorB_levels) {
    DataForAnal2 <- DataForAnal %>% filter(get(FactorB) == factorb)
    plan(multisession, workers = availableCores() - 1)
    results_list <- future_map(variables, function(var) {
      if (all(is.na(DataForAnal2[[var]]))) return(NULL)
      analyze_variable_one(
        data = DataForAnal2, response_var = var,
        formula_effects = FactorA, random_factors = formula_Randoms)},
      .progress = TRUE, .options = furrr_options(seed = TRUE)) %>% set_names(variables)
    results_list <- results_list %>% compact()
    plan(sequential)
    all_tables <- map(0:3, function(i) {
      table_key <- paste0("table", i)
      map_dfr(results_list, ~ .x[[table_key]])
    }) %>% set_names(paste0("table", 0:3))
    three_line_A_original <- create_three_line(all_tables$table3 %>% filter(
      Variable %in% variables_original), FactorA, factorA_levels)
    wb <- createWorkbook()
    sheet_names_original <-
      c("ModelEval", "FixedANOVA", "RandomTest", "Tukey")
    walk2(seq_along(sheet_names_original), sheet_names_original, ~ {
      table_data <- all_tables[[.x]] %>%
        filter(Variable %in% variables_original) %>%
        arrange(match(Variable, variables_original))
      if (nrow(table_data) > 0) {
        addWorksheet(wb, .y)
        writeData(wb, .y, table_data)}})
    addWorksheet(wb, "Line3")
    if (!is.null(three_line_A_original)) {writeData(wb, "Line3", three_line_A_original)}
    workbook_list[[factorb]] <- wb
  }
  wb_combined <- createWorkbook()
  all_sheet_names <- unique(unlist(lapply(workbook_list, names)))
  all_sheet_names <- all_sheet_names[!sapply(all_sheet_names, is.null)]
  for (sheet_name in all_sheet_names) {
    all_sheets <- list()
    factorb_names <- c()
    all_columns <- c()
    for (factorb in existing_B) {
      if (factorb %in% names(workbook_list)) {
        wb <- workbook_list[[factorb]]
        if (sheet_name %in% names(wb)) {
          sheet_data <- readWorkbook(wb, sheet = sheet_name, colNames = TRUE)
          if (nrow(sheet_data) > 0) {
            sheet_data_with_factor <- data.frame(
              temp_col = rep(factorb, nrow(sheet_data)),
              sheet_data, stringsAsFactors = FALSE, check.names = FALSE)
            colnames(sheet_data_with_factor)[1] <- FactorB
            all_sheets[[factorb]] <- sheet_data_with_factor
            factorb_names <- c(factorb_names, factorb)
            all_columns <- union(all_columns, colnames(sheet_data_with_factor))
          }}}}
    if (length(all_sheets) == 0) next
    combined_data_list <- list()
    for (i in seq_along(all_sheets)) {
      factorb <- factorb_names[i]
      sheet_data <- all_sheets[[factorb]]
      missing_cols <- setdiff(all_columns, colnames(sheet_data))
      if (length(missing_cols) > 0) {
        for (col in missing_cols) {sheet_data[[col]] <- NA}}
      sheet_data <- sheet_data[, all_columns, drop = FALSE]
      combined_data_list[[i]] <- sheet_data}
    combined_data <- do.call(rbind, combined_data_list)
    addWorksheet(wb_combined, sheet_name)
    writeData(wb_combined, sheet_name, combined_data)
  }
  output_filename <- paste0("Output/Sta_202125_Sep_", FactorB, "_", FactorA, ".xlsx")
  saveWorkbook(wb_combined, output_filename, overwrite = TRUE)
}

#### Crop-specific statistical visualization ---------------------------------------------------------------------
dataSta99Crop99Treat_Is <- read.xlsx(
  "output/Sta_202125_Sep_Crop_Treat_Is.xlsx",
  "FixedANOVA", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  filter(Term != "(Intercept)") %>% mutate(Sigused = SigF)
dataSta99Crop99Treat_All <- read.xlsx(
  "output/Sta_202125_Sep_Crop_Treat_All.xlsx",
  "Tukey", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  mutate(Sigused = SignoA)
variables_original <-
  c("Yield", "SQI", "EMF")
resp_df <- data.frame(
  NameyResp = variables_original,
  y_title_text =
    c("Yield ~ (t ~ ha^-1)", "SQI", "EMF"), stringsAsFactors = FALSE)

NameLevel <- "Crop"
NamexTreat <- "Treat_Is"
# NamexTreat <- "Treat_All"
# if (!dir.exists("Output")) dir.create("Output", recursive = TRUE)
for (i in 1:nrow(resp_df)) {
  NameyResp <- resp_df$NameyResp[i]
  y_title <- parse(text = resp_df$y_title_text[i])
  stat_data <- get(paste0("dataSta99", NameLevel, "99", NamexTreat)) %>% filter(Variable == NameyResp)
  if (nrow(stat_data) == 0) {message("skip: ", NameyResp, " - stat_data is empty"); next}
  sigtext <- stat_data %>% pull(Sigused) %>% ifelse(is.na(.), "", .)
  dataforplot <- dataTotal_202125MonoMaturity %>% mutate(
    xTreat = !!sym(NamexTreat), yResp = !!sym(NameyResp)) %>% filter(!is.na(yResp)) %>%
    filter(!!sym(NameLevel) %in% (stat_data[[NameLevel]] %>% unique())) %>%
    mutate(!!sym(NameLevel) := factor(!!sym(NameLevel)))
  OutputName1 <- paste0("Output/Raincloud-", NameLevel, "-", NameyResp, "-", NamexTreat, ".pdf")
  p1 <- create_raincloud_sep(
    dataforplot = dataforplot, y_title = y_title, sigtext = sigtext, dodge_width = 0.9,
    point_size = 2.5)
  # point_size = 2)
  cairo_pdf(OutputName1, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54)
  print(p1); dev.off()
  OutputName2 <- paste0("Output/Forest-", NameLevel, "-", NameyResp, "-", NamexTreat, ".pdf")
  p2 <- create_forest_sep(
    dataforplot = dataforplot, y_title = y_title, sigtext = sigtext,
    dodge_width = 0.7)
  # dodge_width = 0.9)
  cairo_pdf(OutputName2, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54)
  print(p2); dev.off()
}

### Site-specific statistical analysis -------------------------------------------------------------------
FactorAs <- c("Treat_All", "Treat_Is")
FactorB <- "Site"
formula_Randoms <- c("Year", "Crop", "PlotID")

for (factor in FactorAs) {
  DataForAnal <- dataTotal_202125MonoMaturity
  FactorA <- factor
  ideal_A_order <- switch(
    FactorA,
    "Treat_All" = c("CK", "PP-P", "PP-F", "PLA-P", "PLA-F"),
    "Treat_Is"  = c("Without", "With"))
  existing_A <- ideal_A_order[ideal_A_order %in% unique(DataForAnal[[FactorA]])]
  DataForAnal[[FactorA]] <- factor(DataForAnal[[FactorA]], levels = existing_A)
  factorA_levels <- existing_A
  ideal_B_order <- levels(DataForAnal[[FactorB]])
  existing_B <- ideal_B_order[ideal_B_order %in% unique(DataForAnal[[FactorB]])]
  DataForAnal[[FactorB]] <- factor(DataForAnal[[FactorB]], levels = existing_B)
  factorB_levels <- existing_B
  workbook_list <- list()
  for (factorb in factorB_levels) {
    DataForAnal2 <- DataForAnal %>% filter(get(FactorB) == factorb)
    plan(multisession, workers = availableCores() - 1)
    results_list <- future_map(variables, function(var) {
      if (all(is.na(DataForAnal2[[var]]))) return(NULL)
      analyze_variable_one(
        data = DataForAnal2, response_var = var,
        formula_effects = FactorA, random_factors = formula_Randoms)},
      .progress = TRUE, .options = furrr_options(seed = TRUE)) %>% set_names(variables)
    results_list <- results_list %>% compact()
    plan(sequential)
    all_tables <- map(0:3, function(i) {
      table_key <- paste0("table", i)
      map_dfr(results_list, ~ .x[[table_key]])
    }) %>% set_names(paste0("table", 0:3))
    three_line_A_original <- create_three_line(all_tables$table3 %>% filter(
      Variable %in% variables_original), FactorA, factorA_levels)
    wb <- createWorkbook()
    sheet_names_original <- c("ModelEval", "FixedANOVA", "RandomTest", "Tukey")
    walk2(seq_along(sheet_names_original), sheet_names_original, ~ {
      table_data <- all_tables[[.x]] %>%
        filter(Variable %in% variables_original) %>%
        arrange(match(Variable, variables_original))
      if (nrow(table_data) > 0) {
        addWorksheet(wb, .y)
        writeData(wb, .y, table_data)}})
    addWorksheet(wb, "Line3")
    if (!is.null(three_line_A_original)) {writeData(wb, "Line3", three_line_A_original)}
    workbook_list[[factorb]] <- wb
  }
  wb_combined <- createWorkbook()
  all_sheet_names <- unique(unlist(lapply(workbook_list, names)))
  all_sheet_names <- all_sheet_names[!sapply(all_sheet_names, is.null)]
  for (sheet_name in all_sheet_names) {
    all_sheets <- list()
    factorb_names <- c()
    all_columns <- c()
    for (factorb in existing_B) {
      if (factorb %in% names(workbook_list)) {
        wb <- workbook_list[[factorb]]
        if (sheet_name %in% names(wb)) {
          sheet_data <- readWorkbook(wb, sheet = sheet_name, colNames = TRUE)
          if (nrow(sheet_data) > 0) {
            sheet_data_with_factor <- data.frame(
              temp_col = rep(factorb, nrow(sheet_data)),
              sheet_data, stringsAsFactors = FALSE, check.names = FALSE)
            colnames(sheet_data_with_factor)[1] <- FactorB
            all_sheets[[factorb]] <- sheet_data_with_factor
            factorb_names <- c(factorb_names, factorb)
            all_columns <- union(all_columns, colnames(sheet_data_with_factor))
          }}}}
    if (length(all_sheets) == 0) next
    combined_data_list <- list()
    for (i in seq_along(all_sheets)) {
      factorb <- factorb_names[i]
      sheet_data <- all_sheets[[factorb]]
      missing_cols <- setdiff(all_columns, colnames(sheet_data))
      if (length(missing_cols) > 0) {
        for (col in missing_cols) {sheet_data[[col]] <- NA}}
      sheet_data <- sheet_data[, all_columns, drop = FALSE]
      combined_data_list[[i]] <- sheet_data}
    combined_data <- do.call(rbind, combined_data_list)
    addWorksheet(wb_combined, sheet_name)
    writeData(wb_combined, sheet_name, combined_data)
  }
  output_filename <- paste0("Output/Sta_202125_Sep_", FactorB, "_", FactorA, ".xlsx")
  saveWorkbook(wb_combined, output_filename, overwrite = TRUE)
}

#### Site-specific statistical visualization ---------------------------------------------------------------------
dataSta99Site99Treat_Is <- read.xlsx(
  "output/Sta_202125_Sep_Site_Treat_Is.xlsx",
  "FixedANOVA", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  filter(Term != "(Intercept)") %>% mutate(Sigused = SigF)
dataSta99Site99Treat_All <- read.xlsx(
  "output/Sta_202125_Sep_Site_Treat_All.xlsx",
  "Tukey", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  mutate(Sigused = SignoA)
variables_original <-
  c("Yield", "SQI", "EMF")
resp_df <- data.frame(
  NameyResp = variables_original,
  y_title_text =
    c("Yield ~ (t ~ ha^-1)", "SQI", "EMF"), stringsAsFactors = FALSE)

NameLevel <- "Site"
NamexTreat <- "Treat_Is"
# NamexTreat <- "Treat_All"
for (i in 1:nrow(resp_df)) {
  NameyResp <- resp_df$NameyResp[i]
  y_title <- parse(text = resp_df$y_title_text[i])
  stat_data <- get(paste0("dataSta99", NameLevel, "99", NamexTreat)) %>% filter(Variable == NameyResp)
  if (nrow(stat_data) == 0) {message("skip: ", NameyResp, " - stat_data is empty"); next}
  sigtext <- stat_data %>% pull(Sigused) %>% ifelse(is.na(.), "", .)
  dataforplot <- dataTotal_202125MonoMaturity %>% mutate(
    xTreat = !!sym(NamexTreat), yResp = !!sym(NameyResp)) %>% filter(!is.na(yResp)) %>%
    filter(!!sym(NameLevel) %in% (stat_data[[NameLevel]] %>% unique())) %>%
    mutate(!!sym(NameLevel) := factor(!!sym(NameLevel)))
  OutputName1 <- paste0("Output/Raincloud-", NameLevel, "-", NameyResp, "-", NamexTreat, ".pdf")
  p1 <- create_raincloud_sep(
    dataforplot = dataforplot, y_title = y_title, sigtext = sigtext, dodge_width = 0.9,
    point_size = 2.5)
  # point_size = 2)
  cairo_pdf(OutputName1, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54)
  print(p1); dev.off()
  OutputName2 <- paste0("Output/Forest-", NameLevel, "-", NameyResp, "-", NamexTreat, ".pdf")
  p2 <- create_forest_sep(
    dataforplot = dataforplot, y_title = y_title, sigtext = sigtext,
    dodge_width = 0.7)
  # dodge_width = 0.9)
  cairo_pdf(OutputName2, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54)
  print(p2); dev.off()
}

### Year-specific statistical analysis -------------------------------------------------------------------
FactorAs <- c("Treat_All", "Treat_Is")
# FactorB <- "Year"
FactorB <- "Year2"
formula_Randoms <- c("Site", "Crop", "PlotID")

for (factor in FactorAs) {
  DataForAnal <- dataTotal_202125MonoMaturity
  DataForAnal[[FactorB]] <- factor(
    as.character(DataForAnal[[FactorB]]),
    # levels = c("2021", "2022", "2023", "2024", "2025")) # Year
    levels = c("1", "2", "3", "4", "5")) # Year2
  FactorA <- factor
  ideal_A_order <- switch(
    FactorA,
    "Treat_All" = c("CK", "PP-P", "PP-F", "PLA-P", "PLA-F"),
    "Treat_Is"  = c("Without", "With"))
  existing_A <- ideal_A_order[ideal_A_order %in% unique(DataForAnal[[FactorA]])]
  DataForAnal[[FactorA]] <- factor(DataForAnal[[FactorA]], levels = existing_A)
  factorA_levels <- existing_A
  ideal_B_order <- levels(DataForAnal[[FactorB]])
  existing_B <- ideal_B_order[ideal_B_order %in% unique(DataForAnal[[FactorB]])]
  DataForAnal[[FactorB]] <- factor(DataForAnal[[FactorB]], levels = existing_B)
  factorB_levels <- existing_B
  workbook_list <- list()
  for (factorb in factorB_levels) {
    DataForAnal2 <- DataForAnal %>% filter(get(FactorB) == factorb)
    plan(multisession, workers = availableCores() - 1)
    results_list <- future_map(variables, function(var) {
      if (all(is.na(DataForAnal2[[var]]))) return(NULL)
      analyze_variable_one(
        data = DataForAnal2, response_var = var,
        formula_effects = FactorA, random_factors = formula_Randoms)},
      .progress = TRUE, .options = furrr_options(seed = TRUE)) %>% set_names(variables)
    results_list <- results_list %>% compact()
    plan(sequential)
    all_tables <- map(0:3, function(i) {
      table_key <- paste0("table", i)
      map_dfr(results_list, ~ .x[[table_key]])
    }) %>% set_names(paste0("table", 0:3))
    three_line_A_original <- create_three_line(all_tables$table3 %>% filter(
      Variable %in% variables_original), FactorA, factorA_levels)
    wb <- createWorkbook()
    sheet_names_original <-
      c("ModelEval", "FixedANOVA", "RandomTest", "Tukey")
    walk2(seq_along(sheet_names_original), sheet_names_original, ~ {
      table_data <- all_tables[[.x]] %>%
        filter(Variable %in% variables_original) %>%
        arrange(match(Variable, variables_original))
      if (nrow(table_data) > 0) {
        addWorksheet(wb, .y)
        writeData(wb, .y, table_data)}})
    addWorksheet(wb, "Line3")
    if (!is.null(three_line_A_original)) {writeData(wb, "Line3", three_line_A_original)}
    workbook_list[[factorb]] <- wb
  }
  wb_combined <- createWorkbook()
  all_sheet_names <- unique(unlist(lapply(workbook_list, names)))
  all_sheet_names <- all_sheet_names[!sapply(all_sheet_names, is.null)]
  for (sheet_name in all_sheet_names) {
    all_sheets <- list()
    factorb_names <- c()
    all_columns <- c()
    for (factorb in existing_B) {
      if (factorb %in% names(workbook_list)) {
        wb <- workbook_list[[factorb]]
        if (sheet_name %in% names(wb)) {
          sheet_data <- readWorkbook(wb, sheet = sheet_name, colNames = TRUE)
          if (nrow(sheet_data) > 0) {
            sheet_data_with_factor <- data.frame(
              temp_col = rep(factorb, nrow(sheet_data)),
              sheet_data, stringsAsFactors = FALSE, check.names = FALSE)
            colnames(sheet_data_with_factor)[1] <- FactorB
            all_sheets[[factorb]] <- sheet_data_with_factor
            factorb_names <- c(factorb_names, factorb)
            all_columns <- union(all_columns, colnames(sheet_data_with_factor))
          }}}}
    if (length(all_sheets) == 0) next
    combined_data_list <- list()
    for (i in seq_along(all_sheets)) {
      factorb <- factorb_names[i]
      sheet_data <- all_sheets[[factorb]]
      missing_cols <- setdiff(all_columns, colnames(sheet_data))
      if (length(missing_cols) > 0) {
        for (col in missing_cols) {sheet_data[[col]] <- NA}}
      sheet_data <- sheet_data[, all_columns, drop = FALSE]
      combined_data_list[[i]] <- sheet_data}
    combined_data <- do.call(rbind, combined_data_list)
    addWorksheet(wb_combined, sheet_name)
    writeData(wb_combined, sheet_name, combined_data)
  }
  output_filename <- paste0("Output/Sta_202125_Sep_", FactorB, "_", FactorA, ".xlsx")
  saveWorkbook(wb_combined, output_filename, overwrite = TRUE)
}

#### Year-specific statistical visualization ---------------------------------------------------------------------
dataSta99Year299Treat_Is <- read.xlsx(
  "output/Sta_202125_Sep_Year2_Treat_Is.xlsx",
  "FixedANOVA", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  filter(Term != "(Intercept)") %>% mutate(Sigused = SigF)
dataSta99Year299Treat_All <- read.xlsx(
  "output/Sta_202125_Sep_Year2_Treat_All.xlsx",
  "Tukey", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  mutate(Sigused = SignoA)
variables_original <-
  c("Yield", "SQI", "EMF")
resp_df <- data.frame(
  NameyResp = variables_original,
  y_title_text =
    c("Yield ~ (t ~ ha^-1)", "SQI", "EMF"), stringsAsFactors = FALSE)

NameLevel <- "Year2"
NamexTreat <- "Treat_Is"
# NamexTreat <- "Treat_All"
for (i in 1:nrow(resp_df)) {
  NameyResp <- resp_df$NameyResp[i]
  y_title <- parse(text = resp_df$y_title_text[i])
  stat_data <- get(paste0("dataSta99", NameLevel, "99", NamexTreat)) %>% filter(Variable == NameyResp)
  if (nrow(stat_data) == 0) {message("skip: ", NameyResp, " - stat_data is empty"); next}
  sigtext <- stat_data %>% pull(Sigused) %>% ifelse(is.na(.), "", .)
  dataforplot <- dataTotal_202125MonoMaturity %>% mutate(
    xTreat = !!sym(NamexTreat), yResp = !!sym(NameyResp)) %>% filter(!is.na(yResp)) %>%
    filter(!!sym(NameLevel) %in% (stat_data[[NameLevel]] %>% unique())) %>%
    mutate(!!sym(NameLevel) := factor(!!sym(NameLevel)))
  OutputName1 <- paste0("Output/Raincloud-", NameLevel, "-", NameyResp, "-", NamexTreat, ".pdf")
  p1 <- create_raincloud_sep(
    dataforplot = dataforplot, y_title = y_title, sigtext = sigtext, dodge_width = 0.9,
    point_size = 2.5)
  # point_size = 2)
  cairo_pdf(OutputName1, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54)
  print(p1); dev.off()
  OutputName2 <- paste0("Output/Forest-", NameLevel, "-", NameyResp, "-", NamexTreat, ".pdf")
  p2 <- create_forest_sep(
    dataforplot = dataforplot, y_title = y_title, sigtext = sigtext,
    dodge_width = 0.7)
  # dodge_width = 0.9)
  cairo_pdf(OutputName2, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54)
  print(p2); dev.off()
}

### Regression analysis - year -------------------------------------------------------------------
FactorAs <- c("Treat_All", "Treat_Is")
FactorB <- "Year2"
formula_Randoms <- c("Site", "Crop", "PlotID")

for (factor in FactorAs) {
  DataForAnal <- dataTotal_202125MonoMaturity
  FactorA <- factor
  ideal_A_order <- switch(
    FactorA, "Treat_All" = c("CK", "PP-P", "PP-F", "PLA-P", "PLA-F"),
    "Treat_Is"  = c("Without", "With"))
  existing_A <- ideal_A_order[ideal_A_order %in% unique(DataForAnal[[FactorA]])]
  DataForAnal[[FactorA]] <- factor(DataForAnal[[FactorA]], levels = existing_A)
  DataForAnal[[FactorB]] <- as.numeric(DataForAnal[[FactorB]])
  formula_Fix <- paste(FactorA, "*", FactorB)
  plan(multisession, workers = availableCores() - 1)
  results_list <- future_map(variables, function(var) {
    if (all(is.na(DataForAnal[[var]]))) return(NULL)
    analyze_variable_tworeg(
      data = DataForAnal, response_var = var,
      formula_effects = formula_Fix, random_factors = formula_Randoms)},
    .progress = TRUE, .options = furrr_options(seed = TRUE)) %>% set_names(variables)
  results_list <- results_list %>% compact()
  plan(sequential)
  all_tables <- map(0:4, function(i) {
    table_key <- paste0("table", i)
    map_dfr(results_list, ~ .x[[table_key]])
  }) %>% set_names(paste0("table", 0:4))
  wb <- createWorkbook()
  sheet_names_original <-
    c("ModelEval", "FixedANOVA", "RandomTest", "RegTotal", "RegIndi")
  walk2(seq_along(sheet_names_original), sheet_names_original, ~ {
    table_data <- all_tables[[.x]] %>%
      filter(Variable %in% variables_original) %>%
      arrange(match(Variable, variables_original))
    if (nrow(table_data) > 0) {
      addWorksheet(wb, .y)
      writeData(wb, .y, table_data)}})
  output_filename <- paste0("Output/Sta_202125_Two_", FactorB, "_", FactorA, ".xlsx")
  saveWorkbook(wb, output_filename, overwrite = TRUE)
}

#### Regression analysis - year visualization -------------------------------------------------------------------
dataSta99Year299Treat_Is <- read.xlsx(
  "output/Sta_202125_Two_Year2_Treat_Is.xlsx",
  "RegIndi", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  mutate(Sigused1 = SloPSig, Sigused2 = SloSignoA)
dataSta99Year299Treat_All <- read.xlsx(
  "output/Sta_202125_Two_Year2_Treat_All.xlsx",
  "RegIndi", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  mutate(Sigused1 = SloPSig, Sigused2 = SloSignoA)
variables_original <-
  c("Yield", "SQI", "EMF")
resp_df <- data.frame(
  NameyResp = variables_original,
  y_title_text =
    c("Yield ~ (t ~ ha^-1)", "SQI", "EMF"), stringsAsFactors = FALSE)

NameLevel <- "Year2"
x_title <- "Duration (a)"
NamexTreat <- "Treat_Is"
# NamexTreat <- "Treat_All"
for (i in 1:nrow(resp_df)) {
  NameyResp <- resp_df$NameyResp[i]
  y_title <- parse(text = resp_df$y_title_text[i])
  stat_data <- get(paste0("dataSta99", NameLevel, "99", NamexTreat)) %>% filter(Variable == NameyResp) %>% mutate(
    xTreat = !!sym(NamexTreat))
  if (nrow(stat_data) == 0) {message("skip: ", NameyResp, " - stat_data is empty"); next}
  dataforplot <- dataTotal_202125MonoMaturity %>% mutate(
    xTreat = !!sym(NamexTreat), yResp = !!sym(NameyResp)) %>% filter(!is.na(yResp))
  OutputNameR <- paste0("Output/Regression-", NameLevel, "-", NameyResp, "-", NamexTreat, ".pdf")
  pR <- create_regression(
    dataforplot = dataforplot, stat_data = stat_data, y_title = y_title, x_title = x_title)
  cairo_pdf(OutputNameR, bg = "transparent", width = 34.32/2.54, height = 5.44/2.54)
  print(pR); dev.off()
}


# map Location ------------------------------------------------------------
setwd('D:/R') # set current working directory
rm(list = ls()) # clear the environment
cat("\014") # clear the console
pacman::p_load(openxlsx, tidyverse, sf, terra, tidyterra)

BaseMap <- rast("D:/ArcGIS/GlobalAIPET/ai_v31_yr.tif") / 10000
points_df <- read.xlsx("Input/FieldData.xlsx", "PointLL")
pts <- vect(points_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
xrange <- diff(range(points_df$Longitude))
yrange <- diff(range(points_df$Latitude))
plot_ext <- ext(
  min(points_df$Longitude) - xrange * 0.2,
  max(points_df$Longitude) + xrange * 0.3,
  min(points_df$Latitude)  - yrange * 1.4,
  max(points_df$Latitude)  + yrange * 0.6)

BaseMap_sub <- crop(BaseMap, plot_ext)
BaseMap_sub <- classify(BaseMap_sub, cbind(0, NA))
ggplot() +
  geom_spatraster(data = BaseMap_sub, maxcell = 2e7) +
  scale_fill_gradientn(
    name = NULL, na.value = NA,
    guide = guide_colourbar(barheight = grid::unit(3, "cm")),
    colours = c("#BFA383", "#DDCBB2", "#F4ECDE", "#EAF3EE", "#C0DCD3", "#97C9BC"),
    limits = c(0, 1.5),
    oob = scales::squish,
    values = scales::rescale(c(0, 0.03, 0.2, 0.5, 0.65, 1.5)),
    breaks = c(0, 0.03, 0.2, 0.5, 0.65, 1.5),
    labels = c("", "0.03", "0.20", "0.50", "0.65", "1.5+")) +
  geom_spatvector(data = pts, color = "red", size = 2.5) +
  geom_text(data = points_df, aes(x = Longitude, y = Latitude, label = Site),
            nudge_y = 0.3, size = 16, size.unit = "pt") +
  coord_sf(xlim = c(xmin(plot_ext), xmax(plot_ext)), ylim = c(ymin(plot_ext), ymax(plot_ext)), expand = FALSE) +
  scale_x_continuous("Longitude (°)", breaks = c(115, 120), labels = c("115", "120"), expand = expansion(mult = 0)) +
  scale_y_continuous("Latitude (°)", breaks = c(38, 42), labels = c("38", "42"), expand = expansion(mult = 0)) +
  theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 16),
    plot.margin = margin(0, 0, 0, 0),
    legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
    legend.box = "vertical",
    legend.background = element_rect(fill = scales::alpha("grey90", 0.4), colour = NA),
    legend.key = element_blank(), legend.title = element_text(size = 14),
    legend.text = element_text(size = 11)) -> p_Map0
pdf("Output/MapLocation.pdf", height = 11.5/2.54)
p_Map0
dev.off()

BaseMap_admin_sub <- crop(BaseMap, ext(72, 136, 3, 54))
BaseMap_admin_sub <- classify(BaseMap_admin_sub, cbind(0, NA))
world <- st_read("D:/ArcGIS/GlobalCountry/global_all_country.shp", quiet = TRUE)
tenline <- st_read("D:/ArcGIS/CTAmap/十段线.shp", quiet = TRUE)
world <- st_transform(world, 4326)
tenline <- st_transform(tenline, 4326)
cn <- world %>% filter(CNTRY_NAME == "China")
plot_box <- as.data.frame(as.list(plot_ext))
ggplot() +
  geom_spatraster(data = BaseMap_admin_sub, maxcell = 2e7) +
  scale_fill_gradientn(
    name = NULL, na.value = NA,
    guide = guide_colourbar(barheight = grid::unit(3, "cm")),
    colours = c("#BFA383", "#DDCBB2", "#F4ECDE", "#EAF3EE", "#C0DCD3", "#97C9BC"),
    limits = c(0, 1.5),
    oob = scales::squish,
    values = scales::rescale(c(0, 0.03, 0.2, 0.5, 0.65, 1.5)),
    breaks = c(0, 0.03, 0.2, 0.5, 0.65, 1.5),
    labels = c("", "0.03", "0.20", "0.50", "0.65", "1.5+")) +
  geom_sf(data = world, fill = NA, color = "#C7C7C7", linewidth = 0.1) +
  geom_sf(data = cn, fill = NA, color = "#7A7A7A", linewidth = 0.25) +
  geom_sf(data = tenline, color = "#7A7A7A", linewidth = 0.25) +
  geom_rect(data = plot_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "#8c510a", linewidth = 0.6, linetype = "solid") +
  coord_sf(xlim = c(72, 136), ylim = c(3, 54),
           expand = FALSE, crs = st_crs(4326), datum = st_crs(4326)) +
  scale_x_continuous("Longitude (°)", breaks = c(80, 120), labels = c("80", "120"), expand = expansion(mult = 0)) +
  scale_y_continuous("Latitude (°)", breaks = c(10, 40), labels = c("10", "40"), expand = expansion(mult = 0)) +
  theme_classic() + theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.text = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 16),
    plot.margin = margin(0, 0, 0, 0),
    legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
    legend.box = "vertical",
    legend.background = element_rect(fill = scales::alpha("grey90", 0.4), colour = NA),
    legend.key = element_blank(), legend.title = element_text(size = 14),
    legend.text = element_text(size = 11)
  ) -> p_Map00
pdf("Output/MapRegion.pdf", height = 11.5/2.54)
p_Map00
dev.off()

# temperature and precipitation -----------------------------------------------------------
setwd('D:/R') # set current working directory
rm(list = ls()) # clear the environment
cat("\014") # clear the console
pacman::p_load(openxlsx, tidyverse)

res_out <- read.xlsx("Input/FieldData.xlsx", "monthly_climate") %>% as_tibble()
clim_keep <- res_out %>% filter(
  (Site == "Chifeng"  & Year %in% 2021:2022) |
    (Site == "Ulanqab"  & Year %in% 2021:2023) |
    (Site %in% c("Zhangbei", "Youyu") & Year %in% 2021:2025))
plot_dat <- clim_keep %>% mutate(
  Date = as.Date(sprintf("%d-%02d-01", Year, Month)),
  Site = factor(Site, levels = c("Zhangbei", "Ulanqab", "Youyu", "Chifeng")
  )) %>% arrange(Site, Date)
temp_min <- floor(min(plot_dat$t_mon, na.rm = TRUE)) - 7
temp_max <- ceiling(max(plot_dat$t_mon, na.rm = TRUE)) + 2
precip_max <- ceiling(max(plot_dat$p_mon, na.rm = TRUE) / 10) * 10
scale_factor <- precip_max / (temp_max - temp_min)
plot_dat <- plot_dat %>% mutate(p_plot = p_mon / scale_factor + temp_min)
shade_dat <- plot_dat %>% distinct(Year) %>% mutate(
  xmin = as.Date(sprintf("%d-04-16", Year)), xmax = as.Date(sprintf("%d-09-15", Year)))

col_temp <- "#F6C188"
col_temp_line <- "#E97800"
col_precip_fill <- "#95D0D9"
col_precip_border <- "#068CA3"
ggplot(plot_dat, aes(x = Date)) +
  geom_rect(data = shade_dat, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#F3F3F3", color = NA) +
  geom_tile(aes(y = (temp_min + p_plot) / 2, height = p_plot - temp_min),
            width = 25, fill = col_precip_fill, color = col_precip_border, linewidth = 0.2) +
  geom_line(aes(y = t_mon), color = col_temp_line, linewidth = 0.5) +
  geom_point(aes(y = t_mon), color = col_temp_line, size = 0.8) +
  facet_grid(Site ~ ., switch = "y") + scale_x_date(
    breaks = seq(as.Date("2021-05-01"), max(plot_dat$Date), by = "1 year"),
    date_labels = "%Y/%m", expand = c(0.01, 0.01)) + scale_y_continuous(
      name = "Temperature (°C)", limits = c(temp_min, temp_max),
      breaks = c(-20, 0, 20), expand = expansion(mult = c(0, 0.03)),
      sec.axis = sec_axis(~ (. - temp_min) * scale_factor, name = "Precipitation (mm)")) +
  theme_classic() + theme(
    axis.line = element_blank(), panel.border = element_rect(),
    plot.background = element_blank(), panel.background = element_blank(),
    axis.title = element_text(size = 16), axis.title.x = element_blank(),
    strip.text.y.left = element_blank(), axis.text.x = element_text(size = 13),
    axis.text.y = element_text(color = col_temp_line, size = 14),
    axis.text.y.right = element_text(color = col_precip_border, size = 14),
    legend.position = 'none', plot.margin = margin(0, 0, 0, 0),
    aspect.ratio = 7/24) -> p_clim; p_clim
pdf("Output/Climate.pdf", height = 10.9/2.54)
p_clim
dev.off()
