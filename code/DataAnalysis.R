# Setup ----
rm(list = ls())
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(openxlsx, tidyverse, lmerTest, nortest, car, emmeans, multcomp, multcompView, DHARMa, performance, furrr, smplot2, sdamr, GGally, patchwork, sf, terra, tidyterra)
options(contrasts = c("contr.sum", "contr.poly"))
input_dir <- if (dir.exists("Input")) "Input" else "."
field_file <- file.path(input_dir, "FieldData.xlsx")
if (!file.exists(field_file)) stop("FieldData.xlsx not found. Put it in the working directory or Input/.")
dir.create("Output", showWarnings = FALSE, recursive = TRUE)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

# Custom functions ----
check_normality <- function(data) {
  n <- length(data)
  if (n <= 50) {
    test_result <- shapiro.test(data)
    method <- "Shapiro-Wilk"} else {
      test_result <- lillie.test(data)
      method <- "Lilliefors"}
  list(method = method, p.value = test_result$p.value, normality = test_result$p.value > 0.05)}

analyze_variable_two <- function(
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

    normality_result <- check_normality(residuals_lmer)
    if (length(fixed_effects_terms) == 1) {
      interaction_term <- data_clean[[fixed_effects_terms]]} else {
        interaction_term <- interaction(
          data_clean[[fixed_effects_terms[1]]], data_clean[[fixed_effects_terms[2]]])}
    interaction_term_dropped <- droplevels(interaction_term)
    sim_res <- simulateResiduals(model_lmer)
    levene_result <- testCategorical(sim_res, interaction_term_dropped)
    table0 <- cbind(base_info, data.frame(
      NormalityMethod = normality_result$method,
      NormalityP = normality_result$p.value,
      IsNormal = normality_result$normality,
      LeveneP = levene_result$homogeneity[1, "Pr(>F)"],
      Homogeneous = levene_result$homogeneity[1, "Pr(>F)"] > 0.05))

    loglik <- logLik(model_lmer)
    mse <- mean(residuals_lmer^2)
    table1 <- cbind(base_info, data.frame(
      UsedFormula = gsub("\\s+", " ", paste(deparse(used_formula), collapse = "")),
      AIC = AIC(model_lmer), BIC = BIC(model_lmer),
      logLik = as.numeric(loglik), npar = attr(loglik, "df"),
      Deviance = -2 * as.numeric(loglik), DfResid = df.residual(model_lmer),
      R2marginal <- r2(model_lmer, tolerance = 1e-1000)$R2_marginal,
      R2conditional <- r2(model_lmer, tolerance = 1e-1000)$R2_conditional,
      MSE = mse, RMSE = sqrt(mse)))
    colnames(table1) <- c(
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
    table2 <- cbind(base_info, data.frame(
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
    table3 <- cbind(base_info, data.frame(
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
    table4 <- data.frame(
      base_info,
      summary_stats_A[[fixed_effects_terms[1]]],
      summary_stats_A[, c("Mean", "SD", "SE", "n", "lower", "upper")],
      Sig = cld_A_sorted$.group, Sig_noA = cld_A_sorted$Sig_noA,
      cld_A_sorted[, c("emmean", "SE", "df", "lower.CL", "upper.CL")])
    colnames(table4) <- c(
      names(base_info), fixed_effects_terms[1], "Mean", "SD", "SE", "n",
      "CLlower", "CLupper", "Sig", "SignoA",
      "emmean", "SEemm", "dfemm", "CLloweremm", "CLupperemm")

    summary_stats_B <- data_clean %>%
      group_by(!!sym(fixed_effects_terms[2]), .drop = FALSE) %>% summarise(
        Mean = mean(.data[[response_var]], na.rm = TRUE), SD = sd(.data[[response_var]], na.rm = TRUE),
        n = sum(!is.na(.data[[response_var]])), SE = ifelse(n > 0, SD / sqrt(n), NA_real_),
        lower = ifelse(n > 0, Mean - qt(0.975, n - 1) * SE, NA_real_),
        upper = ifelse(n > 0, Mean + qt(0.975, n - 1) * SE, NA_real_), .groups = "drop")
    emm_B <- suppressMessages(emmeans(model_lmer, reformulate(fixed_effects_terms[2])))
    cld_B_results <- suppressMessages(
      multcomp::cld(emm_B, adjust = "tukey", Letters = letters, sort = TRUE, reverse = TRUE))
    cld_B_results$.group <- gsub(" ", "", cld_B_results$.group)
    cld_B_results$Sig_noA <- if (all(cld_B_results$.group == "a")) "" else cld_B_results$.group
    all_levels_B <- data.frame(x = levels(data_clean[[fixed_effects_terms[2]]]))
    names(all_levels_B) <- fixed_effects_terms[2]
    cld_B_resultsNew <- merge(
      all_levels_B, cld_B_results, by = fixed_effects_terms[2], all.x = TRUE)
    cld_B_sorted <- cld_B_resultsNew[match(levels(data_clean[[fixed_effects_terms[2]]]), cld_B_resultsNew[[fixed_effects_terms[2]]]), ]
    table5 <- data.frame(
      base_info,
      summary_stats_B[[fixed_effects_terms[2]]],
      summary_stats_B[, c("Mean", "SD", "SE", "n", "lower", "upper")],
      Sig = cld_B_sorted$.group, Sig_noA = cld_B_sorted$Sig_noA,
      cld_B_sorted[, c("emmean", "SE", "df", "lower.CL", "upper.CL")])
    colnames(table5) <- c(
      names(base_info), fixed_effects_terms[2], "Mean", "SD", "SE", "n",
      "CLlower", "CLupper", "Sig", "SignoA",
      "emmean", "SEemm", "dfemm", "CLloweremm", "CLupperemm")

    summary_stats_A_B <- data_clean %>% group_by(
      !!sym(fixed_effects_terms[1]), !!sym(fixed_effects_terms[2]), .drop = FALSE) %>% summarise(
        Mean = mean(.data[[response_var]], na.rm = TRUE),
        SD = sd(.data[[response_var]], na.rm = TRUE),
        n = sum(!is.na(.data[[response_var]])), SE = ifelse(n > 0, SD / sqrt(n), NA_real_),
        lower = ifelse(n > 0, Mean - qt(0.975, n - 1) * SE, NA_real_),
        upper = ifelse(n > 0, Mean + qt(0.975, n - 1) * SE, NA_real_), .groups = "drop")
    emm_A_B <- suppressMessages(emmeans(
      model_lmer, reformulate(paste(fixed_effects_terms[2], "|", fixed_effects_terms[1]))))
    cld_A_B_results <- suppressMessages(
      multcomp::cld(emm_A_B, adjust = "tukey", Letters = letters, sort = TRUE, reverse = TRUE))
    cld_A_B_results$.group <- gsub(" ", "", cld_A_B_results$.group)
    cld_A_B_results$Sig_noA <- cld_A_B_results$.group
    cld_A_B_results <- cld_A_B_results %>% group_by(.data[[fixed_effects_terms[1]]]) %>%
      mutate(Sig_noA = if (all(.group == "a")) "" else .group) %>% ungroup()
    cld_A_B_resultsNew <- cld_A_B_results %>% complete(
      !!sym(fixed_effects_terms[2]) := factor(levels(data_clean[[fixed_effects_terms[2]]]), levels = levels(data_clean[[fixed_effects_terms[2]]])),
      !!sym(fixed_effects_terms[1]) := factor(levels(data_clean[[fixed_effects_terms[1]]]), levels = levels(data_clean[[fixed_effects_terms[1]]])),
      fill = list(.group = "", Sig_noA = "")) %>% arrange(!!sym(fixed_effects_terms[1]), !!sym(fixed_effects_terms[2]))
    cld_A_B_sorted <- cld_A_B_resultsNew %>% arrange(
      .data[[fixed_effects_terms[1]]], match(.data[[fixed_effects_terms[2]]], levels(data_clean[[fixed_effects_terms[2]]])))
    table6 <- data.frame(
      base_info,
      summary_stats_A_B[[fixed_effects_terms[1]]],
      summary_stats_A_B[[fixed_effects_terms[2]]],
      summary_stats_A_B[, c("Mean", "SD", "SE", "n", "lower", "upper")],
      Sig = cld_A_B_sorted$.group, Sig_noA = cld_A_B_sorted$Sig_noA,
      cld_A_B_sorted[, c("emmean", "SE", "df", "lower.CL", "upper.CL")])
    colnames(table6) <- c(
      names(base_info), fixed_effects_terms[1], fixed_effects_terms[2], "Mean", "SD", "SE", "n",
      "CLlower", "CLupper", "Sig", "SignoA",
      "emmean", "SEemm", "dfemm", "CLloweremm", "CLupperemm")

    summary_stats_B_A <- data_clean %>% group_by(
      !!sym(fixed_effects_terms[2]), !!sym(fixed_effects_terms[1]), .drop = FALSE) %>% summarise(
        Mean = mean(.data[[response_var]], na.rm = TRUE), SD = sd(.data[[response_var]], na.rm = TRUE),
        n = sum(!is.na(.data[[response_var]])), SE = ifelse(n > 0, SD / sqrt(n), NA_real_),
        lower = ifelse(n > 0, Mean - qt(0.975, n - 1) * SE, NA_real_),
        upper = ifelse(n > 0, Mean + qt(0.975, n - 1) * SE, NA_real_), .groups = "drop")
    emm_B_A <- suppressMessages(emmeans(
      model_lmer, reformulate(paste(fixed_effects_terms[1], "|", fixed_effects_terms[2]))))
    cld_B_A_results <- suppressMessages(
      multcomp::cld(emm_B_A, adjust = "tukey", Letters = letters, sort = TRUE, reverse = TRUE))
    cld_B_A_results$.group <- gsub(" ", "", cld_B_A_results$.group)
    cld_B_A_results$Sig_noA <- cld_B_A_results$.group
    cld_B_A_results <- cld_B_A_results %>% group_by(.data[[fixed_effects_terms[2]]]) %>%
      mutate(Sig_noA = if (all(.group == "a")) "" else .group) %>% ungroup()
    cld_B_A_resultsNew <- cld_B_A_results %>% complete(
      !!sym(fixed_effects_terms[1]) := factor(levels(data_clean[[fixed_effects_terms[1]]]), levels = levels(data_clean[[fixed_effects_terms[1]]])),
      !!sym(fixed_effects_terms[2]) := factor(levels(data_clean[[fixed_effects_terms[2]]]), levels = levels(data_clean[[fixed_effects_terms[2]]])),
      fill = list(.group = "", Sig_noA = "")) %>% arrange(!!sym(fixed_effects_terms[2]), !!sym(fixed_effects_terms[1]))
    cld_B_A_sorted <- cld_B_A_resultsNew %>% arrange(
      .data[[fixed_effects_terms[2]]], match(.data[[fixed_effects_terms[1]]], levels(data_clean[[fixed_effects_terms[1]]])))
    table7 <- data.frame(
      base_info,
      summary_stats_B_A[[fixed_effects_terms[2]]],
      summary_stats_B_A[[fixed_effects_terms[1]]],
      summary_stats_B_A[, c("Mean", "SD", "SE", "n", "lower", "upper")],
      Sig = cld_B_A_sorted$.group, Sig_noA = cld_B_A_sorted$Sig_noA,
      cld_B_A_sorted[, c("emmean", "SE", "df", "lower.CL", "upper.CL")])
    colnames(table7) <- c(
      names(base_info), fixed_effects_terms[2], fixed_effects_terms[1], "Mean", "SD", "SE", "n",
      "CLlower", "CLupper", "Sig", "SignoA",
      "emmean", "SEemm", "dfemm", "CLloweremm", "CLupperemm")

    tables <- list(
      table0 = table0, table1 = table1, table2 = table2, table3 = table3,
      table4 = table4, table5 = table5, table6 = table6, table7 = table7)
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

analyze_variable_tworeg <- function(
  data, response_var, formula_effects, random_factors) {
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

    normality_result <- check_normality(residuals_lmer)
    if (length(fixed_effects_terms) == 1) {
      interaction_term <- data_clean[[fixed_effects_terms]]} else {
        interaction_term <- interaction(
          data_clean[[fixed_effects_terms[1]]], data_clean[[fixed_effects_terms[2]]])}
    interaction_term_dropped <- droplevels(interaction_term)
    sim_res <- simulateResiduals(model_lmer)
    levene_result <- testCategorical(sim_res, interaction_term_dropped)
    table0 <- cbind(base_info, data.frame(
      NormalityMethod = normality_result$method,
      NormalityP = normality_result$p.value,
      IsNormal = normality_result$normality,
      LeveneP = levene_result$homogeneity[1, "Pr(>F)"],
      Homogeneous = levene_result$homogeneity[1, "Pr(>F)"] > 0.05))

    loglik <- logLik(model_lmer)
    mse <- mean(residuals_lmer^2)
    table1 <- cbind(base_info, data.frame(
      UsedFormula = gsub("\\s+", " ", paste(deparse(used_formula), collapse = "")),
      AIC = AIC(model_lmer), BIC = BIC(model_lmer),
      logLik = as.numeric(loglik), npar = attr(loglik, "df"),
      Deviance = -2 * as.numeric(loglik), DfResid = df.residual(model_lmer),
      R2marginal <- r2(model_lmer, tolerance = 1e-1000)$R2_marginal,
      R2conditional <- r2(model_lmer, tolerance = 1e-1000)$R2_conditional,
      MSE = mse, RMSE = sqrt(mse)))
    colnames(table1) <- c(
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
    table2 <- cbind(base_info, data.frame(
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
    table3 <- cbind(base_info, data.frame(
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
    table4 <- data.frame(base_info, overall_slope_df[, 2:9], overall_intercept_df[, 2:8])
    colnames(table4) <- c(
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
    table5 <- data.frame(base_info, slope_df_sorted[, 1:11], intercepts_df_sorted[, 2:8])
    colnames(table5) <- c(
      names(base_info), fixed_effects_terms[1],
      "Slope", "SloSE", "SloDf", "SloCLlower", "SloCLupper", "SloSig", "SloSignoA", "SloTratio", "SloP", "SloPSig",
      "Intercept", "IntSE", "IntDf", "IntCLlower", "IntCLupper", "IntTratio", "IntP")

    tables <- list(
      table0 = table0, table1 = table1, table2 = table2,
      table3 = table3, table4 = table4, table5 = table5)
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

    normality_result <- check_normality(residuals_lmer)
    if (length(fixed_effects_terms) == 1) {
      interaction_term <- data_clean[[fixed_effects_terms]]} else {
        interaction_term <- interaction(
          data_clean[[fixed_effects_terms[1]]], data_clean[[fixed_effects_terms[2]]])}
    interaction_term_dropped <- droplevels(interaction_term)
    sim_res <- simulateResiduals(model_lmer)
    levene_result <- testCategorical(sim_res, interaction_term_dropped)
    table0 <- cbind(base_info, data.frame(
      NormalityMethod = normality_result$method,
      NormalityP = normality_result$p.value,
      IsNormal = normality_result$normality,
      LeveneP = levene_result$homogeneity[1, "Pr(>F)"],
      Homogeneous = levene_result$homogeneity[1, "Pr(>F)"] > 0.05))

    loglik <- logLik(model_lmer)
    mse <- mean(residuals_lmer^2)
    table1 <- cbind(base_info, data.frame(
      UsedFormula = gsub("\\s+", " ", paste(deparse(used_formula), collapse = "")),
      AIC = AIC(model_lmer), BIC = BIC(model_lmer),
      logLik = as.numeric(loglik), npar = attr(loglik, "df"),
      Deviance = -2 * as.numeric(loglik), DfResid = df.residual(model_lmer),
      R2marginal <- r2(model_lmer, tolerance = 1e-1000)$R2_marginal,
      R2conditional <- r2(model_lmer, tolerance = 1e-1000)$R2_conditional,
      MSE = mse, RMSE = sqrt(mse)))
    colnames(table1) <- c(
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
    table2 <- cbind(base_info, data.frame(
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
    table3 <- cbind(base_info, data.frame(
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
    table4 <- data.frame(
      base_info,
      summary_stats_A[[fixed_effects_terms[1]]],
      summary_stats_A[, c("Mean", "SD", "SE", "n", "lower", "upper")],
      Sig = cld_A_sorted$.group, Sig_noA = cld_A_sorted$Sig_noA,
      cld_A_sorted[, c("emmean", "SE", "df", "lower.CL", "upper.CL")])
    colnames(table4) <- c(
      names(base_info), fixed_effects_terms[1], "Mean", "SD", "SE", "n",
      "CLlower", "CLupper", "Sig", "SignoA",
      "emmean", "SEemm", "dfemm", "CLloweremm", "CLupperemm")

    tables <- list(table0 = table0, table1 = table1, table2 = table2, table3 = table3, table4 = table4)
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
create_three_line2 <- function(
    tukey_table, factor1_var, factor1_levels, factor2_var, factor2_levels) {
  if (is.null(tukey_table)) return(NULL)
  three_line_table <- tukey_table %>% mutate(
    Mean = sprintf("%.2f", as.numeric(Mean)), SE = sprintf("%.2f", as.numeric(SE)),
    value = if_else(
      is.na(SignoA) | SignoA == "", paste0(Mean, " ± ", SE), paste0(Mean, " ± ", SE, " ", SignoA))) %>%
    dplyr::select(all_of(factor1_var), all_of(factor2_var), Variable, value) %>%
    pivot_wider(names_from = Variable, values_from = value) %>%
    arrange(factor(!!sym(factor1_var), levels = factor1_levels),
            factor(!!sym(factor2_var), levels = factor2_levels))
  return(three_line_table)}

analyze_variable_lm <- function(
  data, formula, response_var, site, year, crop, pattern, period) {
  safe_analyze <- safely(function() {
    model_lm <- lm(formula, data = data)
    residuals_lm <- residuals(model_lm)
    summary_lm <- summary(model_lm)
    base_info <- data.frame(
      Site = site, Year2 = year, Crop = crop,
      Pattern = pattern, Period = period, Variable = response_var)

    normality_result <- check_normality(residuals_lm)
    levene_result <- leveneTest(model_lm)
    table0 <- cbind(base_info, data.frame(
      NormalityMethod = normality_result$method,
      NormalityP = normality_result$p.value,
      IsNormal = normality_result$normality,
      LeveneP = levene_result$`Pr(>F)`[1],
      Homogeneous = levene_result$`Pr(>F)`[1] > 0.05))

    loglik <- logLik(model_lm)
    mse <- mean(residuals_lm^2)
    table1 <- cbind(base_info, data.frame(
      AIC = AIC(model_lm), BIC = BIC(model_lm),
      logLik = as.numeric(loglik), npar = attr(loglik, "df"),
      Deviance = -2 * as.numeric(loglik),
      DfResid = df.residual(model_lm),
      R2 = summary_lm$r.squared,
      R2Adjusted = summary_lm$adj.r.squared,
      MSE = mse, RMSE = sqrt(mse)))

    anova_table <- car::Anova(model_lm, type = 3)
    anova_sig <- symnum(
      anova_table$`Pr(>F)`, corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.01, 0.05, 1), symbols = c("**", "*", ""))
    table2 <- cbind(base_info, data.frame(
      Term = rownames(anova_table), Df = anova_table$Df,
      SumSq = anova_table$`Sum Sq`,
      MeanSq = anova_table$`Sum Sq` / anova_table$`Df`,
      Fvalue = anova_table$`F value`, P = anova_table$`Pr(>F)`,
      Sig = as.character(anova_sig)))

    summary_stats <- data %>%
      group_by(across(all_of(formula_effects))) %>%
      summarise(
        Mean = mean(.data[[response_var]], na.rm = TRUE),
        SD = sd(.data[[response_var]], na.rm = TRUE),
        n = sum(!is.na(.data[[response_var]])), SE = SD / sqrt(n),
        lower = Mean - qt(0.975, n - 1) * SE,
        upper = Mean + qt(0.975, n - 1) * SE, .groups = "drop")
    emm <- emmeans(model_lm, reformulate(formula_effects))
    cld_results <- suppressMessages(
      multcomp::cld(emm, adjust = "tukey", Letters = letters, sort = TRUE, reverse = TRUE))
    cld_sorted <- cld_results[match(factor_levels, cld_results[[formula_effects]]), ]
    cld_sorted$.group <- gsub(" ", "", cld_sorted$.group)
    cld_sorted$Sig_noA <- if (all(cld_sorted$.group == "a")) "" else cld_sorted$.group
    table3 <- data.frame(
      base_info,
      summary_stats[[formula_effects]],
      summary_stats[, c("Mean", "SD", "SE", "n", "lower", "upper")],
      Sig = cld_sorted$.group, Sig_noA = cld_sorted$Sig_noA,
      cld_sorted[, c("emmean", "SE", "df", "lower.CL", "upper.CL")])
    colnames(table3) <- c(
      names(base_info), formula_effects, "Mean", "SD", "SE", "n",
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
      "Analysis failed for %s %s %s %s %s %s - Error: %s",
      site, year, crop, pattern, period, response_var, result$error$message))
    return(NULL)}
  return(result$result)
}

analyze_variable_lm_2025 <- function(
  data, formula, response_var, site, year, crop, pattern, period) {
  safe_analyze <- safely(function() {
    model_lm <- lm(formula, data = data)
    residuals_lm <- residuals(model_lm)
    summary_lm <- summary(model_lm)
    base_info <- data.frame(
      Site = site, Year = year, Crop = crop,
      Pattern = pattern, Period = period, Variable = response_var)

    normality_result <- check_normality(residuals_lm)
    levene_result <- leveneTest(model_lm)
    table0 <- cbind(base_info, data.frame(
      NormalityMethod = normality_result$method,
      NormalityP = normality_result$p.value,
      IsNormal = normality_result$normality,
      LeveneP = levene_result$`Pr(>F)`[1],
      Homogeneous = levene_result$`Pr(>F)`[1] > 0.05))

    loglik <- logLik(model_lm)
    mse <- mean(residuals_lm^2)
    table1 <- cbind(base_info, data.frame(
      AIC = AIC(model_lm), BIC = BIC(model_lm),
      logLik = as.numeric(loglik), npar = attr(loglik, "df"),
      Deviance = -2 * as.numeric(loglik),
      DfResid = df.residual(model_lm),
      R2 = summary_lm$r.squared,
      R2Adjusted = summary_lm$adj.r.squared,
      MSE = mse, RMSE = sqrt(mse)))

    anova_table <- car::Anova(model_lm, type = 3)
    anova_sig <- symnum(
      anova_table$`Pr(>F)`, corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.01, 0.05, 1), symbols = c("**", "*", ""))
    table2 <- cbind(base_info, data.frame(
      Term = rownames(anova_table), Df = anova_table$Df,
      SumSq = anova_table$`Sum Sq`,
      MeanSq = anova_table$`Sum Sq` / anova_table$`Df`,
      Fvalue = anova_table$`F value`, P = anova_table$`Pr(>F)`,
      Sig = as.character(anova_sig)))

    summary_stats <- data %>%
      group_by(across(all_of(formula_effects))) %>%
      summarise(
        Mean = mean(.data[[response_var]], na.rm = TRUE),
        SD = sd(.data[[response_var]], na.rm = TRUE),
        n = sum(!is.na(.data[[response_var]])), SE = SD / sqrt(n),
        lower = Mean - qt(0.975, n - 1) * SE,
        upper = Mean + qt(0.975, n - 1) * SE, .groups = "drop")
    emm <- emmeans(model_lm, reformulate(formula_effects))
    cld_results <- suppressMessages(
      multcomp::cld(emm, adjust = "tukey", Letters = letters, sort = TRUE, reverse = TRUE))
    cld_sorted <- cld_results[match(factor_levels, cld_results[[formula_effects]]), ]
    cld_sorted$.group <- gsub(" ", "", cld_sorted$.group)
    cld_sorted$Sig_noA <- if (all(cld_sorted$.group == "a")) "" else cld_sorted$.group
    table3 <- data.frame(
      base_info,
      summary_stats[[formula_effects]],
      summary_stats[, c("Mean", "SD", "SE", "n", "lower", "upper")],
      Sig = cld_sorted$.group, Sig_noA = cld_sorted$Sig_noA,
      cld_sorted[, c("emmean", "SE", "df", "lower.CL", "upper.CL")])
    colnames(table3) <- c(
      names(base_info), formula_effects, "Mean", "SD", "SE", "n",
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
      "Analysis failed for %s %s %s %s %s %s - Error: %s",
      site, year, crop, pattern, period, response_var, result$error$message))
    return(NULL)}
  return(result$result)
}

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
  if (!all(c("xTreat", "yResp") %in% names(dataforplot))) {stop("dataforplot must contain xTreat and yResp columns")}
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
  if (!all(c("xTreat", "yResp") %in% names(dataforplot))) {stop("dataforplot must contain xTreat and yResp columns")}
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
  if (!all(c("xTreat", "yResp") %in% names(dataforplot))) {stop("dataforplot must contain xTreat and yResp columns")}
  if (!is.factor(dataforplot$xTreat)) {dataforplot$xTreat <- as.factor(dataforplot$xTreat)}
  if (!as.character(NameLevel) %in% names(dataforplot)) {stop(paste0("dataforplot must contain '", NameLevel, "' column"))}
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
  if (!all(c("xTreat", "yResp") %in% names(dataforplot))) {stop("dataforplot must contain xTreat and yResp columns")}
  if (!is.factor(dataforplot$xTreat)) {dataforplot$xTreat <- as.factor(dataforplot$xTreat)}
  if (!as.character(NameLevel) %in% names(dataforplot)) {stop(paste0("dataforplot must contain '", NameLevel, "' column"))}
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
  if (!all(c("xTreat", "yResp") %in% names(dataforplot))) {stop("dataforplot must contain xTreat and yResp columns")}
  if (!is.factor(dataforplot$xTreat)) {dataforplot$xTreat <- as.factor(dataforplot$xTreat)}
  if (!as.character(NameLevel) %in% names(dataforplot)) {stop(paste0("dataforplot must contain '", NameLevel, "' column"))}
  required_stat_cols <- c("Slope", "Intercept", "xTreat", "Sigused1", "Sigused2")
  if (!all(required_stat_cols %in% names(stat_data))) {stop(paste0("stat_data must contain columns: ", paste(required_stat_cols, collapse = ", ")))}
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

# Load rebuilt field data ----
dataOrigin <- read.xlsx(field_file, "dataall(plain)") %>% mutate(
  Site = factor(Site, levels = c("Zhangbei", "Ulanqab", "Youyu", "Chifeng")),
  Year = as.numeric(Year),
  Crop = factor(Crop, levels = c("Oat", "Soybean")),
  Period = factor(Period, levels = c("Jointing", "FlowerBudDifferentiation", "Filling", "Maturity")),
  Pattern = factor(Pattern, levels = c("monoculture", "rotation", "intercropping")),
  Treat_All = factor(Treat_All, levels = c("CK", "PP-P", "PP-F", "PLA-P", "PLA-F")),
  Block = factor(Block, levels = c("1", "2", "3", "4")),
  Treat_Is = factor(if_else(Treat_All == "CK", "Without", "With"), levels = c("Without", "With")),
  Treat_Type = factor(str_extract(Treat_All, "CK|PP|PLA"), levels = c("CK", "PP", "PLA")),
  Treat_Shape = factor(case_when(Treat_All == "CK" ~ "CK", str_detect(Treat_All, "-P") ~ "Powder", str_detect(Treat_All, "-F") ~ "Fiber"), levels = c("CK", "Powder", "Fiber"))) %>%
  rename_with(~ str_replace(., "^(\\d)", "X\\1"))

# Fig. 1, Fig. S2 and Fig. S3: 2021-2025 monoculture maturity analysis ----
# Index calculation ----
dataOrigin_202125MonoMaturity <- dataOrigin %>% filter(
  Year %in% 2021:2025, Period == "Maturity", Pattern == "monoculture")

dataOrigin_202125MonoMaturity_clean <- dataOrigin_202125MonoMaturity %>% dplyr::select(where(~!all(is.na(.))))

dataOrigin_202125MonoMaturity_norm <- dataOrigin_202125MonoMaturity_clean %>%
  mutate(across(SeedWater:BDL60, ~ {
    if(all(is.na(.))) return(.)
    reverse_cols <-
      c("BD", "pH", "POX", "PER", "OX", "VectorLength", "CQI", "BDL40", "BDL60")
    col_name <- cur_column()
    if(col_name %in% reverse_cols) {(max(., na.rm = TRUE) - .) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))
    } else {(. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))}}, .names = "{.col}"))

data_for_pca_202125MonoMaturity <- dataOrigin_202125MonoMaturity_norm %>%
  dplyr::select(SWC:PER, -SOC, -TN, -SOCTNratio, -DOC, -DON, -DOCNratio) %>%
  filter(if_all(everything(), ~ !is.na(.) & !is.infinite(.)))
weights_202125MonoMaturity <- prcomp(data_for_pca_202125MonoMaturity, center = TRUE, scale. = TRUE) |>
  (\(pca) {kept <- pca$sdev^2 > 1
  contributions <- apply(pca$rotation[, kept], 1, \(u) sqrt(sum((u^2) * pca$sdev[kept]^2)))
  contributions / sum(contributions)})()

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
    dplyr::select(SWC:PER, -SOC, -TN, -SOCTNratio, -DOC, -DON, -DOCNratio) %>% as.matrix()
  sapply(1:nrow(norm_data), function(i) {
    row_data <- norm_data[i, ]
    valid <- !is.na(row_data)
    if (sum(valid) > 0) {
      available_weights <- weights_202125MonoMaturity[valid]
      adjusted_weights <- available_weights / sum(available_weights)
      sum(row_data[valid] * adjusted_weights)} else NA})
}) %>% bind_cols(
  dataOrigin_202125MonoMaturity_norm %>% dplyr::select(
    PlantProd, PlantGrow, PlantForBelow,
    SoilNutriTurnover, SoilActi, SoilPhy, SoilChem, SoilCSta,
    Support, Regulation, EMF))

dataTotal_202125MonoMaturity <- dataTotal_202125MonoMaturity %>%
  mutate(Year2 = case_when(
    Year == "2021" ~ 1, Year == "2022" ~ 2,
    Year == "2023" ~ 3, Year == "2024" ~ 4, Year == "2025" ~ 5))

# Statistical tests ----
pacman::p_load(openxlsx, tidyverse, lmerTest, nortest, car, emmeans, multcomp,
               multcompView, DHARMa, performance, furrr, smplot2, sdamr, GGally)

dataTotal_202125MonoMaturity <- dataTotal_202125MonoMaturity %>% mutate(
    Site = factor(Site, levels = c("Zhangbei", "Ulanqab", "Youyu", "Chifeng")),
    Year = as.numeric(Year),
    Year2 = as.numeric(Year2),
    Crop = factor(Crop, levels = c("Oat", "Soybean")),
    Period = factor(Period, levels = c("Jointing", "FlowerBudDifferentiation", "Filling", "Maturity")),
    Pattern = factor(Pattern, levels = c("monoculture", "rotation", "intercropping")),
    Treat_All = factor(Treat_All, levels = c("CK", "PP-P", "PP-F", "PLA-P", 'PLA-F')),
    Block = factor(Block, levels = c('1', '2', '3', '4')),
    Treat_Is = factor(Treat_Is, levels = c("Without", "With")),
    Treat_Type = factor(Treat_Type, levels = c("CK", "PP", "PLA")),
    Treat_Shape = factor(Treat_Shape, levels = c("CK", "Powder", "Fiber"))
  ) %>% rename_with(~ str_replace(., "^(\\d)", "X\\1"))
variables_original <-
  c("Yield", "Biomass", "HarvestIndex",
    "SpikeNumber", "KernelsPerSpike", "Weight1000Grain", "SpikeLength",
    "PodNumberPerPlant", "SeedNumberPerPlant", "SeedNumberPerPod", "Weight100Seed",
    "Height", "SeedWater",
    "BiomassAboveground", "BiomassBelowground", "RootShootRatio",
    "RootLength", "RootDiameter", "RootSurfaceArea", "RootVolume", "SpecificRootLength",
    "SWC", "BD", "SWCL40", "BDL40", "SWCL60", "BDL60",
    "pH", "EC", "SOC", "TN", "SOCTNratio",
    "NH4", "NO3", "AvailP", "DOC", "DON", "DOCNratio",
    "BG", "BX", "CBH", "LAP", "NAG", "ALP", "POX", "PER", "Cacq", "Nacq", "Pacq", "OX",
    "EnCratioCN", "EnCratioCP", "VectorLength", "VectorAngle", "CQI", "EnCNratio", "EnCPratio", "EnNPratio",
    "SQI", "EMF", "PlantProd", "PlantGrow", "PlantForBelow",
    "SoilNutriTurnover", "SoilActi", "SoilPhy", "SoilChem", "SoilCSta", "Support", "Regulation")
variables_log <- paste0(variables_original, "_log")
variables <- c(variables_original, variables_log)

dataTotal_202125MonoMaturity <- dataTotal_202125MonoMaturity %>% mutate(across(
  all_of(variables_original), ~ {
    if (all(is.na(.))) return(rep(NA, length(.)))
    min_val <- min(., na.rm = TRUE)
    if (min_val <= 0) {ifelse(is.finite(.), log(. + abs(min_val) + 1), NA)
    } else {ifelse(. > 0, log(.), NA)}}, .names = "{.col}_log"))

FactorAs <- c("Treat_All", "Treat_Is", "Treat_Type", "Treat_Shape")

formula_Randoms <- c("Site", "Year", "Crop")

start_time <- Sys.time()
for (factor in FactorAs) {
  DataForAnal <- dataTotal_202125MonoMaturity

  FactorA <- factor
  ideal_A_order <- switch(
    FactorA,
    "Treat_All" = c("CK", "PP-P", "PP-F", "PLA-P", "PLA-F"),
    "Treat_Is"  = c("Without", "With"),
    "Treat_Type" = c("CK", "PP", "PLA"),
    "Treat_Shape" = c("CK", "Powder", "Fiber"))
  existing_A <- ideal_A_order[ideal_A_order %in% unique(DataForAnal[[FactorA]])]
  DataForAnal[[FactorA]] <- factor(DataForAnal[[FactorA]], levels = existing_A)
  factorA_levels <- existing_A

  cat("\nAnalyzing factor:", FactorA, "\n")

  plan(multisession, workers = availableCores() - 1)
  results_list <- future_map(variables, function(var) {
    if (all(is.na(DataForAnal[[var]]))) return(NULL)
    analyze_variable_one(
      data = DataForAnal, response_var = var,
      formula_effects = FactorA, random_factors = formula_Randoms)},
    .progress = TRUE, .options = furrr_options(seed = TRUE)) %>% set_names(variables)
  results_list <- results_list %>% compact()
  plan(sequential)

  all_tables <- map(0:4, function(i) {
    table_key <- paste0("table", i)
    map_dfr(results_list, ~ .x[[table_key]])
  }) %>% set_names(paste0("table", 0:4))

  three_line_A_original <- create_three_line(all_tables$table4 %>% filter(
    Variable %in% variables_original), FactorA, factorA_levels)
  three_line_A_log <- create_three_line(all_tables$table4 %>% filter(
    Variable %in% variables_log), FactorA, factorA_levels)

  wb <- createWorkbook()
  sheet_names_original <-
    c("NorHomo", "ModelEval", "FixedANOVA", "RandomTest", "Tukey")
  walk2(seq_along(sheet_names_original), sheet_names_original, ~ {
    table_data <- all_tables[[.x]] %>%
      filter(Variable %in% variables_original) %>%
      arrange(match(Variable, variables_original))
    if (nrow(table_data) > 0) {
      addWorksheet(wb, .y)
      writeData(wb, .y, table_data)}})
  sheet_names_log <- c(
    "NorHomolog", "ModelEvallog", "FixedANOVAlog", "RandomTestlog", "Tukeylog")
  walk2(seq_along(sheet_names_log), sheet_names_log, ~ {
    table_data <- all_tables[[.x]] %>%
      filter(Variable %in% variables_log) %>%
      arrange(match(Variable, variables_log))
    if (nrow(table_data) > 0) {
      addWorksheet(wb, .y)
      writeData(wb, .y, table_data)}})
  addWorksheet(wb, "Line3")
  addWorksheet(wb, "Line3log")
  if (!is.null(three_line_A_original)) {writeData(wb, "Line3", three_line_A_original)}
  if (!is.null(three_line_A_log)) {writeData(wb, "Line3log", three_line_A_log)}

  output_filename <- paste0("Output/Sta_202125_Overview_", FactorA, "20260410.xlsx")
  saveWorkbook(wb, output_filename, overwrite = TRUE)
}
end_time <- Sys.time()
end_time - start_time

dataSta99Overview99Treat_Is <- read.xlsx(
  "Output/Sta_202125_Overview_Treat_Is20260410.xlsx",
  "FixedANOVA", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  filter(Term != "(Intercept)") %>% mutate(Sigused = SigF)

dataSta99Overview99Treat_All <- read.xlsx(
  "Output/Sta_202125_Overview_Treat_All20260410.xlsx",
  "Tukey", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  mutate(Sigused = SignoA)

variables_original <-
  c("Yield", "SQI", "EMF")

resp_df <- data.frame(
  NameyResp = variables_original,
  y_title_text =
    c("Yield ~ (t ~ ha^-1)", "SQI", "EMF"
    ), stringsAsFactors = FALSE)

# Fig. 1: overall treatment effects ----
NameLevel <- "Overview"
NamexTreat <- "Treat_Is"

for (i in 1:nrow(resp_df)) {
  NameyResp <- resp_df$NameyResp[i]
  y_title <- parse(text = resp_df$y_title_text[i])
  stat_data <- get(paste0("dataSta99", NameLevel, "99", NamexTreat)) %>% filter(Variable == NameyResp)
  if (nrow(stat_data) == 0) {message("Skipped: ", NameyResp, " - stat_data is empty"); next}
  sigtext <- stat_data %>% pull(Sigused) %>% ifelse(is.na(.), "", .)
  dataforplot <- dataTotal_202125MonoMaturity %>% mutate(
    xTreat = !!sym(NamexTreat), yResp = !!sym(NameyResp)) %>% filter(!is.na(yResp))
  OutputName1 <- paste0("Output/Raincloud-", NameLevel, "-", NameyResp, "-", NamexTreat, "20260624.pdf")
  p1 <- create_raincloud(
    dataforplot = dataforplot, y_title = y_title, sigtext = sigtext,
    point_size = 2.5)
  cairo_pdf(OutputName1, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54
      )
  print(p1); dev.off()
  OutputName2 <- paste0("Output/Forest-", NameLevel, "-", NameyResp, "-", NamexTreat, "20260624.pdf")
  p2 <- create_forest(dataforplot = dataforplot, y_title = y_title, sigtext = sigtext)
  cairo_pdf(OutputName2, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54
  )
  print(p2); dev.off()
}

factors_to_analyze <- c("Treat_All", "Treat_Is", "Treat_Type", "Treat_Shape")

start_time <- Sys.time()
for (Factor_Only in factors_to_analyze) {
  cat("\nAnalyzing factor:", Factor_Only, "\n")
  factor_levels <- switch(
    Factor_Only,
    "Treat_All" = c("CK", "PP-P", "PP-F", "PLA-P", "PLA-F"),
    "Treat_Is" = c("Without", "With"),
    "Treat_Type" = c("CK", "PP", "PLA"),
    "Treat_Shape" = c("CK", "Powder", "Fiber"))
  formula_effects <- Factor_Only
  filter_data <- function(data, siteX, yearX, cropX, patternX, periodX) {
    data %>% filter(
      Site == siteX, Year2 == yearX, Crop == cropX,
      Pattern == patternX, Period == periodX)}

  plan(multisession, workers = parallel::detectCores() - 1)
  valid_combinations <- expand_grid(
    var = variables,
    siteX = unique(dataTotal_202125MonoMaturity$Site),
    yearX = unique(dataTotal_202125MonoMaturity$Year2),
    cropX = unique(dataTotal_202125MonoMaturity$Crop),
    patternX = unique(dataTotal_202125MonoMaturity$Pattern),
    periodX = unique(dataTotal_202125MonoMaturity$Period)
  ) %>% mutate(has_data = future_pmap_lgl(
    list(siteX, yearX, cropX, patternX, periodX, var),
    function(siteX, yearX, cropX, patternX, periodX, var) {
      DataForAnal <- filter_data(dataTotal_202125MonoMaturity, siteX, yearX, cropX, patternX, periodX)
      !all(is.na(DataForAnal[[var]]))}, .progress = TRUE, .options = furrr_options(seed = TRUE))) %>%
    filter(has_data) %>% dplyr::select(-has_data)
  results_list <- future_pmap(valid_combinations, function(var, siteX, yearX, cropX, patternX, periodX) {
    DataForAnal <- filter_data(dataTotal_202125MonoMaturity, siteX, yearX, cropX, patternX, periodX)
    response_var <- var
    formula <- as.formula(paste(response_var, "~", formula_effects))
    analyze_variable_lm(
      data = DataForAnal, formula, response_var, siteX, yearX, cropX, patternX, periodX)
  }, .progress = TRUE, .options = furrr_options(seed = TRUE))
  names(results_list) <- with(
    valid_combinations, paste(siteX, yearX, cropX, patternX, periodX, var, sep = "_"))
  results_list <- compact(results_list)
  plan(sequential)

  all_tables <- map(0:3, function(i) {
    table_list <- map(results_list, ~ .x[[paste0("table", i)]])
    table_list <- compact(table_list)
    if (length(table_list) > 0) {bind_rows(table_list)} else NULL})
  names(all_tables) <- paste0("table", 0:3)

  arrange_data <- function(data, factor_var = NULL, var_levels = NULL) {
    data %>% arrange(
      factor(Site, levels = c("Zhangbei", "Ulanqab", "Youyu", "Chifeng")),
      Year2,
      factor(Crop, levels = c("Oat", "Soybean")),
      factor(Pattern, levels = c("monoculture", "rotation", "intercropping")),
      factor(Period, levels = c("Jointing", "FlowerBudDifferentiation", "Filling", "Maturity")),
      across(any_of(c(factor_var, "Variable")),
             ~ factor(., levels = c(factor_levels, var_levels) %||% levels(.))))}
  create_three_line_table <- function(tukey_table, factor_var) {
    if (is.null(tukey_table)) return(NULL)
    three_line_table <- tukey_table %>%
      mutate(Mean = sprintf("%.2f", as.numeric(Mean)),
             SE = sprintf("%.2f", as.numeric(SE)),
             value = if_else(is.na(SignoA) | SignoA == "", paste0(Mean, " ± ", SE),
                             paste0(Mean, " ± ", SE, " ", SignoA))) %>%
      dplyr::select(Site, Year2, Crop, Pattern, Period, all_of(factor_var), Variable, value) %>%
      pivot_wider(names_from = Variable, values_from = value) %>%
      arrange_data(factor_var = factor_var)
    return(three_line_table)}

  three_line_original <- create_three_line_table(
    all_tables$table3 %>% filter(Variable %in% variables_original), Factor_Only)
  three_line_log <- create_three_line_table(
    all_tables$table3 %>% filter(Variable %in% variables_log), Factor_Only)

  wb <- createWorkbook()
  sheet_names_original <- c("NorHomo", "ModelEval", "ANOVA", "Tukey")
  for (i in seq_along(sheet_names_original)) {
    addWorksheet(wb, sheet_names_original[i])
    table_data <- all_tables[[paste0("table", i - 1)]] %>%
      filter(Variable %in% variables_original)
    if (!is.null(table_data)) {
      table_data <- arrange_data(table_data, var_levels = variables_original)
      writeData(wb, sheet_names_original[i], table_data)}}

  sheet_names_log <- c("NorHomolog", "ModelEvallog", "ANOVAlog", "Tukeylog")
  for (i in seq_along(sheet_names_log)) {
    addWorksheet(wb, sheet_names_log[i])
    table_data <- all_tables[[paste0("table", i - 1)]] %>%
      filter(Variable %in% variables_log)
    if (!is.null(table_data)) {
      table_data <- arrange_data(table_data, var_levels = variables_log)
      writeData(wb, sheet_names_log[i], table_data)}}

  addWorksheet(wb, "Line3")
  addWorksheet(wb, "Line3log")
  if (!is.null(three_line_original)) {writeData(wb, "Line3", three_line_original)}
  if (!is.null(three_line_log)) {writeData(wb, "Line3log", three_line_log)}

  output_filename <- paste0("Output/Sta_202125_Specific_", Factor_Only, "20260410.xlsx")
  saveWorkbook(wb, output_filename, overwrite = TRUE)
}
end_time <- Sys.time()
end_time - start_time

FactorAs <- c("Treat_All", "Treat_Is", "Treat_Type", "Treat_Shape")
FactorB <- "Crop"

formula_Randoms <- c("Site", "Year")

start_time <- Sys.time()
for (factor in FactorAs) {
  DataForAnal <- dataTotal_202125MonoMaturity

  FactorA <- factor
  ideal_A_order <- switch(
    FactorA,
    "Treat_All" = c("CK", "PP-P", "PP-F", "PLA-P", "PLA-F"),
    "Treat_Is"  = c("Without", "With"),
    "Treat_Type" = c("CK", "PP", "PLA"),
    "Treat_Shape" = c("CK", "Powder", "Fiber"))
  existing_A <- ideal_A_order[ideal_A_order %in% unique(DataForAnal[[FactorA]])]
  DataForAnal[[FactorA]] <- factor(DataForAnal[[FactorA]], levels = existing_A)
  factorA_levels <- existing_A

  ideal_B_order <- levels(DataForAnal[[FactorB]])
  existing_B <- ideal_B_order[ideal_B_order %in% unique(DataForAnal[[FactorB]])]
  DataForAnal[[FactorB]] <- factor(DataForAnal[[FactorB]], levels = existing_B)
  factorB_levels <- existing_B

  workbook_list <- list()
  for (factorb in factorB_levels) {
    cat("\nAnalyzing factor:", FactorA, factorb, "\n")
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

    all_tables <- map(0:4, function(i) {
      table_key <- paste0("table", i)
      map_dfr(results_list, ~ .x[[table_key]])
    }) %>% set_names(paste0("table", 0:4))

    three_line_A_original <- create_three_line(all_tables$table4 %>% filter(
      Variable %in% variables_original), FactorA, factorA_levels)
    three_line_A_log <- create_three_line(all_tables$table4 %>% filter(
      Variable %in% variables_log), FactorA, factorA_levels)

    wb <- createWorkbook()
    sheet_names_original <- c(
      "NorHomo", "ModelEval", "FixedANOVA", "RandomTest", "Tukey")
    walk2(seq_along(sheet_names_original), sheet_names_original, ~ {
      table_data <- all_tables[[.x]] %>%
        filter(Variable %in% variables_original) %>%
        arrange(match(Variable, variables_original))
      if (nrow(table_data) > 0) {
        addWorksheet(wb, .y)
        writeData(wb, .y, table_data)}})
    sheet_names_log <- c(
      "NorHomolog", "ModelEvallog", "FixedANOVAlog", "RandomTestlog", "Tukeylog")
    walk2(seq_along(sheet_names_log), sheet_names_log, ~ {
      table_data <- all_tables[[.x]] %>%
        filter(Variable %in% variables_log) %>%
        arrange(match(Variable, variables_log))
      if (nrow(table_data) > 0) {
        addWorksheet(wb, .y)
        writeData(wb, .y, table_data)}})
    addWorksheet(wb, "Line3")
    addWorksheet(wb, "Line3log")
    if (!is.null(three_line_A_original)) {writeData(wb, "Line3", three_line_A_original)}
    if (!is.null(three_line_A_log)) {writeData(wb, "Line3log", three_line_A_log)}
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
    if (length(all_sheets) == 0) {
      cat("Sheet ", sheet_name, " has no data to merge\n")
      next}
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
  output_filename <- paste0("Output/Sta_202125_Sep_", FactorB, "_", FactorA, "20260410.xlsx")
  saveWorkbook(wb_combined, output_filename, overwrite = TRUE)
}
end_time <- Sys.time()
end_time - start_time

dataSta99Crop99Treat_Is <- read.xlsx(
  "Output/Sta_202125_Sep_Crop_Treat_Is20260410.xlsx",
  "FixedANOVA", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  filter(Term != "(Intercept)") %>% mutate(Sigused = SigF)

dataSta99Crop99Treat_All <- read.xlsx(
  "Output/Sta_202125_Sep_Crop_Treat_All20260410.xlsx",
  "Tukey", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  mutate(Sigused = SignoA)

variables_original <-
  c("Yield", "SQI", "EMF")

resp_df <- data.frame(
  NameyResp = variables_original,
  y_title_text =
    c("Yield ~ (t ~ ha^-1)", "SQI", "EMF"
    ), stringsAsFactors = FALSE)

# Fig. S2: crop-specific treatment effects ----
NameLevel <- "Crop"
NamexTreat <- "Treat_Is"

for (i in 1:nrow(resp_df)) {
  NameyResp <- resp_df$NameyResp[i]
  y_title <- parse(text = resp_df$y_title_text[i])
  stat_data <- get(paste0("dataSta99", NameLevel, "99", NamexTreat)) %>% filter(Variable == NameyResp)
  if (nrow(stat_data) == 0) {message("Skipped: ", NameyResp, " - stat_data is empty"); next}
  sigtext <- stat_data %>% pull(Sigused) %>% ifelse(is.na(.), "", .)
  dataforplot <- dataTotal_202125MonoMaturity %>% mutate(
    xTreat = !!sym(NamexTreat), yResp = !!sym(NameyResp)) %>% filter(!is.na(yResp)) %>%
    filter(!!sym(NameLevel) %in% (stat_data[[NameLevel]] %>% unique())) %>%
    mutate(!!sym(NameLevel) := factor(!!sym(NameLevel)))
  OutputName1 <- paste0("Output/Raincloud-", NameLevel, "-", NameyResp, "-", NamexTreat, "20260624.pdf")
  p1 <- create_raincloud_sep(
    dataforplot = dataforplot, y_title = y_title, sigtext = sigtext, dodge_width = 0.9,
    point_size = 2.5)
  cairo_pdf(OutputName1, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54
  )
  print(p1); dev.off()
  OutputName2 <- paste0("Output/Forest-", NameLevel, "-", NameyResp, "-", NamexTreat, "20260624.pdf")
  p2 <- create_forest_sep(
    dataforplot = dataforplot, y_title = y_title, sigtext = sigtext,
    dodge_width = 0.7)
  cairo_pdf(OutputName2, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54
  )
  print(p2); dev.off()
}

FactorAs <- c("Treat_All", "Treat_Is", "Treat_Type", "Treat_Shape")
FactorB <- "Site"

formula_Randoms <- c("Year", "Crop")

start_time <- Sys.time()
for (factor in FactorAs) {
  DataForAnal <- dataTotal_202125MonoMaturity

  FactorA <- factor
  ideal_A_order <- switch(
    FactorA,
    "Treat_All" = c("CK", "PP-P", "PP-F", "PLA-P", "PLA-F"),
    "Treat_Is"  = c("Without", "With"),
    "Treat_Type" = c("CK", "PP", "PLA"),
    "Treat_Shape" = c("CK", "Powder", "Fiber"))
  existing_A <- ideal_A_order[ideal_A_order %in% unique(DataForAnal[[FactorA]])]
  DataForAnal[[FactorA]] <- factor(DataForAnal[[FactorA]], levels = existing_A)
  factorA_levels <- existing_A

  ideal_B_order <- levels(DataForAnal[[FactorB]])
  existing_B <- ideal_B_order[ideal_B_order %in% unique(DataForAnal[[FactorB]])]
  DataForAnal[[FactorB]] <- factor(DataForAnal[[FactorB]], levels = existing_B)
  factorB_levels <- existing_B

  workbook_list <- list()
  for (factorb in factorB_levels) {
    cat("\nAnalyzing factor:", FactorA, factorb, "\n")
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

    all_tables <- map(0:4, function(i) {
      table_key <- paste0("table", i)
      map_dfr(results_list, ~ .x[[table_key]])
    }) %>% set_names(paste0("table", 0:4))

    three_line_A_original <- create_three_line(all_tables$table4 %>% filter(
      Variable %in% variables_original), FactorA, factorA_levels)
    three_line_A_log <- create_three_line(all_tables$table4 %>% filter(
      Variable %in% variables_log), FactorA, factorA_levels)

    wb <- createWorkbook()
    sheet_names_original <- c(
      "NorHomo", "ModelEval", "FixedANOVA", "RandomTest", "Tukey")
    walk2(seq_along(sheet_names_original), sheet_names_original, ~ {
      table_data <- all_tables[[.x]] %>%
        filter(Variable %in% variables_original) %>%
        arrange(match(Variable, variables_original))
      if (nrow(table_data) > 0) {
        addWorksheet(wb, .y)
        writeData(wb, .y, table_data)}})
    sheet_names_log <- c(
      "NorHomolog", "ModelEvallog", "FixedANOVAlog", "RandomTestlog", "Tukeylog")
    walk2(seq_along(sheet_names_log), sheet_names_log, ~ {
      table_data <- all_tables[[.x]] %>%
        filter(Variable %in% variables_log) %>%
        arrange(match(Variable, variables_log))
      if (nrow(table_data) > 0) {
        addWorksheet(wb, .y)
        writeData(wb, .y, table_data)}})
    addWorksheet(wb, "Line3")
    addWorksheet(wb, "Line3log")
    if (!is.null(three_line_A_original)) {writeData(wb, "Line3", three_line_A_original)}
    if (!is.null(three_line_A_log)) {writeData(wb, "Line3log", three_line_A_log)}
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
    if (length(all_sheets) == 0) {
      cat("Sheet ", sheet_name, " has no data to merge\n")
      next}
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
  output_filename <- paste0("Output/Sta_202125_Sep_", FactorB, "_", FactorA, "20260410.xlsx")
  saveWorkbook(wb_combined, output_filename, overwrite = TRUE)
}
end_time <- Sys.time()
end_time - start_time

dataSta99Site99Treat_Is <- read.xlsx(
  "Output/Sta_202125_Sep_Site_Treat_Is20260410.xlsx",
  "FixedANOVA", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  filter(Term != "(Intercept)") %>% mutate(Sigused = SigF)

dataSta99Site99Treat_All <- read.xlsx(
  "Output/Sta_202125_Sep_Site_Treat_All20260410.xlsx",
  "Tukey", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  mutate(Sigused = SignoA)

variables_original <-
  c("Yield", "SQI", "EMF")

resp_df <- data.frame(
  NameyResp = variables_original,
  y_title_text =
    c("Yield ~ (t ~ ha^-1)", "SQI", "EMF"
    ), stringsAsFactors = FALSE)

# Fig. S2: site-specific treatment effects ----
NameLevel <- "Site"
NamexTreat <- "Treat_Is"

for (i in 1:nrow(resp_df)) {
  NameyResp <- resp_df$NameyResp[i]
  y_title <- parse(text = resp_df$y_title_text[i])
  stat_data <- get(paste0("dataSta99", NameLevel, "99", NamexTreat)) %>% filter(Variable == NameyResp)
  if (nrow(stat_data) == 0) {message("Skipped: ", NameyResp, " - stat_data is empty"); next}
  sigtext <- stat_data %>% pull(Sigused) %>% ifelse(is.na(.), "", .)
  dataforplot <- dataTotal_202125MonoMaturity %>% mutate(
    xTreat = !!sym(NamexTreat), yResp = !!sym(NameyResp)) %>% filter(!is.na(yResp)) %>%
    filter(!!sym(NameLevel) %in% (stat_data[[NameLevel]] %>% unique())) %>%
    mutate(!!sym(NameLevel) := factor(!!sym(NameLevel)))
  OutputName1 <- paste0("Output/Raincloud-", NameLevel, "-", NameyResp, "-", NamexTreat, "20260624.pdf")
  p1 <- create_raincloud_sep(
    dataforplot = dataforplot, y_title = y_title, sigtext = sigtext, dodge_width = 0.9,
    point_size = 2.5)
  cairo_pdf(OutputName1, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54
  )
  print(p1); dev.off()
  OutputName2 <- paste0("Output/Forest-", NameLevel, "-", NameyResp, "-", NamexTreat, "20260624.pdf")
  p2 <- create_forest_sep(
    dataforplot = dataforplot, y_title = y_title, sigtext = sigtext,
    dodge_width = 0.7)
  cairo_pdf(OutputName2, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54
  )
  print(p2); dev.off()
}

FactorAs <- c("Treat_All", "Treat_Is", "Treat_Type", "Treat_Shape")
FactorB <- "Year2"

formula_Randoms <- c("Site", "Crop")

start_time <- Sys.time()
for (factor in FactorAs) {
  DataForAnal <- dataTotal_202125MonoMaturity

  DataForAnal[[FactorB]] <- factor(
    as.character(DataForAnal[[FactorB]]),
    levels = c("1", "2", "3", "4", "5"))

  FactorA <- factor
  ideal_A_order <- switch(
    FactorA,
    "Treat_All" = c("CK", "PP-P", "PP-F", "PLA-P", "PLA-F"),
    "Treat_Is"  = c("Without", "With"),
    "Treat_Type" = c("CK", "PP", "PLA"),
    "Treat_Shape" = c("CK", "Powder", "Fiber"))
  existing_A <- ideal_A_order[ideal_A_order %in% unique(DataForAnal[[FactorA]])]
  DataForAnal[[FactorA]] <- factor(DataForAnal[[FactorA]], levels = existing_A)
  factorA_levels <- existing_A

  ideal_B_order <- levels(DataForAnal[[FactorB]])
  existing_B <- ideal_B_order[ideal_B_order %in% unique(DataForAnal[[FactorB]])]
  DataForAnal[[FactorB]] <- factor(DataForAnal[[FactorB]], levels = existing_B)
  factorB_levels <- existing_B

  workbook_list <- list()
  for (factorb in factorB_levels) {
    cat("\nAnalyzing factor:", FactorA, factorb, "\n")
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

    all_tables <- map(0:4, function(i) {
      table_key <- paste0("table", i)
      map_dfr(results_list, ~ .x[[table_key]])
    }) %>% set_names(paste0("table", 0:4))

    three_line_A_original <- create_three_line(all_tables$table4 %>% filter(
      Variable %in% variables_original), FactorA, factorA_levels)
    three_line_A_log <- create_three_line(all_tables$table4 %>% filter(
      Variable %in% variables_log), FactorA, factorA_levels)

    wb <- createWorkbook()
    sheet_names_original <- c(
      "NorHomo", "ModelEval", "FixedANOVA", "RandomTest", "Tukey")
    walk2(seq_along(sheet_names_original), sheet_names_original, ~ {
      table_data <- all_tables[[.x]] %>%
        filter(Variable %in% variables_original) %>%
        arrange(match(Variable, variables_original))
      if (nrow(table_data) > 0) {
        addWorksheet(wb, .y)
        writeData(wb, .y, table_data)}})
    sheet_names_log <- c(
      "NorHomolog", "ModelEvallog", "FixedANOVAlog", "RandomTestlog", "Tukeylog")
    walk2(seq_along(sheet_names_log), sheet_names_log, ~ {
      table_data <- all_tables[[.x]] %>%
        filter(Variable %in% variables_log) %>%
        arrange(match(Variable, variables_log))
      if (nrow(table_data) > 0) {
        addWorksheet(wb, .y)
        writeData(wb, .y, table_data)}})
    addWorksheet(wb, "Line3")
    addWorksheet(wb, "Line3log")
    if (!is.null(three_line_A_original)) {writeData(wb, "Line3", three_line_A_original)}
    if (!is.null(three_line_A_log)) {writeData(wb, "Line3log", three_line_A_log)}
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
    if (length(all_sheets) == 0) {
      cat("Sheet ", sheet_name, " has no data to merge\n")
      next}
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
  output_filename <- paste0("Output/Sta_202125_Sep_", FactorB, "_", FactorA, "20260410.xlsx")
  saveWorkbook(wb_combined, output_filename, overwrite = TRUE)
}
end_time <- Sys.time()
end_time - start_time

dataSta99Year299Treat_Is <- read.xlsx(
  "Output/Sta_202125_Sep_Year2_Treat_Is20260410.xlsx",
  "FixedANOVA", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  filter(Term != "(Intercept)") %>% mutate(Sigused = SigF)

dataSta99Year299Treat_All <- read.xlsx(
  "Output/Sta_202125_Sep_Year2_Treat_All20260410.xlsx",
  "Tukey", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  mutate(Sigused = SignoA)

variables_original <-
  c("Yield", "SQI", "EMF")

resp_df <- data.frame(
  NameyResp = variables_original,
  y_title_text =
    c("Yield ~ (t ~ ha^-1)", "SQI", "EMF"
    ), stringsAsFactors = FALSE)

# Fig. S3: year-specific treatment effects ----
NameLevel <- "Year2"
NamexTreat <- "Treat_Is"

for (i in 1:nrow(resp_df)) {
  NameyResp <- resp_df$NameyResp[i]
  y_title <- parse(text = resp_df$y_title_text[i])
  stat_data <- get(paste0("dataSta99", NameLevel, "99", NamexTreat)) %>% filter(Variable == NameyResp)
  if (nrow(stat_data) == 0) {message("Skipped: ", NameyResp, " - stat_data is empty"); next}
  sigtext <- stat_data %>% pull(Sigused) %>% ifelse(is.na(.), "", .)
  dataforplot <- dataTotal_202125MonoMaturity %>% mutate(
    xTreat = !!sym(NamexTreat), yResp = !!sym(NameyResp)) %>% filter(!is.na(yResp)) %>%
    filter(!!sym(NameLevel) %in% (stat_data[[NameLevel]] %>% unique())) %>%
    mutate(!!sym(NameLevel) := factor(!!sym(NameLevel)))
  OutputName1 <- paste0("Output/Raincloud-", NameLevel, "-", NameyResp, "-", NamexTreat, "20260624.pdf")
  p1 <- create_raincloud_sep(
    dataforplot = dataforplot, y_title = y_title, sigtext = sigtext, dodge_width = 0.9,
    point_size = 2.5)
  cairo_pdf(OutputName1, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54
  )
  print(p1); dev.off()
  OutputName2 <- paste0("Output/Forest-", NameLevel, "-", NameyResp, "-", NamexTreat, "20260624.pdf")
  p2 <- create_forest_sep(
    dataforplot = dataforplot, y_title = y_title, sigtext = sigtext,
    dodge_width = 0.7)
  cairo_pdf(OutputName2, bg = "transparent", width = 34.32/2.54, height = 4.9/2.54
  )
  print(p2); dev.off()
}

FactorAs <- c("Treat_All", "Treat_Is", "Treat_Type", "Treat_Shape")
FactorB <- "Year2"

formula_Randoms <- c("Site", "Crop")

start_time <- Sys.time()
for (factor in FactorAs) {
  DataForAnal <- dataTotal_202125MonoMaturity

  FactorA <- factor
  ideal_A_order <- switch(
    FactorA, "Treat_All" = c("CK", "PP-P", "PP-F", "PLA-P", "PLA-F"),
    "Treat_Is"  = c("Without", "With"), "Treat_Type" = c("CK", "PP", "PLA"),
    "Treat_Shape" = c("CK", "Powder", "Fiber"))
  existing_A <- ideal_A_order[ideal_A_order %in% unique(DataForAnal[[FactorA]])]
  DataForAnal[[FactorA]] <- factor(DataForAnal[[FactorA]], levels = existing_A)
  DataForAnal[[FactorB]] <- as.numeric(DataForAnal[[FactorB]])

  formula_Fix <- paste(FactorA, "*", FactorB)

  cat("\nAnalyzing factor:", FactorA, "*", FactorB, "\n")

  plan(multisession, workers = availableCores() - 1)
  results_list <- future_map(variables, function(var) {
    if (all(is.na(DataForAnal[[var]]))) return(NULL)
    analyze_variable_tworeg(
      data = DataForAnal, response_var = var,
      formula_effects = formula_Fix, random_factors = formula_Randoms)},
    .progress = TRUE, .options = furrr_options(seed = TRUE)) %>% set_names(variables)
  results_list <- results_list %>% compact()
  plan(sequential)

  all_tables <- map(0:5, function(i) {
    table_key <- paste0("table", i)
    map_dfr(results_list, ~ .x[[table_key]])
  }) %>% set_names(paste0("table", 0:5))

  wb <- createWorkbook()
  sheet_names_original <- c(
    "NorHomo", "ModelEval", "FixedANOVA", "RandomTest", "RegTotal", "RegIndi")
  walk2(seq_along(sheet_names_original), sheet_names_original, ~ {
    table_data <- all_tables[[.x]] %>%
      filter(Variable %in% variables_original) %>%
      arrange(match(Variable, variables_original))
    if (nrow(table_data) > 0) {
      addWorksheet(wb, .y)
      writeData(wb, .y, table_data)}})
  sheet_names_log <- c(
    "NorHomolog", "ModelEvallog", "FixedANOVAlog", "RandomTestlog", "RegTotallog", "RegIndilog")
  walk2(seq_along(sheet_names_log), sheet_names_log, ~ {
    table_data <- all_tables[[.x]] %>%
      filter(Variable %in% variables_log) %>%
      arrange(match(Variable, variables_log))
    if (nrow(table_data) > 0) {
      addWorksheet(wb, .y)
      writeData(wb, .y, table_data)}})

  output_filename <- paste0("Output/Sta_202125_Two_", FactorB, "_", FactorA, "20260410.xlsx")
  saveWorkbook(wb, output_filename, overwrite = TRUE)
}
end_time <- Sys.time()
end_time - start_time

dataSta99Year299Treat_Is <- read.xlsx(
  "Output/Sta_202125_Two_Year2_Treat_Is20260410.xlsx",
  "RegIndi", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  mutate(Sigused1 = SloPSig, Sigused2 = SloSignoA)

dataSta99Year299Treat_All <- read.xlsx(
  "Output/Sta_202125_Two_Year2_Treat_All20260410.xlsx",
  "RegIndi", na.strings = "") %>% type.convert(as.is = TRUE) %>%
  mutate(Sigused1 = SloPSig, Sigused2 = SloSignoA)

variables_original <-
  c("Yield", "SQI", "EMF")

resp_df <- data.frame(
  NameyResp = variables_original,
  y_title_text =
    c("Yield ~ (t ~ ha^-1)", "SQI", "EMF"
    ), stringsAsFactors = FALSE)

NameLevel <- "Year2"
x_title <- "Duration (a)"

NamexTreat <- "Treat_Is"

for (i in 1:nrow(resp_df)) {
  NameyResp <- resp_df$NameyResp[i]
  y_title <- parse(text = resp_df$y_title_text[i])
  stat_data <- get(paste0("dataSta99", NameLevel, "99", NamexTreat)) %>% filter(Variable == NameyResp) %>% mutate(
    xTreat = !!sym(NamexTreat))
  if (nrow(stat_data) == 0) {message("Skipped: ", NameyResp, " - stat_data is empty"); next}
  dataforplot <- dataTotal_202125MonoMaturity %>% mutate(
    xTreat = !!sym(NamexTreat), yResp = !!sym(NameyResp)) %>% filter(!is.na(yResp))
  OutputNameR <- paste0("Output/Regression-", NameLevel, "-", NameyResp, "-", NamexTreat, "20260624.pdf")
  pR <- create_regression(
    dataforplot = dataforplot, stat_data = stat_data, y_title = y_title, x_title = x_title)
  cairo_pdf(OutputNameR, bg = "transparent", width = 34.32/2.54, height = 5.44/2.54)
  print(pR); dev.off()
}

# Fig. S1a: field-site map ----
pacman::p_load(openxlsx, tidyverse, sf, terra, tidyterra)

BaseMap <- rast("D:/ArcGIS/GlobalAIPET/ai_v31_yr.tif") / 10000
points_df <- read.xlsx(field_file, "PointLL20260327")
pts <- vect(points_df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")

xrange <- diff(range(points_df$Longitude))
yrange <- diff(range(points_df$Latitude))
plot_ext <- ext(
  min(points_df$Longitude) - xrange * 0.2,
  max(points_df$Longitude) + xrange * 0.5,
  min(points_df$Latitude)  - yrange * 1.4,
  max(points_df$Latitude)  + yrange * 0.6)

BaseMap_sub <- crop(BaseMap, plot_ext)
BaseMap_sub <- classify(BaseMap_sub, rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

ggplot() +
  geom_spatraster(data = BaseMap_sub) +
  scale_fill_gradientn(
    name = NULL, limits  = c(0, 1), na.value = NA, guide = "colourbar",
    colours = c("#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#80cdc1", "#35978f"),
    values  = scales::rescale(c(0, 0.03, 0.2, 0.5, 0.65, 1)),
    breaks  = c(0.03, 0.2, 0.5, 0.65, 1)) +
  geom_spatvector(data = pts, color = "red", size = 2.5) +
  geom_text(data = points_df, aes(x = Longitude, y = Latitude, label = Site),
            nudge_y = 0.3, size = 16, size.unit = "pt") +
  coord_sf(xlim = c(xmin(plot_ext), xmax(plot_ext)), ylim = c(ymin(plot_ext), ymax(plot_ext)), expand = FALSE) +
  scale_x_continuous("Longitude(°)", breaks = c(115, 120), labels = c("115", "120"), expand = expansion(mult = 0)) +
  scale_y_continuous("Latitude(°)", breaks = c(38, 42), labels = c("38", "42"), expand = expansion(mult = 0)) +
  theme_classic() + theme(
    axis.line = element_blank(),
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 16),
    plot.margin = margin(0, 0, 0, 0),
    legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
    legend.box = "vertical",
    legend.background = element_rect(fill = scales::alpha("grey90", 0.4), colour = NA),
    legend.key = element_blank(), legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)) -> p_Map0; p_Map0
pdf("Output/MapLocation20260328.pdf", width = 17/2.54)
p_Map0
dev.off()

world <- st_read("D:/ArcGIS/GlobalCountry/global_all_country.shp", quiet = TRUE)
tenline <- st_read("D:/ArcGIS/CTAmap/十段线.shp", quiet = TRUE)
world <- st_transform(world, 4326)
tenline <- st_transform(tenline, 4326)
cn <- world %>% filter(CNTRY_NAME == "China")

plot_box <- as.data.frame(as.list(plot_ext))

ggplot() +
  geom_sf(data = world, fill = scales::alpha("#F5F5F5", 0.2), color = "#C7C7C7", linewidth = 0.1) +
  geom_sf(data = cn, fill = scales::alpha("#F5F5F5", 0.2), color = "#7A7A7A", linewidth = 0.1) +
  geom_sf(data = tenline, color = "#7A7A7A", linewidth = 0.1) +
  geom_rect(data = plot_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "#8c510a", linewidth = 0.6, linetype = "solid") +
  coord_sf(xlim = c(72, 138), ylim = c(0, 56),
           expand = FALSE, crs = st_crs(4326), datum = st_crs(4326)) +
  scale_x_continuous(breaks = c(80, 120), labels = function(x) paste0(x)) +
  scale_y_continuous(breaks = c(10, 40), labels = function(y) paste0(y)) +
  theme_classic() + theme(
    axis.line = element_blank(),
    panel.background = element_rect(fill = scales::alpha("grey90", 0.2), color = NA),
    plot.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(color = "#7A7A7A", size = 9),
    axis.ticks = element_line(color = "#7A7A7A", linewidth = 0.25),
    panel.border = element_rect(color = "#7A7A7A", fill = NA, linewidth = 0.25),
    plot.margin = margin(0, 0, 0, 0)
  ) -> p_Map00; p_Map00
pdf("Output/MapChina20260328.pdf", width = 5.24/2.54)
p_Map00
dev.off()


# Fig. S1b: monthly climate panel ----
monthly_climate <- read.xlsx(field_file, "monthly_climate") %>% as_tibble()
if (nrow(monthly_climate) == 0) {
  warning("Fig. S1b was not generated: FieldData.xlsx/monthly_climate has no rows.")
} else {
  monthly_climate <- monthly_climate %>% select(any_of(c("Station", "Site", "Year", "Month", "Temperature", "Precipitation", "GrowingSeason"))) %>% mutate(
    Site = factor(Site, levels = c("Zhangbei", "Ulanqab", "Youyu", "Chifeng")),
    Year = as.numeric(Year), Month = as.numeric(Month),
    Temperature = as.numeric(Temperature), Precipitation = as.numeric(Precipitation))
  p_Climate <- ggplot(monthly_climate, aes(x = Month)) +
    geom_col(aes(y = Precipitation), fill = "#95B7DA", width = 0.75) +
    geom_line(aes(y = Temperature * 10), color = "#D6AB85", linewidth = 0.45) +
    geom_rect(data = distinct(monthly_climate, Site, Year), aes(xmin = 5, xmax = 9, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey80", alpha = 0.25) +
    facet_grid(Site ~ Year) +
    scale_x_continuous(breaks = c(1, 6, 12)) +
    scale_y_continuous("Precipitation", sec.axis = sec_axis(~ . / 10, name = "Temperature")) +
    theme_classic() + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 10), strip.background = element_blank())
  pdf("Output/Climate20260410.pdf", height = 10.9/2.54)
  print(p_Climate)
  dev.off()
}
