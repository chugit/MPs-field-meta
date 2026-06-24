# Setup ----
rm(list = ls())
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(openxlsx, tidyverse, metafor, metaforest, glmulti, MuMIn, future.apply, patchwork, terra, tidyterra, plotbiomes)
input_dir <- if (dir.exists("Input")) "Input" else "."
meta_file <- file.path(input_dir, "metaData.xlsx")
if (!file.exists(meta_file)) stop("metaData.xlsx not found. Put it in the working directory or Input/.")
dir.create("Output", showWarnings = FALSE, recursive = TRUE)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

# Custom functions ----
cv_avg <- function(x, sd, n, group, data, label = NULL, sub_b = TRUE, cv2 = FALSE) {
  if (is.null(label)) {label <- purrr::map_chr(enquos(x), rlang::as_label)}

  b_grp_cv_data <- data %>%
    dplyr::group_by({{group}}) %>%
    dplyr::mutate(w_CV2 = weighted_CV({{sd}}, {{x}}, {{n}}, cv2 = cv2),
                  n_mean = mean({{n}}, na.rm = TRUE)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(b_CV2 = weighted.mean(w_CV2, n_mean, na.rm = TRUE), .keep = "used")

  names(b_grp_cv_data) <- paste0(names(b_grp_cv_data), "_", label)

  if (sub_b) {
    b_grp_cv_data <- b_grp_cv_data %>% dplyr::select(grep("b_", names(b_grp_cv_data)))
    dat_new <- cbind(data, b_grp_cv_data)} else {dat_new <- cbind(data, b_grp_cv_data)}

  return(data.frame(dat_new))
}

weighted_CV <- function(sd, x, n, cv2 = FALSE) {
  if (cv2) {weighted.mean(na_if((sd / abs(x))^2, Inf), n, na.rm = TRUE)
  } else {weighted.mean(na_if((sd / abs(x)), Inf), n, na.rm = TRUE)^2}
}

get_est <- function(model) {
  est <- coef(model)
  ci.lb <- model$ci.lb
  ci.ub <- model$ci.ub
  se <- model$se
  return(data.frame(Est. = est, SE = se, "95% LCI" = ci.lb, "95% UCI" = ci.ub, check.names = FALSE))
}

lnrr_laj <- function(m1, m2, cv1_2, cv2_2, n1, n2, taylor = TRUE) {
  if (taylor) {log(m2 / m1) + 0.5*((cv2_2 / n2) - (cv1_2 / n1))
  } else {log(m2 / m1)}
}

v_lnrr_laj <- function(cv1_2, cv2_2, n1, n2, taylor = TRUE) {
  if (taylor) {((cv1_2) / n1) + ((cv2_2) / n2) +
      ((cv1_2)^2 / (2*n1^2)) + ((cv2_2)^2 / (2*n2^2))
  } else {((cv1_2) / n1) + ((cv2_2) / n2)}
}

calc.v <- function(data, m1, sd1, n1, vi) {
  v <- matrix((data[[sd1]][1]^2 / (data[[n1]][1] * data[[m1]][1]^2)),
              nrow = nrow(data), ncol = nrow(data))
  diag(v) <- data[[vi]]
  v
}
calc.vbd <- function(data, CommonID, m1, sd1, n1, vi) {
  vars <- sapply(match.call()[-1], deparse)[2:6]
  CommonID <- vars[1]; m1 <- vars[2]; sd1 <- vars[3]; n1 <- vars[4]; vi <- vars[5]
  full_V <- matrix(0, nrow(data), nrow(data))
  for (g in unique(data[[CommonID]])) {
    idx <- which(data[[CommonID]] == g)
    sub <- data[idx, ]
    cov_val <- sub[[sd1]][1]^2 / (sub[[n1]][1] * sub[[m1]][1]^2)
    V_g <- matrix(cov_val, length(idx), length(idx))
    diag(V_g) <- sub[[vi]]
    full_V[idx, idx] <- V_g}
  full_V
}

# Load rebuilt meta-analysis base table ----
data_combined_wide3 <- read.xlsx(meta_file, "data_combined_wide") %>%
  group_by(VariableNew, SameVarID) %>% mutate(SameVarIDN = cur_group_id()) %>% ungroup()

# Fig. 2, Fig. S5-S9 and Table S6 meta-analysis ----
# Effect-size calculation ----
data_combined_with_con <- data_combined_wide3 %>% filter(!is.na(ConMean))

data_combined_with_con2 <- data_combined_with_con %>%
  group_by(VariableNew) %>% group_modify(~ {
    cv_avg(x = ConMean, sd = ConSD, n = ConN,
           group = ExpID,
           label = "1",
           data = .x)}) %>% ungroup() %>%
  group_by(VariableNew) %>% group_modify(~ {
    cv_avg(x = TreatMean, sd = TreatSD, n = TreatN,
           group = ExpID, label = "2",
           data = .x)}) %>% ungroup() %>%
  mutate(across(c(b_CV2_1, b_CV2_2), ~ ifelse(is.nan(.), NA, .)))

data_combined_with_con3 <- data_combined_with_con2 %>% mutate(
  ConSD_new = if_else(is.na(ConSD), abs(ConMean) * sqrt(b_CV2_1), ConSD),
  TreatSD_new = if_else(is.na(TreatSD), abs(TreatMean) * sqrt(b_CV2_2), TreatSD))

data_combined_with_con4 <- data_combined_with_con3 %>% mutate(
  ConCV2_new = na_if(ConSD_new / abs(ConMean), Inf)^2,
  TreatCV2_new = na_if(TreatSD_new / abs(TreatMean), Inf)^2)

data_combined_with_con5 <- data_combined_with_con4 %>% group_by(CommonIDN) %>%
  filter(all(sign(TreatMean) == sign(ConMean))) %>% ungroup() %>% mutate(
    yi = lnrr_laj(m1 = ConMean, m2 = TreatMean, cv1_2 = ConCV2_new,
                  cv2_2 = TreatCV2_new, n1 = ConN, n2 = TreatN),
    vi = v_lnrr_laj(cv1_2 = ConCV2_new, n1 = ConN,
                    cv2_2 = TreatCV2_new, n2 = TreatN))

data_combined_with_con6 <- data_combined_with_con5 %>% filter(!is.na(EMFCat)) %>%
  filter(!(is.na(ConSD_new) & is.na(TreatSD_new))) %>%
  mutate(ess.var = 1/ConN + 1/TreatN, ess.se = sqrt(ess.var))

col_names <- c("ExpID", "SampleTime2", "SoilDepth", "Niche", "SameVarIDN")
sapply(data_combined_with_con6[col_names], function(x) length(unique(x)))

# Fig. S5: plant and soil variable-level effects ----
data_plant <- data_combined_with_con6 %>% filter(Niche == "plant")

Yield_data <- data_plant %>% filter(
  VariableNew == "Yield", Irrigation %in% c("No"))

sapply(Yield_data[col_names], function(x) length(unique(x)))
V_Yield <- calc.vbd(data = Yield_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_Yield <- rma.mv(
  yi, V_Yield, data = Yield_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  test = "t", dfs = "contain",
  control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_Yield)

model_Yield_ess.se <- rma.mv(
  yi, V_Yield, data = Yield_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_Yield_ess.se)

BiomAbove_data <- data_plant %>% filter(
  VariableNew %in% c("BiomAbove", "BiomStraw"), Irrigation %in% c("No"))
sapply(BiomAbove_data[col_names], function(x) length(unique(x)))
V_BiomAbove <- calc.vbd(data = BiomAbove_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_BiomAbove <- rma.mv(
  yi, V_BiomAbove, data = BiomAbove_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_BiomAbove)
model_BiomAbove_ess.se <- rma.mv(
  yi, V_BiomAbove, data = BiomAbove_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_BiomAbove_ess.se)

PlantHeightAbove_data <- data_plant %>% filter(
  VariableNew == "PlantHeightAbove", Irrigation %in% c("No"))
sapply(PlantHeightAbove_data[col_names], function(x) length(unique(x)))
V_PlantHeightAbove <- calc.vbd(data = PlantHeightAbove_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_PlantHeightAbove <- rma.mv(
  yi, V_PlantHeightAbove, data = PlantHeightAbove_data,
  random = list(~ 1 | ExpID),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_PlantHeightAbove)
model_PlantHeightAbove_ess.se <- rma.mv(
  yi, V_PlantHeightAbove, data = PlantHeightAbove_data,
  random = list(~ 1 | ExpID),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_PlantHeightAbove_ess.se)

Chlorophyll_data <- data_plant %>% filter(
  VariableNew == "Chlorophyll", Irrigation %in% c("No"))
sapply(Chlorophyll_data[col_names], function(x) length(unique(x)))
V_Chlorophyll <- calc.vbd(data = Chlorophyll_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_Chlorophyll <- rma.mv(
  yi, V_Chlorophyll, data = Chlorophyll_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_Chlorophyll)
model_Chlorophyll_ess.se <- rma.mv(
  yi, V_Chlorophyll, data = Chlorophyll_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_Chlorophyll_ess.se)

RootLength_data <- data_plant %>% filter(
  VariableNew == "RootLength", Irrigation %in% c("No"))
sapply(RootLength_data[col_names], function(x) length(unique(x)))
V_RootLength <- calc.vbd(data = RootLength_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_RootLength <- rma.mv(
  yi, V_RootLength, data = RootLength_data,
  random = list(~ 1 | ExpID),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_RootLength)
model_RootLength_ess.se <- rma.mv(
  yi, V_RootLength, data = RootLength_data,
  random = list(~ 1 | ExpID),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_RootLength_ess.se)

model_list_Plant <- list(
  Yield = model_Yield,
  Biomass = model_BiomAbove,
  Height = model_PlantHeightAbove,
  Chlorophyll = model_Chlorophyll,
  RootLength = model_RootLength)
results_Plant <- data.frame(model = factor(
  names(model_list_Plant), levels = c("Chlorophyll", "RootLength", "Height", "Biomass", "Yield")),
  estimate = sapply(model_list_Plant, coef),
  ci.lb = sapply(model_list_Plant, function(m) m$ci.lb),
  ci.ub = sapply(model_list_Plant, function(m) m$ci.ub),
  n = sapply(model_list_Plant, function(m) m$k)) %>%
  mutate(contains_zero = ci.lb <= 0 & ci.ub >= 0)

data_soil <- data_combined_with_con6 %>% filter(Niche != "plant" | is.na(PlantAppear))

BD_data <- data_soil %>% filter(
  VariableNew == "BD", Irrigation %in% c("No"))
sapply(BD_data[col_names], function(x) length(unique(x)))
V_BD <- calc.vbd(data = BD_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_BD <- rma.mv(
  yi, V_BD, data = BD_data,
  random = list(~ 1 | ExpID),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_BD)
model_BD_ess.se <- rma.mv(
  yi, V_BD, data = BD_data,
  random = list(~ 1 | ExpID),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_BD_ess.se)

pH_data <- data_soil %>% filter(
  VariableNew == "pH", Irrigation %in% c("No"))
sapply(pH_data[col_names], function(x) length(unique(x)))
V_pH <- calc.vbd(data = pH_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_pH <- rma.mv(
  yi, V_pH, data = pH_data,
  random = list(~ 1 | ExpID),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_pH)
model_pH_ess.se <- rma.mv(
  yi, V_pH, data = pH_data,
  random = list(~ 1 | ExpID),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_pH_ess.se)

EC_data <- data_soil %>% filter(
  VariableNew == "EC", Irrigation %in% c("No"))
sapply(EC_data[col_names], function(x) length(unique(x)))
V_EC <- calc.vbd(data = EC_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_EC <- rma.mv(
  yi, V_EC, data = EC_data,
  random = list(~ 1 | ExpID),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_EC)
model_EC_ess.se <- rma.mv(
  yi, V_EC, data = EC_data,
  random = list(~ 1 | ExpID),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_EC_ess.se)

SOC_data <- data_soil %>% filter(
  VariableNew == "SOC", Irrigation %in% c("No"))
sapply(SOC_data[col_names], function(x) length(unique(x)))
V_SOC <- calc.vbd(data = SOC_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_SOC <- rma.mv(
  yi, V_SOC, data = SOC_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_SOC)
model_SOC_ess.se <- rma.mv(
  yi, V_SOC, data = SOC_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_SOC_ess.se)

TN_data <- data_soil %>% filter(
  VariableNew == "TN", Irrigation %in% c("No"))
sapply(TN_data[col_names], function(x) length(unique(x)))
V_TN <- calc.vbd(data = TN_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_TN <- rma.mv(
  yi, V_TN, data = TN_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_TN)
model_TN_ess.se <- rma.mv(
  yi, V_TN, data = TN_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_TN_ess.se)

DOC_data <- data_soil %>% filter(
  VariableNew == "DOC", Irrigation %in% c("No"))
sapply(DOC_data[col_names], function(x) length(unique(x)))
V_DOC <- calc.vbd(data = DOC_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_DOC <- rma.mv(
  yi, V_DOC, data = DOC_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_DOC)
model_DOC_ess.se <- rma.mv(
  yi, V_DOC, data = DOC_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_DOC_ess.se)

AvailN_data <- data_soil %>% filter(
  VariableNew == "AvailN", Irrigation %in% c("No"))
sapply(AvailN_data[col_names], function(x) length(unique(x)))
V_AvailN <- calc.vbd(data = AvailN_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_AvailN <- rma.mv(
  yi, V_AvailN, data = AvailN_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_AvailN)
model_AvailN_ess.se <- rma.mv(
  yi, V_AvailN, data = AvailN_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_AvailN_ess.se)

AvailP_data <- data_soil %>% filter(
  VariableNew == "AvailP", Irrigation %in% c("No"))
sapply(AvailP_data[col_names], function(x) length(unique(x)))
V_AvailP <- calc.vbd(data = AvailP_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_AvailP <- rma.mv(
  yi, V_AvailP, data = AvailP_data,
  random = list(~ 1 | ExpID),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_AvailP)
model_AvailP_ess.se <- rma.mv(
  yi, V_AvailP, data = AvailP_data,
  random = list(~ 1 | ExpID),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_AvailP_ess.se)

CAcqActi_data <- data_soil %>% filter(
  VariableNew == "CAcqActi", Irrigation %in% c("No"))
sapply(CAcqActi_data[col_names], function(x) length(unique(x)))
V_CAcqActi <- calc.vbd(data = CAcqActi_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_CAcqActi <- rma.mv(
  yi, V_CAcqActi, data = CAcqActi_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_CAcqActi)
model_CAcqActi_ess.se <- rma.mv(
  yi, V_CAcqActi, data = CAcqActi_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_CAcqActi_ess.se)

NAcqActi_data <- data_soil %>% filter(
  VariableNew == "NAcqActi", Irrigation %in% c("No"))
sapply(NAcqActi_data[col_names], function(x) length(unique(x)))
V_NAcqActi <- calc.vbd(data = NAcqActi_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_NAcqActi <- rma.mv(
  yi, V_NAcqActi, data = NAcqActi_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_NAcqActi)
model_NAcqActi_ess.se <- rma.mv(
  yi, V_NAcqActi, data = NAcqActi_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_NAcqActi_ess.se)

PAcqActi_data <- data_soil %>% filter(
  VariableNew == "PAcqActi", Irrigation %in% c("No"))
sapply(PAcqActi_data[col_names], function(x) length(unique(x)))
V_PAcqActi <- calc.vbd(data = PAcqActi_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_PAcqActi <- rma.mv(
  yi, V_PAcqActi, data = PAcqActi_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_PAcqActi)
model_PAcqActi_ess.se <- rma.mv(
  yi, V_PAcqActi, data = PAcqActi_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_PAcqActi_ess.se)

EnzyActi_data <- data_soil %>% filter(
  VariableNew %in% c("CAcqActi", "NAcqActi", "PAcqActi"), Irrigation %in% c("No"))
sapply(EnzyActi_data[col_names], function(x) length(unique(x)))
V_EnzyActi <- calc.vbd(data = EnzyActi_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_EnzyActi <- rma.mv(
  yi, V_EnzyActi, data = EnzyActi_data,
  random = list(~ 1 | ExpID, ~ 1 | VariableNew/SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_EnzyActi)
model_EnzyActi_ess.se <- rma.mv(
  yi, V_EnzyActi, data = EnzyActi_data,
  random = list(~ 1 | ExpID, ~ 1 | VariableNew/SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_EnzyActi_ess.se)

OXActi_data <- data_soil %>% filter(
  VariableNew %in% c("OXActi"), Irrigation %in% c("No"))
sapply(OXActi_data[col_names], function(x) length(unique(x)))
V_OXActi <- calc.vbd(data = OXActi_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_OXActi <- rma.mv(
  yi, V_OXActi, data = OXActi_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_OXActi)
model_OXActi_ess.se <- rma.mv(
  yi, V_OXActi, data = OXActi_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_OXActi_ess.se)

MicroBiomass_data <- data_soil %>% filter(
  VariableNew == "MicroBiomass", Irrigation %in% c("No"))
sapply(MicroBiomass_data[col_names], function(x) length(unique(x)))
V_MicroBiomass <- calc.vbd(data = MicroBiomass_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_MicroBiomass <- rma.mv(
  yi, V_MicroBiomass, data = MicroBiomass_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_MicroBiomass)
model_MicroBiomass_ess.se <- rma.mv(
  yi, V_MicroBiomass, data = MicroBiomass_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_MicroBiomass_ess.se)

MicroRichness_data <- data_soil %>% filter(
  VariableNew == "MicroRichness", Irrigation %in% c("No"))
sapply(MicroRichness_data[col_names], function(x) length(unique(x)))
V_MicroRichness <- calc.vbd(data = MicroRichness_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_MicroRichness <- rma.mv(
  yi, V_MicroRichness, data = MicroRichness_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_MicroRichness)
model_MicroRichness_ess.se <- rma.mv(
  yi, V_MicroRichness, data = MicroRichness_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_MicroRichness_ess.se)

MicroDiversity_data <- data_soil %>% filter(
  VariableNew == "MicroDiversity", Irrigation %in% c("No")) %>%
  mutate(yi = if_else(Variable %in% c("Simpson"), -yi, yi))
sapply(MicroDiversity_data[col_names], function(x) length(unique(x)))
V_MicroDiversity <- calc.vbd(data = MicroDiversity_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_MicroDiversity <- rma.mv(
  yi, V_MicroDiversity, data = MicroDiversity_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_MicroDiversity)
model_MicroDiversity_ess.se <- rma.mv(
  yi, V_MicroDiversity, data = MicroDiversity_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_MicroDiversity_ess.se)

model_list_Soil <- list(
  BD = model_BD,
  pH = model_pH,
  EC = model_EC,
  SOC = model_SOC,
  TN = model_TN,
  DOC = model_DOC,
  AvailN = model_AvailN,
  AvailP = model_AvailP,
  EnzyActi = model_EnzyActi,
  MicroBiomass = model_MicroBiomass,
  MicroDiversity = model_MicroDiversity)
results_Soil <- data.frame(model = factor(
  names(model_list_Soil), levels = rev(names(model_list_Soil))),
  estimate = sapply(model_list_Soil, coef),
  ci.lb = sapply(model_list_Soil, function(m) m$ci.lb),
  ci.ub = sapply(model_list_Soil, function(m) m$ci.ub),
  n = sapply(model_list_Soil, function(m) m$k)) %>%
  mutate(contains_zero = ci.lb <= 0 & ci.ub >= 0)

ggplot(results_Soil, aes(x = estimate, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), color = "#BBB74E", width = 0.15, linewidth = 1) +
  geom_point(aes(fill = contains_zero), shape = 21, size = 3, color = "#BBB74E", stroke = 1.5) +
  scale_fill_manual(values = c("FALSE" = "#BBB74E", "TRUE" = "white"), guide = "none", drop = FALSE) +
  geom_text(aes(label = paste0("(", n, ")")),
            x = max(abs(c(results_Soil$ci.lb, results_Soil$ci.ub))) * 1.4*0.995,
            hjust = 1, size = 12/2.8346, color = "black") +
  labs(x = "Effect Size", y = NULL) + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0)) +
  coord_cartesian(xlim = c(-max(abs(c(results_Soil$ci.lb, results_Soil$ci.ub))) * 1.4,
                           max(abs(c(results_Soil$ci.lb, results_Soil$ci.ub))) * 1.4)) +
  theme(aspect.ratio = 2.5/2) -> p_ES_Soil; p_ES_Soil
cairo_pdf("Output/ES-Soil20260623.pdf", bg = "transparent", width = 30/2.54, height = 10/2.54)
p_ES_Soil
dev.off()

neg_vars <- c("IGeffect", "IGeffectS", "BD", "CAT", "POX", "PER", "OX",
              "BFratio", "Q10", "Simpson", "nirK", "nirS", "N2O", "CO2")
data_combined_with_con7 <- data_combined_with_con6 %>%
  mutate(yi = if_else(Variable %in% neg_vars, -yi, yi)) %>%
  filter(VariableNew != "pH")

# Composite PGI/SQI/EMF effects ----
data_PGI <- data_combined_with_con7 %>%
  filter(Niche == "plant", Irrigation %in% c("No"))

data_SQI <- data_combined_with_con7 %>%
  filter(Niche != "plant" | is.na(PlantAppear)) %>%
  filter(Irrigation %in% c("No"))

data_EMF <- data_combined_with_con7 %>%
  filter(Irrigation %in% c("No"))

sapply(data_PGI[col_names], function(x) length(unique(x)))
V_PGI <- calc.vbd(data = data_PGI, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
model_PGI <- rma.mv(
  yi, V_PGI, data = data_PGI,
  random = list(~ 1 | ExpID, ~ 1 | VariableNew/SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_PGI)
model_PGI_ess.se <- rma.mv(
  yi, V_PGI, data = data_PGI,
  random = list(~ 1 | ExpID, ~ 1 | VariableNew/SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
summary(model_PGI_ess.se)

sapply(data_SQI[col_names], function(x) length(unique(x)))
V_SQI <- calc.vbd(data = data_SQI, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
start1 <- Sys.time()
model_SQI <- try(rma.mv(
  yi, V_SQI, data = data_SQI,
  random = list(~ 1 | ExpID, ~ 1 | VariableNew/SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000)))
end1 <- Sys.time()
end1 - start1
summary(model_SQI)
start2 <- Sys.time()
model_SQI_ess.se <- try(rma.mv(
  yi, V_SQI, data = data_SQI,
  random = list(~ 1 | ExpID, ~ 1 | VariableNew/SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000)))
end2 <- Sys.time()
end2 - start2
summary(model_SQI_ess.se)

sapply(data_EMF[col_names], function(x) length(unique(x)))
V_EMF <- calc.vbd(data = data_EMF, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
start3 <- Sys.time()
model_EMF <- try(rma.mv(
  yi, V_EMF, data = data_EMF,
  random = list(~ 1 | ExpID, ~ 1 | EMFCat/EMFFun/VariableNew/SameVarIDN),
  test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000)))
end3 <- Sys.time()
end3 - start3
summary(model_EMF)
start4 <- Sys.time()
model_EMF_ess.se <- try(rma.mv(
  yi, V_EMF, data = data_EMF,
  random = list(~ 1 | ExpID, ~ 1 | EMFCat/EMFFun/VariableNew/SameVarIDN),
  mods = ~ ess.se, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000)))
end4 <- Sys.time()
end4 - start4
summary(model_EMF_ess.se)

model_list_Index <- list(
  PGI = model_PGI,
  SQI = model_SQI,
  EMF = model_EMF)
results_Index <- data.frame(model = factor(
  names(model_list_Index), levels = rev(names(model_list_Index))),
  estimate = sapply(model_list_Index, coef),
  ci.lb = sapply(model_list_Index, function(m) m$ci.lb),
  ci.ub = sapply(model_list_Index, function(m) m$ci.ub),
  n = sapply(model_list_Index, function(m) m$k)) %>%
  mutate(contains_zero = ci.lb <= 0 & ci.ub >= 0)

new1 <- rbind(
  results_Plant[results_Plant$model == "Yield", ],
  results_Index[results_Index$model %in% c("SQI", "EMF"), ])

orch_cols <- c("Yield" = "#65C2AD", "SQI" = "#BBB74E", "EMF" = "#95B7DA")
new1$fill_color <- ifelse(new1$contains_zero, "white", orch_cols[new1$model])
ggplot(new1, aes(x = estimate, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub, color = model), width = 0.15, linewidth = 1) +
  geom_point(aes(fill = fill_color, color = model), shape = 21, size = 3, stroke = 1.5) +
  scale_fill_identity(guide = "none") +
  geom_text(aes(label = paste0("(", n, ")")),
            x = max(abs(c(new1$ci.lb, new1$ci.ub))) * 1.2,
            hjust = 0.5, size = 13/2.8346, color = "black") +
  scale_y_discrete(limits = c("EMF", "SQI", "Yield")) +
  scale_color_manual(values = orch_cols) +
  labs(x = "Effect Size", y = NULL) + theme_classic() + theme(
    axis.line = element_blank(), panel.border = element_rect(),
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0)) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(-max(abs(c(new1$ci.lb, new1$ci.ub))) * 1.4,
                           max(abs(c(new1$ci.lb, new1$ci.ub))) * 1.4)) +
  theme(aspect.ratio = 3.5/2) -> p_ES_Index; p_ES_Index
pdf("Output/ES-Index20260415.pdf",
    height = 12/2.54)
p_ES_Index
dev.off()

new2 <- rbind(
  results_Index[results_Index$model == "PGI", ],
  results_Plant[results_Plant$model != "Yield", ])

ggplot(new2, aes(x = estimate, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), color = "#65C2AD", width = 0.15, linewidth = 1) +
  geom_point(aes(fill = contains_zero), shape = 21, size = 3, color = "#65C2AD", stroke = 1.5) +
  scale_fill_manual(values = c("FALSE" = "#65C2AD", "TRUE" = "white"), guide = "none", drop = FALSE) +
  geom_text(aes(label = paste0("(", n, ")")),
            x = max(abs(c(new2$ci.lb, new2$ci.ub))) * 1.4*0.995,
            hjust = 1, size = 12/2.8346, color = "black") +
  scale_y_discrete(limits = c("RootLength", "Chlorophyll", "Height", "Biomass", "PGI")) +
  labs(x = "Effect Size", y = NULL) + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0)) +
  coord_cartesian(xlim = c(-max(abs(c(new2$ci.lb, new2$ci.ub))) * 1.4,
                           max(abs(c(new2$ci.lb, new2$ci.ub))) * 1.4)) +
  theme(aspect.ratio = 2.5/2) -> p_Plant_Index; p_Plant_Index
cairo_pdf("Output/ES-Plant20260623.pdf", bg = "transparent", width = 30/2.54, height = 10/2.54)
p_Plant_Index
dev.off()

# Fig. S6: publication-bias funnel plots ----
generate_funnel_plot <- function(
    model, model_ess.se,
    filename_prefix,
    effect_label = "Effect size",
    height_cm = 7.32, width_cm = 7.32,
    yaxis = "sei",
    levels = c(90, 95, 99), shades = c("white", "gray75", "gray55"),
    refline = 0, digits = 3) {

  pdf_filename <- paste0("Funnel-", filename_prefix, ".pdf")
  if (!dir.exists("Output")) {dir.create("Output")}
  pdf_filename <- file.path("Output", paste0("Funnel-", filename_prefix, ".pdf"))

  pdf(pdf_filename, height = height_cm/2.54, width = width_cm/2.54)
  par(mar = c(2, 2, 0.5, 0.5), mgp = c(1.2, 0.3, 0), cex = 1.1)
  funnel(model, xlab = effect_label,
         yaxis = yaxis, level = levels,
         shade = shades, refline = refline)
  ess.se_p <- summary(model_ess.se)$pval[2]
  legend("topright",
         legend = bquote(italic(P)~"= "~.(formatC(ess.se_p, format = "f", digits = digits))),
         bty = "n", text.col = "red", cex = 1, inset = 0.00)
  dev.off()

  message("Generated: ", pdf_filename)
  return(pdf_filename)
}

batch_generate_funnel <- function(indicator_list, suffix = NULL) {
  generated_files <- c()
  for (indicator in indicator_list) {
    model_name <- paste0("model_", indicator)
    ess.se_name <- paste0("model_", indicator, "_ess.se")
    if (!exists(model_name) || !exists(ess.se_name)) {
      warning(paste("model ", model_name, " or ", ess.se_name, " not found; skipped ", indicator))
      next}
    model_obj <- get(model_name)
    model_ess.se_obj <- get(ess.se_name)
    filename_prefix <- paste0(indicator, suffix)
    filename <- generate_funnel_plot(
      model = model_obj,
      model_ess.se = model_ess.se_obj,
      filename_prefix = filename_prefix)
    generated_files <- c(generated_files, filename)
  }
  message("\nBatch generation completed. Files generated: ", length(generated_files), "  files")
  return(generated_files)
}

my_indicators <-
  c("Chlorophyll", "PlantHeightAbove", "BiomAbove", "Yield", "RootLength",
    "MicroDiversity", "MicroBiomass", "EnzyActi",
    "AvailP", "AvailN", "DOC", "TN", "SOC", "EC",
    "pH", "BD", "EMF", "SQI", "PGI",
    "CAcqActi", "NAcqActi", "PAcqActi", "OXActi", "MicroRichness")
batch_generate_funnel(my_indicators, suffix = "20260316")

pacman::p_load(metafor, future.apply)

data_ExpID <- read.xlsx(meta_file, "基础1整理")
Yield_data$`Exp. ID` <- NA
Yield_data$`Exp. ID` <- data_ExpID$`Exp..ID`[match(Yield_data$ExpID, data_ExpID$`ExpID原`)]
data_SQI$`Exp. ID` <- NA
data_SQI$`Exp. ID` <- data_ExpID$`Exp..ID`[match(data_SQI$ExpID, data_ExpID$`ExpID原`)]
data_EMF$`Exp. ID` <- NA
data_EMF$`Exp. ID` <- data_ExpID$`Exp..ID`[match(data_EMF$ExpID, data_ExpID$`ExpID原`)]

plan(multisession, workers = availableCores() - 1)
start_Yield99 <- Sys.time()
exp_ids_Yield <- unique(Yield_data$`Exp. ID`)
# Fig. S7: leave-one-out sensitivity ----
resultsLOO_Yield <- future_lapply(seq_along(exp_ids_Yield), function(i) {
  data_loo <- Yield_data[Yield_data$`Exp. ID` != exp_ids_Yield[i], ]
  V_loo <- calc.vbd(
    data = data_loo, CommonID = CommonIDN, m1 = ConMean,
    sd1 = ConSD_new, n1 = ConN, vi = vi)
  model_loo <- rma.mv(
    yi, V_loo, data = data_loo,
    random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
    test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
  list(removed_exp = exp_ids_Yield[i],
       estimate = model_loo$beta[1],
       se = model_loo$se, tval = model_loo$zval,
       df = model_loo$ddf, pval = model_loo$pval,
       ci.lb = model_loo$ci.lb, ci.ub = model_loo$ci.ub,
       k = model_loo$k)
})
sensitivity_results_Yield <- do.call(rbind, lapply(resultsLOO_Yield, as.data.frame))
end_Yield99 <- Sys.time()
end_Yield99 - start_Yield99

sensitivity_results_Yield <- sensitivity_results_Yield %>% mutate(
  model = factor(removed_exp, levels = sort(removed_exp, decreasing = TRUE)),
  contains_zero = ci.lb <= 0 & ci.ub >= 0, n_label = paste0("(", k, ")"))

ggplot(sensitivity_results_Yield, aes(x = estimate, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), color = "#65C2AD", width = 0.15, linewidth = 1) +
  geom_point(aes(fill = contains_zero), shape = 21, size = 3, color = "#65C2AD", stroke = 1.5) +
  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = "#65C2AD"), guide = "none", drop = FALSE) +
  geom_text(aes(label = n_label),
            x = max(abs(c(sensitivity_results_Yield$ci.lb, sensitivity_results_Yield$ci.ub))) * 1.4 * 0.995,
            hjust = 1, size = 12/2.8346, color = "black") +
  labs(x = "Effect Size", y = "Removed Exp. ID") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0)) +
  coord_cartesian(xlim = c(-max(abs(c(sensitivity_results_Yield$ci.lb, sensitivity_results_Yield$ci.ub))) * 1.4,
                           max(abs(c(sensitivity_results_Yield$ci.lb, sensitivity_results_Yield$ci.ub))) * 1.4)) +
  theme(aspect.ratio = 3.5/2) -> p_LOO_Yield; p_LOO_Yield
cairo_pdf("Output/LOO-Yield20260623.pdf", bg = "transparent", width = 30/2.54, height = 15/2.54)
p_LOO_Yield
dev.off()

plan(multisession, workers = availableCores() - 1)
start_SQI99 <- Sys.time()
exp_ids_SQI <- unique(data_SQI$`Exp. ID`)
resultsLOO_SQI <- future_lapply(seq_along(exp_ids_SQI), function(i) {
  data_loo <- data_SQI[data_SQI$`Exp. ID` != exp_ids_SQI[i], ]
  V_loo <- calc.vbd(
    data = data_loo, CommonID = CommonIDN, m1 = ConMean,
    sd1 = ConSD_new, n1 = ConN, vi = vi)
  model_loo <- rma.mv(
    yi, V_loo, data = data_loo,
    random = list(~ 1 | ExpID, ~ 1 | VariableNew/SameVarIDN),
    test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
  list(removed_exp = exp_ids_SQI[i],
       estimate = model_loo$beta[1],
       se = model_loo$se, tval = model_loo$zval,
       df = model_loo$ddf, pval = model_loo$pval,
       ci.lb = model_loo$ci.lb, ci.ub = model_loo$ci.ub,
       k = model_loo$k)
})
sensitivity_results_SQI <- do.call(rbind, lapply(resultsLOO_SQI, as.data.frame))
end_SQI99 <- Sys.time()
end_SQI99 - start_SQI99

sensitivity_results_SQI <- sensitivity_results_SQI %>% mutate(
  model = factor(removed_exp, levels = sort(removed_exp, decreasing = TRUE)),
  contains_zero = ci.lb <= 0 & ci.ub >= 0, n_label = paste0("(", k, ")"))

ggplot(sensitivity_results_SQI, aes(x = estimate, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), color = "#BBB74E", width = 0.15, linewidth = 1) +
  geom_point(aes(fill = contains_zero), shape = 21, size = 3, color = "#BBB74E", stroke = 1.5) +
  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = "#BBB74E"), guide = "none", drop = FALSE) +
  geom_text(aes(label = n_label),
            x = max(abs(c(sensitivity_results_SQI$ci.lb, sensitivity_results_SQI$ci.ub))) * 1.4 * 0.995,
            hjust = 1, size = 12/2.8346, color = "black") +
  labs(x = "Effect Size", y = "Removed Exp. ID") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0)) +
  coord_cartesian(xlim = c(-max(abs(c(sensitivity_results_SQI$ci.lb, sensitivity_results_SQI$ci.ub))) * 1.4,
                           max(abs(c(sensitivity_results_SQI$ci.lb, sensitivity_results_SQI$ci.ub))) * 1.4)) +
  theme(aspect.ratio = 3.5/2) -> p_LOO_SQI; p_LOO_SQI
cairo_pdf("Output/LOO-SQI20260623.pdf", bg = "transparent", width = 30/2.54, height = 15/2.54)
p_LOO_SQI
dev.off()

plan(multisession, workers = availableCores() - 1)
start_EMF99 <- Sys.time()
exp_ids_EMF <- unique(data_EMF$`Exp. ID`)
resultsLOO_EMF <- future_lapply(seq_along(exp_ids_EMF), function(i) {
  data_loo <- data_EMF[data_EMF$`Exp. ID` != exp_ids_EMF[i], ]
  V_loo <- calc.vbd(
    data = data_loo, CommonID = CommonIDN, m1 = ConMean,
    sd1 = ConSD_new, n1 = ConN, vi = vi)
  model_loo <- rma.mv(
    yi, V_loo, data = data_loo,
    random = list(~ 1 | ExpID, ~ 1 | EMFCat/EMFFun/VariableNew/SameVarIDN),
    test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
  list(removed_exp = exp_ids_EMF[i],
       estimate = model_loo$beta[1],
       se = model_loo$se, tval = model_loo$zval,
       df = model_loo$ddf, pval = model_loo$pval,
       ci.lb = model_loo$ci.lb, ci.ub = model_loo$ci.ub,
       k = model_loo$k)
})
sensitivity_results_EMF <- do.call(rbind, lapply(resultsLOO_EMF, as.data.frame))
end_EMF99 <- Sys.time()
end_EMF99 - start_EMF99

sensitivity_results_EMF <- sensitivity_results_EMF %>% mutate(
  model = factor(removed_exp, levels = sort(removed_exp, decreasing = TRUE)),
  contains_zero = ci.lb <= 0 & ci.ub >= 0, n_label = paste0("(", k, ")"))

ggplot(sensitivity_results_EMF, aes(x = estimate, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), color = "#95B7DA", width = 0.15, linewidth = 1) +
  geom_point(aes(fill = contains_zero), shape = 21, size = 3, color = "#95B7DA", stroke = 1.5) +
  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = "#95B7DA"), guide = "none", drop = FALSE) +
  geom_text(aes(label = n_label),
            x = max(abs(c(sensitivity_results_EMF$ci.lb, sensitivity_results_EMF$ci.ub))) * 1.4 * 0.995,
            hjust = 1, size = 12/2.8346, color = "black") +
  labs(x = "Effect Size", y = "Removed Exp. ID") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0)) +
  coord_cartesian(xlim = c(-max(abs(c(sensitivity_results_EMF$ci.lb, sensitivity_results_EMF$ci.ub))) * 1.4,
                           max(abs(c(sensitivity_results_EMF$ci.lb, sensitivity_results_EMF$ci.ub))) * 1.4)) +
  theme(aspect.ratio = 3.5/2) -> p_LOO_EMF; p_LOO_EMF
cairo_pdf("Output/LOO-EMF20260623.pdf", bg = "transparent", width = 30/2.54, height = 15/2.54)
p_LOO_EMF
dev.off()

check_vars_Yield <-
  c("SampleTime2", "Niche", "SoilDepth", "PlantAppear",
    "PlasType2", "PlasShape", "PlasSize2", "PlasDose",
    "Ecosystem", "AridityIndex", "MAT", "MAP",
    "FAO90", "AWC", "TEXTURE_USDA",
    "BULK", "ORG_CARBON", "PH_WATER", "TOTAL_N",
    "ELEC_COND")

valid_vars_Yield <- c()
for (var in check_vars_Yield) {
  if (var %in% names(Yield_data)) {
    unique_vals <- unique(Yield_data[[var]])
    valid_vals <- unique_vals[!is.na(unique_vals)]
    if (length(valid_vals) >= 2) {valid_vars_Yield <- c(valid_vars_Yield, var)}}
}

Yield_data$ExpID <- as.factor(Yield_data$ExpID)
Yield_data[valid_vars_Yield] <- lapply(Yield_data[valid_vars_Yield], function(x) {
  if (is.character(x)) as.factor(x) else x})

valid_vars_Yield
start_Yield1 <- Sys.time()
set.seed(99)
# Fig. S8: predictor importance and continuous moderators ----
mf_rep_Yield <- MetaForest(
  yi ~ ., data = Yield_data[, c("ExpID", "yi", "vi", valid_vars_Yield)],
  study = "ExpID", num.trees = 5000)
plot(mf_rep_Yield)
VarImpPlot(mf_rep_Yield$forest)
var_imp_Yield <- mf_rep_Yield$forest$variable.importance
var_imp_Yield_sorted <- sort(var_imp_Yield, decreasing = TRUE)
print(head(var_imp_Yield_sorted, 20))
preselected_Yield0 <- preselect(mf_rep_Yield, replications = 100, algorithm = "replicate")
plot(preselected_Yield0)
retain_mods_Yield0 <- preselect_vars(preselected_Yield0, cutoff = 0)
retain_mods_Yield0
end_Yield1 <- Sys.time()
end_Yield1 - start_Yield1

retain_mods_Yield0
start_Yield2 <- Sys.time()
set.seed(999)
mf_step_Yield <- MetaForest(
  yi ~ ., data = Yield_data[, c("ExpID", "yi", "vi", retain_mods_Yield0)],
  study = "ExpID", num.trees = 5000)
plot(mf_step_Yield)
VarImpPlot(mf_step_Yield$forest)
var_imp_step_Yield <- mf_step_Yield$forest$variable.importance
var_imp_step_Yield_sorted <- sort(var_imp_step_Yield, decreasing = TRUE)
print(head(var_imp_step_Yield_sorted, 20))
preselected_Yield <- preselect(mf_step_Yield, replications = 100, algorithm = "recursive")
plot(preselected_Yield)
retain_mods_Yield <- preselect_vars(preselected_Yield, cutoff = .3)
retain_mods_Yield
end_Yield2 <- Sys.time()
end_Yield2 - start_Yield2

retain_mods_Yield
all_models_Yield <- list()
for (k in 0:length(retain_mods_Yield)) {
  combos_Yield <- combn(retain_mods_Yield, k, simplify = FALSE)
  for (combo in combos_Yield) {
    if (length(combo) == 0) {formula_str <- "~ 1"} else {
      formula_str <- paste("~", paste(combo, collapse = " + "))}
    all_models_Yield[[length(all_models_Yield) + 1]] <- formula_str}
}
V_Yield <- calc.vbd(data = Yield_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
fit_model_Yield <- function(formula_str, data, V) {
  formula <- as.formula(paste("yi", formula_str))
  rma.mv(formula, V, data = data,
         random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
         test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
}
plan(multisession, workers = availableCores() - 1)
start_Yield39 <- Sys.time()
model_results_Yield <- future_lapply(all_models_Yield, function(formula_str) {
  tryCatch({
    mod <- fit_model_Yield(formula_str, Yield_data, V_Yield)
    data.frame(
      formula = formula_str, AICc = AICc(mod),
      logLik = logLik(mod), converged = TRUE)
  }, error = function(e) {
    data.frame(
      formula = formula_str, AICc = NA,
      logLik = NA, converged = FALSE)})
}, future.seed = TRUE)
model_results_Yield <- do.call(rbind, model_results_Yield)
end_Yield39 <- Sys.time()
end_Yield39 - start_Yield39
plan(sequential)

model_results_Yield <- model_results_Yield[order(model_results_Yield$AICc, na.last = TRUE), ]
model_results_Yield$delta_AICc <- model_results_Yield$AICc - min(model_results_Yield$AICc, na.rm = TRUE)
model_results_Yield$weight <- exp(-0.5 * model_results_Yield$delta_AICc)
model_results_Yield$weight <- model_results_Yield$weight / sum(model_results_Yield$weight, na.rm = TRUE)

model_results_Yield_top_100 <- head(model_results_Yield, 100)
model_results_Yield_top_100$weight <- exp(-0.5 * model_results_Yield_top_100$delta_AICc)
model_results_Yield_top_100$weight <- model_results_Yield_top_100$weight / sum(model_results_Yield_top_100$weight, na.rm = TRUE)
importance_df_Yield_top_100 <- data.frame(
  variable = retain_mods_Yield,
  importance = sapply(retain_mods_Yield, function(var) {
    idx <- grepl(var, model_results_Yield_top_100$formula)
    sum(model_results_Yield_top_100$weight[idx], na.rm = TRUE)}))
importance_df_Yield_top_100 <- importance_df_Yield_top_100[order(importance_df_Yield_top_100$importance, decreasing = TRUE), ]

full_importance_Yield <- data.frame(
  Variable = retain_mods_Yield,
  Importance = sapply(retain_mods_Yield, function(var) {
    if (var %in% c("PlasShape", "PlasSize2", "WRB_PHASES", "FAO90")) {
      rows <- grep(paste0("^", var), importance_df_Yield_top_100$variable)
      if (length(rows) > 0) max(importance_df_Yield_top_100$importance[rows]) else 0
    } else {if (var %in% importance_df_Yield_top_100$variable) {
      importance_df_Yield_top_100$importance[importance_df_Yield_top_100$variable == var]} else 0}})
)
full_importance_Yield <- full_importance_Yield %>% arrange(desc(Importance))

full_importance_Yield$Variable
full_importance_Yield$Variable_New <-
  c("Dose", "MAP", "AI", "AWC", "MAT", "Shape", "Size")
full_importance_Yield$Variable_New <- factor(
  full_importance_Yield$Variable_New, levels = rev(unique(full_importance_Yield$Variable_New)))

ggplot(full_importance_Yield, aes(y = Variable_New, x = Importance)) +
  geom_bar(aes(fill = Importance), fill = "#65C2AD", width = 0.7, stat = "identity", color = NA) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 1.1), expand = c(0, 0), breaks = c(0, 0.5, 1)) +
  labs(y = NULL, x = "Importance") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0)) +
  theme(aspect.ratio = 0.8/1) -> p_VarImp_Yield; p_VarImp_Yield
cairo_pdf("Output/VarImp-Yield20260623.pdf", bg = "transparent", width = 30/2.54, height = 6.7/2.54)
p_VarImp_Yield
dev.off()

check_vars_SQI <-
  c("SampleTime2", "Niche", "SoilDepth", "PlantAppear",
    "PlasType2", "PlasShape", "PlasSize2", "PlasDose",
    "Ecosystem", "AridityIndex", "MAT", "MAP",
    "FAO90", "AWC", "TEXTURE_USDA",
    "BULK", "ORG_CARBON", "PH_WATER", "TOTAL_N",
    "ELEC_COND")

valid_vars_SQI <- c()
for (var in check_vars_SQI) {
  if (var %in% names(data_SQI)) {
    unique_vals <- unique(data_SQI[[var]])
    valid_vals <- unique_vals[!is.na(unique_vals)]
    if (length(valid_vals) >= 2) {valid_vars_SQI <- c(valid_vars_SQI, var)}}
}

data_SQI$ExpID <- as.factor(data_SQI$ExpID)
data_SQI[valid_vars_SQI] <- lapply(data_SQI[valid_vars_SQI], function(x) {
  if (is.character(x)) as.factor(x) else x})

valid_vars_SQI
start_SQI1 <- Sys.time()
set.seed(99)
mf_rep_SQI <- MetaForest(
  yi ~ ., data = data_SQI[, c("ExpID", "yi", "vi", valid_vars_SQI)],
  study = "ExpID", num.trees = 5000)
plot(mf_rep_SQI)
VarImpPlot(mf_rep_SQI$forest)
var_imp_SQI <- mf_rep_SQI$forest$variable.importance
var_imp_SQI_sorted <- sort(var_imp_SQI, decreasing = TRUE)
print(head(var_imp_SQI_sorted, 20))
preselected_SQI0 <- preselect(mf_rep_SQI, replications = 100, algorithm = "replicate")
plot(preselected_SQI0)
retain_mods_SQI0 <- preselect_vars(preselected_SQI0, cutoff = 0)
retain_mods_SQI0
end_SQI1 <- Sys.time()
end_SQI1 - start_SQI1

retain_mods_SQI0
start_SQI2 <- Sys.time()
set.seed(999)
mf_step_SQI <- MetaForest(
  yi ~ ., data = data_SQI[, c("ExpID", "yi", "vi", retain_mods_SQI0)],
  study = "ExpID", num.trees = 5000)
plot(mf_step_SQI)
VarImpPlot(mf_step_SQI$forest)
var_imp_step_SQI <- mf_step_SQI$forest$variable.importance
var_imp_step_SQI_sorted <- sort(var_imp_step_SQI, decreasing = TRUE)
print(head(var_imp_step_SQI_sorted, 20))
preselected_SQI <- preselect(mf_step_SQI, replications = 100, algorithm = "recursive")
plot(preselected_SQI)
retain_mods_SQI <- preselect_vars(preselected_SQI, cutoff = .3)
retain_mods_SQI
end_SQI2 <- Sys.time()
end_SQI2 - start_SQI2

retain_mods_SQI
all_models_SQI <- list()
for (k in 0:length(retain_mods_SQI)) {
  combos_SQI <- combn(retain_mods_SQI, k, simplify = FALSE)
  for (combo in combos_SQI) {
    if (length(combo) == 0) {formula_str <- "~ 1"} else {
      formula_str <- paste("~", paste(combo, collapse = " + "))}
    all_models_SQI[[length(all_models_SQI) + 1]] <- formula_str}
}
V_SQI <- calc.vbd(data = data_SQI, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
fit_model_SQI <- function(formula_str, data, V) {
  formula <- as.formula(paste("yi", formula_str))
  rma.mv(formula, V, data = data,
         random = list(~ 1 | ExpID, ~ 1 | VariableNew/SameVarIDN),
         test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
}
plan(multisession, workers = availableCores() - 1)
start_SQI39 <- Sys.time()
model_results_SQI <- future_lapply(all_models_SQI, function(formula_str) {
  tryCatch({
    mod <- fit_model_SQI(formula_str, data_SQI, V_SQI)
    data.frame(
      formula = formula_str, AICc = AICc(mod),
      logLik = logLik(mod), converged = TRUE)
  }, error = function(e) {
    data.frame(
      formula = formula_str, AICc = NA,
      logLik = NA, converged = FALSE)})
}, future.seed = TRUE)
model_results_SQI <- do.call(rbind, model_results_SQI)
end_SQI39 <- Sys.time()
end_SQI39 - start_SQI39
plan(sequential)

model_results_SQI <- model_results_SQI[order(model_results_SQI$AICc, na.last = TRUE), ]
model_results_SQI$delta_AICc <- model_results_SQI$AICc - min(model_results_SQI$AICc, na.rm = TRUE)
model_results_SQI$weight <- exp(-0.5 * model_results_SQI$delta_AICc)
model_results_SQI$weight <- model_results_SQI$weight / sum(model_results_SQI$weight, na.rm = TRUE)

model_results_SQI_top_100 <- head(model_results_SQI, 100)
model_results_SQI_top_100$weight <- exp(-0.5 * model_results_SQI_top_100$delta_AICc)
model_results_SQI_top_100$weight <- model_results_SQI_top_100$weight / sum(model_results_SQI_top_100$weight, na.rm = TRUE)
importance_df_SQI_top_100 <- data.frame(
  variable = retain_mods_SQI,
  importance = sapply(retain_mods_SQI, function(var) {
    idx <- grepl(var, model_results_SQI_top_100$formula)
    sum(model_results_SQI_top_100$weight[idx], na.rm = TRUE)}))
importance_df_SQI_top_100 <- importance_df_SQI_top_100[order(importance_df_SQI_top_100$importance, decreasing = TRUE), ]

full_importance_SQI <- data.frame(
  Variable = retain_mods_SQI,
  Importance = sapply(retain_mods_SQI, function(var) {
    if (var %in% c("PlasShape", "PlasSize2", "WRB_PHASES", "FAO90")) {
      rows <- grep(paste0("^", var), importance_df_SQI_top_100$variable)
      if (length(rows) > 0) max(importance_df_SQI_top_100$importance[rows]) else 0
    } else {if (var %in% importance_df_SQI_top_100$variable) {
      importance_df_SQI_top_100$importance[importance_df_SQI_top_100$variable == var]} else 0}})
)
full_importance_SQI <- full_importance_SQI %>% arrange(desc(Importance))

full_importance_SQI$Variable
full_importance_SQI$Variable_New <-
  c("Type", "Dose", "Duration", "MAT", "MAP", "Niche")
full_importance_SQI$Variable_New <- factor(
  full_importance_SQI$Variable_New, levels = rev(unique(full_importance_SQI$Variable_New)))

ggplot(full_importance_SQI, aes(y = Variable_New, x = Importance)) +
  geom_bar(aes(fill = Importance), fill = "#BBB74E", width = 0.7, stat = "identity", color = NA) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 1.1), expand = c(0, 0), breaks = c(0, 0.5, 1)) +
  labs(y = NULL, x = "Importance") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0)) +
  theme(aspect.ratio = 0.8/1) -> p_VarImp_SQI; p_VarImp_SQI
cairo_pdf("Output/VarImp-SQI20260623.pdf", bg = "transparent", width = 30/2.54, height = 6.7/2.54)
p_VarImp_SQI
dev.off()

check_vars_EMF <-
  c("SampleTime2", "Niche", "SoilDepth", "PlantAppear",
    "PlasType2", "PlasShape", "PlasSize2", "PlasDose",
    "Ecosystem", "AridityIndex", "MAT", "MAP",
    "FAO90", "AWC", "TEXTURE_USDA",
    "BULK", "ORG_CARBON", "PH_WATER", "TOTAL_N",
    "ELEC_COND")

valid_vars_EMF <- c()
for (var in check_vars_EMF) {
  if (var %in% names(data_EMF)) {
    unique_vals <- unique(data_EMF[[var]])
    valid_vals <- unique_vals[!is.na(unique_vals)]
    if (length(valid_vals) >= 2) {valid_vars_EMF <- c(valid_vars_EMF, var)}}
}

data_EMF$ExpID <- as.factor(data_EMF$ExpID)
data_EMF[valid_vars_EMF] <- lapply(data_EMF[valid_vars_EMF], function(x) {
  if (is.character(x)) as.factor(x) else x})

valid_vars_EMF
start_EMF1 <- Sys.time()
set.seed(99)
mf_rep_EMF <- MetaForest(
  yi ~ ., data = data_EMF[, c("ExpID", "yi", "vi", valid_vars_EMF)],
  study = "ExpID", num.trees = 5000)
plot(mf_rep_EMF)
VarImpPlot(mf_rep_EMF$forest)
var_imp_EMF <- mf_rep_EMF$forest$variable.importance
var_imp_EMF_sorted <- sort(var_imp_EMF, decreasing = TRUE)
print(head(var_imp_EMF_sorted, 20))
preselected_EMF0 <- preselect(mf_rep_EMF, replications = 100, algorithm = "replicate")
plot(preselected_EMF0)
retain_mods_EMF0 <- preselect_vars(preselected_EMF0, cutoff = 0)
retain_mods_EMF0
end_EMF1 <- Sys.time()
end_EMF1 - start_EMF1

retain_mods_EMF0
start_EMF2 <- Sys.time()
set.seed(99)
mf_step_EMF <- MetaForest(
  yi ~ ., data = data_EMF[, c("ExpID", "yi", "vi", retain_mods_EMF0)],
  study = "ExpID", num.trees = 5000)
plot(mf_step_EMF)
VarImpPlot(mf_step_EMF$forest)
var_imp_step_EMF <- mf_step_EMF$forest$variable.importance
var_imp_step_EMF_sorted <- sort(var_imp_step_EMF, decreasing = TRUE)
print(head(var_imp_step_EMF_sorted, 20))
preselected_EMF <- preselect(mf_step_EMF, replications = 100, algorithm = "recursive")
plot(preselected_EMF)
retain_mods_EMF <- preselect_vars(preselected_EMF, cutoff = .3)
retain_mods_EMF
end_EMF2 <- Sys.time()
end_EMF2 - start_EMF2

retain_mods_EMF
all_models_EMF <- list()
for (k in 0:length(retain_mods_EMF)) {
  combos_EMF <- combn(retain_mods_EMF, k, simplify = FALSE)
  for (combo in combos_EMF) {
    if (length(combo) == 0) {formula_str <- "~ 1"} else {
      formula_str <- paste("~", paste(combo, collapse = " + "))}
    all_models_EMF[[length(all_models_EMF) + 1]] <- formula_str}
}
V_EMF <- calc.vbd(data = data_EMF, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)
fit_model_EMF <- function(formula_str, data, V) {
  formula <- as.formula(paste("yi", formula_str))
  rma.mv(formula, V, data = data,
         random = list(~ 1 | ExpID, ~ 1 | EMFCat/EMFFun/VariableNew/SameVarIDN),
         test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))
}
plan(multisession, workers = availableCores() - 1)
start_EMF39 <- Sys.time()
model_results_EMF <- future_lapply(all_models_EMF, function(formula_str) {
  tryCatch({
    mod <- fit_model_EMF(formula_str, data_EMF, V_EMF)
    data.frame(
      formula = formula_str, AICc = AICc(mod),
      logLik = logLik(mod), converged = TRUE)
  }, error = function(e) {
    data.frame(
      formula = formula_str, AICc = NA,
      logLik = NA, converged = FALSE)})
}, future.seed = TRUE)
model_results_EMF <- do.call(rbind, model_results_EMF)
end_EMF39 <- Sys.time()
end_EMF39 - start_EMF39
plan(sequential)

model_results_EMF <- model_results_EMF[order(model_results_EMF$AICc, na.last = TRUE), ]
model_results_EMF$delta_AICc <- model_results_EMF$AICc - min(model_results_EMF$AICc, na.rm = TRUE)
model_results_EMF$weight <- exp(-0.5 * model_results_EMF$delta_AICc)
model_results_EMF$weight <- model_results_EMF$weight / sum(model_results_EMF$weight, na.rm = TRUE)

model_results_EMF_top_100 <- head(model_results_EMF, 100)
model_results_EMF_top_100$weight <- exp(-0.5 * model_results_EMF_top_100$delta_AICc)
model_results_EMF_top_100$weight <- model_results_EMF_top_100$weight / sum(model_results_EMF_top_100$weight, na.rm = TRUE)
importance_df_EMF_top_100 <- data.frame(
  variable = retain_mods_EMF,
  importance = sapply(retain_mods_EMF, function(var) {
    idx <- grepl(var, model_results_EMF_top_100$formula)
    sum(model_results_EMF_top_100$weight[idx], na.rm = TRUE)}))
importance_df_EMF_top_100 <- importance_df_EMF_top_100[order(importance_df_EMF_top_100$importance, decreasing = TRUE), ]

full_importance_EMF <- data.frame(
  Variable = retain_mods_EMF,
  Importance = sapply(retain_mods_EMF, function(var) {
    if (var %in% c("PlasShape", "PlasSize2", "WRB_PHASES", "FAO90")) {
      rows <- grep(paste0("^", var), importance_df_EMF_top_100$variable)
      if (length(rows) > 0) max(importance_df_EMF_top_100$importance[rows]) else 0
    } else {if (var %in% importance_df_EMF_top_100$variable) {
      importance_df_EMF_top_100$importance[importance_df_EMF_top_100$variable == var]} else 0}})
)
full_importance_EMF <- full_importance_EMF %>% arrange(desc(Importance))

full_importance_EMF$Variable
full_importance_EMF$Variable_New <-
  c("Plant", "Type", "Size", "Dose", "Duration", "MAT", "MAP", "Ecosystem")
full_importance_EMF$Variable_New <- factor(
  full_importance_EMF$Variable_New, levels = rev(unique(full_importance_EMF$Variable_New)))

ggplot(full_importance_EMF, aes(y = Variable_New, x = Importance)) +
  geom_bar(aes(fill = Importance), fill = "#95B7DA", width = 0.7, stat = "identity", color = NA) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 1.1), expand = c(0, 0), breaks = c(0, 0.5, 1)) +
  labs(y = NULL, x = "Importance") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0)) +
  theme(aspect.ratio = 0.8/1) -> p_VarImp_EMF; p_VarImp_EMF
cairo_pdf("Output/VarImp-EMF20260623.pdf", bg = "transparent", width = 30/2.54, height = 6.7/2.54)
p_VarImp_EMF
dev.off()

sapply(Yield_data[col_names], function(x) length(unique(x)))
V_Yield <- calc.vbd(data = Yield_data, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)

# Fig. 2 and Fig. S9: subgroup and polymer-specific analyses ----
model_Yield_PlasType2 <- rma.mv(
  yi, V_Yield, data = Yield_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ PlasType2 - 1, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))

model_Yield_PlasSize2 <- rma.mv(
  yi, V_Yield, data = Yield_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ PlasSize2 - 1, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))

model_Yield_PlasDose <- rma.mv(
  yi, V_Yield, data = Yield_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ PlasDose, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))

model_Yield_SampleTime2 <- rma.mv(
  yi, V_Yield, data = Yield_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ SampleTime2, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))

model_Yield_MAT <- rma.mv(
  yi, V_Yield, data = Yield_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ MAT, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))

model_Yield_PlasShape <- rma.mv(
  yi, V_Yield, data = Yield_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ PlasShape - 1, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))

model_Yield_PlasType <- rma.mv(
  yi, V_Yield, data = Yield_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ PlasType - 1, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))

model_Yield_PlasDose2 <- rma.mv(
  yi, V_Yield, data = Yield_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ PlasDose2 - 1, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))

model_Yield_PlasSize <- rma.mv(
  yi, V_Yield, data = Yield_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ PlasSize, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))

model_Yield_MAP <- rma.mv(
  yi, V_Yield, data = Yield_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ MAP, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))

model_Yield_AridityIndex <- rma.mv(
  yi, V_Yield, data = Yield_data,
  random = list(~ 1 | ExpID, ~ 1 | SameVarIDN),
  mods = ~ AridityIndex, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000))

extract_model_results <- function(model, group_names) {
  summ <- summary(model)
  n_values <- if (!is.null(model$X)) {
    colSums(model$X)} else {rep(NA, length(group_names))}
  data.frame(
    Group = group_names, estimate = as.numeric(summ$b),
    ci.lb = summ$ci.lb, ci.ub = summ$ci.ub, n = n_values,
    contains_zero = (summ$ci.lb <= 0 & summ$ci.ub >= 0),
    QMp = summ$QMp)
}

summary(model_Yield_PlasType2)
summary(model_Yield_PlasSize2)
summary(model_Yield_PlasDose2)
summary(model_Yield_PlasShape)
summary(model_Yield_PlasType)

unique(data_combined_with_con7$PlasType2)
unique(data_combined_with_con7$PlasSize2)
unique(data_combined_with_con7$PlasDose2)
unique(data_combined_with_con7$PlasShape)
unique(data_combined_with_con7$PlasType)

df_Yield_PlasType2 <- extract_model_results(
  model_Yield_PlasType2, c("Biodegradable", "NonBiodegradable"))
df_Yield_PlasSize2 <- extract_model_results(
  model_Yield_PlasSize2, c("<0.1", ">5", "0.1-1", "unknown"))
df_Yield_PlasDose2 <- extract_model_results(
  model_Yield_PlasDose2, c("<0.02", "0.02-0.2", "0.2-1"))
df_Yield_PlasShape <- extract_model_results(
  model_Yield_PlasShape, c("fiber", "fragment", "granule"))
df_Yield_PlasType <- extract_model_results(
  model_Yield_PlasType,
  c("BioMix", "PE", "PLA", "PP"))

result_Sub_Yield <- rbind(
  df_Yield_PlasType2[order(factor(
    df_Yield_PlasType2$Group, levels = c("NonBiodegradable", "Biodegradable", "Mixture"))), ],
  df_Yield_PlasSize2[order(factor(
    df_Yield_PlasSize2$Group, levels = c("<0.1", "0.1-1", "1-5", ">5", "unknown"))), ],
  df_Yield_PlasDose2[order(factor(
    df_Yield_PlasDose2$Group, levels = c("<0.02", "0.02-0.2", "0.2-1", ">1", "unknown"))), ],
  df_Yield_PlasShape[order(factor(
    df_Yield_PlasShape$Group, levels = c("granule", "fiber", "fragment", "unknown"))), ]
); rownames(result_Sub_Yield) <- NULL
result_Sub_Yield

blank_row <- function(group) {
  data.frame(Group = group, estimate = NA, ci.lb = NA, ci.ub = NA, n = NA, contains_zero = NA, QMp = NA)}
result_Sub_Yield <- rbind(
  blank_row("PlasType"), result_Sub_Yield[1:2, ], blank_row("Mixture"),
  blank_row("PlasSize"), result_Sub_Yield[3:4, ], blank_row("1-5"), result_Sub_Yield[5:6, ],
  blank_row("PlasDose"), result_Sub_Yield[7:9, ], blank_row(">1"), blank_row("unknownD"),
  blank_row("PlasShape"), result_Sub_Yield[10:12, ], blank_row("unknownS"),
  blank_row(c("Plant", "noplant", "withplant"))
); rownames(result_Sub_Yield) <- NULL
result_Sub_Yield

result_Sub_Yield <- result_Sub_Yield %>%
  mutate(Group1 = factor(
    c("PlasType", "NonBio", "Bio", "Mixture",
      "PlasSize", "≤0.1", "0.1-1", "1-5", ">5", "unknown",
      "PlasDose", "≤0.02", "0.02-0.2", "0.2-1", ">1", "unknownD",
      "PlasShape", "Granule", "Fiber", "Fragment", "unknownS",
      "Plant", "NoPlant", "WithPlant")))

result_Sub_Yield$Group1 <- factor(result_Sub_Yield$Group1, levels = rev(unique(result_Sub_Yield$Group1)))

summ <- summary(model_Yield)
overall_est <- summ$b[1]
overall_ci.lb <- summ$ci.lb[1]
overall_ci.ub <- summ$ci.ub[1]
total_n <- nrow(Yield_data)
overall_contains_zero <- (overall_ci.lb <= 0 & overall_ci.ub >= 0)
overall_row <- data.frame(
  Group = "Overall", estimate = overall_est,
  ci.lb = overall_ci.lb, ci.ub = overall_ci.ub, n = total_n,
  contains_zero = overall_contains_zero, QMp = NA_real_,
  Group1 = "Overall", stringsAsFactors = FALSE)

orig_levels <- as.character(unique(result_Sub_Yield$Group1))
result_Sub_Yield <- rbind(overall_row, result_Sub_Yield)
result_Sub_Yield$Group1 <- factor(result_Sub_Yield$Group1, levels = c(rev(orig_levels), "Overall"))

x_max_Yield <- with(result_Sub_Yield[!is.na(result_Sub_Yield$estimate), ], max(abs(c(ci.lb, ci.ub)), na.rm = TRUE))
x_limit_Yield <- x_max_Yield * 1.4

ggplot(result_Sub_Yield, aes(x = estimate, y = Group1)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3.5, ymax = 8.5, fill = "grey98") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 14.5, ymax = 20.5, fill = "grey98") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 24.5, ymax = Inf, fill = "grey98") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), color = "#65C2AD", width = 0.15, linewidth = 1) +
  geom_point(aes(fill = contains_zero), shape = 21, size = 3, color = "#65C2AD", stroke = 1.5) +
  scale_fill_manual(values = c("FALSE" = "#65C2AD", "TRUE" = "white"), guide = "none", drop = FALSE) +
  geom_text(aes(label = ifelse(!is.na(n), paste0("(", n, ")"), "")),
            x = x_limit_Yield * 0.995, hjust = 1, size = 12/2.8346, color = "black") +
  geom_text(data = subset(result_Sub_Yield, Group1 %in% c("NonBio", "≤0.1", "≤0.02", "Granule", "NoPlant")),
            aes(label = ifelse(QMp < 0.001, "italic('P') ~ '<' ~ 0.001", sprintf("italic('P') ~ '= %.3f'", QMp))),
            x = -x_limit_Yield * 0.995, nudge_y = 1, hjust = 0, vjust = 0.5, size = 13/2.8346, color = "black", parse = TRUE) +
  scale_y_discrete(labels = c(
    "PlasType" = "Type", "PlasSize" = "Size", "PlasDose" = "Dose",
    "PlasShape" = "Shape", "unknownD" = "unknown", "unknownS" = "unknown",
    "WithPlant" = "Presence", "NoPlant" = "Absence")) +
  labs(x = "Effect Size", y = NULL) + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    axis.text.y = element_text(size = 15, face = ifelse(levels(result_Sub_Yield$Group1) %in% c(
      "Overall", "PlasType", "PlasSize", "PlasDose", "PlasShape", "Plant"), "bold", "plain")),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 5/2.5) +
  coord_cartesian(xlim = c(-x_limit_Yield, x_limit_Yield)) -> p_Sub_Yield; p_Sub_Yield
cairo_pdf("Output/Sub-Yield20260623.pdf", bg = "transparent", width = 10.85/2.54, height = 30/2.54)
p_Sub_Yield
dev.off()

Yield_data$ci.lb <- Yield_data$yi - 1.96 * sqrt(Yield_data$vi)
Yield_data$ci.ub <- Yield_data$yi + 1.96 * sqrt(Yield_data$vi)
Yield_data <- Yield_data[order(Yield_data$yi), ]
Yield_data$StudyID <- factor(1:nrow(Yield_data), levels = 1:nrow(Yield_data))

x_max_Yield <- max(abs(c(Yield_data$ci.lb, Yield_data$ci.ub)), na.rm = TRUE)
x_limit_Yield <- x_max_Yield * 1.4

ggplot(Yield_data, aes(x = yi, y = StudyID)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), width = 0, color = "#65C2AD", alpha = 0.3, linewidth = 0.2) +
  geom_point(shape = 16, size = 0.7, color = "#65C2AD") +
  labs(x = "Effect Size", y = "Observations") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1.1/2.5) +
  coord_cartesian(xlim = c(-x_limit_Yield, x_limit_Yield)) -> p_Indi_Yield; p_Indi_Yield
cairo_pdf("Output/Indi-Yield20260623.pdf", bg = "transparent", width = 8.88/2.54, height = 30/2.54)
p_Indi_Yield
dev.off()

result_Sub_Yield_PlasType <- df_Yield_PlasType %>% arrange(factor(
  Group, levels = c("PE", "PP", "PS", "PVC", "PES", "PAN", "POM", "TP",
                    "PLA", "PHA", "PBS", "PBAT", "PCL", "BioMix", "Mixture"))) %>% `rownames<-`(NULL)
result_Sub_Yield_PlasType$Group
result_Sub_Yield_PlasType <- result_Sub_Yield_PlasType %>% mutate(Group1 = factor(
  c("PE", "PP", "PLA", "BioMix")))

result_Sub_Yield_PlasType$Group1 <- with(result_Sub_Yield_PlasType, factor(Group1, levels = rev(unique(Group1))))
x_max_Yield_PlasType <- with(result_Sub_Yield_PlasType[!is.na(result_Sub_Yield_PlasType$estimate), ], max(abs(c(ci.lb, ci.ub)), na.rm = TRUE))
x_limit_Yield_PlasType <- x_max_Yield_PlasType * 1.4

ggplot(result_Sub_Yield_PlasType, aes(x = estimate, y = Group1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), color = "#65C2AD", width = 0.15, linewidth = 1) +
  geom_point(aes(fill = contains_zero), shape = 21, size = 3, color = "#65C2AD", stroke = 1.5) +
  scale_fill_manual(values = c("FALSE" = "#65C2AD", "TRUE" = "white"), guide = "none", drop = FALSE) +
  geom_text(aes(label = ifelse(!is.na(n), paste0("(", n, ")"), "")),
            x = x_max_Yield_PlasType * 1.4*0.995, hjust = 1, size = 12/2.8346, color = "black") +
  coord_cartesian(xlim = c(-x_limit_Yield_PlasType, x_limit_Yield_PlasType)) +
  labs(x = "Effect Size", y = NULL) + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 3/2) -> p_Sub_Yield_PlasType; p_Sub_Yield_PlasType
cairo_pdf("Output/Sub-Yield-PlasType20260623.pdf", bg = "transparent", width = 30/2.54, height = 12/2.54)
p_Sub_Yield_PlasType
dev.off()

summary(model_Yield_PlasSize)
summary(model_Yield_PlasDose)
summary(model_Yield_SampleTime2)
summary(model_Yield_MAT)
summary(model_Yield_MAP)
summary(model_Yield_AridityIndex)

unique(data_combined_with_con7$PlasSize)
unique(data_combined_with_con7$PlasDose)
unique(data_combined_with_con7$SampleTime2)
unique(data_combined_with_con7$MAT)
unique(data_combined_with_con7$MAP)
unique(data_combined_with_con7$AridityIndex)

pred_Yield_data_PlasDose <- with(Yield_data, data.frame(
  PlasDose = seq(min(PlasDose, na.rm = TRUE), max(PlasDose, na.rm = TRUE), length.out = 400)))
pred_Yield_data_PlasDose <- cbind(pred_Yield_data_PlasDose, predict(model_Yield_PlasDose, newmods = pred_Yield_data_PlasDose$PlasDose))
r2_val_Yield_PlasDose <- r2_ml(model_Yield_PlasDose)[1]
p_val_Yield_PlasDose <- model_Yield_PlasDose$pval[2]
r2_text_Yield_PlasDose <- ifelse(r2_val_Yield_PlasDose < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_Yield_PlasDose))
p_formatted_Yield_PlasDose <- ifelse(p_val_Yield_PlasDose < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_Yield_PlasDose))

ggplot() +
  geom_point(data = Yield_data, aes(x = PlasDose, y = yi), size = 3, alpha = 0.1, color = "#65C2AD") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_Yield_data_PlasDose, fill = "#65C2AD", alpha = 0.25,
              aes(x = PlasDose, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_Yield_data_PlasDose, linetype = "solid",
            aes(x = PlasDose, y = pred), color = "#65C2AD", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_Yield_PlasDose, p_formatted_Yield_PlasDose, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "Plastic dose (%)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_Yield_PlasDose; p_Sub_Yield_PlasDose
cairo_pdf("Output/Sub-Yield-PlasDose20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_Yield_PlasDose
dev.off()

pred_Yield_data_SampleTime2 <- with(Yield_data, data.frame(
  SampleTime2 = seq(min(SampleTime2, na.rm = TRUE), max(SampleTime2, na.rm = TRUE), length.out = 400)))
pred_Yield_data_SampleTime2 <- cbind(pred_Yield_data_SampleTime2, predict(model_Yield_SampleTime2, newmods = pred_Yield_data_SampleTime2$SampleTime2))
r2_val_Yield_SampleTime2 <- r2_ml(model_Yield_SampleTime2)[1]
p_val_Yield_SampleTime2 <- model_Yield_SampleTime2$pval[2]
r2_text_Yield_SampleTime2 <- ifelse(r2_val_Yield_SampleTime2 < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_Yield_SampleTime2))
p_formatted_Yield_SampleTime2 <- ifelse(p_val_Yield_SampleTime2 < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_Yield_SampleTime2))

ggplot() +
  geom_point(data = Yield_data, aes(x = SampleTime2, y = yi), size = 3, alpha = 0.1, color = "#65C2AD") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_Yield_data_SampleTime2, fill = "#65C2AD", alpha = 0.25,
              aes(x = SampleTime2, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_Yield_data_SampleTime2, linetype = "twodash",
            aes(x = SampleTime2, y = pred), color = "#65C2AD", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_Yield_SampleTime2, p_formatted_Yield_SampleTime2, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "Duration (month)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_Yield_SampleTime2; p_Sub_Yield_SampleTime2
cairo_pdf("Output/Sub-Yield-SampleTime2 20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_Yield_SampleTime2
dev.off()

pred_Yield_data_MAT <- with(Yield_data, data.frame(
  MAT = seq(min(MAT, na.rm = TRUE), max(MAT, na.rm = TRUE), length.out = 400)))
pred_Yield_data_MAT <- cbind(pred_Yield_data_MAT, predict(model_Yield_MAT, newmods = pred_Yield_data_MAT$MAT))
r2_val_Yield_MAT <- r2_ml(model_Yield_MAT)[1]
p_val_Yield_MAT <- model_Yield_MAT$pval[2]
r2_text_Yield_MAT <- ifelse(r2_val_Yield_MAT < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_Yield_MAT))
p_formatted_Yield_MAT <- ifelse(p_val_Yield_MAT < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_Yield_MAT))

ggplot() +
  geom_point(data = Yield_data, aes(x = MAT, y = yi), size = 3, alpha = 0.1, color = "#65C2AD") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_Yield_data_MAT, fill = "#65C2AD", alpha = 0.25,
              aes(x = MAT, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_Yield_data_MAT, linetype = "solid",
            aes(x = MAT, y = pred), color = "#65C2AD", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_Yield_MAT, p_formatted_Yield_MAT, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "MAT (°C)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_Yield_MAT; p_Sub_Yield_MAT
cairo_pdf("Output/Sub-Yield-MAT20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_Yield_MAT
dev.off()

pacman::p_load(openxlsx, tidyverse, patchwork)
pred_Yield_data_PlasSize <- with(Yield_data, data.frame(
  PlasSize = seq(min(PlasSize, na.rm = TRUE), max(PlasSize, na.rm = TRUE), length.out = 400)))
pred_Yield_data_PlasSize <- cbind(pred_Yield_data_PlasSize, predict(model_Yield_PlasSize, newmods = pred_Yield_data_PlasSize$PlasSize))
r2_val_Yield_PlasSize <- r2_ml(model_Yield_PlasSize)[1]
p_val_Yield_PlasSize <- model_Yield_PlasSize$pval[2]
r2_text_Yield_PlasSize <- ifelse(r2_val_Yield_PlasSize < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_Yield_PlasSize))
p_formatted_Yield_PlasSize <- ifelse(p_val_Yield_PlasSize < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_Yield_PlasSize))

ggplot() +
  geom_point(data = Yield_data, aes(x = PlasSize, y = yi), size = 3, alpha = 0.1, color = "#65C2AD") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_Yield_data_PlasSize, fill = "#65C2AD", alpha = 0.25,
              aes(x = PlasSize, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_Yield_data_PlasSize, linetype = "twodash",
            aes(x = PlasSize, y = pred), color = "#65C2AD", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_Yield_PlasSize, p_formatted_Yield_PlasSize, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "Plastic size (mm)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_Yield_PlasSize
ggplot() +
  geom_point(data = subset(Yield_data, PlasSize < 6), aes(x = PlasSize, y = yi), size = 2, alpha = 0.1, color = "#65C2AD") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_line(data = subset(pred_Yield_data_PlasSize, PlasSize < 2), linetype = "twodash",
            aes(x = PlasSize, y = pred), color = "#65C2AD", linewidth = 0.8) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = NULL, y = NULL) + theme_classic() + theme(
    axis.line = element_blank(), panel.border = element_rect(),
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 9),
    plot.margin = margin(0, 0, 0, 0)) -> p_inset_Yield_PlasSize
p_Sub_Yield_PlasSize + inset_element(
  p_inset_Yield_PlasSize, left = 0.6, right = 0.98,
  bottom = 0.01, top = 0.43) -> p_Sub_Yield_PlasSize2; p_Sub_Yield_PlasSize2
cairo_pdf("Output/Sub-Yield-PlasSize20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_Yield_PlasSize
dev.off()
cairo_pdf("Output/Sub-Yield-PlasSize2 20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_Yield_PlasSize2
dev.off()

pred_Yield_data_MAP <- with(Yield_data, data.frame(
  MAP = seq(min(MAP, na.rm = TRUE), max(MAP, na.rm = TRUE), length.out = 400)))
pred_Yield_data_MAP <- cbind(pred_Yield_data_MAP, predict(model_Yield_MAP, newmods = pred_Yield_data_MAP$MAP))
r2_val_Yield_MAP <- r2_ml(model_Yield_MAP)[1]
p_val_Yield_MAP <- model_Yield_MAP$pval[2]
r2_text_Yield_MAP <- ifelse(r2_val_Yield_MAP < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_Yield_MAP))
p_formatted_Yield_MAP <- ifelse(p_val_Yield_MAP < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_Yield_MAP))

ggplot() +
  geom_point(data = Yield_data, aes(x = MAP, y = yi), size = 3, alpha = 0.1, color = "#65C2AD") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_Yield_data_MAP, fill = "#65C2AD", alpha = 0.25,
              aes(x = MAP, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_Yield_data_MAP, linetype = "solid",
            aes(x = MAP, y = pred), color = "#65C2AD", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_Yield_MAP, p_formatted_Yield_MAP, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "MAP (mm)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_Yield_MAP; p_Sub_Yield_MAP
cairo_pdf("Output/Sub-Yield-MAP20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_Yield_MAP
dev.off()

pred_Yield_data_AridityIndex <- with(Yield_data, data.frame(
  AridityIndex = seq(min(AridityIndex, na.rm = TRUE), max(AridityIndex, na.rm = TRUE), length.out = 400)))
pred_Yield_data_AridityIndex <- cbind(pred_Yield_data_AridityIndex, predict(model_Yield_AridityIndex, newmods = pred_Yield_data_AridityIndex$AridityIndex))
r2_val_Yield_AridityIndex <- r2_ml(model_Yield_AridityIndex)[1]
p_val_Yield_AridityIndex <- model_Yield_AridityIndex$pval[2]
r2_text_Yield_AridityIndex <- ifelse(r2_val_Yield_AridityIndex < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_Yield_AridityIndex))
p_formatted_Yield_AridityIndex <- ifelse(p_val_Yield_AridityIndex < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_Yield_AridityIndex))

ggplot() +
  geom_point(data = Yield_data, aes(x = AridityIndex, y = yi), size = 3, alpha = 0.1, color = "#65C2AD") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_Yield_data_AridityIndex, fill = "#65C2AD", alpha = 0.25,
              aes(x = AridityIndex, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_Yield_data_AridityIndex, linetype = "solid",
            aes(x = AridityIndex, y = pred), color = "#65C2AD", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_Yield_AridityIndex, p_formatted_Yield_AridityIndex, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "Aridity index", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_Yield_AridityIndex; p_Sub_Yield_AridityIndex
cairo_pdf("Output/Sub-Yield-AridityIndex20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_Yield_AridityIndex
dev.off()

sapply(data_SQI[col_names], function(x) length(unique(x)))
V_SQI <- calc.vbd(data = data_SQI, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)

plan(multisession, workers = availableCores() - 1)
models_list_SQI_Sub <- list(
  PlasType2 = list(formula = "~ PlasType2 - 1", mods = ~ PlasType2 - 1),
  PlasSize2 = list(formula = "~ PlasSize2 - 1", mods = ~ PlasSize2 - 1),
  PlasDose = list(formula = "~ PlasDose", mods = ~ PlasDose),
  SampleTime2 = list(formula = "~ SampleTime2", mods = ~ SampleTime2),
  MAT = list(formula = "~ MAT", mods = ~ MAT),
  PlantAppear = list(formula = "~ PlantAppear - 1", mods = ~ PlantAppear - 1))
run_model_SQI_Sub <- function(model_name, model_info) {
  model <- try(rma.mv(
    yi, V_SQI, data = data_SQI,
    random = list(~ 1 | ExpID, ~ 1 | VariableNew/SameVarIDN),
    mods = model_info$mods, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000)))
  return(list(name = model_name,
              formula = model_info$formula,
              model = model))
}
start_SQI_Sub <- Sys.time()
results_SQI_Sub <- future_lapply(
  names(models_list_SQI_Sub),
  function(model_name) {run_model_SQI_Sub(model_name, models_list_SQI_Sub[[model_name]])},
  future.seed = TRUE)
model_results_SQI_Sub <- list()
for (result in results_SQI_Sub) {
  model_name <- paste0("model_SQI_", result$name)
  model_results_SQI_Sub[[model_name]] <- result$model
}
end_SQI_Sub <- Sys.time()
end_SQI_Sub - start_SQI_Sub
plan(sequential)

model_SQI_PlasType2 <- model_results_SQI_Sub$model_SQI_PlasType2
model_SQI_PlasSize2 <- model_results_SQI_Sub$model_SQI_PlasSize2
model_SQI_PlasDose <- model_results_SQI_Sub$model_SQI_PlasDose
model_SQI_SampleTime2 <- model_results_SQI_Sub$model_SQI_SampleTime2
model_SQI_MAT <- model_results_SQI_Sub$model_SQI_MAT
model_SQI_PlantAppear <- model_results_SQI_Sub$model_SQI_PlantAppear

plan(multisession, workers = availableCores() - 1)
models_list_SQI_Sub2 <- list(
  PlasShape = list(formula = "~ PlasShape - 1", mods = ~ PlasShape - 1),
  PlasType = list(formula = "~ PlasType - 1", mods = ~ PlasType - 1),
  PlasSize = list(formula = "~ PlasSize", mods = ~ PlasSize),
  MAP = list(formula = "~ MAP", mods = ~ MAP),
  AridityIndex = list(formula = "~ AridityIndex", mods = ~ AridityIndex),
  PlasDose2 = list(formula = "~ PlasDose2 - 1", mods = ~ PlasDose2 - 1))
run_model_SQI_Sub <- function(model_name, model_info) {
  model <- try(rma.mv(
    yi, V_SQI, data = data_SQI,
    random = list(~ 1 | ExpID, ~ 1 | VariableNew/SameVarIDN),
    mods = model_info$mods, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000)))
  return(list(name = model_name,
              formula = model_info$formula,
              model = model))
}
start_SQI_Sub2 <- Sys.time()
results_SQI_Sub2 <- future_lapply(
  names(models_list_SQI_Sub2),
  function(model_name) {run_model_SQI_Sub(model_name, models_list_SQI_Sub2[[model_name]])},
  future.seed = TRUE)
model_results_SQI_Sub2 <- list()
for (result in results_SQI_Sub2) {
  model_name <- paste0("model_SQI_", result$name)
  model_results_SQI_Sub2[[model_name]] <- result$model
}
end_SQI_Sub2 <- Sys.time()
end_SQI_Sub2 - start_SQI_Sub2
plan(sequential)

model_SQI_PlasShape <- model_results_SQI_Sub2$model_SQI_PlasShape
model_SQI_PlasType <- model_results_SQI_Sub2$model_SQI_PlasType
model_SQI_PlasDose2 <- model_results_SQI_Sub2$model_SQI_PlasDose2
model_SQI_PlasSize <- model_results_SQI_Sub2$model_SQI_PlasSize
model_SQI_MAP <- model_results_SQI_Sub2$model_SQI_MAP
model_SQI_AridityIndex <- model_results_SQI_Sub2$model_SQI_AridityIndex

extract_model_results <- function(model, group_names) {
  summ <- summary(model)
  n_values <- if (!is.null(model$X)) {
    colSums(model$X)} else {rep(NA, length(group_names))}
  data.frame(
    Group = group_names, estimate = as.numeric(summ$b),
    ci.lb = summ$ci.lb, ci.ub = summ$ci.ub, n = n_values,
    contains_zero = (summ$ci.lb <= 0 & summ$ci.ub >= 0),
    QMp = summ$QMp)
}

summary(model_SQI_PlasType2)
summary(model_SQI_PlasSize2)
summary(model_SQI_PlasDose2)
summary(model_SQI_PlasShape)
summary(model_SQI_PlantAppear)
summary(model_SQI_PlasType)

unique(data_combined_with_con7$PlasType2)
unique(data_combined_with_con7$PlasSize2)
unique(data_combined_with_con7$PlasDose2)
unique(data_combined_with_con7$PlasShape)
unique(data_combined_with_con7$PlantAppear)
unique(data_combined_with_con7$PlasType)

df_SQI_PlasType2 <- extract_model_results(
  model_SQI_PlasType2, c("Biodegradable", "Mixture", "NonBiodegradable"))
df_SQI_PlasSize2 <- extract_model_results(
  model_SQI_PlasSize2, c("<0.1", ">5", "0.1-1", "1-5", "unknown"))
df_SQI_PlasDose2 <- extract_model_results(
  model_SQI_PlasDose2, c("<0.02", ">1", "0.02-0.2", "0.2-1", "unknown"))
df_SQI_PlasShape <- extract_model_results(
  model_SQI_PlasShape, c("fiber", "fragment", "granule", "unknown"))
df_SQI_PlantAppear <- extract_model_results(
  model_SQI_PlantAppear, c("noplant", "withplant"))
df_SQI_PlasType <- extract_model_results(
  model_SQI_PlasType,
  c("BioMix", "Mixture", "PBAT", "PBS", "PE",
    "PES", "PHA", "PLA", "POM", "PP", "PS", "PVC"))

result_Sub_SQI <- rbind(
  df_SQI_PlasType2[order(factor(
    df_SQI_PlasType2$Group, levels = c("NonBiodegradable", "Biodegradable", "Mixture"))), ],
  df_SQI_PlasSize2[order(factor(
    df_SQI_PlasSize2$Group, levels = c("<0.1", "0.1-1", "1-5", ">5", "unknown"))), ],
  df_SQI_PlasDose2[order(factor(
    df_SQI_PlasDose2$Group, levels = c("<0.02", "0.02-0.2", "0.2-1", ">1", "unknown"))), ],
  df_SQI_PlasShape[order(factor(
    df_SQI_PlasShape$Group, levels = c("granule", "fiber", "fragment", "unknown"))), ],
  df_SQI_PlantAppear[order(factor(
    df_SQI_PlantAppear$Group, levels = c("noplant", "withplant"))), ]
); rownames(result_Sub_SQI) <- NULL
result_Sub_SQI

blank_row <- function(group) {
  data.frame(Group = group, estimate = NA, ci.lb = NA, ci.ub = NA, n = NA, contains_zero = NA, QMp = NA)}
result_Sub_SQI <- rbind(
  blank_row("PlasType"), result_Sub_SQI[1:3, ],
  blank_row("PlasSize"), result_Sub_SQI[4:8, ],
  blank_row("PlasDose"), result_Sub_SQI[9:13, ],
  blank_row("PlasShape"), result_Sub_SQI[14:17, ],
  blank_row("Plant"), result_Sub_SQI[18:19, ]
); rownames(result_Sub_SQI) <- NULL
result_Sub_SQI

result_Sub_SQI <- result_Sub_SQI %>%
  mutate(Group1 = factor(
    c("PlasType", "NonBio", "Bio", "Mixture",
      "PlasSize", "≤0.1", "0.1-1", "1-5", ">5", "unknown",
      "PlasDose", "≤0.02", "0.02-0.2", "0.2-1", ">1", "unknownD",
      "PlasShape", "Granule", "Fiber", "Fragment", "unknownS",
      "Plant", "NoPlant", "WithPlant")))

result_Sub_SQI$Group1 <- factor(result_Sub_SQI$Group1, levels = rev(unique(result_Sub_SQI$Group1)))

summ <- summary(model_SQI)
overall_est <- summ$b[1]
overall_ci.lb <- summ$ci.lb[1]
overall_ci.ub <- summ$ci.ub[1]
total_n <- nrow(data_SQI)
overall_contains_zero <- (overall_ci.lb <= 0 & overall_ci.ub >= 0)
overall_row <- data.frame(
  Group = "Overall", estimate = overall_est,
  ci.lb = overall_ci.lb, ci.ub = overall_ci.ub, n = total_n,
  contains_zero = overall_contains_zero, QMp = NA_real_,
  Group1 = "Overall", stringsAsFactors = FALSE)

orig_levels <- as.character(unique(result_Sub_SQI$Group1))
result_Sub_SQI <- rbind(overall_row, result_Sub_SQI)
result_Sub_SQI$Group1 <- factor(result_Sub_SQI$Group1, levels = c(rev(orig_levels), "Overall"))

x_max_SQI <- with(result_Sub_SQI[!is.na(result_Sub_SQI$estimate), ], max(abs(c(ci.lb, ci.ub)), na.rm = TRUE))
x_limit_SQI <- x_max_SQI * 1.4

ggplot(result_Sub_SQI, aes(x = estimate, y = Group1)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3.5, ymax = 8.5, fill = "grey98") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 14.5, ymax = 20.5, fill = "grey98") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 24.5, ymax = Inf, fill = "grey98") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), color = "#BBB74E", width = 0.15, linewidth = 1) +
  geom_point(aes(fill = contains_zero), shape = 21, size = 3, color = "#BBB74E", stroke = 1.5) +
  scale_fill_manual(values = c("FALSE" = "#BBB74E", "TRUE" = "white"), guide = "none", drop = FALSE) +
  geom_text(aes(label = ifelse(!is.na(n), paste0("(", n, ")"), "")),
            x = x_limit_SQI * 0.995, hjust = 1, size = 12/2.8346, color = "black") +
  geom_text(data = subset(result_Sub_SQI, Group1 %in% c("NonBio", "≤0.1", "≤0.02", "Granule", "NoPlant")),
            aes(label = ifelse(QMp < 0.001, "italic('P') ~ '<' ~ 0.001", sprintf("italic('P') ~ '= %.3f'", QMp))),
            x = -x_limit_SQI * 0.995, nudge_y = 1, hjust = 0, vjust = 0.5, size = 13/2.8346, color = "black", parse = TRUE) +
  scale_y_discrete(labels = c(
    "PlasType" = "Type", "PlasSize" = "Size", "PlasDose" = "Dose",
    "PlasShape" = "Shape", "unknownD" = "unknown", "unknownS" = "unknown",
    "WithPlant" = "Presence", "NoPlant" = "Absence")) +
  labs(x = "Effect Size", y = NULL) + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    axis.text.y = element_text(size = 15, face = ifelse(levels(result_Sub_SQI$Group1) %in% c(
      "Overall", "PlasType", "PlasSize", "PlasDose", "PlasShape", "Plant"), "bold", "plain")),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 5/2.5) +
  coord_cartesian(xlim = c(-x_limit_SQI, x_limit_SQI)) -> p_Sub_SQI; p_Sub_SQI
cairo_pdf("Output/Sub-SQI20260623.pdf", bg = "transparent", width = 10.85/2.54, height = 30/2.54)
p_Sub_SQI
dev.off()

data_SQI$ci.lb <- data_SQI$yi - 1.96 * sqrt(data_SQI$vi)
data_SQI$ci.ub <- data_SQI$yi + 1.96 * sqrt(data_SQI$vi)
data_SQI <- data_SQI[order(data_SQI$yi), ]
data_SQI$StudyID <- factor(1:nrow(data_SQI), levels = 1:nrow(data_SQI))

x_max_SQI <- max(abs(c(data_SQI$ci.lb, data_SQI$ci.ub)), na.rm = TRUE)
x_limit_SQI <- x_max_SQI * 1.4

ggplot(data_SQI, aes(x = yi, y = StudyID)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), width = 0, color = "#BBB74E", alpha = 0.3, linewidth = 0.2) +
  geom_point(shape = 16, size = 0.7, color = "#BBB74E") +
  labs(x = "Effect Size", y = "Observations") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1.1/2.5) +
  coord_cartesian(xlim = c(-x_limit_SQI, x_limit_SQI)) -> p_Indi_SQI; p_Indi_SQI
cairo_pdf("Output/Indi-SQI20260623.pdf", bg = "transparent", width = 8.88/2.54, height = 30/2.54)
p_Indi_SQI
dev.off()

result_Sub_SQI_PlasType <- df_SQI_PlasType %>% arrange(factor(
  Group, levels = c("PE", "PP", "PS", "PVC", "PES", "PAN", "POM", "TP",
                    "PLA", "PHA", "PBS", "PBAT", "PCL", "BioMix", "Mixture"))) %>% `rownames<-`(NULL)
result_Sub_SQI_PlasType$Group
result_Sub_SQI_PlasType <- result_Sub_SQI_PlasType %>% mutate(Group1 = factor(
  c("PE", "PP", "PS", "PVC", "PES", "POM",
    "PLA", "PHA", "PBS", "PBAT", "BioMix", "Mixture")))

result_Sub_SQI_PlasType$Group1 <- with(result_Sub_SQI_PlasType, factor(Group1, levels = rev(unique(Group1))))
x_max_SQI_PlasType <- with(result_Sub_SQI_PlasType[!is.na(result_Sub_SQI_PlasType$estimate), ], max(abs(c(ci.lb, ci.ub)), na.rm = TRUE))
x_limit_SQI_PlasType <- x_max_SQI_PlasType * 1.4

ggplot(result_Sub_SQI_PlasType, aes(x = estimate, y = Group1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), color = "#BBB74E", width = 0.15, linewidth = 1) +
  geom_point(aes(fill = contains_zero), shape = 21, size = 3, color = "#BBB74E", stroke = 1.5) +
  scale_fill_manual(values = c("FALSE" = "#BBB74E", "TRUE" = "white"), guide = "none", drop = FALSE) +
  geom_text(aes(label = ifelse(!is.na(n), paste0("(", n, ")"), "")),
            x = x_max_SQI_PlasType * 1.4*0.995, hjust = 1, size = 13/2.8346, color = "black") +
  coord_cartesian(xlim = c(-x_limit_SQI_PlasType, x_limit_SQI_PlasType)) +
  labs(x = "Effect Size", y = NULL) + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 3/2) -> p_Sub_SQI_PlasType; p_Sub_SQI_PlasType
cairo_pdf("Output/Sub-SQI-PlasType20260623.pdf", bg = "transparent", width = 30/2.54, height = 12/2.54)
p_Sub_SQI_PlasType
dev.off()

summary(model_SQI_PlasSize)
summary(model_SQI_PlasDose)
summary(model_SQI_SampleTime2)
summary(model_SQI_MAT)
summary(model_SQI_MAP)
summary(model_SQI_AridityIndex)

unique(data_combined_with_con7$PlasSize)
unique(data_combined_with_con7$PlasDose)
unique(data_combined_with_con7$SampleTime2)
unique(data_combined_with_con7$MAT)
unique(data_combined_with_con7$MAP)
unique(data_combined_with_con7$AridityIndex)

pred_data_SQI_PlasDose <- with(data_SQI, data.frame(
  PlasDose = seq(min(PlasDose, na.rm = TRUE), max(PlasDose, na.rm = TRUE), length.out = 400)))
pred_data_SQI_PlasDose <- cbind(pred_data_SQI_PlasDose, predict(model_SQI_PlasDose, newmods = pred_data_SQI_PlasDose$PlasDose))
r2_val_SQI_PlasDose <- r2_ml(model_SQI_PlasDose)[1]
p_val_SQI_PlasDose <- model_SQI_PlasDose$pval[2]
r2_text_SQI_PlasDose <- ifelse(r2_val_SQI_PlasDose < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_SQI_PlasDose))
p_formatted_SQI_PlasDose <- ifelse(p_val_SQI_PlasDose < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_SQI_PlasDose))

ggplot() +
  geom_point(data = data_SQI, aes(x = PlasDose, y = yi), size = 3, alpha = 0.1, color = "#BBB74E") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_data_SQI_PlasDose, fill = "#BBB74E", alpha = 0.25,
              aes(x = PlasDose, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_data_SQI_PlasDose, linetype = "solid",
            aes(x = PlasDose, y = pred), color = "#BBB74E", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_SQI_PlasDose, p_formatted_SQI_PlasDose, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "Plastic dose (%)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_SQI_PlasDose; p_Sub_SQI_PlasDose
cairo_pdf("Output/Sub-SQI-PlasDose20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_SQI_PlasDose
dev.off()

pred_data_SQI_SampleTime2 <- with(data_SQI, data.frame(
  SampleTime2 = seq(min(SampleTime2, na.rm = TRUE), max(SampleTime2, na.rm = TRUE), length.out = 400)))
pred_data_SQI_SampleTime2 <- cbind(pred_data_SQI_SampleTime2, predict(model_SQI_SampleTime2, newmods = pred_data_SQI_SampleTime2$SampleTime2))
r2_val_SQI_SampleTime2 <- r2_ml(model_SQI_SampleTime2)[1]
p_val_SQI_SampleTime2 <- model_SQI_SampleTime2$pval[2]
r2_text_SQI_SampleTime2 <- ifelse(r2_val_SQI_SampleTime2 < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_SQI_SampleTime2))
p_formatted_SQI_SampleTime2 <- ifelse(p_val_SQI_SampleTime2 < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_SQI_SampleTime2))

ggplot() +
  geom_point(data = data_SQI, aes(x = SampleTime2, y = yi), size = 3, alpha = 0.1, color = "#BBB74E") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_data_SQI_SampleTime2, fill = "#BBB74E", alpha = 0.25,
              aes(x = SampleTime2, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_data_SQI_SampleTime2, linetype = "solid",
            aes(x = SampleTime2, y = pred), color = "#BBB74E", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_SQI_SampleTime2, p_formatted_SQI_SampleTime2, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "Duration (month)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_SQI_SampleTime2; p_Sub_SQI_SampleTime2
cairo_pdf("Output/Sub-SQI-SampleTime2 20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_SQI_SampleTime2
dev.off()

pred_data_SQI_MAT <- with(data_SQI, data.frame(
  MAT = seq(min(MAT, na.rm = TRUE), max(MAT, na.rm = TRUE), length.out = 400)))
pred_data_SQI_MAT <- cbind(pred_data_SQI_MAT, predict(model_SQI_MAT, newmods = pred_data_SQI_MAT$MAT))
r2_val_SQI_MAT <- r2_ml(model_SQI_MAT)[1]
p_val_SQI_MAT <- model_SQI_MAT$pval[2]
r2_text_SQI_MAT <- ifelse(r2_val_SQI_MAT < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_SQI_MAT))
p_formatted_SQI_MAT <- ifelse(p_val_SQI_MAT < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_SQI_MAT))

ggplot() +
  geom_point(data = data_SQI, aes(x = MAT, y = yi), size = 3, alpha = 0.1, color = "#BBB74E") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_data_SQI_MAT, fill = "#BBB74E", alpha = 0.25,
              aes(x = MAT, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_data_SQI_MAT, linetype = "solid",
            aes(x = MAT, y = pred), color = "#BBB74E", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_SQI_MAT, p_formatted_SQI_MAT, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "MAT (°C)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_SQI_MAT; p_Sub_SQI_MAT
cairo_pdf("Output/Sub-SQI-MAT20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_SQI_MAT
dev.off()

pred_data_SQI_PlasSize <- with(data_SQI, data.frame(
  PlasSize = seq(min(PlasSize, na.rm = TRUE), max(PlasSize, na.rm = TRUE), length.out = 400)))
pred_data_SQI_PlasSize <- cbind(pred_data_SQI_PlasSize, predict(model_SQI_PlasSize, newmods = pred_data_SQI_PlasSize$PlasSize))
r2_val_SQI_PlasSize <- r2_ml(model_SQI_PlasSize)[1]
p_val_SQI_PlasSize <- model_SQI_PlasSize$pval[2]
r2_text_SQI_PlasSize <- ifelse(r2_val_SQI_PlasSize < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_SQI_PlasSize))
p_formatted_SQI_PlasSize <- ifelse(p_val_SQI_PlasSize < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_SQI_PlasSize))

ggplot() +
  geom_point(data = data_SQI, aes(x = PlasSize, y = yi), size = 3, alpha = 0.1, color = "#BBB74E") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_data_SQI_PlasSize, fill = "#BBB74E", alpha = 0.25,
              aes(x = PlasSize, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_data_SQI_PlasSize, linetype = "solid",
            aes(x = PlasSize, y = pred), color = "#BBB74E", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_SQI_PlasSize, p_formatted_SQI_PlasSize, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "Plastic size (mm)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_SQI_PlasSize
ggplot() +
  geom_point(data = subset(data_SQI, PlasSize < 6), aes(x = PlasSize, y = yi), size = 2, alpha = 0.1, color = "#BBB74E") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = subset(pred_data_SQI_PlasSize, PlasSize < 6), fill = "#BBB74E", alpha = 0.25,
              aes(x = PlasSize, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = subset(pred_data_SQI_PlasSize, PlasSize < 6), linetype = "solid",
            aes(x = PlasSize, y = pred), color = "#BBB74E", linewidth = 0.8) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = NULL, y = NULL) + theme_classic() + theme(
    axis.line = element_blank(), panel.border = element_rect(),
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 9),
    plot.margin = margin(0, 0, 0, 0)) -> p_inset_SQI_PlasSize
p_Sub_SQI_PlasSize + inset_element(
  p_inset_SQI_PlasSize, left = 0.45, right = 0.98,
  bottom = 0.01, top = 0.43) -> p_Sub_SQI_PlasSize2; p_Sub_SQI_PlasSize2
cairo_pdf("Output/Sub-SQI-PlasSize20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_SQI_PlasSize
dev.off()
cairo_pdf("Output/Sub-SQI-PlasSize2 20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_SQI_PlasSize2
dev.off()

pred_data_SQI_MAP <- with(data_SQI, data.frame(
  MAP = seq(min(MAP, na.rm = TRUE), max(MAP, na.rm = TRUE), length.out = 400)))
pred_data_SQI_MAP <- cbind(pred_data_SQI_MAP, predict(model_SQI_MAP, newmods = pred_data_SQI_MAP$MAP))
r2_val_SQI_MAP <- r2_ml(model_SQI_MAP)[1]
p_val_SQI_MAP <- model_SQI_MAP$pval[2]
r2_text_SQI_MAP <- ifelse(r2_val_SQI_MAP < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_SQI_MAP))
p_formatted_SQI_MAP <- ifelse(p_val_SQI_MAP < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_SQI_MAP))

ggplot() +
  geom_point(data = data_SQI, aes(x = MAP, y = yi), size = 3, alpha = 0.1, color = "#BBB74E") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_data_SQI_MAP, fill = "#BBB74E", alpha = 0.25,
              aes(x = MAP, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_data_SQI_MAP, linetype = "solid",
            aes(x = MAP, y = pred), color = "#BBB74E", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_SQI_MAP, p_formatted_SQI_MAP, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "MAP (mm)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_SQI_MAP; p_Sub_SQI_MAP
cairo_pdf("Output/Sub-SQI-MAP20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_SQI_MAP
dev.off()

pred_data_SQI_AridityIndex <- with(data_SQI, data.frame(
  AridityIndex = seq(min(AridityIndex, na.rm = TRUE), max(AridityIndex, na.rm = TRUE), length.out = 400)))
pred_data_SQI_AridityIndex <- cbind(pred_data_SQI_AridityIndex, predict(model_SQI_AridityIndex, newmods = pred_data_SQI_AridityIndex$AridityIndex))
r2_val_SQI_AridityIndex <- r2_ml(model_SQI_AridityIndex)[1]
p_val_SQI_AridityIndex <- model_SQI_AridityIndex$pval[2]
r2_text_SQI_AridityIndex <- ifelse(r2_val_SQI_AridityIndex < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_SQI_AridityIndex))
p_formatted_SQI_AridityIndex <- ifelse(p_val_SQI_AridityIndex < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_SQI_AridityIndex))

ggplot() +
  geom_point(data = data_SQI, aes(x = AridityIndex, y = yi), size = 3, alpha = 0.1, color = "#BBB74E") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_data_SQI_AridityIndex, fill = "#BBB74E", alpha = 0.25,
              aes(x = AridityIndex, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_data_SQI_AridityIndex, linetype = "twodash",
            aes(x = AridityIndex, y = pred), color = "#BBB74E", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_SQI_AridityIndex, p_formatted_SQI_AridityIndex, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "Aridity index", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_SQI_AridityIndex; p_Sub_SQI_AridityIndex
cairo_pdf("Output/Sub-SQI-AridityIndex20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_SQI_AridityIndex
dev.off()

sapply(data_EMF[col_names], function(x) length(unique(x)))
V_EMF <- calc.vbd(data = data_EMF, CommonID = CommonIDN, m1 = ConMean, sd1 = ConSD_new, n1 = ConN, vi = vi)

plan(multisession, workers = availableCores() - 1)
models_list_EMF_Sub <- list(
  PlasType2 = list(formula = "~ PlasType2 - 1", mods = ~ PlasType2 - 1),
  PlasSize2 = list(formula = "~ PlasSize2 - 1", mods = ~ PlasSize2 - 1),
  PlasDose = list(formula = "~ PlasDose", mods = ~ PlasDose),
  SampleTime2 = list(formula = "~ SampleTime2", mods = ~ SampleTime2),
  MAT = list(formula = "~ MAT", mods = ~ MAT),
  PlantAppear = list(formula = "~ PlantAppear - 1", mods = ~ PlantAppear - 1))
run_model_EMF_Sub <- function(model_name, model_info) {
  model <- try(rma.mv(
    yi, V_EMF, data = data_EMF,
    random = list(~ 1 | ExpID, ~ 1 | EMFCat/EMFFun/VariableNew/SameVarIDN),
    mods = model_info$mods, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000)))
  return(list(name = model_name,
              formula = model_info$formula,
              model = model))
}
start_EMF_Sub <- Sys.time()
results_EMF_Sub <- future_lapply(
  names(models_list_EMF_Sub),
  function(model_name) {run_model_EMF_Sub(model_name, models_list_EMF_Sub[[model_name]])},
  future.seed = TRUE)
model_results_EMF_Sub <- list()
for (result in results_EMF_Sub) {
  model_name <- paste0("model_EMF_", result$name)
  model_results_EMF_Sub[[model_name]] <- result$model
}
end_EMF_Sub <- Sys.time()
end_EMF_Sub - start_EMF_Sub
plan(sequential)

model_EMF_PlasType2 <- model_results_EMF_Sub$model_EMF_PlasType2
model_EMF_PlasSize2 <- model_results_EMF_Sub$model_EMF_PlasSize2
model_EMF_PlasDose <- model_results_EMF_Sub$model_EMF_PlasDose
model_EMF_SampleTime2 <- model_results_EMF_Sub$model_EMF_SampleTime2
model_EMF_MAT <- model_results_EMF_Sub$model_EMF_MAT
model_EMF_PlantAppear <- model_results_EMF_Sub$model_EMF_PlantAppear

plan(multisession, workers = availableCores() - 1)
models_list_EMF_Sub2 <- list(
  PlasShape = list(formula = "~ PlasShape - 1", mods = ~ PlasShape - 1),
  PlasType = list(formula = "~ PlasType - 1", mods = ~ PlasType - 1),
  PlasSize = list(formula = "~ PlasSize", mods = ~ PlasSize),
  MAP = list(formula = "~ MAP", mods = ~ MAP),
  AridityIndex = list(formula = "~ AridityIndex", mods = ~ AridityIndex),
  PlasDose2 = list(formula = "~ PlasDose2 - 1", mods = ~ PlasDose2 - 1))
run_model_EMF_Sub <- function(model_name, model_info) {
  model <- try(rma.mv(
    yi, V_EMF, data = data_EMF,
    random = list(~ 1 | ExpID, ~ 1 | EMFCat/EMFFun/VariableNew/SameVarIDN),
    mods = model_info$mods, test = "t", dfs = "contain", control = list(optimizer = "Nelder-Mead", maxit = 10000)))
  return(list(name = model_name,
              formula = model_info$formula,
              model = model))
}
start_EMF_Sub2 <- Sys.time()
results_EMF_Sub2 <- future_lapply(
  names(models_list_EMF_Sub2),
  function(model_name) {run_model_EMF_Sub(model_name, models_list_EMF_Sub2[[model_name]])},
  future.seed = TRUE)
model_results_EMF_Sub2 <- list()
for (result in results_EMF_Sub2) {
  model_name <- paste0("model_EMF_", result$name)
  model_results_EMF_Sub2[[model_name]] <- result$model
}
end_EMF_Sub2 <- Sys.time()
end_EMF_Sub2 - start_EMF_Sub2
plan(sequential)

model_EMF_PlasShape <- model_results_EMF_Sub2$model_EMF_PlasShape
model_EMF_PlasType <- model_results_EMF_Sub2$model_EMF_PlasType
model_EMF_PlasDose2 <- model_results_EMF_Sub2$model_EMF_PlasDose2
model_EMF_PlasSize <- model_results_EMF_Sub2$model_EMF_PlasSize
model_EMF_MAP <- model_results_EMF_Sub2$model_EMF_MAP
model_EMF_AridityIndex <- model_results_EMF_Sub2$model_EMF_AridityIndex

extract_model_results <- function(model, group_names) {
  summ <- summary(model)
  n_values <- if (!is.null(model$X)) {
    colSums(model$X)} else {rep(NA, length(group_names))}
  data.frame(
    Group = group_names, estimate = as.numeric(summ$b),
    ci.lb = summ$ci.lb, ci.ub = summ$ci.ub, n = n_values,
    contains_zero = (summ$ci.lb <= 0 & summ$ci.ub >= 0),
    QMp = summ$QMp)
}

summary(model_EMF_PlasType2)
summary(model_EMF_PlasSize2)
summary(model_EMF_PlasDose2)
summary(model_EMF_PlasShape)
summary(model_EMF_PlantAppear)
summary(model_EMF_PlasType)

unique(data_combined_with_con7$PlasType2)
unique(data_combined_with_con7$PlasSize2)
unique(data_combined_with_con7$PlasDose2)
unique(data_combined_with_con7$PlasShape)
unique(data_combined_with_con7$PlantAppear)
unique(data_combined_with_con7$PlasType)

df_EMF_PlasType2 <- extract_model_results(
  model_EMF_PlasType2, c("Biodegradable", "Mixture", "NonBiodegradable"))
df_EMF_PlasSize2 <- extract_model_results(
  model_EMF_PlasSize2, c("<0.1", ">5", "0.1-1", "1-5", "unknown"))
df_EMF_PlasDose2 <- extract_model_results(
  model_EMF_PlasDose2, c("<0.02", ">1", "0.02-0.2", "0.2-1", "unknown"))
df_EMF_PlasShape <- extract_model_results(
  model_EMF_PlasShape, c("fiber", "fragment", "granule", "unknown"))
df_EMF_PlantAppear <- extract_model_results(
  model_EMF_PlantAppear, c("noplant", "withplant"))
df_EMF_PlasType <- extract_model_results(
  model_EMF_PlasType,
  c("BioMix", "Mixture", "PBAT", "PBS", "PE",
    "PES", "PHA", "PLA", "POM", "PP", "PS", "PVC"))

result_Sub_EMF <- rbind(
  df_EMF_PlasType2[order(factor(
    df_EMF_PlasType2$Group, levels = c("NonBiodegradable", "Biodegradable", "Mixture"))), ],
  df_EMF_PlasSize2[order(factor(
    df_EMF_PlasSize2$Group, levels = c("<0.1", "0.1-1", "1-5", ">5", "unknown"))), ],
  df_EMF_PlasDose2[order(factor(
    df_EMF_PlasDose2$Group, levels = c("<0.02", "0.02-0.2", "0.2-1", ">1", "unknown"))), ],
  df_EMF_PlasShape[order(factor(
    df_EMF_PlasShape$Group, levels = c("granule", "fiber", "fragment", "unknown"))), ],
  df_EMF_PlantAppear[order(factor(
    df_EMF_PlantAppear$Group, levels = c("noplant", "withplant"))), ]
); rownames(result_Sub_EMF) <- NULL
result_Sub_EMF

blank_row <- function(group) {
  data.frame(Group = group, estimate = NA, ci.lb = NA, ci.ub = NA, n = NA, contains_zero = NA, QMp = NA)}
result_Sub_EMF <- rbind(
  blank_row("PlasType"), result_Sub_EMF[1:3, ],
  blank_row("PlasSize"), result_Sub_EMF[4:8, ],
  blank_row("PlasDose"), result_Sub_EMF[9:13, ],
  blank_row("PlasShape"), result_Sub_EMF[14:17, ],
  blank_row("Plant"), result_Sub_EMF[18:19, ]
); rownames(result_Sub_EMF) <- NULL
result_Sub_EMF

result_Sub_EMF <- result_Sub_EMF %>%
  mutate(Group1 = factor(
    c("PlasType", "NonBio", "Bio", "Mixture",
      "PlasSize", "≤0.1", "0.1-1", "1-5", ">5", "unknown",
      "PlasDose", "≤0.02", "0.02-0.2", "0.2-1", ">1", "unknownD",
      "PlasShape", "Granule", "Fiber", "Fragment", "unknownS",
      "Plant", "NoPlant", "WithPlant")))

result_Sub_EMF$Group1 <- factor(result_Sub_EMF$Group1, levels = rev(unique(result_Sub_EMF$Group1)))

summ <- summary(model_EMF)
overall_est <- summ$b[1]
overall_ci.lb <- summ$ci.lb[1]
overall_ci.ub <- summ$ci.ub[1]
total_n <- nrow(data_EMF)
overall_contains_zero <- (overall_ci.lb <= 0 & overall_ci.ub >= 0)
overall_row <- data.frame(
  Group = "Overall", estimate = overall_est,
  ci.lb = overall_ci.lb, ci.ub = overall_ci.ub, n = total_n,
  contains_zero = overall_contains_zero, QMp = NA_real_,
  Group1 = "Overall", stringsAsFactors = FALSE)

orig_levels <- as.character(unique(result_Sub_EMF$Group1))
result_Sub_EMF <- rbind(overall_row, result_Sub_EMF)
result_Sub_EMF$Group1 <- factor(result_Sub_EMF$Group1, levels = c(rev(orig_levels), "Overall"))

x_max_EMF <- with(result_Sub_EMF[!is.na(result_Sub_EMF$estimate), ], max(abs(c(ci.lb, ci.ub)), na.rm = TRUE))
x_limit_EMF <- x_max_EMF * 1.4

ggplot(result_Sub_EMF, aes(x = estimate, y = Group1)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3.5, ymax = 8.5, fill = "grey98") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 14.5, ymax = 20.5, fill = "grey98") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 24.5, ymax = Inf, fill = "grey98") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), color = "#95B7DA", width = 0.15, linewidth = 1) +
  geom_point(aes(fill = contains_zero), shape = 21, size = 3, color = "#95B7DA", stroke = 1.5) +
  scale_fill_manual(values = c("FALSE" = "#95B7DA", "TRUE" = "white"), guide = "none", drop = FALSE) +
  geom_text(aes(label = ifelse(!is.na(n), paste0("(", n, ")"), "")),
            x = x_limit_EMF * 0.995, hjust = 1, size = 12/2.8346, color = "black") +
  geom_text(data = subset(result_Sub_EMF, Group1 %in% c("NonBio", "≤0.1", "≤0.02", "Granule", "NoPlant")),
            aes(label = ifelse(QMp < 0.001, "italic('P') ~ '<' ~ 0.001", sprintf("italic('P') ~ '= %.3f'", QMp))),
            x = -x_limit_EMF * 0.995, nudge_y = 1, hjust = 0, vjust = 0.5, size = 13/2.8346, color = "black", parse = TRUE) +
  scale_y_discrete(labels = c(
    "PlasType" = "Type", "PlasSize" = "Size", "PlasDose" = "Dose",
    "PlasShape" = "Shape", "unknownD" = "unknown", "unknownS" = "unknown",
    "WithPlant" = "Presence", "NoPlant" = "Absence")) +
  labs(x = "Effect Size", y = NULL) + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    axis.text.y = element_text(size = 15, face = ifelse(levels(result_Sub_EMF$Group1) %in% c(
      "Overall", "PlasType", "PlasSize", "PlasDose", "PlasShape", "Plant"), "bold", "plain")),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 5/2.5) +
  coord_cartesian(xlim = c(-x_limit_EMF, x_limit_EMF)) -> p_Sub_EMF; p_Sub_EMF
cairo_pdf("Output/Sub-EMF20260623.pdf", bg = "transparent", width = 10.85/2.54, height = 30/2.54)
p_Sub_EMF
dev.off()

data_EMF$ci.lb <- data_EMF$yi - 1.96 * sqrt(data_EMF$vi)
data_EMF$ci.ub <- data_EMF$yi + 1.96 * sqrt(data_EMF$vi)
data_EMF <- data_EMF[order(data_EMF$yi), ]
data_EMF$StudyID <- factor(1:nrow(data_EMF), levels = 1:nrow(data_EMF))

x_max_EMF <- max(abs(c(data_EMF$ci.lb, data_EMF$ci.ub)), na.rm = TRUE)
x_limit_EMF <- x_max_EMF * 1.4

ggplot(data_EMF, aes(x = yi, y = StudyID)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), width = 0, color = "#95B7DA", alpha = 0.3, linewidth = 0.2) +
  geom_point(shape = 16, size = 0.7, color = "#95B7DA") +
  labs(x = "Effect Size", y = "Observations") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1.1/2.5) +
  coord_cartesian(xlim = c(-x_limit_EMF, x_limit_EMF)) -> p_Indi_EMF; p_Indi_EMF
cairo_pdf("Output/Indi-EMF20260623.pdf", bg = "transparent", width = 8.88/2.54, height = 30/2.54)
p_Indi_EMF
dev.off()

result_Sub_EMF_PlasType <- df_EMF_PlasType %>% arrange(factor(
  Group, levels = c("PE", "PP", "PS", "PVC", "PES", "PAN", "POM", "TP",
                    "PLA", "PHA", "PBS", "PBAT", "PCL", "BioMix", "Mixture"))) %>% `rownames<-`(NULL)
result_Sub_EMF_PlasType$Group
result_Sub_EMF_PlasType <- result_Sub_EMF_PlasType %>% mutate(Group1 = factor(
  c("PE", "PP", "PS", "PVC", "PES", "POM",
    "PLA", "PHA", "PBS", "PBAT", "BioMix", "Mixture")))

result_Sub_EMF_PlasType$Group1 <- with(result_Sub_EMF_PlasType, factor(Group1, levels = rev(unique(Group1))))
x_max_EMF_PlasType <- with(result_Sub_EMF_PlasType[!is.na(result_Sub_EMF_PlasType$estimate), ], max(abs(c(ci.lb, ci.ub)), na.rm = TRUE))
x_limit_EMF_PlasType <- x_max_EMF_PlasType * 1.4

ggplot(result_Sub_EMF_PlasType, aes(x = estimate, y = Group1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), color = "#95B7DA", width = 0.15, linewidth = 1) +
  geom_point(aes(fill = contains_zero), shape = 21, size = 3, color = "#95B7DA", stroke = 1.5) +
  scale_fill_manual(values = c("FALSE" = "#95B7DA", "TRUE" = "white"), guide = "none", drop = FALSE) +
  geom_text(aes(label = ifelse(!is.na(n), paste0("(", n, ")"), "")),
            x = x_max_EMF_PlasType * 1.4*0.995, hjust = 1, size = 12/2.8346, color = "black") +
  coord_cartesian(xlim = c(-x_limit_EMF_PlasType, x_limit_EMF_PlasType)) +
  labs(x = "Effect Size", y = NULL) + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 3/2) -> p_Sub_EMF_PlasType; p_Sub_EMF_PlasType
cairo_pdf("Output/Sub-EMF-PlasType20260623.pdf", bg = "transparent", width = 30/2.54, height = 12/2.54)
p_Sub_EMF_PlasType
dev.off()

summary(model_EMF_PlasSize)
summary(model_EMF_PlasDose)
summary(model_EMF_SampleTime2)
summary(model_EMF_MAT)
summary(model_EMF_MAP)
summary(model_EMF_AridityIndex)

unique(data_combined_with_con7$PlasSize)
unique(data_combined_with_con7$PlasDose)
unique(data_combined_with_con7$SampleTime2)
unique(data_combined_with_con7$MAT)
unique(data_combined_with_con7$MAP)
unique(data_combined_with_con7$AridityIndex)

pred_data_EMF_PlasDose <- with(data_EMF, data.frame(
  PlasDose = seq(min(PlasDose, na.rm = TRUE), max(PlasDose, na.rm = TRUE), length.out = 400)))
pred_data_EMF_PlasDose <- cbind(pred_data_EMF_PlasDose, predict(model_EMF_PlasDose, newmods = pred_data_EMF_PlasDose$PlasDose))
r2_val_EMF_PlasDose <- r2_ml(model_EMF_PlasDose)[1]
p_val_EMF_PlasDose <- model_EMF_PlasDose$pval[2]
r2_text_EMF_PlasDose <- ifelse(r2_val_EMF_PlasDose < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_EMF_PlasDose))
p_formatted_EMF_PlasDose <- ifelse(p_val_EMF_PlasDose < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_EMF_PlasDose))

ggplot() +
  geom_point(data = data_EMF, aes(x = PlasDose, y = yi), size = 3, alpha = 0.1, color = "#95B7DA") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_data_EMF_PlasDose, fill = "#95B7DA", alpha = 0.25,
              aes(x = PlasDose, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_data_EMF_PlasDose, linetype = "solid",
            aes(x = PlasDose, y = pred), color = "#95B7DA", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_EMF_PlasDose, p_formatted_EMF_PlasDose, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "Plastic dose (%)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_EMF_PlasDose; p_Sub_EMF_PlasDose
cairo_pdf("Output/Sub-EMF-PlasDose20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_EMF_PlasDose
dev.off()

pred_data_EMF_SampleTime2 <- with(data_EMF, data.frame(
  SampleTime2 = seq(min(SampleTime2, na.rm = TRUE), max(SampleTime2, na.rm = TRUE), length.out = 400)))
pred_data_EMF_SampleTime2 <- cbind(pred_data_EMF_SampleTime2, predict(model_EMF_SampleTime2, newmods = pred_data_EMF_SampleTime2$SampleTime2))
r2_val_EMF_SampleTime2 <- r2_ml(model_EMF_SampleTime2)[1]
p_val_EMF_SampleTime2 <- model_EMF_SampleTime2$pval[2]
r2_text_EMF_SampleTime2 <- ifelse(r2_val_EMF_SampleTime2 < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_EMF_SampleTime2))
p_formatted_EMF_SampleTime2 <- ifelse(p_val_EMF_SampleTime2 < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_EMF_SampleTime2))

ggplot() +
  geom_point(data = data_EMF, aes(x = SampleTime2, y = yi), size = 3, alpha = 0.1, color = "#95B7DA") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_data_EMF_SampleTime2, fill = "#95B7DA", alpha = 0.25,
              aes(x = SampleTime2, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_data_EMF_SampleTime2, linetype = "solid",
            aes(x = SampleTime2, y = pred), color = "#95B7DA", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_EMF_SampleTime2, p_formatted_EMF_SampleTime2, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "Duration (month)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_EMF_SampleTime2; p_Sub_EMF_SampleTime2
cairo_pdf("Output/Sub-EMF-SampleTime2 20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_EMF_SampleTime2
dev.off()

pred_data_EMF_MAT <- with(data_EMF, data.frame(
  MAT = seq(min(MAT, na.rm = TRUE), max(MAT, na.rm = TRUE), length.out = 400)))
pred_data_EMF_MAT <- cbind(pred_data_EMF_MAT, predict(model_EMF_MAT, newmods = pred_data_EMF_MAT$MAT))
r2_val_EMF_MAT <- r2_ml(model_EMF_MAT)[1]
p_val_EMF_MAT <- model_EMF_MAT$pval[2]
r2_text_EMF_MAT <- ifelse(r2_val_EMF_MAT < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_EMF_MAT))
p_formatted_EMF_MAT <- ifelse(p_val_EMF_MAT < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_EMF_MAT))

ggplot() +
  geom_point(data = data_EMF, aes(x = MAT, y = yi), size = 3, alpha = 0.1, color = "#95B7DA") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_data_EMF_MAT, fill = "#95B7DA", alpha = 0.25,
              aes(x = MAT, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_data_EMF_MAT, linetype = "solid",
            aes(x = MAT, y = pred), color = "#95B7DA", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_EMF_MAT, p_formatted_EMF_MAT, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "MAT (°C)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_EMF_MAT; p_Sub_EMF_MAT
cairo_pdf("Output/Sub-EMF-MAT20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_EMF_MAT
dev.off()

pred_data_EMF_PlasSize <- with(data_EMF, data.frame(
  PlasSize = seq(min(PlasSize, na.rm = TRUE), max(PlasSize, na.rm = TRUE), length.out = 400)))
pred_data_EMF_PlasSize <- cbind(pred_data_EMF_PlasSize, predict(model_EMF_PlasSize, newmods = pred_data_EMF_PlasSize$PlasSize))
r2_val_EMF_PlasSize <- r2_ml(model_EMF_PlasSize)[1]
p_val_EMF_PlasSize <- model_EMF_PlasSize$pval[2]
r2_text_EMF_PlasSize <- ifelse(r2_val_EMF_PlasSize < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_EMF_PlasSize))
p_formatted_EMF_PlasSize <- ifelse(p_val_EMF_PlasSize < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_EMF_PlasSize))

ggplot() +
  geom_point(data = data_EMF, aes(x = PlasSize, y = yi), size = 3, alpha = 0.1, color = "#95B7DA") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_data_EMF_PlasSize, fill = "#95B7DA", alpha = 0.25,
              aes(x = PlasSize, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_data_EMF_PlasSize, linetype = "solid",
            aes(x = PlasSize, y = pred), color = "#95B7DA", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_EMF_PlasSize, p_formatted_EMF_PlasSize, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "Plastic size (mm)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_EMF_PlasSize
ggplot() +
  geom_point(data = subset(data_EMF, PlasSize < 6), aes(x = PlasSize, y = yi), size = 2, alpha = 0.1, color = "#95B7DA") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = subset(pred_data_EMF_PlasSize, PlasSize < 6), fill = "#95B7DA", alpha = 0.25,
              aes(x = PlasSize, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = subset(pred_data_EMF_PlasSize, PlasSize < 6), linetype = "solid",
            aes(x = PlasSize, y = pred), color = "#95B7DA", linewidth = 0.8) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = NULL, y = NULL) + theme_classic() + theme(
    axis.line = element_blank(), panel.border = element_rect(),
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 9),
    plot.margin = margin(0, 0, 0, 0)) -> p_inset_EMF_PlasSize
p_Sub_EMF_PlasSize + inset_element(
  p_inset_EMF_PlasSize, left = 0.45, right = 0.98,
  bottom = 0.01, top = 0.43) -> p_Sub_EMF_PlasSize2; p_Sub_EMF_PlasSize2
cairo_pdf("Output/Sub-EMF-PlasSize20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_EMF_PlasSize
dev.off()
cairo_pdf("Output/Sub-EMF-PlasSize2 20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_EMF_PlasSize2
dev.off()

pred_data_EMF_MAP <- with(data_EMF, data.frame(
  MAP = seq(min(MAP, na.rm = TRUE), max(MAP, na.rm = TRUE), length.out = 400)))
pred_data_EMF_MAP <- cbind(pred_data_EMF_MAP, predict(model_EMF_MAP, newmods = pred_data_EMF_MAP$MAP))
r2_val_EMF_MAP <- r2_ml(model_EMF_MAP)[1]
p_val_EMF_MAP <- model_EMF_MAP$pval[2]
r2_text_EMF_MAP <- ifelse(r2_val_EMF_MAP < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_EMF_MAP))
p_formatted_EMF_MAP <- ifelse(p_val_EMF_MAP < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_EMF_MAP))

ggplot() +
  geom_point(data = data_EMF, aes(x = MAP, y = yi), size = 3, alpha = 0.1, color = "#95B7DA") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_data_EMF_MAP, fill = "#95B7DA", alpha = 0.25,
              aes(x = MAP, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_data_EMF_MAP, linetype = "solid",
            aes(x = MAP, y = pred), color = "#95B7DA", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_EMF_MAP, p_formatted_EMF_MAP, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "MAP (mm)", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_EMF_MAP; p_Sub_EMF_MAP
cairo_pdf("Output/Sub-EMF-MAP20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_EMF_MAP
dev.off()

pred_data_EMF_AridityIndex <- with(data_EMF, data.frame(
  AridityIndex = seq(min(AridityIndex, na.rm = TRUE), max(AridityIndex, na.rm = TRUE), length.out = 400)))
pred_data_EMF_AridityIndex <- cbind(pred_data_EMF_AridityIndex, predict(model_EMF_AridityIndex, newmods = pred_data_EMF_AridityIndex$AridityIndex))
r2_val_EMF_AridityIndex <- r2_ml(model_EMF_AridityIndex)[1]
p_val_EMF_AridityIndex <- model_EMF_AridityIndex$pval[2]
r2_text_EMF_AridityIndex <- ifelse(r2_val_EMF_AridityIndex < 0.001, "italic(R)^{2} < 0.001", sprintf("italic(R)^{2} == %.3f", r2_val_EMF_AridityIndex))
p_formatted_EMF_AridityIndex <- ifelse(p_val_EMF_AridityIndex < 0.001, "italic(P) < 0.001", sprintf("italic(P) == %.3f", p_val_EMF_AridityIndex))

ggplot() +
  geom_point(data = data_EMF, aes(x = AridityIndex, y = yi), size = 3, alpha = 0.1, color = "#95B7DA") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_ribbon(data = pred_data_EMF_AridityIndex, fill = "#95B7DA", alpha = 0.25,
              aes(x = AridityIndex, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(data = pred_data_EMF_AridityIndex, linetype = "twodash",
            aes(x = AridityIndex, y = pred), color = "#95B7DA", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste(r2_text_EMF_AridityIndex, p_formatted_EMF_AridityIndex, sep = "*\", \"~"),
           parse = TRUE, size = 13/2.8346, hjust = 1.02, vjust = 1.2) +
  scale_y_continuous(labels = ~ sprintf("%.1f", .x)) +
  labs(x = "Aridity index", y = "Effect Size") + theme_classic() + theme(
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/2) -> p_Sub_EMF_AridityIndex; p_Sub_EMF_AridityIndex
cairo_pdf("Output/Sub-EMF-AridityIndex20260623.pdf", bg = "transparent", height = 5.3/2.54, width = 30/2.54)
p_Sub_EMF_AridityIndex
dev.off()

# Table S6: moderator tests ----
extract_mod_test <- function(model, moderator_name) {
  s <- summary(model)
  qm_val <- round(s$QM, 2)
  p_val_raw <- s$QMp
  p_val_formatted <- ifelse(p_val_raw < 0.001, "<0.001", sprintf("%.3f", round(p_val_raw, 3)))
  qe_val <- round(s$QE, 2)
  qe_p_raw <- s$QEp
  qe_p_formatted <- ifelse(qe_p_raw < 0.001, "<0.001", sprintf("%.3f", round(qe_p_raw, 3)))
  data.frame(
    Moderator = moderator_name, QM = qm_val, QM_p = p_val_formatted,
    QE = qe_val, QE_p = qe_p_formatted, stringsAsFactors = FALSE)
}
QM_results_Yield <- bind_rows(
  extract_mod_test(model_Yield_PlasType2, "PlasType2"),
  extract_mod_test(model_Yield_PlasSize2, "PlasSize2"),
  extract_mod_test(model_Yield_PlasDose2, "PlasDose2"),
  extract_mod_test(model_Yield_PlasShape, "PlasShape"),
  extract_mod_test(model_Yield_PlasType, "PlasType"),
  extract_mod_test(model_Yield_PlasDose, "PlasDose"),
  extract_mod_test(model_Yield_PlasSize, "PlasSize"),
  extract_mod_test(model_Yield_SampleTime2, "SampleTime2"),
  extract_mod_test(model_Yield_MAT, "MAT"),
  extract_mod_test(model_Yield_MAP, "MAP"),
  extract_mod_test(model_Yield_AridityIndex, "AridityIndex"))
QM_results_SQI <- bind_rows(
  extract_mod_test(model_SQI_PlasType2, "PlasType2"),
  extract_mod_test(model_SQI_PlasSize2, "PlasSize2"),
  extract_mod_test(model_SQI_PlasDose2, "PlasDose2"),
  extract_mod_test(model_SQI_PlasShape, "PlasShape"),
  extract_mod_test(model_SQI_PlantAppear, "PlantAppear"),
  extract_mod_test(model_SQI_PlasType, "PlasType"),
  extract_mod_test(model_SQI_PlasDose, "PlasDose"),
  extract_mod_test(model_SQI_PlasSize, "PlasSize"),
  extract_mod_test(model_SQI_SampleTime2, "SampleTime2"),
  extract_mod_test(model_SQI_MAT, "MAT"),
  extract_mod_test(model_SQI_MAP, "MAP"),
  extract_mod_test(model_SQI_AridityIndex, "AridityIndex"))
QM_results_EMF <- bind_rows(
  extract_mod_test(model_EMF_PlasType2, "PlasType2"),
  extract_mod_test(model_EMF_PlasSize2, "PlasSize2"),
  extract_mod_test(model_EMF_PlasDose2, "PlasDose2"),
  extract_mod_test(model_EMF_PlasShape, "PlasShape"),
  extract_mod_test(model_EMF_PlantAppear, "PlantAppear"),
  extract_mod_test(model_EMF_PlasType, "PlasType"),
  extract_mod_test(model_EMF_PlasDose, "PlasDose"),
  extract_mod_test(model_EMF_PlasSize, "PlasSize"),
  extract_mod_test(model_EMF_SampleTime2, "SampleTime2"),
  extract_mod_test(model_EMF_MAT, "MAT"),
  extract_mod_test(model_EMF_MAP, "MAP"),
  extract_mod_test(model_EMF_AridityIndex, "AridityIndex"))
QM_results_all <- bind_rows(
  mutate(QM_results_Yield, Response = "Yield"),
  mutate(QM_results_SQI, Response = "SQI"),
  mutate(QM_results_EMF, Response = "EMF")
) %>% select(Response, everything())
write.xlsx(QM_results_all, "Output/QM_results20260410.xlsx", rowNames = FALSE)

cat("\014")

pacman::p_load(openxlsx, tidyverse, terra, tidyterra)

# Fig. S4c: global study-site map ----
BaseMap <- rast("D:/ArcGIS/NaturalEarth/HYP_HR_SR_OB_DR.tif")

points_df <- read.xlsx(meta_file, "PointLL20260312") %>%
  filter(Irrigation == "No")

points_df$location <- paste(points_df$Longitude, points_df$Latitude, sep = "_")

point_count <- points_df %>%
  group_by(Ecosystem,
           Group1, location = tolower(location)) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  separate(location, into = c("Longitude", "Latitude"), sep = "_") %>%
  mutate(across(c(Longitude, Latitude, Count), as.numeric))

point_literature <- point_count %>% filter(Group1 == "Literature")

point_literature <- point_literature %>%
  mutate(Count_size = case_when(
    Count == 1 ~ "1",
    Count == 2 ~ "2",
    Count >= 3 ~ "3+"
  ) %>% factor(levels = c("1", "2", "3+")),
  Ecosystem = factor(Ecosystem, levels = c("Agricultural", "Natural"))
  ) %>% arrange(desc(Ecosystem))

ggplot() +
  geom_spatraster_rgb(data = BaseMap) +
  geom_point(data = point_literature, color = "black", stroke = 0.8,
             shape = 21,
             aes(x = Longitude, y = Latitude,
                 fill = Ecosystem,
                 size = Count_size)) +
  scale_size_manual(values = c("1" = 3, "2" = 3.5, "3+" = 4)) +
  scale_x_continuous(
    "Longitude(°)", expand = expansion(mult = 0),
    breaks = c(-120, 0, 120), labels = c("-120", "0", "120")) +
  scale_y_continuous(
    "Latitude(°)", expand = expansion(mult = 0),
    breaks = c(-60, 0, 60), labels = c("-60", "0", "60")) +
  theme_classic() + theme(
    axis.line = element_blank(),
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 16),
    plot.margin = margin(0, 0, 0, 0),
    legend.position = c(0.05, 0.05), legend.justification = c(0, 0),
    legend.box = "vertical",
    legend.background = element_blank(), legend.key = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)) +
  labs(fill = NULL, shape = NULL, size = NULL) + guides(
    fill = guide_legend(
      order = 2, title = "Ecosystem",
      override.aes = list(shape = 21, color = "black", size = 3)),
    size = guide_legend(
      order = 1, title = "Count",
      override.aes = list(shape = 1, color = "black", fill = NA))) -> p_Map1; p_Map1

ggplot() +
  geom_spatraster_rgb(
    data = BaseMap,
    maxcell = 5e+06) +
  geom_point(data = point_literature, color = "black", stroke = 0.8,
             shape = 21,
             aes(x = Longitude, y = Latitude,
                 fill = Ecosystem,
                 size = Count_size)) +
  scale_size_manual(values = c("1" = 3, "2" = 3.5, "3+" = 4)) +
  scale_x_continuous(
    "Longitude(°)", expand = expansion(mult = 0),
    breaks = c(-120, 0, 120), labels = c("-120", "0", "120")) +
  scale_y_continuous(
    "Latitude(°)", expand = expansion(mult = 0),
    breaks = c(-60, 0, 60), labels = c("-60", "0", "60")) +
  theme_classic() + theme(
    axis.line = element_blank(),
    plot.background = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 15),
    axis.title = element_text(color = "black", size = 16),
    plot.margin = margin(0, 0, 0, 0),
    legend.position = c(0.05, 0.05), legend.justification = c(0, 0),
    legend.box = "vertical",
    legend.background = element_blank(), legend.key = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)) +
  labs(fill = NULL, shape = NULL, size = NULL) + guides(
    fill = guide_legend(
      order = 2, title = "Ecosystem",
      override.aes = list(shape = 21, color = "black", size = 3)),
    size = guide_legend(
      order = 1, title = "Count",
      override.aes = list(shape = 1, color = "black", fill = NA))) -> p_Map2
pdf("Output/MapLocation20260315.pdf", width = 27/2.54)
p_Map2
dev.off()

pacman::p_load_gh("valentinitnelav/plotbiomes")
pacman::p_load(openxlsx, tidyverse, plotbiomes)

# Fig. S4b: Whittaker climate-space panel ----
points <- read.xlsx(meta_file, "Points_with_Climate_Elevation") %>%
  filter(Irrigation == "No")

points_literature <- points %>% filter(Group1 == "Literature")

points_literature <- points_literature %>% mutate(
  Ecosystem = factor(Ecosystem, levels = c("Agricultural", "Natural"))
) %>% arrange(desc(Ecosystem))

whittaker_base_plot() +
  geom_point(data = points_literature,
             aes(x = Annual_Mean_Temp,
                 y = Annual_Precip/10
             ), shape = 21,
             fill = "white", color = "black", stroke = 0.8, size = 3) +
  scale_shape_manual(
    values = c("No" = 21, "Irrigation" = 22, "Control" = 23, "Flood" = 24)) +
  theme_classic() + theme(
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 15),
    axis.line = element_blank(), panel.border = element_rect(),
    plot.background = element_blank(), panel.background = element_blank(),
    legend.position = "right", legend.justification = "bottom",
    legend.box = "vertical",
    legend.background = element_blank(), legend.key = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    plot.margin = margin(0, 0, 0, 0), aspect.ratio = 1/1) +
  labs(fill = NULL, shape = NULL) + guides(
    fill = guide_legend(
      order = 1, title = "Whittaker biomes", override.aes = list(size = 3))
  ) -> p_Whittaker; p_Whittaker
pdf("Output/MapWhittaker20260315.pdf", width = 20/2.54)
p_Whittaker
dev.off()
