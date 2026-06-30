# =============================================================================
# R SCRIPT 1: skNIR Haploid Sorting Analysis (Two Oil-based Sorting Methods & Three Multivariate Methods)
# Paper: Single-kernel near-infrared spectroscopy enables haploid kernel sorting in field and sweet corn using high-oil haploid inducers across diverse donor–inducer combinations
# =============================================================================
# ORGANIZATION:
#   PART 0: Libraries, data loading, and shared setup
#   PART 1: Multivariate functions for the three methods (PLS-DA, SVM, RF)
#   PART 2: Oil-based haploid sorting methods
#   PART 3: Multivariate haploid sorting methods
#   PART 4: Defective kernel sorting 
#   PART 5: Make Figure 3 (based on the results from PART 2)
# =============================================================================

library(dplyr)
library(tidyr)
library(waves)
library(caret)
library(ggplot2)
library(ggtext)

MASTER_SEED <- 2026

RNGkind("Mersenne-Twister", "Inversion")

set.seed(MASTER_SEED)

cat("Global reproducibility seed:", MASTER_SEED, "\n")

# Data loading 
dat <- read.csv("/home/sumans/26June_skNIR_finalcheck/Data_skNIR_Manuscript_June27.2026.csv")

# Spectral columns (940-1640 nm, 1 nm step)
spec_cols <- paste0("X", 940:1640)

# Apply SNV pretreatment to all spectra
dat[, spec_cols] <- pretreat_spectra(
  df = dat[, spec_cols],
  pretreatment = 2
)

# Haploid induction rates (HIRs) for each population
induction_rates <- data.frame(
  population = as.character(1:12),
  IR = c(0.08, 0.13, 0.10, 0.11, 0.10, 0.04, 0.03, 0.08, 0.08, 0.09, 0.08, 0.12),
  stringsAsFactors = FALSE)

# Population labels
pop_labels <- dat %>%
  distinct(population, Donor, Inducer) %>%
  mutate(
    pop_label     = paste0(Donor, " × ", Inducer, " (", population, ")")
  )

# Population IDs in order 1-12
populations <- sort(unique(dat$population))

# Model data: kernels with confirmed gDNA score for haploid and hybrids (dna.score.final = 1 and 4, respectively) 
# Defectives using visual scoring (dna.score.final = 5)
dat_model <- dat %>%
  filter(!is.na(dna.score.final)) %>%
  filter(dna.score.final %in% c(1, 4, 5))

# Training data: haploid and hybrid only (no defectives)
dat_train_only <- dat_model %>%
  filter(dna.score.final %in% c(1, 4))

# Operational benchmark thresholds for FDR and FNR
FDR_TARGET <- 0.25
FNR_MAX    <- 0.50

cat("Data loaded:", nrow(dat), "kernels across", length(populations), "populations\n")
cat("Kernel class breakdown:\n")
print(table(dat_model$dna.score.final))


# =============================================================================
# SHARED HELPER FUNCTIONS
# =============================================================================

# --- Resample population to target haploid induction rate (HIR) ---
# Includes defectives (dna.score.final = 5), downsampled at same rate as haploids when required
# FDR/FNR should be calculated on classes 1 (haploid) and 4 (hybrid) only
resample_to_IR <- function(data, ir, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  haploids   <- data %>% filter(dna.score.final == 1)
  hybrids    <- data %>% filter(dna.score.final == 4)
  defectives <- data %>% filter(dna.score.final == 5)
  
  n_hap        <- nrow(haploids)
  n_hyb        <- nrow(hybrids)
  n_def        <- nrow(defectives)
  n_hyb_target <- round(n_hap * ((1 - ir) / ir))
  
  if (n_hyb >= n_hyb_target) {
    # Enough hybrids: keep all haploids and defectives, undersample hybrids
    hybrids <- hybrids %>% slice_sample(n = n_hyb_target)
  } else {
    # Not enough hybrids: keep all hybrids, undersample haploids
    # Downsample defectives at same rate as haploids
    n_hap_target    <- round(n_hyb * (ir / (1 - ir)))
    downsample_rate <- n_hap_target / n_hap
    haploids        <- haploids %>% slice_sample(n = n_hap_target)
    if (n_def > 0) {
      n_def_target <- max(1, round(n_def * downsample_rate))
      defectives   <- defectives %>% slice_sample(n = n_def_target)
    }
  }
  
  bind_rows(haploids, hybrids, defectives)
}

# --- Compute FDR and FNR (classes 1 (haploid) and 4 (hybrid) only) ---
# NaN FDR -> 0 (no kernels predicted haploid = very conservative = FDR of 0)
# NaN FNR -> 1 (no kernels predicted haploid = all haploids missed)
compute_metrics <- function(predicted, actual) {
  if (is.list(predicted)) predicted <- predicted$predictions
  predicted <- as.integer(predicted)
  actual    <- as.integer(actual)
  idx  <- actual %in% c(1, 4)
  pred <- predicted[idx]
  act  <- actual[idx]
  TP <- sum(pred == 1 & act == 1)
  FP <- sum(pred == 1 & act == 4)
  FN <- sum(pred == 4 & act == 1)
  total_pred_haploid <- TP + FP
  total_true_haploid <- TP + FN
  FDR <- ifelse(total_pred_haploid > 0, FP / total_pred_haploid, 0)
  FNR <- ifelse(total_true_haploid > 0, FN / total_true_haploid, 1)
  data.frame(
    FDR            = round(FDR, 4),
    FNR            = round(FNR, 4),
    combined       = round(FDR + FNR, 4),
    n_pred_haploid = total_pred_haploid,
    n_true_haploid = total_true_haploid
  )
}

# --- Defective kernel metrics ---
# FDR = defectives misclassified as non-defective / all predicted non-defective
# FNR = non-defectives misclassified as defective / all true non-defective
compute_def_metrics <- function(predicted, actual) {
  predicted <- as.integer(predicted)
  actual    <- as.integer(actual)
  # 1 = defective, 4 = non-defective (haploid or hybrid relabeled)
  TP <- sum(predicted == 4 & actual == 4)   # correctly retained non-defective
  FP <- sum(predicted == 4 & actual == 1)   # defective in non-defective pool (FDR)
  FN <- sum(predicted == 1 & actual == 4)   # non-defective sent to defective pool (FNR)
  total_pred_nondef <- TP + FP
  total_true_nondef <- TP + FN
  FDR <- ifelse(total_pred_nondef > 0, FP / total_pred_nondef, 0)
  FNR <- ifelse(total_true_nondef > 0, FN / total_true_nondef, 1)
  data.frame(FDR = round(FDR, 4), FNR = round(FNR, 4))
}


# =============================================================================
# PART 1: MULTIVARIATE FUNCTIONS FOR THE THREE METHODS
# =============================================================================

# --- 5-fold cross-validation framework ---
CV_5 <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "final",
  allowParallel = FALSE
)

# --- PLS-DA classifier ---
# ncomp = 12
run_plsda <- function(X_train, y_train, X_test,
                      ncomp = 12,
                      seed = MASTER_SEED) {
  
  set.seed(seed)
  
  train_df <- data.frame(
    Class = factor(
      ifelse(y_train == 1, "Class1", "Class2"),
      levels = c("Class1", "Class2")
    ),
    X_train
  )
  
  pls_grid <- expand.grid(
    ncomp = ncomp
  )
  
  fit <- train(
    Class ~ .,
    data = train_df,
    method = "pls",
    tuneGrid = pls_grid,
    metric = "Kappa",
    trControl = CV_5
  )
  
  pred <- predict(
    fit,
    newdata = as.data.frame(X_test)
  )
  
  ifelse(pred == "Class1", 1, 4)
}

# --- SVM classifier ---
# Radial kernel
run_svm <- function(X_train, y_train, X_test,
                    seed = MASTER_SEED) {
  
  set.seed(seed)
  
  train_df <- data.frame(
    Class = factor(
      ifelse(y_train == 1, "Class1", "Class2"),
      levels = c("Class1", "Class2")
    ),
    X_train
  )
  
  svm_grid <- expand.grid(
    C = c(1, 10, 100),
    sigma = c(0.0001, 0.001, 0.01)
  )
  
  fit <- train(
    Class ~ .,
    data = train_df,
    method = "svmRadial",
    tuneGrid = svm_grid,
    metric = "Kappa",
    trControl = CV_5
  )
  
  pred <- predict(
    fit,
    newdata = as.data.frame(X_test)
  )
  
  list(
    predictions = ifelse(pred == "Class1", 1, 4),
    best_c = fit$bestTune$C,
    best_sigma = fit$bestTune$sigma,
    cv_results = fit$results
  )
}

# --- Random Forest classifier ---
# 200 trees, mtry tuned over grid using 50 preliminary trees
run_rf <- function(X_train, y_train, X_test,
                   ntree = 200,
                   tune_ntree = 50,
                   seed = MASTER_SEED) {

  set.seed(seed)

  train_df <- data.frame(
    Class = factor(
      ifelse(y_train == 1, "Class1", "Class2"),
      levels = c("Class1", "Class2")
    ),
    X_train
  )

  n_feat <- ncol(X_train)

  mtry_vals <- unique(c(
    10,
    50,
    100,
    round(sqrt(n_feat))
  ))

  mtry_vals <- mtry_vals[mtry_vals <= n_feat]

  rf_grid <- expand.grid(
    mtry = mtry_vals
  )

  # ---- Tune mtry using 50-tree forests ----
  fit_tune <- train(
    Class ~ .,
    data = train_df,
    method = "rf",
    tuneGrid = rf_grid,
    ntree = tune_ntree,
    metric = "Kappa",
    trControl = CV_5
  )

  best_mtry <- fit_tune$bestTune$mtry

  # ---- Final model using all training data ----
  set.seed(seed + 999)

  fit_final <- train(
    Class ~ .,
    data = train_df,
    method = "rf",
    tuneGrid = data.frame(mtry = best_mtry),
    ntree = ntree,
    trControl = trainControl(method = "none")
)

  pred <- predict(
    fit_final,
    newdata = as.data.frame(X_test)
  )

  list(
    predictions = ifelse(pred == "Class1", 1, 4),
    best_mtry = best_mtry,
    tuning_table = fit_tune$results,
    final_model = fit_final
  )
}

# =============================================================================
# PART 2: OIL-BASED HAPLOID SORTING METHODS
# =============================================================================

# --- 2A: Per-population optimal threshold ---
# Sweep 500 thresholds; at each, average FDR/FNR over n_iter resamplings
# Select threshold minimizing FNR subject to FDR <= 0.25, FNR <= 0.50
optimize_population_resampled <- function(pop_data, pop_ir,
                                          fdr_target = FDR_TARGET,
                                          fnr_max    = FNR_MAX,
                                          n_thresh   = 500,
                                          n_iter     = 5, 
                                          base_seed = MASTER_SEED) {
  valid_full <- pop_data %>% filter(dna.score.final %in% c(1, 4, 5))
  oil_range  <- range(valid_full$pred.oil)
  thresholds <- seq(oil_range[1], oil_range[2], length.out = n_thresh)
  
  curve <- lapply(seq_along(thresholds), function(i) {
    
    t <- thresholds[i]
    iter_results <- lapply(1:n_iter, function(j) {
      iter_seed <- base_seed + i * 1000 + j
      resampled <- resample_to_IR(
        valid_full,
        ir = pop_ir,
        seed = iter_seed
      )
      pred_haploid <- resampled$pred.oil < t
      TP <- sum(pred_haploid & resampled$dna.score.final == 1)
      FP <- sum(pred_haploid & resampled$dna.score.final == 4)
      FN <- sum(!pred_haploid & resampled$dna.score.final == 1)
      total_pred <- TP + FP; total_true <- TP + FN
      FDR <- ifelse(total_pred > 0, FP / total_pred, 0)
      FNR <- ifelse(total_true > 0, FN / total_true, 1)
      data.frame(FDR = FDR, FNR = FNR)
    })
    iter_df <- do.call(rbind, iter_results)
    data.frame(threshold = t,
               FDR = round(mean(iter_df$FDR), 4),
               FNR = round(mean(iter_df$FNR), 4))
  })
  curve <- do.call(rbind, curve)
  curve <- curve[complete.cases(curve), ]
  
  feasible <- curve[curve$FDR <= fdr_target & curve$FNR <= fnr_max, ]
  if (nrow(feasible) > 0) {
    opt <- feasible[which.min(feasible$FNR), ]; fdr_achieved <- TRUE
  } else {
    fnr_feasible <- curve[curve$FNR <= fnr_max, ]
    opt <- if (nrow(fnr_feasible) > 0) fnr_feasible[which.min(fnr_feasible$FDR), ] else
      curve[which.min(curve$FDR + curve$FNR), ]
    fdr_achieved <- FALSE
  }
  list(threshold = round(opt$threshold, 3), FDR = round(opt$FDR, 3),
       FNR = round(opt$FNR, 3), combined = round(opt$FDR + opt$FNR, 3),
       fdr_achieved = fdr_achieved, curve = curve)
}

# --- 2B: Lowest percentage method ---
eval_lowest_pct_resampled <- function(pop_data, pop_ir, pct,
                                      n_iter = 20,
                                      base_seed = MASTER_SEED) {
  valid_full <- pop_data %>% filter(dna.score.final %in% c(1, 4, 5))
  iter_results <- lapply(1:n_iter, function(j) {
    resampled <- resample_to_IR(
      valid_full,
      ir = pop_ir,
      seed = base_seed + j
    )
    n_total          <- nrow(resampled)
    n_hap_pred       <- floor(n_total * pct / 100)
    resampled_sorted <- resampled %>% arrange(pred.oil)
    threshold        <- resampled_sorted$pred.oil[n_hap_pred]
    pred_haploid     <- resampled$pred.oil < threshold
    TP <- sum(pred_haploid & resampled$dna.score.final == 1)
    FP <- sum(pred_haploid & resampled$dna.score.final == 4)
    FN <- sum(!pred_haploid & resampled$dna.score.final == 1)
    total_pred <- TP + FP; total_true <- TP + FN
    FDR <- ifelse(total_pred > 0, FP / total_pred, 0)
    FNR <- ifelse(total_true > 0, FN / total_true, 1)
    data.frame(FDR = FDR, FNR = FNR)
  })
  iter_df <- do.call(rbind, iter_results)
  data.frame(pct = pct,
             FDR = round(mean(iter_df$FDR), 3),
             FNR = round(mean(iter_df$FNR), 3),
             combined = round(mean(iter_df$FDR + iter_df$FNR), 3))
}

# --- 2C: Dynamic threshold method ---
eval_dynamic_resampled <- function(pop_data, pop_ir, offset, n_iter = 20, base_seed = MASTER_SEED) {
  valid_full <- pop_data %>% filter(dna.score.final %in% c(1, 4, 5))
  iter_results <- lapply(1:n_iter, function(j) {
    resampled <- resample_to_IR(
      valid_full,
      ir = pop_ir,
      seed = base_seed + j
    )
    threshold    <- mean(resampled$pred.oil) - offset
    pred_haploid <- resampled$pred.oil < threshold
    TP <- sum(pred_haploid & resampled$dna.score.final == 1)
    FP <- sum(pred_haploid & resampled$dna.score.final == 4)
    FN <- sum(!pred_haploid & resampled$dna.score.final == 1)
    total_pred <- TP + FP; total_true <- TP + FN
    FDR <- ifelse(total_pred > 0, FP / total_pred, 0)
    FNR <- ifelse(total_true > 0, FN / total_true, 1)
    data.frame(FDR = FDR, FNR = FNR)
  })
  iter_df <- do.call(rbind, iter_results)
  data.frame(offset = offset,
             FDR = round(mean(iter_df$FDR), 3),
             FNR = round(mean(iter_df$FNR), 3),
             combined = round(mean(iter_df$FDR + iter_df$FNR), 3))
}

# --- 2D: Run full oil-based analysis ---
run_oil_analysis <- function(pct_values   = seq(4, 20, by = 2),
                             offset_values = seq(0.8, 2.4, by = 0.2),
                             save_file    = "oil_results_v5.RData") {
  cat("\n--- PART 2: Oil-based haploid sorting ---\n")
  set.seed(MASTER_SEED)
  
  # Step 1: Per-population optimal
  cat("Running per-population optimal threshold...\n")
  optimal_list <- lapply(populations, function(pop) {
    pop_ir   <- induction_rates$IR[induction_rates$population == pop]
    cat("  Population:", pop, "| IR:", round(pop_ir * 100, 1), "%\n")
    res      <- optimize_population_resampled(dat[dat$population == pop, ], pop_ir)
    data.frame(population = pop, IR = pop_ir,
               threshold_optimal = res$threshold,
               FDR_optimal = res$FDR, FNR_optimal = res$FNR,
               combined_optimal = res$combined, fdr_achieved = res$fdr_achieved)
  })
  optimal_table <- do.call(rbind, optimal_list) %>% left_join(pop_labels, by = "population")
  cat("Per-population optimal results:\n"); print(optimal_table, row.names = FALSE)
  write.csv(optimal_table,
          "oil_optimal_table.csv",
          row.names = FALSE)
  
  # Step 2: Identify benchmark populations
  unsortable <- optimal_table %>% filter(!fdr_achieved) %>% pull(population)
  sortable   <- optimal_table %>% filter(fdr_achieved)  %>% pull(population)
  cat("\nPopulations meeting benchmark:", paste(sortable, collapse = ", "), "\n")
  cat("Populations NOT meeting benchmark:", paste(unsortable, collapse = ", "), "\n")
  
  # Step 3: Lowest percentage method
  cat("\nRunning lowest percentage method...\n")
  pct_results <- lapply(populations, function(pop) {
    pop_ir <- induction_rates$IR[induction_rates$population == pop]
    res    <- lapply(pct_values, function(p)
      eval_lowest_pct_resampled(dat[dat$population == pop, ], pop_ir, p))
    do.call(rbind, res) %>% mutate(population = pop)
  }) %>% do.call(rbind, .) %>% mutate(sortable = !population %in% unsortable)
  
  pct_selection <- pct_results %>%
    group_by(pct) %>%
    summarise(mean_FDR    = round(mean(FDR), 3),
              mean_FNR    = round(mean(FNR), 3),
              max_FDR     = round(max(FDR[sortable], na.rm = TRUE), 3),
              max_FNR     = round(max(FNR[sortable], na.rm = TRUE), 3),
              n_both_pass = sum(FDR[sortable] <= FDR_TARGET &
                                  FNR[sortable] <= FNR_MAX, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(pct)
  cat("Lowest percentage selection summary:\n"); print(pct_selection, row.names = FALSE)
  write.csv(pct_selection,
          "oil_pct_selection.csv",
          row.names = FALSE)
  
  optimal_pct <- pct_selection %>% filter(max_FDR <= FDR_TARGET & max_FNR <= FNR_MAX) %>%
    slice_min(mean_FDR, n = 1)
  if (nrow(optimal_pct) == 0)
    optimal_pct <- pct_selection %>% filter(n_both_pass == max(n_both_pass)) %>%
    slice_min(mean_FDR, n = 1)
  cat("Optimal lowest percentage:", optimal_pct$pct, "%\n")
  
  pct_final <- pct_results %>% filter(pct == optimal_pct$pct) %>%
    left_join(pop_labels, by = "population") %>%
    dplyr::select(population, Donor, Inducer, pop_label,
                  FDR_pct = FDR, FNR_pct = FNR, combined_pct = combined, sortable)
  
  # Step 4: Dynamic threshold method
  cat("\nRunning dynamic threshold method...\n")
  dynamic_results <- lapply(populations, function(pop) {
    pop_ir <- induction_rates$IR[induction_rates$population == pop]
    res    <- lapply(offset_values, function(o)
      eval_dynamic_resampled(dat[dat$population == pop, ], pop_ir, o))
    do.call(rbind, res) %>% mutate(population = pop)
  }) %>% do.call(rbind, .) %>% mutate(sortable = !population %in% unsortable)
  
  dynamic_selection <- dynamic_results %>%
    group_by(offset) %>%
    summarise(mean_FDR    = round(mean(FDR), 3),
              mean_FNR    = round(mean(FNR), 3),
              max_FDR     = round(max(FDR[sortable], na.rm = TRUE), 3),
              max_FNR     = round(max(FNR[sortable], na.rm = TRUE), 3),
              n_both_pass = sum(FDR[sortable] <= FDR_TARGET &
                                  FNR[sortable] <= FNR_MAX, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(offset)
  cat("Dynamic threshold selection summary:\n"); print(dynamic_selection, row.names = FALSE)
  write.csv(dynamic_selection,
          "oil_dynamic_selection.csv",
          row.names = FALSE)
  
  optimal_offset <- dynamic_selection %>% filter(max_FDR <= FDR_TARGET & max_FNR <= FNR_MAX) %>%
    slice_min(mean_FDR, n = 1)
  if (nrow(optimal_offset) == 0)
    optimal_offset <- dynamic_selection %>% filter(n_both_pass == max(n_both_pass)) %>%
    slice_min(mean_FDR, n = 1)
  cat("Optimal dynamic threshold offset:", optimal_offset$offset, "\n")
  
  dynamic_final <- dynamic_results %>% filter(offset == optimal_offset$offset) %>%
    left_join(pop_labels, by = "population") %>%
    dplyr::select(population, Donor, Inducer, pop_label,
                  FDR_dyn = FDR, FNR_dyn = FNR, combined_dyn = combined, sortable)
  
  # Step 5: Consolidated comparison table
  comparison_table <- optimal_table %>%
    dplyr::select(population, Donor, Inducer, pop_label, IR,
                  FDR_optimal, FNR_optimal, combined_optimal, fdr_achieved) %>%
    left_join(pct_final    %>% dplyr::select(population, FDR_pct, FNR_pct, combined_pct),
              by = "population") %>%
    left_join(dynamic_final %>% dplyr::select(population, FDR_dyn, FNR_dyn, combined_dyn),
              by = "population") %>%
    mutate(sortable = !population %in% unsortable)
  cat("\nConsolidated comparison table:\n"); print(comparison_table, row.names = FALSE)
  write.csv(comparison_table,
          "oil_comparison_table.csv",
          row.names = FALSE)
  
  # Summary statistics
  cat("\nSummary - populations meeting benchmark:\n")
  summary_sortable <- comparison_table %>%
  filter(sortable) %>%
  summarise(
    mean_FDR_optimal = round(mean(FDR_optimal), 3),
    mean_FNR_optimal = round(mean(FNR_optimal), 3),
    mean_FDR_pct     = round(mean(FDR_pct), 3),
    mean_FNR_pct     = round(mean(FNR_pct), 3),
    mean_FDR_dyn     = round(mean(FDR_dyn), 3),
    mean_FNR_dyn     = round(mean(FNR_dyn), 3)
  )

print(summary_sortable)

write.csv(summary_sortable,
          "oil_summary_sortable.csv",
          row.names = FALSE)
  
  cat("\nSummary - all 12 populations:\n")
  summary_all <- comparison_table %>%
  summarise(
    mean_FDR_optimal = round(mean(FDR_optimal), 3),
    mean_FNR_optimal = round(mean(FNR_optimal), 3),
    mean_FDR_pct     = round(mean(FDR_pct), 3),
    mean_FNR_pct     = round(mean(FNR_pct), 3),
    mean_FDR_dyn     = round(mean(FDR_dyn), 3),
    mean_FNR_dyn     = round(mean(FNR_dyn), 3)
  )

print(summary_all)

write.csv(summary_all,
          "oil_summary_all.csv",
          row.names = FALSE)
  
  # Save
  oil_non_benchmark <- unsortable
  save(comparison_table, optimal_table, pct_final, dynamic_final,
       pct_results, dynamic_results, pct_selection, dynamic_selection,
       optimal_pct, optimal_offset, sortable, unsortable,
       oil_non_benchmark, induction_rates,
       file = save_file)
  cat("Oil results saved to", save_file, "\n")
  
  invisible(list(comparison_table = comparison_table, optimal_table = optimal_table,
                 pct_selection = pct_selection, dynamic_selection = dynamic_selection,
                 sortable = sortable, unsortable = unsortable,
                 optimal_pct = optimal_pct, optimal_offset = optimal_offset))
}


# =============================================================================
# PART 3: MULTIVARIATE HAPLOID SORTING METHODS
# =============================================================================

# --- 3A: General model - Hold-One-Population-Out CV ---
# Training: all 11 non-held-out populations using all haploid and hybrid kernels
# Test: held-out population resampled to its own IR, with defectives included
run_hopo_cv <- function(pops, label = "", n_rep = 20) {
  cat("\n--- HOPO-CV:", label, "---\n")
  cat("Training populations:", paste(pops, collapse = ", "), "\n")
  cat("Training: all haploid/hybrid kernels, Test IR: population-specific\n")
  cat("Replications:", n_rep, "\n")
  
  all_rep_results <- list()
  
  for (rep in 1:n_rep) {
    set.seed(MASTER_SEED + rep)
    cat("\n  Replication:", rep, "/", n_rep, "\n")
    results <- list()
    
    for (held_out in pops) {
      train_pops <- pops[pops != held_out]
      
      # Training: ALL haploid and hybrid kernels from 11 non-held-out populations
      # No IR resampling for training
      train_data <- dat_train_only %>% filter(population %in% train_pops)
      
      # Test: resample held-out to its own IR, add defectives
      test_ir   <- induction_rates$IR[induction_rates$population == held_out]
      test_pool <- dat_train_only %>% filter(population == held_out)
      local_seed <- MASTER_SEED + rep * 1000
      test_data_resampled <- resample_to_IR(
        test_pool,
        ir = test_ir,
        seed = local_seed
      )
      
      defectives <- dat_model %>% filter(population == held_out, dna.score.final == 5)
      if (nrow(defectives) > 0) {
        n_def_test <- max(1, floor(nrow(defectives) * 0.3))
        def_sample <- defectives %>% slice_sample(n = n_def_test)
        test_data  <- bind_rows(test_data_resampled, def_sample)
      } else {
        test_data <- test_data_resampled
      }
      
      X_train <- train_data[, spec_cols]
      y_train <- train_data$dna.score.final
      X_test  <- test_data[, spec_cols]
      y_test  <- test_data$dna.score.final
      
      pred_plsda <- tryCatch(run_plsda(X_train, y_train, X_test),
                             error = function(e) rep(NA_integer_, nrow(test_data)))
      svm_out    <- tryCatch(run_svm(X_train, y_train, X_test, seed = local_seed + 100),
                             error = function(e) list(predictions = rep(NA_integer_, nrow(test_data))))
      rf_out     <- tryCatch(run_rf(X_train, y_train, X_test, seed = local_seed + 200),
                             error = function(e) list(predictions = rep(NA_integer_, nrow(test_data))))
      
      m_plsda <- compute_metrics(pred_plsda,          y_test)
      m_svm   <- compute_metrics(svm_out$predictions, y_test)
      m_rf    <- compute_metrics(rf_out$predictions,  y_test)
      
      results[[held_out]] <- data.frame(
        population  = held_out, replication = rep, model_config = label,
        validation = "HOPO-CV", test_IR = test_ir,
        FDR_plsda = m_plsda$FDR, FNR_plsda = m_plsda$FNR,
        FDR_svm   = m_svm$FDR,   FNR_svm   = m_svm$FNR,
        FDR_rf    = m_rf$FDR,    FNR_rf    = m_rf$FNR
      )
    }
    all_rep_results[[rep]] <- do.call(rbind, results)
  }
  
  all_reps_df  <- do.call(rbind, all_rep_results)
  results_avg  <- all_reps_df %>%
    group_by(population, model_config, validation, test_IR) %>%
    summarise(FDR_plsda = round(mean(FDR_plsda, na.rm = TRUE), 3),
              FNR_plsda = round(mean(FNR_plsda, na.rm = TRUE), 3),
              FDR_svm   = round(mean(FDR_svm,   na.rm = TRUE), 3),
              FNR_svm   = round(mean(FNR_svm,   na.rm = TRUE), 3),
              FDR_rf    = round(mean(FDR_rf,    na.rm = TRUE), 3),
              FNR_rf    = round(mean(FNR_rf,    na.rm = TRUE), 3),
              .groups = "drop")
  cat("\nHOPO-CV complete. Results:\n"); print(results_avg, row.names = FALSE)
  
  write.csv(results_avg,
          paste0(label, "_hopo_results.csv"),
          row.names = FALSE)
          
   cat("\nMean across all populations:\n")      
    hopo_summary <- results_avg %>%
                    summarise(
                    mean_FDR_plsda = round(mean(FDR_plsda, na.rm = TRUE), 3),
                    mean_FNR_plsda = round(mean(FNR_plsda, na.rm = TRUE), 3),
                    mean_FDR_svm   = round(mean(FDR_svm,   na.rm = TRUE), 3),
                    mean_FNR_svm   = round(mean(FNR_svm,   na.rm = TRUE), 3),
                    mean_FDR_rf    = round(mean(FDR_rf,    na.rm = TRUE), 3),
                    mean_FNR_rf    = round(mean(FNR_rf,    na.rm = TRUE), 3)
  )

cat("\nHOPO-CV mean across populations:\n")
print(hopo_summary)
          
invisible(list(results = results_avg, all_reps = all_reps_df))
}

# --- 3B: Population-specific models ---
# For each of 100 iterations: resample to pop IR, 70/30 split for train and test
run_pop_specific <- function(n_iter = 100) {
  cat("\n--- Population-specific models ---\n")
  cat("Resampling: population-specific IR within each iteration\n")
  cat("Iterations:", n_iter, "\n")
  
  all_pops <- sort(unique(dat_train_only$population))
  results  <- list()
  set.seed(MASTER_SEED)
  
  for (pop in all_pops) {
    cat("  Population:", pop, "\n")
    pop_ir    <- induction_rates$IR[induction_rates$population == pop]
    pop_data  <- dat_model %>% filter(population == pop)
    pop_train <- pop_data  %>% filter(dna.score.final %in% c(1, 4))
    
    if (sum(pop_train$dna.score.final == 1) < 4) {
      cat("    Too few haploids, skipping.\n"); next
    }
    
    iter_results <- list()
    for (j in 1:n_iter) {
      local_seed <- MASTER_SEED + j
      set.seed(local_seed)
      pop_resampled <- resample_to_IR(
        pop_train,
        ir = pop_ir,
        seed = local_seed
      )
      hap_idx <- which(pop_resampled$dna.score.final == 1)
      hyb_idx <- which(pop_resampled$dna.score.final == 4)
      train_hap_idx <- sample(hap_idx, floor(length(hap_idx) * 0.7))
      train_hyb_idx <- sample(hyb_idx, floor(length(hyb_idx) * 0.7))
      train_idx <- c(train_hap_idx, train_hyb_idx)
      test_idx  <- setdiff(seq_len(nrow(pop_resampled)), train_idx)
      
      X_train     <- pop_resampled[train_idx, spec_cols]
      y_train     <- pop_resampled$dna.score.final[train_idx]
      X_test_base <- pop_resampled[test_idx, spec_cols]
      y_test_base <- pop_resampled$dna.score.final[test_idx]
      
      # Add defectives to test set (30% of available defectives)
      class5_idx <- which(pop_data$dna.score.final == 5)
      if (length(class5_idx) > 0) {
        n_def_test   <- max(1, floor(length(class5_idx) * 0.3))
        def_test_idx <- sample(class5_idx, n_def_test)
        X_test <- rbind(X_test_base, pop_data[def_test_idx, spec_cols])
        rownames(X_test) <- NULL
        y_test       <- c(y_test_base, pop_data$dna.score.final[def_test_idx])
      } else {
        X_test <- X_test_base; rownames(X_test) <- NULL; y_test <- y_test_base
      }
      
      pred_plsda <- tryCatch(run_plsda(X_train, y_train, X_test),
                             error = function(e) rep(NA_integer_, nrow(X_test)))
      
      svm_out <- tryCatch(run_svm(X_train, y_train, X_test, seed = local_seed + 100),
                             error = function(e) list(predictions = rep(NA_integer_, nrow(X_test))))
      
      rf_out <- tryCatch(run_rf(X_train, y_train, X_test, seed = local_seed + 200),
                             error = function(e) list(predictions = rep(NA_integer_, nrow(X_test))))
      
      m_plsda <- compute_metrics(pred_plsda, y_test)
      m_svm   <- compute_metrics(svm_out$predictions, y_test)
      m_rf    <- compute_metrics(rf_out$predictions,  y_test)
      
      iter_results[[j]] <- data.frame(
        population = pop, iteration = j,
        FDR_plsda = m_plsda$FDR, FNR_plsda = m_plsda$FNR,
        FDR_svm   = m_svm$FDR,   FNR_svm   = m_svm$FNR,
        FDR_rf    = m_rf$FDR,    FNR_rf    = m_rf$FNR
      )
    }
    
    iter_df <- do.call(rbind, iter_results)
    results[[pop]] <- iter_df %>%
      group_by(population) %>%
      summarise(model_config = "population-specific",
                validation   = "repeated-split",
                pop_IR       = pop_ir,
                FDR_plsda = round(mean(FDR_plsda, na.rm = TRUE), 3),
                FNR_plsda = round(mean(FNR_plsda, na.rm = TRUE), 3),
                FDR_svm   = round(mean(FDR_svm,   na.rm = TRUE), 3),
                FNR_svm   = round(mean(FNR_svm,   na.rm = TRUE), 3),
                FDR_rf    = round(mean(FDR_rf,    na.rm = TRUE), 3),
                FNR_rf    = round(mean(FNR_rf,    na.rm = TRUE), 3),
                .groups = "drop")
    cat("    FDR: plsda =", results[[pop]]$FDR_plsda,
        "svm =", results[[pop]]$FDR_svm,
        "rf =", results[[pop]]$FDR_rf, "\n")
  }
  
  results_df <- do.call(rbind, results)
  cat("\nPopulation-specific results:\n"); print(results_df, row.names = FALSE)
  
  write.csv(results_df,
          "population_specific_results.csv",
          row.names = FALSE)
          
  cat("\nMean across all populations:\n")    
  
  popspec_summary <- results_df %>%
                     summarise(
                     mean_FDR_plsda = round(mean(FDR_plsda, na.rm = TRUE), 3),
                     mean_FNR_plsda = round(mean(FNR_plsda, na.rm = TRUE), 3),
                     mean_FDR_svm   = round(mean(FDR_svm,   na.rm = TRUE), 3),
                     mean_FNR_svm   = round(mean(FNR_svm,   na.rm = TRUE), 3),
                     mean_FDR_rf    = round(mean(FDR_rf,    na.rm = TRUE), 3),
                     mean_FNR_rf    = round(mean(FNR_rf,    na.rm = TRUE), 3)
  )

print(popspec_summary)

  invisible(list(results = results_df))
}

# --- 3C: Run full multivariate analysis ---
run_multivariate_analysis <- function(save_file = "multivariate_results_v5.RData") {
  cat("\n--- PART 4: Multivariate haploid sorting ---\n")
  set.seed(MASTER_SEED)
  
  # General model - all 12 populations, 20 replications
  cat("Starting general model - all 12 populations\n")
  out_all12 <- run_hopo_cv(pops = populations, label = "general-all12", n_rep = 20)
  saveRDS(out_all12, "results_all12_v5.rds")
  
  # Population-specific models - 100 iterations
  cat("\nStarting population-specific models\n")
  pop_specific_out     <- run_pop_specific(n_iter = 100)
  results_pop_specific <- pop_specific_out$results
  saveRDS(pop_specific_out, "results_pop_specific_v5.rds")
  
  # Combine results
  multivariate_results <- bind_rows(out_all12$results, results_pop_specific) %>%
    left_join(pop_labels, by = "population") %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))
  
  # Benchmark: populations meeting FDR <= 0.25 AND FNR <= 0.50 under pop-specific PLS-DA
  mv_benchmark <- multivariate_results %>%
    filter(model_config == "population-specific") %>%
    dplyr::select(population, FDR_plsda, FNR_plsda) %>%
    mutate(meets_benchmark = !is.na(FDR_plsda) & !is.na(FNR_plsda) &
             FDR_plsda <= FDR_TARGET & FNR_plsda <= FNR_MAX)
  
  multivariate_results <- multivariate_results %>%
    left_join(mv_benchmark %>% dplyr::select(population, meets_benchmark),
              by = "population")
  
  benchmark_pops     <- mv_benchmark %>% filter(meets_benchmark)  %>% pull(population)
  non_benchmark_pops <- mv_benchmark %>% filter(!meets_benchmark) %>% pull(population)
  
  cat("Populations meeting benchmark:", paste(benchmark_pops, collapse = ", "), "\n")
  cat("Populations NOT meeting benchmark:", paste(non_benchmark_pops, collapse = ", "), "\n")
  
  save(multivariate_results, benchmark_pops, non_benchmark_pops, file = save_file)
  cat("Multivariate results saved to", save_file, "\n")
  
  invisible(list(multivariate_results = multivariate_results,
                 benchmark_pops = benchmark_pops,
                 non_benchmark_pops = non_benchmark_pops))
}


# =============================================================================
# PART 4: DEFECTIVE KERNEL SORTING GENERAL MODEL
# =============================================================================

# --- 4: General Model - Defective kernel HOPO-CV ---
# Binary classification: defective (1) vs non-defective (4)
# Test set: all haploids + 25% hybrids + all defectives from held-out population
run_defective_hopo_cv <- function(n_rep = 20, save_file = "results_defective_sorter_v1.rds") {
  cat("\n--- PART 5: Defective kernel sorting ---\n")
  
  # Populations with defectives (exclude 3 and 9 - zero defectives)
  pops_with_def <- populations[!populations %in% c("3", "9")]
  cat("populations evaluated:", paste(pops_with_def, collapse = ", "), "\n")
  
  all_rep_results <- list()
  set.seed(MASTER_SEED)
  
  for (rep in 1:n_rep) {
    
    local_seed <- MASTER_SEED + rep * 1000
    set.seed(local_seed)
    
    cat("\n  Replication:", rep, "/", n_rep, "\n")
    results <- list()
    
    for (held_out in pops_with_def) {
      train_pops <- pops_with_def[pops_with_def != held_out]
      
      # Training: haploid + hybrid = non-defective, defective = defective
      train_nondef <- dat_train_only %>%
        filter(population %in% train_pops) %>%
        mutate(def_class = 4)   # non-defective
      train_def <- dat_model %>%
        filter(population %in% train_pops, dna.score.final == 5) %>%
        mutate(def_class = 1)   # defective
      train_data <- bind_rows(train_nondef, train_def)
      
      # Test: all haploids + 25% hybrids (75/25 hap/hyb) + all defectives
      hap_test  <- dat_model %>% filter(population == held_out, dna.score.final == 1)
      hyb_test  <- dat_model %>% filter(population == held_out, dna.score.final == 4)
      def_test  <- dat_model %>% filter(population == held_out, dna.score.final == 5)
      
      n_hap     <- nrow(hap_test)
      n_hyb_target <- round(n_hap * 0.25 / 0.75)
      if (nrow(hyb_test) > n_hyb_target)
        hyb_test <- hyb_test %>% slice_sample(n = n_hyb_target)
      
      test_data <- bind_rows(
        hap_test %>% mutate(def_class = 4),
        hyb_test %>% mutate(def_class = 4),
        def_test %>% mutate(def_class = 1)
      )
      
      X_train <- train_data[, spec_cols]
      y_train <- train_data$def_class
      X_test  <- test_data[, spec_cols]
      y_test  <- test_data$def_class
      
      pred_plsda <- tryCatch(run_plsda(X_train, y_train, X_test),
                             error = function(e) rep(NA_integer_, nrow(test_data)))
      svm_out    <- tryCatch(run_svm(X_train, y_train, X_test, seed = local_seed + 100),
                             error = function(e) list(predictions = rep(NA_integer_, nrow(test_data))))
      rf_out     <- tryCatch(run_rf(X_train, y_train, X_test, seed = local_seed + 200),
                             error = function(e) list(predictions = rep(NA_integer_, nrow(test_data))))
      
      m_plsda <- compute_def_metrics(pred_plsda,          y_test)
      m_svm   <- compute_def_metrics(svm_out$predictions, y_test)
      m_rf    <- compute_def_metrics(rf_out$predictions,  y_test)
      
      results[[held_out]] <- data.frame(
        population = held_out, replication = rep,
        FDR_plsda = m_plsda$FDR, FNR_plsda = m_plsda$FNR,
        FDR_svm   = m_svm$FDR,   FNR_svm   = m_svm$FNR,
        FDR_rf    = m_rf$FDR,    FNR_rf    = m_rf$FNR
      )
    }
    all_rep_results[[rep]] <- do.call(rbind, results)
  }
  
  all_reps_df <- do.call(rbind, all_rep_results)
  results_avg <- all_reps_df %>%
    group_by(population) %>%
    summarise(FDR_plsda = round(mean(FDR_plsda, na.rm = TRUE), 3),
              FNR_plsda = round(mean(FNR_plsda, na.rm = TRUE), 3),
              FDR_svm   = round(mean(FDR_svm,   na.rm = TRUE), 3),
              FNR_svm   = round(mean(FNR_svm,   na.rm = TRUE), 3),
              FDR_rf    = round(mean(FDR_rf,    na.rm = TRUE), 3),
              FNR_rf    = round(mean(FNR_rf,    na.rm = TRUE), 3),
              .groups = "drop") %>%
    left_join(pop_labels, by = "population")
  
  cat("\nDefective kernel sorting results:\n"); print(results_avg, row.names = FALSE)
  
  write.csv(results_avg,
          "defective_kernel_results.csv",
          row.names = FALSE)
          
  cat("\nMean across all populations:\n")
  results_avg %>%
    summarise(mean_FDR_plsda = round(mean(FDR_plsda), 3),
              mean_FNR_plsda = round(mean(FNR_plsda), 3),
              mean_FDR_svm   = round(mean(FDR_svm), 3),
              mean_FNR_svm   = round(mean(FNR_svm), 3),
              mean_FDR_rf    = round(mean(FDR_rf), 3),
              mean_FNR_rf    = round(mean(FNR_rf), 3)) %>% print()
  
  saveRDS(list(results = results_avg, all_reps = all_reps_df), save_file)
  cat("Defective sorter results saved to", save_file, "\n")
  invisible(list(results = results_avg, all_reps = all_reps_df))
}


# =============================================================================
# PART 5: Figure 3: Predicted Oil Content Boxplots with Thresholds
# =============================================================================

make_figure3 <- function(oil_results_file = "oil_results_v5.RData") {
  
  cat("\n--- Figure 3: Oil content boxplot with thresholds ---\n")
  
  if (!exists("comparison_table"))
    load(oil_results_file)

  pop_order <- populations
  
  x_labels <- ifelse(
    pop_order %in% sortable,
    paste0(pop_order, "*"),
    pop_order
  )
  
  label_data <- comparison_table %>%
    mutate(
      x_label = paste0(
        x_labels[match(population, pop_order)],
        "<br><span style='font-size:11pt;'><b>", Donor, "</b></span>",
        "<br>", Inducer
      )
    )
  
  x_labels <- setNames(label_data$x_label, label_data$population)
  dat_plot <- dat_model %>%
    filter(
      !is.na(pred.oil),
      dna.score.final %in% c(1,4,5)
    ) %>%
    mutate(
      population = factor(population,
                          levels = pop_order),
      dna.score.final = factor(
        dna.score.final,
        levels = c(1,4,5),
        labels = c(
          "Haploid",
          "Hybrid",
          "Defective"
        )
      )
    )
  
  # Population-specific optimal thresholds
  thresh_opt <- optimal_table %>%
    mutate(
      population = factor(population,
                          levels = pop_order)
    ) %>%
    select(
      population,
      threshold = threshold_optimal
    )
  
  # Lowest percentage threshold
  thresh_pct <- dat_model %>%
    filter(!is.na(pred.oil)) %>%
    group_by(population) %>%
    summarise(
      threshold = quantile(
        pred.oil,
        probs = optimal_pct$pct / 100
      ),
      .groups = "drop"
    ) %>%
    mutate(
      population = factor(population,
                          levels = pop_order)
    )
  
  # Dynamic threshold
  thresh_dyn <- dat_model %>%
    filter(!is.na(pred.oil)) %>%
    group_by(population) %>%
    summarise(
      threshold = mean(pred.oil) -
        optimal_offset$offset,
      .groups = "drop"
    ) %>%
    mutate(
      population = factor(population,
                          levels = pop_order)
    )
  
  fig3 <- ggplot(
    dat_plot,
    aes(
      x = population,
      y = pred.oil,
      fill = dna.score.final,
      colour = dna.score.final
    )
  ) +
    
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.40,
      colour = "black",
      linewidth = 0.35,
      width = 0.60,
      position = position_dodge(width = 0.8)
    ) +
    
    geom_jitter(
      position = position_jitterdodge(
        jitter.width = 0.20,
        dodge.width = 0.80
      ),
      size = 0.40,
      alpha = 0.60
    ) +
    
    # Optimal threshold
    geom_crossbar(
      data = thresh_opt,
      aes(
        x = population,
        y = threshold,
        ymin = threshold,
        ymax = threshold
      ),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.45,
      width = 0.90
    ) +
    
    # Lowest percentage threshold
    geom_crossbar(
      data = thresh_pct,
      aes(
        x = population,
        y = threshold,
        ymin = threshold,
        ymax = threshold
      ),
      inherit.aes = FALSE,
      colour = "grey40",
      linewidth = 0.45,
      width = 0.80
    ) +
    
    # Dynamic threshold
    geom_crossbar(
      data = thresh_dyn,
      aes(
        x = population,
        y = threshold,
        ymin = threshold,
        ymax = threshold
      ),
      inherit.aes = FALSE,
      colour = "red",
      linewidth = 0.45,
      width = 0.80
    ) +
    
    geom_vline(
      xintercept = 5.5,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.6
    ) +
    
    scale_fill_manual(
      values = c(
        Haploid = "#56B4E9",
        Hybrid = "#E69F00",
        Defective = "#999999"
      )
    ) +
    
    scale_colour_manual(
      values = c(
        Haploid = "#56B4E9",
        Hybrid = "#E69F00",
        Defective = "#999999"
      )
    ) +
    
    scale_x_discrete(labels = x_labels) +
    
    labs(
      x = "Population",
      y = "Predicted Oil Content (%f.w.)",
      fill = "Kernel class",
      colour = "Kernel class"
    ) +
    
    theme_bw(base_family = "Arial") +
    
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA
      ),
      axis.text.x = ggtext::element_markdown(
        family = "Arial",
        size = 12,
        colour = "black",
        lineheight = 1.0
      ),
      axis.text.y = element_text(
        size = 13,
        colour = "black"
      ),
      axis.title = element_text(
        size = 15,
        face = "bold"
      ),
      legend.title = element_text(
        size = 13,
        face = "bold"
      ),
      legend.text = element_text(
        size = 11
      ),
      legend.position = "right"
    )
  
  ggsave("fig3_boxplot_thresholds.png", fig3, width = 12, height = 8, dpi = 600)
  
  cat("Figure 3 saved as fig3_boxplot_thresholds.png\n")
  
  return(fig3)
}

# =========================================================
# RUN ALL ANALYSES
# =========================================================


cat("\nRunning oil sorting analysis...\n")
oil_results <- run_oil_analysis()

cat("\nRunning multivariate analysis...\n")
multivariate_results <- run_multivariate_analysis()

cat("\nRunning defective kernel analysis...\n")
defective_results <- run_defective_hopo_cv()

cat("\nFigure3...\n")
Fig3 <- make_figure3() 

cat("\nALL ANALYSES COMPLETED.\n")

writeLines(
  capture.output(sessionInfo()),
  "sessionInfo.txt"
)
