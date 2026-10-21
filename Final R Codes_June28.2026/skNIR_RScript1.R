# =============================================================================
# R SCRIPT 1: skNIR Haploid Sorting Analysis (Two Oil-based Sorting Methods, Three Multivariate Methods and Figure 3 & Supplementary Figure 3)
# Paper: Single-kernel near-infrared spectroscopy enables haploid kernel sorting in field and sweet corn using high-oil haploid inducers across diverse donor–inducer combinations
# =============================================================================
# ORGANIZATION:
#   PART 0: Libraries, data loading, and shared setup
#   PART 1: Multivariate functions for the three methods (PLS-DA, SVM, RF)
#   PART 2: Oil-based haploid sorting methods
#   PART 3: Multivariate haploid sorting methods
#   PART 4: Defective Sorting + End-to-end Metrics
#   PART 5: Figure 3 (based on the results from the PART 2)
#   PART 6: End-to-end 3 x 3 Confusion Matrix (Supplementary Figure 3)
# =============================================================================

# =============================================================================
#   PART 0: Libraries, data loading, and shared setup
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
dat <- read.csv("/home/27June_skNIR_finalcheck/Data_skNIR_Manuscript.csv")

dat$.row_id <- seq_len(nrow(dat))

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
# --- Defective kernel sorting metrics ---
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
stage1_test_rows <- list()

  for (rep in 1:n_rep) {
    set.seed(MASTER_SEED + rep)
    cat("\n  Replication:", rep, "/", n_rep, "\n")
    results <- list()
    
    for (held_out in pops) {
      train_pops <- pops[pops != held_out]
      
      # Training: All haploid and hybrid kernels from 11 non-held-out populations
      train_data <- dat_train_only %>% filter(population %in% train_pops)
      
      # Test: resample held-out to its own HIR, and add defectives
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
      
      pred_plsda <- tryCatch(run_plsda(X_train, y_train, X_test, seed = local_seed),
                             error = function(e) rep(NA_integer_, nrow(test_data)))
      svm_out    <- tryCatch(run_svm(X_train, y_train, X_test, seed = local_seed + 100),
                             error = function(e) list(predictions = rep(NA_integer_, nrow(test_data))))
      rf_out     <- tryCatch(run_rf(X_train, y_train, X_test, seed = local_seed + 200),
                             error = function(e) list(predictions = rep(NA_integer_, nrow(test_data))))
    
  # Save exact Stage-1 test-set rows and predictions for Stage 2
   stage1_test_rows[[paste(rep, held_out, sep = "_")]] <- test_data %>%
   mutate(
    stage1_model = "svm",
    stage1_replication = rep,
    stage1_prediction = svm_out$predictions
  )

      m_plsda <- compute_metrics(pred_plsda,          y_test)
      m_svm   <- compute_metrics(svm_out$predictions, y_test)
      m_rf    <- compute_metrics(rf_out$predictions,  y_test)
      
      results[[held_out]] <- data.frame(
        population  = held_out, replication = rep, model_config = label,
        validation = "HOPO-CV", test_IR = test_ir,
        FDR_plsda = m_plsda$FDR, FNR_plsda = m_plsda$FNR,
        FDR_svm   = m_svm$FDR,   FNR_svm   = m_svm$FNR,
        FDR_rf    = m_rf$FDR,    FNR_rf    = m_rf$FNR,
        
        check.names = FALSE
        
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
          
invisible(list(
  results = results_avg,
  all_reps = all_reps_df,
  stage1_test_rows = bind_rows(stage1_test_rows)
))
}


# --- 3B: Population-specific models ---
# For each of 100 iterations: resample to population HIR, then 70/30 split for train and test sets
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
      
      pred_plsda <- tryCatch(run_plsda(X_train, y_train, X_test, seed = local_seed),
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

run_multivariate_analysis <- function(
    save_file = "multivariate_results_v5.RData"
) {
  
  cat("\n--- PART 3: Multivariate haploid sorting ---\n")
  
  set.seed(MASTER_SEED)
  
  cat("\nStarting general model - all 12 populations\n")
  
  out_all12 <- run_hopo_cv(
    pops = populations,
    label = "general-all12",
    n_rep = 20
  )

  stage1_test_rows <- out_all12$stage1_test_rows
  
  saveRDS(
    out_all12,
    "results_all12_v5.rds"
  )
  
  cat("\nStarting population-specific models\n")
  
  pop_specific_out <- run_pop_specific(
    n_iter = 100
  )
  
  results_pop_specific <- pop_specific_out$results
  
  saveRDS(
    pop_specific_out,
    "results_pop_specific_v5.rds"
  )
  
  cat("\nCombining multivariate results\n")
  
  multivariate_results <- bind_rows(
    out_all12$results,
    results_pop_specific
  ) %>%
    left_join(
      pop_labels,
      by = "population"
    ) %>%
    mutate(
      across(
        where(is.numeric),
        ~ ifelse(is.nan(.), NA, .)
      )
    )
  
  cat("\nCalculating population-specific benchmarks under PLS-DA\n")
  
  mv_benchmark <- multivariate_results %>%
    filter(
      model_config == "population-specific"
    ) %>%
    dplyr::select(
      population,
      FDR_plsda,
      FNR_plsda
    ) %>%
    mutate(
      meets_benchmark =
        !is.na(FDR_plsda) &
        !is.na(FNR_plsda) &
        FDR_plsda <= FDR_TARGET &
        FNR_plsda <= FNR_MAX
    )
  
  multivariate_results <- multivariate_results %>%
    left_join(
      mv_benchmark %>%
        dplyr::select(
          population,
          meets_benchmark
        ),
      by = "population"
    )
  
  benchmark_pops <- mv_benchmark %>%
    filter(meets_benchmark) %>%
    pull(population)
  
  non_benchmark_pops <- mv_benchmark %>%
    filter(!meets_benchmark) %>%
    pull(population)
  
  cat(
    "\nPopulations meeting benchmark:",
    paste(benchmark_pops, collapse = ", "),
    "\n"
  )
  
  cat(
    "Populations NOT meeting benchmark:",
    paste(non_benchmark_pops, collapse = ", "),
    "\n"
  )
  
  save(
    multivariate_results,
    stage1_test_rows,
    benchmark_pops,
    non_benchmark_pops,
    file = save_file
  )
  
  cat(
    "\nMultivariate results saved to",
    save_file,
    "\n"
  )
  
  invisible(
    list(
      multivariate_results = multivariate_results,
      stage1_test_rows = stage1_test_rows,
      benchmark_pops = benchmark_pops,
      non_benchmark_pops = non_benchmark_pops
    )
  )
  
}

# ============================================================
# PART 4: DEFECTIVE SORTING (Stage 2) + END-TO-END METRICS
# ============================================================

run_defective_hopo_cv <- function(
    stage1_results_file = "results_all12_v5.rds",
    stage1_model_name = "svm",
    save_file = "results_defective_sorter_v2.rds"
) {

  stage1_results <- readRDS(stage1_results_file)

  # Exact Stage-1 test rows, including actual class and
  # Stage-1 SVM predictions
  stage1_test_rows <- stage1_results$stage1_test_rows

  cat("\n============================================================\n")
  cat("STAGE 2: DEFECTIVE SORTING + END-TO-END ANALYSIS\n")
  cat("Stage 1:", stage1_model_name, "\n")
  cat("Stage 2: PLS-DA\n")
  cat("Using exact Stage-1 predicted-haploid pool\n")
  cat("All 12 populations included\n")
  cat("============================================================\n")

  stage1_reps <- sort(unique(stage1_test_rows$stage1_replication))

  defective_pops <- populations

  all_results <- list()
  all_end_to_end <- list()

  for (rep in stage1_reps) {

    cat("\nReplication:", rep, "\n")

    for (held_out in defective_pops) {

      cat("  Held-out population:", held_out, "\n")
      
      train_data <- dat_model %>%
        filter(population != held_out) %>%
        mutate(
          def_class = ifelse(
            dna.score.final == 5,
            1L,
            4L
          )
        )

      test_data <- stage1_test_rows %>%
        filter(
          stage1_replication == rep,
          as.character(population) == as.character(held_out),
          stage1_model == stage1_model_name,
          stage1_prediction == 1
        ) %>%
        mutate(
          def_class = ifelse(
            dna.score.final == 5,
            1L,
            4L
          )
        )

      cat(
        "    Stage-1 predicted haploid rows:",
        nrow(test_data),
        "| Defective:",
        sum(test_data$def_class == 1, na.rm = TRUE),
        "| Non-defective:",
        sum(test_data$def_class == 4, na.rm = TRUE),
        "\n"
      )

      if (nrow(test_data) > 0) {

        local_seed <- MASTER_SEED + rep * 1000
        
        X_train <- train_data[, spec_cols]
        y_train <- train_data$def_class

        X_test <- test_data[, spec_cols]
        y_test <- test_data$def_class

        plsda_out <- run_plsda(
          X_train,
          y_train,
          X_test,
          seed = local_seed 
        )
        
        svm_out <- run_svm(
          X_train,
          y_train,
          X_test,
          seed = local_seed + 100
        )
        
        rf_out <- run_rf(
          X_train,
          y_train,
          X_test,
          seed = local_seed + 200
        )

        plsda_metrics <- compute_def_metrics(
          plsda_out,
          y_test
        )

        svm_metrics <- compute_def_metrics(
          svm_out$predictions,
          y_test
        )

        rf_metrics <- compute_def_metrics(
          rf_out$predictions,
          y_test
        )

      } else {

        y_test <- integer(0)
        plsda_out <- integer(0)
        svm_out <- list(predictions = integer(0))
        rf_out <- list(predictions = integer(0))

        plsda_metrics <- data.frame(FDR = 0, FNR = 1)
        svm_metrics   <- data.frame(FDR = 0, FNR = 1)
        rf_metrics    <- data.frame(FDR = 0, FNR = 1)
      }

      all_results[[paste(rep, held_out, sep = "_")]] <- tibble(
        stage1_model = stage1_model_name,
        replication = rep,
        population = held_out,

        n_test = length(y_test),
        n_defective = sum(y_test == 1, na.rm = TRUE),
        n_non_defective = sum(y_test == 4, na.rm = TRUE),

        plsda_FDR = plsda_metrics$FDR,
        plsda_FNR = plsda_metrics$FNR,

        svm_FDR = svm_metrics$FDR,
        svm_FNR = svm_metrics$FNR,

        rf_FDR = rf_metrics$FDR,
        rf_FNR = rf_metrics$FNR
      )

      stage1_all <- stage1_test_rows %>%
        filter(
          stage1_replication == rep,
          as.character(population) == as.character(held_out),
          stage1_model == stage1_model_name
        ) %>%
        mutate(
          stage1_prediction = as.integer(stage1_prediction)
        )

      stage2_predictions <- test_data %>%
        dplyr::select(.row_id) %>%
        mutate(
          stage2_plsda_prediction = as.integer(plsda_out)
        )

      end_to_end_rep <- stage1_all %>%
        left_join(
          stage2_predictions,
          by = ".row_id"
        ) %>%
        mutate(
          final_prediction = case_when(

            # Stage 1 says hybrid -> final hybrid
            stage1_prediction == 4L ~ 4L,

            # Stage 1 says haploid + Stage 2 says non-defective
            # -> final haploid
            stage1_prediction == 1L &
              stage2_plsda_prediction == 4L ~ 1L,

            # Stage 1 says haploid + Stage 2 says defective
            # -> final defective
            stage1_prediction == 1L &
              stage2_plsda_prediction == 1L ~ 5L,

            TRUE ~ NA_integer_
          ),

          actual_class = as.integer(dna.score.final),

          population = held_out,
          replication = rep,
          stage1_model = stage1_model_name,
          stage2_model = "plsda"
        )


      # ------------------------------------------------------
      # 1. Fraction of Stage-1 predicted haploids that are
      #    actually defective
      # ------------------------------------------------------

      n_stage1_pred_haploid <- sum(
        end_to_end_rep$stage1_prediction == 1L,
        na.rm = TRUE
      )

      n_defective_in_stage1_haploid_pool <- sum(
        end_to_end_rep$stage1_prediction == 1L &
        end_to_end_rep$actual_class == 5L,
        na.rm = TRUE
      )

      fraction_defective_stage1_haploid_pool <-
        ifelse(
          n_stage1_pred_haploid > 0,
          n_defective_in_stage1_haploid_pool /
            n_stage1_pred_haploid,
          0
        )

      # ------------------------------------------------------
      # 2. End-to-end haploid purity
      #
      # True haploids / all final predicted haploids
      # ------------------------------------------------------

      n_final_pred_haploid <- sum(
        end_to_end_rep$final_prediction == 1L,
        na.rm = TRUE
      )

      n_true_haploid_final_pool <- sum(
        end_to_end_rep$final_prediction == 1L &
        end_to_end_rep$actual_class == 1L,
        na.rm = TRUE
      )

      end_to_end_purity <-
        ifelse(
          n_final_pred_haploid > 0,
          n_true_haploid_final_pool /
            n_final_pred_haploid,
          0
        )

      # ------------------------------------------------------
      # 3. End-to-end haploid recovery
      #
      # True haploids recovered / all true haploids
      # ------------------------------------------------------

      n_actual_haploid <- sum(
        end_to_end_rep$actual_class == 1L,
        na.rm = TRUE
      )

      end_to_end_recovery <-
        ifelse(
          n_actual_haploid > 0,
          n_true_haploid_final_pool /
            n_actual_haploid,
          0
        )

      end_to_end_rep$population <- held_out
      end_to_end_rep$replication <- rep

      end_to_end_rep$stage1_predicted_haploid_n <-
        n_stage1_pred_haploid

      end_to_end_rep$defective_in_stage1_predicted_haploid_n <-
        n_defective_in_stage1_haploid_pool

      end_to_end_rep$fraction_defective_stage1_predicted_haploid <-
        fraction_defective_stage1_haploid_pool

      end_to_end_rep$final_predicted_haploid_n <-
        n_final_pred_haploid

      end_to_end_rep$true_haploid_in_final_haploid_pool_n <-
        n_true_haploid_final_pool

      end_to_end_rep$actual_haploid_n <-
        n_actual_haploid

      end_to_end_rep$end_to_end_purity <-
        end_to_end_purity

      end_to_end_rep$end_to_end_recovery <-
        end_to_end_recovery

      all_end_to_end[[paste(rep, held_out, sep = "_")]] <-
        end_to_end_rep
    }
  }

  results_df <- bind_rows(all_results) %>%
    left_join(
      pop_labels,
      by = "population"
    )

  write.csv(
    results_df,
    "defective_kernel_results.csv",
    row.names = FALSE
  )

  summary_df <- results_df %>%
    group_by(
      stage1_model,
      population,
      pop_label
    ) %>%
    summarise(
      n_replications = n(),

      mean_FDR_plsda = round(
        mean(plsda_FDR, na.rm = TRUE), 3
      ),
      mean_FNR_plsda = round(
        mean(plsda_FNR, na.rm = TRUE), 3
      ),

      mean_FDR_svm = round(
        mean(svm_FDR, na.rm = TRUE), 3
      ),
      mean_FNR_svm = round(
        mean(svm_FNR, na.rm = TRUE), 3
      ),

      mean_FDR_rf = round(
        mean(rf_FDR, na.rm = TRUE), 3
      ),
      mean_FNR_rf = round(
        mean(rf_FNR, na.rm = TRUE), 3
      ),

      .groups = "drop"
    )

  write.csv(
    summary_df,
    "defective_kernel_summary.csv",
    row.names = FALSE
  )

  end_to_end_predictions_df <- bind_rows(all_end_to_end)

  write.csv(
    end_to_end_predictions_df,
    "end_to_end_predictions_svm_plsda.csv",
    row.names = FALSE
  )

  # ==========================================================
  # END-TO-END METRICS PER REPLICATION
  # ==========================================================

  end_to_end_replications <- end_to_end_predictions_df %>%
    group_by(
      population,
      replication
    ) %>%
    summarise(
      stage1_predicted_haploid_n =
        first(stage1_predicted_haploid_n),

      defective_in_stage1_predicted_haploid_n =
        first(defective_in_stage1_predicted_haploid_n),

      fraction_defective_stage1_predicted_haploid =
        first(
          fraction_defective_stage1_predicted_haploid
        ),

      final_predicted_haploid_n =
        first(final_predicted_haploid_n),

      true_haploid_in_final_haploid_pool_n =
        first(true_haploid_in_final_haploid_pool_n),

      actual_haploid_n =
        first(actual_haploid_n),

      end_to_end_purity =
        first(end_to_end_purity),

      end_to_end_recovery =
        first(end_to_end_recovery),

      .groups = "drop"
    ) %>%
    left_join(
      pop_labels,
      by = "population"
    )

  write.csv(
    end_to_end_replications,
    "end_to_end_metrics_replications_svm_plsda.csv",
    row.names = FALSE
  )

  # ==========================================================
  # END-TO-END METRICS AVERAGED ACROSS 20 REPLICATIONS
  # ==========================================================

  end_to_end_metrics <- end_to_end_replications %>%
    group_by(
      population,
      pop_label
    ) %>%
    summarise(
      n_replications = n(),

      mean_fraction_defective_stage1_predicted_haploid =
        round(
          mean(
            fraction_defective_stage1_predicted_haploid,
            na.rm = TRUE
          ),
          4
        ),

      mean_end_to_end_purity =
        round(
          mean(end_to_end_purity, na.rm = TRUE),
          4
        ),

      mean_end_to_end_recovery =
        round(
          mean(end_to_end_recovery, na.rm = TRUE),
          4
        ),

      .groups = "drop"
    )

  write.csv(
    end_to_end_metrics,
    "end_to_end_metrics_svm_plsda.csv",
    row.names = FALSE
  )

  # ==========================================================
  # 3 x 3 CONFUSION MATRIX
  # ==========================================================

  confusion_replications <- end_to_end_predictions_df %>%
    filter(
      actual_class %in% c(1L, 4L, 5L),
      final_prediction %in% c(1L, 4L, 5L)
    ) %>%
    count(
      population,
      replication,
      actual_class,
      final_prediction,
      name = "count"
    ) %>%
    complete(
      population,
      replication,
      actual_class = c(1L, 4L, 5L),
      final_prediction = c(1L, 4L, 5L),
      fill = list(count = 0)
    ) %>%
    arrange(
      population,
      replication,
      actual_class,
      final_prediction
    )

  write.csv(
    confusion_replications,
    "end_to_end_confusion_replications_svm_plsda.csv",
    row.names = FALSE
  )


  average_confusion_matrices <- confusion_replications %>%
    group_by(
      population,
      actual_class,
      final_prediction
    ) %>%
    summarise(
      mean_count = mean(count),
      .groups = "drop"
    ) %>%
    left_join(
      pop_labels,
      by = "population"
    ) %>%
    arrange(
      population,
      actual_class,
      final_prediction
    )

  write.csv(
    average_confusion_matrices,
    "average_confusion_matrices_svm_plsda.csv",
    row.names = FALSE
  )

  output <- list(
    results = results_df,
    summary = summary_df,

    end_to_end_predictions =
      end_to_end_predictions_df,

    end_to_end_replications =
      end_to_end_replications,

    end_to_end_metrics =
      end_to_end_metrics,

    confusion_replications =
      confusion_replications,

    average_confusion_matrices =
      average_confusion_matrices,

    stage1_model = stage1_model_name,
    stage2_model = "plsda"
  )

  saveRDS(
    output,
    save_file
  )

  cat("\n============================================================\n")
  cat("STAGE-2 + END-TO-END ANALYSIS COMPLETED\n")
  cat("============================================================\n")

  cat("\nEnd-to-end metrics:\n")
  print(end_to_end_metrics, row.names = FALSE)

  cat("\nAverage 3 x 3 confusion matrices:\n")
  print(average_confusion_matrices, row.names = FALSE)

  cat("\nResults saved to:", save_file, "\n")

  return(output)
}


# =============================================================================
# PART 5: Figure 3 (Predicted Oil Content Boxplots with Thresholds)
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

# ============================================================
# PART 6: FINAL END-TO-END 3 x 3 CONFUSION MATRIX FIGURE
# Stage 1 SVM Output + Stage 2 PLS-DA Output
# ============================================================

confusion_plot_data <- defective_results$average_confusion_matrices %>%
  mutate(
    
    actual_label = factor(
      actual_class,
      levels = c(1, 4, 5),
      labels = c("Haploid", "Hybrid", "Defective")
    ),
    
    predicted_label = factor(
      final_prediction,
      levels = c(1, 4, 5),
      labels = c("Haploid", "Hybrid", "Defective")
    ),
    
    cell_type = case_when(
      
      # ========================================================
      # GREEN:
      # All correct classifications
      # ========================================================
      
      actual_class == 1 &
        final_prediction == 1 ~ "Correct Classification",
      
      actual_class == 4 &
        final_prediction == 4 ~ "Correct Classification",
      
      actual_class == 5 &
        final_prediction == 5 ~ "Correct Classification",
      
      # ========================================================
      # RED:
      # Contamination of the final haploid pool
      # ========================================================
      
      actual_class == 4 &
        final_prediction == 1 ~ "Haploid Contamination",
      
      actual_class == 5 &
        final_prediction == 1 ~ "Haploid Contamination",
      
      # ========================================================
      # YELLOW:
      # Loss of true haploids to another pool
      # ========================================================
      
      actual_class == 1 &
        final_prediction == 4 ~ "Haploid Loss",
      
      actual_class == 1 &
        final_prediction == 5 ~ "Haploid Loss",
      
      # ========================================================
      # GRAY:
      # Other misclassification
      # ========================================================
      
      TRUE ~ "Other Misclassification"
    )
  )

confusion_plot_data$pop_label <- pop_labels$pop_label[
  match(
    as.character(confusion_plot_data$population),
    as.character(pop_labels$population)
  )
]

confusion_plot_data$pop_label <- factor(
  confusion_plot_data$pop_label,
  levels = unique(
    pop_labels$pop_label[
      match(
        as.character(populations),
        as.character(pop_labels$population)
      )
    ]
  )
)

p_end_to_end_confusion <- ggplot(
  confusion_plot_data,
  aes(
    x = predicted_label,
    y = actual_label,
    fill = cell_type
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.8
  ) +
  
  geom_text(
    aes(
      label = sprintf("%.2f", mean_count)
    ),
    size = 4.5,
    colour = "black",
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "Correct Classification" = "#66BD63",
      "Haploid Contamination" = "#EF553B",
      "Haploid Loss" = "#FFD966",
      "Other Misclassification" = "#BDBDBD"
    ),
    
    breaks = c(
      "Correct Classification",
      "Haploid Contamination",
      "Haploid Loss",
      "Other Misclassification"
    )
  ) +
  
  scale_x_discrete(
    drop = FALSE
  ) +

  scale_y_discrete(
    drop = FALSE
  ) +
  
  labs(
    x = "Predicted class",
    y = "Actual class",
    fill = NULL
  ) +
  
  facet_wrap(
    ~ pop_label,
    ncol = 3
  ) +
  
  theme_minimal(
    base_size = 14
  ) +

  theme(
    
    panel.grid = element_blank(),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    axis.text.x = element_text(
      angle = 0,
      hjust = 1,
      face = "bold",
      size = 13
    ),
    
    axis.text.y = element_text(
      face = "bold",
      size = 13
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 15
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    legend.text = element_text(
      face = "bold",
      size = 12
    ),
    
    legend.position = "bottom"
  )


print(p_end_to_end_confusion)

ggsave(
  "end_to_end_confusion_matrix_svm_plsda.png",
  p_end_to_end_confusion,
  width = 14,
  height = 16,
  units = "in",
  dpi = 600,
  bg = "white"
)


cat(
  "\nFinal end-to-end 3 x 3 confusion matrix figure saved to:\n",
  "  end_to_end_confusion_matrix_svm_plsda.png\n"
)

cat("\nALL ANALYSES COMPLETED.\n")


writeLines(
  capture.output(sessionInfo()),
  "sessionInfo.txt"
)
