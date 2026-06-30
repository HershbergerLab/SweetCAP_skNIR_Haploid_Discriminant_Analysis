# =================================================================================================================================================================================
# R SCRIPT 2
# Paper: Single-kernel near-infrared spectroscopy enables haploid kernel sorting in field and sweet corn using high-oil haploid inducers across diverse donor–inducer combinations
# =================================================================================================================================================================================
# ORGANIZATION:

# PART 1: Analysis for ANOVA and LS Means for Delta Oil and Haploid Induction Rate (HIR), Kernel Counts, Visual Classification Accuracy and Other Analyses 
       # 1A. Delta Oil and Haploid Induction Rate (HIR) ANOVA + LS Means (for Supplementary Table 1 and Table 2)
       # 1B. Kernel Counts (Table 1)
       # 1C: Visual Classification Accuracy (Table 2)
       # 1D. Shapiro-Wilk Test, Kruskal-Wallis and Dunn's Post-hoc Test for the Predicted Oil Content (for Manual Addition of Significant Difference Letters to the Figure 3)
       # 1E: Correlation of Delta Oil and Haploid Induction Rate (HIR)
       # 1F: NMR-measured Oil Content and skNIR-predicted Oil Content (R2 and RMSE) (for Supplementary Fig. 2)

# PART 2: Figures (Fig. 2, Supplementary Fig. 1, Fig. 4, Fig. 5, and Supplementary Fig. 2)
# ===============================================================================================================================

library(dplyr)
library(tidyr)
library(emmeans)
library(car)
library(FSA)
library(multcompView)
library(rstatix)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggtext)

# Upload data
dat <- read.csv("/home/sumans/26June_skNIR_finalcheck/Data_skNIR_Manuscript_June27.2026.csv")

# Keeping only kernels with final classification scores used in analysis
dat_model <- dat %>%
  filter(!is.na(dna.score.final)) %>%
  filter(dna.score.final %in% c(1, 4, 5))

# Labels for dna.score.final = 1 (haploid), 4 (hybrid) and 5 (defective)
seed_class <- dat_model %>%
  mutate(
    dna.score.final = factor(
      dna.score.final,
      levels = c(1, 4, 5),
      labels = c("Haploid", "Hybrid", "Defective")
    )
  )

# Population-specific Haploid Induction Rates (HIRs) (Corrected by gDNA scoring)
induction_rates <- data.frame(
  population = 1:12,
  IR = c(0.08, 0.13, 0.10, 0.11, 0.10, 0.04, 0.03, 0.08, 0.08, 0.09, 0.08, 0.12),
  stringsAsFactors = FALSE)

# Spectral columns (940-1640 nm, 1 nm step)
spec_cols <- paste0("X", 940:1640)

#---------------------------------------------------------------------------------
# PART 1A: Delta Oil and HIR ANOVA + LS Means (for Supplementary Table 1 and Table 2)
# Delta oil for each population = Mean Hybrid Oil - Mean Haploid Oil
#---------------------------------------------------------------------------------
run_oil_anova <- function() {
  cat("\n--- PART 1A: Oil content ANOVA and LS means ---\n")
  # ---  Compute delta oil per population ---
  cat("\n-- Delta oil computation --\n")
  delta_oil <- dat_model %>%
    filter(dna.score.final %in% c(1, 4)) %>%
    group_by(population, Donor, Inducer, dna.score.final) %>%
    summarise(mean_oil = mean(pred.oil, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = dna.score.final,
      values_from = mean_oil,
      names_prefix = "class_"
    ) %>%
    mutate(delta_oil = class_4 - class_1)
  cat("Delta oil per population:\n")
  print(delta_oil %>% dplyr::select(population, Donor, Inducer,
                                    class_1, class_4, delta_oil),
        row.names = FALSE)
  
  # --- Two-way ANOVA on delta oil ---
  cat("\n--Two-way ANOVA on delta oil (Donor + Inducer) --\n")
  anova_delta <- lm(delta_oil ~ Donor + Inducer, data = delta_oil)
  anova_delta_typeIII <- Anova(anova_delta, type = "III")
  cat("ANOVA Type III SS for delta oil:\n")
  print(anova_delta_typeIII)
  
  # Eta-squared
  ss_donor   <- anova_delta_typeIII["Donor",   "Sum Sq"]
  ss_inducer <- anova_delta_typeIII["Inducer",  "Sum Sq"]
  ss_resid   <- anova_delta_typeIII["Residuals","Sum Sq"]
  ss_total   <- ss_donor + ss_inducer + ss_resid
  
  supp_table1_deltaoil <- data.frame(
    Source   = c("Donor", "Inducer", "Residual"),
    SS       = round(c(ss_donor, ss_inducer, ss_resid), 4),
    df       = c(anova_delta_typeIII["Donor",    "Df"],
                 anova_delta_typeIII["Inducer",   "Df"],
                 anova_delta_typeIII["Residuals", "Df"]),
    F        = round(c(anova_delta_typeIII["Donor",    "F value"],
                       anova_delta_typeIII["Inducer",   "F value"],
                       NA), 3),
    p_value  = round(c(anova_delta_typeIII["Donor",    "Pr(>F)"],
                       anova_delta_typeIII["Inducer",   "Pr(>F)"],
                       NA), 3),
    eta2     = round(c(ss_donor, ss_inducer, ss_resid) / ss_total, 3),
    eta2_pct = paste0(round(c(ss_donor, ss_inducer, ss_resid) / ss_total * 100, 1), "%")
  )
  cat("\nSupplementary Table 1 - Delta oil ANOVA:\n")
  print(supp_table1_deltaoil, row.names = FALSE)
  
  # ---Two-way ANOVA on induction rate (HIR) ---
  cat("\n--Two-way ANOVA on induction rate (HIR) --\n")
  ir_data <- induction_rates %>%
    left_join(
      dat %>%
        dplyr::select(population, Donor, Inducer) %>%
        distinct(),
      by = "population"
    )
  
  anova_ir <- lm(IR ~ Donor + Inducer, data = ir_data)
  anova_ir_typeIII <- Anova(anova_ir, type = "III")
  cat("ANOVA Type III SS for IR:\n")
  print(anova_ir_typeIII)
  
  ss_donor_ir   <- anova_ir_typeIII["Donor",    "Sum Sq"]
  ss_inducer_ir <- anova_ir_typeIII["Inducer",   "Sum Sq"]
  ss_resid_ir   <- anova_ir_typeIII["Residuals", "Sum Sq"]
  ss_total_ir   <- ss_donor_ir + ss_inducer_ir + ss_resid_ir
  
  supp_table1_ir <- data.frame(
    Source   = c("Donor", "Inducer", "Residual"),
    SS       = round(c(ss_donor_ir, ss_inducer_ir, ss_resid_ir), 6),
    df       = c(anova_ir_typeIII["Donor",    "Df"],
                 anova_ir_typeIII["Inducer",   "Df"],
                 anova_ir_typeIII["Residuals", "Df"]),
    F        = round(c(anova_ir_typeIII["Donor",    "F value"],
                       anova_ir_typeIII["Inducer",   "F value"],
                       NA), 3),
    p_value  = round(c(anova_ir_typeIII["Donor",    "Pr(>F)"],
                       anova_ir_typeIII["Inducer",   "Pr(>F)"],
                       NA), 3),
    eta2     = round(c(ss_donor_ir, ss_inducer_ir, ss_resid_ir) / ss_total_ir, 3),
    eta2_pct = paste0(round(c(ss_donor_ir, ss_inducer_ir, ss_resid_ir) / ss_total_ir * 100, 1), "%")
  )
  cat("\nSupplementary Table 1 - HIR ANOVA:\n")
  print(supp_table1_ir, row.names = FALSE)
  
  # --- LS means for delta oil ---
  cat("\n-- LS means for delta oil --\n")
  emm_donor_delta   <- emmeans(anova_delta, ~ Donor)
  emm_inducer_delta <- emmeans(anova_delta, ~ Inducer)
  
  emm_donor_delta_df   <- as.data.frame(summary(emm_donor_delta)) %>%
    rename(LS_mean = emmean, SE = SE, lower_CI = lower.CL, upper_CI = upper.CL) %>%
    mutate(trait = "delta_oil")
  emm_inducer_delta_df <- as.data.frame(summary(emm_inducer_delta)) %>%
    rename(LS_mean = emmean, SE = SE, lower_CI = lower.CL, upper_CI = upper.CL) %>%
    mutate(trait = "delta_oil")
  
  # Bonferroni pairwise comparisons
  pairs_donor_delta   <- pairs(emm_donor_delta,   adjust = "bonferroni")
  pairs_inducer_delta <- pairs(emm_inducer_delta, adjust = "bonferroni")
  
  cat("LS means - Donor (delta oil):\n")
  print(emm_donor_delta_df, row.names = FALSE)
  cat("\nLS means - Inducer (delta oil):\n")
  print(emm_inducer_delta_df, row.names = FALSE)
  cat("\nPairwise comparisons - Donor (delta oil):\n")
  print(pairs_donor_delta)
  cat("\nPairwise comparisons - Inducer (delta oil):\n")
  print(pairs_inducer_delta)
  
  # --- LS means for HIR ---
  cat("\n-- LS means for HIR --\n")
  emm_donor_ir   <- emmeans(anova_ir, ~ Donor)
  emm_inducer_ir <- emmeans(anova_ir, ~ Inducer)
  
  emm_donor_ir_df   <- as.data.frame(summary(emm_donor_ir)) %>%
    rename(LS_mean = emmean, SE = SE, lower_CI = lower.CL, upper_CI = upper.CL) %>%
    mutate(trait = "IR")
  emm_inducer_ir_df <- as.data.frame(summary(emm_inducer_ir)) %>%
    rename(LS_mean = emmean, SE = SE, lower_CI = lower.CL, upper_CI = upper.CL) %>%
    mutate(trait = "IR")
  
  pairs_donor_ir   <- pairs(emm_donor_ir,   adjust = "bonferroni")
  pairs_inducer_ir <- pairs(emm_inducer_ir, adjust = "bonferroni")
  
  cat("LS means - Donor (HIR):\n")
  print(emm_donor_ir_df, row.names = FALSE)
  cat("\nLS means - Inducer (HIR):\n")
  print(emm_inducer_ir_df, row.names = FALSE)
  cat("\nPairwise comparisons - Donor (HIR):\n")
  print(pairs_donor_ir)
  cat("\nPairwise comparisons - Inducer (HIR):\n")
  print(pairs_inducer_ir)
  
  # --- Compile Supplementary Table 2 ---
  cat("\n-- Supplementary Table 2: LS means for delta oil and HIR --\n")
  
  # Compact letter display for grouping
  cld_donor_delta   <- multcomp::cld(emm_donor_delta,   adjust = "bonferroni", Letters = letters)
  cld_inducer_delta <- multcomp::cld(emm_inducer_delta, adjust = "bonferroni", Letters = letters)
  cld_donor_ir      <- multcomp::cld(emm_donor_ir,      adjust = "bonferroni", Letters = letters)
  cld_inducer_ir    <- multcomp::cld(emm_inducer_ir,    adjust = "bonferroni", Letters = letters)
  
  supp_table2 <- list(
    donor_delta   = as.data.frame(cld_donor_delta),
    inducer_delta = as.data.frame(cld_inducer_delta),
    donor_ir      = as.data.frame(cld_donor_ir),
    inducer_ir    = as.data.frame(cld_inducer_ir)
  )
  cat("Supplementary Table 2 data compiled (see supp_table2 list).\n")
  
  print(supp_table2, row.names = FALSE)
  
  # Save Supplementary Table 1
  supp_table1_all <- dplyr::bind_rows(
    supp_table1_deltaoil %>% mutate(Trait = "Delta oil"),
    supp_table1_ir %>% mutate(Trait = "HIR")
  )
  
  write.csv(
    supp_table1_all,
    "Supplementary_Table1.csv",
    row.names = FALSE
  )
  
  # Save Supplementary Table 2
  supp_table2_all <- dplyr::bind_rows(
    supp_table2$donor_delta %>%
      mutate(Factor = "Donor", Trait = "Delta oil"),
    
    supp_table2$inducer_delta %>%
      mutate(Factor = "Inducer", Trait = "Delta oil"),
    
    supp_table2$donor_ir %>%
      mutate(Factor = "Donor", Trait = "HIR"),
    
    supp_table2$inducer_ir %>%
      mutate(Factor = "Inducer", Trait = "HIR")
  )
  
  write.csv(
    supp_table2_all,
    "Supplementary_Table2.csv",
    row.names = FALSE
  )
  
  cat("\nSaved:\n")
  cat("  - Supplementary_Table1.csv\n")
  cat("  - Supplementary_Table2.csv\n")
  
  invisible(list(
    delta_oil           = delta_oil,
    anova_delta         = anova_delta_typeIII,
    anova_ir            = anova_ir_typeIII,
    supp_table1_deltaoil = supp_table1_deltaoil,
    supp_table1_ir      = supp_table1_ir,
    supp_table2         = supp_table2,
    emm_donor_delta     = emm_donor_delta_df,
    emm_inducer_delta   = emm_inducer_delta_df,
    emm_donor_ir        = emm_donor_ir_df,
    emm_inducer_ir      = emm_inducer_ir_df
  ))
}

oil_anova_results <- run_oil_anova()

#---------------------------------------------------
# PART 1B: Kernel Counts (Table 1)
#---------------------------------------------------
kernel_counts <- seed_class %>%
  group_by(population, dna.score.final) %>%
  summarise(
    kernel_count = n(),
    .groups = "drop"
  )

kernel_counts

write.csv(kernel_counts, "kernel_counts_Table1.csv", row.names = FALSE)

#---------------------------------------------------
# PART 1C: Visual Classification Accuracy (Table 2)
#---------------------------------------------------
data_clean <- seed_class %>%
  filter(dna.score.final %in% c("Haploid", "Hybrid")) %>%
  mutate(
    prediction = case_when(
      vis.score == 1 ~ "Haploid",
      vis.score == 4 ~ "Hybrid",
      TRUE ~ NA_character_
    )
  )

# False Discovery Rate (FDR) and False Negative Rate (FNR) by population
results_by_population <- data_clean %>%
  group_by(population) %>%
  summarise(
    TP = sum(prediction == "Haploid" & dna.score.final == "Haploid", na.rm = TRUE),
    FP = sum(prediction == "Haploid" & dna.score.final == "Hybrid", na.rm = TRUE),
    FN = sum(prediction == "Hybrid" & dna.score.final == "Haploid", na.rm = TRUE),
    FDR = round(ifelse((TP + FP) > 0, FP / (TP + FP), 0), 3),
    FNR = round(ifelse((TP + FN) > 0, FN / (TP + FN), 1), 3),
    n = n(),
    .groups = "drop"
  )

results_by_population

write.csv(results_by_population, "visual_scoring_accuracy_Table2.csv", row.names = FALSE)

# Mean FDR and FNR across 12 populations
avg_fdr_fnr <- results_by_population %>%
  summarise(
    avg_FDR = mean(FDR, na.rm = TRUE),
    avg_FNR = mean(FNR, na.rm = TRUE)
  )

avg_fdr_fnr

#------------------------------------------------------------------------------------------------------------------------------
# PART 1D: # Shapiro-Wilk Test (Normality Check), 
           # Kruskal-Wallis and Dunn's Post-hoc Test for the Predicted Oil Content (For Manually Adding the Significant Difference Letters to the Figure 3) 
#------------------------------------------------------------------------------------------------------------------------------
# Step 1: Shapiro-Wilk Test (Normality Check)
normality_by_class <- seed_class %>%
  group_by(dna.score.final) %>%
  shapiro_test(pred.oil)

# Predicted oil content data were Non-normal
normality_by_class

write.csv(normality_by_class, "normality_by_class_predicted_oil.csv", row.names = FALSE)

# Step 2: Kruskal-Wallis and Dunn's Post-hoc Test
kw_results <- data.frame()
all_letters <- data.frame()

for (pop in sort(unique(seed_class$population))) {
  
  temp_data <- filter(seed_class, population == pop)
  
  # Kruskal-Wallis test
  kw <- kruskal.test(pred.oil ~ dna.score.final, data = temp_data)
  
  kw_results <- rbind(
    kw_results,
    data.frame(
      population = pop,
      chi_sq = as.numeric(kw$statistic),
      df = as.numeric(kw$parameter),
      p_value = kw$p.value
    )
  )
  
  # Dunn's test
  if (kw$p.value < 0.05) {
    
    dunn <- dunnTest(
      pred.oil ~ dna.score.final,
      data = temp_data,
      method = "bonferroni"
    )
    
    pvals <- setNames(
      dunn$res$P.adj,
      gsub(" - ", "-", dunn$res$Comparison)
    )
    
    letters <- multcompLetters(pvals)
    
    letter_df <- data.frame(
      population = pop,
      dna.score.final = names(letters$Letters),
      letters = letters$Letters
    )
    
  } else {
    
    letter_df <- data.frame(
      population = pop,
      dna.score.final = c("Haploid", "Hybrid", "Defective"),
      letters = "a"
    )
  }
  
  all_letters <- rbind(all_letters, letter_df)
}

kw_results
all_letters

write.csv(kw_results, "Kruskal_Wallis_Results.csv", row.names = FALSE)
write.csv(all_letters, "Dunn's_Test_Results.csv", row.names = FALSE)


#---------------------------------------------------------------------
# PART 1E: Correlation of Delta Oil and Haploid Induction Rate (HIR)
#---------------------------------------------------------------------

delta_oil <- oil_anova_results$delta_oil

corr_table <- delta_oil %>%
  left_join(
    induction_rates %>%
      dplyr::select(population, IR),
    by = "population"
)

corr_table

pearson_result <- cor.test(
  corr_table$IR,
  corr_table$delta_oil,
  method = "pearson"
)

pearson_result

#------------------------------------------------------------------------------------------------------------
# PART 1F: NMR-measured Oil Content and skNIR-predicted Oil Content (R2 and RMSE) (for Supplementary Fig. 2)
#------------------------------------------------------------------------------------------------------------
stats_table <- seed_class %>%
  filter(!is.na(nmr.oil)) %>%
  group_by(dna.score.final) %>%
  summarise(
    r2 = summary(lm(nmr.oil ~ pred.oil ))$r.squared,
    rmse = sqrt(mean((nmr.oil - pred.oil)^2)),
    .groups = "drop"
  )

stats_table

# =============================================================================
# PART 2: All Figures for the skNIR Manuscript except Figure 3 (Figure 3 is in R Script1)
# =============================================================================

theme_set(
  theme_bw(base_family = "Arial")
)

#---------------------------------------------------
# Figure 2 and Supplementary Figure 1 (Mean Spectra)
#---------------------------------------------------
dat_spec <- dat_model %>%
  mutate(
    dna.score.final = factor(
      dna.score.final,
      levels = c(1, 4, 5),
      labels = c("Haploid", "Hybrid", "Defective")
    )
  )

# Panel A: All populations combined
mean_spectra <- dat_spec %>%
  dplyr::select(dna.score.final, all_of(spec_cols)) %>%
  group_by(dna.score.final) %>%
  summarise(
    across(everything(), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(spec_cols),
    names_to = "wavelength",
    values_to = "absorbance"
  ) %>%
  mutate(
    wavelength = as.numeric(sub("X", "", wavelength))
  )

# Panel B: Population 1 only
mean_spectra_A <- dat_spec %>%
  filter(population == "1") %>%
  dplyr::select(dna.score.final, all_of(spec_cols)) %>%
  group_by(dna.score.final) %>%
  summarise(
    across(everything(), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(spec_cols),
    names_to = "wavelength",
    values_to = "absorbance"
  ) %>%
  mutate(
    wavelength = as.numeric(sub("X", "", wavelength))
  )

plot_dat <- bind_rows(
  mean_spectra %>%
    mutate(facet_label = "All Populations Combined"),
  
  mean_spectra_A %>%
    mutate(facet_label = "Viking 51 × IHO-1 (1)")
)

plot_dat$facet_label <- factor(
  plot_dat$facet_label,
  levels = c(
    "All Populations Combined",
    "Viking 51 × IHO-1 (1)"
  )
)

fig2 <- ggplot(
  plot_dat,
  aes(
    x = wavelength,
    y = absorbance,
    color = dna.score.final,
    linetype = dna.score.final
  )
) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~facet_label, ncol = 2) +
  scale_color_manual(
    values = c(
      "Haploid" = "black",
      "Hybrid" = "darkgrey",
      "Defective" = "red"
    ),
    name = "Kernel class"
  ) +
  scale_linetype_manual(
    values = c(
      "Haploid" = "solid",
      "Hybrid" = "dashed",
      "Defective" = "dotted"
    ),
    name = "Kernel class"
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "Absorbance (log(1/R))"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.position = "bottom"
  )

fig2

ggsave("Figure2.png", fig2, width = 12, height = 5, dpi = 600, bg = "white")


# Supplementary Fig. 1: Per-population MEAN spectra
mean_pop_all <- dat_spec %>%
  dplyr::select(population, dna.score.final, all_of(spec_cols)) %>%
  group_by(population, dna.score.final) %>%
  summarise(
    across(everything(), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(spec_cols),
    names_to = "wavelength",
    values_to = "absorbance"
  ) %>%
  mutate(
    wavelength = as.numeric(sub("X", "", wavelength)),
    facet_label = factor(population, levels = 1:12)
  )

supp_fig1 <- ggplot(
  mean_pop_all,
  aes(
    x = wavelength,
    y = absorbance,
    color = dna.score.final,
    linetype = dna.score.final
  )
) +
  geom_line(linewidth = 0.7) +
  
  scale_color_manual(
    values = c(
      "Haploid" = "black",
      "Hybrid" = "darkgrey",
      "Defective" = "red"
    ),
    name = "Kernel class"
  ) +
  
  scale_linetype_manual(
    values = c(
      "Haploid" = "solid",
      "Hybrid" = "dashed",
      "Defective" = "dotted"
    ),
    name = "Kernel class"
  ) +
  
  facet_wrap(~facet_label, ncol = 4) +
  
  labs(
    x = "Wavelength (nm)",
    y = "Absorbance (log(1/R))"
  ) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.position = "bottom"
  )


supp_fig1

ggsave("SupplementaryFig1.png", supp_fig1, width = 16, height = 14, dpi = 600)

#---------------------------------------------------
# Figure 4: Delta Oil and Inducer Effect
#---------------------------------------------------
# Use the already calculated delta oil
delta_oil_data <- delta_oil %>%
  mutate(
    pop_num = factor(population, levels = 1:12),
    color_group = case_when(
      population %in% c(1, 6, 9)  ~ "G1",
      population %in% c(2, 5, 10) ~ "G2",
      population %in% c(3, 7, 11) ~ "G3",
      population %in% c(4, 8, 12) ~ "G4"
    )
  )

delta_oil_data <- delta_oil_data %>%
  mutate(
    x_label = paste0(
      population,
      "<br><b>", Inducer, "</b>",
      "<br><span style='font-size:11pt;'>", Donor, "</span>"
    )
  )

Inducer_plot <- ggplot(
  delta_oil_data,
  aes(
    x = pop_num,
    y = delta_oil,
    fill = color_group
  )
) +
  geom_col(
    color = "black",
    width = 0.90,
    linewidth = 0.7
  ) +
  
  scale_fill_manual(
    values = c(
      G1 = "#A6CEE3",
      G2 = "#FDBF6F",
      G3 = "#B2DF8A",
      G4 = "#CAB2D6"
    ),
    guide = "none"
  ) +
  
  scale_x_discrete(
    labels = setNames(delta_oil_data$x_label, delta_oil_data$pop_num)
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(
    x = "Population",
    y = expression(bold(Delta~Oil~"(%)"))
  ) +
  
  theme_classic(base_family = "Arial") +
  
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    ),
    
    axis.line = element_blank(),
    
    axis.text.x = ggtext::element_markdown(
      size = 12,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 14,
      colour = "black"
    ),
    
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.title.x = element_text(
      margin = margin(t = 12)
    ),
    
    axis.title.y = element_text(
      margin = margin(r = 12)
    ),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.8
    ),
    
    axis.ticks.length = unit(0.2, "cm"),
    
    plot.margin = margin(10, 10, 10, 10)
  )

Inducer_plot

ggsave("Figure4.png", plot = Inducer_plot, dpi = 600, width = 10, height = 8, units = "in")

#-----------------------------------------------------------------------------
# Figure 5: Correlation of Delta Oil with FDR and FNR from General SVM Model
#-----------------------------------------------------------------------------
# Adding FDR and FNR values manually from Table 5 (General SVM Multivariate Model)
manual_data <- data.frame(
  
  population = 1:12,
  
  FDR = c(0.15, 0.18, 0.02, 0.43, 0.48, 0.50, 0.00, 0.11, 0.64, 0.77, 0.55, 0.02),
  
  FNR = c(0.11, 0.39, 0.44, 0.21, 0.14, 0.37, 0.21, 0.67, 0.67, 0.30, 0.40, 0.69)
)

fig5_data <- delta_oil %>%
  left_join(manual_data, by = "population")

r_fdr <- cor(fig5_data$delta_oil, fig5_data$FDR)

p_fdr <- cor.test(
  fig5_data$delta_oil,
  fig5_data$FDR
)$p.value

r_fnr <- cor(fig5_data$delta_oil, fig5_data$FNR)

p_fnr <- cor.test(
  fig5_data$delta_oil,
  fig5_data$FNR
)$p.value


# FDR Plot
p1 <- ggplot(
  fig5_data,
  aes(
    x = delta_oil,
    y = FDR
  )
) +
  
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "lightgrey"
  ) +
  
  geom_point(
    size = 5,
    color = "darkorange"
  ) +

geom_text_repel(
  aes(label = population),
  size = 5,
  box.padding = 0.4,
  point.padding = 0.3,
  max.overlaps = Inf
) +
  
  geom_hline(
    yintercept = 0.25,
    linetype = "dashed",
    color = "black"
  ) +
  
  annotate(
    "text",
    x = max(fig5_data$delta_oil),
    y = 0.95,
    label = paste0(
      "r = ", round(r_fdr, 2),
      "\np = ", round(p_fdr, 3)
    ),
    hjust = 1,
    size = 4
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  
  labs(
    x = NULL,
    y = "FDR"
  ) +
  
  theme_classic(base_family = "Arial") +
  
  theme(
    legend.position = "none",
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    ),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    axis.title = element_text(
      size = 14,
      face = "bold"
    )
  )

# FNR Plot
p2 <- ggplot(
  fig5_data,
  aes(
    x = delta_oil,
    y = FNR
  )
) +
  
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "lightgrey"
  ) +
  
  geom_point(
    size = 5,
    color = "darkorange"
  ) +

geom_text_repel(
  aes(label = population),
  size = 5,
  box.padding = 0.4,
  point.padding = 0.3,
  max.overlaps = Inf
) +
  geom_hline(
    yintercept = 0.50,
    linetype = "dashed",
    color = "black"
  ) +
  
  annotate(
    "text",
    x = max(fig5_data$delta_oil),
    y = 0.95,
    label = paste0(
      "r = ", round(r_fnr, 2),
      "\np = ", round(p_fnr, 3)
    ),
    hjust = 1,
    size = 4
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  
  labs(
    x = expression(bold(Delta~Oil~"(%)")),
    y = "FNR"
  ) +
  
  theme_classic(base_family = "Arial") +
  
  theme(
    legend.position = "none",
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    ),
    axis.line = element_blank(),
    axis.text.x = element_text(
      size = 12,
      color = "black"
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    axis.title = element_text(
      size = 14,
      face = "bold"
    )
  )

# Combine plots using cowplot
fig5 <- plot_grid(p1, p2, ncol = 1, align = "v")

fig5

ggsave("Figure5.png", plot = fig5, dpi = 600, width = 10, height = 8, units = "in")

# --------------------------------------------------
# Supplementary Figure 2 (NMR vs. skNIR Oil Content)
# --------------------------------------------------
# Create annotation labels
stats_table <- stats_table %>%
  mutate(
    label = paste0(
      dna.score.final,
      ":  R² = ", round(r2, 2),
      "\nRMSE = ", round(rmse, 2)
    )
  )

seed_class <- seed_class %>% filter(!is.na(nmr.oil))

x_min <- min(seed_class$nmr.oil, na.rm = TRUE)
y_max <- max(seed_class$pred.oil, na.rm = TRUE)
y_step <- diff(range(seed_class$pred.oil, na.rm = TRUE)) * 0.08

stats_table <- stats_table %>%
  mutate(
    x = x_min,
    y = y_max - (row_number() - 1) * y_step
  )

oil_content_r2 <- ggplot(
  seed_class,
  aes(
    x = nmr.oil,
    y = pred.oil,
    color = dna.score.final
  )
) +
  geom_point(size = 3, alpha = 0.8) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1.2
  ) +
  
  scale_color_manual(values = c(
    "Haploid" = "#56B4E9",
    "Hybrid" = "#E69F00",
    "Defective" = "#999999"
  )) +
  
  scale_x_continuous(
    limits = c(0, 10),
    breaks = seq(2, 10, 2),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 10),
    breaks = seq(2, 10, 2),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "NMR-measured Oil Content (%f.w.)",
    y = "skNIR-predicted Oil Content (%f.w.)",
    color = "Kernel Class"
  ) +
  
  geom_text(
    data = stats_table,
    aes(
      x = x,
      y = y,
      label = label,
      color = dna.score.final
    ),
    hjust = 0,
    size = 5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  theme_classic(base_family = "Arial") +
  
  theme(
    axis.title = element_text(
      size = 20,
      face = "bold",
      color = "black"
    ),
    
    axis.text = element_text(
      size = 16,
      color = "black"
    ),
    
    legend.title = element_text(
      size = 18,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 16
    ),
    
    legend.position = "right",
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    
    axis.line = element_blank(),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.8
    ),
    
    axis.ticks.length = unit(0.25, "cm"),
    
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    plot.background = element_rect(
      fill = "white",
      color = NA
    )
  )

oil_content_r2

ggsave("SupplementaryFig2.png", plot = oil_content_r2, dpi = 600, width = 10, height = 8, units = "in")

writeLines(
  capture.output(sessionInfo()),
  "sessionInfo.txt"
)
