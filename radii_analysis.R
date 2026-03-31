# ============================================================
#  Bony Endplate Radii — Statistical Analysis
#  Design: Zone (central vs. peripheral) x Vertebra x Cow
#  Model: Linear Mixed-Effects (LMM)
#  Zone = fixed effect; Cow = random effect; Vertebra nested in Cow
# ============================================================

# ── 0. Load packages ───────────────────────────────
library(readxl)
library(lme4)
library(lmerTest)   # p-values
library(emmeans)    # post-hoc comparisons
library(performance)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)

# ── 1. Import Excel file and define your dependant variable ───────────────────────────────
setwd("C:/Users/Ahmad Alminnawi/Desktop/pipeline_publication/modeling-part/submitted_version/second submission/media-guidlines-rev/media/prism/pore_diameter")
df <- read_excel("Radii_R.xlsx", sheet = "Sheet1")

dep_variable <- "Radii"

# Ensure correct data types
df$Cow <- as.factor(df$Cow)
df$Vertebra <- as.factor(df$Vertebra)
df$Zone <- relevel(factor(df$Zone), ref = "Central")

# ── 2. Fit nested mixed-effects model ───────────────────────────────
model <- lmer(Value ~ Zone + (1 | Cow/Vertebra), data = df, REML = TRUE)
m0 <- lmer(Value ~ Zone + (1 | Cow), data = df, REML = TRUE) # take only cow as random effect but not vertebra

# ── 3. Model summary and ANOVA ───────────────────────────────
summary(model)       # main coefficients and random effects
anova(model)  # Type III ANOVA using Satterthwaite's method
ranef(model)         # random effects estimates per group (reflect the vertebra-level variability within a cow and cow-level within all cows)
VarCorr(model)       # variance components

anova(m0, model) # compare model to m0 to see if nested vertebra random effect significantly improves the model fit 



# Make sure the folder "analysis" exists
dir.create("analysis", showWarnings = FALSE)
writeLines(capture.output(summary(model)), "analysis/1-model_summary.txt")
writeLines(capture.output(anova(model)), "analysis/2-model_anova.txt")
writeLines(capture.output(ranef(model)), "analysis/3-model_random_effects.txt")
writeLines(capture.output(VarCorr(model)), "analysis/4-model_variance_components.txt")
writeLines(capture.output(anova(m0, model)), "analysis/5-model-vs-m0_anova.txt")



# ── 4. Estimated marginal means and contrasts ───────────────────────────────
emm <- emmeans(model, ~ Zone)
emm_df <- as.data.frame(emm)
write.csv(emm_df, "analysis/emm_df.csv", row.names = FALSE)


# Pairwise contrasts with 95% CI
contrasts_ci <- confint(pairs(emm), level = 0.95)
write.csv(as.data.frame(contrasts_ci), "analysis/contrasts_ci.csv", row.names = TRUE)

diff_df <- as.data.frame(contrasts_ci)

# ── 5. Model diagnostics ───────────────────────────────
# Make sure the folder "graphs" exists
dir.create("graphs", showWarnings = FALSE)
plot(model)                # Residuals vs fitted
# save figure
png("graphs/residuals_vs_fitted.png", width = 800, height = 600)
plot(model, which = 1)   # 'which = 1' ensures residuals vs fitted
dev.off()

png("graphs/qqplot_residuals.png", width = 800, height = 600)
qqnorm(resid(model))
qqline(resid(model))
dev.off()

# ── 6. Effect size (Cohen's d) ───────────────────────────────
resid_sd <- sigma(model)  # residual SD

residuals_df <- data.frame(fitted = fitted(model), residuals = resid(model))
write.csv(residuals_df, "analysis/residuals_vs_fitted.csv", row.names = FALSE)

# Save sigma (residual SD)
resid_sd <- sigma(model)
writeLines(paste("Residual SD:", resid_sd), "analysis/residual_sd.txt")


ci <- confint(model, parm = "ZonePeripheral", method = "boot", nsim = 1000)

d_est <- fixef(model)["ZonePeripheral"] / resid_sd
d_lower <- ci["ZonePeripheral", "2.5 %"] / resid_sd
d_upper <- ci["ZonePeripheral", "97.5 %"] / resid_sd

cat("Cohen's d =", round(d_est,2),
    "95% CI [", round(d_lower,2), ",", round(d_upper,2), "]\n",
    file = "analysis/cohens_d.txt")

# ── 7. R² (marginal and conditional) with 95% CI ───────────────────────────────
r2_semiparametric <- r2(model, type = "semiparametric", ci = 0.95, bootstrap = 1000)
r2_parametric <- r2(model, type = "parametric", ci = 0.95, bootstrap = 1000)

writeLines(capture.output(r2_semiparametric), "analysis/R2_semiparametric.txt")
writeLines(capture.output(r2_parametric), "analysis/R2_parametric.txt")

# ── 8. Figure (Nature/Elsevier style) ───────────────────────────────

# Panel A: Raw data
p1 <- ggplot(df, aes(x = Zone, y = Value)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +
  geom_quasirandom(width = 0.2, alpha = 0.7, color = "blue") +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  theme_classic() +
  labs(title = "A. Raw Data", y = dep_variable, x = "Zone") +
  theme(plot.title = element_text(face = "bold", size = 12))
# Save the plot
ggsave(filename = "graphs/panelA_raw_data.png", plot = p1,
       width = 6, height = 4, dpi = 300)

# Panel B: Estimated marginal means ± 95% CI
p2 <- ggplot(emm_df, aes(x = Zone, y = emmean)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "darkgreen") +
  theme_classic() +
  labs(title = "B. Estimated Marginal Means", y = "Estimated Value", x = "Zone") +
  theme(plot.title = element_text(face = "bold", size = 12))
# Save as PNG
ggsave("graphs/panelB_emmeans.png", plot = p2,
       width = 6, height = 4, dpi = 300)

# Panel C: Contrast (Peripheral − Central) ± 95% CI
p3 <- ggplot(diff_df, aes(x = contrast, y = estimate)) +
  geom_point(size = 3, color = "purple") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "purple") +
  theme_classic() +
  labs(title = "C. Zone Difference", y = "Difference (Peripheral − Central)", x = "") +
  theme(plot.title = element_text(face = "bold", size = 12))
# Save as PNG
ggsave("graphs/panelC_zone_difference.png", plot = p3,
       width = 6, height = 4, dpi = 300)


# Combine panels (A | B on top, C below)
multi_panel <- (p1 | p2) / p3 + plot_layout(heights = c(2,1))
print(multi_panel)

# Save figure as TIFF, 600 dpi (publication quality)
ggsave("graphs/figure_multi_panel.tiff", multi_panel, width = 6, height = 6, dpi = 600)


# General interpretation file for mixed-model results
interpretation_text <- c(
  "==================== Mixed Model Analysis Interpretation ====================",
  "",
  "1. Fixed Effects:",
  "   - Coefficients represent the estimated effect of each factor while accounting for random effects.",
  "   - Positive coefficient: factor increases the outcome; negative coefficient: factor decreases the outcome.",
  "",
  "2. Random Effects:",
  "   - Variance associated with grouping factors (e.g., subjects, clusters, nested measurements).",
  "   - Larger variance indicates more variability between groups.",
  "",
  "3. ANOVA (Type III or similar):",
  "   - Tests whether a factor significantly affects the outcome while controlling for other factors.",
  "   - p-value < 0.05 generally indicates a significant effect.",
  "",
  "4. Estimated Marginal Means (EMMs):",
  "   - Model-adjusted means for each level of a factor, accounting for random effects.",
  "   - Useful for comparing group means after adjusting for the model structure.",
  "",
  "5. Contrasts:",
  "   - Pairwise comparisons between factor levels (differences in EMMs).",
  "   - If 95% CI does NOT include 0, the difference is considered statistically significant.",
  "",
  "6. Residuals:",
  "   - Residuals = observed value - fitted value.",
  "   - Residuals vs Fitted plot checks homoscedasticity (equal variance).",
  "   - QQ plot of residuals checks normality; points close to the line indicate approximate normal distribution.",
  "",
  "7. Effect Size (Cohen's d):",
  "   - Standardized measure of difference between two groups (estimate / residual SD).",
  "   - Interpretation guidelines:",
  "       d ~ 0.2 → small effect",
  "       d ~ 0.5 → moderate effect",
  "       d ~ 0.8 → large effect",
  "   - Sign indicates direction of effect (positive or negative).",
  "",
  "8. R² for Mixed Models:",
  "   - Marginal R²: proportion of variance explained by fixed effects only.",
  "   - Conditional R²: proportion of variance explained by fixed + random effects.",
  "   - Useful for understanding how much of the outcome variability is captured by the model.",
  "",
  "9. Graphs:",
  "   - Panel A: Raw data with points and summary statistics.",
  "   - Panel B: Estimated marginal means ± 95% CI.",
  "   - Panel C: Contrasts between groups ± 95% CI.",
  "",
  "============================================================================="
)

# Save to text file
dir.create("analysis", showWarnings = FALSE)
writeLines(interpretation_text, "analysis/results_interpretation_general.txt")


# ── EXTRA. Outliers and robustness check ───────────────────────────────
shapiro.test(resid(model))

std_res <- resid(model) / sigma(model)
range(std_res)

which(abs(resid(model) / sigma(model)) > 3)

df[which(abs(resid(model) / sigma(model)) > 3), ]

model_no_outliers <- lmer(Value ~ Zone + (1 | Cow/Vertebra),
                          data = df[-c(27, 31), ])

summary(model)
summary(model_no_outliers)



