# ============================================================
# Bony Endplate Tortuosity — MANUAL Power Analysis
# ============================================================
# This script avoids simr::powerSim / doTest entirely.
#
# Strategy:
# 1) Fit pilot model with REML to estimate variance components
# 2) Simulate new datasets manually from the fitted hierarchical model
# 3) Refit full and reduced models with ML
# 4) Use stats::anova(full, reduced) likelihood-ratio test for Zone
# 5) Estimate:
#    - power for the CURRENT design
#    - power curve by number of cows
#    - optional balanced power grid for cows x vertebrae per cow
# ============================================================

# -----------------------------
# 0. Packages
# -----------------------------
packages <- c("readxl", "dplyr", "lme4")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!p %in% installed) install.packages(p)
  library(p, character.only = TRUE)
}

set.seed(123)

# -----------------------------
# 1. User settings
# -----------------------------
setwd("C:/Users/Ahmad Alminnawi/Desktop/pipeline_publication/modeling-part/submitted_version/second submission/media-guidlines-rev/media/prism/tortuosity/power_analysis")

file_path    <- "tortuosity_R.xlsx"
sheet_name   <- "Sheet1"
response_var <- "Value"      # change to "Tortuosity" if needed
cow_var      <- "Cow"
vertebra_var <- "Vertebra"
zone_var     <- "Zone"

output_dir <- "power_Tortuosity_manual"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

nsim_current <- 500   # increase to 1000 for final reporting
nsim_curve   <- 300   # increase later if needed
target_power <- 0.80
alpha_level  <- 0.05

# Optional:
# overwrite the observed zone effect with a biologically meaningful one
# Example: effect_override <- -0.05
effect_override <- NULL

# -----------------------------
# 2. Read and clean data
# -----------------------------
df0 <- read_excel(file_path, sheet = sheet_name)

df <- df0 %>%
  dplyr::select(
    Cow      = all_of(cow_var),
    Vertebra = all_of(vertebra_var),
    Zone     = all_of(zone_var),
    Tortuosity = all_of(response_var)
  ) %>%
  mutate(
    Cow = factor(Cow),
    Vertebra = factor(Vertebra),
    Zone = as.character(Zone)
  )

df$Zone <- dplyr::case_when(
  df$Zone %in% c("Central", "central", "CENTER", "Center", "C") ~ "Central",
  df$Zone %in% c("Peripheral", "peripheral", "Periphery", "P", "P1", "P2") ~ "Peripheral",
  TRUE ~ df$Zone
)

df$Zone <- factor(df$Zone, levels = c("Central", "Peripheral"))

df <- df %>% filter(!is.na(Cow), !is.na(Vertebra), !is.na(Zone), !is.na(Tortuosity))

# IMPORTANT:
# Create explicit nested grouping column so nothing depends on on-the-fly interaction parsing
df$CowVertebra <- interaction(df$Cow, df$Vertebra, drop = TRUE)

sink(file.path(output_dir, "00_design_overview.txt"))
cat("========== Tortuosity MANUAL POWER ANALYSIS ==========\n\n")
cat("Rows used:", nrow(df), "\n")
cat("Cows:", nlevels(df$Cow), "\n")
cat("Vertebrae:", nlevels(df$CowVertebra), "\n\n")
cat("Zone counts:\n")
print(table(df$Zone))
cat("\nObservations per vertebra:\n")
print(table(df$CowVertebra))
sink()

# -----------------------------
# 3. Fit pilot model with REML
# -----------------------------
# Use REML for parameter estimation
pilot_reml <- lme4::lmer(
  Tortuosity ~ Zone + (1 | Cow) + (1 | CowVertebra),
  data = df,
  REML = TRUE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

writeLines(capture.output(summary(pilot_reml)),
           file.path(output_dir, "01_pilot_model_summary_REML.txt"))

# Extract parameters used for simulation
beta <- fixef(pilot_reml)

if (!is.null(effect_override)) {
  beta["ZonePeripheral"] <- effect_override
}

vc <- as.data.frame(VarCorr(pilot_reml))
sd_cow   <- vc$sdcor[vc$grp == "Cow"]
sd_cv    <- vc$sdcor[vc$grp == "CowVertebra"]
sd_resid <- sigma(pilot_reml)

param_text <- c(
  paste("Intercept used for simulation:", beta["(Intercept)"]),
  paste("ZonePeripheral effect used for simulation:", beta["ZonePeripheral"]),
  paste("Cow SD:", sd_cow),
  paste("CowVertebra SD:", sd_cv),
  paste("Residual SD:", sd_resid)
)

writeLines(param_text,
           file.path(output_dir, "02_simulation_parameters.txt"))

# -----------------------------
# 4. Helper functions
# -----------------------------

# Simulate one response vector from the hierarchical LMM
simulate_Tortuosity <- function(dat, beta, sd_cow, sd_cv, sd_resid) {
  dat$Zone <- factor(dat$Zone, levels = c("Central", "Peripheral"))
  dat$Cow <- factor(dat$Cow)
  dat$CowVertebra <- factor(dat$CowVertebra)
  
  X <- model.matrix(~ Zone, data = dat)
  
  cow_levels <- levels(dat$Cow)
  cv_levels  <- levels(dat$CowVertebra)
  
  b_cow <- stats::rnorm(length(cow_levels), mean = 0, sd = sd_cow)
  names(b_cow) <- cow_levels
  
  b_cv <- stats::rnorm(length(cv_levels), mean = 0, sd = sd_cv)
  names(b_cv) <- cv_levels
  
  mu <- as.numeric(X %*% beta) +
    b_cow[as.character(dat$Cow)] +
    b_cv[as.character(dat$CowVertebra)]
  
  y <- mu + stats::rnorm(nrow(dat), mean = 0, sd = sd_resid)
  y
}

# One simulation run:
# fit full and reduced models with ML and test Zone by LRT
one_power_run <- function(dat, beta, sd_cow, sd_cv, sd_resid, alpha = 0.05) {
  dat <- dat
  dat$y_sim <- simulate_Tortuosity(dat, beta, sd_cow, sd_cv, sd_resid)
  
  full_fit <- tryCatch(
    suppressWarnings(
      suppressMessages(
        lme4::lmer(
          y_sim ~ Zone + (1 | Cow) + (1 | CowVertebra),
          data = dat,
          REML = FALSE,
          control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
        )
      )
    ),
    error = function(e) NULL
  )
  
  red_fit <- tryCatch(
    suppressWarnings(
      suppressMessages(
        lme4::lmer(
          y_sim ~ (1 | Cow) + (1 | CowVertebra),
          data = dat,
          REML = FALSE,
          control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
        )
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(full_fit) || is.null(red_fit)) return(NA)
  
  lrt <- tryCatch(
    suppressWarnings(
      suppressMessages(
        stats::anova(red_fit, full_fit)
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(lrt)) return(NA)
  if (!("Pr(>Chisq)" %in% colnames(lrt))) return(NA)
  if (nrow(lrt) < 2) return(NA)
  
  pval <- lrt$`Pr(>Chisq)`[2]
  if (is.na(pval)) return(NA)
  
  pval < alpha
}

# Summarize logical/NA simulation results
summarize_power <- function(x) {
  valid <- !is.na(x)
  trials <- sum(valid)
  successes <- sum(x[valid])
  
  if (trials == 0) {
    return(data.frame(
      successes = 0,
      trials = 0,
      power = NA,
      lower = NA,
      upper = NA,
      errors = sum(is.na(x))
    ))
  }
  
  ci <- binom.test(successes, trials)$conf.int
  
  data.frame(
    successes = successes,
    trials = trials,
    power = successes / trials,
    lower = ci[1],
    upper = ci[2],
    errors = sum(is.na(x))
  )
}

# Extend current design by duplicating existing cow-level row patterns
# This preserves the observed unbalanced structure approximately.
extend_design_by_cows <- function(df, n_cows) {
  base_cows <- unique(as.character(df$Cow))
  out <- vector("list", n_cows)
  
  for (i in seq_len(n_cows)) {
    src_cow <- base_cows[((i - 1) %% length(base_cows)) + 1]
    tmp <- df[df$Cow == src_cow, c("Vertebra", "Zone"), drop = FALSE]
    
    old_verts <- unique(as.character(tmp$Vertebra))
    new_verts <- paste0("V", seq_along(old_verts))
    vmap <- setNames(new_verts, old_verts)
    
    tmp$Cow <- paste0("C", i)
    tmp$Vertebra <- vmap[as.character(tmp$Vertebra)]
    tmp$CowVertebra <- interaction(tmp$Cow, tmp$Vertebra, drop = TRUE)
    
    out[[i]] <- tmp[, c("Cow", "Vertebra", "CowVertebra", "Zone")]
  }
  
  newdat <- do.call(rbind, out)
  newdat$Cow <- factor(newdat$Cow)
  newdat$Vertebra <- factor(newdat$Vertebra)
  newdat$CowVertebra <- factor(newdat$CowVertebra)
  newdat$Zone <- factor(newdat$Zone, levels = c("Central", "Peripheral"))
  newdat
}

# Optional balanced design generator:
# n_cows cows, n_vert_per_cow vertebrae per cow, 3 zone rows per vertebra
# using the observed 1 central + 2 peripheral pattern
make_balanced_design <- function(n_cows, n_vert_per_cow) {
  one_vert <- data.frame(
    Zone = factor(c("Central", "Peripheral", "Peripheral"),
                  levels = c("Central", "Peripheral"))
  )
  
  out <- list()
  idx <- 1
  
  for (i in seq_len(n_cows)) {
    for (j in seq_len(n_vert_per_cow)) {
      tmp <- one_vert
      tmp$Cow <- paste0("C", i)
      tmp$Vertebra <- paste0("V", j)
      tmp$CowVertebra <- interaction(tmp$Cow, tmp$Vertebra, drop = TRUE)
      out[[idx]] <- tmp[, c("Cow", "Vertebra", "CowVertebra", "Zone")]
      idx <- idx + 1
    }
  }
  
  dat <- do.call(rbind, out)
  dat$Cow <- factor(dat$Cow)
  dat$Vertebra <- factor(dat$Vertebra)
  dat$CowVertebra <- factor(dat$CowVertebra)
  dat$Zone <- factor(dat$Zone, levels = c("Central", "Peripheral"))
  dat
}

# Run many simulations on a given design
estimate_power_manual <- function(dat, nsim, beta, sd_cow, sd_cv, sd_resid, alpha = 0.05) {
  out <- rep(NA, nsim)
  
  for (i in seq_len(nsim)) {
    out[i] <- one_power_run(
      dat = dat,
      beta = beta,
      sd_cow = sd_cow,
      sd_cv = sd_cv,
      sd_resid = sd_resid,
      alpha = alpha
    )
  }
  
  summarize_power(out)
}

# -----------------------------
# 5. Current-design power
# -----------------------------
power_current <- estimate_power_manual(
  dat = df[, c("Cow", "Vertebra", "CowVertebra", "Zone")],
  nsim = nsim_current,
  beta = beta,
  sd_cow = sd_cow,
  sd_cv = sd_cv,
  sd_resid = sd_resid,
  alpha = alpha_level
)

write.csv(power_current,
          file.path(output_dir, "03_power_current_design.csv"),
          row.names = FALSE)

print(power_current)

# -----------------------------
# 6. Power curve by number of cows
# -----------------------------
cow_breaks <- seq(5, 60, by = 5)

curve_cows <- lapply(cow_breaks, function(nc) {
  dat_nc <- extend_design_by_cows(df, n_cows = nc)
  
  res <- estimate_power_manual(
    dat = dat_nc,
    nsim = nsim_curve,
    beta = beta,
    sd_cow = sd_cow,
    sd_cv = sd_cv,
    sd_resid = sd_resid,
    alpha = alpha_level
  )
  
  cbind(
    n_cows = nc,
    n_rows = nrow(dat_nc),
    n_vertebrae = nlevels(dat_nc$CowVertebra),
    res
  )
})

curve_cows <- do.call(rbind, curve_cows)

write.csv(curve_cows,
          file.path(output_dir, "04_power_curve_by_cows.csv"),
          row.names = FALSE)

png(file.path(output_dir, "04_power_curve_by_cows.png"), width = 900, height = 650)
plot(curve_cows$n_cows, curve_cows$power,
     type = "b", pch = 16, ylim = c(0, 1),
     xlab = "Number of cows", ylab = "Power",
     main = "Manual power curve for Tortuosity by number of cows")
arrows(curve_cows$n_cows, curve_cows$lower,
       curve_cows$n_cows, curve_cows$upper,
       angle = 90, code = 3, length = 0.05)
abline(h = target_power, lty = 2)
dev.off()

# First number of cows reaching target power
idx80_cows <- which(curve_cows$power >= target_power)
n_cows_80 <- if (length(idx80_cows) == 0) NA else curve_cows$n_cows[min(idx80_cows)]

writeLines(
  c(
    paste("Target power:", target_power),
    paste("Estimated cows needed for 80% power:", n_cows_80)
  ),
  file.path(output_dir, "05_required_cows_for_80pct_power.txt")
)

# -----------------------------
# 7. Optional balanced design grid:
#    cows x vertebrae per cow
# -----------------------------
cow_grid <- seq(5, 60, by = 5)
vert_grid <- c(6, 7, 8, 9, 10)

grid_results <- list()
k <- 1

for (nc in cow_grid) {
  for (nv in vert_grid) {
    dat_grid <- make_balanced_design(n_cows = nc, n_vert_per_cow = nv)
    
    res <- estimate_power_manual(
      dat = dat_grid,
      nsim = 150,   # keep modest; increase later if needed
      beta = beta,
      sd_cow = sd_cow,
      sd_cv = sd_cv,
      sd_resid = sd_resid,
      alpha = alpha_level
    )
    
    grid_results[[k]] <- cbind(
      n_cows = nc,
      vertebrae_per_cow = nv,
      n_rows = nrow(dat_grid),
      n_vertebrae = nlevels(dat_grid$CowVertebra),
      res
    )
    k <- k + 1
  }
}

grid_results <- do.call(rbind, grid_results)

write.csv(grid_results,
          file.path(output_dir, "06_power_grid_cows_by_vertebrae.csv"),
          row.names = FALSE)

grid80 <- grid_results[grid_results$power >= target_power, , drop = FALSE]
if (nrow(grid80) > 0) {
  grid80 <- grid80[order(grid80$n_cows, grid80$vertebrae_per_cow), ]
  best_grid <- grid80[1, , drop = FALSE]
  
  write.csv(best_grid,
            file.path(output_dir, "07_first_balanced_design_reaching_80pct_power.csv"),
            row.names = FALSE)
}

