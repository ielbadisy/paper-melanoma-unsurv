# ============================================================
# Code for Survival Trajectory Phenotypes in Advanced Melanoma Under Systemic Therapy
# ============================================================

## setup 
pacman::p_load(
  survival, survdnn, unsurv, tvrmst, ggplot2,
  dplyr, tidyr, table1, patchwork
)

set.seed(01042026)

# plotting style

theme_pub <- function() {
  theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.box = "horizontal",
      legend.justification = "center",
      strip.text = element_text(face = "bold")
      )
      }

theme_shared_legend <- function() {
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.title = element_text(face = "bold")
  )
}

guide_pub <- guides(color = guide_legend(nrow = 1, byrow = TRUE))

treatment_type_cols <- c(
  "Dacarbazine" = "#d95f02",
  "Pembrolizumab" = "#1b9e77"
  )

therapy_line_cols <- c(
  "First-line"  = "#1b9e77",
  "Second-line" = "#d95f02"
  )

risk_cols <- c(
  "Low risk"    = "#1b9e77",
  "Medium risk" = "#7570b3",
  "High risk"   = "#d95f02"
  )

# 0) DATA
melanoma2 <- read.csv("melanoma2.csv")

melanoma2 <- melanoma2 |>
  mutate(
    PFS = as.numeric(PFS),
    Event_of_PFS = as.integer(Event_of_PFS)) |>
  mutate(across(c(
    Sex, Histology, Comorbidity, PS,
    Lung_Metastasis, GG_metastasis, Bone_metastasis, Other,
    Therapeutic_line, therapy_type, Response_type), 
    as.factor))

# Therapy labels
melanoma2 <- melanoma2 |>
  mutate(
    therapy_type = factor(
      recode(
        therapy_type,
        dacarbazine = "Dacarbazine",
        pembrolizumab = "Pembrolizumab"
      ),
      levels = c("Dacarbazine", "Pembrolizumab")
    ),
    Therapeutic_line = factor(
      recode(
        Therapeutic_line,
        first_line  = "First-line",
        second_line = "Second-line"
      ),
      levels = c("First-line", "Second-line")
    )
  )

# 1) deep survival model (survdnn with AFT loss)

formula <- Surv(PFS, Event_of_PFS) ~
  Sex + Age +
  Histology +
  Comorbidity +
  PS +
  Lung_Metastasis +
  GG_metastasis +
  Bone_metastasis +
  Other +
  therapy_type +
  Therapeutic_line +
  Total_doses

mod <- survdnn(
  formula,
  data = melanoma2,
  hidden = c(32, 64, 16),
  epochs = 500,
  loss = "aft",
  verbose = TRUE,
  .seed = 01042024
  )

summary(mod)

# 2) Model performance metrics (survdnn)
tmax_pfs_obs <- max(melanoma2$PFS, na.rm = TRUE)
eval_times <- sort(unique(c(12, 24, 36, 48, tmax_pfs_obs)))
eval_times <- eval_times[eval_times <= tmax_pfs_obs]

perf_survdnn <- evaluate_survdnn(
  model   = mod,
  metrics = c("cindex", "brier", "ibs"),
  times   = eval_times
  )

perf_survdnn


# 4) Predict individualized survival curves

tmax <- floor(tmax_pfs_obs)
time_grid <- seq(0, tmax, by = 1)

S <- predict(
  mod,
  melanoma2,
  type = "survival",
  time = time_grid
  )

S <- as.matrix(S)

# Primary clustering horizon
cluster_horizon <- 24
cluster_idx <- which(time_grid <= cluster_horizon)
time_grid_cluster <- time_grid[cluster_idx]
S_cluster <- S[, cluster_idx, drop = FALSE]


# 5) Unsupervised clustering of survival trajectories

fit_unsurv <- unsurv(
  S_cluster,
  time_grid_cluster,
  K = 3
  )

summary(fit_unsurv)

melanoma2$cluster_raw <- factor(fit_unsurv$clusters)

# cluster stability
stab <- unsurv_stability(
  S_cluster,
  time_grid_cluster,
  fit_unsurv,
  B = 50,
  frac = 0.7,
  mode = "subsample"
  )

print(stab$mean)

# Sensitivity: cluster stability over truncated horizons
compute_stability_at_horizon <- function(S, time_grid, horizon, K = 3, B = 50, frac = 0.7) {
  keep_idx <- which(time_grid <= horizon)
  S_h <- S[, keep_idx, drop = FALSE]
  t_h <- time_grid[keep_idx]

  fit_h <- unsurv(
    S_h,
    t_h,
    K = K
  )

  stab_h <- unsurv_stability(
    S_h,
    t_h,
    fit_h,
    B = B,
    frac = frac,
    mode = "subsample"
  )

  data.frame(
    horizon = paste0("0-", horizon, " months"),
    mean_stability = stab_h$mean
  )
}

stability_summary <- bind_rows(
  data.frame(
    horizon = paste0("0-", cluster_horizon, " months (primary)"),
    mean_stability = stab$mean
  ),
  compute_stability_at_horizon(S, time_grid, horizon = 36),
  compute_stability_at_horizon(S, time_grid, horizon = max(time_grid))
)

print(stability_summary)

# ============================================================
# 6) Dynamic RMST computation
# ============================================================

survmat_raw <- as_survmat(
  S,
  time_grid,
  group = melanoma2$cluster_raw
)

rmst_raw <- rmst_dynamic(survmat_raw)

# ------------------------------------------------------------
# Relabel clusters by mean RMST at tau = 24 months
# highest RMST = Low risk ; lowest RMST = High risk
# ------------------------------------------------------------
tau_star <- 24
i_tau <- which.min(abs(rmst_raw$time - tau_star))

cluster_rmst24 <- data.frame(
  cluster_raw = levels(melanoma2$cluster_raw),
  rmst24 = sapply(levels(melanoma2$cluster_raw), function(cl) {
    idx <- which(melanoma2$cluster_raw == cl)
    mean(rmst_raw$individual[idx, i_tau], na.rm = TRUE)
  })
) |>
  arrange(desc(rmst24)) |>
  mutate(
    risk_group = c("Low risk", "Medium risk", "High risk")
  )

print(cluster_rmst24)

risk_map <- setNames(cluster_rmst24$risk_group, cluster_rmst24$cluster_raw)

melanoma2 <- melanoma2 |>
  mutate(
    risk_group = factor(
      risk_map[as.character(cluster_raw)],
      levels = c("Low risk", "Medium risk", "High risk")
    )
  )

# Recompute survmat / RMST with relabeled groups
survmat <- as_survmat(
  S,
  time_grid,
  group = melanoma2$risk_group
)

rmst_res <- rmst_dynamic(survmat)

# ============================================================
# 7) Table 2 by risk phenotype
# ============================================================

melanoma_table2 <- melanoma2 |>
  transmute(
    `Risk group` = risk_group,
    Sex = factor(Sex, levels = c("F", "M"), labels = c("Female", "Male")),
    Age = Age,
    Histology = Histology,
    Comorbidity = factor(Comorbidity, levels = c("no", "yes"), labels = c("No", "Yes")),
    `Performance status` = factor(as.character(PS), levels = c("0", "1", "2", "3")),
    `Lung metastasis` = factor(Lung_Metastasis, levels = c("no", "yes"), labels = c("No", "Yes")),
    `Nodal metastasis` = factor(GG_metastasis, levels = c("no", "yes"), labels = c("No", "Yes")),
    `Bone metastasis` = factor(Bone_metastasis, levels = c("no", "yes"), labels = c("No", "Yes")),
    `Other metastasis` = factor(Other, levels = c("no", "yes"), labels = c("No", "Yes")),
    Treatment = therapy_type,
    `Therapy line` = Therapeutic_line,
    `Total doses` = Total_doses,
    Response = factor(
      Response_type,
      levels = c("NE", "PD", "RC", "stable"),
      labels = c("Non-evaluable", "Progressive disease", "Complete response", "Stable disease")
    )
  )

tab2_risk <- table1(
  ~ Sex + Age + Histology + Comorbidity + `Performance status` +
    `Lung metastasis` + `Nodal metastasis` + `Bone metastasis` +
    `Other metastasis` + Treatment + `Therapy line` + `Total doses` + Response |
    `Risk group`,
  data = melanoma_table2,
  overall = c(left = "Overall"),
  caption = "Clinical characteristics by discovered survival risk phenotype"
)

print(tab2_risk)

tab2_export <- as.data.frame(tab2_risk, stringsAsFactors = FALSE)
names(tab2_export)[1] <- "Characteristic"

write.csv(
  tab2_export,
  "table2_by_risk_group.csv",
  row.names = FALSE,
  na = ""
)

# ============================================================
# 8) Build long-format data for ggplot
# ============================================================

# Predicted survival curves (all individuals)
surv_long <- data.frame(
  id         = rep(seq_len(nrow(S)), times = length(time_grid)),
  time       = rep(time_grid, each = nrow(S)),
  surv       = as.vector(S),
  risk_group = rep(melanoma2$risk_group, times = length(time_grid)),
  treatment_type = rep(melanoma2$therapy_type, times = length(time_grid)),
  therapy_line = rep(melanoma2$Therapeutic_line, times = length(time_grid))
)

# Individual RMST curves
RMST <- rmst_res$individual

rmst_long <- data.frame(
  id         = rep(seq_len(nrow(RMST)), times = ncol(RMST)),
  tau        = rep(time_grid, each = nrow(RMST)),
  rmst       = as.vector(RMST),
  risk_group = rep(melanoma2$risk_group, times = ncol(RMST)),
  treatment_type = rep(melanoma2$therapy_type, times = ncol(RMST)),
  therapy_line = rep(melanoma2$Therapeutic_line, times = ncol(RMST))
)

# ============================================================
# 9) Publication-ready plots (all ggplot, same style)
# ============================================================

# 9A) Mean predicted survival by risk phenotype
# Suggested main-draft caption: Mean predicted PFS trajectories by discovered risk group.
surv_mean_risk <- surv_long |>
  group_by(risk_group, time) |>
  summarise(mean_surv = mean(surv, na.rm = TRUE), .groups = "drop")

p_surv_risk <- ggplot(
  surv_mean_risk,
  aes(x = time, y = mean_surv, color = risk_group)
) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = risk_cols) +
  labs(
    x = "Time (months)",
    y = "Predicted survival probability",
    color = "Risk group"
  ) +
  theme_pub()

print(p_surv_risk)

# 9B) Mean predicted survival by treatment type
# Suggested main-draft caption: Mean predicted PFS trajectories by treatment type.
surv_mean_tx <- surv_long |>
  group_by(treatment_type, time) |>
  summarise(mean_surv = mean(surv, na.rm = TRUE), .groups = "drop")

p_surv_tx <- ggplot(
  surv_mean_tx,
  aes(x = time, y = mean_surv, color = treatment_type)
) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = treatment_type_cols) +
  labs(
    x = "Time (months)",
    y = "Predicted survival probability",
    color = "Treatment type"
  ) +
  theme_pub()

print(p_surv_tx)

# 9C) Secondary result: mean predicted survival by therapy line
# Suggested main-draft caption: Mean predicted PFS trajectories by therapy line.
surv_mean_line <- surv_long |>
  group_by(therapy_line, time) |>
  summarise(mean_surv = mean(surv, na.rm = TRUE), .groups = "drop")

p_surv_line <- ggplot(
  surv_mean_line,
  aes(x = time, y = mean_surv, color = therapy_line)
) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = therapy_line_cols) +
  labs(
    x = "Time (months)",
    y = "Predicted survival probability",
    color = "Therapy line"
  ) +
  theme_pub()

print(p_surv_line)

# 9D) Representative survival curves of unsurv clusters (medoids if available)
# Suggested main-draft caption: Representative predicted survival curves for the discovered risk groups.
get_medoid_indices <- function(fit_obj, n_clusters, cluster_membership) {
  if (!is.null(fit_obj$medoids)) {
    med <- fit_obj$medoids
    if (is.numeric(med) && length(med) == n_clusters) return(as.integer(med))
  }
  # fallback: first subject from each cluster
  sapply(seq_len(n_clusters), function(k) which(cluster_membership == k)[1])
}

medoid_idx <- get_medoid_indices(
  fit_obj = fit_unsurv,
  n_clusters = length(levels(melanoma2$cluster_raw)),
  cluster_membership = as.integer(melanoma2$cluster_raw)
)

medoid_df <- data.frame(
  time = rep(time_grid, times = length(medoid_idx)),
  surv = as.vector(t(S[medoid_idx, , drop = FALSE])),
  cluster_raw = factor(rep(levels(melanoma2$cluster_raw), each = length(time_grid)),
                       levels = levels(melanoma2$cluster_raw))
) |>
  mutate(
    risk_group = factor(
      risk_map[as.character(cluster_raw)],
      levels = c("Low risk", "Medium risk", "High risk")
    )
  )

p_medoids <- ggplot(
  medoid_df,
  aes(x = time, y = surv, color = risk_group)
) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = risk_cols) +
  labs(
    x = "Time (months)",
    y = "Survival probability",
    color = "Risk group"
  ) +
  theme_pub()

print(p_medoids)

# 9E) Mean RMST by risk group
# Suggested main-draft caption: Mean dynamic RMST trajectories by discovered risk group.
rmst_mean_risk <- rmst_long |>
  group_by(risk_group, tau) |>
  summarise(mean_rmst = mean(rmst, na.rm = TRUE), .groups = "drop")

p_rmst_risk <- ggplot(
  rmst_mean_risk,
  aes(x = tau, y = mean_rmst, color = risk_group)
) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = risk_cols) +
  labs(
    x = "Time horizon tau (months)",
    y = "Mean RMST(tau)",
    color = "Risk group"
  ) +
  theme_pub()

print(p_rmst_risk)

# 9F) RMST by risk group and treatment type
# Suggested main-draft caption: Mean dynamic RMST trajectories by treatment type within each discovered risk group.
rmst_mean_rt <- rmst_long |>
  group_by(risk_group, treatment_type, tau) |>
  summarise(mean_rmst = mean(rmst, na.rm = TRUE), .groups = "drop")

p_rmst_tx <- ggplot(
  rmst_mean_rt,
  aes(x = tau, y = mean_rmst, color = treatment_type)
) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ risk_group, nrow = 1) +
  scale_color_manual(values = treatment_type_cols, drop = TRUE) +
  labs(
    x = "Time horizon tau (months)",
    y = "Mean RMST(tau)",
    color = "Treatment type"
  ) +
  theme_pub()

print(p_rmst_tx)

# 9G) Secondary result: RMST by risk group and therapy line
# Suggested main-draft caption: Mean dynamic RMST trajectories by therapy line within each discovered risk group.
rmst_mean_line <- rmst_long |>
  group_by(risk_group, therapy_line, tau) |>
  summarise(mean_rmst = mean(rmst, na.rm = TRUE), .groups = "drop")

p_rmst_line <- ggplot(
  rmst_mean_line,
  aes(x = tau, y = mean_rmst, color = therapy_line)
) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ risk_group, nrow = 1) +
  scale_color_manual(values = therapy_line_cols, drop = TRUE) +
  labs(
    x = "Time horizon tau (months)",
    y = "Mean RMST(tau)",
    color = "Therapy line"
  ) +
  theme_pub()

print(p_rmst_line)

# 9H) Kaplan-Meier validation (ggplot)
# Suggested main-draft caption: Observed Kaplan-Meier PFS curves stratified by discovered risk group.
fit_km <- survfit(
  Surv(PFS, Event_of_PFS) ~ risk_group,
  data = melanoma2
)

km_df <- data.frame(
  time   = fit_km$time,
  surv   = fit_km$surv,
  strata = rep(names(fit_km$strata), fit_km$strata)
) |>
  mutate(
    risk_group = factor(
      gsub("^risk_group=", "", strata),
      levels = c("Low risk", "Medium risk", "High risk")
    )
  ) |>
  bind_rows(
    tibble(
      time = 0,
      surv = 1,
      strata = paste0("risk_group=", levels(melanoma2$risk_group)),
      risk_group = factor(
        levels(melanoma2$risk_group),
        levels = levels(melanoma2$risk_group)
      )
    )
  ) |>
  arrange(risk_group, time)

p_km <- ggplot(
  km_df,
  aes(x = time, y = surv, color = risk_group)
) +
  geom_step(linewidth = 1.2) +
  scale_color_manual(values = risk_cols) +
  labs(
    x = "Time (months)",
    y = "Observed survival probability",
    color = "Risk group"
  ) +
  theme_pub()

print(p_km)

# 9I) Follow-up / observed time distribution
# Suggested main-draft caption: Distribution of observed PFS times across discovered risk groups.
p_follow <- melanoma2 |>
  ggplot(aes(x = risk_group, y = PFS, color = risk_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.15) +
  geom_jitter(width = 0.12, alpha = 0.6) +
  scale_color_manual(values = risk_cols) +
  guides(color = "none") +
  labs(
    x = "Risk group",
    y = "Observed PFS time (months)",
    color = "Risk group"
  ) +
  theme_pub() +
  theme(legend.position = "none")

print(p_follow)

# Suggested main-draft caption: Individual dynamic RMST trajectories and group-level means by discovered risk group.
p_rmst_individual <- ggplot(
  rmst_long,
  aes(x = tau, y = rmst, group = id, color = risk_group)
) +
  geom_line(alpha = 0.3, linewidth = 0.4, show.legend = FALSE) +
  stat_summary(
    aes(group = risk_group),
    fun = mean,
    geom = "line",
    linewidth = 1.8,
    show.legend = FALSE
  ) +
  facet_wrap(~ risk_group, nrow = 1) +   # <- remove scales = "free_y"
  scale_color_manual(values = risk_cols) +
  labs(
    x = "Time horizon t (months)",
    y = "RMST_i(t)",
    color = "Risk group"
  ) +
  theme_pub()

print(p_rmst_individual)

# ============================================================
# 9J) Composite figures for the main manuscript
# ============================================================

figure_1_treatment_rmst <- (
  (
    p_rmst_tx +
      guide_pub +
      theme(legend.position = "bottom") +
      labs(color = "Treatment type")
  ) |
  (
    p_rmst_line +
      guide_pub +
      theme(legend.position = "bottom") +
      labs(color = "Therapy line")
  )
) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme_shared_legend()

figure_2_followup_km <- (
  (
    p_follow +
      theme(legend.position = "none")
  ) |
  (
    p_km +
      guide_pub +
      theme(legend.position = "bottom")
  )
) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme_shared_legend()

figure_3_rmst_risk <- (
  (
    p_rmst_individual +
      theme(legend.position = "none")
  ) |
  (
    p_rmst_risk +
      guide_pub +
      theme(legend.position = "bottom")
  )
) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme_shared_legend()

print(figure_1_treatment_rmst)
print(figure_2_followup_km)
print(figure_3_rmst_risk)
# 10) Diagnostics

print(table(melanoma2$risk_group, melanoma2$Event_of_PFS))

print(aggregate(PFS ~ risk_group, melanoma2, summary))

print(tapply(
  melanoma2$PFS,
  melanoma2$risk_group,
  quantile,
  probs = c(.1, .25, .5, .75, .9),
  na.rm = TRUE
))

# Reverse Kaplan–Meier median follow-up
fit_fu <- survfit(
  Surv(PFS, 1 - Event_of_PFS) ~ 1,
  data = melanoma2
)

median_followup <- summary(fit_fu)$table["median"]
print(median_followup)

# ============================================================
# 11) Draft table exports
# ============================================================

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

write_html_table <- function(df, file, title = NULL) {
  header_cells <- paste0("<th>", html_escape(names(df)), "</th>", collapse = "")
  body_rows <- apply(df, 1, function(row) {
    paste0(
      "<tr>",
      paste0("<td>", html_escape(as.character(row)), "</td>", collapse = ""),
      "</tr>"
    )
  })

  html_lines <- c(
    "<html>",
    "<head>",
    "<meta charset='utf-8' />",
    "<style>",
    "body { font-family: Arial, sans-serif; margin: 24px; }",
    "table { border-collapse: collapse; width: 100%; }",
    "th, td { border: 1px solid #999; padding: 6px 8px; text-align: left; vertical-align: top; }",
    "th { background: #f2f2f2; }",
    "caption { caption-side: top; text-align: left; font-weight: bold; margin-bottom: 8px; }",
    "</style>",
    "</head>",
    "<body>",
    "<table>",
    if (!is.null(title)) paste0("<caption>", html_escape(title), "</caption>"),
    paste0("<thead><tr>", header_cells, "</tr></thead>"),
    "<tbody>",
    body_rows,
    "</tbody>",
    "</table>",
    "</body>",
    "</html>"
  )

  writeLines(html_lines, con = file)
}

writeLines(
  c(
    "<html>",
    "<head><meta charset='utf-8' /></head>",
    "<body>",
    as.character(tab2_risk),
    "</body>",
    "</html>"
  ),
  con = "table2_by_risk_group.html"
)

cluster_rmst24_export <- cluster_rmst24 |>
  rename(
    Raw_cluster = cluster_raw,
    Mean_RMST_at_24_months = rmst24,
    Assigned_risk_group = risk_group
  )

write.csv(
  cluster_rmst24_export,
  "cluster_rmst24_relabeling.csv",
  row.names = FALSE
)

write_html_table(
  transform(
    cluster_rmst24_export,
    Mean_RMST_at_24_months = sprintf("%.3f", Mean_RMST_at_24_months)
  ),
  file = "cluster_rmst24_relabeling.html",
  title = "RMST-based relabeling of discovered clusters at 24 months"
)

rmst_export_times <- c(12, 24, 36, 48)
rmst_selected_horizons <- rmst_long |>
  filter(tau %in% rmst_export_times) |>
  group_by(risk_group, tau) |>
  summarise(
    mean_rmst = mean(rmst, na.rm = TRUE),
    sd_rmst = sd(rmst, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    rmst_mean_sd = sprintf("%.2f (%.2f)", mean_rmst, sd_rmst),
    tau_label = paste0("RMST at ", tau, " months, mean (SD)")
  ) |>
  select(risk_group, tau_label, rmst_mean_sd) |>
  pivot_wider(names_from = tau_label, values_from = rmst_mean_sd) |>
  rename(Risk_group = risk_group)

write.csv(
  rmst_selected_horizons,
  "table1_rmst_by_risk_group_selected_horizons.csv",
  row.names = FALSE
)

write.csv(
  rmst_selected_horizons,
  "rmst_by_risk_group_selected_horizons.csv",
  row.names = FALSE
)

write_html_table(
  rmst_selected_horizons,
  file = "rmst_by_risk_group_selected_horizons.html",
  title = "Dynamic RMST by discovered risk group at selected horizons"
)

# ============================================================
# 11) Save publication figures
# ============================================================

ggsave("predicted_survival_by_risk.png", p_surv_risk, width = 7, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("predicted_survival_by_therapy_type.png", p_surv_tx, width = 7, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("predicted_survival_by_therapy_line.png", p_surv_line, width = 7, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("unsurv_medoid_curves.png", p_medoids, width = 7, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("rmst_by_risk_group.png", p_rmst_risk, width = 7, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("rmst_by_risk_and_therapy_type.png", p_rmst_tx, width = 7, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("rmst_by_risk_and_therapy_line.png", p_rmst_line, width = 7, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("rmst_individual_trajectories_by_risk_group.png", p_rmst_individual, width = 7, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("km_by_risk_group.png", p_km, width = 7, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("followup_by_risk_group.png", p_follow, width = 7, height = 7, units = "in", dpi = 300, bg = "white")

ggsave("figure_1_treatment_rmst.png", figure_1_treatment_rmst, width = 14, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("figure_2_followup_km.png", figure_2_followup_km, width = 14, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("figure_3_rmst_risk.png", figure_3_rmst_risk, width = 14, height = 7, units = "in", dpi = 300, bg = "white")

# Draft-ready multi-page PDF with the three paired manuscript figures
pdf("draft_patched_figures_1_3.pdf", width = 14, height = 7, onefile = TRUE)
print(figure_1_treatment_rmst)
print(figure_2_followup_km)
print(figure_3_rmst_risk)
dev.off()
