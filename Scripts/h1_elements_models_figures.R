suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(performance)
  library(patchwork)
  library(vegan)
})

options(stringsAsFactors = FALSE)
options(emmeans.check.rank = FALSE)
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
arg_get <- function(key, default = NULL) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}

in_file  <- arg_get("in",  "data/processed/All_elements_wide.csv")
out_base <- arg_get("out", "outputs")
control_level <- arg_get("control", "T0:P0")
run_pca_permanova <- as.logical(tolower(arg_get("pca", "true")))
permutations_perm <- as.integer(arg_get("perm", "999"))

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(out_base, paste0("h1_elements_", stamp))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures", "main"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures", "supp"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures", "diagnostics"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures", "pca"), recursive = TRUE, showWarnings = FALSE)

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_max  <- function(x) if (all(is.na(x))) NA_real_ else suppressWarnings(max(x, na.rm = TRUE))

p_fmt <- function(p) ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
sig_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ ".",
    TRUE ~ ""
  )
}

get_ci_cols <- function(nms) {
  lwr <- intersect(c("lower.CL","asymp.LCL","lower","LCL"), nms)[1]
  upr <- intersect(c("upper.CL","asymp.UCL","upper","UCL"), nms)[1]
  if (is.na(lwr)) lwr <- grep("LCL|lower", nms, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(upr)) upr <- grep("UCL|upper", nms, ignore.case = TRUE, value = TRUE)[1]
  list(
    lwr = ifelse(length(lwr) == 0, NA_character_, lwr),
    upr = ifelse(length(upr) == 0, NA_character_, upr)
  )
}

safe_r2 <- function(mod) {
  tryCatch({
    rr <- performance::r2(mod)
    tibble(R2_marginal = rr$R2_marginal, R2_conditional = rr$R2_conditional)
  }, error = function(e) tibble(R2_marginal = NA_real_, R2_conditional = NA_real_))
}

save_plot <- function(p, filename_base, w = 12, h = 8, dpi = 600) {
  png_file <- paste0(filename_base, ".png")
  pdf_file <- paste0(filename_base, ".pdf")

  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(png_file, width = w, height = h, units = "in", res = dpi, background = "white")
    print(p)
    dev.off()
  } else {
    ggsave(png_file, p, width = w, height = h, dpi = dpi, bg = "white")
  }

  if (capabilities("cairo")) {
    ggsave(pdf_file, p, width = w, height = h, device = cairo_pdf)
  } else {
    ggsave(pdf_file, p, width = w, height = h, device = "pdf")
  }

  invisible(list(png = png_file, pdf = pdf_file))
}

save_diag <- function(mod, out_png) {
  png(out_png, width = 2400, height = 1100, res = 260)
  par(mfrow = c(1,2), mar = c(4.5,4.5,1,1), cex.axis = 1.2, cex.lab = 1.4, font.lab = 2)
  plot(fitted(mod), resid(mod), xlab = "Fitted values", ylab = "Residuals", pch = 16, col = "black")
  abline(h = 0, lty = 2, lwd = 2, col = "black")
  qqnorm(resid(mod), pch = 16, col = "black", main = "")
  qqline(resid(mod), lwd = 2, col = "black")
  par(mfrow = c(1,1))
  dev.off()
}

fit_with_fallback <- function(fml, d, singular_tol = 1e-5) {
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

  mod_lmm <- tryCatch(
    suppressMessages(lmerTest::lmer(fml, data = d, control = ctrl, REML = TRUE)),
    error = function(e) NULL
  )

  if (!is.null(mod_lmm)) {
    sing <- lme4::isSingular(mod_lmm, tol = singular_tol)
    vc <- as.data.frame(VarCorr(mod_lmm))
    re_var <- suppressWarnings(vc$vcov[vc$grp == "unit"][1])
    msg <- mod_lmm@optinfo$conv$lme4$messages
    ok <- !isTRUE(sing) && is.finite(re_var) && re_var > 1e-10 && is.null(msg)
    if (ok) return(list(model = mod_lmm, model_class = "LMM"))
  }

  mod_lm <- lm(lme4::nobars(fml), data = d)
  list(model = mod_lm, model_class = "LM")
}

dat0 <- read_csv(in_file, show_col_types = FALSE)

dat0 <- dat0 %>%
  mutate(
    age_f   = factor(age),
    depth_f = factor(depth, levels = c("0-10cm","10-20cm","20-40cm")),
    thin_f  = factor(t_code, levels = 0:4, labels = paste0("T", 0:4)),
    prune_f = factor(p_code, levels = 0:2, labels = paste0("P", 0:2)),
    TP      = interaction(thin_f, prune_f, sep = ":", drop = TRUE),
    unit    = if ("plot_id" %in% names(.)) factor(plot_id) else interaction(age_f, thin_f, prune_f, drop = TRUE)
  )

stopifnot(control_level %in% levels(dat0$TP))

dat <- dat0 %>%
  group_by(age_f, depth_f, thin_f, prune_f, TP, unit) %>%
  summarise(
    Fe     = safe_mean(wi_Fe),
    Cu     = safe_mean(wi_Cu),
    Mn     = safe_mean(wi_Mn),
    Zn     = safe_mean(wi_Zn),
    Mn_imp = safe_max(wi_Mn_imputed_flag),
    Zn_imp = safe_max(wi_Zn_imputed_flag),
    .groups = "drop"
  )

elements <- c("Fe","Cu","Mn","Zn")

pal_age  <- c("#0072B2", "#D55E00", "#009E73")
pal_elem <- c(Fe = "#0072B2", Cu = "#D55E00", Mn = "#009E73", Zn = "#CC79A7")

age_levels <- levels(dat$age_f)
age_cols <- setNames(rep(pal_age, length.out = length(age_levels)), age_levels)

theme_pub <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(face = "bold", colour = "black"),
      axis.title = element_text(size = base_size + 2, face = "bold"),
      axis.text  = element_text(size = base_size + 1, colour = "black"),
      strip.text = element_text(size = base_size + 2, face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.35),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.45),
      legend.title = element_text(size = base_size + 1, face = "bold"),
      legend.text  = element_text(size = base_size),
      legend.position = "bottom",
      plot.title = element_blank(),
      plot.margin = margin(10, 12, 10, 10)
    )
}

fit_element_bundle <- function(el, data = dat, control_level = "T0:P0") {
  y <- data[[el]]
  if (any(is.na(y))) stop("Missing values after aggregation for: ", el)

  eps <- if (any(y <= 0, na.rm = TRUE)) min(y[y > 0], na.rm = TRUE) / 2 else 0

  d <- data %>%
    mutate(resp_log = log(.data[[el]] + eps)) %>%
    filter(is.finite(resp_log))

  fml_main <- resp_log ~ age_f * depth_f + thin_f * prune_f + (1 | unit)
  fit_main <- fit_with_fallback(fml_main, d)
  mod_main <- fit_main$model

  fml_tp <- resp_log ~ age_f * depth_f + TP + (1 | unit)
  fit_tp <- fit_with_fallback(fml_tp, d)
  mod_tp <- fit_tp$model

  r2_tbl <- if (inherits(mod_main, "lmerMod")) {
    safe_r2(mod_main)
  } else {
    tibble(R2_marginal = summary(mod_main)$r.squared, R2_conditional = NA_real_)
  } %>%
    mutate(element = el, model_used = fit_main$model_class)

  emm_ad <- emmeans(mod_main, ~ depth_f | age_f)
  s_ad <- as_tibble(summary(emm_ad, infer = c(TRUE, TRUE), adjust = "none"))
  ci_ad <- get_ci_cols(names(s_ad))
  lwr <- if (!is.na(ci_ad$lwr)) s_ad[[ci_ad$lwr]] else rep(NA_real_, nrow(s_ad))
  upr <- if (!is.na(ci_ad$upr)) s_ad[[ci_ad$upr]] else rep(NA_real_, nrow(s_ad))

  emm_ad_tbl <- s_ad %>%
    transmute(
      element = el,
      age_f, depth_f,
      mean_bt = exp(emmean) - eps,
      lwr_bt  = exp(lwr) - eps,
      upr_bt  = exp(upr) - eps,
      eps_used = eps
    )

  emm0 <- emmeans(mod_tp, ~ TP)
  ref_idx0 <- which(levels(emm0@grid$TP) == control_level)
  ctr0 <- contrast(emm0, method = "trt.vs.ctrl", ref = ref_idx0)
  s0 <- as_tibble(summary(ctr0, infer = c(TRUE, TRUE), adjust = "none"))
  ci0 <- get_ci_cols(names(s0))
  l0 <- if (!is.na(ci0$lwr)) s0[[ci0$lwr]] else rep(NA_real_, nrow(s0))
  u0 <- if (!is.na(ci0$upr)) s0[[ci0$upr]] else rep(NA_real_, nrow(s0))

  tp_overall <- s0 %>%
    mutate(
      element = el,
      TP_comp = str_trim(str_remove(contrast, paste0("\\s*-\\s*", control_level, "\\s*$"))),
      pct_change = 100 * (exp(estimate) - 1),
      pct_lwr    = 100 * (exp(l0) - 1),
      pct_upr    = 100 * (exp(u0) - 1),
      p_adj_fdr  = p.adjust(p.value, method = "BH"),
      sig_adj    = sig_stars(p.adjust(p.value, method = "BH"))
    ) %>%
    filter(TP_comp != control_level) %>%
    select(element, TP_comp, pct_change, pct_lwr, pct_upr, p.value, p_adj_fdr, sig_adj)

  emm1 <- emmeans(mod_tp, ~ TP | age_f * depth_f)
  ref_idx1 <- which(levels(emm1@grid$TP) == control_level)
  ctr1 <- contrast(emm1, method = "trt.vs.ctrl", ref = ref_idx1)
  s1 <- as_tibble(summary(ctr1, infer = c(TRUE, TRUE), adjust = "none"))
  ci1 <- get_ci_cols(names(s1))
  l1 <- if (!is.na(ci1$lwr)) s1[[ci1$lwr]] else rep(NA_real_, nrow(s1))
  u1 <- if (!is.na(ci1$upr)) s1[[ci1$upr]] else rep(NA_real_, nrow(s1))

  tp_strat <- s1 %>%
    mutate(
      element = el,
      TP_comp = str_trim(str_remove(contrast, paste0("\\s*-\\s*", control_level, "\\s*$"))),
      pct_change = 100 * (exp(estimate) - 1),
      pct_lwr    = 100 * (exp(l1) - 1),
      pct_upr    = 100 * (exp(u1) - 1),
      p_adj_fdr  = p.adjust(p.value, method = "BH"),
      sig_adj    = sig_stars(p.adjust(p.value, method = "BH"))
    ) %>%
    filter(TP_comp != control_level) %>%
    select(element, age_f, depth_f, TP_comp, pct_change, pct_lwr, pct_upr, p.value, p_adj_fdr, sig_adj)

  list(
    element = el,
    mod_main = mod_main,
    r2 = r2_tbl,
    emm_ad = emm_ad_tbl,
    tp_overall = tp_overall,
    tp_strat = tp_strat
  )
}

bundles <- lapply(elements, fit_element_bundle, data = dat, control_level = control_level)
names(bundles) <- elements

r2_all         <- bind_rows(lapply(bundles, `[[`, "r2"))
emm_ad_all     <- bind_rows(lapply(bundles, `[[`, "emm_ad"))
tp_overall_all  <- bind_rows(lapply(bundles, `[[`, "tp_overall"))
tp_strat_all    <- bind_rows(lapply(bundles, `[[`, "tp_strat"))

write_csv(r2_all,        file.path(out_dir, "tables", "r2_all_elements.csv"))
write_csv(emm_ad_all,    file.path(out_dir, "tables", "emm_age_depth_backtransformed.csv"))
write_csv(tp_overall_all, file.path(out_dir, "tables", "tp_vs_control_overall.csv"))
write_csv(tp_strat_all,  file.path(out_dir, "tables", "tp_vs_control_by_age_depth.csv"))

for (el in elements) {
  save_diag(bundles[[el]]$mod_main, file.path(out_dir, "figures", "diagnostics", paste0("diag_", el, ".png")))
}

plot_age_depth <- function(el) {
  dd <- emm_ad_all %>%
    filter(element == el) %>%
    mutate(age_f = factor(age_f, levels = age_levels))

  ggplot(dd, aes(x = depth_f, y = mean_bt, group = age_f, colour = age_f)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.8, shape = 21, fill = "white", stroke = 1.0) +
    geom_errorbar(aes(ymin = lwr_bt, ymax = upr_bt), width = 0.10, linewidth = 0.8, colour = "grey35") +
    scale_colour_manual(values = age_cols, name = "Stand age (years)") +
    labs(x = "Soil depth (cm)", y = paste0(el, " (mg kg\u207B\u00B9)")) +
    theme_pub(14)
}

p1 <- (plot_age_depth("Fe") + plot_age_depth("Cu")) /
  (plot_age_depth("Mn") + plot_age_depth("Zn")) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_plot(p1, file.path(out_dir, "figures", "main", "fig1_age_depth_profiles_4elements"), w = 12.8, h = 8.8)

plot_tp_overall <- function(el) {
  dd <- tp_overall_all %>%
    filter(element == el) %>%
    arrange(pct_change) %>%
    mutate(TP_comp = factor(TP_comp, levels = unique(TP_comp))) %>%
    mutate(
      x_star = ifelse(pct_change >= 0, pct_upr + 2.5, pct_lwr - 2.5),
      hjust_star = ifelse(pct_change >= 0, 0, 1)
    )

  ggplot(dd, aes(x = pct_change, y = TP_comp)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.9, colour = "black") +
    geom_errorbarh(aes(xmin = pct_lwr, xmax = pct_upr), height = 0.20, linewidth = 0.9, colour = "grey25") +
    geom_point(size = 2.8, colour = pal_elem[[el]]) +
    geom_text(aes(x = x_star, label = sig_adj, hjust = hjust_star), size = 5.2, colour = pal_elem[[el]], fontface = "bold") +
    labs(x = paste0("% change vs ", control_level), y = "Treatment (TP)") +
    theme_pub(14) +
    theme(legend.position = "none")
}

p2 <- (plot_tp_overall("Fe") + plot_tp_overall("Cu")) /
  (plot_tp_overall("Mn") + plot_tp_overall("Zn")) +
  plot_annotation(tag_levels = "A")

save_plot(p2, file.path(out_dir, "figures", "main", "fig2_tp_vs_control_overall_4elements"), w = 13.2, h = 10.2)

plot_tp_strat <- function(el) {
  dd <- tp_strat_all %>%
    filter(element == el) %>%
    mutate(
      age_f = factor(age_f, levels = age_levels),
      TP_comp = fct_infreq(TP_comp)
    ) %>%
    mutate(
      x_star = ifelse(pct_change >= 0, pct_upr + 2.0, pct_lwr - 2.0),
      hjust_star = ifelse(pct_change >= 0, 0, 1)
    )

  ggplot(dd, aes(x = pct_change, y = TP_comp)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.7, colour = "black") +
    geom_errorbarh(aes(xmin = pct_lwr, xmax = pct_upr), height = 0.16, linewidth = 0.7, colour = "grey30") +
    geom_point(size = 2.0, colour = pal_elem[[el]]) +
    geom_text(aes(x = x_star, label = sig_adj, hjust = hjust_star), size = 4.2, colour = pal_elem[[el]], fontface = "bold") +
    facet_grid(age_f ~ depth_f) +
    labs(x = paste0("% change vs ", control_level), y = "Treatment (TP)") +
    theme_pub(12) +
    theme(legend.position = "none", strip.text = element_text(size = 13))
}

for (el in elements) {
  g_el <- plot_tp_strat(el)
  save_plot(g_el, file.path(out_dir, "figures", "supp", paste0("supp_tp_vs_control_by_age_depth_", el)), w = 14.5, h = 9.5)
}

r2_long <- r2_all %>%
  pivot_longer(cols = c(R2_marginal, R2_conditional), names_to = "R2_type", values_to = "R2") %>%
  mutate(
    R2_type = recode(R2_type, R2_marginal = "Marginal R\u00B2", R2_conditional = "Conditional R\u00B2"),
    element = factor(element, levels = elements)
  )

g_r2 <- ggplot(r2_long, aes(x = element, y = R2, fill = R2_type)) +
  geom_col(position = position_dodge(width = 0.72), colour = "black", linewidth = 0.7) +
  scale_fill_manual(values = c("Marginal R\u00B2" = "grey70", "Conditional R\u00B2" = "grey40")) +
  labs(x = "Element", y = "R\u00B2") +
  theme_pub(14) +
  theme(legend.position = "bottom")

save_plot(g_r2, file.path(out_dir, "figures", "supp", "supp_r2_summary"), w = 8.8, h = 5.6)

if (isTRUE(run_pca_permanova)) {

  pca_eps <- function(x) {
    xp <- x[x > 0 & is.finite(x)]
    if (length(xp) == 0) return(1e-6)
    min(xp) / 2
  }

  run_pca_perm <- function(depth_level, data = dat, permutations = 999) {
    d <- data %>%
      filter(depth_f == depth_level) %>%
      select(Fe, Cu, Mn, Zn, age_f, thin_f, prune_f)

    e <- c(Fe = pca_eps(d$Fe), Cu = pca_eps(d$Cu), Mn = pca_eps(d$Mn), Zn = pca_eps(d$Zn))
    d <- d %>%
      mutate(
        Fe = log(Fe + e["Fe"]),
        Cu = log(Cu + e["Cu"]),
        Mn = log(Mn + e["Mn"]),
        Zn = log(Zn + e["Zn"])
      )

    X <- d %>% select(Fe, Cu, Mn, Zn) %>% scale(center = TRUE, scale = TRUE)
    pca <- prcomp(X, center = FALSE, scale. = FALSE)

    scores <- as_tibble(pca$x) %>%
      bind_cols(d %>% select(age_f, thin_f, prune_f))

    distX <- dist(X, method = "euclidean")
    perm <- adonis2(distX ~ age_f + thin_f * prune_f, data = d, permutations = permutations, by = "margin")

    list(depth = depth_level, pca = pca, scores = scores, perm = perm)
  }

  depths <- levels(dat$depth_f)
  pca_list <- lapply(depths, run_pca_perm, data = dat, permutations = permutations_perm)

  thin_levels <- levels(dat$thin_f)
  thin_cols <- setNames(c("#0072B2","#E69F00","#009E73","#56B4E9","#D55E00"), thin_levels)

  for (pp in pca_list) {
    dep <- pp$depth

    perm_df <- as.data.frame(pp$perm) %>%
      tibble::rownames_to_column("term") %>%
      as_tibble()
    write_csv(perm_df, file.path(out_dir, "tables", paste0("permanova_", dep, ".csv")))

    sc <- pp$scores %>%
      mutate(
        age_f = factor(age_f, levels = age_levels),
        thin_f = factor(thin_f, levels = thin_levels),
        prune_f = factor(prune_f)
      )

    ve <- (pp$pca$sdev^2) / sum(pp$pca$sdev^2)
    xlab <- paste0("PC1 (", sprintf("%.1f", 100 * ve[1]), "%)")
    ylab <- paste0("PC2 (", sprintf("%.1f", 100 * ve[2]), "%)")

    g_pca <- ggplot(sc, aes(x = PC1, y = PC2)) +
      stat_ellipse(aes(colour = thin_f, group = thin_f), type = "norm", level = 0.68, linewidth = 0.7, alpha = 0.25) +
      geom_point(aes(colour = thin_f, shape = prune_f), size = 2.6, alpha = 0.9) +
      facet_wrap(~ age_f, nrow = 1) +
      scale_colour_manual(values = thin_cols, name = "Thinning") +
      labs(x = xlab, y = ylab) +
      coord_equal() +
      theme_pub(13) +
      theme(legend.position = "bottom")

    save_plot(g_pca, file.path(out_dir, "figures", "pca", paste0("pca_", dep)), w = 14.2, h = 5.2)
  }
}

readme <- c(
  "Outputs: H1 elements",
  "",
  "Figures:",
  "  figures/main/ fig1_age_depth_profiles_4elements.(png/pdf)",
  "  figures/main/ fig2_tp_vs_control_overall_4elements.(png/pdf)",
  "  figures/supp/ supp_tp_vs_control_by_age_depth_<Element>.(png/pdf)",
  "  figures/supp/ supp_r2_summary.(png/pdf)",
  "  figures/diagnostics/ diag_<Element>.png",
  "  figures/pca/ pca_<depth>.(png/pdf) (optional)",
  "",
  "Tables:",
  "  tables/ r2_all_elements.csv",
  "  tables/ emm_age_depth_backtransformed.csv",
  "  tables/ tp_vs_control_overall.csv",
  "  tables/ tp_vs_control_by_age_depth.csv",
  "  tables/ permanova_<depth>.csv (optional)"
)
writeLines(readme, file.path(out_dir, "README.txt"))

message("Done: ", out_dir)
