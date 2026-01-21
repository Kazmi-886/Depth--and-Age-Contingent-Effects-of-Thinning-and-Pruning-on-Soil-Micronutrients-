suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(performance)
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
runA     <- as.logical(tolower(arg_get("runA", "true")))
runB     <- as.logical(tolower(arg_get("runB", "true")))
perm_n   <- as.integer(arg_get("perm", "999"))

elements <- c("Fe","Cu","Mn","Zn")

need_pkgs <- c("dplyr","tidyr","tibble","readr","stringr","ggplot2","lme4","lmerTest","emmeans","performance")
missing <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nInstall them (e.g., install.packages(...)) or restore your renv environment.")
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(out_base, paste0("h2_h3_objective3_moderation_", stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures", "optionA_factor"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures", "optionB_intensity"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures", "diagnostics"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

writeLines(
  c(
    paste0("run_started: ", Sys.time()),
    paste0("wd: ", getwd()),
    paste0("in_file: ", in_file),
    paste0("out_dir: ", out_dir),
    paste0("runA: ", runA),
    paste0("runB: ", runB),
    paste0("perm: ", perm_n)
  ),
  con = file.path(out_dir, "logs", "run_metadata.txt")
)

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

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_max  <- function(x) if (all(is.na(x))) NA_real_ else suppressWarnings(max(x, na.rm = TRUE))

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

get_trend_col <- function(nms) {
  cand <- intersect(c("trend","estimate","emtrend","slope"), nms)
  if (length(cand) > 0) return(cand[1])
  g <- grep("trend", nms, ignore.case = TRUE, value = TRUE)
  if (length(g) > 0) return(g[1])
  stop("No slope column detected in emtrends summary. Columns: ", paste(nms, collapse = ", "))
}

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(face = "bold", colour = "black"),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      axis.text  = element_text(size = base_size, colour = "black"),
      strip.text = element_text(size = base_size + 1, face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.35),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.45),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_blank(),
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_plot <- function(p, filename_base, w = 11.5, h = 4.4, dpi = 600) {
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
  png(out_png, width = 1800, height = 800, res = 250)
  par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
  plot(fitted(mod), resid(mod), xlab = "Fitted", ylab = "Residuals", pch = 16, col = "black")
  abline(h = 0, lty = 2, col = "black")
  qqnorm(resid(mod), pch = 16, col = "black", main = "")
  qqline(resid(mod), col = "black", lwd = 2)
  par(mfrow = c(1, 1))
  dev.off()
}

safe_r2_tbl <- function(mod) {
  if (inherits(mod, "lmerMod")) {
    tryCatch({
      rr <- performance::r2(mod)
      tibble(R2_marginal = rr$R2_marginal, R2_conditional = rr$R2_conditional)
    }, error = function(e) tibble(R2_marginal = NA_real_, R2_conditional = NA_real_))
  } else {
    tibble(R2_marginal = summary(mod)$r.squared, R2_conditional = NA_real_)
  }
}

fit_with_fallback <- function(fml, d, singular_tol = 1e-5) {
  ctrl1 <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  mod_lmm <- tryCatch(
    suppressMessages(lmerTest::lmer(fml, data = d, control = ctrl1, REML = TRUE)),
    error = function(e) NULL
  )

  if (!is.null(mod_lmm)) {
    singular_flag <- lme4::isSingular(mod_lmm, tol = singular_tol)
    vc <- as.data.frame(lme4::VarCorr(mod_lmm))
    re_var <- suppressWarnings(vc$vcov[vc$grp == "unit"][1])
    msg <- mod_lmm@optinfo$conv$lme4$messages
    ok <- !isTRUE(singular_flag) && is.finite(re_var) && re_var > 1e-10 && is.null(msg)
    if (ok) return(list(model = mod_lmm, model_class = "LMM"))
  }

  mod_lm <- lm(lme4::nobars(fml), data = d)
  list(model = mod_lm, model_class = "LM")
}

extract_type3 <- function(mod, element, model_tag) {
  has_car <- requireNamespace("car", quietly = TRUE)

  if (inherits(mod, "lmerMod")) {
    a <- as.data.frame(anova(mod, type = 3, ddf = "Satterthwaite")) |>
      tibble::rownames_to_column("term") |>
      dplyr::filter(!term %in% c("(Intercept)", "Residuals"))

    out <- tibble(
      element = element,
      model_tag = model_tag,
      anova_type = "TypeIII",
      term = a$term,
      NumDF = a$NumDF,
      DenDF = a$DenDF,
      F_value = a$`F value`,
      p_value = a$`Pr(>F)`
    )
  } else {
    if (has_car) {
      a <- as.data.frame(car::Anova(mod, type = 3)) |>
        tibble::rownames_to_column("term") |>
        dplyr::filter(!term %in% c("(Intercept)", "Residuals"))

      out <- tibble(
        element = element,
        model_tag = model_tag,
        anova_type = "TypeIII",
        term = a$term,
        NumDF = a$Df,
        DenDF = df.residual(mod),
        F_value = a$`F value`,
        p_value = a$`Pr(>F)`
      )
    } else {
      a <- as.data.frame(anova(mod)) |>
        tibble::rownames_to_column("term") |>
        dplyr::filter(!term %in% c("(Intercept)", "Residuals"))

      out <- tibble(
        element = element,
        model_tag = model_tag,
        anova_type = "TypeI_fallback",
        term = a$term,
        NumDF = a$Df,
        DenDF = df.residual(mod),
        F_value = a$`F value`,
        p_value = a$`Pr(>F)`
      )
    }
  }

  out |>
    mutate(p_fmt = p_fmt(p_value), sig = sig_stars(p_value))
}

prep_log <- function(data, resp) {
  y <- data[[resp]]
  if (any(is.na(y))) stop("Missing values after aggregation for: ", resp)
  eps <- if (any(y <= 0, na.rm = TRUE)) min(y[y > 0], na.rm = TRUE) / 2 else 0
  d <- data |>
    mutate(resp_log = log(.data[[resp]] + eps)) |>
    filter(is.finite(resp_log))
  list(d = d, eps = eps)
}

dat0 <- readr::read_csv(in_file, show_col_types = FALSE)

id_col <- dplyr::case_when(
  "plot_id" %in% names(dat0) ~ "plot_id",
  "plot"    %in% names(dat0) ~ "plot",
  TRUE ~ NA_character_
)

if (is.na(id_col)) {
  dat0 <- dat0 |>
    mutate(unit = interaction(age, t_code, p_code, drop = TRUE))
} else {
  dat0 <- dat0 |>
    mutate(unit = factor(.data[[id_col]]))
}

dat0 <- dat0 |>
  mutate(
    age_f   = factor(age),
    depth_f = factor(depth, levels = c("0-10cm","10-20cm","20-40cm")),
    thin_f  = factor(t_code, levels = 0:4, labels = paste0("T", 0:4)),
    prune_f = factor(p_code, levels = 0:2, labels = paste0("P", 0:2)),
    TP      = interaction(thin_f, prune_f, sep = ":", drop = TRUE),
    thin_pct  = as.numeric(thinning_pct),
    prune_pct = as.numeric(pruning_pct)
  )

dat <- dat0 |>
  group_by(age_f, depth_f, thin_f, prune_f, TP, unit) |>
  summarise(
    Fe = safe_mean(wi_Fe),
    Cu = safe_mean(wi_Cu),
    Mn = safe_mean(wi_Mn),
    Zn = safe_mean(wi_Zn),
    Mn_imp = safe_max(wi_Mn_imputed_flag),
    Zn_imp = safe_max(wi_Zn_imputed_flag),
    thin_pct  = safe_mean(thin_pct),
    prune_pct = safe_mean(prune_pct),
    .groups = "drop"
  ) |>
  mutate(
    thin10  = (thin_pct  - mean(thin_pct,  na.rm = TRUE)) / 10,
    prune10 = (prune_pct - mean(prune_pct, na.rm = TRUE)) / 10
  )

writeLines(
  c(
    paste0("n_rows: ", nrow(dat)),
    paste0("n_units: ", n_distinct(dat$unit)),
    paste0("age_levels: ", paste(levels(dat$age_f), collapse = ", ")),
    paste0("depth_levels: ", paste(levels(dat$depth_f), collapse = ", ")),
    paste0("tp_levels: ", nlevels(dat$TP))
  ),
  con = file.path(out_dir, "logs", "data_summary.txt")
)

emm_bt <- function(mod, eps, spec, label, element, model_tag) {
  em <- if (inherits(mod, "lmerMod")) emmeans(mod, spec, lmer.df = "satterthwaite") else emmeans(mod, spec)
  s <- as_tibble(summary(em, infer = c(TRUE, TRUE), adjust = "none"))
  ci <- get_ci_cols(names(s))
  lwr <- if (!is.na(ci$lwr)) s[[ci$lwr]] else rep(NA_real_, nrow(s))
  upr <- if (!is.na(ci$upr)) s[[ci$upr]] else rep(NA_real_, nrow(s))

  s |>
    mutate(
      element = element,
      model_tag = model_tag,
      label = label,
      eps_used = eps,
      mean_bt = exp(emmean) - eps,
      lwr_bt  = exp(lwr)   - eps,
      upr_bt  = exp(upr)   - eps
    )
}

run_optionA <- function(resp, data) {
  pl <- prep_log(data, resp)
  d <- pl$d
  eps <- pl$eps

  fml <- resp_log ~
    age_f * depth_f +
    thin_f * prune_f +
    age_f:thin_f + age_f:prune_f +
    depth_f:thin_f + depth_f:prune_f +
    (1 | unit)

  fit <- fit_with_fallback(fml, d)
  mod <- fit$model
  model_tag <- paste0("OptionA_factor_", fit$model_class)

  anov <- extract_type3(mod, element = resp, model_tag = model_tag)
  r2v  <- safe_r2_tbl(mod) |> mutate(element = resp, model_tag = model_tag)

  emms <- bind_rows(
    emm_bt(mod, eps, ~ thin_f  | depth_f, "thin_by_depth",  resp, model_tag),
    emm_bt(mod, eps, ~ prune_f | depth_f, "prune_by_depth", resp, model_tag),
    emm_bt(mod, eps, ~ thin_f  | age_f,   "thin_by_age",    resp, model_tag),
    emm_bt(mod, eps, ~ prune_f | age_f,   "prune_by_age",   resp, model_tag)
  )

  plot_emm <- function(df, xvar, facetvar, fname) {
    g <- ggplot(df, aes(x = .data[[xvar]], y = mean_bt, group = 1)) +
      geom_line(colour = "black", linewidth = 0.7) +
      geom_point(shape = 21, fill = "white", colour = "black", size = 2.3, stroke = 0.8) +
      geom_errorbar(aes(ymin = lwr_bt, ymax = upr_bt), width = 0.12, colour = "black", linewidth = 0.6) +
      facet_wrap(as.formula(paste0("~", facetvar)), nrow = 1) +
      labs(x = xvar, y = paste0(resp, " (back-transformed)")) +
      theme_pub(11)

    save_plot(g, fname, w = 11.5, h = 4.2)
    invisible(g)
  }

  figs_dir <- file.path(out_dir, "figures", "optionA_factor")
  plot_emm(emms |> filter(label == "thin_by_depth"),  "thin_f",  "depth_f", file.path(figs_dir, paste0("A_", resp, "_thin_by_depth")))
  plot_emm(emms |> filter(label == "prune_by_depth"), "prune_f", "depth_f", file.path(figs_dir, paste0("A_", resp, "_prune_by_depth")))
  plot_emm(emms |> filter(label == "thin_by_age"),    "thin_f",  "age_f",   file.path(figs_dir, paste0("A_", resp, "_thin_by_age")))
  plot_emm(emms |> filter(label == "prune_by_age"),   "prune_f", "age_f",   file.path(figs_dir, paste0("A_", resp, "_prune_by_age")))

  save_diag(mod, file.path(out_dir, "figures", "diagnostics", paste0("diag_optionA_", resp, ".png")))

  list(anova = anov, r2 = r2v, emms = emms)
}

run_optionB <- function(resp, data) {
  pl <- prep_log(data, resp)
  d <- pl$d
  eps <- pl$eps

  fml <- resp_log ~
    age_f * depth_f +
    thin10 * prune10 +
    age_f:thin10 + age_f:prune10 +
    depth_f:thin10 + depth_f:prune10 +
    (1 | unit)

  fit <- fit_with_fallback(fml, d)
  mod <- fit$model
  model_tag <- paste0("OptionB_intensity_", fit$model_class)

  anov <- extract_type3(mod, element = resp, model_tag = model_tag)
  r2v  <- safe_r2_tbl(mod) |> mutate(element = resp, model_tag = model_tag)

  trends_one <- function(by, var, label) {
    tr <- if (inherits(mod, "lmerMod")) {
      emtrends(mod, specs = by, var = var, lmer.df = "satterthwaite")
    } else {
      emtrends(mod, specs = by, var = var)
    }

    s <- as_tibble(summary(tr, infer = c(TRUE, TRUE), adjust = "none"))
    slope_col <- get_trend_col(names(s))
    ci <- get_ci_cols(names(s))

    lwr <- if (!is.na(ci$lwr)) s[[ci$lwr]] else rep(NA_real_, nrow(s))
    upr <- if (!is.na(ci$upr)) s[[ci$upr]] else rep(NA_real_, nrow(s))
    slope <- s[[slope_col]]

    s |>
      mutate(
        element = resp,
        model_tag = model_tag,
        label = label,
        effect_per10pct = 100 * (exp(slope) - 1),
        lwr_per10pct    = 100 * (exp(lwr)   - 1),
        upr_per10pct    = 100 * (exp(upr)   - 1),
        p_fmt = p_fmt(p.value),
        sig = sig_stars(p.value)
      )
  }

  trends <- bind_rows(
    trends_one(~ age_f,   "thin10",  "slope_thin_by_age"),
    trends_one(~ depth_f, "thin10",  "slope_thin_by_depth"),
    trends_one(~ age_f,   "prune10", "slope_prune_by_age"),
    trends_one(~ depth_f, "prune10", "slope_prune_by_depth")
  )

  make_pred <- function(vary = c("thin10","prune10"), hold_other = 0) {
    vary <- match.arg(vary)
    seqx <- seq(min(d[[vary]], na.rm = TRUE), max(d[[vary]], na.rm = TRUE), length.out = 60)

    grid <- expand.grid(
      age_f = levels(d$age_f),
      depth_f = levels(d$depth_f),
      thin10 = 0,
      prune10 = 0
    )

    grid <- grid[rep(seq_len(nrow(grid)), each = length(seqx)), ]
    grid[[vary]] <- rep(seqx, times = n_distinct(grid$age_f) * n_distinct(grid$depth_f))
    if (vary == "thin10")  grid$prune10 <- hold_other
    if (vary == "prune10") grid$thin10  <- hold_other

    pred <- if (inherits(mod, "lmerMod")) {
      predict(mod, newdata = grid, re.form = NA, allow.new.levels = TRUE)
    } else {
      predict(mod, newdata = grid)
    }

    grid |>
      mutate(
        element = resp,
        model_tag = model_tag,
        vary = vary,
        pred_log = pred,
        pred_bt = exp(pred_log) - eps,
        x = .data[[vary]]
      )
  }

  preds <- bind_rows(make_pred("thin10", 0), make_pred("prune10", 0))

  lt_pool <- c("solid","dashed","dotted","dotdash","twodash","longdash")
  age_levels <- levels(d$age_f)
  lt_vals <- setNames(rep(lt_pool, length.out = length(age_levels)), age_levels)

  plot_pred <- function(df, vary_var, fname) {
    dd <- df |> filter(vary == vary_var) |> arrange(depth_f, age_f, x)

    xlab <- if (vary_var == "thin10") {
      "Thinning intensity (centred; 1 unit = +10%)"
    } else {
      "Pruning intensity (centred; 1 unit = +10%)"
    }

    g <- ggplot(dd, aes(x = x, y = pred_bt, linetype = age_f, group = age_f)) +
      geom_line(colour = "black", linewidth = 0.85) +
      facet_grid(age_f ~ depth_f) +
      scale_linetype_manual(values = lt_vals, name = "Age class") +
      labs(x = xlab, y = paste0(resp, " (back-transformed)")) +
      theme_pub(10) +
      theme(legend.position = "top")

    save_plot(g, fname, w = 11.5, h = 7.8)
    invisible(g)
  }

  figs_dir <- file.path(out_dir, "figures", "optionB_intensity")
  plot_pred(preds, "thin10",  file.path(figs_dir, paste0("B_", resp, "_thin_ageDepth")))
  plot_pred(preds, "prune10", file.path(figs_dir, paste0("B_", resp, "_prune_ageDepth")))

  save_diag(mod, file.path(out_dir, "figures", "diagnostics", paste0("diag_optionB_", resp, ".png")))

  list(anova = anov, r2 = r2v, trends = trends, preds = preds)
}

resA <- list()
resB <- list()

if (isTRUE(runA)) {
  resA <- lapply(elements, run_optionA, data = dat)
  names(resA) <- elements
}

if (isTRUE(runB)) {
  resB <- lapply(elements, run_optionB, data = dat)
  names(resB) <- elements
}

if (isTRUE(runA)) {
  anovaA <- bind_rows(lapply(resA, `[[`, "anova"))
  r2A    <- bind_rows(lapply(resA, `[[`, "r2"))
  emmsA  <- bind_rows(lapply(resA, `[[`, "emms"))

  write_csv(anovaA, file.path(out_dir, "tables", "optionA_anova.csv"))
  write_csv(r2A,    file.path(out_dir, "tables", "optionA_r2.csv"))
  write_csv(emmsA,  file.path(out_dir, "tables", "optionA_emm_backtransformed.csv"))
}

if (isTRUE(runB)) {
  anovaB  <- bind_rows(lapply(resB, `[[`, "anova"))
  r2B     <- bind_rows(lapply(resB, `[[`, "r2"))
  trendsB <- bind_rows(lapply(resB, `[[`, "trends"))
  predsB  <- bind_rows(lapply(resB, `[[`, "preds"))

  write_csv(anovaB,  file.path(out_dir, "tables", "optionB_anova.csv"))
  write_csv(r2B,     file.path(out_dir, "tables", "optionB_r2.csv"))
  write_csv(trendsB, file.path(out_dir, "tables", "optionB_emtrends_per10pct.csv"))
  write_csv(predsB,  file.path(out_dir, "tables", "optionB_predictions_for_plots.csv"))
}

readme <- c(
  "Objective 3 moderation (H2–H3)",
  "",
  "Option A (factor moderation):",
  "  log(y+eps) ~ age*depth + thin*prune + age:thin + age:prune + depth:thin + depth:prune + (1|unit)",
  "  tables/optionA_anova.csv, tables/optionA_r2.csv, tables/optionA_emm_backtransformed.csv",
  "  figures/optionA_factor/ per-element EMM plots",
  "",
  "Option B (intensity moderation):",
  "  thin10/prune10 are centred and scaled so 1 unit = +10% intensity change.",
  "  log(y+eps) ~ age*depth + thin10*prune10 + age:thin10 + age:prune10 + depth:thin10 + depth:prune10 + (1|unit)",
  "  tables/optionB_anova.csv, tables/optionB_r2.csv, tables/optionB_emtrends_per10pct.csv, tables/optionB_predictions_for_plots.csv",
  "  figures/optionB_intensity/ per-element prediction plots",
  "",
  "Diagnostics:",
  "  figures/diagnostics/ residual-vs-fitted + QQ plots (Option A and Option B)."
)
writeLines(readme, con = file.path(out_dir, "README.txt"))

writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "logs", "sessionInfo.txt"))
message("Done: ", out_dir)
