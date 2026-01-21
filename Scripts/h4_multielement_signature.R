suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(vegan)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(123)

args <- commandArgs(trailingOnly = TRUE)
arg_get <- function(key, default = NULL) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}

in_file  <- arg_get("in",  "data/processed/All_elements_wide.csv")
out_base <- arg_get("out", "outputs")
perm_n   <- as.integer(arg_get("perm", "9999"))

need_pkgs <- c("dplyr","tidyr","tibble","readr","stringr","ggplot2","vegan","patchwork")
missing <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nInstall them or restore your renv environment.")
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(out_base, paste0("h4_multielement_signature_", stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

writeLines(
  c(
    paste0("run_started: ", Sys.time()),
    paste0("wd: ", getwd()),
    paste0("in_file: ", in_file),
    paste0("out_dir: ", out_dir),
    paste0("perm: ", perm_n)
  ),
  con = file.path(out_dir, "logs", "run_metadata.txt")
)

save_png <- function(filename, plot, w = 12, h = 7, dpi = 600) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename, width = w, height = h, units = "in", res = dpi, background = "white")
    print(plot)
    dev.off()
  } else {
    ggsave(filename, plot, width = w, height = h, dpi = dpi, bg = "white")
  }
}

save_pdf <- function(filename, plot, w = 12, h = 7) {
  if (capabilities("cairo")) {
    ggsave(filename, plot, width = w, height = h, device = cairo_pdf)
  } else {
    ggsave(filename, plot, width = w, height = h, device = "pdf")
  }
}

theme_pub_strong <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      text            = element_text(colour = "black", face = "bold"),
      axis.title      = element_text(size = base_size + 2, face = "bold"),
      axis.text       = element_text(size = base_size,     face = "bold"),
      strip.text      = element_text(size = base_size + 1, face = "bold"),
      legend.title    = element_text(size = base_size + 1, face = "bold"),
      legend.text     = element_text(size = base_size,     face = "bold"),
      panel.grid      = element_blank(),
      plot.margin     = margin(8, 8, 8, 8)
    )
}

dat0 <- readr::read_csv(in_file, show_col_types = FALSE)

pick_first <- function(cands, nms) {
  hit <- cands[cands %in% nms]
  if (length(hit) == 0) NA_character_ else hit[1]
}

unit_col  <- pick_first(c("unit","plot_id","plot","exp_unit","experimental_unit","eu","Unit","UNIT"), names(dat0))
age_col   <- pick_first(c("age_f","age"), names(dat0))
depth_col <- pick_first(c("depth_f","depth"), names(dat0))

thin_col  <- pick_first(c("thin","thin_f","thinning","thin_level","t_code","T"), names(dat0))
prune_col <- pick_first(c("prune","prune_f","pruning","prune_level","p_code","P"), names(dat0))

elements <- c("Fe","Cu","Mn","Zn")
has_elements <- all(elements %in% names(dat0))

has_raw <- all(c("wi_Fe","wi_Cu","wi_Mn","wi_Zn") %in% names(dat0))

if (is.na(age_col) || is.na(depth_col)) {
  stop("Input must contain age (age or age_f) and depth (depth or depth_f).")
}
if (is.na(thin_col) || is.na(prune_col)) {
  stop("Input must contain thinning/pruning identifiers (e.g., thin_f & prune_f, or t_code & p_code).")
}
if (is.na(unit_col)) {
  if (!all(c("age","t_code","p_code") %in% names(dat0))) {
    stop("No unit column found and cannot build one (need age, t_code, p_code). Provide unit/plot_id.")
  }
  dat0 <- dat0 %>% mutate(unit = interaction(age, t_code, p_code, drop = TRUE))
  unit_col <- "unit"
}

d <- dat0 %>%
  transmute(
    unit    = factor(.data[[unit_col]]),
    age_f   = factor(.data[[age_col]]),
    depth_f = factor(.data[[depth_col]]),
    thin_raw  = .data[[thin_col]],
    prune_raw = .data[[prune_col]],
    Fe_in = if ("Fe" %in% names(dat0)) Fe else NA_real_,
    Cu_in = if ("Cu" %in% names(dat0)) Cu else NA_real_,
    Mn_in = if ("Mn" %in% names(dat0)) Mn else NA_real_,
    Zn_in = if ("Zn" %in% names(dat0)) Zn else NA_real_,
    wi_Fe = if ("wi_Fe" %in% names(dat0)) wi_Fe else NA_real_,
    wi_Cu = if ("wi_Cu" %in% names(dat0)) wi_Cu else NA_real_,
    wi_Mn = if ("wi_Mn" %in% names(dat0)) wi_Mn else NA_real_,
    wi_Zn = if ("wi_Zn" %in% names(dat0)) wi_Zn else NA_real_
  )

as_tp_factor <- function(x, prefix, max_level) {
  if (is.numeric(x) || is.integer(x)) {
    return(factor(x, levels = 0:max_level, labels = paste0(prefix, 0:max_level)))
  }
  x <- as.character(x)
  if (all(grepl(paste0("^", prefix, "\\d+$"), x))) return(factor(x))
  if (suppressWarnings(all(!is.na(as.numeric(x))))) {
    xx <- as.numeric(x)
    return(factor(xx, levels = 0:max_level, labels = paste0(prefix, 0:max_level)))
  }
  factor(x)
}

thin  <- as_tp_factor(d$thin_raw,  "T", 4)
prune <- as_tp_factor(d$prune_raw, "P", 2)

d <- d %>%
  mutate(thin = thin, prune = prune) %>%
  select(-thin_raw, -prune_raw)

if (has_elements) {
  dd <- d %>%
    mutate(
      Fe = Fe_in, Cu = Cu_in, Mn = Mn_in, Zn = Zn_in
    ) %>%
    select(unit, age_f, depth_f, thin, prune, Fe, Cu, Mn, Zn)
} else if (has_raw) {
  dd <- d %>%
    group_by(unit, age_f, depth_f, thin, prune) %>%
    summarise(
      Fe = mean(wi_Fe, na.rm = TRUE),
      Cu = mean(wi_Cu, na.rm = TRUE),
      Mn = mean(wi_Mn, na.rm = TRUE),
      Zn = mean(wi_Zn, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  stop("Input must contain either Fe/Cu/Mn/Zn columns or wi_Fe/wi_Cu/wi_Mn/wi_Zn columns.")
}

dd <- dd %>%
  mutate(
    depth_f = if (all(c("0-10cm","10-20cm","20-40cm") %in% levels(factor(depth_f)))) {
      factor(depth_f, levels = c("0-10cm","10-20cm","20-40cm"))
    } else {
      factor(depth_f)
    },
    TP = interaction(thin, prune, sep = ":", drop = TRUE)
  )

stopifnot(all(elements %in% names(dd)))

writeLines(
  c(
    paste0("n_rows: ", nrow(dd)),
    paste0("n_units: ", n_distinct(dd$unit)),
    paste0("age_levels: ", paste(levels(dd$age_f), collapse = ", ")),
    paste0("depth_levels: ", paste(levels(dd$depth_f), collapse = ", ")),
    paste0("thin_levels: ", paste(levels(dd$thin), collapse = ", ")),
    paste0("prune_levels: ", paste(levels(dd$prune), collapse = ", ")),
    paste0("TP_levels: ", nlevels(dd$TP))
  ),
  con = file.path(out_dir, "logs", "data_summary.txt")
)

eps_tbl <- sapply(elements, function(v) {
  x <- dd[[v]]
  if (any(x <= 0, na.rm = TRUE)) {
    min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
    if (is.finite(min_pos)) min_pos * 0.5 else 1e-6
  } else {
    0
  }
})
write_csv(
  tibble(element = names(eps_tbl), eps = as.numeric(eps_tbl)),
  file.path(out_dir, "tables", "H4_eps_by_element.csv")
)

Xlog <- sapply(elements, function(v) log(dd[[v]] + eps_tbl[[v]]))
colnames(Xlog) <- paste0(elements, "_log")

Xz <- scale(Xlog)
colnames(Xz) <- paste0(elements, "_z")

meta <- dd %>% select(unit, age_f, depth_f, thin, prune, TP)

dist_euc <- dist(Xz, method = "euclidean")

cap_h4 <- capscale(dist_euc ~ thin * prune + Condition(age_f + depth_f), data = meta)

anova_global <- anova.cca(cap_h4, permutations = perm_n)
anova_terms  <- anova.cca(cap_h4, by = "terms", permutations = perm_n)

tab_global <- as.data.frame(anova_global) %>% tibble::rownames_to_column("term")
tab_terms  <- as.data.frame(anova_terms)  %>% tibble::rownames_to_column("term")

write_csv(tab_global, file.path(out_dir, "tables", "Table_H4-1_capscale_global.csv"))
write_csv(tab_terms,  file.path(out_dir, "tables", "Table_H4-2_capscale_terms.csv"))

perm_h4 <- adonis2(dist_euc ~ age_f + depth_f + thin * prune, data = meta, permutations = perm_n)
tab_perm <- as.data.frame(perm_h4) %>% tibble::rownames_to_column("term")
write_csv(tab_perm, file.path(out_dir, "tables", "Table_H4-3_PERMANOVA.csv"))

bd <- betadisper(dist_euc, group = meta$TP)
bd_test <- permutest(bd, permutations = perm_n)
tab_bd <- as.data.frame(bd_test$tab) %>% tibble::rownames_to_column("term")
write_csv(tab_bd, file.path(out_dir, "tables", "Table_H4-4_betadisper_TP.csv"))

depth_levels <- levels(meta$depth_f)
tab_by_depth <- bind_rows(lapply(depth_levels, function(dep) {
  idx <- which(meta$depth_f == dep)
  dist_d <- dist(Xz[idx, , drop = FALSE], method = "euclidean")
  perm_d <- adonis2(dist_d ~ age_f + thin * prune, data = meta[idx, , drop = FALSE], permutations = perm_n)
  as.data.frame(perm_d) %>%
    tibble::rownames_to_column("term") %>%
    mutate(depth_f = dep)
}))
write_csv(tab_by_depth, file.path(out_dir, "tables", "Table_H4-S1_PERMANOVA_by_depth.csv"))

sc_sites <- as.data.frame(scores(cap_h4, display = "sites", choices = 1:2))
colnames(sc_sites) <- c("CAP1","CAP2")
sc_sites <- bind_cols(meta, sc_sites)

centroids <- sc_sites %>%
  group_by(TP, thin, prune) %>%
  summarise(CAP1 = mean(CAP1), CAP2 = mean(CAP2), .groups = "drop")

p_ord <- ggplot(sc_sites, aes(CAP1, CAP2)) +
  geom_point(aes(shape = thin, colour = prune), size = 2.0, alpha = 0.65) +
  geom_point(data = centroids, aes(CAP1, CAP2), shape = 4, size = 3.2, stroke = 1.0, colour = "black") +
  labs(
    x = "Partial dbRDA axis 1 (Condition: age + depth)",
    y = "Partial dbRDA axis 2 (Condition: age + depth)",
    shape = "Thinning",
    colour = "Pruning"
  ) +
  theme_pub_strong(13) +
  theme(legend.position = "right")

save_png(file.path(out_dir, "figures", "Fig_H4-1_dbRDA_partial_management.png"), p_ord, w = 8.2, h = 5.8, dpi = 650)
save_pdf(file.path(out_dir, "figures", "Fig_H4-1_dbRDA_partial_management.pdf"), p_ord, w = 8.2, h = 5.8)

sig_tbl <- as.data.frame(Xz) %>%
  mutate(
    unit    = meta$unit,
    age_f   = meta$age_f,
    depth_f = meta$depth_f,
    thin    = meta$thin,
    prune   = meta$prune,
    TP      = meta$TP
  ) %>%
  group_by(depth_f, TP, thin, prune) %>%
  summarise(across(ends_with("_z"), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(cols = ends_with("_z"), names_to = "element", values_to = "z_mean") %>%
  mutate(
    element = sub("_z$", "", element),
    element = factor(element, levels = c("Fe","Cu","Mn","Zn"))
  )

p_heat <- ggplot(sig_tbl, aes(x = TP, y = element, fill = z_mean)) +
  geom_tile(colour = "grey40", linewidth = 0.25) +
  facet_wrap(~ depth_f, nrow = 1) +
  labs(x = "Treatment (Thin:Prune)", y = "Element", fill = "Mean z-score\n(log scale)") +
  theme_pub_strong(12) +
  theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1))

save_png(file.path(out_dir, "figures", "Fig_H4-2_signature_heatmap_TP_by_depth.png"), p_heat, w = 13.0, h = 4.0, dpi = 650)
save_pdf(file.path(out_dir, "figures", "Fig_H4-2_signature_heatmap_TP_by_depth.pdf"), p_heat, w = 13.0, h = 4.0)

p_combo <- p_ord / p_heat + plot_layout(heights = c(2.2, 1.4))
save_png(file.path(out_dir, "figures", "Fig_H4-3_combined_dbRDA_plus_heatmap.png"), p_combo, w = 13.5, h = 9.5, dpi = 700)
save_pdf(file.path(out_dir, "figures", "Fig_H4-3_combined_dbRDA_plus_heatmap.pdf"), p_combo, w = 13.5, h = 9.5)

term_p <- function(tab, term_name) {
  i <- which(tab$term == term_name)
  if (length(i) == 0) return(NA_character_)
  p <- tab$`Pr(>F)`[i[1]]
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) "<0.001" else formatC(p, format = "f", digits = 3)
}

glob_F <- tab_global$F[1]
glob_p <- tab_global$`Pr(>F)`[1]
glob_p_txt <- if (is.na(glob_p)) NA_character_ else if (glob_p < 0.001) "<0.001" else formatC(glob_p, format = "f", digits = 3)

summary_txt <- c(
  "H4: Multielement signature (Fe–Cu–Mn–Zn)",
  "",
  "Primary test: partial dbRDA (capscale) using Euclidean distance on z-standardised log-transformed elements,",
  "conditioning on stand age and soil depth (Condition(age_f + depth_f)).",
  paste0("Global capscale test: F = ", round(glob_F, 3), ", p = ", glob_p_txt, "."),
  paste0("Term tests (capscale): thin p = ", term_p(tab_terms, "thin"),
         "; prune p = ", term_p(tab_terms, "prune"),
         "; thin:prune p = ", term_p(tab_terms, "thin:prune"), "."),
  "",
  "Cross-check: PERMANOVA (adonis2) with age and depth as covariates (Table_H4-3_PERMANOVA.csv).",
  "Dispersion: betadisper by TP (Table_H4-4_betadisper_TP.csv); significant dispersion implies potential heterogeneity effects.",
  "Sensitivity: PERMANOVA repeated within each depth level (Table_H4-S1_PERMANOVA_by_depth.csv).",
  "",
  "Figures:",
  "  Fig_H4-1_dbRDA_partial_management.(png/pdf)",
  "  Fig_H4-2_signature_heatmap_TP_by_depth.(png/pdf)",
  "  Fig_H4-3_combined_dbRDA_plus_heatmap.(png/pdf)"
)

writeLines(summary_txt, con = file.path(out_dir, "H4_results_auto_summary.txt"))
writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "logs", "sessionInfo.txt"))

message("Done: ", out_dir)
