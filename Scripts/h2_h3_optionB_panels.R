suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
arg_get <- function(key, default = NULL) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}

results_dir <- arg_get("results", NULL)

find_pred <- function(root) {
  hits <- list.files(root, pattern = "^optionB_predictions_for_plots\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("optionB_predictions_for_plots.csv not found under: ", root)
  hits[1]
}

pred_path <- if (!is.null(results_dir)) {
  find_pred(results_dir)
} else {
  find_pred(getwd())
}

root_dir <- dirname(dirname(pred_path))
panel_dir <- file.path(root_dir, "figures", "optionB_panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(pred_path, show_col_types = FALSE)

need <- c("element","vary","x","pred_bt","age_f","depth_f")
miss <- setdiff(need, names(df))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

df <- df %>%
  mutate(
    age_f   = factor(age_f, levels = sort(unique(age_f))),
    depth_f = factor(depth_f, levels = c("0-10cm","10-20cm","20-40cm"))
  )

elements <- c("Fe","Cu","Mn","Zn")

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(face = "bold", colour = "black"),
      axis.title = element_text(size = base_size + 2),
      axis.text  = element_text(size = base_size, colour = "black"),
      strip.text = element_text(size = base_size + 1, face = "bold", colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.35),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.45),
      legend.position = "bottom",
      legend.title = element_text(size = base_size + 1, face = "bold"),
      legend.text  = element_text(size = base_size),
      plot.tag = element_text(size = base_size + 3, face = "bold"),
      plot.margin = margin(10, 12, 10, 10)
    )
}

lt_pool <- c("solid","dashed","dotted","dotdash","twodash","longdash")
age_levels <- levels(df$age_f)
lt_vals <- setNames(rep(lt_pool, length.out = length(age_levels)), age_levels)

cb_pal <- c("#0072B2", "#D55E00", "#009E73")
col_vals <- setNames(rep(cb_pal, length.out = length(age_levels)), age_levels)

pct_label <- function(z) sprintf("%+g", z * 10)

plot_one <- function(dat, el, vary_var) {
  dd <- dat %>%
    filter(element == el, vary == vary_var) %>%
    arrange(depth_f, age_f, x)

  xlab <- if (vary_var == "prune10") {
    "Pruning intensity relative to mean (percentage points)"
  } else {
    "Thinning intensity relative to mean (percentage points)"
  }

  ggplot(dd, aes(x = x, y = pred_bt, colour = age_f, linetype = age_f, group = age_f)) +
    geom_line(linewidth = 0.95) +
    facet_wrap(~ depth_f, nrow = 1) +
    scale_colour_manual(values = col_vals, name = "Age class") +
    scale_linetype_manual(values = lt_vals, name = "Age class") +
    scale_x_continuous(labels = pct_label) +
    labs(x = xlab, y = paste0(el, " (back-transformed)")) +
    theme_pub(12)
}

p_prune <- lapply(elements, function(el) plot_one(df, el, "prune10"))
names(p_prune) <- elements
Fig3 <- (p_prune$Fe + p_prune$Cu) / (p_prune$Mn + p_prune$Zn) +
  plot_annotation(tag_levels = "A")

p_thin <- lapply(elements, function(el) plot_one(df, el, "thin10"))
names(p_thin) <- elements
Fig4 <- (p_thin$Fe + p_thin$Cu) / (p_thin$Mn + p_thin$Zn) +
  plot_annotation(tag_levels = "A")

save_png <- function(filename, plot, w = 16, h = 9.5, dpi = 600) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename, width = w, height = h, units = "in", res = dpi, background = "white")
    print(plot)
    dev.off()
  } else {
    ggsave(filename, plot, width = w, height = h, dpi = dpi, bg = "white")
  }
}

save_png(file.path(panel_dir, "fig3_optionB_prune_4elements.png"), Fig3)
ggsave(file.path(panel_dir, "fig3_optionB_prune_4elements.pdf"), Fig3, width = 16, height = 9.5)

save_png(file.path(panel_dir, "fig4_optionB_thin_4elements.png"), Fig4)
ggsave(file.path(panel_dir, "fig4_optionB_thin_4elements.pdf"), Fig4, width = 16, height = 9.5)

message("Saved: ", panel_dir)
message("Predictions: ", pred_path)
