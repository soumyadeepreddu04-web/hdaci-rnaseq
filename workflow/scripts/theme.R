# Shared look for every figure: colours, ggplot theme and a save helper.
suppressPackageStartupMessages(library(ggplot2))

# Treatment colours are colour-blind checked as a set; treatment is also drawn as point shape.
# Up/down use a blue-red diverging pair with a neutral grey midpoint.
col <- list(
  treatment = c(DMSO = "#2a78d6", SAHA = "#eb6834", RMD = "#1baf7a"),
  up = "#e34948",
  down = "#2a78d6",
  neutral = "#c3c2b7",
  mid = "#f0efec",
  series = "#2a78d6",
  sequential = c("#cde2fb", "#6da7ec", "#256abf", "#0d366b"),
  ink = "#0b0b0b",
  ink2 = "#52514e",
  muted = "#898781",
  grid = "#e1e0d9",
  axis = "#c3c2b7"
)
shape_treatment <- c(DMSO = 16, SAHA = 17, RMD = 15)

theme_report <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(colour = col$ink),
      plot.title = element_text(face = "bold", size = rel(1.1)),
      plot.subtitle = element_text(colour = col$ink2, size = rel(0.95)),
      plot.caption = element_text(colour = col$muted, hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      axis.title = element_text(colour = col$ink2),
      axis.text = element_text(colour = col$ink2),
      axis.line = element_line(colour = col$axis, linewidth = 0.3),
      axis.ticks = element_line(colour = col$axis, linewidth = 0.3),
      panel.grid.major = element_line(colour = col$grid, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0, colour = col$ink),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_text(colour = col$ink2),
      legend.text = element_text(colour = col$ink2),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

save_png <- function(plot, path, width, height) {
  ggsave(path, plot, width = width, height = height, dpi = 200, bg = "white")
}

# Display names for the contrasts
contrast_label <- c(
  SAHA_vs_DMSO = "Vorinostat (SAHA) 1 µM vs DMSO",
  RMD_vs_DMSO = "Romidepsin (RMD) 15 nM vs DMSO",
  RMD_vs_SAHA = "Romidepsin vs vorinostat"
)
