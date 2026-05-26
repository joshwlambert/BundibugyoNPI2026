library(data.table)
library(ggplot2)

# read simulation results saved by running inst/scripts/run_analysis.R
BVD_results <- readRDS(file.path("inst", "extdata", "BVD_parameter_sweep.rds"))

npi_start_day <- 42

# 3 month outbreak control horizon
control_3_month <- ceiling((npi_start_day + 92) / 7)

# 6 month outbreak control horizon
control_6_month <- ceiling((npi_start_day + 184) / 7)

# 12 month outbreak control horizon
control_12_month <- ceiling((npi_start_day + 365) / 7)

BVD_results[, pext_3_month := ringbp::extinct_prob(
  sims[[1]], extinction_week = control_3_month),
  by = scenario
]

BVD_results[, pext_6_month := ringbp::extinct_prob(
  sims[[1]], extinction_week = control_6_month),
  by = scenario
]

BVD_results[, pext_12_month := ringbp::extinct_prob(
  sims[[1]], extinction_week = control_12_month),
  by = scenario
]

BVD_results

BVD_data <- rbindlist(BVD_results$data)
BVD_data[, `:=`(
  scenario = BVD_results$scenario,
  pext_3_month = BVD_results$pext_3_month,
  pext_6_month = BVD_results$pext_6_month,
  pext_12_month = BVD_results$pext_12_month
)]

rm(BVD_results)

# convert to percentages for plotting
BVD_data[, symptomatic_traced_num := as.numeric(gsub(
  pattern = "CT", replacement = "", x = symptomatic_traced_name
))]
pext_cols <- grep("pext", names(BVD_data), value = TRUE)
BVD_data[, (pext_cols) := lapply(.SD, `*`, 100), .SDcols = pext_cols]

# tidy data to facet by outbreak control horizon
BVD_data <- melt(
  BVD_data,
  measure.vars = patterns("^pext_"),
  variable.name = "horizon",
  value.name = "pext"
)

facet_labels <- c(
  pext_3_month = "3 months",
  pext_6_month = "6 months",
  pext_12_month = "12 months",
  slow = "Mean 5 days",
  medium = "Mean 3 days",
  fast = "Mean 1 day"
)

# pull colours from Spectral palette
spectral_colours <- RColorBrewer::brewer.pal(11, "Spectral")[c(10, 4, 2)]

outbreak_control_plot <- ggplot(data = BVD_data[delay == "medium"]) +
  geom_line(
    mapping = aes(
      x = symptomatic_traced_num,
      y = pext,
      colour = as.factor(r0_community)
    ),
    linewidth = 0.75
  ) +
  geom_point(
    mapping = aes(
      x = symptomatic_traced_num,
      y = pext,
      fill = as.factor(r0_community)
    ),
    shape = 21,
    colour = "black",
    size = 3,
    stroke = 0.75
  ) +
  facet_wrap(
    facets = vars(horizon),
    labeller = as_labeller(facet_labels)
  ) +
  scale_x_continuous(
    name = "Contacts traced (%)",
    limits = c(0, 100)
  ) +
  scale_y_continuous(
    name = "Simulated outbreaks controlled (%)",
    limits = c(0, 100)
  ) +
  scale_colour_manual(values = spectral_colours) +
  scale_fill_manual(values = spectral_colours) +
  labs(
    colour = expression("Basic Reproduction Number (" * R[0] * ")"),
    fill = expression("Basic Reproduction Number (" * R[0] * ")"),
    title = "Proportion of BVD outbreaks controlled",
    subtitle = paste0(
      "Proportion of outbreaks controlled within 3, ",
      "6 or 12 months from the start of PHSMs (42 days after first case)"
    )
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21)
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box="vertical",
    strip.background = element_blank(),
    strip.text = element_text(size = 12, hjust = 0.5)
  )

ggsave(
  file.path("inst", "plots", "outbreak_control.png"),
  plot = outbreak_control_plot,
  device = "png",
  width = 150,
  height = 125,
  units = "mm",
  dpi = 300
)

outbreak_control_iso_speed_plot <- ggplot(data = BVD_data) +
  geom_line(
    mapping = aes(
      x = symptomatic_traced_num,
      y = pext,
      colour = as.factor(r0_community)
    ),
    linewidth = 0.75
  ) +
  geom_point(
    mapping = aes(
      x = symptomatic_traced_num,
      y = pext,
      fill = as.factor(r0_community)
    ),
    shape = 21,
    colour = "black",
    size = 3,
    stroke = 0.75
  ) +
  facet_grid(
    rows = vars(delay),
    cols = vars(horizon),
    labeller = as_labeller(facet_labels)
  ) +
  scale_x_continuous(
    name = "Contacts traced (%)",
    limits = c(0, 100)
  ) +
  scale_y_continuous(
    name = "Simulated outbreaks controlled (%)",
    limits = c(0, 100)
  ) +
  scale_colour_manual(values = spectral_colours) +
  scale_fill_manual(values = spectral_colours) +
  labs(
    colour = expression("Basic Reproduction Number (" * R[0] * ")"),
    fill = expression("Basic Reproduction Number (" * R[0] * ")"),
    title = "Proportion of BVD outbreaks controlled",
    subtitle = paste0(
      "Proportion of outbreaks controlled within 3, ",
      "6 or 12 months from the start of PHSMs (42 days after first case)"
    )
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21)
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box="vertical",
    strip.background = element_blank(),
    strip.text = element_text(size = 12, hjust = 0.5)
  )

ggsave(
  file.path("inst", "plots", "outbreak_control_iso_speed.png"),
  plot = outbreak_control_iso_speed_plot,
  device = "png",
  width = 250,
  height = 225,
  units = "mm",
  dpi = 300
)
