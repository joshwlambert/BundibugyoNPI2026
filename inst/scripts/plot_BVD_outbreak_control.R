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

# Probability of control over time, used for the heatmap plot below.
# `extinct_prob` is evaluated weekly from week 7 (just after NPIs activate on
# day 42) to the simulation cap (week 104 = 730 days).
weeks_seq <- 7:floor(730 / 7)
heatmap_data <- rbindlist(lapply(weeks_seq, function(w) {
  res <- BVD_results[
    , .(pext = ringbp::extinct_prob(sims[[1]], extinction_week = w)),
    by = scenario
  ]
  res[, week := w]
  res
}))
scenario_params <- rbindlist(BVD_results$data)
scenario_params[, scenario := BVD_results$scenario]
heatmap_data <- merge(heatmap_data, scenario_params, by = "scenario")
heatmap_data[, symptomatic_traced_num := as.numeric(gsub(
  pattern = "CT", replacement = "", x = symptomatic_traced_name
))]
heatmap_data[, pext := pext * 100]

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

# Wilson 95% binomial confidence intervals for the proportion of conditioned
# simulations meeting the control criterion (n = 200 simulations per scenario).
n_sims <- 200
z_975 <- qnorm(0.975)
BVD_data[, c("pext_low", "pext_high") := {
  p <- pext / 100
  centre <- (p + z_975^2 / (2 * n_sims)) / (1 + z_975^2 / n_sims)
  margin <- (z_975 * sqrt((p * (1 - p) + z_975^2 / (4 * n_sims)) / n_sims)) /
    (1 + z_975^2 / n_sims)
  list((centre - margin) * 100, (centre + margin) * 100)
}]

facet_labels <- c(
  pext_3_month = "3 months",
  pext_6_month = "6 months",
  pext_12_month = "12 months",
  slow = "Mean 5 days",
  medium = "Mean 3 days",
  fast = "Mean 1 day"
)

# pull colours from Spectral palette
spectral_colours <- RColorBrewer::brewer.pal(11, "Spectral")[c(10, 2)]

outbreak_control_plot <- ggplot(
  data = BVD_data[delay == "medium" & quarantine == FALSE & test_sensitivity == 0.75]
) +
  geom_ribbon(
    mapping = aes(
      x = symptomatic_traced_num,
      ymin = pext_low,
      ymax = pext_high,
      fill = as.factor(r0_community)
    ),
    alpha = 0.2,
    colour = NA
  ) +
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
    fill = guide_legend(override.aes = list(shape = 21)),
    shape = guide_legend(override.aes = list(fill = "grey50"))
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box="vertical",
    strip.background = element_blank(),
    strip.text = element_text(size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 8)
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

outbreak_control_iso_speed_plot <- ggplot(
  data = BVD_data[quarantine == FALSE & test_sensitivity == 0.75]
) +
  geom_ribbon(
    mapping = aes(
      x = symptomatic_traced_num,
      ymin = pext_low,
      ymax = pext_high,
      fill = as.factor(r0_community)
    ),
    alpha = 0.2,
    colour = NA
  ) +
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

outbreak_control_iso_speed_sensitivity_plot <- ggplot(
  data = BVD_data[quarantine == FALSE]
) +
  geom_ribbon(
    mapping = aes(
      x = symptomatic_traced_num,
      ymin = pext_low,
      ymax = pext_high,
      fill = as.factor(r0_community),
      group = interaction(r0_community, test_sensitivity)
    ),
    alpha = 0.2,
    colour = NA
  ) +
  geom_line(
    mapping = aes(
      x = symptomatic_traced_num,
      y = pext,
      colour = as.factor(r0_community),
      linetype = as.factor(test_sensitivity),
      group = interaction(r0_community, test_sensitivity)
    ),
    linewidth = 0.75
  ) +
  geom_point(
    mapping = aes(
      x = symptomatic_traced_num,
      y = pext,
      fill = as.factor(r0_community),
      shape = as.factor(test_sensitivity)
    ),
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
  scale_shape_manual(
    values = c("0.5" = 21, "0.75" = 22, "1" = 23),
    labels = c("50%", "75%", "100%")
  ) +
  scale_linetype_manual(
    values = c("0.5" = "dotted", "0.75" = "dashed", "1" = "solid"),
    labels = c("50%", "75%", "100%")
  ) +
  labs(
    colour = expression("Basic Reproduction Number (" * R[0] * ")"),
    fill = expression("Basic Reproduction Number (" * R[0] * ")"),
    shape = "Proportion of untraced infections isolated",
    linetype = "Proportion of untraced infections isolated",
    title = "Proportion of BVD outbreaks controlled",
    subtitle = paste0(
      "Proportion of outbreaks controlled within 3, ",
      "6 or 12 months from the start of PHSMs (42 days after first case)"
    )
  ) +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21)),
    shape = guide_legend(override.aes = list(fill = "grey50"))
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box="vertical",
    strip.background = element_blank(),
    strip.text = element_text(size = 12, hjust = 0.5)
  )

ggsave(
  file.path("inst", "plots", "outbreak_control_iso_speed_sensitivity.png"),
  plot = outbreak_control_iso_speed_sensitivity_plot,
  device = "png",
  width = 250,
  height = 225,
  units = "mm",
  dpi = 300
)

# Heatmap of control probability over contact-tracing ascertainment and time
# (medium isolation speed, no quarantine; rows per community R0, columns per
# test sensitivity)
outbreak_control_heatmap_plot <- ggplot(
  heatmap_data[delay == "medium" & quarantine == FALSE]
) +
  geom_tile(
    mapping = aes(
      x = symptomatic_traced_num,
      y = week,
      fill = pext
    ),
    width = 20,
    height = 1
  ) +
  facet_grid(
    rows = vars(r0_community),
    cols = vars(test_sensitivity),
    labeller = labeller(
      r0_community = as_labeller(
        c("2" = "R[0] == 2", "3" = "R[0] == 3"),
        default = label_parsed
      ),
      test_sensitivity = as_labeller(c(
        "0.5"  = "Traced contacts isolated: 50%",
        "0.75" = "Traced contacts isolated: 75%",
        "1"    = "Traced contacts isolated: 100%"
      ))
    )
  ) +
  scale_x_continuous(
    name = "Contacts traced (%)",
    breaks = seq(0, 100, 20),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Weeks since outbreak start",
    expand = c(0, 0)
  ) +
  scale_fill_viridis_c(
    name = "Outbreaks controlled (%)",
    limits = c(0, 100)
  ) +
  labs(
    title = "Probability of BVD outbreak control over time",
    subtitle = "Mean 3-day onset-to-isolation; no quarantine; NPIs activate on day 42 (week 6)"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 8)
  )

ggsave(
  file.path("inst", "plots", "outbreak_control_heatmap.png"),
  plot = outbreak_control_heatmap_plot,
  device = "png",
  width = 250,
  height = 220,
  units = "mm",
  dpi = 300
)
