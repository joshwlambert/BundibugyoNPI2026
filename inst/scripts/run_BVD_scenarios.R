library(ringbp)
library(data.table)
library(epiparameter)
library(ggplot2)

# --- Conditioning targets ----------------------------------------------------
# information on timing of outbreak start from:
# https://virological.org/t/initial-genomes-from-may-2026-bundibugyo-virus-disease-outbreak-in-the-democratic-republic-of-the-congo-and-uganda/1032/4

outbreak_start_date <- as.Date("2026-04-11")
time_since_outbreak_start <- as.numeric(Sys.Date() - outbreak_start_date)


day_lower  <- as.numeric(Sys.Date() - as.Date("2026-04-18"))    # window start (days)
day_upper  <- as.numeric(Sys.Date() - as.Date("2026-04-4"))    # window end   (days)

day_lower
day_upper

# information on estimates of outbreak size from:
# https://epiforecasts.io/BVDOutbreakSize/stable/analysis#Joint-model-estimates

case_lower <- 440   # cumulative-case lower bound (inclusive)
case_upper <- 2230  # cumulative-case upper bound (inclusive)

# Weeks whose [7*w, 7*w + 6] day span lies entirely within [day_lower, day_upper].
target_weeks <- seq(floor(day_lower / 7), floor((day_upper - 6) / 7))

# --- Loop controls -----------------------------------------------------------
n_scenarios  <- 10L      # number of qualifying outbreaks to collect
max_attempts <- 10000L   # safety cap so the loop cannot run forever
base_seed    <- 1L       # first seed; incremented each attempt for reproducibility
report_every <- 100L     # progress print frequency

# --- Branching-process parameters --------------------------------------------
k  <- 0.27   # community offspring dispersion (Polonsky et al. 2021)
r0 <- 2     # community basic reproduction number (mu of nbinom offspring)

# --- ringbp model parameters -------------------------------------------------
initial_cases <- 1

ebola_incubation_period <- epiparameter::epiparameter_db(
  disease = "Ebola",
  epi_name = "incubation",
  single_epiparameter = TRUE
)
plot(ebola_incubation_period)

onset_to_isolation_params <- epiparameter::convert_summary_stats_to_params(
  "gamma",
  mean = 3,
  sd   = 2
)
onset_to_isolation <- \(n) rgamma(
  n     = n,
  shape = onset_to_isolation_params$shape,
  scale = onset_to_isolation_params$scale
)
hist(onset_to_isolation(1e5))

delays <- delay_opts(
  incubation_period  = as.function(ebola_incubation_period, func_type = "generate"),
  onset_to_isolation = onset_to_isolation,
  latent_period      = 0
)

# "This suggests that either asymptomatic infection is occurring at very low
# numbers or not at all" from https://royalsocietypublishing.org/rstb/article/372/1721/20160303/23182

event_probs <- event_prob_opts(
  asymptomatic = 0.01,
  presymptomatic_transmission = 0.01,
  symptomatic_traced = \(t) ifelse(t > day_upper, yes = 1, no = 0)
)

interventions <- intervention_opts(
  quarantine = FALSE,
  test_sensitivity = \(t) ifelse(t > day_upper, yes = 1, no = 0)
)

sim <- sim_opts(
  cap_max_days = 730,
  cap_cases    = 10000
)

offspring <- offspring_opts(
  community = \(n) rnbinom(n = n, mu = r0, size = k),
  isolated  = \(n) rep(0, n)
)

# --- Conditioning loop -------------------------------------------------------
accepted_runs  <- vector("list", n_scenarios)
accepted_seeds <- integer(0)
n_kept         <- 0L

for (attempt in seq_len(max_attempts)) {

  seed <- base_seed + attempt - 1L
  set.seed(seed)

  run <- suppressWarnings(
    outbreak_model(
      initial_cases = initial_cases,
      offspring     = offspring,
      delays        = delays,
      event_probs   = event_probs,
      interventions = interventions,
      sim           = sim
    )
  )

  window_cumulative <- run[week %in% target_weeks, cumulative]
  condition_met <- any(
    window_cumulative >= case_lower & window_cumulative <= case_upper
  )

  if (condition_met) {
    n_kept <- n_kept + 1L
    run[, `:=`(scenario = n_kept, seed = seed)]
    accepted_runs[[n_kept]] <- run
    accepted_seeds <- c(accepted_seeds, seed)
    message(sprintf(
      "  kept scenario %d/%d on attempt %d (seed = %d, total cases = %d).",
      n_kept, n_scenarios, attempt, seed,
      as.integer(run[, max(cumulative)])
    ))
    if (n_kept == n_scenarios) break
  }

  if (attempt %% report_every == 0L) {
    message(sprintf(
      "  ...%d attempts, %d/%d scenarios kept so far.",
      attempt, n_kept, n_scenarios
    ))
  }
}

total_attempts  <- attempt
acceptance_prop <- n_kept / total_attempts
message(sprintf(
  "Conditioning acceptance: %d/%d runs met the target (%.2f%%).",
  n_kept, total_attempts, 100 * acceptance_prop
))

if (n_kept < n_scenarios) {
  stop(sprintf(
    "Only %d/%d qualifying outbreaks in %d attempts.",
    n_kept, n_scenarios, max_attempts
  ))
}

all_runs_npi <- data.table::rbindlist(accepted_runs)

# Trim each trajectory to its last week with non-zero weekly cases (drops the
# zero-padded tail that outbreak_model() always extends to floor(cap_max_days / 7)),
# and flag whether the outbreak went extinct (neither cap was binding).
cap_max_week <- floor(sim$cap_max_days / 7)
all_runs_npi <- all_runs_npi[, {
  last_idx <- max(which(weekly_cases > 0))
  out      <- .SD[seq_len(last_idx)]
  out[, extinct := cumulative[.N] < sim$cap_cases & week[.N] < cap_max_week]
  out
}, by = scenario]

# Endpoint of each extinct trajectory, for marker dots on the plots.
extinct_endpoints <- all_runs_npi[extinct == TRUE, .SD[.N], by = scenario]

# --- Per-scenario summary ----------------------------------------------------
print(all_runs_npi[, .(seed = seed[1L], total_cases = max(cumulative)),
                   by = scenario])

# --- Plots: NPI scenarios ----------------------------------------------------
plot_title <- sprintf(
  "NPIs activate after day %d", day_upper
)

# Cumulative trajectories with NPIs
ggplot(
  all_runs_npi,
  aes(x = week, y = cumulative, group = scenario, colour = factor(scenario))
) +
  geom_vline(xintercept = day_upper / 7, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(data = extinct_endpoints, size = 2.5) +
  labs(
    x     = "Week",
    y     = "Cumulative number of cases",
    title = plot_title
  ) +
  scale_y_continuous(limits = c(0, 10000)) +
  scale_x_continuous(limits = c(0, 22)) +
  theme_bw() +
  theme(legend.position = "none")

# Weekly incidence with NPIs
ggplot(
  all_runs_npi,
  aes(x = week, y = weekly_cases, group = scenario, colour = factor(scenario))
) +
  geom_vline(xintercept = day_upper / 7, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(data = extinct_endpoints, size = 2.5) +
  labs(
    x     = "Week",
    y     = "Weekly number of cases",
    title = plot_title
  ) +
  theme_minimal() +
  theme(legend.position = "none")
