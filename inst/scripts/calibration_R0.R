# calibration_R0.R
# -----------------------------------------------------------------------------
# Sweep community R0 across a grid and, for each value, run a fixed number of
# outbreak simulations (no early stopping) to estimate the proportion that meet
# the conditioning criteria used in analysis.R: cumulative cases in
# [case_lower, case_upper] at some week overlapping [day_lower, day_upper].
#
# Output: `calibration_results`, a data.table with columns r0, n_sims, n_met,
# proportion. Each R0 row is independent; seeds are reset per R0 so runs are
# reproducible.
# -----------------------------------------------------------------------------

library(ringbp)
library(data.table)
library(epiparameter)

# --- Conditioning targets (must match analysis.R) ----------------------------
day_lower  <- 21
day_upper  <- 35
case_lower <- 300
case_upper <- 1300

target_weeks <- seq(floor(day_lower / 7), floor(day_upper / 7))

# --- Calibration controls ----------------------------------------------------
r0_values    <- c(1, 1.5, 2, 2.5, 3, 3.5, 4)
n_sims       <- 10000L   # number of simulations per R0 value (no early stop)
base_seed    <- 1L       # first seed; incremented each attempt
report_every <- 500L     # progress print frequency within an R0 value

# --- Fixed ringbp model components (R0 set per-iteration below) --------------
initial_cases <- 1

ebola_incubation_period <- epiparameter::epiparameter_db(
  disease             = "Ebola",
  epi_name            = "incubation",
  single_epiparameter = TRUE
)

onset_to_isolation_params <- epiparameter::convert_summary_stats_to_params(
  "gamma",
  mean = 5,
  sd   = 2
)
onset_to_isolation <- \(n) rgamma(
  n     = n,
  shape = onset_to_isolation_params$shape,
  scale = onset_to_isolation_params$scale
)

delays <- delay_opts(
  incubation_period  = as.function(ebola_incubation_period, func_type = "generate"),
  onset_to_isolation = onset_to_isolation,
  latent_period      = 0
)

event_probs <- event_prob_opts(
  asymptomatic                = 0.25,
  presymptomatic_transmission = 0.01,
  symptomatic_ascertained     = 0.5
)

interventions <- intervention_opts(
  quarantine       = FALSE,
  test_sensitivity = 1,
  test_capacity    = \(t) if (t > day_upper) 1e5 else 0
)

sim <- sim_opts(
  cap_max_days = 100,
  cap_cases    = 10000
)

# --- Calibration loop --------------------------------------------------------
calibration_results <- data.table::data.table(
  r0         = numeric(0),
  n_sims     = integer(0),
  n_met      = integer(0),
  proportion = numeric(0)
)

for (r0 in r0_values) {

  message(sprintf("=== R0 = %.2f: running %d simulations ===", r0, n_sims))

  offspring <- offspring_opts(
    community                        = local({
      mu <- r0
      \(n) rnbinom(n = n, mu = mu, size = 0.27)
    }),
    isolated                         = \(n) rep(0, n),
    community_contact_prob_infect    = 1,
    isolated_contact_prob_infect     = 1,
    asymptomatic_contact_prob_infect = 1
  )

  n_met <- 0L

  for (attempt in seq_len(n_sims)) {

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
    if (condition_met) n_met <- n_met + 1L

    if (attempt %% report_every == 0L) {
      message(sprintf(
        "  R0=%.2f ... %d/%d attempts, %d met conditioning so far (%.2f%%).",
        r0, attempt, n_sims, n_met, 100 * n_met / attempt
      ))
    }
  }

  calibration_results <- rbind(
    calibration_results,
    data.table::data.table(
      r0         = r0,
      n_sims     = n_sims,
      n_met      = n_met,
      proportion = n_met / n_sims
    )
  )

  message(sprintf(
    "R0 = %.2f: %d/%d met conditioning (%.2f%%).",
    r0, n_met, n_sims, 100 * n_met / n_sims
  ))
}

# --- Result ------------------------------------------------------------------
print(calibration_results)

saveRDS(
  object = calibration_results,
  file = file.path("inst", "extdata", "R0_calibration.rds")
)

plot(
  calibration_results$r0,
  calibration_results$proportion,
  type = "b",
  pch  = 19,
  xlab = "Community R0",
  ylab = "Proportion of runs meeting conditioning",
  main = sprintf(
    "Calibration sweep (%d sims per R0)",
    n_sims
  )
)
