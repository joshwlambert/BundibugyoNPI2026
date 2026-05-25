library(BundibugyoNPI2026)
library(ringbp)
library(data.table)
library(epiparameter)

slow_onset_to_isolation_params <- epiparameter::convert_summary_stats_to_params(
  "gamma",
  mean = 5,
  sd   = 1
)
medium_onset_to_isolation_params <- epiparameter::convert_summary_stats_to_params(
  "gamma",
  mean = 3,
  sd   = 1
)
fast_onset_to_isolation_params <- epiparameter::convert_summary_stats_to_params(
  "gamma",
  mean = 1,
  sd   = 1
)
ebola_incubation_period <- epiparameter::epiparameter_db(
  disease = "Ebola",
  epi_name = "incubation",
  single_epiparameter = TRUE
)

# Put parameters that are grouped by disease into this data.table
scenarios <- data.table(
  expand.grid(
    delay_group = list(data.table(
      delay = c("slow", "medium", "fast"),
      onset_to_isolation = c(
        \(n) rgamma(
          n= n,
          shape = slow_onset_to_isolation_params$shape,
          scale = slow_onset_to_isolation_params$scale
        ),
        \(n) rgamma(
          n= n,
          shape = medium_onset_to_isolation_params$shape,
          scale = medium_onset_to_isolation_params$scale
        ),
        \(n) rgamma(
          n= n,
          shape = fast_onset_to_isolation_params$shape,
          scale = fast_onset_to_isolation_params$scale
        )
      )
    )),
    incubation_period_group = list(data.table(
      subtype = c("Zaire2015"),
      incubation_period = c(
        as.function(ebola_incubation_period, func_type = "generate")
      )
    )),
    r0_community = c(2, 3, 4),
    r0_isolated = 0,
    disp_community = 0.27,
    disp_isolated = 1,
    prop_presymptomatic = c(0.01),
    prop_asymptomatic = c(0.01),
    initial_cases = c(1),
    quarantine = c(FALSE),
    test_sensitivity = 0.75,
    cap_max_days = 730,
    cap_cases = 10000
  )
)

list_cols <- grep("_group", colnames(scenarios), value = TRUE)
non_list_cols <- setdiff(colnames(scenarios), list_cols)

expanded_groups <- scenarios[, rbindlist(delay_group), by = c(non_list_cols)]
expanded_incub <- scenarios[, rbindlist(incubation_period_group), by = c(non_list_cols)]

scenarios <- merge(
  expanded_groups, expanded_incub, by = non_list_cols, allow.cartesian = TRUE
)

# Day on which the NPI package activates: contact tracing (gated on contact
# exposure time) and testing (gated on symptom onset time) both switch on for
# events occurring strictly after this day.
npi_start_day <- 42

symptomatic_traced <- list(
  CT0 = \(t) ifelse(t > npi_start_day, 0, 0),
  CT20 = \(t) ifelse(t > npi_start_day, 0.2, 0),
  CT40 = \(t) ifelse(t > npi_start_day, 0.4, 0),
  CT60 = \(t) ifelse(t > npi_start_day, 0.6, 0),
  CT80 = \(t) ifelse(t > npi_start_day, 0.8, 0),
  CT100 = \(t) ifelse(t > npi_start_day, 1, 0)
)

idx <- CJ(row_idx = 1:nrow(scenarios), fn_idx = 1:length(symptomatic_traced))

scenarios2 <- scenarios[idx$row_idx]
scenarios2[, symptomatic_traced_name := names(symptomatic_traced)[idx$fn_idx]]
scenarios2[, symptomatic_traced := symptomatic_traced[idx$fn_idx]]

scenarios <- scenarios2

scenarios[, scenario :=  1:.N]

scenario_sims <- scenarios[, list(data = list(.SD)), by = scenario]

n <- 100

# Run parameter sweep
scenario_sims[, sims := lapply(data, \(x, n) {
  scenario_sim_cond(
    n = n,
    initial_cases = x$initial_cases,
    offspring = offspring_opts(
      community = \(n) rnbinom(n = n, mu = x$r0_community, size = x$disp_community),
      isolated = \(n) rnbinom(n = n, mu = x$r0_isolated, size = x$disp_isolated)
    ),
    delays = delay_opts(
      incubation_period = x$incubation_period[[1]],
      onset_to_isolation = x$onset_to_isolation[[1]]
    ),
    event_probs = event_prob_opts(
      asymptomatic = x$prop_asymptomatic,
      presymptomatic_transmission = x$prop_presymptomatic,
      symptomatic_traced = x$symptomatic_traced[[1]]
    ),
    interventions = intervention_opts(
      quarantine = x$quarantine,
      test_sensitivity = local({
        sens <- x$test_sensitivity
        \(t) ifelse(t > npi_start_day, sens, 0)
      })
    ),
    sim = sim_opts(
      cap_max_days = x$cap_max_days,
      cap_cases = x$cap_cases
    ),
    cond = list(
      max_attempts = 1000,
      target_weeks = c(5, 7),
      target_cases = c(440, 2230)
    )
  )
},
n = n
)]

saveRDS(
  object = scenario_sims,
  file = file.path(
    "inst", "extdata", "BVD_parameter_sweep.rds"
  )
)

cat("Finished \n")
