#' [ringbp::scenario_sim()] with conditioning on each outbreak
#'
#' @inheritParams ringbp::scenario_sim
#' @param cond A list including:
#'   * `max_attempts`: the maximum number of attempts to simulate each outbreak
#'   that meets conditioning before erroring. The `max_attempts` applies to
#'   each outbreak simulation iteration, so if the `scenario_sim()` runs
#'   100 replicates (`n = 100`) and `max_attempts = 1000` then there may be
#'   up to 100,000 outbreaks simulated.
#'   * `target_weeks`: an `integer` vector of weeks during which the cumulative
#'   outbreak size is checked against `target_cases`
#'   * `target_cases`: a `numeric` vector of length 2 with the lower and upper
#'   bounds (inclusive) of the cumulative outbreak size for the conditioning
#'
#' @return A `data.table`
#' @importFrom data.table rbindlist setattr
#' @export
#' @autoglobal
scenario_sim_cond <- function(n,
                              initial_cases,
                              offspring,
                              delays,
                              event_probs,
                              interventions,
                              sim,
                              cond) {

  checkmate::assert_number(n, lower = 1, finite = TRUE)
  checkmate::assert_number(initial_cases, lower = 1, finite = TRUE)
  checkmate::assert_class(offspring, "ringbp_offspring_opts")
  checkmate::assert_class(delays, "ringbp_delay_opts")
  checkmate::assert_class(event_probs, "ringbp_event_prob_opts")
  checkmate::assert_class(interventions, "ringbp_intervention_opts")
  checkmate::assert_class(sim, "ringbp_sim_opts")
  ringbp:::cross_check_opts(delays, event_probs)

  # cross-check once per scenario not per replicate
  attr(x = delays$onset_to_self_isolation, which = "cross_checked") <- TRUE

  res <- vector("list", n)
  total_attempts <- 0L
  for (sim_idx in seq_len(n)) {

    accepted <- FALSE
    for (attempt in seq_len(cond$max_attempts)) {
      total_attempts <- total_attempts + 1L

      run <- ringbp::outbreak_model(
        initial_cases = initial_cases,
        offspring = offspring,
        delays = delays,
        event_probs = event_probs,
        interventions = interventions,
        sim = sim
      )

      window_cumulative <- run[week %in% cond$target_weeks, cumulative]

      condition_met <- any(
        window_cumulative >= cond$target_cases[1] &
          window_cumulative <= cond$target_cases[2]
      )

      if (condition_met) {
        res[[sim_idx]] <- run
        accepted <- TRUE
        break
      }
    }

    if (!accepted) {
      community_r0 <- mean(offspring$community(10000L))
      n_kept <- sim_idx - 1L
      stop(
        sprintf(
          paste0(
            "Conditioning not met: kept %d of %d replicates; replicate %d ",
            "exhausted its per-replicate budget of %d attempts (cumulative ",
            "%d attempts across replicates, acceptance %.2f%%).\n",
            "  target_weeks: %s\n",
            "  target_cases: [%g, %g] (inclusive cumulative bounds)\n",
            "  community R0 (empirical mean of %d draws): %.3f\n",
            "  initial_cases: %d; cap_max_days: %d; cap_cases: %d\n",
            "Consider widening `cond$target_cases`, shifting ",
            "`cond$target_weeks`, increasing `cond$max_attempts`, or ",
            "revising the offspring / delay / event parameters."
          ),
          n_kept, n, sim_idx, cond$max_attempts, total_attempts,
          100 * n_kept / total_attempts,
          paste(cond$target_weeks, collapse = ", "),
          cond$target_cases[1], cond$target_cases[2],
          10000L, community_r0,
          as.integer(initial_cases), sim$cap_max_days, sim$cap_cases
        ),
        call. = FALSE
      )
    }
  }

  extinct <- vapply(
    res, attr, FUN.VALUE = logical(1), which = "extinct", exact = TRUE
  )

  # bind output together and add simulation index
  res <- rbindlist(res, idcol = "sim")

  setattr(res, name = "cap_cases", value = sim$cap_cases)
  setattr(res, name = "extinct", value = extinct)

  res[]
}
