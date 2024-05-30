expected_out_list <- as.list(seq(1, M))

# Testing with no multicore backend
test_that(
  "colsbm_lapply no multicore backend, default R lapply",
  {
    expect_equal(colsbm_lapply(
      X = seq(1, M), function(i) i,
      backend = "no_mc"
    ), expected_out_list)
  }
)

#  Testing with parallel backend
if (parallel_is_installed) {
  test_that(
    "colsbm_lapply with parallel backend works",
    {
      expect_equal(colsbm_lapply(
        X = seq(1, M), function(i) i,
        backend = "parallel",
        nb_cores = 2L
      ), expected_out_list)
    }
  )
}

#  Testing with future backend
if (future_is_installed) {
  library(future)
  library(future.apply)

  if (os_is_windows) {
    plan(multisession)
  } else {
    plan(multicore)
  }

  test_that(
    "colsbm_lapply with future backend works",
    {
      expect_equal(
        colsbm_lapply(
          X = seq(1, M), function(i) i,
          backend = "future",
          nb_cores = 2L
        ), expected_out_list
      )
    }
  )
}

#  Testing with bettermc backend
if (bettermc_is_installed) {
  library(bettermc)

  test_that(
    "colsbm_lapply with bettermc backend works",
    {
      expect_equal(
        colsbm_lapply(
          X = seq(1, M), function(i) i,
          backend = "bettermc",
          nb_cores = 2L
        ), expected_out_list
      )
    }
  )
}
