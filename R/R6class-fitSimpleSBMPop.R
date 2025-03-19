#' An R6 Class object, a fitted population of netowrks sbm
#' once $optimize() is done
#'
#' @import R6

fitSimpleSBMPop <- R6::R6Class(
  classname = "fitSimpleSBMPop",
  # inherit = "bmpop",
  #
  public = list(
    #' @field n A list of size M with the number of nodes per network
    n = NULL, # a vector of size M
    #' @field M Number of networks
    M = NULL, # Number of networks
    #' @field A List of incidence matrices of size `n \times n`
    A = NULL, # List of n[m]xn[m] Matrix
    #' @field mask List of M masks, indicating NAs in the matrices. 1 for NA, 0 else
    mask = NULL, # List of NA Matrices
    #' @field nb_inter A vector of length M the number of unique non NA entries
    nb_inter = NULL, # A vector of length M the number of unique non NA entries
    #' @field directed A boolean indicating if the networks are directed or not
    directed = NULL, # Boolean for network direction, Constant
    #' @field Q An integer indicating the number of blocks
    Q = NULL, # Number of clusters, Constant
    #' @field tau List of length M, variational parameters `n[m]xQ[m]` matrices
    tau = NULL, # List of length M, variational parameters n[m]xQ[m] matrices
    #' @field alpha Matrix of size QxQ, connection parameters
    alpha = NULL, # Matrix of size QxQ, connection parameters
    #' @field delta Vector of M,  density parameters with `delta[1] = 1`
    delta = NULL, # Vector of M,  density parameters with delta[1] = 1
    #' @field pi List of M vectors of size Q, the mixture parameters
    pi = NULL, # List of M vectors of size Q, the mixture parameters
    #' @field e Vector of size M, the sum of unique entries
    e = NULL, # Vector of size M, the sum of unique entries
    #' @field emqr  List of M QxQ matrix, the sum of edges between q and r in m
    emqr = NULL, # List of M QxQ matrix, the sum of edges between q and r in m
    #' @field nmqr List of M QxQ matrix, the number of entries between q and r in m
    nmqr = NULL, # list of M QxQ matrix, the number of entries between q and r in m
    #' @field pim  List of M vectors of size Q, the mixture parameters (pi_tilde)
    pim = NULL, # List of M vectors of size Q, the mixture parameters (pi_tilde)
    #' @field alpham  list of M QxQ matrix, the classic sbm parameters (alpha_tilde)
    alpham = NULL, # list of M QxQ matrix, the classic sbm parameters (alpha_tilde)
    #' @field free_mixture A boolean indicating if the model is with free
    #' mixture
    free_mixture = NULL, # A boolean
    #' @field free_density A boolean indicating if the model is with free
    #' density
    free_density = NULL, # A boolean
    #' @field weight  A vector of size M for weighted likelihood
    weight = NULL, # A vector of size M for weighted likelihood
    #' @field distribution The emission distribution, either bernoulli or
    #' poisson
    distribution = NULL, # "poisson", "bernoulli"
    #' @field Cpi  A list of matrices of size Q x M containing TRUE (1)
    #' or FALSE (0) if the cluster is represented in the network m
    Cpi = NULL,
    #' @field Calpha The corresponding support on the connectivity parameters
    #' computed with Cpi.
    Calpha = NULL,
    #' @field mloss Loss on the M step of the VEM
    mloss = Inf,
    #' @field vloss Loss on the VE step of the VEM
    vloss = NULL,
    #' @field vbound The variational bound
    vbound = NULL,
    #' @field net_id A vector containing the "ids" or names of the networks
    #' (if none given, they are set to their number in A list)
    net_id = NULL,
    #' @field df_mixture The degrees of freedom for mixture parameters pi,used
    #' to compute penalty
    df_mixture = NULL,
    #' @field df_connect The degrees of freedom for connection parameters
    #' alpha,used to compute penalty
    df_connect = NULL,
    #' @field df_density The degrees of freedom for density parameters delta,
    #' used to compute penalty
    df_density = NULL,
    #' @field logfactA A quantity used with the Poisson probability distribution
    logfactA = NULL,
    #  algo_ve = NULL,
    #  minibatch = NULL,
    #' @field init_method The initialization method used for the first clustering
    init_method = NULL,
    #  verbosity = NULL,
    #' @field penalty  The penalty computed based on the number of parameters
    penalty = NULL,
    # approx_pois = NULL,
    #' @field Z  The clusters memberships, a list of size M of two matrices : 1
    #' for rows clusters memberships and 2 for columns clusters memberships
    Z = NULL,
    #' @field map Maximum a posteriori
    map = NULL,
    #' @field map_parameters MAP params
    map_parameters = NULL,
    #' @field ICL Stores the ICL of the model
    ICL = NULL,
    #' @field penalty_clustering Unused attribute
    penalty_clustering = NULL,
    #' @field BICL Stores the BICL of the model
    BICL = NULL,
    #' @field net_clustering Unused parameter
    net_clustering = NULL,
    #' @field counter_merge A counter for the merge (backward) steps
    counter_merge = 0,
    #' @field counter_split A counter for the splitting (forward) steps
    counter_split = 0,
    #' @field fit_opts Fit parameters, used to determine the fitting method/
    fit_opts = NULL,

    #' @description
    #' Initializes the fitBipartiteSBMPop object
    #'
    #' @param A List of incidence Matrix of size `n[[2]][m]xn[[2]][m]`
    #' @param Q The number of blocks
    #' @param Z  The block memberships, a list of size M of two matrices : 1
    #' for rows clusters memberships and 2 for columns clusters memberships
    #' @param mask List of M masks, indicating NAs in the matrices. 1 for NA, 0 else
    #' @param net_id A vector containing the "ids" or names of the networks
    #' (if none given, they are set to their number in A list)
    #' @param distribution Emission distribution either : "poisson" or
    #' "bernoulli"
    #' @param free_mixture A boolean indicating if there is a free mixture
    #' @param free_density A boolean indicating if there is a free_density
    #' @param directed A boolean specifying if the networks are directed or not
    #' @param weight A vector of size M for weighted likelihood
    #' @param Cpi  A list of matrices of size Qd x M containing TRUE (1)
    #' or FALSE (0) if the d-th dimension cluster is represented
    #' in the network m
    #' @param Calpha The corresponding support on the connectivity parameters
    #' computed with Cpi.
    #' @param init_method The initialization method used for the first clustering
    #' @param greedy_exploration_starting_point Stores the coordinates Q1 & Q2
    #' from the greedy exploration to  keep track of the starting_point
    #' @param logfactA A quantity used with the Poisson probability distribution
    #' @param fit_opts Fit parameters, used to determine the fitting method/
    initialize = function(A = NULL,
                          Q = NULL,
                          Z = NULL,
                          mask = NULL,
                          net_id = NULL,
                          distribution = "bernoulli",
                          free_mixture = TRUE,
                          free_density = TRUE,
                          directed = NULL,
                          init_method = "spectral",
                          weight = NULL, # A vector of size M, the weight of each network
                          Cpi = NULL,
                          Calpha = NULL,
                          logfactA = NULL,
                          fit_opts = list(
                            algo_ve = "fp",
                            approx_pois = FALSE,
                            minibatch = TRUE,
                            verbosity = 1
                          )) {
      if (is.null(directed)) {
        directed <- !isSymmetric.matrix(A[[1]])
      }
      self$directed <- directed
      self$A <- A
      self$M <- length(A)
      self$n <- vapply(seq_along(A), function(m) nrow(A[[m]]), FUN.VALUE = .1)
      self$e <- vapply(seq_along(A), function(m) sum(A[[m]], na.rm = TRUE), FUN.VALUE = .1)
      if (!is.null(mask)) {
        self$mask <- mask
      } else {
        self$mask <- lapply(
          seq_along(self$A),
          function(m) {
            mask <- diag(1, self$n[m])
            if (sum(is.na(self$A[[m]])) > 0) {
              mask[is.na(self$A[[m]])] <- 1
            }
            mask
          }
        )
      }
      if (is.null(net_id)) {
        self$net_id <- seq(self$M)
      } else {
        self$net_id <- net_id
      }
      self$Q <- Q
      self$free_mixture <- free_mixture
      self$free_density <- free_density
      if (is.null(Cpi) | is.null(Calpha) | !self$free_mixture) {
        self$Cpi <- matrix(TRUE, self$Q, self$M)
        self$Calpha <- (Reduce("+", lapply(
          seq(self$M),
          function(m) tcrossprod(self$Cpi[, m])
        ))
        > 0)
      } else {
        self$Cpi <- Cpi
        self$Calpha <- Calpha
      }
      self$distribution <- distribution
      if (self$distribution == "poisson") {
        if (is.null(logfactA)) {
          self$logfactA <- vapply(
            seq_along(self$A),
            function(m) {
              sum(lfactorial(self$A[[m]]) * (1 - self$mask[[m]]), na.rm = TRUE)
            },
            FUN.VALUE = .1
          )
        } else {
          self$logfactA <- logfactA
        }
      }
      lapply(seq_along(self$A), function(m) self$A[[m]][is.na(self$A[[m]])] <- -1)
      self$fit_opts <- list(
        algo_ve = "fp",
        approx_pois = FALSE,
        minibatch = TRUE,
        verbosity = 1,
        max_step = 100L,
        nlopt_algo = "NLOPT_LD_MMA"
      )
      self$fit_opts <- modifyList(self$fit_opts, fit_opts)
      self$weight <- weight
      self$pi <- vector("list", self$M)
      self$pim <- vector("list", self$M)
      self$tau <- vector("list", self$M)
      self$alpham <- vector("list", self$M)
      self$vloss <- vector("list", self$M)
      self$emqr <- array(1, dim = c(self$M, self$Q, self$Q))
      self$nmqr <- array(1, dim = c(self$M, self$Q, self$Q))
      self$map <- list()
      self$map$emqr <- self$emqr
      self$map$nmqr <- self$nmqr
      self$map$pi <- self$pi
      if (self$free_mixture) {
        self$df_mixture <- colSums(self$Cpi) - 1
      } else {
        self$df_mixture <- self$Q - 1
      }
      self$df_density <- self$M - 1
      self$df_connect <- ifelse(self$directed,
        sum(self$Calpha),
        sum(self$Calpha[lower.tri(self$Calpha, diag = TRUE)])
      )
      self$Z <- if (is.null(Z)) {
        vector("list", self$M)
      } else {
        Z
      }
      self$delta <- rep(1, self$M)
      if (self$free_density) {
        self$delta <- self$e / ((self$n * (self$n - 1))) /
          (self$e[1] / ((self$n[1] * (self$n[1] - 1))))
      }
      self$alpha <- matrix(.5, Q, Q) * self$Calpha
      self$init_method <- init_method
      self$nb_inter <- self$dircoef *
        vapply(seq(self$M), function(m) self$n[m]**2 - sum(self$mask[[m]]), FUN.VALUE = .1)
      self$vbound <- vector("list", self$M)
    },
    #' Method to compute the maximum a posteriori for Z clustering
    #' @return nothing; stores the values
    compute_map = function() {
      self$Z <- lapply(
        self$tau,
        function(tau) apply(tau, 1, which.max)
      )
      invisible(self$Z)
    },
    #' Objective function
    #' @return The evaluation of the function
    objective = function() {
      sum(vapply(seq_along(self$A), function(m) vb(m)), FUN.VALUE = .1)
    },
    #' Computes the portion of the vbound with tau and alpha
    #' @param m The id of the network for which to compute
    #' @param map Wether to use the MAP parameters or not, a boolean, defaults
    #' to FALSE.
    #' @return The computed quantity.
    vb_tau_alpha = function(m, map = FALSE) {
      if (map) {
        emqr <- outer(self$Cpi[, m], self$Cpi[, m]) * self$map$emqr[m, , ]
        nmqr <- outer(self$Cpi[, m], self$Cpi[, m]) * self$map$nmqr[m, , ]
        alpha <- self$Calpha * self$map$alpha
        delta <- self$map$delta[m]
      } else {
        emqr <- outer(self$Cpi[, m], self$Cpi[, m]) * self$emqr[m, , ]
        nmqr <- outer(self$Cpi[, m], self$Cpi[, m]) * self$nmqr[m, , ]
        alpha <- self$Calpha * self$alpha
        delta <- self$delta[m]
      }
      switch(self$distribution,
        "bernoulli" = {
          self$dircoef * sum(
            self$Calpha * .xlogy(emqr, alpha * delta, eps = 1e-12) +
              self$Calpha * .xlogy(nmqr - emqr, 1 - alpha * delta, eps = 1e-12)
          )
        },
        "poisson" = {
          #          tau_tmp <- self$tau[[m]]
          self$dircoef * (sum(self$Calpha * emqr * log(pmax(alpha * delta, 1e-12))) -
            sum(nmqr * self$Calpha * alpha * delta) -
            self$logfactA[m])
        }
      )
    },
    #' @param m The id of the network for which to compute
    #' @param map Wether to use the MAP parameters or not, a boolean, defaults
    #' to FALSE.
    #' @return The computed quantity.
    vb_tau_pi = function(m, map = FALSE) {
      if (!map) {
        sum(self$tau[[m]][, which(self$Cpi[, m]), drop = FALSE] %*% log(self$pi[[m]][which(self$Cpi[, m])]))
      } else {
        sum(.xlogy(tabulate(self$Z[[m]], self$Q), self$map$pi[[m]]))
      }
    },
    #' Computes the entropy of the model
    #' @param m The id of the network for which to compute
    #' @return The computed quantity.
    entropy_tau = function(m) {
      -sum(.xlogx(self$tau[[m]][, which(self$Cpi[, m])]))
    },
    #' Objective function for the variational bound regarding
    #' the alpha and delta parameters.
    #'
    #' @param par The parameters, alpha and delta combined in one big vector.
    #' @param emqr List of M QxQ matrix, the sum of edges between q and r in m
    #' @param nmqr List of M QxQ matrix, the number of entries between q and r
    #' in m
    #'
    #' @return The evaluation of the function
    fn_vb_alpha_delta = function(par, emqr, nmqr) {
      df_connect <- ifelse(self$directed, self$Q**2, choose(self$Q + 1, 2))
      alpha <- matrix(0, self$Q, self$Q)
      if (self$directed) {
        alpha <- par[1:df_connect]
      } else {
        alpha[lower.tri(alpha, diag = TRUE)] <- par[1:df_connect]
        alpha[upper.tri(alpha)] <- alpha[lower.tri(alpha)]
      }
      delta <- c(1, par[1:self$df_density + df_connect])
      res <-
        vapply(seq_along(self$A),
          function(m) {
            -.5 * sum(emqr[m, , ] * .log(alpha * delta[m], eps = 1e-9) +
              (nmqr[m, , ] - emqr[m, , ]) * .log(1 - alpha * delta[m], eps = 1e-9))
          },
          FUN.VALUE = .1
        )
      sum(res)
    },
    #' Gradient of the objective function for the variational bound regarding
    #' the alpha and delta parameters.
    #'
    #' @param par The parameters, alpha and delta combined in one big vector.
    #' @param emqr List of M QxQ matrix, the sum of edges between q and r in m
    #' @param nmqr List of M QxQ matrix, the number of entries between q and r
    #' in m
    #'
    #' @return The evaluation of the function
    gr_vb_alpha_delta = function(par, emqr, nmqr) {
      # browser()
      df_connect <- ifelse(self$directed, self$Q**2, choose(self$Q + 1, 2))
      alpha <- matrix(0, self$Q, self$Q)
      if (self$directed) {
        alpha <- par[1:df_connect]
      } else {
        alpha[lower.tri(alpha, diag = TRUE)] <- par[1:df_connect]
        alpha[upper.tri(alpha)] <- alpha[lower.tri(alpha)]
      }
      delta <- c(1, par[1:self$df_density + df_connect])
      res_alpha <-
        lapply(
          seq_along(self$A),
          function(m) {
            emqr[m, , ] / pmax(alpha, 1e-9) -
              (nmqr[m, , ] - emqr[m, , ]) * (delta[m] / pmax(1 - alpha * delta[m], 1e-9))
          }
        )
      res_alpha <- Reduce("+", res_alpha)
      res_delta <-
        vapply(seq_along(self$A),
          function(m) {
            sum(emqr[m, , ] / (pmax(delta[m], 1e-9)) -
              (nmqr[m, , ] - emqr[m, , ]) * (delta[m] / pmax(1 - alpha * delta[m], 1e-9)))
          },
          FUN.VALUE = .1
        )
      if (self$directed) {
        res <- -c(as.vector(res_alpha), res_delta[2:self$M])
      } else {
        res <- -c(as.vector(res_alpha[lower.tri(res_alpha, diag = TRUE)]), res_delta[2:self$M])
      }
      invisible(res)
    },
    #' Constraint
    #'
    #' @param par The parameters, alpha and delta combined in one big vector.
    #' @param emqr List of M QxQ matrix, the sum of edges between q and r in m
    #' @param nmqr List of M QxQ matrix, the number of entries between q and r
    #' in m
    #'
    #' @return The evaluation of the function
    eval_g0_vb_alpha_delta = function(par, emqr, nmqr) {
      df_connect <- ifelse(self$directed, self$Q**2, choose(self$Q + 1, 2))
      as.vector(
        vapply(
          seq(self$M - 1),
          function(m) {
            c(par[1:self$df_density + df_connect])[m] *
              par[1:df_connect] - 1 + 1e-9
          },
          FUN.VALUE = rep(.1, df_connect)
        )
      )
    },
    #' Jacobian of the constraint
    #'
    #' @param par The parameters, alpha and delta combined in one big vector.
    #' @param emqr List of M QxQ matrix, the sum of edges between q and r in m
    #' @param nmqr List of M QxQ matrix, the number of entries between q and r
    #' in m
    #'
    #' @return The evaluation of the function
    eval_jac_g0_vb_alpha_delta = function(par, emqr, nmqr) {
      df_connect <- ifelse(self$directed, self$Q**2, choose(self$Q + 1, 2))
      jac_d <- aperm(
        vapply(
          seq(self$df_density),
          function(m) diag(par[m + df_connect], df_connect),
          FUN.VALUE = matrix(.1, df_connect, df_connect)
        ),
        perm = c(2, 3, 1)
      )
      dim(jac_d) <- c((self$df_density) * df_connect, df_connect)
      jac_a <- matrix(0, (self$df_density) * df_connect, self$df_density)
      for (m in seq(self$df_density)) {
        jac_a[1:(df_connect) + (m - 1) * (df_connect), m] <-
          par[1:(df_connect)]
      }
      cbind(jac_d, jac_a)
    },
    #' Updates the alpha and delta parameters
    #'
    #' @param map Wether to use the MAP parameters or not, a boolean, defaults
    #' to FALSE.
    #'
    #' @return nothing; but stores the values
    update_alpha_delta = function(map = FALSE) {
      # browser()
      # Only undirected
      df_connect <- ifelse(self$directed, self$Q**2, choose(self$Q + 1, 2))
      d <- self$delta
      a <- self$alpha
      if (!self$directed) {
        a <- a[lower.tri(a, diag = TRUE)]
      }
      a <- pmin(a * d[1], 1 - 1e-12)
      a <- pmax(a * d[1], 1e-12)
      d <- d / d[1]
      if (map) {
        emqr <- self$map$emqr
        nmqr <- self$map$nmqr
      } else {
        emqr <- self$emqr
        nmqr <- self$nmqr
      }
      hat <- nloptr::nloptr(
        x0 = c(a, d[2:self$M]),
        eval_f = self$fn_vb_alpha_delta,
        eval_grad_f = self$gr_vb_alpha_delta,
        lb = c(rep(
          10 * .Machine$double.eps,
          df_connect + self$df_density
        )),
        ub = c(
          rep(
            1 - 10 * .Machine$double.eps,
            df_connect
          ),
          rep(Inf, self$df_density)
        ),
        eval_g_ineq = self$eval_g0_vb_alpha_delta,
        eval_jac_g_ineq = self$eval_jac_g0_vb_alpha_delta,
        opts = list( # "algorithm" = "NLOPT_LD_MMA",
          "algorithm" = self$fit_opts$nlopt_algo,
          # "local_opts" = list("algorithm" = "NLOPT_LD_LBFGS",
          #                     "xtol_rel" = 1.0e-4),
          "xtol_rel" = 1.0e-4
        ),
        emqr = emqr, nmqr = nmqr
      )
      a <- matrix(0, self$Q, self$Q)
      if (self$directed) {
        a <- matrix(hat$solution[1:df_connect], self$Q, self$Q)
      } else {
        a[lower.tri(a, diag = TRUE)] <- hat$solution[1:(df_connect)]
        a[upper.tri(a)] <- a[lower.tri(a)]
      }
      d <- c(1, hat$solution[1:(self$df_density) + df_connect])
      if (map) {
        self$map$alpha <- a
        self$map$delta <- d
      } else {
        self$alpha <- a
        self$delta <- d
      }
    },
    #' Computes the variational bound (vbound)
    #'
    #' @return The variational bound for the model.
    compute_vbound = function() {
      sum(vapply(
        seq_along(self$A),
        function(m) self$entropy_tau(m) + self$vb_tau_pi(m) + self$vb_tau_alpha(m),
        FUN.VALUE = .1
      ))
    },
    #' Computes the penalty for the model
    #'
    #' @return the computed penalty using the formulae.
    compute_penalty = function() {
      if (self$free_mixture) {
        self$df_mixture <- colSums(self$Cpi) - 1
      } else {
        self$df_mixture <- self$Q - 1
      }
      self$df_density <- self$M - 1
      self$df_connect <- ifelse(self$directed,
        sum(self$Calpha),
        sum(self$Calpha[lower.tri(self$Calpha, diag = TRUE)])
      )
      df_connect <- self$df_connect
      if (self$free_density) df_connect <- self$df_connect + self$df_density
      self$penalty <- .5 * (df_connect * log(sum(self$nb_inter)) +
        ifelse(self$free_mixture,
          sum(self$df_mixture * log(self$n)),
          self$df_mixture * log(sum(self$n))
        ))
      invisible(self$penalty)
    },

    # compute_penalty_clustering = function() {
    #   df_mixture <- colSums(self$Cpi)
    #   df_connect <- sum(self$Calpha)
    #   if (self$free_density) df_connect <- self$df_connect + self$df_density
    #   self$penalty  <- .5*(df_connect*log(sum(self$nb_inter)) +
    #                          sum(self$df_mixture *log(self$n)))
    #   invisible(self$penalty)
    # },
    #' Computes the ICL criterion
    #' @param map Wether to use the MAP parameters or not, a boolean, defaults
    #' to FALSE.
    #' @return The ICL for the model.
    compute_icl = function(map = FALSE) {
      ICL <-
        sum(vapply(
          seq_along(self$A),
          function(m) self$vb_tau_pi(m, map = map) + self$vb_tau_alpha(m, map = map),
          FUN.VALUE = .1
        )) - self$compute_penalty()
      if (map) {
        self$map$ICL <- ICL
      } else {
        self$ICL <- ICL
      }
      return(ICL)
    },
    #' Computes the BICL criterion
    #' @param map Wether to use the MAP parameters or not, a boolean, defaults
    #' to FALSE.
    #' @return The ICL for the model.
    compute_BICL = function(map = TRUE) {
      # browser()
      #    Z_unique <- lapply(self$Z, function(Z) sort(unique(Z)))
      #    uz <- unique(Z_unique)
      #    self$net_clustering <- match(Z_unique, uz)
      # #   nb_netclusters <- max(self$net_clustering)
      #    pen <- 0
      #    for (g in seq_along(uz)) {
      #      Q <- length(uz[[g]]) # + self$Q - length(unique(unlist(uz)))
      #      M <- sum(self$net_clustering == g)
      #      df_mixture <- Q-1
      #      df_connect <-
      #        ifelse(self$directed, Q**2, Q*(Q+1)/2)
      #      if (self$free_density) df_connect <- df_connect + M - 1
      #      pen <- pen + .5*(df_connect*log(sum(self$nb_inter[self$net_clustering == g])) +
      #                          sum(df_mixture *log(self$n[self$net_clustering == g])))
      #    }
      #    self$penalty_clustering <- pen
      self$BICL <- self$compute_vbound() - self$compute_penalty() -
        ifelse(self$free_mixture, sum(log(choose(self$Q, colSums(self$Cpi)))) + self$M * log(self$Q), 0) #-
      #        ifelse(self$free_mixture, self$Q*self$M*log(2), 0)
      # nqr <- colSums(self$map$nmqr, dim = 1)
      # df_mixture <- vapply(seq(self$M), function(m) sum(self$map$pi[[m]] > 0),
      #                     FUN.VALUE = 1)
      # if(self$directed) {
      #   df_connect <- sum(nqr > 0)
      # } else {
      #   df_connect <- sum(nqr[upper.tri(nqr, diag = TRUE)] > 0)
      # }
      # if (self$free_density) df_connect <- df_connect + self$M - 1
      # pen <- .5*(df_connect * log(sum(self$nb_inter)) +
      #              sum(df_mixture*log(self$n)) ) + (self$M*self$Q*log(2))
      # self$BICL <-
      #   sum(vapply(
      #   seq_along(self$A),
      #   function(m) self$vb_tau_pi(m, map = map) + self$vb_tau_alpha(m, map = map),
      #   FUN.VALUE = .1)) - pen
      invisible(self$BICL)
    },
    #' Computes the exact ICL criterion
    #' @return The exact ICL for the model.
    compute_exact_icl = function() {
      ## directed not implemented yet
      df_mixture <- ifelse(self$free_mixture, (self$Q - 1), self$Q - 1)
      df_connect <-
        ifelse(self$directed, self$Q**2, self$Q * (self$Q + 1) / 2)
      if (self$free_density) df_connect <- df_connect + self$M - 1
      eqr <- colSums(self$emqr, 1)
      nqr <- colSums(self$nmqr, 1)
      if (!self$directed) {
        diag(eqr) <- diag(eqr) / 2
        diag(nqr) <- diag(nqr) / 2
        exicl <- -df_mixture * lbeta(.5, .5) +
          self$M * (lgamma(.5 * self$Q) - self$Q * lgamma(.5)) +
          sum(lbeta(.5 + eqr, .5 + nqr - eqr)[lower.tri(eqr, diag = TRUE)]) +
          Reduce(
            sum,
            lapply(
              seq_along(self$A),
              function(m) {
                sum(lgamma(self$n[m] * self$pi[[m]] + .5)) -
                  lgamma(sum(self$n[m] * self$pi[[m]] + .5))
              }
            )
          )
      } else {
        exicl <- df_mixture * (1 / lbeta(.5, .5)) +
          self$M * (lgamma(.5 * self$Q) - self$Q * lgamma(.5)) +
          sum(lbeta(.5 + eqr, .5 + nqr - eqr)) +
          Reduce(
            sum,
            lapply(
              seq_along(self$A),
              function(m) {
                sum(lgamma(self$n[m] * self$pi[[m]] + .5)) -
                  lgamma(sum(self$n[m] * self$pi[[m]] + .5))
              }
            )
          )
      }
      return(exicl)
    },
    #' Computes the exact ICL criterion for iid
    #' @return The exact ICL for the iid model.
    compute_exact_icl_iid = function() {
      ## directed not implemented yet
      df_mixture <- ifelse(self$free_mixture, self$M * (self$Q - 1), self$Q - 1)
      emqr <- self$emqr
      nmqr <- self$nmqr
      for (m in seq_along(self$A)) {
        diag(emqr[m, , ]) <- diag(emqr[m, , ]) / 2
        diag(nmqr[m, , ]) <- diag(nmqr[m, , ]) / 2
      }
      Reduce(
        sum,
        lapply(
          seq_along(self$A),
          function(m) {
            -df_mixture * lbeta(.5, .5) +
              (lgamma(.5 * self$Q) - self$Q * lgamma(.5)) +
              sum(lbeta(.5 + emqr[m, , ], .5 + nmqr[m, , ] - emqr[m, , ])[lower.tri(emqr[m, , ], diag = TRUE)]) +
              sum(lgamma(self$n[m] * self$pi[[m]] + .5)) -
              lgamma(sum(self$n[m] * self$pi[[m]] + .5))
          }
        )
      )
    },
    #' Updates the MAP parameters
    #'
    #' @return nothing; but stores the values
    update_map_parameters = function() {
      Z <- self$compute_map()
      Z <- lapply(
        seq_along(Z),
        function(m) {
          .one_hot(Z[[m]], self$Q)
        }
      )
      lapply(
        seq_along(Z),
        function(m) {
          self$map$emqr[m, , ] <- .tquadform(
            Z[[m]],
            as.matrix(self$A[[m]]) * (1 - as.matrix(self$mask[[m]]))
          )
        }
      )
      lapply(
        seq_along(Z),
        function(m) {
          self$map$nmqr[m, , ] <- .tquadform(
            Z[[m]],
            1 - as.matrix(self$mask[[m]])
          )
        }
      )
      self$map$Z <- Z
      self$map$alpha <- self$alpha
      self$map$delta <- self$delta
      self$m_step(map = TRUE)
      if (!self$free_density) {
        self$map$delta <- self$delta
      }
      lapply(seq.int(self$M), function(m) self$update_alpham(m, map = TRUE))
      invisible(Z)
    },
    #' Method to update tau values with a fixed-point algorithm
    #'
    #' @param m The number of the network in the netlist
    #' @param max_iter The maximum number of iterations to perform, defaults
    #' to 1
    #' @param tol The tolerance for which to stop iterating defaults to 1e-2
    #'
    #' @return The new tau values
    fixed_point_tau = function(m, max_iter = 1, tol = 1e-2) {
      condition <- TRUE
      it <- 0
      # reup_counter <- 0
      self$vloss[[m]] <-
        c(self$vloss[[m]], self$vb_tau_alpha(m) + self$vb_tau_pi(m) + self$entropy_tau(m))
      tau_old <- t(self$Cpi[, m] * t(self$tau[[m]]))
      Um <- 1 - as.matrix(self$mask[[m]])
      Am <- Um * as.matrix(self$A[[m]])
      while (condition) {
        tau_new <- switch(self$distribution,
          "bernoulli" = {
            tau_new <-
              t(matrix(.xlogy(self$Cpi[, m], self$pi[[m]]), self$Q, self$n[m])) +
              Am %*%
              tau_old %*%
              t(.logit(self$delta[m] * self$alpha, eps = 1e-9)) +
              Um %*%
              tau_old %*%
              t(.log(1 - self$alpha * self$delta[m], eps = 1e-9))
            if (self$directed) {
              tau_new <- tau_new +
                crossprod(
                  Am,
                  tau_old %*%
                    .logit(self$delta[m] * self$alpha, eps = 1e-9)
                ) +
                crossprod(
                  Um,
                  tau_old %*%
                    .log(1 - self$alpha * self$delta[m], eps = 1e-9)
                )
            }
            invisible(tau_new)
          },
          "poisson" = {
            tau_new <-
              t(matrix(.xlogy(self$Cpi[, m], self$pi[[m]]), self$Q, self$n[m])) +
              Am %*%
              tau_old %*%
              t(log(pmax(self$delta[m] * self$alpha, 1e-12))) -
              Um %*%
              tau_old %*%
              t(self$alpha * self$delta[m])
            if (self$directed) {
              tau_new <- tau_new +
                crossprod(
                  Am,
                  tau_old %*%
                    log(pmax(self$delta[m] * self$alpha, 1e-12))
                ) -
                crossprod(
                  Um,
                  tau_old %*%
                    self$alpha * self$delta[m]
                )
            }
            invisible(tau_new)
          }
        )
        # tau_new <- torch::nnf_softmax(tau_new, dim = 2)
        tau_new[, which(!self$Cpi[, m])] <- 0
        tau_new[, which(self$Cpi[, m])] <-
          .softmax(tau_new[, which(self$Cpi[, m]), drop = FALSE])
        tau_new[, which(self$Cpi[, m])][
          tau_new[, which(self$Cpi[, m])] < 1e-9
        ] <- 1e-9
        tau_new[, which(self$Cpi[, m])][
          tau_new[, which(self$Cpi[, m])] > 1 - 1e-9
        ] <- 1 - 1e-9
        tau_new[, which(self$Cpi[, m])] <-
          tau_new[, which(self$Cpi[, m]), drop = FALSE] /
            .rowSums(tau_new[, which(self$Cpi[, m]), drop = FALSE], self$n[m], sum(self$Cpi[, m]))
        tau_new[, which(!self$Cpi[, m])] <- 0
        it <- it + 1
        emqr <- .tquadform(tau_new, Am)
        nmqr <- .tquadform(tau_new, Um)
        emqr[is.nan(emqr)] <- 0
        nmqr[is.nan(nmqr)] <- 0
        if (it >= 1) {
          if (self$distribution == "bernoulli") {
            vl <- self$dircoef * sum(
              outer(self$Cpi[, m], self$Cpi[, m]) * (
                .xlogy(emqr, self$alpha * self$delta[m], eps = 1e-12) +
                  .xlogy(nmqr - emqr, 1 - self$alpha * self$delta[m], eps = 1e-12))
            ) +
              sum(tau_new[, which(self$Cpi[, m]), drop = FALSE] %*%
                log(self$pi[[m]][which(self$Cpi[, m])])) -
              sum(.xlogx(tau_new[which(self$Cpi[, m])]))
          }
          if (self$distribution == "poisson") {
            vl <- self$dircoef * sum(
              outer(self$Cpi[, m], self$Cpi[, m]) * (
                emqr * log(pmax(self$alpha * self$delta[m], 1e-12)) -
                  nmqr * self$alpha * self$delta[m])
            ) +
              sum(tau_new[, which(self$Cpi[, m]), drop = FALSE] %*%
                log(self$pi[[m]][which(self$Cpi[, m])])) -
              sum(.xlogx(tau_new[which(self$Cpi[, m])])) - self$logfactA[m]
          }
          self$vloss[[m]] <- c(self$vloss[[m]], vl)
          up <- self$vloss[[m]][length(self$vloss[[m]])] >
            self$vloss[[m]][length(self$vloss[[m]]) - 1]
        } else {
          up <- TRUE
        }
        # if (self$vloss[[m]][length(self$vloss[[m]])] <
        #     self$vloss[[m]][length(self$vloss[[m]]) - 1] & it == 1 & reup_counter == 0) {
        #     tau_old <- t(vapply(seq(self$n[m]),
        #                     function(i) gtools::rdirichlet(1, exp(self$tau[[m]][i])/sum(exp(self$tau[[m]][i]))),
        #                     FUN.VALUE = rep(.1, self$Q)))
        #     it == 0
        # }
        if (up) {
          self$tau[[m]] <- tau_new
          tau_old <- tau_new
        }
        condition <- mean((tau_new - tau_old)**2) > tol &
          it < max_iter & up
      }
      invisible(self$tau[[m]])
    },
    #' Fixed point to update alpha and delta
    #'
    #' @param map A boolean wether to use MAP parameters or not, defaults to
    #' FALSE
    #' @param max_iter The maximum number of iterations, default to 50
    #' @param tol The tolerance for which to stop iterating
    #'
    #' @return nothing; stores the values
    fixed_point_alpha_delta = function(map = FALSE, max_iter = 50, tol = 1e-6) {
      # switch(
      #   self$distribution,
      #   "poisson" = {
      condition <- TRUE
      d <- self$delta
      a <- self$alpha
      d_old <- d
      a_old <- a
      it <- 0
      if (map) {
        emqr <- self$map$emqr
        nmqr <- self$map$nmqr
      } else {
        emqr <- self$emqr
        nmqr <- self$nmqr
      }
      while (condition) {
        d <- rowSums(emqr, dim = 1) /
          rowSums(aperm(array(a, c(self$Q, self$Q, self$M))) * nmqr, dim = 1)
        d[1] <- 1
        a <- colSums(emqr, dim = 1) /
          colSums(array(d, c(self$M, self$Q, self$Q)) * nmqr, dim = 1)
        a[is.nan(a)] <- 0
        it <- it + 1
        condition <- mean((a - a_old)**2) + mean((d - d_old)**2) > 2 * tol &
          it < max_iter
        d_old <- d
        a_old <- a
      }
      if (map) {
        self$map$alpha <- a
        self$map$delta <- d
      } else {
        self$alpha <- a
        self$delta <- d
      }
    },
    #' Computes the pi for the whole model
    #'
    #' @param m The number of the network in the netlist
    #' @param map A boolean wether to use the MAP parameters or not, defaults to
    #' FALSE
    #'
    #' @return the pi and stores the values
    update_pi = function(m, map = FALSE) {
      self$update_pim(m, map)
      if (!map) {
        if (self$free_mixture) {
          #          pi <- self$Cpi[,m]*.colMeans(self$tau[[m]], self$n[m], self$Q)
          self$pi[[m]] <- self$pim[[m]]
        } else {
          # A changer, update juste pim a chaque m puis dire pi = sum(n*pim)/sum(n)
          pi1 <- vapply(seq(self$M),
            function(x) self$n[x] * self$pim[[x]],
            FUN.VALUE = rep(.1, self$Q)
          )
          #                        function(x) .colSums(self$tau[[x]], self$n[x], self$Q),
          #                       FUN.VALUE = rep(.1, self$Q)) #QxM matrix
          pi <- rowSums(pi1) / sum(pi1)
          #          self$pim[[m]] <- self$Cpi[,m]*.colMeans(self$tau[[m]], self$n[m], self$Q)
          self$pi <- rep(list(pi), self$M)
        }
      } else {
        if (self$free_mixture) {
          self$map$pi[[m]] <- self$map$pim[[m]]
          #          pi <- tabulate(self$Z[[m]], self$Q)/self$n[m]
        } else {
          pi1 <- vapply(seq(self$M),
            function(x) self$n[x] * self$map$pim[[x]],
            FUN.VALUE = rep(.1, self$Q)
          )
          #                        function(x) .colSums(self$tau[[x]], self$n[x], self$Q),
          #                       FUN.VALUE = rep(.1, self$Q)) #QxM matrix
          pi <- rowSums(pi1) / sum(pi1)
          #         pi <- tabulate(unlist(self$Z), self$Q)/sum(self$n)
          self$map$pi <- rep(list(pi), self$M)
        }
      }
      invisible(pi)
    },
    #' Computes the pi per network, known as the pim
    #'
    #' @param m The number of the network in the netlist
    #' @param map A boolean wether to use the MAP parameters or not, defaults to
    #' FALSE
    #'
    #' @return nothing; stores the values
    update_pim = function(m, map = FALSE) {
      if (!map) {
        pi <- self$Cpi[, m] * .colMeans(self$tau[[m]], self$n[m], self$Q)
        self$pim[[m]] <- pi
      } else {
        pi <- tabulate(self$Z[[m]], self$Q) / self$n[m]
        self$map$pim[[m]] <- pi
      }
      invisible(pi)
    },
    #' Computes the alpha per network, known as the alpham
    #'
    #' @param m The number of the network in the netlist
    #' @param map A boolean wether to use the MAP parameters or not, defaults to
    #' FALSE
    #'
    #' @return the alpham and stores the values
    update_alpham = function(m, map = FALSE) {
      if (!map) {
        alpham <- outer(self$Cpi[, m], self$Cpi[, m]) * self$emqr[m, , ] / self$nmqr[m, , ]
        alpham[is.nan(alpham)] <- 0
        self$alpham[[m]] <- alpham
      } else {
        alpham <- outer(self$Cpi[, m], self$Cpi[, m]) * self$map$emqr[m, , ] / self$map$nmqr[m, , ]
        alpham[is.nan(alpham)] <- 0
        self$map$alpham[[m]] <- alpham
      }
      invisible(alpham)
    },
    #' Computes the alpha for the whole model
    #'
    #' @param map A boolean wether to use the MAP parameters or not, defaults to
    #' FALSE
    #'
    #' @return the alpha and stores the values
    update_alpha = function(map = FALSE) {
      if (!map) {
        alpha <- self$Calpha * (colSums(self$emqr, dims = 1) / colSums(self$nmqr, dims = 1))
        alpha[is.nan(alpha)] <- 0
        self$alpha <- alpha
      } else {
        alpha <- self$Calpha * (colSums(self$map$emqr, dims = 1) / colSums(self$map$nmqr, dims = 1))
        alpha[is.nan(alpha)] <- 0
        self$map$alpha <- alpha
      }
      invisible(alpha)
    },

    #' @description Initialize clusters
    #'
    #' @importFrom gtools rdirichlet
    init_clust = function() {
      self$tau <-
        switch(self$init_method,
          "random" = lapply(
            X = seq_along(self$A),
            FUN = function(m) {
              Um <- 1 - as.matrix(self$mask[[m]])
              Am <- as.matrix(self$A[[m]]) * Um
              tau_tmp <- gtools::rdirichlet(self$n[m], rep(1, self$Q))
              self$emqr[m, , ] <- .tquadform(tau_tmp, Am)
              self$nmqr[m, , ] <- .tquadform(tau_tmp, Um)
              a <- self$emqr[m, , ] / self$nmqr[m, , ]
              prob <- self$Q * diag(a) #+ rowSums(a)
              p <- sample.int(self$Q, prob = prob)
              tau_tmp <- tau_tmp[, p]
              self$emqr[m, , ] <- self$emqr[m, p, p]
              self$nmqr[m, , ] <- self$nmqr[m, p, p]
              tau_tmp
            }
          ),
          "spectral" = lapply(
            X = seq_along(self$A),
            FUN = function(m) {
              Um <- 1 - as.matrix(self$mask[[m]])
              Am <- as.matrix(self$A[[m]]) * Um
              tau_tmp <- .one_hot(spectral_clustering(
                as.matrix(self$A[[m]]), self$Q
              ), self$Q)
              tau_tmp[tau_tmp < 1e-6] <- 1e-6
              tau_tmp[tau_tmp > 1 - 1e-6] <- 1 - 1e-6
              tau_tmp <- tau_tmp / .rowSums(tau_tmp, self$n[m], self$Q)
              self$emqr[m, , ] <- .tquadform(tau_tmp, Am)
              self$nmqr[m, , ] <- .tquadform(tau_tmp, Um)
              a <- self$emqr[m, , ] / self$nmqr[m, , ]
              prob <- self$Q * diag(a) # + rowSums(a)
              p <- sample.int(self$Q, prob = prob)
              tau_tmp <- tau_tmp[, p]
              self$emqr[m, , ] <- self$emqr[m, p, p]
              self$nmqr[m, , ] <- self$nmqr[m, p, p]
              tau_tmp
            }
          ),
          "given" = lapply(
            X = seq_along(self$Z),
            FUN = function(m) {
              Um <- 1 - as.matrix(self$mask[[m]])
              Am <- as.matrix(self$A[[m]]) * Um
              if (is.matrix(self$Z[[m]])) {
                tau_tmp <- self$Z[[m]]
                tau_tmp[, self$Cpi[, m]] <-
                  .threshold(tau_tmp[, self$Cpi[, m], drop = FALSE])
                self$emqr[m, , ] <- .tquadform(tau_tmp, Am)
                self$nmqr[m, , ] <- .tquadform(tau_tmp, Um)
                tau_tmp
              } else {
                tau_tmp <- .one_hot(self$Z[[m]], self$Q)
                tau_tmp[, self$Cpi[, m]] <-
                  .threshold(tau_tmp[, self$Cpi[, m], drop = FALSE])
                # tau_tmp[,self$Cpi[,m]][tau_tmp[,self$Cpi[,m]] < 1e-9] <- 1e-9
                # tau_tmp[,self$Cpi[,m]][tau_tmp[,self$Cpi[,m]] > 1-1e-9] <- 1-1e-9
                # tau_tmp[,self$Cpi[,m]] <- tau_tmp[,self$Cpi[,m], drop = FALSE] /
                #   .rowSums(tau_tmp[,self$Cpi[,m], drop = FALSE], self$n[m], sum(self$Cpi[,m]))
                self$emqr[m, , ] <- .tquadform(tau_tmp, Am)
                self$nmqr[m, , ] <- .tquadform(tau_tmp, Um)
                # a <- self$emqr[m,,]/self$nmqr[m,,]
                # prob <- self$Q*diag(a) #+ rowSums(a)
                # p <- sample(self$Q, prob = prob)
                # tau_tmp <- tau_tmp[,p]
                # self$emqr[m,,] <- self$emqr[m,p,p]
                # self$nmqr[m,,] <- self$nmqr[m,p,p]
                tau_tmp
              }
              tau_tmp
            }
          ),
          "empty" = lapply(
            X = seq_along(self$tau),
            FUN = function(m) {
              Um <- 1 - as.matrix(self$mask[[m]])
              Am <- as.matrix(self$A[[m]]) * Um
              tau_tmp <- self$tau[[m]]
              tau_tmp[, which(!self$Cpi[, m])] <- 0
              tau_tmp[, self$Cpi[, m]] <-
                .threshold(tau_tmp[, self$Cpi[, m], drop = FALSE])
              self$emqr[m, , ] <- .tquadform(tau_tmp, Am)
              self$nmqr[m, , ] <- .tquadform(tau_tmp, Um)
              tau_tmp
            }
          )
        )
      lapply(seq(self$M), function(m) self$update_alpham(m))
      # To have a good initialization point for gradient descent
      if (self$distribution == "bernoulli" & self$free_density &
        !self$fit_opts$approx_pois) {
        self$fixed_point_alpha_delta()
      }
    },
    #' @description TODO Remove
    make_permutation = function() {

    },
    #' @description
    #' The M step of the VEM
    #'
    #' @param map A boolean wether to use the MAP parameters or not, defaults to
    #' FALSE
    #' @param max_iter The maximum number of iterations, default to 2
    #' @param tol The tolerance for which to stop iterating defaults to 1e-3
    #' @param ... Other parameters
    #'
    #' @return nothing; stores values
    m_step = function(map = FALSE, max_iter = 100, tol = 1e-3, ...) {
      ## A faire: Separer l'update de pi du m_step, on le fait beaucou trop.
      ## Completement sous optimal lorsque M est grand.
      # browser()
      #      lapply(seq_along(self$pi), function(m) self$update_pi(m, map = map))
      if (self$free_density == FALSE) {
        self$update_alpha(map = map)
      } else {
        switch(self$distribution,
          "poisson" = self$fixed_point_alpha_delta(map = map),
          "bernoulli" =
            ifelse(self$fit_opts$approx_pois,
              self$fixed_point_alpha_delta(map = map),
              self$update_alpha_delta(map = map)
            )
        )
      }
    },
    #' An optimization version for the VE step of the VEM but currently a
    #' placeholder
    #'
    #' @param m The number of the network in the netlist
    #' @param max_iter The maximum number of iterations, default to 2
    #' @param tol The tolerance for which to stop iterating defaults to 1e-3
    #' @param ... Other parameters
    #'
    #' @return nothing; stores values
    ve_step = function(m, max_iter = 20, tol = 1e-3, ...) {
      # Place holder for gradient ascent or other optimization methods
    },
    #' Updates the mqr quantities
    #'
    #' @description
    #' Namely, it updates the emqr and nmqr.
    #'
    #' @param m The number of the network in the netlist
    update_mqr = function(m) {
      tau_tmp <- self$tau[[m]]
      Um <- 1 - as.matrix(self$mask[[m]])
      Am <- as.matrix(self$A[[m]]) * Um
      self$emqr[m, , ] <- outer(self$Cpi[, m], self$Cpi[, m]) * .tquadform(tau_tmp, Am)
      self$nmqr[m, , ] <- outer(self$Cpi[, m], self$Cpi[, m]) * .tquadform(tau_tmp, Um)
    },
    #' Perform the whole initialization and VEM algorithm
    #'
    #' @param max_step The maximum number of steps to perform optimization
    #' @param tol The tolerance for which to stop iterating, default to 1e-3
    #' @param ... Other parameters
    #'
    #' @return nothing; stores values
    optimize = function(max_step = self$fit_opts$max_step, tol = 1e-3, ...) {
      if (self$Q == 1) {
        self$tau <- lapply(seq(self$M), function(m) matrix(1, self$n[m], 1))
        self$Z <- lapply(seq(self$M), function(m) rep(1, self$n[m]))
        self$pi <- lapply(seq(self$M), function(m) 1)
        self$pim <- lapply(seq(self$M), function(m) 1)
        lapply(seq(self$M), function(m) self$update_mqr(m))
        if (self$free_density) {
          self$alpha <- matrix(sum(self$emqr[1, , ]) / sum(self$nmqr[1, , ]), 1, 1)
          self$delta <- vapply(
            seq(self$M),
            function(m) {
              sum(self$emqr[m, , ]) / sum(self$alpha * self$nmqr[m, , ])
            },
            FUN.VALUE = .1
          )
        } else {
          self$alpha <- matrix(sum(self$emqr) / sum(self$nmqr), 1, 1)
          self$delta <- rep(1, self$M)
        }
        self$alpham <- lapply(
          seq(self$M),
          function(m) matrix(sum(self$emqr[m, , ]) / sum(self$nmqr[m, , ]), 1, 1)
        )
        self$map <- list(
          Z = self$tau,
          pi = self$pi,
          alpha = self$alpha,
          delta = self$delta,
          alpham = self$alpham,
          emqr = self$emqr,
          nmqr = self$nmqr
        )
        self$vbound <- self$compute_vbound()
      } else {
        #  browser()
        self$init_clust()
        # lapply(seq(self$M), function(m) self$update_mqr(m))
        self$m_step(...)
        lapply(seq(self$M), function(m) self$update_pim(m, map = FALSE))
        if (self$free_mixture) {
          self$pi <- self$pim
        } else {
          self$update_pi(1, map = FALSE)
        }
        step <- 0
        vb <- self$compute_vbound()
        self$vbound <- vb
        step_condition <- TRUE
        while (step_condition) {
          if (!self$fit_opts$minibatch) {
            lapply(
              seq(self$M),
              function(m) {
                switch(self$fit_opts$algo_ve,
                  "fp" = self$fixed_point_tau(m),
                  self$ve_step(m, ...)
                )
                self$update_mqr(m)
              }
            )
            lapply(seq_along(self$pi), function(m) self$update_pi(m, map = FALSE))
            self$m_step(...)
          } else {
            seq_m <- sample.int(self$M)
            lapply(
              seq(self$M),
              function(m) {
                switch(self$fit_opts$algo_ve,
                  "fp" = self$fixed_point_tau(seq_m[m]),
                  self$ve_step(seq_m[m], ...)
                )
                self$update_mqr(seq_m[m])
                self$update_pi(seq_m[m], map = FALSE)
                self$m_step(...)
              }
            )
          }
          vb <- self$compute_vbound()
          self$vbound <- c(self$vbound, vb)
          step <- step + 1
          step_condition <- step < max_step &
            self$vbound[length(self$vbound)] -
              self$vbound[length(self$vbound) - 1] > tol
          if (step %% 5 == 0) {
            if (self$fit_opts$verbosity >= 1) {
              print(paste0(step, ": ", vb))
              print(self$alpha)
            }
          }
        }
        self$compute_map()
        lapply(seq(self$M), function(m) self$update_alpham(m))
        lapply(seq(self$M), function(m) self$update_pim(m, map = FALSE))
        self$update_map_parameters()
        lapply(seq(self$M), function(m) self$update_pim(m, map = TRUE))
        if (self$free_mixture) {
          self$map$pi <- self$map$pim
        } else {
          self$update_pi(1, map = TRUE)
        }
        #  self$compute_parameters()
        #    self$compute_icl()
      }
      # if (self$fit_opts$approx_pois & self$distribution == "bernoulli" & self$free_density) {
      #   fit$delta[d] <- fit$delta[d] / max(fit$delta[d]*fit$alpha[tcrossprod(Cpi[d])])
      #   vb <- self$compute_vbound()
      #   self$vbound <- c(self$vbound, vb)
      # }
      self$compute_icl()

      self$compute_icl(map = TRUE)
      self$compute_BICL()
    },
    #' The message printed when one prints the object
    #' @param type The title above the message.
    show = function(type = "Fitted Collection of Simple SBM") {
      cat(type, "--", self$distribution, "variant for", self$M, "networks \n")
      cat("=====================================================================\n")
      cat("Dimension = (", self$n, ") - (", self$Q, ") blocks.\n")
      cat("BICL = ", self$BICL, " -- #Empty blocks : ", sum(!self$Cpi), " \n")
      cat("=====================================================================\n")
      cat("* Useful fields \n")
      cat("  $distribution, $nb_nodes, $nb_clusters, $support, $Z \n")
      cat("  $memberships, $parameters, $BICL, $vbound, $pred_dyads \n")
    },
    #' The print method
    #'
    #' @return nothing; print to console
    print = function() {
      self$show()
    }
  ),
  active = list(
    #' @field dircoef The coefficients used change if the network is directed
    #' or not
    dircoef = function() ifelse(self$directed, 1, .5),
    #' @field nb_nodes Returns n a list of the number of nodes per network
    nb_nodes = function(value) self$n,
    #' @field nb_clusters Returns Q an integer with the number of blocks
    nb_clusters = function(value) self$Q,
    #' @field support Returns the Cpi, a list of M boolean matrices indicating
    #' which blocks are populated
    support = function(value) self$Cpi,
    #' @field memberships Returns the tau, the probabilities of memberships
    #' "a posteriori", after seeing the data
    memberships = function(value) self$tau,
    #' @field parameters Returns the list of parameters of the model, alpha, pi
    #' and delta
    parameters = function(value) {
      list(
        alpha = self$alpha,
        pi = self$pi,
        delta = self$delta
      )
    },
    #' @field pred_dyads Predicted dyads from the estimated probabilities and
    #' parameters
    pred_dyads = function(value) {
      lapply(
        seq(self$M),
        function(m) {
          A_hat <- .quadform(self$tau[[m]], self$delta[m] * self$alpha)
          diag(A_hat) <- 0
          A_hat
        }
      )
    }
  )
)
