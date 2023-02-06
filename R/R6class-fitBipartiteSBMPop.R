#' An R6 Class object, a fitted population of netowrks sbm
#' once $optimize() is done
#' @noMd
#' @noRd

fitBipartiteSBMPop <- R6::R6Class(
  classname = "fitBipartiteSBMPop",
  # inherit = "bmpop",
  #
  # TODO : fix the comments to adapt to bipartite
  public = list(
    nr = NULL, # a vector of size M counting the number of rows for each matrix
    nc = NULL, # a vector of size M counting the number of columns for each matrix
    M = NULL, # Number of networks
    A = NULL, # List of incidence Matrix of size nr[m]xnc[m]
    mask = NULL, # List of NA Matrices
    nb_inter = NULL, # A vector of length M the number of unique non NA entries
    directed = NULL, # Boolean for network direction, Constant
    Q = NULL, # Number of clusters, vectors of size2
    tau = NULL, # List of size M variational parameters n[m]xQ matrices
    alpha = NULL, # Matrix of size QxQ, connection parameters
    delta = NULL, # Vector of M,  density parameters with delta[1] = 1
    pi = NULL, # List of M vectors of size Q, the mixture parameters
    e = NULL, # Vector of size M, the sum of unique entries
    emqr = NULL, # List of M QxQ matrix, the sum of edges between q and r in m, ie the edges that are observed
    nmqr = NULL, # list of M QxQ matrix, the number of entries between q and r in m, ie all the possible edges
    alpham = NULL, # list of M QxQ matrix, the classic sbm parameters
    free_mixture = NULL, # A boolean
    free_density = NULL, # A boolean
    weight = NULL, # A vector of size M for weighted likelihood
    probablityDistribution = NULL, # "poisson", "bernoulli"
    mloss = Inf, # Loss on the M step of the VEM
    vloss = NULL, # Loss on the VE step of the VEM
    vbound = NULL, # The Variational bound
    net_id = NULL, # The "ids" or names of the networks (if none given, they are set to their number in A list)
    df_mixture = NULL, # The degrees of freedom for mixture parameters pi,used to compute penalty
    df_connect = NULL, # The degrees of freedom for connection parameters alpha,used to compute penalty
    df_density = NULL, # The degrees of freedom for density parameters delta,used to compute penalty
    logfactA = NULL, # used with the Poisson probability distribution
    #  algo_ve = NULL,
    #  minibatch = NULL,
    init_method = NULL, # The initialization method used for the first clustering
    #  verbosity = NULL,
    penalty = NULL, # The penalty computed based on the number of parameters
    # approx_pois = NULL,
    Z = NULL,
    MAP = NULL, # Maximum a posteriori
    MAP_parameters = NULL,
    parameters = NULL,
    ICL = NULL,
    penalty_clustering = NULL,
    ICL_clustering = NULL,
    net_clustering = NULL,
    counter_merge = 0,
    counter_split = 0,
    fit_opts = NULL,
    initialize = function(A = NULL,
                          Q = NULL,
                          Z = NULL,
                          mask = NULL,
                          net_id = NULL,
                          probablityDistribution = "bernoulli",
                          free_mixture = TRUE,
                          free_density = TRUE,
                          directed = NULL,
                          init_method = "spectral",
                          weight = NULL, # A vector of size M, the weight of each network
                          fit_opts = list(
                            algo_ve = "fp",
                            approx_pois = TRUE,
                            minibatch = TRUE,
                            verbosity = 1
                          )) {
      # if (is.null(directed)) {
      # Bipartite networks have rectangular matrices so no symmetry
      #   directed <- ! isSymmetric.matrix(A[[1]])
      # }
      # TODO : let the user specify is the directed network represented, ie Yij = {-1,0,1}
      self$directed <- directed

      # TODO : implement sanity checks ?


      self$A <- A
      self$M <- length(A)
      self$nr <- vapply(seq_along(A), function(m) nrow(A[[m]]), FUN.VALUE = .1)
      self$nc <- vapply(seq_along(A), function(m) ncol(A[[m]]), FUN.VALUE = .1)
      self$e <- vapply(seq_along(A), function(m) sum(A[[m]], na.rm = TRUE),
        FUN.VALUE = .1
      )

      # Getting or computing the NA Mask
      if (!is.null(mask)) {
        self$mask <- mask
      } else {
        self$mask <- lapply(
          seq_along(self$A),
          function(m) {
            mask <- matrix(1, nrow(self$A), ncol(self$A))
            if (sum(is.na(self$A[[m]] > 0))) {
              mask[is.na(self$A[[m]])] <- 0
            }
            mask
          }
        )
      }

      # Setting the net_id, if none is given, it's the position in the list that's used
      if (is.null(net_id)) {
        self$net_id <- seq(self$M)
      } else {
        self$net_id <- net_id
      }

      self$Q <- Q
      self$probablityDistribution <- probablityDistribution
      if (self$probablityDistribution == "poisson") {
        self$logfactA <- vapply(
          seq_along(self$A),
          function(m) {
            sum(lfactorial(self$A[[m]]), na.rm = TRUE)
          },
          FUN.VALUE = .1
        )
      }

      # TODO: check if needed
      # Replacing NAs by -1 in the incidence matrices
      lapply(
        seq_along(self$A),
        function(m) self$A[[m]][is.na(self$A[[m]])] <- -1
      )

      # Setting default fit options
      self$fit_opts <- list(
        algo_ve = "fp",
        approx_pois = TRUE,
        minibatch = TRUE,
        verbosity = 1
      )
      # If the user provided custom fit options they are applied here
      self$fit_opts <- utils::modifyList(self$fit_opts, fit_opts)

      # Registering the colSBM model to fit
      # iid : free_mixture = F, free_density = F
      # pi colSBM : free_mixture = T, free_density = F
      # delta colSBM : free_mixture = F, free_density = T
      # pi-delta colSBM : free_mixture = T, free_density = T
      self$free_mixture <- free_mixture
      self$free_density <- free_density

      # TODO: check if needed
      self$weight <- weight

      self$pi <- vector("list", self$M)
      self$tau <- vector("list", self$M)
      self$alpham <- vector("list", self$M)
      self$vloss <- vector("list", self$M)
      self$emqr <- array(1, dim = c(self$M, self$Q[1], self$Q[2]))
      self$nmqr <- array(1, dim = c(self$M, self$Q[1], self$Q[2]))

      # Maximum a posteriori
      self$MAP <- list()
      self$MAP$emqr <- self$emqr
      self$MAP$nmqr <- self$nmqr
      self$MAP$pi <- self$pi

      # Degrees of freedom
      self$df_mixture <- self$Q - 1
      self$df_density <- self$M - 1
      self$df_connect <- self$Q[1] * self$Q[2]

      self$Z <- if (is.null(Z)) {
        vector("list", self$M)
      } else {
        Z
      }

      # Computing the density of the m networks
      self$delta <- rep(1, self$M)
      if (self$free_density) {
        self$delta <- (self$e / (self$nr * self$nc)) /
          (self$e[1] / ((self$nr[1] * self$nc[1])))
      }

      self$alpha <- matrix(.5, Q[1], Q[2])
      self$init_method <- init_method
      self$nb_inter <- self$dircoef *
        vapply(seq(self$M), function(m) sum(self$mask[[m]]), FUN.VALUE = .1)
      self$vbound <- vector("list", self$M)
    },
    compute_MAP = function() {
      self$Z <- lapply(
        self$tau,
        function(tau) {
          list(
            row = apply(tau[[1]], 1, which.max),
            col = apply(tau[[2]], 1, which.max)
          )
        }
      )
      invisible(self$Z)
    },
    objective = function() {
      sum(vapply(seq_along(self$A), function(m) vb(m)), FUN.VALUE = .1)
    },
    vb_tau_alpha = function(m, MAP = FALSE) {
      if (MAP) {
        # If we are calculating the maximum a posteriori
        emqr <- self$MAP$emqr[m, , ]
        nmqr <- self$MAP$nmqr[m, , ]
        alpha <- self$MAP$alpha
        delta <- self$MAP$delta[m]
      } else {
        emqr <- self$emqr[m, , ]
        nmqr <- self$nmqr[m, , ]
        alpha <- self$alpha
        delta <- self$delta[m]
      }
      switch(self$probablityDistribution,
        "bernoulli" = {
          # tau_tmp <- self$tau[[m]]
          # self$dircoef*sum(
          #   self$M[[m]] *
          #     (self$A[[m]] *
          #        tcrossprod(tau_tmp %*% .logit(self$alpha*self$delta[m], eps=1e-6),
          #                   tau_tmp) +
          #        tcrossprod(tau_tmp %*% .log(1-self$alpha*self$delta[m], eps = 1e-6),
          #                   tau_tmp)))
          sum(
            .xlogy(emqr, alpha * delta, eps = 1e-12) +
              .xlogy(nmqr - emqr, 1 - alpha * delta, eps = 1e-12)
          )
        },
        "poisson" = {
          #          tau_tmp <- self$tau[[m]]
          sum(.xlogy(emqr, alpha * delta, eps = 1e-12) -
            nmqr * alpha * delta -
            self$logfactA[m])
          # self$dircoef*sum(
          #   self$M[[m]] * (self$A[[m]] *
          #                    tcrossprod(tau_tmp %*% .log(self$alpha*self$delta[m], eps=1e-9),
          #                               tau_tmp) -
          #                    tcrossprod(tau_tmp %*% (self$alpha*self$delta[m]),
          #                               tau_tmp)) -
          #   lgamma(A[[m]]+1))
        }
      )
    },
    vb_tau_pi = function(m, MAP = FALSE) {
      if (!MAP) {
        sum(self$tau[[m]][[1]] %*% log(self$pi[[m]][[1]])) +
          sum(self$tau[[m]][[2]] %*% log(self$pi[[m]][[2]]))
      } else {
        sum(.xlogy(tabulate(self$Z[[m]][[1]], self$Q[1]), self$MAP$pi[[m]][[1]])) +
          sum(.xlogy(tabulate(self$Z[[m]][[2]], self$Q[2]), self$MAP$pi[[m]][[2]]))
      }
    },
    entropy_tau = function(m) {
      -sum(.xlogx(self$tau[[m]][[1]])) - sum(.xlogx(self$tau[[m]][[2]]))
    },
    fn_vb_alpha_delta = function(par, emqr, nmqr) {
      alpha <- par[1:self$df_connect]
      delta <- c(1, par[1:self$df_density + self$df_connect])
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
    gr_vb_alpha_delta = function(par, emqr, nmqr) {
      # browser()
      alpha <- par[1:self$df_connect]
      delta <- c(1, par[1:self$df_density + self$df_connect])
      res_alpha <-
        lapply(
          seq_along(self$A),
          function(m) {
            emqr[m, , ] / pmax(alpha, 1e-9) -
              (nmqr[m, , ] - emqr[m, , ]) *
                (delta[m] / pmax(1 - alpha * delta[m], 1e-9))
          }
        )
      res_alpha <- Reduce("+", res_alpha)
      res_delta <-
        vapply(seq_along(self$A),
          function(m) {
            sum(emqr[m, , ] / (pmax(delta[m], 1e-9)) -
              (nmqr[m, , ] - emqr[m, , ]) *
                (delta[m] / pmax(1 - alpha * delta[m], 1e-9)))
          },
          FUN.VALUE = .1
        )
      if (self$directed) {
        res <- -c(as.vector(res_alpha), res_delta[2:self$M])
      } else {
        res <- -c(
          as.vector(res_alpha[lower.tri(res_alpha, diag = TRUE)]),
          res_delta[2:self$M]
        )
      }
      invisible(res)
    },
    eval_g0_vb_alpha_delta = function(par, emqr, nmqr) {
      as.vector(
        vapply(
          seq(self$M - 1),
          function(m) {
            c(par[1:self$df_density + self$df_connect])[m] *
              par[1:self$df_connect] - 1 + 1e-9
          },
          FUN.VALUE = rep(.1, self$df_connect)
        )
      )
    },
    eval_jac_g0_vb_alpha_delta = function(par, emqr, nmqr) {
      jac_d <- aperm(
        vapply(
          seq(self$df_density),
          function(m) diag(par[m + self$df_connect], self$df_connect),
          FUN.VALUE = matrix(.1, self$df_connect, self$df_connect)
        ),
        perm = c(2, 3, 1)
      )
      dim(jac_d) <- c((self$df_density) * self$df_connect, self$df_connect)
      jac_a <- matrix(0, (self$df_density) * self$df_connect, self$df_density)
      for (m in seq(self$df_density)) {
        jac_a[1:(self$df_connect) + (m - 1) * (self$df_connect), m] <-
          par[1:(self$df_connect)]
      }
      cbind(jac_d, jac_a)
    },
    update_alpha_delta = function(MAP = FALSE) {
      # browser()
      d <- self$delta
      a <- self$alpha
      a <- pmin(a * d[1], 1 - 1e-12)
      d <- d / d[1]
      if (MAP) {
        emqr <- self$MAP$emqr
        nmqr <- self$MAP$nmqr
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
          self$df_connect + self$df_density
        )),
        ub = c(
          rep(
            1 - 10 * .Machine$double.eps,
            self$df_connect
          ),
          rep(Inf, self$df_density)
        ),
        eval_g_ineq = self$eval_g0_vb_alpha_delta,
        eval_jac_g_ineq = self$eval_jac_g0_vb_alpha_delta,
        opts = list(
          "algorithm" = "NLOPT_LD_MMA",
          # "local_opts" = list("algorithm" = "NLOPT_LD_LBFGS",
          #                     "xtol_rel" = 1.0e-4),
          "xtol_rel" = 1.0e-4
        ),
        emqr = emqr, nmqr = nmqr
      )
      a <- matrix(hat$solution[1:self$df_connect], self$Q[1], self$Q[2])
      d <- c(1, hat$solution[1:(self$df_density) + self$df_connect])
      if (MAP) {
        self$MAP$alpha <- a
        self$MAP$delta <- d
      } else {
        self$alpha <- a
        self$delta <- d
      }
    },
    compute_vbound = function() {
      sum(vapply(
        seq_along(self$A),
        function(m) self$entropy_tau(m) + self$vb_tau_pi(m) + self$vb_tau_alpha(m),
        FUN.VALUE = .1
      ))
    },
    compute_penalty = function() {
      df_connect <- self$df_connect
      if (self$free_density) df_connect <- self$df_connect + self$df_density
      self$penalty <- .5 * (df_connect * log(sum(self$nb_inter)) +
        sum(self$df_mixture * log(self$n)))
      invisible(self$penalty)
    },
    compute_icl = function(MAP = FALSE) {
      ICL <-
        sum(vapply(
          seq_along(self$A),
          function(m) self$vb_tau_pi(m, MAP = MAP) + self$vb_tau_alpha(m, MAP = MAP),
          FUN.VALUE = .1
        )) - self$compute_penalty()
      if (MAP) {
        self$MAP$ICL <- ICL
      } else {
        self$ICL <- ICL
      }
      return(ICL)
    },
    ################################################################################
    ## A modifier
    compute_icl_clustering = function(MAP = TRUE) {
      # browser()
      Z_unique <- lapply(self$Z, function(Z) sort(unique(Z[[1]])))
      uz <- unique(Z_unique)
      self$net_clustering <- match(Z_unique, uz)
      #   nb_netclusters <- max(self$net_clustering)
      pen <- 0
      for (g in seq_along(uz)) {
        Q <- length(uz[[g]]) # + self$Q - length(unique(unlist(uz)))
        M <- sum(self$net_clustering == g)
        df_mixture <- Q - 1
        df_connect <- Q * Q
        if (self$free_density) df_connect <- df_connect + M - 1
        pen <- pen + .5 * (df_connect * log(sum(self$nb_inter[self$net_clustering == g])) +
          sum(df_mixture * log(self$n[self$net_clustering == g])))
      }
      self$penalty_clustering <- pen
      self$ICL_clustering <-
        sum(vapply(
          seq_along(self$A),
          function(m) self$vb_tau_pi(m, MAP = MAP) + self$vb_tau_alpha(m, MAP = MAP),
          FUN.VALUE = .1
        )) - pen
      invisible(self$ICL_clustering)
    },
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
    ################################################################################
    update_MAP_parameters = function() {
      Z <- self$compute_MAP()
      ZR <- lapply(
        seq_along(Z),
        function(m) {
          .one_hot(Z[[m]][[1]], self$Q[1])
        }
      )
      ZC <- lapply(
        seq_along(Z),
        function(m) {
          .one_hot(Z[[m]][[2]], self$Q[2])
        }
      )
      lapply(
        seq_along(Z),
        function(m) {
          self$MAP$emqr[m, , ] <-
            crossprod(ZR[[m]], (self$A[[m]] * self$mask[[m]]) %*% ZC[[m]])
        }
      )
      lapply(
        seq_along(Z),
        function(m) {
          self$MAP$nmqr[m, , ] <-
            crossprod(ZR[[m]], self$mask[[m]] %*% ZC[[m]])
        }
      )
      self$MAP$Z <- Z
      self$MAP$alpha <- self$alpha
      self$MAP$delta <- self$delta
      self$m_step(MAP = TRUE)
      if (!self$free_density) {
        self$MAP$delta <- self$delta
      }
      lapply(seq.int(self$M), function(m) self$update_alpham(m, MAP = TRUE))
      invisible(Z)
    },

    # The fixed point algorithm to update the tau
    fixed_point_tau = function(m, d, max_iter = 10, tol = 1e-3) {
      condition <- TRUE
      it <- 0
      # reup_counter <- 0
      self$vloss[[m]] <-
        c(self$vloss[[m]], self$vb_tau_alpha(m) + self$vb_tau_pi(m) +
          self$entropy_tau(m))
      tau_old <- self$tau[[m]][[d]]
      tau_new <- switch(self$probablityDistribution,
        "bernoulli" = {
          tau_new <-
            if (d == 1) {
              tau_new <-
                t(matrix(log(self$pi[[m]][[d]]), self$Q[d], self$nr[m])) +
                (self$mask[[m]] * self$A[[m]]) %*%
                self$tau[[m]][[2]] %*%
                t(.logit(self$delta[m] * self$alpha, eps = 1e-9)) +
                self$mask[[m]] %*%
                self$tau[[m]][[2]] %*%
                t(.log(1 - self$alpha * self$delta[m], eps = 1e-9))
            }
          if (d == 2) {
            tau_new <-
              t(matrix(log(self$pi[[m]][[d]]), self$Q[d], self$nc[m])) +
              (self$mask[[m]] * self$A[[m]]) %*% # FIXME : check t(self$mask[[m]] * self$A[[m]])
              self$tau[[m]][[1]] %*%
              .logit(self$delta[m] * self$alpha, eps = 1e-9) +
              self$mask[[m]] %*% # FIXME : check t(self$mask[[m]])
              self$tau[[m]][[1]] %*%
              .log(1 - self$alpha * self$delta[m], eps = 1e-9)
          }
          invisible(tau_new)
        },
        "poisson" = {
          if (d == 1) {
            tau_new <-
              t(matrix(log(self$pi[[m]][[d]]), self$Q[d], self$nr[m])) +
              (self$mask[[m]] * self$A[[m]]) %*%
              self$tau[[m]][[2]] %*%
              t(log(self$delta[m] * self$alpha)) -
              self$mask[[m]] %*%
              self$tau[[m]][[2]] %*%
              t(self$alpha * self$delta[m])
          }
          if (d == 2) {
            tau_new <-
              t(matrix(log(self$pi[[m]][[d]]), self$Q[d], self$nc[m])) +
              (self$mask[[m]] * self$A[[m]]) %*% # FIXME : check shouldn't it be t(self$mask[[m]] * self$A[[m]]) %*%
              self$tau[[m]][[1]] %*%
              t(log(self$delta[m] * self$alpha)) - # FIXME : check, just log(self$delta[m] * self$alpha)
              self$mask[[m]] %*% # FIXME : check, t(self$mask[[m]])
              self$tau[[m]][[1]] %*%
              t(self$alpha * self$delta[m]) # FIXME : check, (self$alpha$ * self$delta[m])
          }
          invisible(tau_new)
        }
      )
      tau_new <- .softmax(tau_new)
      tau_new[tau_new < 1e-9] <- 1e-9
      tau_new[tau_new > 1 - 1e-9] <- 1 - 1e-9
      tau_new <- tau_new / rowSums(tau_new)
      self$tau[[m]][[d]] <- tau_new
      invisible(tau_new)
    },
    fixed_point_alpha_delta = function(MAP = FALSE, max_iter = 50, tol = 1e-6) {
      # switch(
      #   self$probablityDistribution,
      #   "poisson" = {
      condition <- TRUE
      d <- self$delta
      a <- self$alpha
      d_old <- d
      a_old <- a
      it <- 0
      if (MAP) {
        emqr <- self$MAP$emqr
        nmqr <- self$MAP$nmqr
      } else {
        emqr <- self$emqr
        nmqr <- self$nmqr
      }
      while (condition) {
        d <- rowSums(emqr, dim = 1) /
          rowSums(aperm(array(a, c(self$Q[1], self$Q[2], self$M))) * nmqr, dim = 1)
        d[1] <- 1
        a <- colSums(emqr, dim = 1) /
          colSums(array(d, c(self$M, self$Q[1], self$Q[2])) * nmqr, dim = 1)
        a[is.nan(a)] <- 0
        it <- it + 1
        condition <- mean((a - a_old)**2) + mean((d - d_old)**2) > 2 * tol &
          it < max_iter
        d_old <- d
        a_old <- a
      }
      if (MAP) {
        self$MAP$alpha <- a
        self$MAP$delta <- d
      } else {
        self$alpha <- a
        self$delta <- d
      }
    },
    update_pi = function(m, MAP = FALSE) {
      if (!MAP) {
        pi1 <- .colMeans(self$tau[[m]][[1]], self$nr[m], self$Q[1])
        pi2 <- .colMeans(self$tau[[m]][[2]], self$nr[m], self$Q[2])
        self$pi[[m]] <- list(row = pi1, col = pi2)
      } else {
        pi1 <- tabulate(self$Z[[m]][[1]], self$Q[1]) / self$nr[m]
        pi2 <- tabulate(self$Z[[m]][[2]], self$Q[2]) / self$nc[m]
        self$MAP$pi[[m]] <- list(row = pi1, col = pi2)
      }
      invisible(pi)
    },
    update_alpham = function(m, MAP = FALSE) {
      if (!MAP) {
        alpham <- self$emqr[m, , ] / self$nmqr[m, , ]
        alpham[is.nan(alpham)] <- 0
        self$alpham[[m]] <- alpham
      } else {
        alpham <- self$MAP$emqr[m, , ] / self$MAP$nmqr[m, , ]
        alpham[is.nan(alpham)] <- 0
        self$MAP$alpham[[m]] <- alpham
      }
      invisible(alpham)
    },
    update_alpha = function(MAP = FALSE) {
      if (!MAP) {
        alpha <- colSums(self$emqr, dims = 1) / colSums(self$nmqr, dims = 1)
        alpha[is.nan(alpha)] <- 0
        self$alpha <- alpha
      } else {
        alpha <- colSums(self$MAP$emqr, dims = 1) / colSums(self$MAP$nmqr, dims = 1)
        alpha[is.nan(alpha)] <- 0
        self$MAP$alpha <- alpha
      }
      invisible(alpha)
    },
    init_clust = function() {
      self$tau <-
        switch(self$init_method,
          "random" = lapply(
            X = seq_along(self$A),
            FUN = function(m) {
              tau_tmp <- gtools::rdirichlet(self$n[m], rep(1, self$Q))
              self$emqr[m, , ] <- .tquadform(tau_tmp, self$A[[m]] * self$mask[[m]])
              self$nmqr[m, , ] <- .tquadform(tau_tmp, self$mask[[m]])
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
              tau_tmp <- .one_hot(spectral_clustering(self$A[[m]], self$Q), self$Q)
              tau_tmp[tau_tmp < 1e-6] <- 1e-6
              tau_tmp[tau_tmp > 1 - 1e-6] <- 1 - 1e-6
              tau_tmp <- tau_tmp / .rowSums(tau_tmp, self$n[m], self$Q)
              self$emqr[m, , ] <- .tquadform(tau_tmp, self$A[[m]] * self$mask[[m]])
              self$nmqr[m, , ] <- .tquadform(tau_tmp, self$mask[[m]])
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
              if (is.matrix(self$Z[[m]])) {
                tau_tmp <- self$Z[[m]]
              } else {
                tau_tmp <- .one_hot(self$Z[[m]], self$Q)
                tau_tmp[tau_tmp < 1e-6] <- 1e-6
                tau_tmp[tau_tmp > 1 - 1e-6] <- 1 - 1e-6
                tau_tmp <- tau_tmp / .rowSums(tau_tmp, self$n[m], self$Q)
                self$emqr[m, , ] <- .tquadform(tau_tmp, self$A[[m]] * self$mask[[m]])
                self$nmqr[m, , ] <- .tquadform(tau_tmp, self$mask[[m]])
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
          )
        )
      lapply(seq(self$M), function(m) self$update_alpham(m))
      if (self$probablityDistribution == "bernoulli" & self$free_density &
        !self$fit_opts$approx_pois) {
        self$fixed_point_alpha_delta()
      }
    },
    make_permutation = function() {

    },
    m_step = function(MAP = FALSE, max_iter = 100, tol = 1e-3, ...) {
      # browser()
      lapply(seq_along(self$pi), function(m) self$update_pi(m, MAP = MAP))
      if (self$free_density == FALSE) {
        self$update_alpha(MAP = MAP)
      } else {
        switch(self$probablityDistribution,
          "poisson" = self$fixed_point_alpha_delta(MAP = MAP),
          "bernoulli" =
            ifelse(self$fit_opts$approx_pois,
              self$fixed_point_alpha_delta(MAP = MAP),
              self$update_alpha_delta(MAP = MAP)
            )
        )
      }
    },
    ve_step = function(m, max_iter = 20, tol = 1e-3, ...) {
      # Place holder for gradient ascent or other optimization methods
    },
    update_mqr = function(m) {
      # Compute the "mqr" quantities
      # emqr : the realised edges observed
      # nmqr : all the possible edges
      tau_tmp <- self$tau[[m]]
      self$emqr[m, , ] <- .tquadform(tau_tmp, self$A[[m]] * self$mask[[m]])
      self$nmqr[m, , ] <- .tquadform(tau_tmp, self$mask[[m]])
    },
    optimize = function(max_step = 100, tol = 1e-3, ...) {
      if (self$Q == 1) {
        self$tau <- lapply(seq(self$M), function(m) matrix(1, self$n[m], 1))
        self$Z <- lapply(seq(self$M), function(m) rep(1, self$n[m]))
        self$pi <- lapply(seq(self$M), function(m) 1)
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
        self$MAP <- list(
          Z = self$tau,
          pi = self$pi,
          alpha = self$alpha,
          delta = self$delta,
          alpham = self$alpham,
          emqr = self$emqr,
          nmqr = self$nmqr
        )
      } else {
        #  browser()
        self$init_clust()
        # lapply(seq(self$M), function(m) self$update_mqr(m))
        self$m_step(...)
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
        self$compute_MAP()
        lapply(seq(self$M), function(m) self$update_alpham(m))
        # self$compute_parameters()
        self$compute_icl()
        self$update_MAP_parameters()
      }
      self$compute_icl()
      self$compute_icl(MAP = TRUE)
      self$compute_icl_clustering()
    }
  ),
  active = list(
    dircoef = function() ifelse(self$directed, 1, .5)
  )
)
