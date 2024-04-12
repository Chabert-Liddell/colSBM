#' An R6 Class object, a fitted population of netowrks sbm
#' once $optimize() is done
#' 
#' @importFrom R6 R6Class
#'
#' @export
fitBipartiteSBMPop <- R6::R6Class(
  classname = "fitBipartiteSBMPop",
  # inherit = "bmpop",
  #
  public = list(
    # For user function add helpful `n[[m]]$row and n[[m]][[]]`
    #' @field n A list with two dimensions, each of size M for the rows and cols
    n = NULL,
    #' @field M Number of networks
    M = NULL,
    #' @field A List of incidence Matrix of size `n[[1]][m]xn[[2]][m]`
    A = NULL, 
    #' @field mask List of M masks, indicating NAs in the matrices. 1 for NA, 0 else
    mask = NULL,
    #'@field nonNAs List of M masks, indicating non NAs in the matrices. 1 - mask, so 0 for NA, 1 for non NA
    nonNAs = NULL, 
    #'@field nb_inter A vector of length M the number of unique non NA entries
    nb_inter = NULL,
    #'@field Q Number of clusters, vectors of size2
    Q = NULL, # 
    #'@field tau List of size M of list of two variational parameters.
    #' `n[[1]][m]xQ` matrices and `n[[2]][m]xQ` matrices
    tau = NULL,
    #'@field alpha Matrix of size QxQ, connection parameters
    alpha = NULL, # 
    #'@field delta  Vector of M,  density parameters with `delta[1] = 1`
    delta = NULL, #
    #'@field pi List of M vectors of size Q, the mixture parameters
    pi = NULL, #
    #'@field pim List of M vectors of size Q, the mixture parameters in case
    #' of free_mixture
    pim = NULL, 
    #'@field e Vector of size M, the sum of unique entries
    e = NULL, # 
    #'@field emqr List of M QxQ matrix, the sum of edges between q and r in m,
    #' ie the edges that are observed.
    emqr = NULL,
    #'@field nmqr list of M QxQ matrix, the number of entries between q and r 
    #' in m, ie all the possible edges.
    nmqr = NULL,
    #'@field alpham list of M QxQ matrix, the classic sbm parameters.
    alpham = NULL,
    #'@field free_mixture_row A boolean indicating if there is a free mixture
    #' on the rows
    free_mixture_row = NULL,
    #'@field free_mixture_col A boolean indicating if there is a free mixture
    #' on the columns
    free_mixture_col = NULL,
    #'@field free_density A boolean TODO delete
    free_density = NULL,
    #'@field weight A vector of size M for weighted likelihood
    weight = NULL,
    #'@field distribution Emission distribution either : "poisson" or 
    #' "bernoulli"
    distribution = NULL,
    #'@field mloss Loss on the M step of the VEM
    mloss = Inf,
    #'@field vloss Loss on the VE step of the VEM
    vloss = NULL,
    #'@field vbound The variational bound
    vbound = NULL,
    #'@field net_id A vector containing the "ids" or names of the networks 
    #' (if none given, they are set to their number in A list)
    net_id = NULL,
    #'@field df_mixture The degrees of freedom for mixture parameters pi,used 
    #' to compute penalty
    df_mixture = NULL,
    #'@field df_connect The degrees of freedom for connection parameters 
    #' alpha,used to compute penalty
    df_connect = NULL,
    #'@field df_density The degrees of freedom for density parameters delta,
    #' used to compute penalty
    df_density = NULL,
    #'@field Cpi  A list of matrices of size Qd x M containing TRUE (1)
    #' or FALSE (0) if the d-th dimension cluster is represented
    #' in the network m
    Cpi = NULL,
    #'@field Calpha The corresponding support on the connectivity parameters
    #' computed with Cpi.
    Calpha = NULL,
    #'@field logfactA A quantity used with the Poisson probability distribution
    logfactA = NULL,
    #'@field init_method The initialization method used for the first clustering
    init_method = NULL, # 
    #  verbosity = NULL,
    #'@field penalty  The penalty computed based on the number of parameters
    penalty = NULL,
    # approx_pois = NULL,
    #'@field Z  The clusters memberships, a list of size M of two matrices : 1 
    #' for rows clusters memberships and 2 for columns clusters memberships
    Z = NULL,
    #'@field MAP Maximum a posteriori
    MAP = NULL,
    #'@field MAP_parameters MAP params
    MAP_parameters = NULL,
    #'@field ICL Stores the ICL of the model
    ICL = NULL,
    #'@field BICL Stores the BICL of the model
    BICL = NULL,
    #'@field fit_opts Fit parameters, used to determine the fitting method/
    fit_opts = NULL,
    #'@field step_counter Counts the number of passes
    step_counter = 0,
    #'@field greedy_exploration_starting_point Stores the coordinates Q1 & Q2 
    #' from the greedy exploration to  keep track of the starting_point
    greedy_exploration_starting_point = NULL, 
    #' @field effective_clustering_list A list of size M storing the number 
    #' of the clusters that contains at least one point. Used for safety checks.
    effective_clustering_list = NULL, 
    #'@field clustering_is_complete A boolean used to know if the model real
    #' blocks match the expected blocks.
    clustering_is_complete = TRUE,

    #' @description
    #' Initializes the fitBipartiteSBMPop object
    #' 
    #' @param A List of incidence Matrix of size `n[[2]][m]xn[[2]][m]`
    #' @param Q A vector of size 2 with the number of row blocks and column 
    #' blocks
    #' @param Z  The clusters memberships, a list of size M of two matrices : 1 
    #' for rows clusters memberships and 2 for columns clusters memberships
    #' @param mask List of M masks, indicating NAs in the matrices. 1 for NA, 0 else
    #' @param net_id A vector containing the "ids" or names of the networks 
    #' (if none given, they are set to their number in A list)
    #' @param distribution Emission distribution either : "poisson" or 
    #' "bernoulli"
    #' @param free_mixture_row A boolean indicating if there is a free mixture
    #' on the rows
    #' @param free_mixture_col A boolean indicating if there is a free mixture
    #' on the columns
    #' @param free_density A boolean not used TODO delete
    #' @param Cpi  A list of matrices of size Qd x M containing TRUE (1)
    #' or FALSE (0) if the d-th dimension cluster is represented
    #' in the network m
    #' @param Calpha The corresponding support on the connectivity parameters
    #' computed with Cpi.
    #' @param init_method The initialization method used for the first clustering
    #' @param weight A vector of size M for weighted likelihood
    #' @param greedy_exploration_starting_point Stores the coordinates Q1 & Q2 
    #' from the greedy exploration to  keep track of the starting_point
    #' @param fit_opts Fit parameters, used to determine the fitting method/
    initialize = function(A = NULL,
                          Q = NULL,
                          Z = NULL,
                          mask = NULL,
                          net_id = NULL,
                          distribution = "bernoulli",
                          free_mixture_row = TRUE,
                          free_mixture_col = TRUE,
                          free_density = TRUE,
                          Cpi = NULL,
                          Calpha = NULL,
                          init_method = "spectral",
                          weight = NULL, # A vector of size M, the weight of each network
                          greedy_exploration_starting_point = NULL,
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


      # TODO : implement sanity checks

      self$greedy_exploration_starting_point <- greedy_exploration_starting_point

      self$A <- A
      self$M <- length(A)
      self$n <- vector("list", 2)
      self$n[[1]] <- vapply(seq_along(A),
        function(m) nrow(A[[m]]),
        FUN.VALUE = .1
      )
      self$n[[2]] <- vapply(seq_along(A),
        function(m) ncol(A[[m]]),
        FUN.VALUE = .1
      )
      self$e <- vapply(seq_along(A), function(m) sum(A[[m]], na.rm = TRUE),
        FUN.VALUE = .1
      )

      # Getting or computing the NA mask
      if (!is.null(mask)) {
        self$mask <- mask
      } else {
        self$mask <- lapply(
          seq_along(self$A),
          function(m) {
            mask <- matrix(0, nrow(self$A[[m]]), ncol(self$A[[m]]))
            if (sum(is.na(self$A[[m]] > 0))) {
              mask[is.na(self$A[[m]])] <- 1
            }
            mask
          }
        )
      }

      # Computing the non NA mask from the NA mask
      self$nonNAs <- lapply(
        seq_along(self$A),
        function(m) {
          nonNAs <- 1 - self$mask[[m]]
          nonNAs
        }
      )

      # Setting the net_id, if none is given, it's the position in the list that's used
      if (is.null(net_id)) {
        self$net_id <- seq(self$M)
      } else {
        self$net_id <- net_id
      }

      self$Q <- Q
      self$distribution <- distribution
      if (self$distribution == "poisson") {
        self$logfactA <- vapply(
          seq_along(self$A),
          function(m) {
            sum(lfactorial(self$A[[m]]), na.rm = TRUE)
          },
          FUN.VALUE = .1
        )
      }

      # Replacing NAs by -1 in the incidence matrices
      # To have no problems with the spectral clustering
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
      self$free_mixture_row <- free_mixture_row
      self$free_mixture_col <- free_mixture_col
      self$free_density <- free_density

      # Setting the needed matrices for free mixture models
      if (is.null(Cpi[[1]]) | is.null(Cpi[[2]]) |
        is.null(Calpha) | !self$free_mixture_row | !self$free_mixture_col) {
        if (is.null(Cpi[[1]])) {
          self$Cpi[[1]] <- matrix(TRUE, self$Q[1], self$M)
        }
        if (is.null(Cpi[[2]])) {
          self$Cpi[[2]] <- matrix(TRUE, self$Q[2], self$M)
        }
        self$Calpha <- tcrossprod(self$Cpi[[1]], self$Cpi[[2]]) > 0
      } else {
        # Sanity checks
        if (nrow(Cpi[[1]]) != self$Q[1] | ncol(Cpi[[1]]) != self$M) {
          stop(
            paste(
              "Wrong dimensions for Cpi[[1]] ! Should be nrow = ",
              self$Q[1],
              ", ncol = ", self$M
            )
          )
        }
        if (nrow(Cpi[[2]]) != self$Q[2] | ncol(Cpi[[2]]) != self$M) {
          stop(
            paste(
              "Wrong dimensions for Cpi[[2]] ! Should be nrow = ",
              self$Q[2],
              ", ncol = ", self$M
            )
          )
        }
        self$Cpi[[1]] <- matrix(Cpi[[1]],
          nrow = nrow(Cpi[[1]]), ncol = ncol(Cpi[[1]])
        )
        self$Cpi[[2]] <- matrix(Cpi[[2]],
          nrow = nrow(Cpi[[2]]), ncol = ncol(Cpi[[2]])
        )

        if (all(dim(Calpha) != self$Q)) {
          warning(
            "The provided Calpha has wrong dimensions ! Will recompute from Cpi"
          )
          Calpha <- tcrossprod(Cpi[[1]], Cpi[[2]])
        }

        self$Calpha <- Calpha
      }

      self$pi <- vector("list", self$M)
      self$pim <- vector("list", self$M)
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
      self$MAP$pim <- self$pim

      # Degrees of freedom
      # for iid
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
        self$delta <- (self$e / (self$n[[1]] * self$n[[2]])) /
          (self$e[1] / ((self$n[[1]][1] * self$n[[2]][1])))
      }

      self$alpha <- matrix(.5, Q[1], Q[2])
      self$init_method <- init_method
      self$nb_inter <- vapply(seq(self$M), function(m) sum(self$nonNAs[[m]]), FUN.VALUE = .1)
      self$vbound <- vector("list", self$M)
    },
    #' Method to compute the maximum a posteriori for Z clustering
    #' @return nothing; stores the values
    compute_MAP = function() {
      self$Z <- lapply(
        self$tau,
        function(tau) {
          list(
            # Retrieves the block membership for row
            row = apply(tau[[1]], 1, which.max),
            # Retrieves the block membership for col
            col = apply(tau[[2]], 1, which.max)
          )
        }
      )
      invisible(self$Z)
    },
    #' Computes the portion of the vbound with tau and alpha
    #' @param m The id of the network for which to compute
    #' @param MAP Wether to use the MAP parameters or not, a boolean, defaults
    #' to FALSE.
    #' @return The computed quantity.
    vb_tau_alpha = function(m, MAP = FALSE) {
      # Loading all the quantities useful for the computation

      # matrix_mqr_to_use <- outer(self$Cpi[[1]][, m], self$Cpi[[2]][, m])
      # OPTIM
      matrix_mqr_to_use <- self$Cpi[[1]][, m] %*% t(self$Cpi[[2]][, m])
      if (MAP) {
        # If we are calculating the maximum a posteriori
        emqr <- matrix_mqr_to_use * self$MAP$emqr[m, , ]
        nmqr <- matrix_mqr_to_use * self$MAP$nmqr[m, , ]
        alpha <- self$Calpha * self$MAP$alpha
        delta <- self$MAP$delta[m]
      } else {
        emqr <- matrix_mqr_to_use * self$emqr[m, , ]
        nmqr <- matrix_mqr_to_use * self$nmqr[m, , ]
        alpha <- self$Calpha * self$alpha
        delta <- self$delta[m]
      }
      switch(self$distribution,
        "bernoulli" = {
          sum(
            self$Calpha * .xlogy(emqr, alpha * delta, eps = 1e-12) +
              self$Calpha * .xlogy(nmqr - emqr, 1 - alpha * delta, eps = 1e-12)
          )
        },
        "poisson" = {
          sum(.xlogy(emqr, alpha * delta, eps = 1e-12) -
            nmqr * alpha * delta -
            self$logfactA[m])
        }
      )
    },
    #' Computes the portion of the vbound with tau pi and rho
    #' @param m The id of the network for which to compute
    #' @param MAP Wether to use the MAP parameters or not, a boolean, defaults
    #' to FALSE.
    #' @return The computed quantity.
    vb_tau_pi = function(m, MAP = FALSE) {
      if (!MAP) {
        sum(self$tau[[m]][[1]][, which(self$Cpi[[1]][, m])] %*%
          matrix(log(self$pi[[m]][[1]][which(self$Cpi[[1]][, m])]))) +
          sum(self$tau[[m]][[2]][, which(self$Cpi[[2]][, m])] %*%
            matrix(log(self$pi[[m]][[2]][which(self$Cpi[[2]][, m])])))
      } else {
        sum(.xlogy(
          tabulate(self$Z[[m]][[1]], self$Q[1]),
          self$MAP$pi[[m]][[1]]
        )) +
          sum(.xlogy(
            tabulate(self$Z[[m]][[2]], self$Q[2]),
            self$MAP$pi[[m]][[2]]
          ))
      }
    },
    #' Computes the entropy of the model
    #' @param m The id of the network for which to compute
    #' @return The computed quantity.
    entropy_tau = function(m) {
      -sum(.xlogx(self$tau[[m]][[1]][, which(self$Cpi[[1]][, m])])) -
        sum(.xlogx(self$tau[[m]][[2]][, which(self$Cpi[[2]][, m])]))
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
    #' Jacobian of the constraint
    #'
    #' @param par The parameters, alpha and delta combined in one big vector.
    #' @param emqr List of M QxQ matrix, the sum of edges between q and r in m
    #' @param nmqr List of M QxQ matrix, the number of entries between q and r 
    #' in m
    #'
    #' @return The evaluation of the function
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
    #' Updates the alpha and delta parameters
    #' 
    #' @param MAP Wether to use the MAP parameters or not, a boolean, defaults
    #' to FALSE.
    #' 
    #' @return nothing; but stores the values
    update_alpha_delta = function(MAP = FALSE) {
      # browser()
      # TODO Modify and adapt to bipartite case
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
      Cpi <- list()
      if (self$free_mixture_row) {
        Cpi[[1]] <- self$Cpi[[1]]
        pi1_penalty <- sum((colSums(Cpi[[1]]) - 1) * log(self$n[[1]]))
        log_p_Q1 <- -self$M * log(self$Q[1]) -
          sum(log(choose(
            rep(self$Q[1], self$M), colSums(Cpi[[1]])
          )))
        S1_penalty <- -2 * log_p_Q1
      } else {
        # To account for the possibility of the other free_mixture we store a
        # temporary support full of TRUE
        Cpi[[1]] <- matrix(TRUE, self$Q[1], self$M) # Cpi must be Q x M !
        # If there is no free mixture on the cols
        pi1_penalty <- (self$Q[1] - 1) * log(sum(self$n[[1]]))
        S1_penalty <- 0
      }

      if (self$free_mixture_col) {
        Cpi[[2]] <- self$Cpi[[2]]
        pi2_penalty <- sum((colSums(Cpi[[2]]) - 1) * log(self$n[[2]]))
        log_p_Q2 <- -self$M * log(self$Q[2]) - sum(log(choose(
          rep(self$Q[2], self$M), colSums(Cpi[[2]])
        )))
        S2_penalty <- -2 * log_p_Q2
      } else {
        # If there is no free mixture on the cols
        Cpi[[2]] <- matrix(TRUE, self$Q[2], self$M)
        pi2_penalty <- (self$Q[2] - 1) * log(sum(self$n[[2]]))
        S2_penalty <- 0
      }

      N_M <- sum(self$n[[1]] * self$n[[2]])

      if (self$free_mixture_row || self$free_mixture_col) {
        # If there is at least one free_mixture, the alpha penalty
        # will be computed using supports
        alpha_penalty <- sum(self$Calpha) * log(N_M)
      } else {
        # iid or delta
        alpha_penalty <- self$Q[1] * self$Q[2] * log(N_M)
      }
      self$penalty <- 0.5 * (pi1_penalty + pi2_penalty +
        alpha_penalty +
        S1_penalty + S2_penalty)
      if (self$penalty == -Inf | self$penalty == Inf) {
        stop("Infinite penalty !")
      }
      return(self$penalty)
    },
    #' Computes the ICL criterion
    #' @param MAP Wether to use the MAP parameters or not, a boolean, defaults
    #' to FALSE.
    #' @return The ICL for the model.
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
    #' Computes the BICL criterion
    #' @param MAP Wether to use the MAP parameters or not, a boolean, defaults
    #' to FALSE.
    #' @return The BICL for the model.
    compute_BICL = function(MAP = TRUE) {
      self$BICL <- self$compute_vbound() -
        self$compute_penalty()
      invisible(self$BICL)
    },

    #' Updates the MAP parameters
    #'
    #' @return nothing; but stores the values
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
            crossprod(ZR[[m]], (self$A[[m]] * (self$nonNAs[[m]])) %*% ZC[[m]])
        }
      )
      lapply(
        seq_along(Z),
        function(m) {
          self$MAP$nmqr[m, , ] <-
            crossprod(ZR[[m]], (self$nonNAs[[m]]) %*% ZC[[m]])
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

    #' Method to update tau values
    #' 
    #' @description
    #' Not really a fixed point as tau^1 depends only tau^2.
    #' 
    #' @param m The number of the network in the netlist
    #' @param d The dimension to update
    #' @param max_iter The maximum number of iterations to perform, defaults
    #' to 1
    #' @param tol The tolerance for which to stop iterating
    #' 
    #' @return The new tau values
    fixed_point_tau = function(m, d, max_iter = 1, tol = 1e-3) {
      # Just 1 step is necessary because tau1 depends only on tau2
      condition <- TRUE
      it <- 0
      # reup_counter <- 0
      self$vloss[[m]] <-
        c(self$vloss[[m]], self$vb_tau_alpha(m) + self$vb_tau_pi(m) +
          self$entropy_tau(m))
      tau_new <- switch(self$distribution,
        "bernoulli" = {
          tau_new <-
            if (d == 1) {
              # n[[1]] * Q1
              tau_new <-
                t(matrix(
                  .xlogy(self$Cpi[[1]][, m],
                    self$pi[[m]][[d]],
                    eps = NULL
                  ),
                  self$Q[d], self$n[[1]][m]
                )) +
                ((self$nonNAs[[m]]) * self$A[[m]]) %*%
                t(self$Cpi[[2]][, m] * t(self$tau[[m]][[2]])) %*%
                t(.logit(self$Calpha * self$delta[m] * self$alpha,
                  eps = 1e-9
                )) +
                (self$nonNAs[[m]]) %*%
                t(self$Cpi[[2]][, m] * t(self$tau[[m]][[2]])) %*%
                t(.log(1 - self$Calpha * self$alpha * self$delta[m],
                  eps = 1e-9
                ))
              # In order to fix NaN appearing in the formula (log(Pi) when Pi
              # = 0), the .xlogy function is used with eps = 1e-9
              # POSSIBLE POINT OF FAILURE
            }
          if (d == 2) {
            # n[[2]] * Q2
            tau_new <-
              t(matrix(
                .xlogy(self$Cpi[[2]][, m],
                  self$pi[[m]][[d]],
                  eps = NULL
                ),
                self$Q[d], self$n[[2]][m]
              )) +
              t((self$nonNAs[[m]]) * self$A[[m]]) %*%
              t(self$Cpi[[1]][, m] * t(self$tau[[m]][[1]])) %*%
              .logit(self$Calpha * self$delta[m] * self$alpha, eps = 1e-9) +
              t(self$nonNAs[[m]]) %*%
              t(self$Cpi[[1]][, m] * t(self$tau[[m]][[1]])) %*%
              .log(1 - self$Calpha * self$alpha * self$delta[m], eps = 1e-9)
            # In order to fix NaN appearing in the formula (log(Pi) when Pi
            # = 0), the .xlogy function is used with eps = 1e-9
            # POSSIBLE POINT OF FAILURE
          }
          invisible(tau_new)
        },
        "poisson" = {
          if (d == 1) {
            tau_new <-
              t(matrix(log(self$pi[[m]][[d]]), self$Q[d], self$n[[1]][m])) +
              ((self$nonNAs[[m]]) * self$A[[m]]) %*%
              self$tau[[m]][[2]] %*%
              t(log(self$delta[m] * self$alpha)) -
              (self$nonNAs[[m]]) %*%
              self$tau[[m]][[2]] %*%
              t(self$alpha * self$delta[m])
          }
          if (d == 2) {
            tau_new <-
              t(matrix(log(self$pi[[m]][[d]]), self$Q[d], self$n[[2]][m])) +
              t((self$nonNAs[[m]]) * self$A[[m]]) %*%
              self$tau[[m]][[1]] %*%
              log(self$delta[m] * self$alpha) -
              t(self$nonNAs[[m]]) %*%
              self$tau[[m]][[1]] %*%
              (self$alpha * self$delta[m])
          }
          invisible(tau_new)
        }
      )
      # To force the invalid taus to be zero
      # which(self$Cpi[[d]][,m]) returns a vector with the index of the support
      # for the network m that are TRUE
      # Select the tau for the individuals
      tau_new[, which(!self$Cpi[[d]][, m])] <- 0

      # tau_new[, which(self$Cpi[[d]][, m])] <- .softmax(
      #   matrix(tau_new[, which(self$Cpi[[d]][, m])], nrow = ifelse(d == 1,
      #     self$n[[1]][m], # This can be improved by using self$n[[d]][m] # TODO implement
      #     self$n[[2]][m]
      #   ), ncol = sum(self$Cpi[[d]][,m]))
      # )
      # tau_new[, which(self$Cpi[[d]][, m])][tau_new[, which(self$Cpi[[d]][, m])] < 1e-9] <- 1e-9
      # tau_new[, which(self$Cpi[[d]][, m])][tau_new[, which(self$Cpi[[d]][, m])] > 1 - 1e-9] <- 1 - 1e-9
      # OPTIM
      tau_cols_to_use <- which(self$Cpi[[d]][, m])
      tau_new[, tau_cols_to_use] <- .softmax(
        matrix(tau_new[, tau_cols_to_use],
          nrow = self$n[[d]][m],
          ncol = sum(self$Cpi[[d]][, m])
        )
      )
      tau_new[, tau_cols_to_use][tau_new[, tau_cols_to_use] < 1e-9] <- 1e-9
      tau_new[, tau_cols_to_use][tau_new[, tau_cols_to_use] > 1 - 1e-9] <- 1 - 1e-9

      tau_new <- tau_new / rowSums(tau_new)
      self$tau[[m]][[d]] <- tau_new
      invisible(tau_new)
    },

    #' Fixed point to update alpha and delta 
    #' TODO Check and fix
    #' 
    #' @param MAP A boolean wether to use MAP parameters or not, defaults to
    #' FALSE
    #' @param max_iter The maximum number of iterations, default to 50
    #' @param tol The tolerance for which to stop iterating
    #' 
    #' @return nothing; stores the values
    fixed_point_alpha_delta = function(MAP = FALSE, max_iter = 50, tol = 1e-6) {
      # switch(
      #   self$distribution,
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
    #' Computes the pi per network, known as the pim
    #' 
    #' @param m The number of the network in the netlist
    #' @param MAP A boolean wether to use the MAP parameters or not, defaults to
    #' FALSE
    #' 
    #' @return nothing; stores the values
    update_pim = function(m, MAP = FALSE) {
      # Here we compute the pi for each of the m networks
      # So those are the "pim"
      if (!MAP) {
        # The pim are the mean of all the tau for each column, ie the mean tau per block
        pim1 <- self$Cpi[[1]][, m] * matrix(.colMeans(
          self$tau[[m]][[1]],
          self$n[[1]][m], self$Q[1]
        ), nrow = 1)
        pim2 <- self$Cpi[[2]][, m] * matrix(.colMeans(
          self$tau[[m]][[2]],
          self$n[[2]][m], self$Q[2]
        ), nrow = 1)
        self$pim[[m]] <- list(pim1, pim2)
      } else {
        # The MAP Pi are obtained by counting the number of nodes in each block and dividing by the total (n[[1]] or n[[2]])
        pim1 <- tabulate(self$Z[[m]][[1]], self$Q[1]) / self$n[[1]][m]
        pim2 <- tabulate(self$Z[[m]][[2]], self$Q[2]) / self$n[[2]][m]
        self$MAP$pim[[m]] <- list(pim1, pim2)
      }
    },
    #' Computes the pi for the whole model
    #'
    #' @param MAP A boolean wether to use the MAP parameters or not, defaults to
    #' FALSE
    #'
    #' @return the pi and stores the values
    update_pi = function(MAP = FALSE) {
      # The function is m subscripted but it might not need it to work
      for (m in seq.int(self$M)) {
        self$update_pim(m, MAP)
      }
      if (!MAP) {
        if (self$free_mixture_row | self$free_mixture_col) {
          # TODO : Can remove  the outer OR if and else to use the below
          # If free_mixture the pi are the pim
          if (self$free_mixture_row) {
            for (m in seq.int(self$M)) {
              self$pi[[m]][[1]] <- self$pim[[m]][[1]]
            }
          } else {
            # No free_mixture on row so averaging on pi[[m]][[1]]
            pi1 <- matrix(self$n[[1]] * vapply(seq(self$M), function(m) {
              self$pim[[m]][[1]]
            }, FUN.VALUE = rep(.1, self$Q[1])), ncol = self$M, nrow = self$Q[1])
            pi1 <- matrix(rowSums(pi1) / sum(pi1), nrow = 1)
            for (m in seq.int(self$M)) {
              self$pi[[m]][[1]] <- pi1
            }
          }

          if (self$free_mixture_col) {
            for (m in seq.int(self$M)) {
              self$pi[[m]][[2]] <- self$pim[[m]][[2]]
            }
          } else {
            # No free_mixture on col so averaging on pi[[m]][[2]]
            pi2 <- matrix(self$n[[2]] * vapply(seq(self$M), function(m) {
              self$pim[[m]][[2]]
            }, FUN.VALUE = rep(.1, self$Q[2])), ncol = self$M, nrow = self$Q[2])
            pi2 <- matrix(rowSums(pi2) / sum(pi2), nrow = 1)
            for (m in seq.int(self$M)) {
              self$pi[[m]][[2]] <- pi2
            }
          }
        } else {
          # Otherwise we need to ponder based on the size for
          # Rows
          pi1 <- matrix(self$n[[1]] * vapply(seq(self$M), function(m) {
            self$pim[[m]][[1]]
          }, FUN.VALUE = rep(.1, self$Q[1])), ncol = self$M, nrow = self$Q[1])
          pi1 <- matrix(rowSums(pi1) / sum(pi1), nrow = 1)
          # Need matrix() because there is a dot product in the vb_tau_pi

          # And columns
          pi2 <- matrix(self$n[[2]] * vapply(seq(self$M), function(m) {
            self$pim[[m]][[2]]
          }, FUN.VALUE = rep(.1, self$Q[2])), ncol = self$M, nrow = self$Q[2])
          pi2 <- matrix(rowSums(pi2) / sum(pi2), nrow = 1)
          # Need matrix() because there is a dot product in the vb_tau_pi

          self$pi <- lapply(seq.int(self$M), function(m) list(pi1, pi2))
        }
      } else {
        if (self$free_mixture_row | self$free_mixture_col) {
          # If free_mixture the pi are the pim
          self$MAP$pi <- self$MAP$pim
        } else {
          # Otherwise we need to ponder based on the size for
          # Rows
          pi1 <- matrix(self$n[[1]] * vapply(seq(self$M), function(m) {
            self$MAP$pim[[m]][[1]]
          }, FUN.VALUE = rep(.1, self$Q[1])), ncol = self$M, nrow = self$Q[1])

          pi1 <- rowSums(pi1) / sum(pi1)

          # And columns
          pi2 <- matrix(self$n[[2]] * vapply(seq(self$M), function(m) {
            self$MAP$pim[[m]][[2]]
          }, FUN.VALUE = rep(.1, self$Q[2])), ncol = self$M, nrow = self$Q[2])

          pi2 <- rowSums(pi2) / sum(pi2)

          self$MAP$pi <- lapply(seq.int(self$M), function(m) list(pi1, pi2))
        }
      }
      invisible(pi)
    },
    #' Computes the alpha per network, known as the alpĥam
    #' 
    #' @param m The number of the network in the netlist
    #' @param MAP A boolean wether to use the MAP parameters or not, defaults to
    #' FALSE
    #' 
    #' @return the alpham and stores the values
    update_alpham = function(m, MAP = FALSE) {
      if (!MAP) {
        alpham <- self$emqr[m, , ] / self$nmqr[m, , ]
        alpham[is.nan(alpham)] <- 0
        alpham <- matrix(alpham, self$Q[1], self$Q[2])
        self$alpham[[m]] <- alpham
      } else {
        alpham <- self$MAP$emqr[m, , ] / self$MAP$nmqr[m, , ]
        alpham[is.nan(alpham)] <- 0
        alpham <- matrix(alpham, self$Q[1], self$Q[2])
        self$MAP$alpham[[m]] <- alpham
      }
      invisible(alpham)
    },
    #' Computes the alpha for the whole model
    #'
    #' @param MAP A boolean wether to use the MAP parameters or not, defaults to
    #' FALSE
    #'
    #' @return the alpha and stores the values
    update_alpha = function(MAP = FALSE) {
      # For iid and pi-colSBM
      if (!MAP) {
        alpha <- self$Calpha * colSums(self$emqr, dims = 1) / colSums(self$nmqr, dims = 1)
        alpha[is.nan(alpha)] <- 0
        self$alpha <- alpha
      } else {
        alpha <- self$Calpha * colSums(self$MAP$emqr, dims = 1) / colSums(self$MAP$nmqr, dims = 1)
        alpha[is.nan(alpha)] <- 0
        self$MAP$alpha <- alpha
      }
      invisible(alpha)
    },

    #' @description Initialize clusters
    #'
    #' @importFrom gtools rdirichlet
    #' @return nothing; stores
    init_clust = function() {
      self$tau <-
        switch(self$init_method,
          # TODO later : adapt init_clust "random" to bipartite
          "random" = lapply(
            X = seq_along(self$A),
            FUN = function(m) {
              tau_1 <- gtools::rdirichlet(self$n[[1]][m], rep(1, self$Q[[1]]))
              tau_2 <- gtools::rdirichlet(self$n[[2]][m], rep(1, self$Q[[2]]))
              self$emqr[m, , ] <- .tquadform(tau_tmp, self$A[[m]] * (self$nonNAs[[m]]))
              self$nmqr[m, , ] <- .tquadform(tau_tmp, (self$nonNAs[[m]]))
              a <- self$emqr[m, , ] / self$nmqr[m, , ]
              prob <- self$Q * diag(a) #+ rowSums(a)
              p <- sample.int(self$Q, prob = prob)
              tau_tmp <- tau_tmp[, p]
              self$emqr[m, , ] <- self$emqr[m, p, p]
              self$nmqr[m, , ] <- self$nmqr[m, p, p]
              tau_tmp
            }
          ),
          "hca" = lapply(
            X = seq_along(self$A),
            FUN = function(m) {
              # The .one_hot performs a one hot encoding of the spectral clustering performed
              # DONE : Adapt this step to handle two taus
              biclustering <- bipartite_hierarchic_clustering(self$nonNAs[[m]] * self$A[[m]], self$Q)
              row_clustering <- biclustering$row_clustering
              col_clustering <- biclustering$col_clustering

              # Tau for the rows
              tau_1 <- .one_hot(row_clustering, self$Q[[1]])
              tau_1[tau_1 < 1e-6] <- 1e-6
              tau_1[tau_1 > 1 - 1e-6] <- 1 - 1e-6
              tau_1 <- tau_1 / .rowSums(tau_1, self$n[[1]][m], self$Q[[1]])

              # Tau for the columns
              tau_2 <- .one_hot(col_clustering, self$Q[[2]])
              tau_2[tau_2 < 1e-6] <- 1e-6
              tau_2[tau_2 > 1 - 1e-6] <- 1 - 1e-6
              tau_2 <- tau_2 / .rowSums(tau_2, self$n[[2]][m], self$Q[[2]])

              self$emqr[m, , ] <- t(tau_1) %*% (self$A[[m]] * (self$nonNAs[[m]])) %*% tau_2
              self$nmqr[m, , ] <- t(tau_1) %*% (self$nonNAs[[m]]) %*% tau_2

              # Permuter avec un ordre favorisant plus fortement
              # les plus gros clusters mais de manière stochastique
              # colSums sur les tau1 et tau2 pour obtenir les pi1 et pi2
              a <- self$emqr[m, , ] / self$nmqr[m, , ] # ? estimate of the alpha
              pi1 <- .colSums(tau_1, self$n[[1]][m], self$Q[[1]]) / sum(tau_1)
              pi2 <- .colSums(tau_2, self$n[[2]][m], self$Q[[2]]) / sum(tau_2)

              prob1 <- as.vector(pi2 %*% t(a))
              prob2 <- as.vector(pi1 %*% a)
              # p1 and p2 contain the clustering ordered using the highest probabilities first

              p1 <- order(prob1)
              p2 <- order(prob2)

              # Tau, emqr and nmqr are reordered accordingly
              if (ncol(tau_1) != 1) {
                tau_1 <- tau_1[, p1]
                self$emqr[m, , ] <- self$emqr[m, p1, ]
                self$nmqr[m, , ] <- self$nmqr[m, p1, ]
              }
              if (ncol(tau_2) != 1) {
                tau_2 <- tau_2[, p2]
                self$emqr[m, , ] <- self$emqr[m, , p2]
                self$nmqr[m, , ] <- self$nmqr[m, , p2]
              }

              # The output is tau[[m]][[1]] for tau_1 and tau[[m]][[2]] for tau_2
              list(tau_1, tau_2)
            }
          ),
          # DONE : adapt "spectral" clustering to bipartite with two clustering on the rows and columns
          "spectral" = lapply(
            X = seq_along(self$A),
            FUN = function(m) {
              # The .one_hot performs a one hot encoding of the spectral clustering performed

              biclustering <- spectral_biclustering(self$nonNAs[[m]] * self$A[[m]], self$Q)
              row_clustering <- biclustering$row_clustering
              col_clustering <- biclustering$col_clustering

              # Tau for the rows
              tau_1 <- .one_hot(row_clustering, self$Q[[1]])
              tau_1[tau_1 < 1e-6] <- 1e-6
              tau_1[tau_1 > 1 - 1e-6] <- 1 - 1e-6
              tau_1 <- tau_1 / .rowSums(tau_1, self$n[[1]][m], self$Q[[1]])

              # Tau for the columns
              tau_2 <- .one_hot(col_clustering, self$Q[[2]])
              tau_2[tau_2 < 1e-6] <- 1e-6
              tau_2[tau_2 > 1 - 1e-6] <- 1 - 1e-6
              tau_2 <- tau_2 / .rowSums(tau_2, self$n[[2]][m], self$Q[[2]])

              self$emqr[m, , ] <- t(tau_1) %*% (self$A[[m]] * (self$nonNAs[[m]])) %*% tau_2
              self$nmqr[m, , ] <- t(tau_1) %*% (self$nonNAs[[m]]) %*% tau_2

              # Permuter avec un ordre favorisant plus fortement
              # les plus gros clusters mais de manière stochastique
              # colSums sur les tau1 et tau2 pour obtenir les pi1 et pi2
              a <- self$emqr[m, , ] / self$nmqr[m, , ] # ? estimate of the alpha
              pi1 <- .colSums(tau_1, self$n[[1]][m], self$Q[[1]]) / sum(tau_1)
              pi2 <- .colSums(tau_2, self$n[[2]][m], self$Q[[2]]) / sum(tau_2)

              prob1 <- as.vector(pi2 %*% t(a))
              prob2 <- as.vector(pi1 %*% a)
              # p1 and p2 contain the clustering ordered using the highest probabilities first
              # But it is still sampled and random according to the probabilities
              # p1 <- sample.int(self$Q[[1]], prob = prob1)
              # p2 <- sample.int(self$Q[[2]], prob = prob2)

              if (ncol(tau_1) != 1) {
                p1 <- order(prob1)
              }

              if (ncol(tau_2)) {
                p2 <- order(prob2)
              }

              # Tau, emqr and nmqr are reordered accordingly
              if (ncol(tau_1) != 1) {
                tau_1 <- tau_1[, p1]
                self$emqr[m, , ] <- self$emqr[m, p1, ]
                self$nmqr[m, , ] <- self$nmqr[m, p1, ]
              }
              if (ncol(tau_2) != 1) {
                tau_2 <- tau_2[, p2]
                self$emqr[m, , ] <- self$emqr[m, , p2]
                self$nmqr[m, , ] <- self$nmqr[m, , p2]
              }

              list(tau_1, tau_2)
            }
          ),
          "empty" = lapply(
            seq_along(self$Z),
            function(m) {
              tau_1 <- matrix(self$tau[[m]][[1]],
                nrow = self$n[[1]][m],
                ncol = self$Q[1]
              )
              tau_1[, which(!self$Cpi[[1]][, m])] <- 0
              tau_1[, which(self$Cpi[[1]][, m])] <-
                .threshold(matrix(tau_1[, which(self$Cpi[[1]][, m])],
                  nrow = self$n[[1]][m]
                ))

              tau_2 <- matrix(self$tau[[m]][[2]],
                nrow = self$n[[2]][m],
                ncol = self$Q[2]
              )
              tau_2[, which(!self$Cpi[[2]][, m])] <- 0
              tau_2[, which(self$Cpi[[2]][, m])] <-
                .threshold(matrix(tau_2[, which(self$Cpi[[2]][, m])],
                  nrow = self$n[[2]][m]
                ))

              self$emqr[m, , ] <- t(tau_1) %*%
                (self$A[[m]] * (self$nonNAs[[m]])) %*% tau_2
              self$nmqr[m, , ] <- t(tau_1) %*% (self$nonNAs[[m]]) %*% tau_2
              list(tau_1, tau_2)
            }
          ),
          "given" = lapply(
            X = seq_along(self$Z),
            FUN = function(m) {
              if (is.matrix(self$Z[[m]][[1]]) &&
                is.matrix(self$Z[[m]][[2]]) &&
                all(dim(self$Z[[m]][[1]]) == c(self$n[[1]][m], self$Q[1])) &&
                all(dim(self$Z[[m]][[2]]) == c(self$n[[2]][m], self$Q[2]))) {
                # If Z was already provided as a list of two matrices
                tau_1 <- self$Z[[m]][[1]]
                tau_2 <- self$Z[[m]][[2]]
              } else {
                # Tau for the rows
                tau_1 <- .one_hot(self$Z[[m]][[1]], self$Q[[1]])
                tau_1[tau_1 < 1e-6] <- 1e-6
                tau_1[tau_1 > 1 - 1e-6] <- 1 - 1e-6
                tau_1 <- tau_1 / .rowSums(tau_1, self$n[[1]][m], self$Q[[1]])

                # Tau for the columns
                tau_2 <- .one_hot(self$Z[[m]][[2]], self$Q[[2]])
                tau_2[tau_2 < 1e-6] <- 1e-6
                tau_2[tau_2 > 1 - 1e-6] <- 1 - 1e-6
                tau_2 <- tau_2 / .rowSums(tau_2, self$n[[2]][m], self$Q[[2]])

                # update_mqr(m)
                self$emqr[m, , ] <- t(tau_1) %*% (self$A[[m]] * (self$nonNAs[[m]])) %*% tau_2
                self$nmqr[m, , ] <- t(tau_1) %*% (self$nonNAs[[m]]) %*% tau_2
              }
              return(list(tau_1, tau_2))
            }
          )
          # TODO add a "given_tau" init method
        )
      lapply(seq(self$M), function(m) self$update_alpham(m))
      if (self$distribution == "bernoulli" & self$free_density &
        !self$fit_opts$approx_pois) {
        self$fixed_point_alpha_delta()
      }
    },
    #' The M step of the VEM
    #'
    #' @param MAP A boolean wether to use the MAP parameters or not, defaults to
    #' FALSE
    #' @param max_iter The maximum number of iterations, default to 2
    #' @param tol The tolerance for which to stop iterating defaults to 1e-3
    #' @param ... Other parameters
    #' 
    #' @return nothing; stores values
    m_step = function(MAP = FALSE, max_iter = 2, tol = 1e-3, ...) {
      # browser()
      # lapply(seq_along(self$pi), function(m) self$update_pi(m, MAP = MAP))
      self$update_pi(MAP = MAP)
      if (self$free_density == FALSE) {
        # deltas are all equals to 1
        self$update_alpha(MAP = MAP)
      } else {
        switch(self$distribution,
          "poisson" = self$fixed_point_alpha_delta(MAP = MAP),
          "bernoulli" =
            ifelse(self$fit_opts$approx_pois,
              self$fixed_point_alpha_delta(MAP = MAP),
              self$update_alpha_delta(MAP = MAP)
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
    ve_step = function(m, max_iter = 2, tol = 1e-3, ...) {
      # Place holder for gradient ascent or other optimization methods
    },
    #' Updates the mqr quantities
    #' 
    #' @description
    #' Namely, it updates the emqr and nmqr.
    #' 
    #' @param m The number of the network in the netlist
    update_mqr = function(m) {
      # Compute the "mqr" quantities
      # emqr : the observed edges
      # nmqr : all the possible edges

      tau_m_1 <- self$tau[[m]][[1]]
      tau_m_2 <- self$tau[[m]][[2]]

      self$emqr[m, , ] <- t(tau_m_1) %*%
        (self$A[[m]] * (self$nonNAs[[m]])) %*% tau_m_2
      self$nmqr[m, , ] <- t(tau_m_1) %*% (self$nonNAs[[m]]) %*% tau_m_2
    },
    #' TODO Investigate what its supposed to do
    #' @return nothing
    make_permutation = function() {
      # TODO : ask @Chabert-Liddell or @demiperimetre to see if such permutation are useful at the end
      if (self$free_mixture_row | !self$free_mixture_col) {
        # Use pim
        pi1 <- lapply(seq.int(self$M), function(m) {
          self$pim[[m]][[1]]
        })
        pi2 <- lapply(seq.int(self$M), function(m) {
          self$pim[[m]][[2]]
        })
      } else {
        # Use pi
        pi1 <- lapply(seq.int(self$M), function(m) {
          self$pi[[m]][[1]]
        })
        pi2 <- lapply(seq.int(self$M), function(m) {
          self$pi[[m]][[2]]
        })
      }

      # Compute the marginal laws
      prob_row <- pi2 %*% t(self$alpha)
      prob_col <- pi1 %*% self$alpha

      # Find the order
      row_order <- order(prob_row, decreasing = TRUE)
      col_order <- order(col_row, decreasing = TRUE)

      # Tau, emqr and nmqr, alpha & pi are reordered accordingly
      # if (ncol(tau_1) != 1) {
      #   next
      # }
      # if (ncol(tau_2) != 1) {
      #   next
      # }
    },
    #' Computes the number of blocks that are effectively populated
    #' 
    #' @return nothing; but stores the value
    compute_effective_clustering = function() {
      self$effective_clustering_list <- lapply(seq.int(self$M), function(m) {
        list(
          length(unique(self$MAP$Z[[m]][[1]])),
          length(unique(self$MAP$Z[[m]][[2]]))
        )
      })
    },
    #' Perform the whole initialization and VEM algorithm
    #' 
    #' @param max_step The maximum number of steps to perform optimization
    #' @param tol The tolerance for which to stop iterating, default to 1e-3
    #' @param ... Other parameters
    #' 
    #' @return nothing
    optimize = function(max_step = 20, tol = 1e-3, ...) {
      if (all(self$Q == c(1, 1))) {
        self$tau <- lapply(
          seq(self$M),
          function(m) list(matrix(1, self$n[[1]][m], 1), matrix(1, self$n[[2]][m], 1))
        )
        self$Z <- lapply(
          seq(self$M),
          function(m) list(rep(1, self$n[[1]][m]), rep(1, self$n[[2]][m]))
        )
        self$pi <- lapply(
          seq(self$M),
          function(m) list(matrix(1), matrix(1))
        )
        lapply(
          seq(self$M),
          function(m) self$update_mqr(m)
        )
        if (self$free_density) {
          # TODO later : implement the free density
          self$alpha <- matrix(
            sum(self$emqr[1, , ]) / sum(self$nmqr[1, , ]), 1, 1
          )
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
          function(m) {
            matrix(
              sum(self$emqr[m, , ]) / sum(self$nmqr[m, , ]), 1, 1
            )
          }
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
        # If there is more than one cluster
        self$init_clust()
        # lapply(seq(self$M), function(m) self$update_mqr(m))
        self$m_step(...)
        step <- 0
        vb <- self$compute_vbound()
        self$vbound <- vb
        step_condition <- TRUE
        while (step_condition) {
          seq_m_minibatch <- sample.int(self$M)
          lapply(
            # Minibatch or not
            ifelse(self$fit_opts$minibatch, seq_m_minibatch, seq(self$M)),
            function(m) {
              switch(self$fit_opts$algo_ve,
                "fp" = {
                  self$fixed_point_tau(m, d = 1)
                  self$update_mqr(m)
                  self$m_step(...)

                  self$fixed_point_tau(m, d = 2)
                  self$update_mqr(m)
                  self$m_step(...)
                },
                # If we're not using the previous methods default to gradient ascent
                self$ve_step(m, ...)
              )
            }
          )

          vb <- self$compute_vbound()
          self$vbound <- c(self$vbound, vb)
          step <- step + 1
          step_condition <- step < max_step &
            (self$vbound[length(self$vbound)] -
              self$vbound[length(self$vbound) - 1] > tol)
          if (step %% 5 == 0) {
            if (self$fit_opts$verbosity >= 1) {
              print(paste0(step, ": ", vb))
              print(self$alpha)
            }
          }
        }

        # Here we stock the number of steps needed to converge
        self$step_counter <- step

        self$compute_MAP()
        lapply(seq(self$M), function(m) self$update_alpham(m))
        # self$compute_parameters()
        self$compute_icl()
        self$update_MAP_parameters()
      }
      self$compute_icl()
      self$compute_icl(MAP = TRUE)

      # After all the effective clustering is computed
      self$compute_effective_clustering()

      # Here the matching of the effective clustering is checked
      if (!any(self$free_mixture_row, self$free_mixture_col)) {
        non_correct_clusterings <- sapply(seq.int(self$M), function(m) {
          self$effective_clustering_list[[m]][[1]] != self$Q[1] &&
            self$effective_clustering_list[[m]][[2]] != self$Q[2]
        })
        if (any(non_correct_clusterings)) {
          self$clustering_is_complete <- FALSE
        }
      }

      self$compute_BICL()
      if (self$BICL == Inf | self$BICL == -Inf) {
        stop("Infinite BICL !")
      }
      self$reorder_parameters()
    },
    #' Reorder the blocks putting the "strongest" ones first in order to have
    #' a coherent ordering of blocks with SBM and LBM for visualisation.
    #'
    #' @return nothing; stores the new ordering
    reorder_parameters = function() {
      Z_label_switch <- function(Z, new_order) {
        # Create a mapping of old labels to new labels
        old_names <- names(Z)
        label_map <- setNames(new_order, unique(Z))

        # Use the mapping to replace labels in the vector
        switched_labels <- label_map[Z]
        names(switched_labels) <- old_names
        return(switched_labels)
      }
      if (all(self$Q == c(1, 1))) {
        return(self)
      }
      mean_pi <- sapply(self$pi, function(pi) pi[[1]])
      if (self$Q[1] > 1) {
        mean_pi <- matrixStats::rowMeans2(mean_pi)
      } else {
        mean_pi <- matrix(1, 1, 1)
      }
      mean_rho <- sapply(self$pi, function(pi) pi[[2]])
      if (self$Q[2] > 1) {
        mean_rho <- matrixStats::rowMeans2(mean_rho)
      } else {
        mean_rho <- matrix(1, 1, 1)
      }
      # The row clustering are reordered according to their marginal distribution
      prob1 <- as.vector(mean_rho %*% t(self$MAP$alpha))
      p1 <- order(prob1, decreasing = TRUE)

      # The col clustering are reordered according to their marginal distribution
      prob2 <- as.vector(mean_pi %*% self$MAP$alpha)
      p2 <- order(prob2, decreasing = TRUE)

      # m independent
      self$MAP$alpha <- self$MAP$alpha[p1, p2, drop = FALSE]
      self$Calpha <- self$Calpha[p1, p2, drop = FALSE]
      self$alpha <- self$alpha[p1, p2, drop = FALSE]

      # m dependent
      lapply(seq.int(self$M), function(m) {
        # Reordering the parameters
        self$Cpi[[1]][, m] <- self$Cpi[[1]][p1, m, drop = FALSE]
        self$Cpi[[2]][, m] <- self$Cpi[[2]][p2, m, drop = FALSE]

        self$pim[[m]][[1]] <- self$pim[[m]][[1]][p1, drop = FALSE]
        self$pim[[m]][[2]] <- self$pim[[m]][[2]][p2, drop = FALSE]

        self$pi[[m]][[1]] <- self$pi[[m]][[1]][p1, drop = FALSE]
        self$pi[[m]][[2]] <- self$pi[[m]][[2]][p2, drop = FALSE]

        self$emqr[m, , ] <- self$emqr[m, p1, p2, drop = FALSE]
        self$nmqr[m, , ] <- self$nmqr[m, p1, p2, drop = FALSE]
        self$alpham[[m]] <- matrix(
          self$alpham[[m]][p1, p2, drop = FALSE],
          self$Q[1], self$Q[2]
        )
        self$tau[[m]][[1]] <- self$tau[[m]][[1]][, p1, drop = FALSE]
        self$tau[[m]][[2]] <- self$tau[[m]][[2]][, p2, drop = FALSE]
        self$Z[[m]][[1]] <- Z_label_switch(self$Z[[m]][[1]], p1)
        self$Z[[m]][[2]] <- Z_label_switch(self$Z[[m]][[2]], p2)

        # MAP parameters
        # Work needed to relabel correctly!
        self$MAP$Z[[m]][[1]] <- Z_label_switch(self$MAP$Z[[m]][[1]], p1)
        self$MAP$Z[[m]][[2]] <- Z_label_switch(self$MAP$Z[[m]][[2]], p2)
        self$MAP$emqr[m, , ] <- self$MAP$emqr[m, p1, p2, drop = FALSE]
        self$MAP$nmqr[m, , ] <- self$MAP$nmqr[m, p1, p2, drop = FALSE]
        self$MAP$alpham[[m]] <- matrix(
          self$MAP$alpham[[m]][p1, p2, drop = FALSE],
          self$Q[1], self$Q[2]
        )
        self$MAP$pim[[m]][[1]] <- self$MAP$pim[[m]][[1]][p1, drop = FALSE]
        self$MAP$pim[[m]][[2]] <- self$MAP$pim[[m]][[2]][p2, drop = FALSE]
        self$MAP$pi[[m]][[1]] <- self$MAP$pi[[m]][[1]][p1, drop = FALSE]
        self$MAP$pi[[m]][[2]] <- self$MAP$pi[[m]][[2]][p2, drop = FALSE]
      })
    },
    #' The message printed when one prints the object
    #' @param type The title above the message.
    show = function(type = "Fitted Collection of Bipartite SBM") {
      cat(type, "--", self$distribution, "variant for", self$M, "networks \n")
      cat("=====================================================================\n")
      cat("net_id = (", self$net_id, ")\n")
      cat(
        "Dimensions = (", toString(lapply(seq.int(self$M), function(m) {
          c(self$n[[1]][[m]], self$n[[2]][[m]])
        })), ") - (",
        toString(self$Q), ") blocks.\n"
      )
      cat(
        "BICL = ", self$BICL,
        "\n#Empty row blocks on all networks: ", sum(!self$Cpi[[1]]),
        " -- #Empty columns blocks on all networks: ", sum(!self$Cpi[[2]]), " \n"
      )
      cat("* Useful fields \n")
      cat("  $distribution, $nb_nodes, $nb_blocks, $support, $prob_memberships \n")
      cat("  $memberships, $parameters, $BICL, $vbound, $pred_dyads \n")
      cat("=====================================================================")
    },
    #' The print method
    #'
    #' @return nothing; print to console
    print = function() self$show(),
    #' @title The method to plot the fitBipartite objects
    #'
    #' @import ggplot2
    #' @import dplyr
    #' @importFrom matrixStats rowMeans2
    #' @importFrom patchwork plot_layout plot_annotation wrap_plots
    #' @importFrom reshape2 melt
    #' @importFrom purrr map_dfc
    #' 
    #' @param type The type of the plot. Could be "graphon", "meso" or "block".
    #' @param oRow A reordering of the row blocks.
    #' @param oCol A reordering of the column blocks.
    #' @param mixture Should the block proportions of each network be plotted as
    #' well?
    #' @param net_id Use to plot only on network in "graphon" view.
    #' @param ... Further argument to be passed
    #' @return A plot, a ggplot2 object.
    plot = function(type = "graphon", oRow = NULL, oCol = NULL, mixture = FALSE, net_id = NULL, ...) {
      # The below order use mean over all networks to have a consistent display
      if (is.null(oRow)) {
        if (self$Q[2] == 1) {
          mean_rho <- 1
        } else {
          mean_rho <- matrixStats::rowMeans2(sapply(self$pi, function(pi) pi[[2]]))
        }
        oRow <- order(self$alpha %*% mean_rho, decreasing = TRUE)
      }
      if (is.null(oCol)) {
        if (self$Q[1] == 1) {
          mean_pi <- 1
        } else {
          mean_pi <- matrixStats::rowMeans2(sapply(self$pi, function(pi) pi[[1]]))
        }
        oCol <- order(mean_pi %*% self$alpha, decreasing = TRUE)
      }
      p <- switch(type,
        graphon = {
          if (self$Q[1] == 1) {
            ymin <- rep(0, each = self$Q[2])
            ymax <- rep(1, each = self$Q[2])
          } else {
            ymin <- rep(c(0, cumsum(self$pi[[net_id]][[1]][oRow][1:(self$Q[1] - 1)])), each = self$Q[2])
            ymax <- rep(c(cumsum(self$pi[[net_id]][[1]][oRow])), each = self$Q[2])
          }
          if (self$Q[2] == 1) {
            xmin <- rep(0, self$Q[1])
            xmax <- rep(1, self$Q[1])
          } else {
            xmin <- rep(c(0, cumsum(self$pi[[net_id]][[2]][oCol][1:(self$Q[2] - 1)])), self$Q[1])
            xmax <- rep(cumsum(self$pi[[net_id]][[2]][oCol]), self$Q[1])
          }
          (self$alpha[oRow, oCol] * mean(self$delta)) %>%
            t() %>%
            reshape2::melt() %>%
            dplyr::mutate(
              xmin = xmin,
              ymin = ymin,
              xmax = xmax,
              ymax = ymax
            ) %>%
            ggplot2::ggplot(ggplot2::aes(
              xmin = xmin, ymin = ymin,
              xmax = xmax, ymax = ymax, fill = value
            )) +
            ggplot2::geom_rect() +
            ggplot2::scale_fill_gradient2("alpha", low = "white", mid = "red", midpoint = 1) +
            ggplot2::geom_hline(yintercept = cumsum(self$pi[[net_id]][[1]][oRow][1:(self$Q[1] - 1)]), linewidth = .2) +
            ggplot2::geom_vline(xintercept = cumsum(self$pi[[net_id]][[2]][oCol][1:(self$Q[2] - 1)]), linewidth = .2) +
            ggplot2::scale_y_reverse() +
            ggplot2::theme_bw(base_size = 15, base_rect_size = 1, base_line_size = 1) +
            ggplot2::xlab("Column Blocks") +
            ggplot2::ylab("Row Blocks") +
            ggplot2::coord_equal(expand = FALSE)
        },
        meso = {
          p_alpha <- self$alpha[oRow, oCol, drop = FALSE] %>%
            t() %>%
            reshape2::melt() %>%
            ggplot2::ggplot(ggplot2::aes(x = Var1, y = Var2, fill = value)) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_gradient2("alpha", low = "white", high = "red") +
            ggplot2::geom_hline(yintercept = seq(self$Q[1]) + .5) +
            ggplot2::geom_vline(xintercept = seq(self$Q[2]) + .5) +
            ggplot2::scale_y_reverse() +
            ggplot2::theme_bw(base_size = 15, base_rect_size = 1, base_line_size = 1) +
            ggplot2::xlab("") +
            ggplot2::ylab("") +
            ggplot2::coord_fixed(expand = FALSE)
          #  scale_y_reverse()
          if (self$free_density) {
            xl <- paste(round(self$delta, 1))
          } else {
            xl <- ""
          }
          df_pi <- purrr::map_dfc(
            seq_along(self$net_id),
            function(m) setNames(data.frame(self$pim[[m]][[1]][oRow]), self$net_id[m])
          )
          df_rho <- purrr::map_dfc(
            seq_along(self$net_id),
            function(m) setNames(data.frame(self$pim[[m]][[2]][oCol]), 
              self$net_id[m])
          )
          # names(df_pi) <- self$net_id
          if (mixture) {
            p_pi <-
              df_pi %>%
              #    rename() %>%
              dplyr::mutate(q = seq(self$Q[1])) %>%
              tidyr::pivot_longer(cols = -c(q)) %>%
              dplyr::mutate(Proportion = value) %>%
              ggplot2::ggplot(ggplot2::aes(
                fill = as.factor(q), y = name,
                x = Proportion
              )) +
              ggplot2::geom_col() +
              ggplot2::coord_flip(expand = FALSE) +
              ggplot2::scale_fill_brewer("Row block",
                type = "qual", palette = "Paired",
                direction = -1
              ) +
              ggplot2::guides(fill = ggplot2::guide_legend(
                ncol = self$Q[1] %/% 3 + 1,
                byrow = TRUE
              )) +
              ggplot2::ylab("") +
              ggplot2::ylab(xl) +
              ggplot2::theme(axis.text.x = ggplot2::element_text(
                angle = 90, vjust = .5,
                hjust = 1
              ))
            ggplot2::theme_bw(base_size = 15)
            p_rho <- df_rho %>%
              #    rename() %>%
              dplyr::mutate(q = seq(self$Q[2])) %>%
              tidyr::pivot_longer(cols = -c(q)) %>%
              dplyr::mutate(Proportion = value) %>%
              ggplot2::ggplot(ggplot2::aes(
                fill = as.factor(q), y = name,
                x = Proportion
              )) +
              ggplot2::geom_col() +
              # ggplot2::coord_flip(expand = FALSE) +
              ggplot2::scale_fill_brewer("Column block",
                type = "qual", palette = "Set2",
                direction = -1
              ) +
              ggplot2::guides(fill = ggplot2::guide_legend(
                ncol = self$Q[2] %/% 3 + 1,
                byrow = TRUE
              )) +
              ggplot2::ylab("") +
              ggplot2::ylab(xl) +
              ggplot2::theme_bw(base_size = 15)
            # Merging the plots with patchwork
            mixture_layout <- "
                              ##CCCC
                              ##CCCC
                              RRAAAA
                              RRAAAA
                              RRAAAA
                              "
            p_alpha <- patchwork::wrap_plots(
              R = p_pi, C = p_rho, A = p_alpha,
              design = mixture_layout
            ) +
              patchwork::plot_layout(
                guides = "collect",
                design = mixture_layout
              )
          }
          return(p_alpha)
        },
        "block" = {
          as.matrix(self$A[[net_id]])[
            order(self$Z[[net_id]][[1]]),
            order(self$Z[[net_id]][[2]])
          ] %>%
            reshape2::melt() %>%
            ggplot2::ggplot(ggplot2::aes(x = Var2, y = rev(Var1), fill = value)) +
            ggplot2::geom_tile(show.legend = FALSE) +
            ggplot2::geom_hline(
              yintercept = cumsum(tabulate(self$Z[[net_id]][[1]])[self$Q[1]:2]) + .5,
              col = "red", size = .5
            ) +
            ggplot2::geom_vline(
              xintercept = cumsum(tabulate(self$Z[[net_id]][[2]])[1:(self$Q[2] - 1)]) + .5,
              col = "red", size = .5
            ) +
            ggplot2::scale_fill_gradient(low = "white", high = "black", na.value = "transparent") +
            ggplot2::ylab("") +
            ggplot2::xlab(self$net_id[net_id]) +
            ggplot2::scale_x_discrete(
              breaks = ""
            ) +
            # ggplot2::scale_y_reverse() +
            ggplot2::scale_y_discrete(
              breaks = "",
              guide = ggplot2::guide_axis(angle = 0)
            ) +
            ggplot2::coord_equal(expand = FALSE) +
            ggplot2::theme_bw(base_size = 15) +
            ggplot2::theme(axis.ticks = ggplot2::element_blank())
        }
      )
      return(p)
    }
  ),
  active = list(
    #' @field nb_nodes Returns n a list of the number of nodes per network
    nb_nodes = function(value) self$n,
    #' @field nb_blocks Returns Q a vector with 2 coordinates, Q1 and Q2 for the
    #' row blocks and the column blocks 
    nb_blocks = function(value) self$Q,
    #' @field support Returns the Cpi, a list of M boolean matrices indicating 
    #' which blocks are populated
    support = function(value) self$Cpi,
    #' @field prob_memberships Returns the tau, the probabilities of memberships
    #' "a posteriori", after seeing the data
    prob_memberships = function(value) self$tau,
    #' @field parameters Returns the list of parameters of the model, alpha, pi 
    #' and rho
    parameters = function(value) {
      list(
        alpha = self$alpha,
        pi = lapply(self$pi, function(list_pi_rho) list_pi_rho[[1]]),
        rho = lapply(self$pi, function(list_pi_rho) list_pi_rho[[2]])
      )
    },
    #' @field pred_dyads Predicted dyads from the estimated probabilities and
    #' parameters
    pred_dyads = function(value) {
      lapply(
        seq(self$M),
        function(m) {
          A_hat <- self$tau[[m]][[1]] %*% (self$delta[m] * self$alpha) %*% t(self$tau[[m]][[2]])
          A_hat
        }
      )
    },
    #' @field memberships The block memberships
    memberships = function(value) {
      self$Z
    }
  )
)
