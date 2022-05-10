#' An R6 Class object, a fitted population of netowrks sbm
#' once $optimize() is done
#'

fitBipartiteSBMPop <- R6::R6Class(
  classname = "fitBipartiteSBMPop",
  # inherit = "bmpop",
  #
  public = list(
    nr = NULL, # a vector of size M
    nc = NULL, # a vector of size M
    M = NULL, # Number of networks
    A = NULL, # List of nr[m]xnc[m] Matrix
    mask = NULL, # List of NA Matrices
    nb_inter = NULL, # A vector of length M the number of unique non NA entries
    directed = NULL, # Boolean for network direction, Constant
    Q = NULL, # Number of clusters, vectors of size2
    tau = NULL, # List of size M variational parameters n[m]xQ matrices
    alpha = NULL, # Matrix of size QxQ, connection parameters
    delta = NULL, # Vector of M,  density parameters with delta[1] = 1
    pi = NULL, # List of M vectors of size Q, the mixture parameters
    e = NULL, # Vector of size M, the sum of unique entries
    emqr = NULL, # List of M QxQ matrix, the sum of edges between q and r in m
    nmqr = NULL, # list of M QxQ matrix, the number of entries between q and r in m
    alpham = NULL, # list of M QxQ matrix, the classic sbm parameters
    free_mixture = NULL, # A boolean
    free_density = NULL, # A boolean
    weight = NULL, # A vector of size M for weighted likelihood
    model = NULL, # "poisson", "bernoulli"
    mloss = Inf,
    vloss = NULL,
    vbound = NULL,
    net_id = NULL,
    df_mixture = NULL,
    df_connect = NULL,
    df_density = NULL,
    logfactA = NULL,
    #  algo_ve = NULL,
    #  minibatch = NULL,
    init_method = NULL,
    #  verbosity = NULL,
    penalty = NULL,
    # approx_pois = NULL,
    Z = NULL,
    map = NULL,
    map_parameters = NULL,
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
                          model = "bernoulli",
                          free_mixture = TRUE,
                          free_density = TRUE,
                          directed = NULL,
                          init_method = "spectral",
                          weight = NULL, # A vector of size M, the weight of each network
                          fit_opts = list(algo_ve = "fp",
                                          approx_pois = TRUE,
                                          minibatch = TRUE,
                                          verbosity = 1)) {
      if (is.null(directed)) {
        directed <- ! isSymmetric.matrix(A[[1]])
      }
      self$directed <- directed
      self$A <- A
      self$M <- length(A)
      self$nr <- vapply(seq_along(A), function(m) nrow(A[[m]]), FUN.VALUE = .1)
      self$nc <- vapply(seq_along(A), function(m) ncol(A[[m]]), FUN.VALUE = .1)
      self$e <- vapply(seq_along(A), function(m) sum(A[[m]], na.rm = TRUE),
                       FUN.VALUE = .1)
      if (! is.null(mask)) {
        self$mask <- mask
      } else {
        self$mask <- lapply(
          seq_along(self$A),
          function(m) {
            mask <-  matrix(1, nrow(self$A), ncol(self$A))
            if(sum(is.na(self$A[[m]] > 0))) {
              mask[is.na(self$A[[m]])] <- 0
            }
            mask
          }
        )
      }
      if(is.null(net_id)) {
        self$net_id <- seq(self$M)
      } else {
        self$net_id <- net_id
      }
      self$Q <- Q
      self$model <- model
      if (self$model == "poisson") {
        self$logfactA <- vapply(
          seq_along(self$A),
          function(m) {
            sum(lfactorial(self$A[[m]],na.rm = TRUE))
          },
          FUN.VALUE = .1)
      }
      lapply(seq_along(self$A),
             function(m) self$A[[m]][is.na(self$A[[m]])] <- -1)
      self$fit_opts <- list(algo_ve = "fp",
                            approx_pois = TRUE,
                            minibatch = TRUE,
                            verbosity = 1)
      self$fit_opts <- modifyList(self$fit_opts, fit_opts)
      self$free_mixture <- free_mixture
      self$free_density <- free_density
      self$weight <- weight
      self$pi <- vector("list", self$M)
      self$tau <- vector("list", self$M)
      self$alpham <- vector("list", self$M)
      self$vloss <-  vector("list", self$M)
      self$emqr <- array(1, dim = c(self$M, self$Q[1], self$Q[2]))
      self$nmqr <- array(1, dim = c(self$M, self$Q[1], self$Q[2]))
      self$map <- list()
      self$map$emqr <- self$emqr
      self$map$nmqr <- self$nmqr
      self$map$pi <- self$pi
      self$df_mixture <- self$Q - 1
      self$df_density <- self$M - 1
      self$df_connect <- self$Q[1]*self$Q[2]
      self$Z <- if (is.null(Z)) {
        vector("list", self$M)
      } else {Z}
      self$delta <- rep(1, self$M)
      if(self$free_density) {
        self$delta <- (self$e/(self$nr*self$nc))/
          (self$e[1]/((self$nr[1]*self$nc[1])))
      }
      self$alpha <- matrix(.5, Q[1], Q[2])
      self$init_method <- init_method
      self$nb_inter <- self$dircoef *
        vapply(seq(self$M), function(m) sum(self$mask[[m]]), FUN.VALUE = .1)
      self$vbound <- vector("list", self$M)
    },


    compute_map = function() {
      self$Z <- lapply(self$tau,
                       function(tau)
                         list(row = apply(tau[[1]], 1, which.max),
                              col = apply(tau[[2]], 1, which.max)))
      invisible(self$Z)
    },


    objective = function() {
      sum(vapply(seq_along(self$A), function(m) vb(m)), FUN.VALUE = .1)
    },


    vb_tau_alpha = function(m, map = FALSE) {
      if (map) {
        emqr <- self$map$emqr[m,,]
        nmqr <- self$map$nmqr[m,,]
        alpha <- self$map$alpha
        delta <- self$map$delta[m]
      } else {
        emqr <- self$emqr[m,,]
        nmqr <- self$nmqr[m,,]
        alpha <- self$alpha
        delta <- self$delta[m]
      }
      switch(
        self$model,
        "bernoulli" = {
          #tau_tmp <- self$tau[[m]]
          # self$dircoef*sum(
          #   self$M[[m]] *
          #     (self$A[[m]] *
          #        tcrossprod(tau_tmp %*% .logit(self$alpha*self$delta[m], eps=1e-6),
          #                   tau_tmp) +
          #        tcrossprod(tau_tmp %*% .log(1-self$alpha*self$delta[m], eps = 1e-6),
          #                   tau_tmp)))
         sum(
            .xlogy(emqr, alpha*delta, eps = 1e-12) +
              .xlogy(nmqr - emqr, 1 - alpha*delta, eps = 1e-12))
        },
        "poisson" = {
          #          tau_tmp <- self$tau[[m]]
          sum(.xlogy(emqr, alpha*delta, eps = 1e-12) -
                               nmqr * alpha*delta -
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


    vb_tau_pi = function(m, map = FALSE) {
      if (! map) {
        sum(self$tau[[m]][[1]] %*% log(self$pi[[m]][[1]])) +
          sum(self$tau[[m]][[2]] %*% log(self$pi[[m]][[2]]))
      } else {
        sum(.xlogy(tabulate(self$Z[[m]][[1]], self$Q[1]), self$map$pi[[m]][[1]])) +
        sum(.xlogy(tabulate(self$Z[[m]][[2]], self$Q[2]), self$map$pi[[m]][[2]]))
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
                 -.5*sum(emqr[m,,]*.log(alpha*delta[m], eps=1e-9) +
                           (nmqr[m,,] - emqr[m,,]) * .log(1-alpha*delta[m], eps = 1e-9))
               },
               FUN.VALUE = .1)
      sum(res)
    },


    gr_vb_alpha_delta = function(par, emqr, nmqr) {
      # browser()
      alpha <- par[1:self$df_connect]
      delta <- c(1, par[1:self$df_density + self$df_connect])
      res_alpha <-
        lapply(seq_along(self$A),
               function(m) {
                 emqr[m,,]/pmax(alpha, 1e-9) -
                   (nmqr[m,,] - emqr[m,,]) *
                   (delta[m]/pmax(1-alpha*delta[m], 1e-9))
               })
      res_alpha <- Reduce( "+",res_alpha)
      res_delta <-
        vapply(seq_along(self$A),
               function(m) {
                 sum(emqr[m,,]/(pmax(delta[m], 1e-9)) -
                       (nmqr[m,,] - emqr[m,,]) *
                       (delta[m]/pmax(1-alpha*delta[m], 1e-9)))
               },
               FUN.VALUE = .1)
      if(self$directed) {
        res <- -c(as.vector(res_alpha), res_delta[2:self$M])
      }
      else {
        res <- -c(as.vector(res_alpha[lower.tri(res_alpha, diag = TRUE)]),
                  res_delta[2:self$M])
      }
      invisible(res)
    },

    eval_g0_vb_alpha_delta = function(par, emqr, nmqr) {
      as.vector(
        vapply(
          seq(self$M-1),
          function(m) {
            c(par[1 : self$df_density + self$df_connect])[m]  *
              par[1 : self$df_connect] - 1+1e-9
          },
          FUN.VALUE = rep(.1,self$df_connect)
        )
      )
    },

    eval_jac_g0_vb_alpha_delta = function(par, emqr, nmqr) {
      jac_d <- aperm(
        vapply(
          seq(self$df_density),
          function(m) diag(par[m + self$df_connect],self$df_connect ),
          FUN.VALUE = matrix(.1, self$df_connect, self$df_connect)),
        perm = c(2,3,1))
      dim(jac_d) <- c((self$df_density)*self$df_connect, self$df_connect)
      jac_a <- matrix(0,(self$df_density)*self$df_connect, self$df_density)
      for (m in seq(self$df_density)) {
        jac_a[1:(self$df_connect) + (m-1)*(self$df_connect), m] <-
          par[1:(self$df_connect)]
      }
      cbind(jac_d, jac_a)
    },

    update_alpha_delta = function(map = FALSE) {
      # browser()
      d <- self$delta
      a <- self$alpha
      a <- pmin(a*d[1], 1-1e-12)
      d <- d/d[1]
      if (map) {
        emqr <- self$map$emqr
        nmqr <- self$map$nmqr
      } else {
        emqr <- self$emqr
        nmqr <- self$nmqr
      }
      hat <- nloptr::nloptr(
        x0 = c(a, d[2:self$M]),
        eval_f =  self$fn_vb_alpha_delta,
        eval_grad_f =  self$gr_vb_alpha_delta,
        lb = c(rep(10*.Machine$double.eps,
                   self$df_connect + self$df_density)),
        ub = c(rep(1-10*.Machine$double.eps,
                   self$df_connect),
               rep(Inf, self$df_density)),
        eval_g_ineq =  self$eval_g0_vb_alpha_delta,
        eval_jac_g_ineq = self$eval_jac_g0_vb_alpha_delta,
        opts = list("algorithm" = "NLOPT_LD_MMA",
                    # "local_opts" = list("algorithm" = "NLOPT_LD_LBFGS",
                    #                     "xtol_rel" = 1.0e-4),
                    "xtol_rel" = 1.0e-4),
        emqr = emqr, nmqr = nmqr)
      a <- matrix(hat$solution[1:self$df_connect], self$Q[1], self$Q[2])
      d <- c(1, hat$solution[1:(self$df_density)+self$df_connect])
      if (map) {
        self$map$alpha <- a
        self$map$delta <- d
      } else {
        self$alpha <- a
        self$delta <- d
      }
    },


    compute_vbound = function() {
      sum(vapply(
        seq_along(self$A),
        function(m) self$entropy_tau(m) + self$vb_tau_pi(m) + self$vb_tau_alpha(m),
        FUN.VALUE = .1))
    },


    compute_penalty = function() {
      df_connect <- self$df_connect
      if (self$free_density) df_connect <- self$df_connect + self$df_density
      self$penalty  <- .5*(df_connect*log(sum(self$nb_inter)) +
                             sum(self$df_mixture *log(self$n)))
      invisible(self$penalty)
    },


    compute_icl = function(map = FALSE) {
      ICL <-
        sum(vapply(
          seq_along(self$A),
          function(m) self$vb_tau_pi(m, map = map) + self$vb_tau_alpha(m, map = map),
          FUN.VALUE = .1)) - self$compute_penalty()
      if(map) {
        self$map$ICL <- ICL
      } else {
        self$ICL <- ICL
      }
      return(ICL)
    },
################################################################################
    ## A modifier
    compute_icl_clustering = function(map = TRUE) {
      # browser()
      Z_unique <- lapply(self$Z, function(Z) sort(unique(Z[[1]])))
      uz <- unique(Z_unique)
      self$net_clustering <- match(Z_unique, uz)
      #   nb_netclusters <- max(self$net_clustering)
      pen <- 0
      for (g in seq_along(uz)) {
        Q <- length(uz[[g]]) # + self$Q - length(unique(unlist(uz)))
        M <- sum(self$net_clustering == g)
        df_mixture <- Q-1
        df_connect <- Q * Q
        if (self$free_density) df_connect <- df_connect + M - 1
        pen <- pen + .5*(df_connect*log(sum(self$nb_inter[self$net_clustering == g])) +
                           sum(df_mixture *log(self$n[self$net_clustering == g])))
      }
      self$penalty_clustering <- pen
      self$ICL_clustering <-
        sum(vapply(
          seq_along(self$A),
          function(m) self$vb_tau_pi(m, map = map) + self$vb_tau_alpha(m, map = map),
          FUN.VALUE = .1)) - pen
      invisible(self$ICL_clustering)
    },


    compute_exact_icl = function() {
      ## directed not implemented yet
      df_mixture <- ifelse(self$free_mixture,(self$Q-1), self$Q-1)
      df_connect <-
        ifelse(self$directed, self$Q**2, self$Q*(self$Q+1)/2)
      if (self$free_density) df_connect <- df_connect + self$M - 1
      eqr <- colSums(self$emqr, 1)
      nqr <- colSums(self$nmqr, 1)
      if (! self$directed) {
        diag(eqr) <- diag(eqr)/2
        diag(nqr) <- diag(nqr)/2
        exicl <-       - df_mixture * lbeta(.5,.5) +
          self$M * (lgamma(.5*self$Q) - self$Q * lgamma(.5)) +
          sum(lbeta(.5 + eqr, .5 + nqr - eqr)[lower.tri(eqr, diag = TRUE)]) +
          Reduce(sum,
                 lapply(
                   seq_along(self$A),
                   function(m) {
                     sum(lgamma(self$n[m]*self$pi[[m]] + .5)) -
                       lgamma(sum(self$n[m]*self$pi[[m]] + .5))
                   }
                 )
          )
      } else {
        exicl <-       df_mixture * (1/lbeta(.5,.5)) +
          self$M * (lgamma(.5*self$Q) - self$Q * lgamma(.5)) +
          sum(lbeta(.5 + eqr, .5 + nqr - eqr)) +
          Reduce(sum,
                 lapply(
                   seq_along(self$A),
                   function(m) {
                     sum(lgamma(self$n[m]*self$pi[[m]] + .5)) -
                       lgamma(sum(self$n[m]*self$pi[[m]] + .5))
                   }
                 )
          )
      }
      return(exicl)
    },


    compute_exact_icl_iid = function() {
      ## directed not implemented yet
      df_mixture <- ifelse(self$free_mixture, self$M*(self$Q-1), self$Q-1)
      emqr <- self$emqr
      nmqr <- self$nmqr
      for( m in seq_along(self$A)) {
        diag(emqr[m,,]) <- diag(emqr[m,,])/2
        diag(nmqr[m,,]) <- diag(nmqr[m,,])/2
      }
      Reduce(
        sum,
        lapply(
          seq_along(self$A),
          function(m)  {
            -df_mixture *lbeta(.5,.5) +
              (lgamma(.5*self$Q) - self$Q * lgamma(.5)) +
              sum(lbeta(.5 + emqr[m,,], .5 + nmqr[m,,] - emqr[m,,])[lower.tri(emqr[m,,], diag = TRUE)]) +
              sum(lgamma(self$n[m]*self$pi[[m]] + .5)) -
              lgamma(sum(self$n[m]*self$pi[[m]] + .5))
          }
        )
      )
    },
################################################################################
    update_map_parameters = function() {
      Z <- self$compute_map()
      ZR <- lapply(
        seq_along(Z),
        function(m) {
          .one_hot(Z[[m]][[1]], self$Q[1])
        })
      ZC <- lapply(
        seq_along(Z),
        function(m) {
          .one_hot(Z[[m]][[2]], self$Q[2])
        })
      lapply(
        seq_along(Z),
        function(m) self$map$emqr[m,,] <-
          crossprod(ZR[[m]], (self$A[[m]] * self$mask[[m]]) %*% ZC[[m]])
      )
      lapply(
        seq_along(Z),
        function(m) self$map$nmqr[m,,] <-
          crossprod(ZR[[m]], self$mask[[m]] %*% ZC[[m]])
      )
      self$map$Z <- Z
      self$map$alpha <- self$alpha
      self$map$delta <- self$delta
      self$m_step(map = TRUE)
      if(!self$free_density) {
        self$map$delta <- self$delta
      }
      lapply(seq.int(self$M), function(m) self$update_alpham(m, map = TRUE))
      invisible(Z)
    },


    fixed_point_tau = function(m, d, max_iter = 10, tol = 1e-3) {
      condition <- TRUE
      it <- 0
      # reup_counter <- 0
      self$vloss[[m]] <-
        c(self$vloss[[m]], self$vb_tau_alpha(m) + self$vb_tau_pi(m) +
            self$entropy_tau(m))
      tau_old <- self$tau[[m]][[d]]
      tau_new <- switch (
        self$model,
        "bernoulli" = {
          tau_new <-
            if (d == 1) {
              tau_new <-
                t(matrix(log(self$pi[[m]][[d]]), self$Q[d], self$nr[m])) +
                (self$mask[[m]] * self$A[[m]]) %*%
                self$tau[[m]][[2]] %*%
                t(.logit(self$delta[m]*self$alpha, eps = 1e-9)) +
                self$mask[[m]] %*%
                self$tau[[m]][[2]] %*%
                t(.log(1-self$alpha*self$delta[m], eps = 1e-9))
            }
          if (d == 2) {
            tau_new <-
              t(matrix(log(self$pi[[m]][[d]]), self$Q[d], self$nc[m])) +
              (self$mask[[m]] * self$A[[m]]) %*%
              self$tau[[m]][[1]] %*%
              .logit(self$delta[m]*self$alpha, eps = 1e-9) +
              self$mask[[m]] %*%
              self$tau[[m]][[1]] %*%
              .log(1-self$alpha*self$delta[m], eps = 1e-9)
          }
          invisible(tau_new)
        },
        "poisson" = {
          if (d == 1) {
            tau_new <-
              t(matrix(log(self$pi[[m]][[d]]), self$Q[d], self$nr[m])) +
              (self$mask[[m]] * self$A[[m]]) %*%
              self$tau[[m]][[2]] %*%
              t(log(self$delta[m]*self$alpha)) -
              self$mask[[m]] %*%
              self$tau[[m]][[2]] %*%
              t(self$alpha*self$delta[m])
          }
          if (d == 2) {
            tau_new <-
              t(matrix(log(self$pi[[m]][[d]]), self$Q[d], self$nc[m])) +
              (self$mask[[m]] * self$A[[m]]) %*%
              self$tau[[m]][[1]] %*%
              t(log(self$delta[m]*self$alpha)) -
              self$mask[[m]] %*%
              self$tau[[m]][[1]] %*%
              t(self$alpha*self$delta[m])
          }
          invisible(tau_new)
        }
      )
      tau_new <- .softmax(tau_new)
      tau_new[tau_new < 1e-9] <- 1e-9
      tau_new[tau_new > 1-1e-9] <- 1-1e-9
      tau_new <- tau_new/rowSums(tau_new)
      self$tau[[m]][[d]] <- tau_new
      invisible(tau_new)
    },


    fixed_point_alpha_delta = function(map = FALSE, max_iter = 50, tol = 1e-6) {
      # switch(
      #   self$model,
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
        d <- rowSums(emqr, dim = 1)/
          rowSums(aperm(array(a, c(self$Q[1], self$Q[2], self$M))) * nmqr, dim = 1)
        d[1] <- 1
        a <- colSums(emqr, dim = 1)/
          colSums(array(d, c(self$M, self$Q[1], self$Q[2])) * nmqr, dim = 1)
        a[is.nan(a)] <- 0
        it <- it+1
        condition <- mean((a - a_old)**2) + mean((d - d_old)**2) > 2*tol &
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


    update_pi = function(m, map = FALSE) {
      if (! map) {
        pi1 <- .colMeans(self$tau[[m]][[1]], self$nr[m], self$Q[1])
        pi2 <- .colMeans(self$tau[[m]][[2]], self$nr[m], self$Q[2])
        self$pi[[m]] <- list(row = pi1, col = pi2)
      } else {
        pi1 <- tabulate(self$Z[[m]][[1]], self$Q[1])/self$nr[m]
        pi2 <- tabulate(self$Z[[m]][[2]], self$Q[2])/self$nc[m]
        self$map$pi[[m]] <- list(row = pi1, col = pi2)
      }
      invisible(pi)
    },


    update_alpham = function(m, map = FALSE) {
      if (! map) {
        alpham <-  self$emqr[m,,]/self$nmqr[m,,]
        alpham[is.nan(alpham)] <- 0
        self$alpham[[m]] <- alpham
      } else {
        alpham <- self$map$emqr[m,,]/self$map$nmqr[m,,]
        alpham[is.nan(alpham)] <- 0
        self$map$alpham[[m]] <- alpham
      }
      invisible(alpham)
    },


    update_alpha = function(map = FALSE) {
      if (! map) {
        alpha <- colSums(self$emqr, dims = 1)/colSums(self$nmqr, dims = 1)
        alpha[is.nan(alpha)] <- 0
        self$alpha <- alpha
      } else {
        alpha <- colSums(self$map$emqr, dims = 1)/colSums(self$map$nmqr, dims = 1)
        alpha[is.nan(alpha)] <- 0
        self$map$alpha <- alpha
      }
      invisible(alpha)
    },


    init_clust = function() {
      self$tau <-
        switch(
          self$init_method,
          "random" = lapply(
            X = seq_along(self$A),
            FUN = function(m) {
              tau_tmp <- gtools::rdirichlet(self$n[m], rep(1, self$Q))
              self$emqr[m,,] <- .tquadform(tau_tmp, self$A[[m]]*self$mask[[m]])
              self$nmqr[m,,] <- .tquadform(tau_tmp, self$mask[[m]])
              a <- self$emqr[m,,]/self$nmqr[m,,]
              prob <- self$Q*diag(a) #+ rowSums(a)
              p <- sample.int(self$Q, prob = prob)
              tau_tmp <- tau_tmp[,p]
              self$emqr[m,,] <- self$emqr[m,p,p]
              self$nmqr[m,,] <- self$nmqr[m,p,p]
              tau_tmp
            }
          ),
          "spectral" = lapply(
            X = seq_along(self$A),
            FUN = function(m) {
              tau_tmp <- .one_hot(spectral_clustering(self$A[[m]], self$Q), self$Q)
              tau_tmp[tau_tmp < 1e-6] <- 1e-6
              tau_tmp[tau_tmp > 1-1e-6] <- 1-1e-6
              tau_tmp <- tau_tmp / .rowSums(tau_tmp, self$n[m], self$Q)
              self$emqr[m,,] <- .tquadform(tau_tmp, self$A[[m]]*self$mask[[m]])
              self$nmqr[m,,] <- .tquadform(tau_tmp, self$mask[[m]])
              a <- self$emqr[m,,]/self$nmqr[m,,]
              prob <- self$Q*diag(a)# + rowSums(a)
              p <- sample.int(self$Q, prob = prob)
              tau_tmp <- tau_tmp[,p]
              self$emqr[m,,] <- self$emqr[m,p,p]
              self$nmqr[m,,] <- self$nmqr[m,p,p]
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
                tau_tmp[tau_tmp > 1-1e-6] <- 1-1e-6
                tau_tmp <- tau_tmp / .rowSums(tau_tmp, self$n[m], self$Q)
                self$emqr[m,,] <- .tquadform(tau_tmp, self$A[[m]]*self$mask[[m]])
                self$nmqr[m,,] <- .tquadform(tau_tmp, self$mask[[m]])
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
      if(self$model == "bernoulli" & self$free_density &
         ! self$fit_opts$approx_pois)
        self$fixed_point_alpha_delta()
    },


    make_permutation = function() {

    },


    m_step = function(map = FALSE, max_iter = 100, tol = 1e-3,...) {
      #browser()
      lapply(seq_along(self$pi), function(m) self$update_pi(m, map = map))
      if (self$free_density == FALSE) {
        self$update_alpha(map = map)
      } else {
        switch(
          self$model,
          "poisson" = self$fixed_point_alpha_delta(map = map),
          "bernoulli" =
            ifelse(self$fit_opts$approx_pois,
                   self$fixed_point_alpha_delta(map = map),
                   self$update_alpha_delta(map = map))
        )
      }
    },



    ve_step = function(m,max_iter = 20, tol = 1e-3,...) {
      # Place holder for gradient ascent or other optimization methods
    },

    update_mqr = function(m) {
      tau_tmp <- self$tau[[m]]
      self$emqr[m,,] <- .tquadform(tau_tmp, self$A[[m]] * self$mask[[m]])
      self$nmqr[m,,] <- .tquadform(tau_tmp, self$mask[[m]])
    },


    optimize = function(max_step = 100, tol = 1e-3, ...) {
      if (self$Q == 1) {
        self$tau <- lapply(seq(self$M), function(m) matrix(1, self$n[m], 1))
        self$Z <- lapply(seq(self$M), function(m) rep(1, self$n[m]))
        self$pi <- lapply(seq(self$M), function(m) 1)
        lapply(seq(self$M), function(m) self$update_mqr(m))
        if(self$free_density) {
          self$alpha <- matrix(sum(self$emqr[1,,])/sum(self$nmqr[1,,]), 1, 1)
          self$delta <- vapply(
            seq(self$M),
            function(m) {
              sum(self$emqr[m,,])/sum(self$alpha * self$nmqr[m,,])
            }, FUN.VALUE = .1)} else {
              self$alpha <- matrix(sum(self$emqr)/sum(self$nmqr), 1, 1)
              self$delta <- rep(1, self$M)
            }
        self$alpham <- lapply(seq(self$M),
                              function(m)  matrix(sum(self$emqr[m,,])/sum(self$nmqr[m,,]), 1, 1))
        self$map <- list(Z = self$tau,
                         pi = self$pi,
                         alpha = self$alpha,
                         delta = self$delta,
                         alpham = self$alpham,
                         emqr = self$emqr,
                         nmqr = self$nmqr)
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
          if (! self$fit_opts$minibatch) {
            lapply(seq(self$M),
                   function(m) {
                     switch (self$fit_opts$algo_ve,
                             "fp" = self$fixed_point_tau(m),
                             self$ve_step(m,...))
                     self$update_mqr(m)
                   }
            )
            self$m_step(...)
          } else {
            seq_m <- sample.int(self$M)
            lapply(seq(self$M),
                   function(m) {
                     switch (self$fit_opts$algo_ve,
                             "fp" = self$fixed_point_tau(seq_m[m]),
                             self$ve_step(seq_m[m],...))
                     self$update_mqr(seq_m[m])
                     self$m_step(...)
                   })
          }
          vb <- self$compute_vbound()
          self$vbound <- c(self$vbound, vb)
          step <- step + 1
          step_condition <- step < max_step &
            self$vbound[length(self$vbound)] -
            self$vbound[length(self$vbound)-1] > tol
          if (step %% 5 == 0) {
            if(self$fit_opts$verbosity >= 1) {
              print(paste0(step, ": ", vb))
              print(self$alpha)
            }
          }
        }
        self$compute_map()
        lapply(seq(self$M), function(m) self$update_alpham(m))
        # self$compute_parameters()
        self$compute_icl()
        self$update_map_parameters()
      }
      self$compute_icl()
      self$compute_icl(map = TRUE)
      self$compute_icl_clustering()
    }
  ),
  active = list(
    dircoef =  function() ifelse (self$directed, 1, .5)
  )
)
