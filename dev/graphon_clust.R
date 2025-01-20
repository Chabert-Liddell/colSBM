devtools::load_all()

data(dorebipartite)
netlist <- dorebipartite[1:3]
colsbm_model <- "iid"
net_id <- NULL
distribution <- "bernoulli"
nb_run <- 3L
global_opts <- list()
fit_opts <- list()

nr <- 75
nc <- 75

current_model <- "iid"
pi <- matrix(c(0.2, 0.3, 0.5), nrow = 1, byrow = TRUE)
rho <- matrix(c(0.2, 0.3, 0.5), nrow = 1, byrow = TRUE)
eps <- 0.3

alpha_assortative <- matrix(0.3, nrow = 3, ncol = 3) +
    matrix(
        c(
            eps, -0.5 * eps, -0.5 * eps,
            -0.5 * eps, eps, -0.5 * eps,
            -0.5 * eps, -0.5 * eps, eps
        ),
        nrow = 3, byrow = TRUE
    )

alpha_core_periphery <- matrix(0.3, nrow = 3, ncol = 3) +
    matrix(
        c(
            1.5 * eps, eps, 0.5 * eps,
            eps, 0.5 * eps, 0,
            0.5 * eps, 0, -0.5 * eps
        ),
        nrow = 3, byrow = TRUE
    )

alpha_disassortative <- matrix(0.3, nrow = 3, ncol = 3) +
    matrix(
        c(
            -0.5 * eps, eps, eps,
            eps, -0.5 * eps, eps,
            eps, eps, -0.5 * eps
        ),
        nrow = 3, byrow = TRUE
    )

assortative_collection <- generate_bipartite_collection(
    nr, nc,
    pi, rho,
    alpha_assortative, 3,
    model = current_model,
    return_memberships = TRUE
)

assortative_incidence <- lapply(
    seq_along(assortative_collection),
    function(m) {
        return(assortative_collection[[m]]$incidence_matrix)
    }
)

assortative_row_clustering <- lapply(
    seq_along(assortative_collection),
    function(m) {
        return(assortative_collection[[m]]$row_clustering)
    }
)

assortative_col_clustering <- lapply(
    seq_along(assortative_collection),
    function(m) {
        return(assortative_collection[[m]]$row_clustering)
    }
)

core_periphery_collection <- generate_bipartite_collection(
    nr, nc,
    pi, rho,
    alpha_core_periphery, 3,
    model = current_model,
    return_memberships = TRUE
)

core_periphery_incidence <- lapply(
    seq_along(core_periphery_collection),
    function(m) {
        return(core_periphery_collection[[m]]$incidence_matrix)
    }
)

core_periphery_row_clustering <- lapply(
    seq_along(core_periphery_collection),
    function(m) {
        return(core_periphery_collection[[m]]$row_clustering)
    }
)

core_periphery_col_clustering <- lapply(
    seq_along(core_periphery_collection),
    function(m) {
        return(core_periphery_collection[[m]]$row_clustering)
    }
)

disassortative_collection <- generate_bipartite_collection(
    nr, nc,
    pi, rho,
    alpha_disassortative, 3,
    model = current_model,
    return_memberships = TRUE
)

disassortative_incidence <- lapply(
    seq_along(disassortative_collection),
    function(m) {
        return(disassortative_collection[[m]]$incidence_matrix)
    }
)

disassortative_row_clustering <- lapply(
    seq_along(disassortative_collection),
    function(m) {
        return(disassortative_collection[[m]]$row_clustering)
    }
)

disassortative_col_clustering <- lapply(
    seq_along(disassortative_collection),
    function(m) {
        return(disassortative_collection[[m]]$row_clustering)
    }
)

real_row_clustering <- append(
    append(
        assortative_row_clustering,
        core_periphery_row_clustering
    ),
    disassortative_row_clustering
)

real_col_clustering <- append(
    append(
        assortative_col_clustering,
        core_periphery_col_clustering
    ),
    disassortative_col_clustering
)

incidence_matrices <- append(
    append(
        assortative_incidence,
        core_periphery_incidence
    ),
    disassortative_incidence
)

netids <- rep(c("as", "cp", "dis"), each = 3)
