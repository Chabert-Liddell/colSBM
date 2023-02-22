# Sourcing all necessary files
require("sbm", quietly = T)
require("aricode", quietly = T)
require("dplyr", quietly = T)
require("tictoc", quietly = T)
# require("visNetwork", quietly = T)
# require("igraph", quietly = T)
source("R/utils.R")
source("R/R6class-fitBipartiteSBMPop.R")

set.seed(1234)

ARImethod <- function(fit_col) {
    out <- list()
    out <- lapply(seq.int(fit_col$M), function(m) {
        list(
            row = aricode::ARI(
                fit_col$MAP$Z[[m]][[1]],
                Z[[m]][[1]]
            ),
            col = aricode::ARI(
                fit_col$MAP$Z[[m]][[2]],
                Z[[m]][[2]]
            )
        )
    })
}

verbose <- TRUE
test_alea <- TRUE

eps <- 0.05
M <- 5
nr <- 100
nc <- 250

pir <- c(0.2, 0.8)
pic <- c(0.2, 0.3, 0.5)

Q <- c(length(pir), length(pic))

alpha <- matrix(
    c(
        0.9, eps, eps,
        eps, 0.8, eps
    ), nrow = Q[1], ncol = Q[2], byrow = TRUE
)

output_file = file("output.txt", "w+")

cat(
    "New run : \nParameters are : \n - M : ", M,
    "\n - Pi row : ", pir,
    "\n - Pi col : ", pic,
    "\n - alpha : ", alpha,
    "\n"
)

cat(
    "\nParameters are : \n - M : ", M,
    "\n - Pi row : ", pir,
    "\n - Pi col : ", pic,
    "\n - alpha : ", alpha,
    "\n",
    file = output_file,
    append = TRUE
)

bipartite_collection <- generate_bipartite_collection(nr, nc, pir, pic, alpha, M)

# This is a list of the M incidence matrices
bipartite_collection_incidence <- lapply(seq.int(M), function(m) {
    bipartite_collection[[m]]$incidence_matrix
})


## Init given with exact membership

Z <- lapply(seq.int(M), function(m) {
    list(bipartite_collection[[m]]$row_clustering, bipartite_collection[[m]]$col_clustering)
})

fitColExactMembership <- fitBipartiteSBMPop$new(
    A = bipartite_collection_incidence,
    Q = Q, free_mixture = FALSE,
    Z = Z,
    free_density = FALSE,
    init_method = "given",
    fit_opts = list(verbosity = ifelse(verbose, 4, 0))
)

fitColExactMembership.init_resultARI <- mutate(
    as.data.frame(do.call("rbind", lapply(
        invisible(seq.int(M)),
        function(m) {
            ARIrow <- aricode::ARI(
                bipartite_collection[[m]]$row_clustering,
                Z[[m]][[1]]
            )
            ARIcol <- aricode::ARI(
                bipartite_collection[[m]]$col_clustering,
                Z[[m]][[2]]
            )
            list(network = m, ARIrow = ARIrow, ARIcol = ARIcol)
        }
    ))),
    network = unlist(network),
    ARIrow = unlist(ARIrow),
    ARIcol = unlist(ARIcol),
    row_error = (ARIrow - 1)^2,
    col_error = (ARIcol - 1)^2
)

fitColExactMembership$optimize()
fitColExactMembership$MAP

fitColExactMembership.end_resultARI <- mutate(
    as.data.frame(do.call("rbind", lapply(
        invisible(seq.int(M)),
        function(m) {
            ARIrow <- ARImethod(fitColExactMembership)[[m]]$row
            ARIcol <- ARImethod(fitColExactMembership)[[m]]$col
            list(network = m, ARIrow = unlist(ARIrow), ARIcol = ARIcol)
        }
    ))),
    network = unlist(network),
    ARIrow = unlist(ARIrow),
    ARIcol = unlist(ARIcol),
    row_error = (ARIrow - 1)^2,
    col_error = (ARIcol - 1)^2
)

cat("Given exact clusters :\n - Real alpha : ", alpha, "\n - Given alpha : ", fitColExactMembership$alpha, "\n")

iter_max <- 1
if (test_alea) {
    for (i in 1:iter_max) {
        if (i %% 5 == 0) {
            cat("Seed ", i, "\n Completed : ", (i / iter_max * 100), "%")
        }
        set.seed(i)

        tic(paste0("Seed ", i))
        ## Init with given LBM
        sepLBM <- lapply(
            seq_along(bipartite_collection),
            function(m) {
                estimateBipartiteSBM(bipartite_collection_incidence[[m]],
                    estimOptions = list(
                        verbosity = ifelse(verbose, 1, 0),
                        plot = FALSE
                    )
                )
            }
        )

        # Extracting the list of the memberships
        Z_LBM <- lapply(seq.int(M), function(m) {
            list(sepLBM[[m]]$memberships$row, sepLBM[[m]]$memberships$col)
        })

        LBM_reorganized_matrices <- lapply(seq.int(M), function(m) {

        })

        # Create the fitBipartite object
        fitColSepLBM <- fitBipartiteSBMPop$new(
            A = bipartite_collection_incidence,
            Q = Q, free_mixture = FALSE,
            free_density = FALSE,
            Z = Z_LBM,
            init_method = "given",
            fit_opts = list(verbosity = ifelse(verbose, 4, 0))
        )


        # Checking the clustering obtained
        fitColSepLBM.init_resultARI <- mutate(
            as.data.frame(do.call("rbind", lapply(
                invisible(seq.int(M)),
                function(m) {
                    ARIrow <- aricode::ARI(
                        bipartite_collection[[m]]$row_clustering,
                        sepLBM[[m]]$memberships$row
                    )
                    ARIcol <- aricode::ARI(
                        bipartite_collection[[m]]$col_clustering,
                        sepLBM[[m]]$memberships$col
                    )
                    list(network = m, ARIrow = ARIrow, ARIcol = ARIcol)
                }
            ))),
            network = unlist(network),
            ARIrow = unlist(ARIrow),
            ARIcol = unlist(ARIcol),
            row_error = (ARIrow - 1)^2,
            col_error = (ARIcol - 1)^2
        )
        if (verbose) {
            cat("SepLBM Clustering\n")
            cat("Init:")
        }
        table.SepLBM.init <- knitr::kable(fitColSepLBM.init_resultARI)
        if (verbose) print(table.SepLBM.init)

        # SepLBM optimize
        fitColSepLBM$optimize()
        if (verbose) cat("End:")
        fitColSepLBM.end_resultARI <- mutate(
            as.data.frame(do.call("rbind", lapply(
                invisible(seq.int(M)),
                function(m) {
                    ARIrow <- ARImethod(fitColSepLBM)[[m]]$row
                    ARIcol <- ARImethod(fitColSepLBM)[[m]]$col
                    list(network = m, ARIrow = unlist(ARIrow), ARIcol = ARIcol)
                }
            ))),
            network = unlist(network),
            ARIrow = unlist(ARIrow),
            ARIcol = unlist(ARIcol),
            row_error = (ARIrow - 1)^2,
            col_error = (ARIcol - 1)^2
        )

        table.SepLBM.end <- knitr::kable(fitColSepLBM.end_resultARI)
        if (verbose) {
            print(table.SepLBM.end)

            cat(
                "SepLBM clustering colSBM :\n - Init :\n\t - Mean row error : ",
                mean(fitColSepLBM.init_resultARI$row_error),
                "\n\t - Mean col error : ", mean(fitColSepLBM.init_resultARI$col_error),
                "\n - After colSBM :\n\t - Mean row error : ",
                mean(fitColSepLBM.end_resultARI$row_error),
                "\n\t - Mean col error : ",
                mean(fitColSepLBM.end_resultARI$col_error)
            )
        }
        if (verbose) cat("\nSepLBM clustering :\n - Real alpha : ", alpha, "\n - SepLBM alpha : ", fitColSepLBM$alpha, "\n")

        # Init with Spectral

        # Create the fitBipartite object
        fitColSpectral <- fitBipartiteSBMPop$new(
            A = bipartite_collection_incidence,
            Q = Q, free_mixture = FALSE,
            free_density = FALSE,
            fit_opts = list(verbosity = ifelse(verbose, 4, 0))
        )

        # Spectral Clustering by row and cols separately
        fitColSpectral.spectralbiclust <- lapply(
            seq_along(fitColSpectral$A),
            function(m) {
                spectral_biclustering(X = fitColSpectral$A[[m]], K = fitColSpectral$Q)
            }
        )

        # Checking the clustering obtained
        fitColSpectral.init_resultARI <- mutate(
            as.data.frame(do.call("rbind", lapply(
                invisible(seq.int(M)),
                function(m) {
                    ARIrow <- aricode::ARI(
                        bipartite_collection[[m]]$row_clustering,
                        fitColSpectral.spectralbiclust[[m]]$row_clustering
                    )
                    ARIcol <- aricode::ARI(
                        bipartite_collection[[m]]$col_clustering,
                        fitColSpectral.spectralbiclust[[m]]$col_clustering
                    )
                    list(network = m, ARIrow = ARIrow, ARIcol = ARIcol)
                }
            ))),
            network = unlist(network),
            ARIrow = unlist(ARIrow),
            ARIcol = unlist(ARIcol),
            row_error = (ARIrow - 1)^2,
            col_error = (ARIcol - 1)^2
        )

        if (verbose) cat("\nSpectral Bi-Clustering\n")
        table.Spectral.init <- knitr::kable(fitColSpectral.init_resultARI)
        if (verbose) print(table.Spectral.init)

        # Spectral optimize
        fitColSpectral$optimize()
        if (verbose) cat("End:")
        fitColSpectral.end_resultARI <- mutate(
            as.data.frame(do.call("rbind", lapply(
                invisible(seq.int(M)),
                function(m) {
                    ARIrow <- ARImethod(fitColSpectral)[[m]]$row
                    ARIcol <- ARImethod(fitColSpectral)[[m]]$col
                    list(network = m, ARIrow = ARIrow, ARIcol = ARIcol)
                }
            ))),
            network = unlist(network),
            ARIrow = unlist(ARIrow),
            ARIcol = unlist(ARIcol),
            row_error = (ARIrow - 1)^2,
            col_error = (ARIcol - 1)^2
        )

        table.Spectral.end <- knitr::kable(fitColSpectral.end_resultARI)
        if (verbose) {
            print(table.Spectral.end)


            cat(
                "Spectral clustering colSBM :\n - Init :\n\t - Mean row error : ",
                mean(fitColSpectral.init_resultARI$row_error),
                "\n\t - Mean col error : ", mean(fitColSpectral.init_resultARI$col_error),
                "\n - After colSBM :\n\t - Mean row error : ",
                mean(fitColSpectral.end_resultARI$row_error),
                "\n\t - Mean col error : ",
                mean(fitColSpectral.end_resultARI$col_error)
            )


            cat("\nSpectral clustering :\n - Real alpha : ", alpha, "\n - Spectral alpha : ", fitColSpectral$alpha, "\n")
        }
        ## Init with HCA
        # Create the fitBipartite object
        fitColHCA <- fitBipartiteSBMPop$new(
            A = bipartite_collection_incidence,
            Q = Q, free_mixture = FALSE,
            free_density = FALSE,
            init_method = "hca",
            fit_opts = list(verbosity = ifelse(verbose, 4, 0))
        )

        # Hierarchical Clustering by row and cols separately
        fitColHCA.hierarchicalbiclust <- lapply(
            seq_along(fitColHCA$A),
            function(m) {
                bipartite_hierarchic_clustering(X = fitColHCA$A[[m]], K = fitColHCA$Q)
            }
        )

        # Checking the clustering obtained
        fitColHCA.init_resultARI <- mutate(
            as.data.frame(do.call("rbind", lapply(
                invisible(seq.int(M)),
                function(m) {
                    ARIrow <- aricode::ARI(
                        bipartite_collection[[m]]$row_clustering,
                        fitColHCA.hierarchicalbiclust[[m]]$row_clustering
                    )
                    ARIcol <- aricode::ARI(
                        bipartite_collection[[m]]$col_clustering,
                        fitColHCA.hierarchicalbiclust[[m]]$col_clustering
                    )
                    list(network = m, ARIrow = ARIrow, ARIcol = ARIcol)
                }
            ))),
            network = unlist(network),
            ARIrow = unlist(ARIrow),
            ARIcol = unlist(ARIcol),
            row_error = (ARIrow - 1)^2,
            col_error = (ARIcol - 1)^2
        )
        if (verbose) {
            cat("\nHierarchical Bi-Clustering\n")
        }
        table.HCA.init <- knitr::kable(fitColHCA.init_resultARI)
        if (verbose) print(table.HCA.init)

        fitColHCA$optimize()

        Z
        fitColHCA$MAP$Z

        fitColHCA.end_resultARI <- mutate(
            as.data.frame(do.call("rbind", lapply(
                invisible(seq.int(M)),
                function(m) {
                    ARIrow <- ARImethod(fitColHCA)[[m]]$row
                    ARIcol <- ARImethod(fitColHCA)[[m]]$col
                    list(network = m, ARIrow = ARIrow, ARIcol = ARIcol)
                }
            ))),
            network = unlist(network),
            ARIrow = unlist(ARIrow),
            ARIcol = unlist(ARIcol),
            row_error = (ARIrow - 1)^2,
            col_error = (ARIcol - 1)^2
        )
        table.HCA.end <- knitr::kable(fitColHCA.end_resultARI)
        if (verbose) {
            print(table.HCA.end)

            cat(
                "HCA clustering colSBM :\n - Init :\n\t - Mean row error : ",
                mean(fitColHCA.init_resultARI$row_error),
                "\n\t - Mean col error : ", mean(fitColHCA.init_resultARI$col_error),
                "\n - After colSBM :\n\t - Mean row error : ",
                mean(fitColHCA.end_resultARI$row_error),
                "\n\t - Mean col error : ",
                mean(fitColHCA.end_resultARI$col_error)
            )

            cat("\nHCA clustering :\n - Real alpha : ", alpha, "\n - HCA alpha : ", fitColHCA$alpha, "\n\n")
        }

        # Writing to a file
        cat("\n---------------------------\nSeed : ", i, "Time : ", toc(quiet = T)$callback_msg,
            file = output_file,
            append = TRUE
        )
        cat("\n--------\nSepLBM\n", file = output_file, append = TRUE)
        cat("Init", file = output_file, append = TRUE)
        capture.output(table.SepLBM.init,
            file = output_file,
            append = TRUE
        )
        cat(
            "Init mean error\nRow : ",
            mean(fitColSepLBM.init_resultARI$row_error),
            "\nCol : ", mean(fitColSepLBM.init_resultARI$col_error),
            file = output_file,
            append = TRUE
        )
        cat("\nAfter colSBM", file = output_file, append = TRUE)
        capture.output(table.SepLBM.end,
            file = output_file,
            append = TRUE
        )
        cat(
            "After colSBM mean error\nRow : ",
            mean(fitColSepLBM.end_resultARI$row_error),
            "\nCol : ", mean(fitColSepLBM.end_resultARI$col_error),
            file = output_file,
            append = TRUE
        )
        cat("\n--------\nSpectral\n", file = output_file, append = TRUE)
        cat("Init", file = output_file, append = TRUE)
        capture.output(table.Spectral.init,
            file = output_file,
            append = TRUE
        )
        cat(
            "Init mean error\nRow : ",
            mean(fitColSpectral.init_resultARI$row_error),
            "\nCol : ", mean(fitColSpectral.init_resultARI$col_error),
            file = output_file,
            append = TRUE
        )
        cat("\nAfter colSBM", file = output_file, append = TRUE)
        capture.output(table.Spectral.end,
            file = output_file,
            append = TRUE
        )
        cat(
            "After colSBM mean error\nRow : ",
            mean(fitColSpectral.end_resultARI$row_error),
            "\nCol : ", mean(fitColSpectral.end_resultARI$col_error),
            file = output_file,
            append = TRUE
        )
        cat("\n--------\nHCA\n", file = output_file, append = TRUE)
        cat("Init", file = output_file, append = TRUE)
        capture.output(table.HCA.init,
            file = output_file,
            append = TRUE
        )
        cat(
            "Init mean error\nRow : ",
            mean(fitColHCA.init_resultARI$row_error),
            "\nCol : ", mean(fitColHCA.init_resultARI$col_error),
            file = output_file,
            append = TRUE
        )
        cat("\nAfter colSBM", file = output_file, append = TRUE)
        capture.output(table.HCA.end,
            file = output_file,
            append = TRUE
        )
        cat(
            "After colSBM mean error\nRow : ",
            mean(fitColHCA.end_resultARI$row_error),
            "\nCol : ", mean(fitColHCA.end_resultARI$col_error),
            file = output_file,
            append = TRUE
        )
    }

    close(output_file)
}
# Add NAs to the collections
# bipartite_collection_NA <- bipartite_collection
# bipartite_collection_NA <- lapply(
#     seq_along(bipartite_collection_NA),
#     function(m) {
#         bipartite_collection_NA[[m]][sample(length(bipartite_collection_NA[[m]]), numberOfNAsPerNetwork)] <- rep(NA, numberOfNAsPerNetwork)
#         bipartite_collection_NA[[m]]
#     }
# )
# fitColWithNA <- fitBipartiteSBMPop$new(A = , Q = Q)
