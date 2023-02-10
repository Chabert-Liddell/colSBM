set.seed(1234)
coupleOfNodesPerClasses <- c(10, 25)
Q <- c(3, 6)
numberOfNetworks <- 5
numberOfNAsPerNetwork <- 10

# Generate a collection
col <- generate_bipartite_collection(coupleOfNodesPerClass = coupleOfNodesPerClasses, coupleOfNumberOfClasses = Q, numberOfNetworks = numberOfNetworks)

# Add NAs to the collections
colWithNA <- col
colWithNA$A <- lapply(
    seq_along(colWithNA$A),
    function(m) {
        colWithNA$A[[m]][sample(length(colWithNA$A[[m]]), numberOfNAsPerNetwork)] <- rep(NA, numberOfNAsPerNetwork)
        colWithNA$A[[m]]
    }
)

# Create the fitBipartite objects
fitCol <- fitBipartiteSBMPop$new(A = col$A, Q = Q)
fitColWithNA <- fitBipartiteSBMPop$new(A = colWithNA$A, Q = Q)

fitCol$optimize()
