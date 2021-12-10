distances <- scan("mice_simple.braycurtis.dist",
                  what=character(),
                  quiet=TRUE,
                  sep="\t")

n_samples <- as.numeric(distances[1])
distances <- distances[-1]
