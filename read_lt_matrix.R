distances <- scan("mice_simple.braycurtis.dist",
                  what=character(),
                  quiet=TRUE,
                  sep="\t")

n_samples <- as.numeric(distances[1])
distances <- distances[-1]

dist_matrix <- matrix(0, nrow=n_samples, ncol=n_samples)
samples <- rep("", n_samples)

samples[1] <- distances[1]
distances <- distances[-1]


samples[2] <- distances[1]
distances <- distances[-1]

dist_matrix[2, 1] <- as.numeric(distances[1])
distances <- distances[-1]


samples[3] <- distances[1]
distances <- distances[-1]

dist_matrix[3, 1:2] <- as.numeric(distances[1:2])
distances <- distances[-c(1,2)]


samples[4] <- distances[1]
distances <- distances[-1]

dist_matrix[4, 1:3] <- as.numeric(distances[1:3])
distances <- distances[-(1:3)]


samples[5] <- distances[1]
distances <- distances[-1]

dist_matrix[5, 1:4] <- as.numeric(distances[1:4])
distances <- distances[-(1:4)]


samples[6] <- distances[1]
distances <- distances[-1]

dist_matrix[6, 1:5] <- as.numeric(distances[1:5])
distances <- distances[-(1:5)]


samples[7] <- distances[1]
distances <- distances[-1]

dist_matrix[7, 1:6] <- as.numeric(distances[1:6])
distances <- distances[-(1:6)]


samples[8] <- distances[1]
distances <- distances[-1]

dist_matrix[8, 1:7] <- as.numeric(distances[1:7])
distances <- distances[-(1:7)]


samples[9] <- distances[1]
distances <- distances[-1]

dist_matrix[9, 1:8] <- as.numeric(distances[1:8])
distances <- distances[-(1:8)]


samples[10] <- distances[1]
distances <- distances[-1]

dist_matrix[10, 1:9] <- as.numeric(distances[1:9])
distances <- distances[-(1:9)]
