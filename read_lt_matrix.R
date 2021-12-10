distances <- scan("mice.braycurtis.dist",
                  what=character(),
                  quiet=TRUE,
                  sep="\t")

n_samples <- as.numeric(distances[1])
distances <- distances[-1]

dist_matrix <- matrix(0, nrow=n_samples, ncol=n_samples)
samples <- rep("", n_samples)

samples[1] <- distances[1]
distances <- distances[-1]

for(i in 2:n_samples){
  samples[i] <- distances[1]
  distances <- distances[-1]
  
  dist_matrix[i, 1:(i-1)] <- as.numeric(distances[1:(i-1)])
  dist_matrix[1:(i-1), i] <- dist_matrix[i, 1:(i-1)]
  distances <- distances[-(1:(i-1))]
}