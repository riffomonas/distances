file <- scan("mice_simple.braycurtis.dist",
                  what=character(),
                  quiet=TRUE,
                  sep="\n")

n_samples <- as.numeric(file[1])
file <- file[-1]

file_split <- strsplit(file, "\t")

fill_in <- function(x, length){
  c(x, rep("0", length - length(x)))
}

filled <- lapply(file_split, fill_in, length=n_samples + 1)

#samples_distance_matrix <- matrix(unlist(filled), nrow=n_samples, byrow=TRUE)
samples_distance_matrix <- do.call(rbind, filled)

samples <- samples_distance_matrix[,1]

dist_matrix <- samples_distance_matrix[,-1]
dist_matrix <- matrix(as.numeric(dist_matrix), nrow=n_samples)
dist_matrix <- dist_matrix + t(dist_matrix)

