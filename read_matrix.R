read_matrix <- function(file_name){
  
  file <- scan(file_name,
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
  
  samples_distance_matrix <- do.call(rbind, filled)
  
  samples <- samples_distance_matrix[,1]
  
  dist_matrix <- samples_distance_matrix[,-1]
  dist_matrix <- matrix(as.numeric(dist_matrix), nrow=n_samples)
  
  if(sum(dist_matrix[upper.tri(dist_matrix)]) == 0){
    dist_matrix <- dist_matrix+t(dist_matrix)
  }
  
  rownames(dist_matrix) <- samples
  colnames(dist_matrix) <- samples

  return(dist_matrix)
}