library(tidyverse)
source("code/read_matrix.R")

dist_matrix <- read_matrix("data/mice.braycurtis.dist")

dist_tbl <- as_tibble(dist_matrix, rownames="samples")

sample_lookup <- dist_tbl %>% 
  select(samples) %>%
  mutate(animal = str_replace(samples, "D\\d*", ""),
         sex = str_replace(samples, "\\d*D\\d*", ""),
         day = as.numeric(str_replace(samples, "\\w\\d*D", "")))
