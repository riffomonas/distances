library(tidyverse)
source("code/read_matrix.R")

dist_matrix <- read_matrix("data/mice.braycurtis.dist")

dist_tbl <- as_tibble(dist_matrix, rownames="samples")

sample_lookup <- dist_tbl %>% 
  select(samples) %>%
  mutate(delimited = str_replace(samples,
                                 "^(([FM])\\d+)D(\\d+)$",
                                 "\\2-\\1-\\3")) %>%
  separate(col=delimited,
           into=c("sex", "animal", "day"), sep="-",
           convert=TRUE)
