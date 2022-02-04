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

days_wanted <- c(0:9, 141:150)

dist_matrix <- dist_tbl %>%
  pivot_longer(cols=-samples, names_to="b", values_to="distances") %>%
  inner_join(., sample_lookup, by="samples") %>%
  inner_join(., sample_lookup, by=c("b" = "samples")) %>%
  filter(day.x %in% days_wanted & day.y %in% days_wanted) %>%
  select(samples, b, distances) %>%
  pivot_wider(names_from="b", values_from="distances") %>%
  select(-samples) %>%
  as.dist()
