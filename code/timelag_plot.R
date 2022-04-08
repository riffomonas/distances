library(tidyverse)
library(vegan)

days_wanted <- c(0:9, 141:150)

shared_tbl <- read_tsv("data/mice.shared") %>%
  mutate(days = as.numeric(str_replace(Group, ".*D", ""))) %>%
  filter(days %in% days_wanted) %>%
  select(Group, starts_with("Otu"))

shared_tbl %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  summarize(n_seqs = sum(value)) %>%
  arrange(n_seqs) %>%
  print(n=15)

shared_df <- shared_tbl %>%
  column_to_rownames("Group")

# generate distance matrix - bray-curtis/rarefied to 18XX
mice_dist <- avgdist(shared_df, sample=1828)
