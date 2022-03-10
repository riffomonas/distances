library(vegan)
library(zCompositions)
library(tidyverse)

days_wanted <- c(0:9, 141:150)

shared <- read_tsv("data/mice.shared") %>%
  select(Group, starts_with("Otu")) %>%
  mutate(day = str_replace(Group, ".*D", "")) %>%
  filter(day %in% days_wanted) %>%
  select(-day) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total = sum(value)) %>%
  filter(total > 1800) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total)

rand <- shared %>%
  uncount(value) %>%
  mutate(rand_name = sample(name)) %>%
  select(-name) %>%
  count(Group, rand_name)

group_count <- rand %>%
  group_by(Group) %>%
  summarize(n_seqs = sum(n))

rand_df <- rand %>%
  pivot_wider(names_from=rand_name, values_from=n, values_fill = 0) %>%
  as.data.frame()

rownames(rand_df) <- rand_df$Group
rand_df <- rand_df[,-1]

rare_dist <- avgdist(rand_df, dmethod="euclidean", sample=min(group_count$n_seqs))
