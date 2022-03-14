# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

library(vegan)
library(zCompositions)
library(DESeq2)
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

rand_matrix <- t(as.matrix(rand_df)) + 1

vst_matrix <- varianceStabilizingTransformation(rand_matrix, fitType = "mean") %>%
  t()

vst_dist <- vegdist(vst_matrix, method="euclidean")

vst_dtbl <- vst_dist %>%
  as.matrix() %>%
  as_tibble(rownames="samples") %>%
  pivot_longer(-samples) %>%
  filter(name < samples)

rare_dtbl <- rare_dist %>%
  as.matrix() %>%
  as_tibble(rownames="samples") %>%
  pivot_longer(-samples) %>%
  filter(name < samples)

inner_join(rare_dtbl, vst_dtbl, by=c("samples", "name")) %>%
  inner_join(., group_count, by=c("samples" = "Group")) %>%
  inner_join(., group_count, by=c("name" = "Group")) %>%
  mutate(diff = abs(n_seqs.x - n_seqs.y)) %>%
  select(samples, name, rare=value.x, vst=value.y, diff) %>%
  pivot_longer(cols=c("rare", "vst"), names_to="method", values_to="dist") %>%
  ggplot(aes(x=diff, y=dist)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~method, scales="free_y", nrow=2)
