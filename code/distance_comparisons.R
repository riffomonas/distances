library(tidyverse)
library(vegan)

source("code/read_matrix.R")

bray_tbl <- read_tsv("data/mice.shared") %>%
  select(Group, starts_with("Otu")) %>%
  column_to_rownames("Group") %>%
  avgdist(sample=1804, dmethod="bray") %>%
  as.matrix() %>%
  as_tibble(rownames="A") %>%
  pivot_longer(-A, names_to="B", values_to="bray")

jaccard_tbl <- read_tsv("data/mice.shared") %>%
  select(Group, starts_with("Otu")) %>%
  column_to_rownames("Group") %>%
  avgdist(sample=1804, dmethod="jaccard") %>%
  as.matrix() %>%
  as_tibble(rownames="A") %>%
  pivot_longer(-A, names_to="B", values_to="jaccard")

unweighted_tbl <- read_matrix("data/mice.unweighted.ave.dist") %>% 
  as_tibble(rownames= "A") %>%
  pivot_longer(-A, names_to="B", values_to="unweighted") %>%
  # mutate(A = str_trim(A),
  #        B = str_trim(B))
  #mutate_at(c("A", "B"), str_trim)
  mutate_if(is.character, str_trim)
  
weighted_tbl <- read_matrix("data/mice.weighted.ave.dist") %>% 
  as_tibble(rownames= "A") %>%
  pivot_longer(-A, names_to="B", values_to="weighted") %>%
  # mutate(A = str_trim(A),
  #        B = str_trim(B))
  #mutate_at(c("A", "B"), str_trim)
  mutate_if(is.character, str_trim)


combined <- inner_join(bray_tbl, jaccard_tbl, by=c("A", "B")) %>%
  inner_join(., unweighted_tbl, by=c("A", "B")) %>%
  inner_join(., weighted_tbl, by=c("A", "B"))

combined %>%
  filter(A < B) %>%
  ggplot(aes(x=bray, y=weighted)) +
  geom_point() +
  geom_smooth(se=FALSE)

ggsave("bray_weighted.png", width=6, height=5)

combined %>%
  filter(A < B) %>%
  ggplot(aes(x=jaccard, y=unweighted)) +
  geom_point() +
  geom_smooth(se=FALSE)

ggsave("jaccard_unweighted.png", width=6, height=5)

combined %>%
  filter(A < B) %>%
  ggplot(aes(x=jaccard, y=bray)) +
  geom_point() +
  geom_smooth(se=FALSE)

ggsave("jaccard_bray.png", width=6, height=5)

combined %>%
  filter(A < B) %>%
  ggplot(aes(x=unweighted, y=weighted)) +
  geom_point() +
  geom_smooth(se=FALSE)

ggsave("unweighted_unweighted.png", width=6, height=5)


bray_dist <- combined %>%
  select(A, B, bray) %>%
  pivot_wider(names_from=B, values_from=bray) %>%
  column_to_rownames("A") %>%
  as.dist()

jaccard_dist <- combined %>%
  select(A, B, jaccard) %>%
  pivot_wider(names_from=B, values_from=jaccard) %>%
  column_to_rownames("A") %>%
  as.dist()

unweighted_dist <- combined %>%
  select(A, B, unweighted) %>%
  pivot_wider(names_from=B, values_from=unweighted) %>%
  column_to_rownames("A") %>%
  as.dist()

weighted_dist <- combined %>%
  select(A, B, weighted) %>%
  pivot_wider(names_from=B, values_from=weighted) %>%
  column_to_rownames("A") %>%
  as.dist()

mantel(bray_dist, weighted_dist, method="spearman")$statistic
mantel(weighted_dist, unweighted_dist, method="spearman")$statistic
mantel(bray_dist, jaccard_dist, method="spearman")$statistic
mantel(jaccard_dist, unweighted_dist, method="spearman")$statistic
