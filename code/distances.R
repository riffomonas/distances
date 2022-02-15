library(tidyverse)
library(vegan)

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
  select(-total) %>%
  pivot_wider(Group) %>%
  as.data.frame()

rownames(shared) <- shared$Group
shared <- shared[, -1]
shared <- as.matrix(shared)

set.seed(19760620)
dist <- avgdist(shared, dmethod="bray", sample=1800)

set.seed(1)
nmds <- metaMDS(dist)

scores(nmds) %>%
  as_tibble(rownames = "Group") %>%
  ggplot(aes(x=NMDS1, y=NMDS2)) + 
  geom_point()
