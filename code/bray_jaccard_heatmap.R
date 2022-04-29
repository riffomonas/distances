library(tidyverse)
library(vegan)

days_wanted <- c(0:9, 141:150)

early_late_tbl <- read_tsv("data/mice.shared") %>%
  select(Group, starts_with("Otu")) %>%
  pivot_longer(-Group) %>%
  separate(Group, into=c("animal", "day"), sep="D",
           remove=FALSE, convert=TRUE) %>%
  filter(day %in% days_wanted) %>%
  group_by(Group) %>%
  mutate(N = sum(value)) %>%
  ungroup() %>%
  filter(N >= 1828) %>%
  select(-N) %>%
  pivot_wider(names_from="name", values_from="value", values_fill=0)
