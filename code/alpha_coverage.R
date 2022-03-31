library(vegan)
library(tidyverse)

set.seed(19760620)

days_wanted <- c(0:9, 141:150)

shared <- read_tsv("data/mice.shared") %>%
  select(Group, starts_with("Otu")) %>%
  mutate(day = str_replace(Group, ".*D", "")) %>%
#  filter(day %in% days_wanted) %>%
  select(-day) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total = sum(value)) %>%
  # filter(total > 1800) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total)

rand <- shared %>%
  uncount(value) %>%
  mutate(name = sample(name)) %>%
  count(Group, name, name="value")

sampling_coverage <- shared %>% 
  group_by(Group) %>%
  summarize(n_seqs = sum(value))

sampling_coverage %>%
  ggplot(aes(x=n_seqs)) +
  geom_histogram(binwidth=500) +
  coord_cartesian(xlim=c(0, 5000))

sampling_coverage %>%
  ggplot(aes(x=1, y=n_seqs)) +
  geom_jitter() +
  scale_y_log10()

sampling_coverage %>%
  arrange(n_seqs) %>%
  ggplot(aes(x=1:nrow(.), y=n_seqs))+
  geom_line() +
  coord_cartesian(xlim=c(0, 50), ylim=c(0, 5000))

sampling_coverage %>%
  arrange(n_seqs) %>%
  print(n=20)


coverage_stats <- shared %>%
  group_by(Group) %>%
  summarize(n_seqs = sum(value),
            n_sings = sum(value == 1),
            goods = 100*(1 - n_sings / n_seqs)) %>%
  filter(n_seqs > 1800)

coverage_stats %>%
  ggplot(aes(x=n_seqs, y=goods)) +
  geom_point()

coverage_stats %>% 
  arrange(goods)
