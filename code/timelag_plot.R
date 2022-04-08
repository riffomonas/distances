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

mice_dist %>%
  as.matrix() %>%
  as_tibble(rownames= "sample") %>%
  pivot_longer(-sample) %>%
  filter(sample < name) %>%
  mutate(animal_a = str_replace(sample, "D.*", ""),
         animal_b = str_replace(name, "D.*", ""),
         day_a = as.numeric(str_replace(sample, ".*D", "")),
         day_b = as.numeric(str_replace(name, ".*D", "")),
         diff = abs(day_a - day_b),
         early = day_a < 10) %>%
  filter(animal_a == animal_b & diff < 10) %>%
  group_by(diff, animal_a, early) %>%
  summarize(median = median(value)) %>%
  ungroup() %>%
  ggplot(aes(x=diff, y=median, color=early, group=paste0(animal_a, early))) +
  geom_line(size=0.25) +
  geom_smooth(aes(group=early), se=FALSE, size=4) +
  labs(x="Days between time points",
       y="Median Bray-Curtis distance") +
  scale_x_continuous(breaks=1:9) +
  scale_color_manual(name=NULL,
                     breaks=c(TRUE, FALSE),
                     values=c("blue", "red"),
                     labels=c("Early", "Late")) +
  guides(color = guide_legend(override.aes = list(size=1))) +
  theme_classic()

ggsave("time_interval.png", width=5, height=3)
