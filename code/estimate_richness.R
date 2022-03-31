library(vegan)
library(breakaway)
library(tidyverse)

set.seed(19760620)

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
  mutate(name = sample(name)) %>%
  count(Group, name, name="value")

get_breakaway <- function(x){
  
  ba <- breakaway(x)
  tibble(est= ba$estimate,
         lci=ba$interval[1], uci=ba$interval[2],
         model=ba$model)
  
}

get_chao <- function(x){
  
  sobs <- sum(x$n)
  sing <- x[x$value == 1, "n"] %>% pull(n)
  doub <- x[x$value == 2, "n"] %>% pull(n)
  
  sobs + sing^2 / (2*doub)
  
}

b_analysis <- rand %>%
  count(Group, value) %>%
  nest(data = -Group) %>%
  mutate(n_seqs = map_dbl(data, ~sum(.x$value * .x$n)),
         sobs = map_dbl(data, ~sum(.x$n)),
         rare = map_dbl(data, ~rarefy(rep(.x$value, .x$n), sample=1828)),
         chao = map_dbl(data, ~get_chao(.x)),
         ba = map(data, ~get_breakaway(.x))) %>%
  select(-data) %>%
  unnest(ba)


b_analysis %>%
  select(Group, n_seqs, sobs, chao, rare, est) %>%
  pivot_longer(-c(Group, n_seqs)) %>%
  ggplot(aes(x=n_seqs, y=value, color=name)) +
  geom_point() +
  geom_smooth() +
  coord_cartesian(ylim=c(0, 2000))
