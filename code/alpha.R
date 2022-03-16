library(vegan)
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

richness <- function(x){
  
  # r <- sum(x > 0)
  # return(r)
  
  sum(x>0)
}

shannon <- function(x){
  
  rabund <- x[x>0]/sum(x)
  -sum(rabund * log(rabund))
  
}

simpson <- function(x){
  
  n <- sum(x)
  
  # sum(x * (x-1) / (n * (n-1)))
  1 - sum((x/n)^2)
}

rand %>%
  group_by(Group) %>%
  summarize(sobs = richness(value),
            shannon = shannon(value),
            simpson = simpson(value),
            invsimpson = 1/simpson,
            n = sum(value)) %>%
  pivot_longer(cols=c(sobs, shannon, invsimpson, simpson),
               names_to="metric") %>%
  ggplot(aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~metric, nrow=4, scales="free_y")


rand %>%
  group_by(Group) %>%
  summarize(sobs = specnumber(value),
            shannon = diversity(value, index="shannon"),
            simpson = diversity(value, index="simpson"),
            invsimpson = 1/simpson,
            n = sum(value)) %>%
  pivot_longer(cols=c(sobs, shannon, invsimpson, simpson),
               names_to="metric") %>%
  ggplot(aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~metric, nrow=4, scales="free_y")
