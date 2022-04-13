library(tidyverse)
library(vegan)

days_wanted <- c(9, 141)

shared_tbl <- read_tsv("data/mice.shared") %>%
  mutate(days = as.numeric(str_replace(Group, ".*D", ""))) %>%
  filter(days %in% days_wanted) %>%
  select(Group, starts_with("Otu"))

shared_df <- shared_tbl %>%
  column_to_rownames("Group")

# generate distance matrix - bray-curtis/rarefied to 18XX
mice_dist <- avgdist(shared_df, sample=1828)

mice_dist %>%
  as.matrix() %>%
  as_tibble(rownames = "samples") %>%
  pivot_longer(-samples) %>%
  filter(samples < name) %>%
  separate(samples, into=c("animal_a", "day_a"), "D", convert=TRUE) %>%
  separate(name, into=c("animal_b", "day_b"), "D", convert=TRUE) %>%
  mutate(comparison = case_when(
    animal_a != animal_b & day_a == 9 & day_a == day_b ~ "early",
    animal_a != animal_b & day_a == 141 & day_a == day_b ~ "late",
    animal_a == animal_b & day_a == 141 & day_b == 9 ~ "same",
    TRUE ~ NA_character_),
    comparison = factor(comparison, levels=c("same", "early", "late"))
    ) %>% 
  drop_na() %>%
  ggplot(aes(x=comparison, y=value)) +
  geom_jitter(width=0.25, color="gray") +
  stat_summary(fun.data=median_hilow, color="red", size=1,
               fun.args = list(conf.int=0.50)) +
  labs(x=NULL, y="Bray-Curtis distances")+
  scale_x_discrete(breaks=c("early", "late", "same"),
                   labels=c("Inter-mouse\ndistances at\n9 dpw",
                            "Inter-mouse\ndistances at\n141 dpw",
                            "Intra-mouse\ndistances between\n9 and 141 dpw")) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  theme_classic()

ggsave("jittered_dist_plot.png", width=4, height=4)
