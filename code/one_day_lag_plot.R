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
  as_tibble(rownames="samples") %>%
  pivot_longer(-samples) %>%
  filter(samples < name) %>%
  mutate(animal_a = str_replace(samples, "D.*", ""),
         animal_b = str_replace(name, "D.*", ""),
         day_a = str_replace(samples, ".*D", "") %>% as.double(),
         day_b = str_replace(name, ".*D", "") %>% as.double(),
         diff = abs(day_a-day_b),
         day = if_else(day_b > day_a, day_b, day_a),
         period = if_else(day < 10, "early", "late"),
         sex = if_else(str_detect(animal_a, "F"), "female", "male")) %>%
  filter(animal_a == animal_b & diff == 1) %>%
  ggplot(aes(x=day, y=value, group=animal_a, color=sex)) +
  geom_line() +
  facet_wrap(~period, scales="free_x", nrow=1) +
  geom_smooth(aes(group=period), se=FALSE, color="black", size=2) +
  scale_x_continuous(breaks=1:150) +
  scale_color_manual(name=NULL,
                     breaks=c("female", "male"),
                     labels=c("Female", "Male"),
                     values=c("purple", "limegreen")) +
  labs(x="Days after weaning",
       y="Bray-Curtis distance to previous day") +
  theme_classic() +
  theme(
    strip.text = element_blank(),
    legend.position = c(0.8, 0.9),
    axis.text.x = element_text(size=7)
  )

ggsave("one_day_lag_plot.png", width=5, height=3)
