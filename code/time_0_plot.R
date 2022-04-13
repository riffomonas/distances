library(tidyverse)
library(vegan)

days_wanted <- c(0:9, 141:150)

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
  as_tibble(rownames="samples") %>%
  pivot_longer(-samples) %>%
  filter(samples < name) %>%
  mutate(animal_a = str_replace(samples, "D.*", ""),
         animal_b = str_replace(name, "D.*", ""),
         day_a = str_replace(samples, ".*D", "") %>% as.numeric(),
         day_b = str_replace(name, ".*D", "") %>% as.numeric(),
         period = if_else(day_b < 10, "early", "late"),
         sex = if_else(str_detect(animal_a, "F"), "female", "male")) %>%
  filter(animal_a == animal_b & day_a == 0) %>%
  ggplot(aes(x=day_b, y=value, group=animal_a, color=sex)) +
  geom_line() +
  geom_smooth(aes(group=period), se=FALSE, color="black", size=2) +
  scale_x_continuous(breaks=1:150) +
  scale_color_manual(name=NULL, 
                     breaks=c("female", "male"),
                     values=c("purple", "limegreen"),
                     labels=c("Female", "Male")) +
  facet_wrap(~period, nrow=1, scales="free_x") +
  labs(x="Days following weaning",
       y="Bray-Curtis distance to\nday of weaning") +
  theme_classic() +
  theme(
    strip.text = element_blank(),
    legend.text = element_text(size=7),
    axis.text.x = element_text(size=7),
    legend.box.spacing = unit(0, "in"),
    panel.border = element_rect(color="black", fill=NA, size=1)
  )

ggsave("time_0_plot.png", width=6, height=3)
