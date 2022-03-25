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
  select(-total)


my_rarefy <- function(x, sample){
  
  x <- x[x>0]
  sum(1-exp(lchoose(sum(x) - x, sample) - lchoose(sum(x), sample)))
  
}

min_n_seqs <- shared %>%
  group_by(Group) %>%
  summarize(n_seqs = sum(value)) %>%
  summarize(min = min(n_seqs)) %>%
  pull(min)

shared_df <- shared %>%
  pivot_wider(names_from="name", values_from="value", values_fill = 0) %>%
  as.data.frame()

rownames(shared_df) <- shared_df$Group
shared_df <- shared_df[,-1]

vegans <- rarefy(shared_df, min_n_seqs) %>% 
  as_tibble(rownames="Group") %>%
  select(Group, vegan=value)

mine <- shared %>%
  group_by(Group) %>%
  summarize(mine = my_rarefy(value, min_n_seqs))

inner_join(vegans, mine, by="Group") %>%
  mutate(diff = vegan-mine) %>%
  summarize(mean(diff), sd(diff))


rrarefy(shared_df, sample=min_n_seqs) %>%
  as_tibble(rownames="Group") %>%
  pivot_longer(-Group)


old <- options(pillar.sigfig = 7)
drarefy(shared_df, sample=min_n_seqs) %>%
  as_tibble(rownames="Group") %>%
  pivot_longer(-Group)

options(old)
drarefy(shared_df, sample=min_n_seqs) %>%
  as_tibble(rownames="Group") %>%
  pivot_longer(-Group)


rarecurve_data <- rarecurve(shared_df,  step=100)

map_dfr(rarecurve_data, bind_rows) %>% 
  bind_cols(Group = rownames(shared_df),.) %>%
  pivot_longer(-Group) %>%
  drop_na() %>%
  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>%
  select(-name) %>%
  ggplot(aes(x=n_seqs, y=value, group=Group)) +
  geom_line(color="gray") +
  geom_vline(xintercept = min_n_seqs, color="red") +
  theme_classic()

ggsave("vegan_rarefaction.tiff")
