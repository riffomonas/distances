library(vegan)
library(zCompositions)
library(tidyverse)

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
  mutate(rand_name = sample(name)) %>%
  select(-name) %>%
  count(Group, rand_name)

group_count <- rand %>%
  group_by(Group) %>%
  summarize(n_seqs = sum(n))

group_count$n_seqs %>% range

rand_df <- rand %>%
  pivot_wider(names_from=rand_name, values_from=n, values_fill = 0) %>%
  as.data.frame()

rownames(rand_df) <- rand_df$Group
rand_df <- rand_df[,-1]

norare_dist <- vegdist(rand_df, method="euclidean")
rare_dist <- avgdist(rand_df, dmethod="euclidean", sample=min(group_count$n_seqs))

gm <- function(x){
  
  exp(mean(log(x[x>0])))
  
}

rclr_df <- rand %>%
  group_by(Group) %>%
  mutate(rclr = log(n/gm(n))) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from=rand_name, values_from=rclr, values_fill=0) %>%
  as.data.frame()

rownames(rclr_df) <- rclr_df$Group
rclr_df <- rclr_df[,-1]

rclr_dist <- vegdist(rclr_df, method="euclidean")

zclr_df <- cmultRepl(rand_df, method="CZM", output="p-count") %>%
  as_tibble(rownames = "Group") %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(zclr = log(value/gm(value))) %>%
  ungroup() %>%
  select(-value) %>%
  pivot_wider(names_from=name, values_from=zclr, values_fill=0) %>%
  column_to_rownames("Group")


zclr_dist <- vegdist(zclr_df, method="euclidean")

norare_dtbl <- norare_dist %>%
  as.matrix %>%
  as_tibble(rownames = "Group") %>%
  pivot_longer(cols= -Group) %>%
  filter(name < Group)

rare_dtbl <- rare_dist %>%
  as.matrix %>%
  as_tibble(rownames = "Group") %>%
  pivot_longer(cols= -Group) %>%
  filter(name < Group)

rclr_dtbl <- rclr_dist %>%
  as.matrix %>%
  as_tibble(rownames = "Group") %>%
  pivot_longer(cols=-Group) %>%
  filter(name < Group)

zclr_dtbl <- zclr_dist %>%
  as.matrix %>%
  as_tibble(rownames = "Group") %>%
  pivot_longer(cols=-Group) %>%
  filter(name < Group)

inner_join(norare_dtbl, rare_dtbl, by=c("Group", "name")) %>%
  inner_join(., rclr_dtbl, by=c("Group", "name")) %>%
  inner_join(., zclr_dtbl, by=c("Group", "name")) %>%
  inner_join(., group_count, by=c("Group" = "Group")) %>%
  inner_join(., group_count, by=c("name" = "Group")) %>%
  mutate(diffs = abs(n_seqs.x - n_seqs.y)) %>%
  select(Group, name, norare=value.x, rare=value.y, rclr=value.x.x, zclr=value.y.y, diffs) %>%
  pivot_longer(cols=c(norare, rare, rclr, zclr), names_to="method", values_to="dist") %>%
  ggplot(aes(x=diffs, y=dist)) +
  geom_point() +
  facet_wrap(~method, nrow=4, scales="free_y") +
  geom_smooth()

