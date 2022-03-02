library(tidyverse)
library(vegan)
library(SRS)

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


shared_group_count <- shared %>%
  group_by(Group) %>%
  summarize(n = sum(value))

rand_group_count <- rand %>%
  group_by(Group) %>%
  summarize(n = sum(n))

inner_join(shared_group_count, rand_group_count,  by="Group")

smallest_group <- min(rand_group_count$n)

shared_otu_count <- shared %>%
  group_by(name) %>%
  summarize(n = sum(value))

rand_otu_count <- rand %>%
  group_by(rand_name) %>%
  summarize(n = sum(n))

inner_join(shared_otu_count, rand_otu_count,  by=c("name" = "rand_name" ))


rand_group_count %>%
  ggplot(aes(x=n)) + geom_histogram()

range(rand_group_count$n)

rand_df <- rand %>%
  pivot_wider(names_from="rand_name", values_from="n", values_fill = 0) %>%
  as.data.frame()

rownames(rand_df) <- rand_df$Group
rand_df <- rand_df[, -1]

rand_matrix <- as.matrix(rand_df)

norare_dist_matrix <- vegdist(rand_matrix, method="bray")
rare_dist_matrix <- avgdist(rand_matrix, dmethod="bray", sample=smallest_group)

norare_dist_tibble <- norare_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample)

rare_dist_tibble <- rare_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample)

relabund <- rand %>%
  group_by(Group) %>%
  mutate(rel_abund = n/sum(n)) %>%
  ungroup() %>%
  select(Group, rand_name, rel_abund) %>%
  pivot_wider(names_from="rand_name", values_from="rel_abund",values_fill=0) %>%
  as.data.frame()

rownames(relabund) <- relabund$Group
relabund <- relabund[, -1]
relabund_matrix <- as.matrix(relabund)

relabund_dist_matrix <- vegdist(relabund_matrix, method="bray")

relabund_dist_tibble <- relabund_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample)

# relabund %>%
#   as_tibble(rownames="Group") %>%
#   pivot_longer(-Group) %>%
#   mutate(scaled = round(value * smallest_group, 0)) %>%
#   group_by(Group) %>%
#   summarize(n_scaled = sum(scaled))

normalized <- rand_df %>%
  t() %>%
  as.data.frame() %>%
  SRS(Cmin=smallest_group) %>%
  t()

normalized_dist_matrix <- vegdist(normalized, method="bray")

normalized_dist_tibble <- normalized_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample)



comparison <- inner_join(norare_dist_tibble, rare_dist_tibble,
                         by=c("sample", "name")) %>%
  inner_join(., relabund_dist_tibble, by=c("sample", "name")) %>%
  inner_join(., normalized_dist_tibble, by=c("sample", "name")) %>%
  select(sample, name, norare=value.x, rare=value.y, relabund=value.x.x,
         normalized=value.y.y) %>%
  inner_join(., rand_group_count, by=c("sample" = "Group")) %>%
  inner_join(., rand_group_count, by=c("name" = "Group")) %>%
  mutate(n_diff = abs(n.x-n.y)) %>%
  select(-n.x, -n.y)

# comparison %>%
#   ggplot(aes(x=norare, y=rare, color=n_diff)) +
#   geom_point(size=0.25, alpha=0.25) +
#   geom_smooth()

comparison %>%
  pivot_longer(cols=c("norare", "rare", "relabund", "normalized"),
               names_to="type", values_to="dist") %>%
  mutate(diff_cat = cut_width(n_diff, width=500, boundary=0),
         diff_cat = str_replace(diff_cat, "[\\[\\(](.*),.*", "\\1"),
         diff_cat = as.numeric(diff_cat)) %>%
  filter(type != "norare") %>%
  group_by(diff_cat, type) %>%
  summarize(mean=mean(dist), sd=sd(dist), n=n(), .groups="drop") %>%
  filter(n > 100) %>%
  pivot_longer(cols=c("mean", "sd", "n")) %>%
  mutate(name = factor(name, levels=c("n", "mean", "sd"))) %>%
  ggplot(aes(x=diff_cat,  y=value, group=type, color=type)) +
  geom_line() +
  facet_wrap(~name, nrow=3, scales="free_y") +
  coord_cartesian(ylim=c(0, NA)) +
  theme_classic()

comparison %>%
  pivot_longer(cols=c("norare", "rare", "relabund", "normalized"),
               names_to="type", values_to="dist") %>%
  ggplot(aes(x=dist, fill=type)) +
  geom_density(alpha=0.5, color=NA) +
  coord_cartesian(ylim=c(0, 50), xlim=c(0,0.25)) +
  theme_classic()
  


rand_df <- rand %>%
  pivot_wider(names_from="rand_name", values_from="n", values_fill=0) %>%
  as.data.frame()

rownames(rand_df) <- rand_df[,1]
rand_df <- rand_df[, -1]


run_rarefy_bray <- function(x){
  
  mean_sd <- avgdist(rand_df, dmethod="bray", sample=x) %>%
    as.matrix() %>%
    as_tibble(rownames="samples") %>%
    pivot_longer(cols=-samples) %>%
    summarize(mean=mean(value), sd=sd(value))
  
  n <- rand_group_count %>% filter(n>=x) %>% nrow()
  
  bind_cols(n_seqs=x, mean_sd, n=n)
}

library(furrr)
plan(multisession)
rarefy_bray_results <- future_map_dfr(seq(1000, 15000, by=1000),
                                      run_rarefy_bray, .progress=TRUE)

rarefy_bray_results %>%
  pivot_longer(cols=c("mean", "sd", "n")) %>%
  mutate(name = factor(name, levels=c("n", "mean", "sd"))) %>%
  ggplot(aes(x=n_seqs, y=value)) +
  geom_line() +
  coord_cartesian(ylim=c(0, NA)) +
  facet_wrap(~name, nrow=3, scales="free_y")


run_rarefy_bray2 <- function(x){
  
  mean_sd <- avgdist(rand_df, dmethod="bray", sample=x) %>%
    as.matrix() %>%
    as_tibble(rownames="samples") %>%
    pivot_longer(cols=-samples) %>%
    filter(name < samples) %>%
    summarize(mean=mean(value), sd=sd(value))
  
  n <- rand_group_count %>% filter(n>=x) %>% nrow()
  
  bind_cols(n_seqs=x, mean_sd, n=n)
}

rarefy_bray_results2 <- future_map_dfr(seq(1000, 15000, by=1000),
                                      run_rarefy_bray2, .progress=TRUE)

rarefy_bray_results2 %>%
  pivot_longer(cols=c("mean", "sd", "n")) %>%
  mutate(name = factor(name, levels=c("n", "mean", "sd"))) %>%
  ggplot(aes(x=n_seqs, y=value)) +
  geom_line() +
  coord_cartesian(ylim=c(0, NA)) +
  facet_wrap(~name, nrow=3, scales="free_y")


