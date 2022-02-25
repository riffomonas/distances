library(tidyverse)
library(vegan)
library(SRS)
library(furrr)

plan(multisession)

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


run_iteration <- function(x){
  rand <- shared %>%
    uncount(value) %>%
    mutate(rand_name = sample(name)) %>%
    select(-name) %>%
    count(Group, rand_name)
  
  rand_group_count <- rand %>%
    group_by(Group) %>%
    summarize(n = sum(n))
  
  smallest_group <- min(rand_group_count$n)
  
  rabund_df <- rand %>%
    group_by(Group) %>%
    mutate(rel_abund = n / sum(n)) %>%
    ungroup() %>%
    select(Group, rand_name, rel_abund) %>%
    pivot_wider(names_from="rand_name", values_from="rel_abund", values_fill=0) %>%
    as.data.frame()
  
  rownames(rabund_df) <- rabund_df[, 1]
  rabund_df <- rabund_df[, -1]
  
  norare_df <- rand %>%
    pivot_wider(names_from="rand_name", values_from="n", values_fill=0) %>%
    as.data.frame()
  
  rownames(norare_df) <- norare_df[, 1]
  norare_df <- norare_df[, -1]
  
  normal_df <- norare_df %>%
    t() %>%
    as.data.frame() %>%
    SRS(., Cmin=smallest_group) %>%
    t() %>%
    as.data.frame()
  
  norare_dist <- vegdist(norare_df, method="bray")
  rabund_dist <- vegdist(rabund_df, method="bray")
  normal_dist <- vegdist(normal_df, method="bray")
  rare_dist <- avgdist(norare_df, dmethod="bray", sample=smallest_group)
  
  treatment <- if_else(rand_group_count$n < median(rand_group_count$n), "A", "B")
  
# run this line for when treatment is not confounded w/ sequencing depth
# treatment <- sample(treatment)
  
  norare_test <- adonis(norare_dist ~ treatment)$aov.tab$`Pr(>F)`[1]
  norare_sig <- norare_test < 0.05
  
  rare_test <- adonis(rare_dist ~ treatment)$aov.tab$`Pr(>F)`[1]
  rare_sig <- rare_test < 0.05
  
  rabund_test <- adonis(rabund_dist ~ treatment)$aov.tab$`Pr(>F)`[1]
  rabund_sig <- rabund_test < 0.05
  
  normal_test <- adonis(normal_dist ~ treatment)$aov.tab$`Pr(>F)`[1]
  normal_sig <- normal_test < 0.05
  
  c(norare=norare_sig, rare=rare_sig, rabund=rabund_sig, normal=normal_sig)
}

sig_results <- future_map_dfr(1:100, run_iteration, .id="seed", .progress=TRUE)

sig_results %>%
  pivot_longer(-seed, names_to="method", values_to="significant") %>%
  group_by(method) %>%
  summarize(false_pos = mean(significant))
