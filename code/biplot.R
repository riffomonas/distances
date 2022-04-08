library(tidyverse)
library(vegan)
library(broom)
library(ggrepel)

set.seed(19760620)

# get shared file - days 0:9 and 141:150
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

# generate nmds ordination
set.seed(2)
nmds <- metaMDS(mice_dist)

nmds_positions <- scores(nmds) %>%
  as_tibble(rownames="Group")
  

nmds_positions %>%
  mutate(days = as.numeric(str_replace(Group, ".*D", "")),
         early = days < 10) %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=early)) +
  geom_point()
  

# generate subsampled shared file
subsample_shared_tbl <- shared_tbl %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total = sum(value)) %>%
  filter(total > 1828) %>%
  uncount(value) %>%
  slice_sample(n=1828) %>%
  count(Group, name) %>%
  pivot_wider(names_from="name", values_from="n", values_fill=0) %>%
  pivot_longer(-Group)

# correlate OTUs with x and y axis positions
nmds_positions
subsample_shared_tbl

nmds_shared <- inner_join(subsample_shared_tbl, nmds_positions)

cor_x <- nmds_shared %>%
  nest(data = -name) %>%
  mutate(cor_x = map(data,
                     ~cor.test(.x$value, .x$NMDS1,
                               method="spearman",
                               exact=FALSE) %>% tidy())) %>%
  unnest(cor_x) %>%
  select(name, estimate, p.value)
  
cor_y <- nmds_shared %>%
  nest(data = -name) %>%
  mutate(cor_y = map(data,
                     ~cor.test(.x$value, .x$NMDS2,
                               method="spearman",
                               exact=FALSE) %>% tidy())) %>%
  unnest(cor_y) %>%
  select(name, estimate, p.value)

correlations <- inner_join(cor_x, cor_y, by="name")

correlations %>%
  filter(p.value.x < 0.001 | p.value.y < 0.001)

high_corr <- correlations %>%
  filter(abs(estimate.x) > 0.75 | abs(estimate.y > 0.75)) %>%
  mutate(name = str_replace(name, "tu0+", "TU"))

# plot segments from (0, 0) to (Rx, Ry)
high_corr %>%
  ggplot(aes(x=0, xend=estimate.x, y=0, yend=estimate.y)) + 
  geom_segment()

nmds_positions %>%
  mutate(days = as.numeric(str_replace(Group, ".*D", "")),
         early = days < 10) %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=early)) +
  geom_point() +
  geom_segment(data=high_corr,
               aes(x=0, xend=estimate.x, y=0, yend=estimate.y), 
               inherit.aes=FALSE)


# plot text at (Rx, Ry)

nmds_positions %>%
  mutate(days = as.numeric(str_replace(Group, ".*D", "")),
         early = days < 10) %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=early)) +
  geom_point() +
  # geom_segment(data=high_corr,
  #              aes(x=0, xend=estimate.x, y=0, yend=estimate.y), 
  #              inherit.aes=FALSE) +
  geom_text_repel(data=high_corr,
            aes(x=estimate.x, y=estimate.y, label=name), 
            min.segment.length = 0.01, xlim=c(-1.2, 0.6),
            inherit.aes=FALSE) +
  geom_point(data=high_corr,
             color="black",
             aes(x=estimate.x, y=estimate.y), 
             inherit.aes = FALSE) +
  scale_color_manual(name=NULL,
                     breaks=c(FALSE, TRUE),
                     values=c("red", "blue"),
                     labels=c("Late", "Early")) +
  coord_cartesian(xlim=c(-1.1, 0.6)) +
  theme_classic()

ggsave("biplot.png", width=5, height=3)
