library(tidyverse)
library(vegan)

days_wanted <- c(0:9, 141:150)

early_late_df <- read_tsv("data/mice.shared") %>%
  select(Group, starts_with("Otu")) %>%
  pivot_longer(-Group) %>%
  separate(Group, into=c("animal", "day"), sep="D",
           remove=FALSE, convert=TRUE) %>%
  filter(day %in% days_wanted) %>%
  group_by(Group) %>%
  mutate(N = sum(value)) %>%
  ungroup() %>%
  filter(N >= 1828 & animal == "F3") %>%
  select(-N, -animal, -day) %>%
  pivot_wider(names_from="name", values_from="value", values_fill=0) %>%
  column_to_rownames("Group")

bray <- avgdist(early_late_df, dmethod="bray", sample=1828) %>%
  as.matrix() %>%
  as_tibble(rownames = "A") %>%
  pivot_longer(-A, names_to="B", values_to="distances")

bray %>%
  separate(A, sep="D", into=c("animal_a", "day_a"), convert=TRUE) %>%
  separate(B, sep="D", into=c("animal_b", "day_b"), convert=TRUE) %>%
  mutate(day_a = fct_reorder(as.character(day_a), day_a),
         day_b = fct_reorder(as.character(day_b), -day_b)) %>%
  select(-starts_with("animal")) %>%
  ggplot(aes(x=day_a, y=day_b, fill=distances)) +
  geom_tile()


jaccard <- avgdist(early_late_df, dmethod="jaccard", sample=1828) %>%
  as.matrix() %>%
  as_tibble(rownames = "A") %>%
  pivot_longer(-A, names_to="B", values_to="distances")

jaccard %>%
  separate(A, sep="D", into=c("animal_a", "day_a"), convert=TRUE) %>%
  separate(B, sep="D", into=c("animal_b", "day_b"), convert=TRUE) %>%
  mutate(day_a = fct_reorder(as.character(day_a), day_a),
         day_b = fct_reorder(as.character(day_b), -day_b)) %>%
  select(-starts_with("animal")) %>%
  ggplot(aes(x=day_a, y=day_b, fill=distances)) +
  geom_tile()

labels <- tibble(
  
  x=c(5, 16),
  y=c(3, 16),
  label=c("Bray-Curtis", "Jaccard")
)

inner_join(bray, jaccard, by=c("A", "B")) %>%
  separate(A, sep="D", into=c("animal_a", "day_a"), convert=TRUE) %>%
  separate(B, sep="D", into=c("animal_b", "day_b"), convert=TRUE) %>%
  select(day_a, day_b, bray=distances.x, jaccard=distances.y) %>%
  mutate(distances = if_else(day_a<day_b, bray, jaccard)) %>%
  mutate(day_a = fct_reorder(as.character(day_a), day_a),
         day_b = fct_reorder(as.character(day_b), -day_b)) %>%
  ggplot(aes(x=day_a, y=day_b, fill=distances)) +
  geom_tile() +
  geom_text(data=labels, aes(x=x, y=y, label=label), inherit.aes=FALSE,
            size=10) +
  scale_fill_gradient(low="#FF0000", high="#FFFFFF", name=NULL) +
  labs(x="Days post weaning", y="Days post weaning") +
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size=8),
        axis.text.y = element_text(hjust= 0.5))

ggsave("bray_jaccard_heatmap.png", width=6, height=5)
