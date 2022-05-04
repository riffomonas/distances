library(tidyverse)
library(vegan)

phylodiv <- read_tsv("data/mice.phylodiv.rarefaction") %>%
  filter(numSampled == 1804) %>%
  pivot_longer(-numSampled, names_to="Group", values_to="phylodiv") %>%
  select(-numSampled)

otu_data <- read_tsv("data/mice.shared") %>%
  select(Group, starts_with("Otu")) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(n = sum(value)) %>%
  ungroup() %>%
  filter(n >= 1804) %>%
  select(-n) %>%
  pivot_wider(names_from="name", values_from="value") %>%
  column_to_rownames("Group")

richness <- otu_data %>%
  rarefy(sample=1804) %>%
  as_tibble(rownames="Group") %>%
  select(Group, richness=value)

shannon_iteration <- function(){

    otu_data %>%
    rrarefy(sample=1804) %>%
    diversity()

}

shannon <- replicate(100, shannon_iteration()) %>%
  as_tibble(rownames="Group", .name_repair = "unique") %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  summarize(shannon = mean(value))

combined <- inner_join(phylodiv, richness, by="Group") %>%
  inner_join(., shannon, by="Group")

combined %>%
  filter(shannon > 2) %>%
  ggplot(aes(x=richness, y=phylodiv)) +
  geom_point() +
  geom_smooth(se=FALSE, size=2, method="lm") +
  labs(x="Richness", y="Phylogenetic Diversity") +
  theme_classic()

ggsave("richness_phylodiv.png", width=6, height=5)

combined %>%
  filter(shannon > 2) %>%
  ggplot(aes(x=shannon, y=phylodiv)) +
  geom_point() +
  geom_smooth()

combined %>%
  filter(shannon > 2) %>%
  ggplot(aes(x=shannon, y=richness)) +
  geom_point() +
  geom_smooth()

cor.test(combined$phylodiv, combined$richness, method="spearman", exact=FALSE)
cor.test(combined$phylodiv, combined$shannon, method="spearman", exact=FALSE)
cor.test(combined$shannon, combined$richness, method="spearman", exact=FALSE)
