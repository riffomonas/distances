library(tidyverse)
library(vegan)

days_wanted <- c(0:9, 141:150)

early_late_tbl <- read_tsv("data/mice.shared") %>%
  select(Group, starts_with("Otu")) %>%
  pivot_longer(-Group) %>%
  separate(Group, into=c("animal", "day"), sep="D", remove=FALSE, convert=TRUE) %>%
  filter(day %in% days_wanted) %>%
  group_by(Group) %>%
  mutate(N = sum(value)) %>%
  ungroup() %>%
  filter(N >= 1828) %>%
  select(-N) %>%
  pivot_wider(names_from="name", values_from="value", values_fill=0)

early_late_metadata <- early_late_tbl %>% 
  select(Group, animal, day) %>%
  mutate(period = if_else(day < 10, "early", "late"),
         sex = if_else(str_detect(animal, "F"), "female", "male"),
         period_sex=paste0(period, sex))

early_late_dist <- early_late_tbl %>%
  select(-animal, -day) %>%
  column_to_rownames("Group") %>%
  avgdist(sample=1828)

set.seed(1)
early_late_nmds <- metaMDS(early_late_dist) %>%
  scores() %>%
  as_tibble(rownames="Group")

metadata_nmds <- inner_join(early_late_metadata, early_late_nmds)

centroid <- metadata_nmds %>%
  group_by(period) %>%
  summarize(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

metadata_nmds %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=period)) +
  geom_point() +
  stat_ellipse(show.legend=FALSE) +
  geom_point(data=centroid, size=5, shape=21, color="black",
             aes(fill=period), show.legend=FALSE)


centroid <- metadata_nmds %>%
  group_by(sex) %>%
  summarize(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

metadata_nmds %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=sex)) +
  geom_point() +
  stat_ellipse(show.legend=FALSE) +
  geom_point(data=centroid, size=5, shape=21, color="black",
             aes(fill=sex), show.legend=FALSE)


test <- adonis(as.dist(early_late_dist)~early_late_metadata$period, permutations = 1e3)
p_value <- test$aov.tab$`Pr(>F)`[1]

adonis(early_late_dist~early_late_metadata$sex)

adonis(early_late_dist~early_late_metadata$period*early_late_metadata$sex)

adonis(early_late_dist~early_late_metadata$period+early_late_metadata$sex)

adonis(early_late_dist~early_late_metadata$period*early_late_metadata$sex,
       strata=early_late_metadata$animal)

adonis(early_late_dist~early_late_metadata$period,
       strata=early_late_metadata$animal)

bd <- betadisper(early_late_dist, early_late_metadata$period)
anova(bd)
permutest(bd)

bd <- betadisper(early_late_dist, early_late_metadata$sex)
anova(bd)
permutest(bd)

bd <- betadisper(early_late_dist, early_late_metadata$period_sex)
anova(bd)
permutest(bd, pairwise = TRUE)
