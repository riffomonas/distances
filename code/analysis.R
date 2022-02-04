library(tidyverse)
library(glue)
source("code/read_matrix.R")

dist_matrix <- read_matrix("data/mice.braycurtis.dist")

dist_tbl <- as_tibble(dist_matrix, rownames="samples")

sample_lookup <- dist_tbl %>% 
  select(samples) %>%
  mutate(delimited = str_replace(samples,
                                 "^(([FM])\\d+)D(\\d+)$",
                                 "\\2-\\1-\\3")) %>%
  separate(col=delimited,
           into=c("sex", "animal", "day"), sep="-",
           convert=TRUE)

days_wanted <- c(0:9, 141:150)

dist_matrix <- dist_tbl %>%
  pivot_longer(cols=-samples, names_to="b", values_to="distances") %>%
  inner_join(., sample_lookup, by="samples") %>%
  inner_join(., sample_lookup, by=c("b" = "samples")) %>%
  filter(day.x %in% days_wanted & day.y %in% days_wanted) %>%
  select(samples, b, distances) %>%
  pivot_wider(names_from="b", values_from="distances") %>%
  select(-samples) %>%
  as.dist()

pcoa <- cmdscale(dist_matrix, eig=TRUE, add=TRUE)

positions <- pcoa$points
colnames(positions) <- c("pcoa1", "pcoa2")

percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)

pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

labels <- c(glue("PCo Axis 1 ({pretty_pe[1]}%)"),
            glue("PCo Axis 2 ({pretty_pe[2]}%)"))

positions %>%
  as_tibble(rownames = "samples") %>%
  ggplot(aes(x=pcoa1, y=pcoa2)) +
  geom_point() +
  labs(x=labels[1], y=labels[2])

tibble(pe = cumsum(percent_explained),
       axis = 1:length(percent_explained)) %>%
  ggplot(aes(x=axis, y=pe)) +
  geom_line() +
  coord_cartesian(xlim = c(1, 10), ylim=c(0, 50)) +
  scale_x_continuous(breaks=1:10)
