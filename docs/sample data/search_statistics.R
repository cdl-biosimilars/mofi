# assume we're in the sample data folder
library(tidyverse)
df <- read.csv("search_statistics.csv", comment.char = "#")
df <- as_tibble(df)
ggplot(df %>%
         mutate(mass = Exp..Mass / 1000) %>%
         filter(Measure != "Search space size"),
       aes(x = mass, y = Value, color = Measure)) +
  geom_point() +
  geom_smooth() +
  ggtitle("Search statistics") +
  xlab("mass (kDa)") +
  ylab("number") +
  scale_color_brewer(type = "qual", palette = 2, name = "") +
  theme(legend.position = "bottom")