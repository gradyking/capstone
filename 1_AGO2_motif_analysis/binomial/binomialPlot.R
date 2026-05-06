library(tidyverse)

data.frame(x = seq(0,100,1), density = dbinom(seq(0,100,1), 171, 18639 / (18639 + 74774))) %>% 
  ggplot(aes(x = x, y = density)) + 
  geom_col() + 
  labs(title = "Binomial(171, 0.1995)", 
       x = "# occurences of motif", 
       y= "Density") +
  theme_minimal()
ggsave("1_AGO2_motif_analysis/binomial/binomialPlot.png",width = 1350, height = 800, units = "px")

qbinom(0.95, 171, 18639 / (18639 + 74774))
