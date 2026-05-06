library(tidyverse)
data.frame(x = seq(0,100,1), density = dbinom(seq(0,100,1), 190, 20710/83058)) %>% 
  ggplot(aes(x = x, y = density)) + geom_col() + labs(title = "binom(190, 0.249)", x = "successes")

data.frame(x = seq(0, 1000, 1), prob = pbinom(seq(0, 1000, 1), 190, 20710/83058)) %>% ggplot(aes(x = x, y = prob)) + geom_col()

18639 / (18639 + 74774)

data.frame(x = seq(0,100,1), density = dbinom(seq(0,100,1), 171, 18639 / (18639 + 74774))) %>% 
  ggplot(aes(x = x, y = density)) + geom_col() + labs(title = "binom(171, 0.1995)", x = "successes")

data.frame(x = seq(0,100,1), density = pbinom(seq(0,100,1), 171, 18639 / 74774)) %>% 
  ggplot(aes(x = x, y = density)) + geom_col() + labs(title = "binom(171, 0.249)", x = "successes")

1 - pbinom(87, 171, 18639 / 74774)
