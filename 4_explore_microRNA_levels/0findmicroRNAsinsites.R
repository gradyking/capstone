library(tidyverse)
library(plyranges)

# find which microRNAs are more common in negatively regulated sites

plotRanges <- function(x, xlim=x, main=deparse(substitute(x)), col="black", sep=0.5, ...)
{
  height <- 1
  if (is(xlim, "IntegerRanges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
  title(main)
  axis(1)
}

clip_data <- read_tsv("2_find_proximal_binding/diff_chimeric_clip.zip") %>%
  filter(Feature == "3'-UTR") #filter to only sequences of interest

clip_data <- clip_data %>%
  dplyr::filter(!(is.na(padj))) %>%
  # label whether the interaction increased, decreased, had no effect, or had a mixed effect
  mutate(
    regulation = case_when(
      (padj < 0.05) & (log2FoldChange < 0) ~ "downregulated",
      (padj < 0.05) & (log2FoldChange > 0) ~ "upregulated",
      (padj > 0.90) & (baseMean > 120) ~ "notregulated",
      TRUE ~ "undetermined" # for all other cases
    ),
    peakName = gsub(" ", "", peak)
  )

miRsByRegulation <- clip_data %>%
  separate_longer_delim(miR, delim = ";") %>%
  summarize(n = n(), .by = c(regulation, miR)) %>%
  pivot_wider(names_from = regulation, values_from = n, values_fill = 0) %>%
  # get rid of infinities or zeroes
  filter(downregulated>0, notregulated>0) %>%
  mutate(
    # ratio, standardized to total # of miRs in the given regulation 
    ratio = downregulated / notregulated,
    rel_downvsnotratio = (downregulated / sum(downregulated)) / (notregulated / sum(notregulated)),
    log2rel_ratio = log2((downregulated / sum(downregulated)) / (notregulated / sum(notregulated))),
    total_n = downregulated + upregulated + undetermined + notregulated
  ) %>%
  arrange(desc(rel_downvsnotratio)) %>%
  mutate(miR = factor(miR, levels = miR))

miRsByRegulation %>% summarize(downSum = sum(downregulated),
                               upSum = sum(upregulated),
                               undeterSum = sum(undetermined),
                               notSum = sum(notregulated))

# highest 20 rel_downvsnotratio
# microRNAs that are much more prevalent in downregulated sites than nonregulated sites
miRsByRegulation %>% 
  filter(total_n > 20) %>%
  tail(20) %>%
  ggplot(aes(x = miR, y = rel_downvsnotratio)) + geom_col() +
  coord_flip() + theme_minimal()

# lowest 20 rel_ratio
miRsByRegulation %>% 
  head(20) %>%
  ggplot(aes(x = miR, y = rel_downvsnotratio)) + geom_col() +
  coord_flip() + theme_minimal()