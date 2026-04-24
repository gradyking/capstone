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

# ok, compare these to ones that have seeds? see if some of them have the motifs?
###############################
# find which microRNAs are more common in negatively regulated seed sites

# extract 3'UTR seeds
seeds <- readRDS("0_2_stoilov_microRNA_seed_mapping_with_genomic_locations/mapped_miRNA_seeds.rds")$UTR3_peaks

# find which seeds are negatively regulated
length(intersect(seeds$peakName, clip_data$peakName))
length(setdiff(seeds$peakName, clip_data$peakName))

# what.... the flip.
sum(seeds$peakName %in% clip_data$peakName)
sum(clip_data$peakName %in% seeds$peakName)

# almost all sequences are 8 long
table(sapply(seeds$seedSeq, nchar))

# many of the top ratios don't have a seed, i guess because there aren't enough examples of binding
miRsByRegulation %>% 
  mutate(hasSeed = miR %in% intersect(miRsByRegulation$miR,seeds$miR)) %>% 
  map_df(rev) %>% # reverse row order https://stackoverflow.com/a/52398627
  head(25)

# make col chart only for miRs with valid seeds
miRsByRegulation %>% 
  mutate(hasSeed = miR %in% intersect(miRsByRegulation$miR,seeds$miR)) %>% 
  filter(hasSeed) %>%
  tail(20) %>%
  ggplot(aes(x = miR, y = rel_downvsnotratio)) + geom_col() +
  coord_flip() + theme_minimal() +
  labs(y = "relative ratio of binding, downregulated / not regulated")

# wait... what does this even mean? 
# i'm trying to answer the question: 
# each seed has a matching microRNA, and I am trying to find which of these microRNAs 
seedsWithRatio <- seeds %>% as_tibble() %>%
  left_join(miRsByRegulation, by = "miR") %>%
  arrange(desc(rel_downvsnotratio))

seedsWithRatio %>% select(seqnames, start, end, strand, miR, seedSeq, downregulated, notregulated, ratio, rel_downvsnotratio) %>% head(20)

table(seedsWithRatio %>% select(miR, seedSeq))
