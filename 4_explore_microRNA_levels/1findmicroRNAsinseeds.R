library(tidyverse)
library(plyranges)

# find which microRNAs are more common in negatively regulated seed sites

# A) find negatively regulated seed sites

# extract 3'UTR seeds, remove duplicate gene info
seeds <- readRDS("0_2_stoilov_microRNA_seed_mapping_with_genomic_locations/mapped_miRNA_seeds.rds")$UTR3_peaks
seedsNoDuplicate <- seeds %>% as_tibble %>% .[!duplicated(.[ , c("peakName", "miR","seedName")]),] %>% makeGRangesFromDataFrame(keep.extra.columns=T)

# count number of miRs in seeds table
length(table(seedsNoDuplicate$miR))
# 156 unique miRs

# almost all sequences are 8 long
table(sapply(seeds$seedSeq, nchar))

# bring in clip data to find negatively regulation
clip_data <- read_tsv("2_find_proximal_binding/diff_chimeric_clip.zip") %>%
  dplyr::filter(Feature == "3'-UTR")

# count number of miRs in clip table
clip_data %>%
  separate_longer_delim(miR, delim = ";") %>%
  summarize(n = n(), .by = c(miR)) %>% dim()
# 728 unique miRs

# compare the sets. all of the seed miRs are in the clip table
length(setdiff(seedsNoDuplicate$miR, (clip_data %>% separate_longer_delim(miR, delim = ";"))$miR))
length(setdiff((clip_data %>% separate_longer_delim(miR, delim = ";"))$miR, seedsNoDuplicate$miR))

# split into regulated and unregulated
clip_data_grouped <- clip_data %>%
  dplyr::filter(!(is.na(padj))) %>%
  # label whether the interaction increased, decreased, had no effect, or had a mixed effect
  # TODO: these are arbitrary cut-off values, must ask for revised numbers
  mutate(
    regulation = case_when(
      (padj < 0.05) & (log2FoldChange < 0) ~ "downregulated",
      (padj < 0.05) & (log2FoldChange > 0) ~ "upregulated",
      (padj > 0.90) & (baseMean > 120) ~ "notregulated",
      TRUE ~ "undetermined" # for all other cases
    ),
    peakName = gsub(" ", "", peak)
  )

clipAGOGR <- makeGRangesFromDataFrame(clip_data_grouped, keep.extra.columns=T)

# investigating peakName to see if I can join using that
length(intersect(seeds$peakName, clip_data_grouped$peakName))
length(setdiff(seeds$peakName, clip_data_grouped$peakName))
sum(seeds$peakName %in% clip_data_grouped$peakName)
sum(clip_data_grouped$peakName %in% seeds$peakName)
# what.... the flip.
# ok these come from different sources. the clip_data comes from the chimeric algorithm that excludes peaks that aren't in both msi-1 knockout and wild type
# seeds has microRNA seeds for all AGO2 binding sites, not just those that come from the chimeric algorithm
# so one will not be a subset of the other
# so instead of using peakName, just do an inner directed join

# about 1300 seeds don't have a corresponding entry in AGO chimeric clip, which is fine
join_overlap_left_directed(seedsNoDuplicate, clipAGOGR) %>% summarize(na = sum(is.na(peakName.y)))

UTR3_seeds  = join_overlap_inner_directed(seedsNoDuplicate, clipAGOGR) %>%
  as.data.frame()

# lost two microRNAs of the 156 for those that happened to not be in the same 
setdiff(intersect(seedsNoDuplicate$miR, (clip_data %>% separate_longer_delim(miR, delim = ";"))$miR), UTR3_seeds$miR.x)

# B) find most common microRNAs in seeds by regulation

miRsinSeedsbyRegulation <- UTR3_seeds %>% 
  summarize(n = n(), .by = c(regulation, miR.x)) %>%
  pivot_wider(names_from = regulation, values_from = n, values_fill = 0) %>%
  mutate(
    ratio = downregulated / notregulated,
    rel_downvsnotratio = (downregulated / sum(downregulated)) / (notregulated / sum(notregulated)),
    n = downregulated + notregulated
  ) %>%
  arrange(desc(ratio))%>%
  mutate(miR.x = factor(miR.x, levels = miR.x))

miRsinSeedsbyRegulation %>%
  filter(n > 10) %>%
  ggplot(aes(x = miR.x, y = rel_downvsnotratio)) + geom_col() +
  coord_flip() + theme_minimal() +
  labs(y = "relative ratio of binding, downregulated / not regulated", x="") + 
  geom_text(
    aes(y = 2.5,label = paste0("n=", n)),
    size = 3, hjust=0
  ) +
  expand_limits(y = 2.6)
ggsave("4_explore_microRNA_levels/differentBindingOfMicroRNAs.png", width = 6, height = 6, unit = "in")

# is this significant? do permutation. take all seeds, make two random groups (match in size)
# repeat 100000x, what ratio is more/less than you'd expect to get? make distribution of ratios, which are too far out in the tails?


# the miR with high ratio tend to have UAG while the others do not
# proving this?
# trimer differential

library("scanMiRData")

data(mmu, package = "scanMiRData")
mmu[["mmu-miR-16-5p"]]$mirseq
