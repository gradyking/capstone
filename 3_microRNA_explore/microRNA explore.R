library(tidyverse)
library(plyranges)

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

clip_data <- read_tsv("3_microRNA_explore/diff_chimeric_clip.zip") %>%
  filter(Feature == "3'-UTR") #filter to only sequences of interest

clip_data_grouped <- clip_data %>%
  dplyr::filter(!(is.na(padj))) %>%
  # label whether the interaction increased, decreased, had no effect, or had a mixed effect
  # TODO: these are arbitrary cut-off values, must ask for revised numbers
  mutate(group = case_when(
    (padj < 0.05) & (log2FoldChange < 0) ~ "negative",
    (padj < 0.05) & (log2FoldChange > 0) ~ "positive",
    (padj > 0.90) & (baseMean > 120) ~ "control",
    TRUE ~ "middle" # for all other cases
  ))

clip_data_grouped %>% group_by(group) %>% reframe(split(miR, ";"))


###############################
seeds <- readRDS("0_stoilov_microRNA_seed_mapping_with_genomic_locations/mapped_miRNA_seeds.rds")

# almost all sequences are 8 long
table(sapply(seeds$UTR3_peaks$seedSeq, nchar))

# some overlap between regions
plotRanges((seeds$UTR3_peaks %>% dplyr::filter(seqnames == 'chr1'))@ranges)

###########################
# finding UAG binding domains within musashi sites
