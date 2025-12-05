# trying to extract gene letters from the ranges
library(tidyverse)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(plyranges)
library(ggdensity)

clip_data <- read_tsv("diff_chimeric_clip.zip") %>%
  filter(Feature == "3'-UTR") #filter to only sequences of interest

# data was pre-analyzed like in this paper https://www.biorxiv.org/content/10.1101/2022.02.13.480296v1.full
# "symbol","miR","Score","description","gene_id" came from there

# first 5 are coordinates in the genome
# baseMean: average normalized count values taken over all samples (raw is number of alignments, normalize)
# log2foldchange: log2(treatment count / control count). ex: log2(50/10) = 2.32
# lfc_SE: standard error (used for standardized statistic)
# stat- difference in deviance in likelihood ratio test
# pvalue- individual p-values
# padj- family-wise adjusted p-values benjamini-hochberg -> control for multiple testing & remove low reads/outliers
# symbol: name gene
# MiR: micro-RNAs on peak (microRNA connected to mRNA in chimeric step)
# Score: more hits for that microRNA at that position (normalized) (reads per million chimeras)
# description: gene description
# gene_id: gene id
# Feature: type of mRNA section
# peak: unique name from coordinates


#how many NAs are in each column
apply(is.na(clip_data), 2, sum)

# visualize groups for changes in expression
ggplot(clip_data, aes(padj, log2FoldChange)) + geom_hdr()
ggplot(clip_data, aes(padj, log2FoldChange)) + geom_point(alpha = 0.05)
ggplot(clip_data, aes(padj)) + geom_histogram()
ggplot(clip_data, aes(log2FoldChange)) + geom_histogram()

clip_data_grouped <- clip_data %>%
  filter(!(is.na(padj))) %>%
  # label whether the interaction increased, decreased, had no effect, or had a mixed effect
  # TODO: these are arbitrary cut-off values, must ask for revised numbers
  mutate(group = case_when(
    (padj < 0.05) & (log2FoldChange < 0) ~ "negative",
    (padj < 0.05) & (log2FoldChange > 0) ~ "positive",
    (padj > 0.90) & (baseMean > 120) ~ "control",
    TRUE ~ "middle" # for all other cases
  ))
table(clip_data_grouped$group)

clip_data_filtered <- clip_data_grouped %>% filter(group == "negative" | group == "control")
clip_granges <- clip_data_filtered %>% 
  as_granges()


# takes into account seqnames, start, end and strand
clip_seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, clip_granges)
clip_granges$sequence <- clip_seqs
names(clip_granges$sequence) <- clip_granges$peak
clip_granges <- clip_granges %>% filter(lapply(clip_granges$sequence, length) > 10)

# now, need to compare the two groups for "words" common within themselves but not common between each other
# (differences between the two groups)

(clip_granges %>% filter(group == 'negative'))$sequence |> writeXStringSet('inputSeq.fasta')
(clip_granges %>% filter(group == 'control'))$sequence |> writeXStringSet('controlSeq.fasta')

# throwing into https://meme-suite.org/meme/tools/streme
# https://meme-suite.org/meme/info/status?service=STREME&id=appSTREME_5.5.817622894229752125478339