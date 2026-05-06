library(Biostrings)
library(tidyverse)

clip_data <- read_tsv("2_find_proximal_binding/diff_chimeric_clip.zip") %>%
  dplyr::filter(Feature == "3'-UTR")

# split into regulated and unregulated
clip_data_grouped <- clip_data %>%
  dplyr::filter(!(is.na(padj))) %>%
  # label whether the interaction increased, decreased, had no effect, or had a mixed effect
  mutate(regulation = case_when(
    (padj < 0.05) & (log2FoldChange < 0) ~ "downregulated",
    (padj < 0.05) & (log2FoldChange > 0) ~ "upregulated",
    (padj > 0.90) & (baseMean > 120) ~ "notregulated",
    TRUE ~ "undetermined" # for all other cases
  ))

table(clip_data_grouped$regulation)

clipAGOGR <- makeGRangesFromDataFrame(clip_data_grouped, seqnames.field="seqnames", start.field='start', end.field='end', strand.field='strand', keep.extra.columns=T) %>%
  mutate(peak = gsub(" ", "", peak))

seeds <- readRDS("0_2_stoilov_microRNA_seed_mapping_with_genomic_locations/mapped_miRNA_seeds.rds")$UTR3_peaks

# there are duplicate miR and seedName
as_tibble(seeds@elementMetadata@listData) %>% dplyr::summarize(n=dplyr::n(), .by = c(miR, seedName)) %>% arrange(desc(n))

# as_tibble(seeds@elementMetadata@listData) %>% filter(miR == "mmu-miR-124-3p", seedName == "chr18:37841491-37841498+")
# ohh the geneID is different. i don't care about that. remove the duplicate rows https://stackoverflow.com/a/46656891
seedsDedup <- seeds %>% as_tibble %>% .[!duplicated(.[ , c("miR","seedName")]),] %>% makeGRangesFromDataFrame(keep.extra.columns=T)

UTR3_seeds_down  = join_overlap_inner_directed(seedsDedup, clipAGOGR) %>%
  as.data.frame() %>% 
  filter(regulation == "downregulated")

UTR3_seeds_down %>% dplyr::summarize(n=dplyr::n(), .by = c(miR.x, seedName)) %>% arrange(desc(n))

seqs <- UTR3_seeds_down$seedSeq
names(seqs) <- paste(UTR3_seeds_down$seedName, UTR3_seeds_down$miR.x, sep="|")

DNAStringSet(seqs) |> writeXStringSet('5_seed_motifs/seedSeqDownregulated.fasta')

###############################3
# put into STREME, somehow the 3-6 and 3-7 got completely different motifs, so i want to look at both
FIMO36 <- read_tsv('5_seed_motifs/3-6 FIMO matches.tsv') %>% drop_na() %>% makeGRangesFromDataFrame(seqnames = "sequence_name", keep.extra.columns=T)
FIMO37 <- read_tsv('5_seed_motifs/3-7 FIMO matches.tsv') %>% drop_na() %>% makeGRangesFromDataFrame(seqnames = "sequence_name", keep.extra.columns=T)

FIMO36f <- FIMO36 %>% filter(`p-value` < 0.05)
FIMO37f <- FIMO37 %>% filter(`p-value` < 0.05)

left36 <- FIMO36f %>% join_overlap_left_directed(FIMO37f)
left37 <- FIMO37f %>% join_overlap_left_directed(FIMO36f)
intersect <- FIMO36f %>% join_overlap_intersect_directed(FIMO37f)

# if you look at low p-value 3-7s, they're all identical to 3-6s essentially other than one letter. so just 3-7 is fine, and that motif is similar for both
left37 %>% as_tibble() %>% View()

UTR3_seeds_with_matches <- UTR3_seeds_down %>% makeGRangesFromDataFrame(keep.extra.columns=T) %>% join_overlap_inner_directed(FIMO37f)



##########################
# motifs found match mmu-mir-124 plus degeneracy where mir-182, 183, 96 are similar sequence except for a few at the beginning
# analysis works, but uhhh not exactly interesting!

# perhaps expand to 
# most UAGs are within 20 (average)-50 nucleotides (75%)
# maybe mask seed itself
# ehh mostly similar to 1_peak analysis anyways