library("AnnotationHub")
library("tidyverse")
library("plyranges")

# essentially merging 2/0findProximalBindingExploring and 0_2/AGO-MSI-CLIP-comaping.R

ah<-AnnotationHub()
qr<-query(ah, c("EnsDb", "GRCm38"))

ensDB<-qr[["AH89211"]]
seqlevelsStyle(ensDB)<-"UCSC"

threeUTRs<-threeUTRsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>%
  unlist() %>% mutate(transcriptID = names(.)) %>%
  reduce_ranges_directed(transcriptID = paste(transcriptID, collapse=";")) %>% # merges overlapping UTRs
  mutate(UTRID = paste(seqnames,":", start,"-",end, strand, sep=""))

# standard chromosomes only
selectSeqLevels<-paste0("chr", c(1:19, "X", "Y"))
threeUTRs<-threeUTRs %>% dplyr::filter(seqnames %in% selectSeqLevels)
seqlevels(threeUTRs)<-selectSeqLevels

#######################3
## read the AGO2 differential (wild-type vs musashi-1 knockout)
clip_data <- read_tsv("2_find_proximal_binding/diff_chimeric_clip.zip") %>%
  dplyr::filter(Feature == "3'-UTR")

# ggplot(clip_data, aes(padj, log2FoldChange)) + geom_point(alpha = 0.05)
# ggplot(clip_data, aes(padj, baseMean)) + geom_point(alpha = 0.05)

# split into regulated and unregulated
clip_data_grouped <- clip_data %>%
  dplyr::filter(!(is.na(padj))) %>%
  # label whether the interaction increased, decreased, had no effect, or had a mixed effect
  # TODO: these are arbitrary cut-off values, must ask for revised numbers
  mutate(regulation = case_when(
    (padj < 0.05) & (log2FoldChange < 0) ~ "downregulated",
    (padj < 0.05) & (log2FoldChange > 0) ~ "upregulated",
    (padj > 0.90) & (baseMean > 120) ~ "notregulated",
    TRUE ~ "undetermined" # for all other cases
  ))

table(clip_data_grouped$regulation)

clipAGOGR <- makeGRangesFromDataFrame(clip_data_grouped, seqnames.field="seqnames", start.field='start', end.field='end', strand.field='strand', keep.extra.columns=T) %>%
  mutate(peak = gsub(" ", "", peak))

## recover from saved RDS
seedsList = readRDS("0_2_stoilov_microRNA_seed_mapping_with_genomic_locations/mapped_miRNA_seeds.rds")

UTR3_seeds  = join_overlap_inner_directed(seedsList[["UTR3_peaks"]] , clipAGOGR) %>%
  #mcols(.) %>%
  as.data.frame() %>% 
  rename(seed.start = "start", seed.end = "end", seed.width = "width")

#########################
# read the Msi-1 cross-linking regions

MSI1_sites = read.table("0_2_stoilov_microRNA_seed_mapping_with_genomic_locations/MSI1-with_input.sites.ucsc.bed", header = F)
colnames(MSI1_sites ) = c("seqnames","start","end","state","score","strand")
MSI1_sites = makeGRangesFromDataFrame(MSI1_sites, keep.extra.columns = T)

## annotate the x-link sites with the UTR data
UTR3_MSI = join_overlap_inner_directed(MSI1_sites,threeUTRs) %>% 
  as.data.frame()  %>%
  rename(MSI.start = "start", MSI.end = "end", MSI.width = "width")

## combine the miRNA seeds with the MSI1 x-link sites
CLIP_x_table = inner_join(UTR3_seeds, UTR3_MSI, by = c("segmentName" = "UTRID", "strand" = "strand", "seqnames" = "seqnames"), relationship = "many-to-many") %>%
  dplyr::mutate(distance = floor(MSI.start - (seed.start+seed.end)/2), absolute.distance = abs(distance))

# this is stoilov's regulation annotation
# instead of filtered "notregulated" with a baseMean > 120, he does abs(log2FoldChange)<0.1 which would include a few that don't have a lot of reads
CLIP_x_table = CLIP_x_table %>% dplyr::mutate(Sregulation = ifelse(padj<0.05 & log2FoldChange<0 , "downregulated", "undetermined"),
                                              Sregulation = ifelse(padj<0.05 & log2FoldChange>0 , "upregulated", Sregulation),
                                              Sregulation = ifelse(padj>0.95 & abs(log2FoldChange)<0.1 , "notregulated", Sregulation),
                                              distance.bin = cut_interval(distance, length=10))

tmp = CLIP_x_table %>% dplyr::group_by(peakName, regulation) %>% dplyr::summarise(absolute.distance = min(absolute.distance)) 

summaryCLIP_x_table = left_join(tmp, CLIP_x_table, by = c("peakName" = "peakName", "regulation" = "regulation", "absolute.distance" = "absolute.distance"))
table(summaryCLIP_x_table$regulation)
table(summaryCLIP_x_table$Sregulation)

summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") ) %>% dplyr::group_by(regulation) %>% 
  # dplyr::slice_sample(n=207) %>%
  ggplot( aes(x = regulation, y = distance)) + geom_violin() + labs(x = "AGO2 Regulation", y="Distance Between Closest AGO2/Msi-1 Pair", title="notregulated = (padj > 0.90) & (baseMean > 120)") + theme_minimal()

summaryCLIP_x_table %>% dplyr::filter(Sregulation %in% c("downregulated","notregulated") ) %>% dplyr::group_by(Sregulation) %>% 
  # dplyr::slice_sample(n=207) %>%
  ggplot( aes(x = Sregulation, y = distance)) + geom_violin() + labs(x = "AGO2 Regulation", y="Distance Between Closest AGO2/Msi-1 Pair", title="notregulated = (padj>0.95) & abs(log2FoldChange)<0.1") + theme_minimal()

# looks the same with his regulation definition
stoilovOriginal = read_tsv("0_2_stoilov_microRNA_seed_mapping_with_genomic_locations/stoilovFinalComparisonTable.tsv")
stoilovOriginal %>% dplyr::filter(regulation %in% c("downregulated","notregulated") ) %>% dplyr::group_by(regulation) %>% 
  # dplyr::slice_sample(n=207) %>%
  ggplot( aes(x = regulation, y = distance)) + geom_violin() + labs(x = "AGO2 Regulation", y="Distance Between Closest AGO2/Msi-1 Pair") + theme_minimal()
                            
CLIP_x_table %>% dplyr::filter(regulation == "notregulated") %>% ggplot(aes(x=baseMean)) + geom_histogram() + xlim(0,500) + labs(title="notregulated = (padj > 0.90) & (baseMean > 120)")
CLIP_x_table %>% dplyr::filter(Sregulation == "notregulated") %>% ggplot(aes(x=baseMean)) + geom_histogram() + xlim(0,500) + labs(title="notregulated = (padj>0.95) & abs(log2FoldChange)<0.1")
