library("BSgenome.Mmusculus.UCSC.mm10")
library('tidyverse')
library("plyranges")

# get threeUTRs
ah<-AnnotationHub()
qr<-query(ah, c("EnsDb", "GRCm38"))

ensDB<-qr[["AH89211"]]
seqlevelsStyle(ensDB)<-"UCSC"

threeUTRs<-threeUTRsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>%
  unlist() %>% mutate(transcriptID = names(.)) %>%
  reduce_ranges_directed(transcriptID = paste(transcriptID, collapse=";")) %>% # merges overlapping UTRs
  mutate(UTRID = paste(seqnames,":", start,"-",end, strand, sep=""))

threeUTRs<-threeUTRs %>% dplyr::filter(seqnames %in% selectSeqLevels)
seqlevels(threeUTRs)<-selectSeqLevels

MSI_peaks = read.table(file="2_find_proximal_binding/MSI1-with_input.regions.ucsc.bed",header=F,stringsAsFactors=F,sep="\t")
colnames(MSI_peaks) = c("Chr","Start","End","Scores","ScoreSum","Strand")
MSI_peaks <- MSI_peaks %>% mutate(peakName = paste(Chr, ":", Start, "-", End, Strand,sep =""))

# keep only the standard chromosomes
selectSeqLevels<-paste0("chr", c(1:19, "X", "Y"))

## convert to genomics ranges so that we can collapse the overlaps and intersect the UTRs
MSI_GR = makeGRangesFromDataFrame(MSI_peaks,seqnames.field="Chr", start.field="Start", end.field="End", strand.field="Strand", keep.extra.columns=T) %>%
  dplyr::filter(seqnames %in% selectSeqLevels) # %>%
# reduce_ranges_directed(ScoreSum = paste(ScoreSum, collapse=";")) # not needed, already done
seqlevels(MSI_GR)<-selectSeqLevels
MSI_threeGR <- MSI_GR %>% join_overlap_intersect_directed(threeUTRs)

genome(MSI_GR) = "mm10"

# expand ranges by 10 either direction
MSI_GR <- MSI_GR + 10
MSI_GR <- MSI_GR %>% mutate(sequence=getSeq(BSgenome.Mmusculus.UCSC.mm10, .))

table(str_count(as.character(MSI_GR$sequence), "[GA]T{1,3}AGT"))

genome(MSI_threeGR) = "mm10"

MSI_threeGR <- MSI_threeGR + 10
MSI_threeGR <- MSI_threeGR %>% mutate(sequence=getSeq(BSgenome.Mmusculus.UCSC.mm10, .))

mean(str_count(as.character(MSI_threeGR$sequence), "[GA]T{1,3}AGT") > 0)
table(str_count(as.character(MSI_threeGR$sequence), "[GA]T{1,3}AGT"))

