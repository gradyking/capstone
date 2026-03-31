library("dplyr")
library("tidyr")
library("stringr")
library("plyranges")
library("AnnotationHub")
library(readr)

library("BSgenome.Mmusculus.UCSC.mm10")
library("Biostrings")
library("BiocParallel")

# BiocManager::install(listOfLibraries)

setwd("2_find_proximal_binding")

## multiporcessing parameters
param <- MulticoreParam(workers = 30, progressbar = TRUE)

## gather anotation and other external data we need
ah<-AnnotationHub()
qr<-query(ah, c("EnsDb", "GRCm38"))

ensDB<-qr[["AH89211"]]
seqlevelsStyle(ensDB)<-"UCSC"

# pull protein coding genes
proteinGenes<-genes(ensDB, columns=c("gene_id","symbol","description"), filter= ~ gene_biotype == "protein_coding")

# keep only the standard chromosomes
selectSeqLevels<-seqlevels(proteinGenes)[str_detect(seqlevels(proteinGenes),regex("^chr.*[0-9XY]"))]
proteinGenes<-proteinGenes %>% dplyr::filter(seqnames %in% selectSeqLevels)
seqlevels(proteinGenes)<-selectSeqLevels

temp <- threeUTRsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>%
  unlist()
sort(table(names(temp)), decreasing=TRUE)[1:20]

# taken from https://www.bioconductor.org/packages/devel/bioc/vignettes/IRanges/inst/doc/IRangesOverview.pdf
# idk why this function isn't just built-in to IRanges
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
# plotRanges(IRanges(c(5,6), c(10,11))) # testing the function

# no overlap on any of these??
plotRanges((temp[names(temp) == "ENSMUST00000000001"]@ranges))
plotRanges((temp[names(temp) == "ENSMUST00000178006"]@ranges))
plotRanges((temp[names(temp) == "ENSMUST00000204457"]@ranges))
plotRanges((temp[names(temp) == "ENSMUST00000179598"]@ranges))

# oh boy some of the UTRs overlap with different UTRs of other genes
sort(countOverlaps(temp@ranges, drop.self=TRUE), decreasing=TRUE)[1:10]

# worst example is 20488 which overlaps with 59 UTRs from 58 other genes
lookup <- findOverlaps(temp@ranges, drop.self=TRUE)
i <- which(names(temp) == "ENSMUST00000020488")
j <- lookup[lookup@from == i]@to
plotRanges((temp[c(i,j)]@ranges), main = "ENSMUST00000020488 overlaps with 59 others")

# but this is fine, right? i'm just filtering the .bams based on whether its within any UTR, doesn't really matter the particular gene

threeUTRs<-threeUTRsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>%
  unlist() %>% mutate(transcriptID = names(.)) %>%
  reduce_ranges_directed(transcriptID = paste(transcriptID, collapse=";")) %>% # merges overlapping UTRs
  mutate(UTRID = paste(seqnames,":", start,"-",end, strand, sep=""))
  
threeUTRs<-threeUTRs %>% dplyr::filter(seqnames %in% selectSeqLevels)
seqlevels(threeUTRs)<-selectSeqLevels

###########################################
## read the AGO2 differential (wild-type vs musashi-1 knockout)

clip_data <- read_tsv("diff_chimeric_clip.zip") %>%
  dplyr::filter(Feature == "3'-UTR")

# split into regulated and unregulated
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

clipAGOGR <- makeGRangesFromDataFrame(clip_data_grouped, seqnames.field="seqnames", start.field='start', end.field='end', strand.field='strand', keep.extra.columns=T)

# negativeregulatedAGO2 <- clip_data_grouped %>% dplyr::filter(group == "negative")
# positiveregulatedAGO2 <- clip_data_grouped %>% dplyr::filter(group == "posi")
# nonregulatedAGO2 <- clip_data_grouped %>% dplyr::filter(group == "control")

#######################################################################

## read the musashi-1 eCLIP 
MSI_peaks = read.table(file="MSI1-with_input.regions.ucsc.bed",header=F,stringsAsFactors=F,sep="\t")
colnames(MSI_peaks) = c("Chr","Start","End","Scores","ScoreSum","Strand")
MSI_peaks <- MSI_peaks %>% tidyr::unite(GeneID, c(Chr,Start,End,Strand), remove=F) 

## convert to genomics ranges so that we can collapse the overlaps and intersect the UTRs

MSI_GR = makeGRangesFromDataFrame(MSI_peaks,seqnames.field="Chr", start.field="Start", end.field="End", strand.field="Strand", keep.extra.columns=T) %>%
  dplyr::filter(seqnames %in% selectSeqLevels) # %>%
  # reduce_ranges_directed(ScoreSum = paste(ScoreSum, collapse=";")) # not needed, already done
seqlevels(MSI_GR)<-selectSeqLevels

# MSI_GR <- MSI_GR %>%
#   mutate(sequence = getSeq(BSgenome.Mmusculus.UCSC.mm10, .),
#          names=paste(seqnames,":",start,"-",end, sep=""))

## not sure if this was only necessary for the other .bed file
# ## filter by number of replicates detecting the peak and the maximum peak score
# ## in the dat the median score is 1
# filteredGR = MSI_GR %>%
#   dplyr::group_by(miR) %>%
#   plyranges::reduce_ranges_directed(Score  = max(Score)) %>%
#   plyranges::ungroup() %>%
#   plyranges::filter(nRep == 3, Score >= 1 ) # Filter top 50% - score >1

##############################################
# filter binding sites to only ones in 3'-UTRs
MSI_GR = MSI_GR %>% join_overlap_intersect_directed(threeUTRs)

MSI_GR_table = as.tibble(MSI_GR) %>% 
  mutate(MSI_center = as.integer((end +start)/2)) %>%
  dplyr::select(UTRID, MSI_center)

Ago_GR_Table = clipAGOGR %>% join_overlap_intersect_directed(threeUTRs) %>%
  as.tibble() %>%
  mutate(AGO_center = as.integer((end +start)/2)) %>%
  dplyr::select(UTRID, AGO_center, group)

Msi_Ago = MSI_GR_table %>% inner_join(Ago_GR_Table, by='UTRID', relationship = "many-to-many")

compared <- Msi_Ago %>% mutate(dist = AGO_center-MSI_center) %>%
  group_by(UTRID, group) %>%
  dplyr::filter(abs(dist) == min(abs(dist)))

compared %>% ggplot(aes(group, dist)) + geom_point()
compared %>% dplyr::filter(abs(dist) < 5000, group %in% c('negative', 'positive')) %>% ggplot(aes(dist)) + geom_histogram() + facet_wrap(~ group)

compared %>% dplyr::filter(abs(dist) < 5000, group %in% c('negative', 'control')) %>% ungroup() %>% summarise(p_val = t.test(dist ~ group)$p.value)

distneg <- compared %>% ungroup() %>% dplyr::filter(group == "negative") %>% dplyr::select(dist) %>% pull()
distcontrol <- compared %>% ungroup() %>% dplyr::filter(group == "control") %>% dplyr::select(dist) %>% pull()

# https://mgimond.github.io/Stats-in-R/F_test.html
var.test(distneg, distcontrol) # neg is lower variance, but less samples

anova = compared %>% ungroup() %>% 
  aov(dist ~ group,.)

summary(anova)
TukeyHSD(anova)
# find distances between nearest regulated AGO2 sites and musashi-1 sites
# find distances between nearest unregulated AGO2 sites and musashi-1 sites

# find if distributions are random
# compare distributions 
# potentially expand to nearest n sites

# more complex analysis for farther interactions? even possible?