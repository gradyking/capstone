library("AnnotationHub")

library("BSgenome.Mmusculus.UCSC.mm10")
library("Biostrings")
library("BiocParallel")

library('tidyverse')
library("plyranges")

# BiocManager::install(listOfLibraries)

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

# temp <- threeUTRsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>%
#   unlist()
# 
# # temp2 <- threeUTRsByTranscript(ensDB, 'tx_biotype') %>% unlist()
# # sort(table(names(temp)), decreasing=TRUE)
# # sum(table(names(temp)) > 2)
# 
# # taken from https://www.bioconductor.org/packages/devel/bioc/vignettes/IRanges/inst/doc/IRangesOverview.pdf
# # idk why this function isn't just built-in to IRanges
# plotRanges <- function(x, xlim=x, main=deparse(substitute(x)), col="black", sep=0.5, ...)
# {
#   height <- 1
#   if (is(xlim, "IntegerRanges"))
#     xlim <- c(min(start(xlim)), max(end(xlim)))
#   bins <- disjointBins(IRanges(start(x), end(x) + 1))
#   plot.new()
#   plot.window(xlim, c(0, max(bins)*(height + sep)))
#   ybottom <- bins * (sep + height) - height
#   rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
#   title(main)
#   axis(1)
# }
# # plotRanges(IRanges(c(5,6), c(10,11))) # testing the function
# 
# # no overlap on any of these??
# plotRanges((temp[names(temp) == "ENSMUST00000000001"]@ranges))
# plotRanges((temp[names(temp) == "ENSMUST00000178006"]@ranges)) # many exons after the end of the transcript -> would be degraded by non-sense mediated decay almost immediately
# plotRanges((temp[names(temp) == "ENSMUST00000204457"]@ranges))
# plotRanges((temp[names(temp) == "ENSMUST00000179598"]@ranges))
# 
# # oh boy some of the UTRs overlap with different UTRs of other genes
# sort(countOverlaps(temp@ranges, drop.self=TRUE), decreasing=TRUE)[1:10]
# 
# # worst example is 20488 which overlaps with 59 UTRs from 58 other genes
# lookup <- findOverlaps(temp@ranges, drop.self=TRUE)
# i <- which(names(temp) == "ENSMUST00000020488")
# j <- lookup[lookup@from == i]@to
# plotRanges((temp[c(i,j)]@ranges), main = "ENSMUST00000020488 overlaps with 59 others")

# but this is fine, right? i'm just filtering the .bams based on whether its within any UTR, doesn't really matter the particular gene

threeUTRs<-threeUTRsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>%
  unlist() %>% mutate(transcriptID = names(.)) %>%
  reduce_ranges_directed(transcriptID = paste(transcriptID, collapse=";")) %>% # merges overlapping UTRs
  mutate(UTRID = paste(seqnames,":", start,"-",end, strand, sep=""))
  
threeUTRs<-threeUTRs %>% dplyr::filter(seqnames %in% selectSeqLevels)
seqlevels(threeUTRs)<-selectSeqLevels

###########################################
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
  mutate(group = case_when(
    (padj < 0.05) & (log2FoldChange < 0) ~ "negative",
    (padj < 0.05) & (log2FoldChange > 0) ~ "positive",
    (padj > 0.90) & (baseMean > 120) ~ "control",
    TRUE ~ "middle" # for all other cases
  ))

table(clip_data_grouped$group)

clipAGOGR <- makeGRangesFromDataFrame(clip_data_grouped, seqnames.field="seqnames", start.field='start', end.field='end', strand.field='strand', keep.extra.columns=T) %>%
  mutate(peak = gsub(" ", "", peak))

# negativeregulatedAGO2 <- clip_data_grouped %>% dplyr::filter(group == "negative")
# positiveregulatedAGO2 <- clip_data_grouped %>% dplyr::filter(group == "posi")
# nonregulatedAGO2 <- clip_data_grouped %>% dplyr::filter(group == "control")

#######################################################################

## read the musashi-1 eCLIP, made by PureCLIP
MSI_peaks = read.table(file="2_find_proximal_binding/MSI1-with_input.regions.ucsc.bed",header=F,stringsAsFactors=F,sep="\t")
colnames(MSI_peaks) = c("Chr","Start","End","Scores","ScoreSum","Strand")
MSI_peaks <- MSI_peaks %>% mutate(peakName = paste(Chr, ":", Start, "-", End, Strand,sep =""))

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
# MSI_GR = MSI_GR %>% join_overlap_intersect_directed(threeUTRs)

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


compared %>% dplyr::filter(abs(dist) < 50, group %in% c('negative', 'control')) %>% ggplot(aes(dist)) + geom_histogram(binwidth=1) + facet_wrap(~ group)

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

############################################################
############################################################
############################################################
# the above analysis has bad resolution-- the peaks are much wider than the true binding sites
# stoilov did a computational method to find "seeds" in the AGO2 peaks for the associated microRNAs
# i need to find the UAG binding sites within the musashi peaks

# mm10 and GRCm38 are the same?
genome(MSI_GR) = "mm10"
MSI_GR <- MSI_GR %>% mutate(sequence=getSeq(BSgenome.Mmusculus.UCSC.mm10, .))

# find UAG binding sites ((G/A)U_n AGU, n=1 to 3) https://www.tandfonline.com/doi/full/10.1128/MCB.21.12.3888-3900.2001
# str_locate_all(as.character(MSI_GR$sequence), "TAG")
mean(str_count(as.character(MSI_GR$sequence), "TAGT") > 0)
mean(str_count(as.character(MSI_GR$sequence), "TAG") > 0)
# hmmm, only 28% of the binding sites contain UAG

# wait, so this sequence is from the original DNA, right? the mRNA would be the reverse complement, right?
# 5' UAG 3'
# 3' AUC 5' -> 5' CTA 3' (DNA)
mean(str_count(as.character(MSI_GR$sequence), "CTA") > 0)
# still only 16% of binding sites have UAG.
# well i can simply continue and fix it later if it's wrong

# some have more than one match. for those, i will expand the search to be the above (G/A) U_n AGU, which would be
# 5' (G/A) U_n AGU 3'
# 3' (C/U) A_n UCA 5' -> 5' ACT A_n C/T 3' DNA
# which can be expressed in regex as ACTA{1,3}[CT]

# about 3.1% have more than one match to the UAG motif
mean(str_count(as.character(MSI_GR$sequence), "CTA") > 1)
table(str_count(as.character(MSI_GR$sequence), "CTA"))

# only 45 have more than one match for the full musashi motif
mean(str_count(as.character(MSI_GR$sequence), "ACTA{1,3}[CT]") > 1)
table(str_count(as.character(MSI_GR$sequence), "ACTA{1,3}[CT]"))

# get UAG motif locations
MSI_GR$matches <- str_locate_all(as.character(MSI_GR$sequence), "TAG")

MSI_motifGR <- MSI_GR %>% 
  as_tibble() %>% 
  rowwise() %>% 
  reframe(matches = split(matches, row(matches)),
          seqnames = seqnames,
          start = start,
          end = end,
          width = width,
          strand = strand) %>%
  as_granges()
# # to remove the duplicates, I would ideally to get the site nearest the high point of the peak on the graph
# # but i don't think that's accessible unfortunately. i'll take the one nearest the center for now
# 
# # find which have more than 1 match
# idxmoreThanOne <- which(sapply(MSI_GR$matches, dim)[1,] > 1)
# idxequalsOne <- which(sapply(MSI_GR$matches, dim)[1,] == 1)
# 
# 
# # extract from MSI_GR to avoid queries & writing to large dataframe
# matches_list <- MSI_GR$matches
# 
# # stupid fix to get rid of nested matrices inside lists for those with only one match
# matches_list[idxequalsOne] <- lapply(idxequalsOne, function(x) matches_list[[x]][1,])
# widths <- width(MSI_GR)
# for(i in idxmoreThanOne){
#   if(i %% 100 == 0) cat(i, "\n")
#   
#   # find which match is closest to the peak center
#   match <- matches_list[[i]]
#   
#   # find centers of motifs (generalized in case motif length changes)
#   centers <- floor((match[,2] + match[,1]) / 2)
#   
#   # find center of peak, find which is closest
#   peakCenter <-  floor((widths[i] + 1) / 2)
#   idxMin <- which.min(abs(centers-peakCenter))
#   
#   # rewrite match to be the best match (list is because for some reason the matches are in lists)
#   matches_list[i] <- list(match[idxMin,])
# }
# 
# # put new de-duplicated matches (as lists) back into column
# MSI_GR$matches <- matches_list
# 
# # no more duplicates!
# table(sapply(MSI_GR$matches, length)/2)
# 
# # ok i will filter to those with a motif
# MSI_motifGR <- MSI_GR %>% 
#   dplyr::filter(sapply(.$matches, length) > 0)

# now matches has the location relative to the beginning of the string, where the first index is 1
# but i need a GRange relative to the original genome
# so, if a sequence starts at 100, and the UAG motif is at 3 to 5, 
# 100, 101, 102, 103, 104
# 1,   2,   3,   4,   5
# the gene location would be 100 + 3 - 1 = 102 to 100 + 5 - 1 = 104

starts <- start(MSI_motifGR)
start(MSI_motifGR) <- starts + sapply(MSI_motifGR$matches, "[[", 1) - 1
end(MSI_motifGR) <- starts + sapply(MSI_motifGR$matches, "[[", 2) - 1

# now the GRanges contain the correct absolute location on the genome!

###################
# import stoilov's AGO2 seeds
threeUTRSeeds <- readRDS("0_stoilov_microRNA_seed_mapping_with_genomic_locations/mapped_miRNA_seeds.rds")$UTR3_peaks

# find groups from chimeric data
clip_data <- read_tsv("2_find_proximal_binding/diff_chimeric_clip.zip")

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

clipAGOGR <- makeGRangesFromDataFrame(clip_data_grouped, seqnames.field="seqnames", start.field='start', end.field='end', strand.field='strand', keep.extra.columns=T) %>%
  mutate(peak = gsub(" ", "", peak))

# obtain the group for each peak from the clipAGOGR dataframe
AGOGRWithGroup <- threeUTRSeeds %>% 
  # as.tibble() %>%
  # mutate(AGO_center = floor((end + start)/2)) %>%
  join_overlap_inner_directed(clipAGOGR) %>%
  as.tibble() %>%
  mutate(AGO_center = floor((end + start)/2))

mean(is.na(AGOGRWithGroup$group)) # cool, no NAs

MSI_motifGR <- MSI_motifGR %>% join_overlap_inner_directed(threeUTRs) %>% as.tibble() %>% mutate(MSI_center = floor((end + start)/2))

# join on the UTRID, which is called segmentName in the RDS
Msi_Agoprecise <- inner_join(MSI_motifGR, AGOGRWithGroup, 
                             by=join_by("UTRID" == "segmentName"), 
                             relationship = "many-to-many")

compared <- Msi_Agoprecise %>% mutate(dist = AGO_center-MSI_center) %>%
  group_by(seedName, group) %>% #MISTAKE: grouping by UTRID instead of peakName or seedName, since multiple AGO2 sites can be in one UTR
  dplyr::filter(abs(dist) == min(abs(dist)))

compared %>% group_by(group) %>% summarize(n=dplyr::n())

compared %>% dplyr::filter(group %in% c('negative', 'control')) %>% ggplot(aes(group, dist)) + geom_violin()

compared %>% dplyr::filter(abs(dist) < 100, group %in% c('negative', 'control')) %>% ggplot(aes(dist)) + geom_histogram(binwidth=1) + facet_wrap(~ group)
ggsave("highPreciseCloseHist.png")

compared %>% dplyr::filter(abs(dist) < 1000, group %in% c('negative', 'control')) %>% ggplot(aes(dist)) + geom_histogram(binwidth=5) + facet_wrap(~ group, scales = "free_y")
ggsave("highPreciseFarHist.png")

compared %>% dplyr::filter(abs(dist) < 3000, group %in% c('negative', 'control')) %>% ungroup() %>% summarise(p_val = t.test(dist ~ group)$p.value)
