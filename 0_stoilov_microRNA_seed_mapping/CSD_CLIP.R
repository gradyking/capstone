library("dplyr")
library("tidyr")
library("stringr")
library("plyranges") 
library("AnnotationHub")
library("ggplot2")
library("openxlsx")
library("mirbase.db")
library("ggforce")

library("scanMiR")
library("scanMiRData")
library("BSgenome.Mmusculus.UCSC.mm10")
library("Biostrings")
library("BiocParallel")
library(miRBaseConverter)

## Wraps the findSeedMatches function
findSeed = function(sequence, seed, name) {
    # Call findSeedMatches on a single pair
    # Note: Ensure arguments match findSeedMatches signature, 
    # you may need to convert sequences to DNAStringSet inside the function.
    noMatchDF$transcript = name
    sequence = DNAStringSet(sequence) %>% setNames(name)
    res = scanMiR::findSeedMatches(seqs = sequence, seeds = mmu[seed], ret = "aggregate", verbose = F)
    if(!is.data.frame(res)){ res = noMatchDF}
    return(list(res))
  }

## multiporcessing parameters
param <- MulticoreParam(workers = 30, progressbar = TRUE)

setwd("/projects/Retina/CLIP/AGO2/difClip/CDS_peaks")

## Prep a list for the results
peaksList = vector(mode = "list", length = 11)
names(peaksList) = c("UTR3_peaks", #peaks in the 3' UTR - genome mapped data
                     "CDS_peaks",  #peaks in the CDS - genome mapped data
                     "UTR5_peaks",  #peaks in the 5' UTR - genome mapped data
                     "intron_peaks",  #peaks in the introns - genome mapped data
                     "intergenic_peaks", # interegnic peaks (not introns, CDS, UTRs) - genome mapped data
                     "UTR3_peaksTX", #peaks in the 3' UTR - mapped from genome to transcript data
                     "CDS_peaksTX",  #peaks in the CDS - mapped from genome to transcript data
                     "UTR5_peaksTX"  #peaks in the 5' UTR - mapped from genome to transcript data
#                      "UTR3_peaksGN", #peaks in the 3' UTR - mapped from genome to transcript data and then back to genome to remove duplicates
#                      "CDS_peaksGN",  #peaks in the CDS - mapped from genome to transcript data and then back to genome to remove duplicates
#                      "UTR5_peaksGN"  #peaks in the 5' UTR - mapped from genome to transcript data and then back to genome to remove duplicates
                     )

## gather anotation and other external data we need
data(mmu, package = "scanMiRData")
ah<-AnnotationHub()
qr<-query(ah, c("EnsDb", "GRCm38"))

ensDB<-qr[["AH89211"]]
seqlevelsStyle(ensDB)<-"UCSC"


proteinGenes<-genes(ensDB, columns=c("gene_id","symbol","description","entrezid"), filter= ~ gene_biotype == "protein_coding") %>% 
          plyranges::filter(!is.na(entrezid), !is.na(gene_id))
selectSeqLevels<-seqlevels(proteinGenes)[str_detect(seqlevels(proteinGenes),regex("^chr.*[0-9XY]"))]

allGenes = genes(ensDB, columns=c("gene_id","symbol","description","entrezid"))

proteinGenes<-proteinGenes %>% plyranges::filter(seqnames %in% selectSeqLevels)
seqlevels(proteinGenes)<-selectSeqLevels

allGenes<-allGenes %>% plyranges::filter(seqnames %in% selectSeqLevels)
seqlevels(allGenes)<-selectSeqLevels

allGenes = allGenes %>% plyranges::mutate(start = start - 1000, end = end + 1000) # extend to capture antisense promoter transcripts and extended UTRs
intergenic = gaps(allGenes, ignore.strand = T) %>%
            reduce_ranges_directed()

CDS = cdsBy(ensDB, by = "tx",filter= ~ tx_biotype == "protein_coding") %>% 
            unlist()%>%
            reduce_ranges_directed()
            
CDS = CDS %>% plyranges::filter(seqnames %in% selectSeqLevels)
seqlevels(CDS)<-selectSeqLevels

threeUTRs<-threeUTRsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>% 
            unlist() %>%
            reduce_ranges_directed()

fiveUTRs<-fiveUTRsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>% 
            unlist() %>%
            reduce_ranges_directed()

threeUTRs<-threeUTRs %>% plyranges::filter(seqnames %in% selectSeqLevels)
seqlevels(threeUTRs)<-selectSeqLevels
#mcols(threeUTRs)$UTRlocation<-paste(seqnames(threeUTRs),":",start(threeUTRs),"-",end(threeUTRs),sep="")

fiveUTRs<-fiveUTRs %>% plyranges::filter(seqnames %in% selectSeqLevels)
seqlevels(fiveUTRs)<-selectSeqLevels
#mcols(fiveUTRs)$UTRlocation<-paste(seqnames(fiveUTRs),":",start(fiveUTRs),"-",end(fiveUTRs),sep="")

introns = intronsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>% 
            unlist()%>%
            reduce_ranges_directed()
            
introns = introns %>% plyranges::filter(seqnames %in% selectSeqLevels)
seqlevels(introns)<-selectSeqLevels


###########################################
## reads the bed files with chimeric reads
## Using CLIP data only from flox animals in which Msi1/2 have not been knocked out

bedFiles = list.files(path = "../../files", pattern= ".*_F.*filtered.chimeric_peaks.bed", recursive =T, full.names=T)

peaksBed = vector(mode = "list", length = length(bedFiles))
for (i in 1:length(bedFiles)){
  peaksBed[[i]] = read.table(file=bedFiles[i],header=F,stringsAsFactors=F,sep="\t")
  peaksBed[[i]][,"Replicate"] = i
}

peaksBed = bind_rows(peaksBed)
colnames(peaksBed)<-c("Chr","Start","End","miR","Score","Strand", "Replicate")



## convert to genomics ranges so that we can collapse the overlaps and intersect the UTRs
peaksGR<-makeGRangesFromDataFrame(peaksBed,seqnames.field="Chr", start.field="Start", end.field="End", strand.field="Strand", keep.extra.columns=T)    
peaksGR<-peaksGR %>% plyranges::filter(seqnames %in% selectSeqLevels)
seqlevels(peaksGR)<-selectSeqLevels

peaksGR

## filter by number of replicates detecting the peak and the maximum peak score 
## in the dat the median score is 1
filteredGR = peaksGR %>%
          plyranges::group_by(miR) %>%
          plyranges::reduce_ranges_directed( nRep = lengths(unique(Replicate)), Score  = max(Score)) %>% 
          plyranges::ungroup() %>%
          plyranges::filter(nRep == 3, Score >= 1 ) # Filter top 50% - score >1
              
filteredGR <- filteredGR %>%
  mutate(sequence = getSeq(BSgenome.Mmusculus.UCSC.mm10, .), 
         names=paste(seqnames,":",start,"-",end, sep=""))

         
#####################################
#### Find miR seeds in the CLIP peaks

## scanMir::findSeedMatches returns null when no matches are found and a table if matches are found
## We need a "noMatch" table to place instead of the NULLs so that bind_rows will work
noMatchDF = data.frame(c("seq1"), c(0), c(0), c(0), c(0), c(0))  
colnames(noMatchDF) = c("transcript", "repression", "8mer", "7mer", "6mer", "non-canonical")


# Find if there is a seed for th miRNA in each CLIP peak
seedMatches = bpmapply(findSeed, filteredGR$sequence, filteredGR$miR, filteredGR$names, BPPARAM = param)


seedMatchesTable = bind_rows(seedMatches)

#fix column names that start with numbers and contain "-" sign
colnames(seedMatchesTable) = c("transcript", "repression", "S.8mer", "S.7mer", "S.6mer", "S.non_canonical")

# check count seeds and summarize
seedMatchesTable = seedMatchesTable %>% 
  rowwise() %>%
  mutate(nSeeds = sum(c_across(c(S.8mer,S.7mer,S.6mer,S.non_canonical))), 
         hasSeed = (nSeeds>0),
         hasCanonicalSeed = ((nSeeds - S.non_canonical)>0)) 

# add the miRNA seed search results to the genomic ranges
mcols(filteredGR) = bind_cols(as.data.frame(mcols(filteredGR)),seedMatchesTable) 


# add gene annotation
filteredGRAnnotated = join_overlap_left_directed(filteredGR,proteinGenes)

# intersect with the different sequence segment types
peaksList[["UTR3_peaks"]] =  join_overlap_inner_directed(filteredGRAnnotated,threeUTRs)
peaksList[["UTR5_peaks"]] =  join_overlap_inner_directed(filteredGRAnnotated, fiveUTRs)
peaksList[["CDS_peaks"]] = join_overlap_inner_directed(filteredGRAnnotated,CDS)
peaksList[["CDS_peaks"]] = subsetByOverlaps(peaksList[["CDS_peaks"]], peaksList[["UTR3_peaks"]], invert = T) # some sequences can be both coding and UTRs depending on transcript

peaksList[["intron_peaks"]] = join_overlap_inner_directed(filteredGRAnnotated,introns)
peaksList[["intron_peaks"]] = subsetByOverlaps(peaksList[["intron_peaks"]], peaksList[["UTR3_peaks"]], invert = T) # some introns may contain exons
peaksList[["intron_peaks"]]= subsetByOverlaps(peaksList[["intron_peaks"]], peaksList[["UTR5_peaks"]], invert = T)
peaksList[["intron_peaks"]]= subsetByOverlaps(peaksList[["intron_peaks"]], peaksList[["CDS_peaks"]], invert = T)
peaksList[["intergenic_peaks"]] = join_overlap_inner(filteredGR,intergenic)

# peaksList[["UTR3_peaks"]] %>% plyranges::filter(nRep ==3) %>% plyranges::summarize(Seed = sum(hasSeed)/n(), Canonical = sum(hasCanonicalSeed)/n())
# 
# peaksList[["CDS_peaks"]] %>% plyranges::filter(nRep ==3) %>% plyranges::summarize(Seed = sum(hasSeed)/n(), Canonical = sum(hasCanonicalSeed)/n())
# 
# peaksList[["intron_peaks"]] %>% plyranges::filter(nRep ==3) %>% plyranges::summarize(Seed = sum(hasSeed)/n(), Canonical = sum(hasCanonicalSeed)/n())
# peaksList[["intergenic_peaks"]] %>% plyranges::filter(nRep ==3) %>% plyranges::summarize(Seed = sum(hasSeed)/n(), Canonical = sum(hasCanonicalSeed)/n())

##################################################################################################
## Do everything by transcript instead of genomic coordinates to account for peaks in the CDS split between exons
##
transcripts = exonsBy(ensDB, by= "tx",  filter= ~ tx_biotype == "protein_coding") 
transcripts = transcripts %>% unlist() %>% plyranges::filter(seqnames %in% selectSeqLevels) %>% split(., names(.))
txSeqs = extractTranscriptSeqs(BSgenome.Mmusculus.UCSC.mm10, transcripts)

peaksTX = mapToTranscripts(peaksGR, transcripts, ignore.strand = F) %>% trim()
mcols(peaksTX) = mcols(peaksGR)[mcols(peaksTX)$xHits, , drop=FALSE]
mcols(peaksTX)$length = seqlengths(peaksTX)[as.character(seqnames(peaksTX))]
 seqlengths(peaksTX) = NA

 
 #######################################################################################
 ### IGNORE
 ### extend the peaks by 25nt on each side to better capture seed positions. 
 ## This only results in less than 1% reduction of peak numbers after peaks are overlapped and reduced
# peaksTX = peaksTX %>%  
#           plyranges::mutate( name = as.character(seqnames(.)), start = pmax(1,start-25), end = pmin(length, end+25))
          
          
### filter by number of replicates detecting the peak and maximum score
selectedPeaksTX = peaksTX %>%
          plyranges::group_by(miR) %>%
          plyranges::reduce_ranges_directed( nRep = lengths(unique(Replicate)), Score  = max(Score)) %>% 
          plyranges::ungroup() %>%
          plyranges::filter(nRep == 3, Score >= 1 )

selectedPeaksTX = selectedPeaksTX %>% plyranges::mutate( name = as.character(seqnames(.)))
txSeqDF = data.frame(row.names = names(txSeqs), sequence = as.character(txSeqs))

## add transcript sequences and cut the sequences to the peak ranges
mcols(selectedPeaksTX)$sequence = txSeqDF[mcols(selectedPeaksTX)$name,"sequence"]
selectedPeaksTX = selectedPeaksTX %>% plyranges::mutate(sequence = substr(sequence,start, end))
 
 

# bptasks(param) = 1000
# noMatchDF = data.frame(c("seq1"), c(0), c(0), c(0), c(0), c(0))  
# colnames(noMatchDF) = c("transcript", "repression", "8mer", "7mer", "6mer", "non-canonical")

## find the seed matches 
seedMatches = bpmapply(findSeed, selectedPeaksTX$sequence, selectedPeaksTX$miR, selectedPeaksTX$name, BPPARAM = param)

seedMatchesTable = bind_rows(seedMatches)
colnames(seedMatchesTable) = c("transcript", "repression", "S.8mer", "S.7mer", "S.6mer", "S.non_canonical")
seedMatchesTable = seedMatchesTable %>% 
  rowwise() %>%
  mutate(nSeeds = sum(c_across(c(S.8mer,S.7mer,S.6mer,S.non_canonical))), 
         hasSeed = (nSeeds>0),
         hasCanonicalSeed = ((nSeeds - S.non_canonical)>0)) 

mcols(selectedPeaksTX) = bind_cols(as.data.frame(mcols(selectedPeaksTX)),seedMatchesTable) 

 

## add transcript annotation
txAnnotation = transcriptsBy(ensDB, by = "gene", columns=c("tx_id","gene_id","symbol","description"), filter = ~ tx_biotype == "protein_coding") %>% 
                unlist() %>% 
                mcols() %>%
                as.data.frame()

mcols(selectedPeaksTX) = left_join(as.data.frame(mcols(selectedPeaksTX)), txAnnotation, by = c("name" = "tx_id"))

selectedPeaksTX = selectedPeaksTX %>% plyranges::mutate(peakID = paste(seqnames, ":", start, "-", end, sep =""))

### intersect with the different sequence segment types
cdsTX = mapToTranscripts(CDS, transcripts, ignore.strand = F) %>% reduce()
threeUTRtx =  mapToTranscripts(threeUTRs, transcripts, ignore.strand = F) %>% reduce()
fiveUTRtx = mapToTranscripts(fiveUTRs, transcripts, ignore.strand = F) %>% reduce()

peaksList[["UTR3_peaksTX"]]  =  join_overlap_inner_directed(selectedPeaksTX,threeUTRtx)
peaksList[["UTR5_peaksTX"]]  =  join_overlap_inner_directed(selectedPeaksTX,fiveUTRtx)
peaksList[["CDS_peaksTX"]] = join_overlap_inner_directed(selectedPeaksTX,cdsTX)
peaksList[["CDS_peaksTX"]] = subsetByOverlaps(peaksList[["CDS_peaksTX"]], peaksList[["UTR3_peaksTX"]], invert = T) # some sequences can be both coding and UTRs depending on transcript


##
# convert to genomic coordinates to collapse duplicate peaks originating from 
# transcript isoforms of the same gene
# NOTE: mapping peaks that span exon-exon jucntions will result in coordinates that include 
# the whole genomic range, including the introns - one-to-one mapping that does not split the peaks into exons
# DO NOT USE THE GENOMIC COORDINATES FROM THIS OPERATION FOR SEQUENCE EXTRACTION


# peaksList[["UTR3_peaksGN"]] = mapFromTranscripts(peaksList[["UTR3_peaksTX"]], transcripts) %>% 
#                   reduce_ranges_directed(txMap= min(xHits))
# peaksList[["UTR5_peaksGN"]] = mapFromTranscripts(peaksList[["UTR5_peaksTX"]], transcripts) %>% 
#                   reduce_ranges_directed(txMap= min(xHits))
# peaksList[["CDS_peaksGN"]] = mapFromTranscripts(peaksList[["CDS_peaksTX"]], transcripts) %>% 
#                   reduce_ranges_directed(txMap= min(xHits))
# # add annotation, including peaks sequence, microRNa etc to the genomic range
# mcols(peaksList[["UTR3_peaksGN"]]) = mcols(peaksList[["UTR3_peaksTX"]])[peaksList[["UTR3_peaksGN"]]$txMap,]
# mcols(peaksList[["UTR5_peaksGN"]]) = mcols(peaksList[["UTR5_peaksTX"]])[peaksList[["UTR5_peaksGN"]]$txMap,]
# mcols(peaksList[["CDS_peaksGN"]]) = mcols(peaksList[["CDS_peaksTX"]])[peaksList[["CDS_peaksGN"]]$txMap,]

# 

peaksList[["UTR3_peaks"]] %>% mcols() %>% as.data.frame() %>% dplyr::select(sequence, miR, hasSeed, hasCanonicalSeed) %>% distinct()  %>% plyranges::summarise(nPeaks=dplyr::n(), Seed = sum(hasSeed)/dplyr::n(), Canonical = sum(hasCanonicalSeed)/dplyr::n())
peaksList[["UTR5_peaks"]] %>% mcols() %>% as.data.frame() %>% dplyr::select(sequence, miR, hasSeed, hasCanonicalSeed) %>% distinct()  %>% plyranges::summarise(nPeaks=dplyr::n(), Seed = sum(hasSeed)/dplyr::n(), Canonical = sum(hasCanonicalSeed)/dplyr::n())
peaksList[["CDS_peaks"]] %>%  mcols() %>% as.data.frame() %>% dplyr::select(sequence, miR, hasSeed, hasCanonicalSeed) %>% distinct()  %>% plyranges::summarise(nPeaks=dplyr::n(), Seed = sum(hasSeed)/dplyr::n(), Canonical = sum(hasCanonicalSeed)/dplyr::n())

peaksList[["intron_peaks"]] %>%  mcols() %>% as.data.frame() %>% dplyr::select(sequence, miR, hasSeed, hasCanonicalSeed) %>% distinct()  %>% plyranges::summarise(nPeaks=dplyr::n(), Seed = sum(hasSeed)/dplyr::n(), Canonical = sum(hasCanonicalSeed)/dplyr::n())
peaksList[["intergenic"]] %>%  mcols() %>% as.data.frame() %>% dplyr::select(sequence, miR, hasSeed, hasCanonicalSeed) %>% distinct()  %>% plyranges::summarise(nPeaks=dplyr::n(), Seed = sum(hasSeed)/dplyr::n(), Canonical = sum(hasCanonicalSeed)/dplyr::n())
 filteredGRAnnotated %>% mcols() %>% as.data.frame() %>% dplyr::select(sequence, miR, hasSeed, hasCanonicalSeed) %>% distinct() %>% dim()

## This produces result very close to the genome mapped peaks
peaksList[["UTR3_peaksTX"]] %>% mcols() %>% as.data.frame() %>% dplyr::select(sequence, miR, hasSeed, hasCanonicalSeed) %>% distinct()  %>% plyranges::summarise(nPeaks=dplyr::n(), Seed = sum(hasSeed)/dplyr::n(), Canonical = sum(hasCanonicalSeed)/dplyr::n())
peaksList[["UTR5_peaksTX"]] %>% mcols() %>% as.data.frame() %>% dplyr::select(sequence, miR, hasSeed, hasCanonicalSeed) %>% distinct()  %>% plyranges::summarise(nPeaks=dplyr::n(), Seed = sum(hasSeed)/dplyr::n(), Canonical = sum(hasCanonicalSeed)/dplyr::n())
peaksList[["CDS_peaksTX"]] %>%  mcols() %>% as.data.frame() %>% dplyr::select(sequence, miR, hasSeed, hasCanonicalSeed) %>% distinct()  %>% plyranges::summarise(nPeaks=dplyr::n(), Seed = sum(hasSeed)/dplyr::n(), Canonical = sum(hasCanonicalSeed)/dplyr::n())

## this reduces the peaks more than it should - need to investigate why
# peaksList[["UTR3_peaksGN"]] %>% mcols() %>% as.data.frame() %>% dplyr::select(sequence, miR, hasSeed, hasCanonicalSeed) %>% distinct()  %>% plyranges::summarise(nPeaks=dplyr::n(), Seed = sum(hasSeed)/dplyr::n(), Canonical = sum(hasCanonicalSeed)/dplyr::n())
# peaksList[["UTR5_peaksGN"]] %>% mcols() %>% as.data.frame() %>% dplyr::select(sequence, miR, hasSeed, hasCanonicalSeed) %>% distinct()  %>% plyranges::summarise(nPeaks=dplyr::n(), Seed = sum(hasSeed)/dplyr::n(), Canonical = sum(hasCanonicalSeed)/dplyr::n())
# peaksList[["CDS_peaksGN"]] %>%  mcols() %>% as.data.frame() %>% dplyr::select(sequence, miR, hasSeed, hasCanonicalSeed) %>% distinct()  %>% plyranges::summarise(nPeaks=dplyr::n(), Seed = sum(hasSeed)/dplyr::n(), Canonical = sum(hasCanonicalSeed)/dplyr::n())

results = vector(mode = "list")
for (n in names(peaksList)){
 if (!is.null(peaksList[[n]])) {
  if ("names" %in% names(mcols(peaksList[[n]]))){
    tmp = peaksList[[n]] %>% 
                          mcols() %>% 
                          as.data.frame() %>%
                          dplyr::group_by(names) %>%
                          summarise(hasSeed = max(hasSeed), hasCanonicalSeed = max(hasCanonicalSeed))
  }
  else{
    tmp = peaksList[[n]] %>% 
                          mcols() %>% 
                          as.data.frame() %>%
                          dplyr::group_by(peakID) %>%
                          summarise(hasSeed = max(hasSeed), hasCanonicalSeed = max(hasCanonicalSeed))
  }
  results[[n]] = tmp %>%
                  dplyr::ungroup() %>%
                  dplyr::summarise(`Total Peaks`=dplyr::n(), `Fraction With Matched Seed` = sum(hasSeed)/dplyr::n(), `Fraction With Canonical Seed` = sum(hasCanonicalSeed)/dplyr::n())
  results[[n]]$Region = n
 }
 # else {results[[n]] = data.frame(n = 0, Seed = 0, Canonical = 0, Region = n)}
}

results = bind_rows(results)
write.csv(file = "AGO2-CLIP_peaks_summary.csv", results, row.names = F)

miRCounts = vector(mode = "list")
for (p in c("UTR3_peaks","UTR5_peaks","intron_peaks","CDS_peaks")){
  miRCounts[[p]] = peaksList[[p]] %>% 
                    plyranges::filter(hasSeed) %>% 
                    plyranges::group_by(miR) %>% 
                    plyranges::summarise(n=n()) %>%
                    as.data.frame() %>%
                    dplyr::mutate(Region = p)#%>% as.data.frame() %>% dplyr::arrange(desc(n)) %>% head(10)
}

miRCounts = bind_rows(miRCounts) 

miRCounts = miRCounts %>% dplyr::group_by(Region) %>% dplyr::mutate(fraction = n/sum(n)) %>% ungroup()


mirAnnotation1 = getMiRNATable(species = "mmu") %>% dplyr::select(Mature1, Mature1_Seq, Mature1_Acc) %>% dplyr::mutate(Arm = "5p")
colnames(mirAnnotation1) = c("miR","Seq","Accession","Arm")
mirAnnotation2 = getMiRNATable(species = "mmu") %>% dplyr::select(Mature2, Mature2_Seq, Mature2_Acc) %>% dplyr::mutate(Arm = "3p")
colnames(mirAnnotation2) = c("miR","Seq","Accession","Arm")
mirAnnotation = bind_rows(mirAnnotation1,mirAnnotation2) %>% tidyr::drop_na() %>% distinct()

miRCounts = left_join(miRCounts, mirAnnotation, by = "miR")

miRCounts = miRCounts %>% mutate(Family = str_extract(miR, regex("(?<=mmu-)[A-z]+-[0-9]+")), Family = paste(Family, Arm, sep ='-'), Seed = str_sub(Seq, 2, 7))

seedFam = miRCounts %>% dplyr::select(Seed, Family, miR) %>% dplyr::group_by(Seed) %>% dplyr::summarize(SeedFam = paste(unique(Family), collapse = "/"), miRNAs = paste(unique(miR), collapse = "/"))

miRCounts = left_join(miRCounts, seedFam, by = "Seed")

miRCounts = miRCounts %>% dplyr::select(miRNAs, Family, SeedFam, Seed, Region, n, fraction) %>% dplyr::group_by(miRNAs, Family, SeedFam, Seed, Region) %>% summarise(Counts = sum(n), Fraction = sum(fraction))



miRCounts %>% dplyr::filter(Fraction > 0.02) %>% mutate(`Seed Family` = paste(Seed, Family)) %>%
  ggplot(aes(`Seed Family`, Fraction)) +
  geom_col() +
  facet_wrap( ~ Region, ncol = 2) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

