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
## Takes a genomic range (peak
findSeed2 = function(GR) {
    # Call findSeedMatches on a single pair
    # Note: Ensure arguments match findSeedMatches signature, 
    # you may need to convert sequences to DNAStringSet inside the function.
   #  noMatchDF = data.frame(type = character(), log.kd = integer(), p3.score = integer(), note = Rle())
    sequence = GR[1]$sequence
    sequence = DNAStringSet(sequence) %>% setNames(GR[1]$names)
    seed = mmu[GR[1]$miR]
    scan = scanMiR::findSeedMatches(seqs = sequence, seeds = seed, ret = "GRange", verbose = F) %>% as.data.frame()
    if(ncol(scan) == 5){scan = bind_cols(scan,noMatchDF)}
    scan = scan %>% mutate(peakName = GR[1]$names)
   # return(scan)
    GRtable = as.data.frame(mcols(GR)) 
    scan = left_join(scan, GRtable, by = c("peakName" = "names"))
 #  return(scan)
    scan = scan %>% mutate(seqnames = as.character(seqnames(GR[1])), 
                           seedSeq = substr(sequence, start, end),
                           start = start + start(GR[1]), 
                           end = end + start(GR[1]),
                           strand = as.character(strand(GR[1])),
                           note = as.character(note),
                           seedName = paste(seqnames, ":", start, "-", end, strand, sep = "")) %>%
                    dplyr::select(-c(sequence))
#return(scan)
    scan = makeGRangesFromDataFrame(scan, keep.extra.columns = T)
    seqlevels(scan) = seqlevels(GR)
    return(scan)
  }
filteredGR[13]

findSeed2(filteredGR[13])

seedMatches = bplapply( split(filteredGR[1:10]), findSeed2,BPPARAM = param)


## multiporcessing parameters
param <- MulticoreParam(workers = 30, progressbar = TRUE)

setwd("/projects/Retina/CLIP/AGO2/difClip/CDS_peaks")

## Prep a list for the results
peaksList = vector(mode = "list", length = 11)
names(peaksList) = c("UTR3_peaks", #peaks in the 3' UTR - genome mapped data
                     "CDS_peaks",  #peaks in the CDS - genome mapped data
                     "UTR5_peaks",  #peaks in the 5' UTR - genome mapped data
                     "intron_peaks"  #peaks in the introns - genome mapped data
#                     "intergenic_peaks", # interegnic peaks (not introns, CDS, UTRs) - genome mapped data
#                     "UTR3_peaksTX", #peaks in the 3' UTR - mapped from genome to transcript data
#                     "CDS_peaksTX",  #peaks in the CDS - mapped from genome to transcript data
#                     "UTR5_peaksTX"  #peaks in the 5' UTR - mapped from genome to transcript data
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
            reduce_ranges_directed() %>%
            mutate(segmentName = paste(seqnames, ":", start, "-", end, strand, sep =""))
            
CDS = CDS %>% plyranges::filter(seqnames %in% selectSeqLevels)
seqlevels(CDS)<-selectSeqLevels

threeUTRs<-threeUTRsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>% 
            unlist() %>%
            reduce_ranges_directed() %>%
            mutate(segmentName = paste(seqnames, ":", start, "-", end, strand, sep =""))

fiveUTRs<-fiveUTRsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>% 
            unlist() %>%
            reduce_ranges_directed()%>%
            mutate(segmentName = paste(seqnames, ":", start, "-", end, strand, sep =""))

threeUTRs<-threeUTRs %>% plyranges::filter(seqnames %in% selectSeqLevels)
seqlevels(threeUTRs)<-selectSeqLevels
#mcols(threeUTRs)$UTRlocation<-paste(seqnames(threeUTRs),":",start(threeUTRs),"-",end(threeUTRs),sep="")

fiveUTRs<-fiveUTRs %>% plyranges::filter(seqnames %in% selectSeqLevels)
seqlevels(fiveUTRs)<-selectSeqLevels
#mcols(fiveUTRs)$UTRlocation<-paste(seqnames(fiveUTRs),":",start(fiveUTRs),"-",end(fiveUTRs),sep="")

introns = intronsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>% 
            unlist()%>%
            reduce_ranges_directed()%>%
            mutate(segmentName = paste(seqnames, ":", start, "-", end, strand, sep =""))
            
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
         names=paste(seqnames,":",start,"-",end, strand, sep=""))

         
                  
#####################################
#### Find miR seeds in the CLIP peaks

## scanMir::findSeedMatches returns null when no matches are found and a table if matches are found
## We need a "noMatch" table to place instead of the NULLs so that bind_rows will work

noMatchDF = data.frame(type = character(), log_kd = integer(), p3.score = integer(), note = Rle())



# Find seeds for the miRNA in each CLIP peak
# UTR3_peaks = join_overlap_inner_directed(filteredGR,threeUTRs)
seedMatches = bplapply( split(filteredGR), findSeed2,BPPARAM = param)

seedMatches = GRangesList(seedMatches) %>% unlist()

## annotate the seedmatches
seedMatches =join_overlap_left_directed(seedMatches,proteinGenes)

# intersect with the different sequence segment types
peaksList[["UTR3_peaks"]] =  join_overlap_inner_directed(seedMatches,threeUTRs)
peaksList[["UTR5_peaks"]] =  join_overlap_inner_directed(seedMatches, fiveUTRs)
peaksList[["CDS_peaks"]] = join_overlap_inner_directed(seedMatches,CDS)
peaksList[["CDS_peaks"]] = subsetByOverlaps(peaksList[["CDS_peaks"]], peaksList[["UTR3_peaks"]], invert = T) # some sequences can be both coding and UTRs depending on transcript

peaksList[["intron_peaks"]] = join_overlap_inner_directed(seedMatches,introns)
peaksList[["intron_peaks"]] = subsetByOverlaps(peaksList[["intron_peaks"]], peaksList[["UTR3_peaks"]], invert = T) # some introns may contain exons
peaksList[["intron_peaks"]]= subsetByOverlaps(peaksList[["intron_peaks"]], peaksList[["UTR5_peaks"]], invert = T)
peaksList[["intron_peaks"]]= subsetByOverlaps(peaksList[["intron_peaks"]], peaksList[["CDS_peaks"]], invert = T)
#peaksList[["intergenic_peaks"]] = join_overlap_inner(seedMatches,intergenic)


## save a table for the 3' UTRs

peaksList[["UTR3_peaks"]] %>% as.data.frame() %>% dplyr::select(-entrezid) %>% write.csv(., "UTR3_mapped_miRNA_seeds.csv", row.names = F)







