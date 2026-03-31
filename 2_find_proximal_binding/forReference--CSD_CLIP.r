library("dplyr")
library("tidyr")
library("stringr")
library("plyranges")
library("AnnotationHub")


library("BSgenome.Mmusculus.UCSC.mm10")
library("Biostrings")
library("BiocParallel")



## multiporcessing parameters
param <- MulticoreParam(workers = 30, progressbar = TRUE)

setwd("/projects/Retina/CLIP/AGO2/difClip/CDS_peaks")

## gather anotation and other external data we need
ah<-AnnotationHub()
qr<-query(ah, c("EnsDb", "GRCm38"))

ensDB<-qr[["AH89211"]]
seqlevelsStyle(ensDB)<-"UCSC"

# pull protein coding genes
proteinGenes<-genes(ensDB, columns=c("gene_id","symbol","description"), filter= ~ gene_biotype == "protein_coding")

# keep only the standard chromosomes
selectSeqLevels<-seqlevels(proteinGenes)[str_detect(seqlevels(proteinGenes),regex("^chr.*[0-9XY]"))]
proteinGenes<-proteinGenes %>% plyranges::filter(seqnames %in% selectSeqLevels)
seqlevels(proteinGenes)<-selectSeqLevels



threeUTRs<-threeUTRsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>%
            unlist() %>%
            reduce_ranges_directed()  # merges overlaping UTRs

threeUTRs<-threeUTRs %>% plyranges::filter(seqnames %in% selectSeqLevels)
seqlevels(threeUTRs)<-selectSeqLevels
#mcols(threeUTRs)$UTRlocation<-paste(seqnames(threeUTRs),":",start(threeUTRs),"-",end(threeUTRs),sep="")

UTR3 = split(threeUTRs)

for (i in 1:100){
  UTR3[[i]] = join_overlap_left_directed(UTR3[[i]], tmp)
}


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


########################################################################### IMPORTANT
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

