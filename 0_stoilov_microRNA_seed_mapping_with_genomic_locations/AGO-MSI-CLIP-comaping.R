library("dplyr")
library("tidyr")
library("stringr")
library("plyranges") 
library("AnnotationHub")
library("ggplot2")
# library("openxlsx")
# library("mirbase.db")
# library("ggforce")

# library("scanMiR")
# library("scanMiRData")
library("BSgenome.Mmusculus.UCSC.mm10")
library("Biostrings")
library("BiocParallel")
# library("miRBaseConverter")
# library("rstatix")

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






################################################################################################################################################
############################################################################################################################################
## temp copy to recreate threeUTRs starting from this point in the code
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


## recover from saved RDS
seedsList = readRDS("0_stoilov_microRNA_seed_mapping_with_genomic_locations/mapped_miRNA_seeds.rds")
## read 
difClip = read_tsv("1_AGO2_motif_analysis/diff_chimeric_clip.zip")

## add differential CLIP data from MSI1/MSI2 knockout
difClip = difClip %>% 
          dplyr::select(seqnames, start, end, width, strand, log2FoldChange, padj) %>%
          dplyr::filter(!is.na(padj))
difClip = makeGRangesFromDataFrame(difClip, keep.extra.columns = T)


UTR3_seeds  = join_overlap_inner_directed(seedsList[["UTR3_peaks"]] , difClip) %>%
              #mcols(.) %>%
              as.data.frame() %>% 
              rename(all_of(c( start = "seed.start", end = "seed.end", width = "seed.width")))


## grab the MSI1 cross-link sites and build peaks by expansing the sites with 10nt in each direction and overlapping them

## read the cross-link sites
MSI1_sites = read.table("0_stoilov_microRNA_seed_mapping_with_genomic_locations/MSI1-with_input.sites.ucsc.bed", header = F)
colnames(MSI1_sites )<-c("seqnames","start","end","state","score","strand")
MSI1_sites = makeGRangesFromDataFrame(MSI1_sites, keep.extra.columns = T)

## generate peaks by grabing a sequences cenetered on the cross-link sites and then reducing the overlaps
MSI1_peaks = MSI1_sites %>% plyranges::stretch(8) %>% reduce_ranges_directed(score = sum(score), n.xlink.sites=n())


## find the UAG sites in the peaks
MSI1_peaks = MSI1_peaks %>% mutate(sequence = getSeq(BSgenome.Mmusculus.UCSC.mm10, .),
                                    names = paste(seqnames,start,end,strand, sep = ":"))
names(MSI1_peaks$sequence) = MSI1_peaks$names


UAG_pos = vmatchPattern("TAG", MSI1_peaks$sequence) %>% unlist() %>% as.data.frame()

##what fraction of the peaks has a UAG?
sum(names(MSI1_peaks$sequence) %in% UAG_pos$names)/length(MSI1_peaks)


## Annotate the UAG poistions with UTR data
UAG_pos = UAG_pos %>% 
            tidyr::separate(names, into = c("seqnames","seq.start","seq.end", "strand"), convert = T, sep =":", remove = F) %>% 
            dplyr::mutate(start = start+seq.start, end = end+seq.start) %>%
            makeGRangesFromDataFrame(., keep.extra.columns = F)

## UAGs in 3' UTRS
UTR3_UAG = join_overlap_inner_directed(UAG_pos,threeUTRs) %>% 
            as.data.frame() %>%
            rename(all_of(c( start = "UAG.start", end = "UAG.end", width = "UAG.width")))
            
### Merge the seeds and UAG tables and calculate the distances
### Keep the distances integers
CLIP_x_table = inner_join(UTR3_seeds, UTR3_UAG, by = c("segmentName" = "UTRID", "strand" = "strand", "seqnames" = "seqnames"), relationship = "many-to-many") %>%
              dplyr::mutate(distance = floor((UAG.start+UAG.end)/2 - (seed.start+seed.end)/2), absolute.distance = abs(distance))
              
### annotate the regulation
CLIP_x_table = CLIP_x_table %>% dplyr::mutate(regulation = ifelse(padj<0.05 & log2FoldChange<0 , "downregulated", "undetermined"), 
                               regulation = ifelse(padj<0.05 & log2FoldChange>0 , "upregulated", regulation), 
                               regulation = ifelse(padj>0.95 & abs(log2FoldChange)<0.1 , "notregulated", regulation))

                               
## Histogram plots of the UAGs relative to the seed
 ggplot(CLIP_x_table, aes(x = distance)) + geom_histogram(binwidth = 5) + xlim(-10000,10000) + facet_wrap(~regulation, scales = "free_y")
 
 
ggplot(CLIP_x_table, aes(x = regulation, y = distance)) + geom_violin() + ylim(-10000,10000) 

 
## Find the UAGs closest to the seed
tmp = CLIP_x_table %>% dplyr::group_by(seedName, regulation) %>% dplyr::summarise(absolute.distance = min(absolute.distance)) 

summaryCLIP_x_table = left_join(tmp, CLIP_x_table, by = c("seedName" = "seedName", "regulation" = "regulation", "absolute.distance" = "absolute.distance"))
summaryCLIP_x_table %>% dplyr::group_by(regulation) %>% dplyr::summarise(n= dplyr::n())

summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") ) %>%
ggplot( aes(x = distance)) + geom_histogram(binwidth = 5) + xlim(-1000,1000) + facet_wrap(~regulation, scales = "free_y")

summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") ) %>% #dplyr::group_by(regulation) %>% dplyr::slice_sample(n=207) %>%
ggplot( aes(x = regulation, y = distance)) + geom_violin() + ylim(-1000,1000)

## Same as above but group by peak instead of seed, as one peak can have more than one seeds
tmp = CLIP_x_table %>% dplyr::group_by(peakName, regulation) %>% dplyr::summarise(absolute.distance = min(absolute.distance)) 

summaryCLIP_x_table = left_join(tmp, CLIP_x_table, by = c("peakName" = "peakName", "regulation" = "regulation", "absolute.distance" = "absolute.distance"))

summaryCLIP_x_table %>% dplyr::group_by(regulation) %>% dplyr::summarise(n= dplyr::n())


summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") ) %>% #dplyr::group_by(regulation) %>% dplyr::slice_sample(n=207) %>%
ggplot( aes(x = distance)) + geom_histogram(binwidth = 5) + xlim(-1000,1000) + facet_wrap(~regulation, scales = "free_y") + theme_clean()

## match the sample sizes
summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") ) %>% dplyr::group_by(regulation) %>% dplyr::slice_sample(n=207) %>%
ggplot( aes(x = distance)) + geom_histogram(binwidth = 5) + xlim(-1000,1000) + facet_wrap(~regulation, scales = "free_y")

summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") ) %>% dplyr::group_by(regulation) %>% dplyr::slice_sample(n=207) %>%
ggplot( aes(x = regulation, y = distance)) + geom_violin() #+ ylim(-500,500)


### stats
summaryCLIP_x_table %>% ungroup() %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>%  rstatix::t_test(distance ~ regulation)

summaryCLIP_x_table %>% ungroup() %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>% aov(distance ~ regulation, data =.) %>% summary()

summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>% ungroup()  %>% var.test(distance ~ regulation, data =.)

summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>% dplyr::group_by(regulation) %>% dplyr::slice_sample(n=207) %>% ungroup()  %>% var.test(distance ~ regulation, data =.)




############################
## As abovem but using Msi1 cross linking sites
## straight cross-link sites

MSI1_sites = read.table("0_stoilov_microRNA_seed_mapping_with_genomic_locations/MSI1-with_input.sites.ucsc.bed", header = F)
colnames(MSI1_sites ) = c("seqnames","start","end","state","score","strand")
MSI1_sites = makeGRangesFromDataFrame(MSI1_sites, keep.extra.columns = T)

## annotate the x-link sites with the UTR data
UTR3_MSI = join_overlap_inner_directed(MSI1_sites,threeUTRs) %>% 
            as.data.frame()  %>%
            rename(all_of(c( start = "MSI.start", end = "MSI.end", width = "MSI.width")))
            
## combine the miRNA seeds with the MSI1 x-link sites
CLIP_x_table = inner_join(UTR3_seeds, UTR3_MSI, by = c("segmentName" = "UTRID", "strand" = "strand", "seqnames" = "seqnames"), relationship = "many-to-many") %>%
              dplyr::mutate(distance = floor(MSI.start - (seed.start+seed.end)/2), absolute.distance = abs(distance))

## annotate the regulation
CLIP_x_table = CLIP_x_table %>% dplyr::mutate(regulation = ifelse(padj<0.05 & log2FoldChange<0 , "downregulated", "undetermined"), 
                               regulation = ifelse(padj<0.05 & log2FoldChange>0 , "upregulated", regulation), 
                               regulation = ifelse(padj>0.95 & abs(log2FoldChange)<0.1 , "notregulated", regulation),
                               distance.bin = cut_interval(distance, length=10))

ggplot(CLIP_x_table, aes(x = distance)) + geom_histogram(binwidth = 5) + xlim(-10000,10000) + facet_wrap(~regulation, scales = "free_y")
ggplot(CLIP_x_table, aes(x = regulation, y = distance)) + geom_violin() + ylim(-1000,1000) 


CLIP_x_table %>% dplyr::group_by(distance, regulation) %>% dplyr::summarise(score = sum(score)/dplyr::n()) %>%
ggplot( aes(x = distance, y= score)) + geom_line() + xlim(-10000,10000) + facet_wrap(~regulation) #, scales = "free_y")

CLIP_x_table %>% dplyr::group_by(distance, regulation) %>% dplyr::summarise(score = sum(score)/dplyr::n()) %>%
ggplot( aes(x = distance, y= score)) + geom_smooth(method = "loess", span=0.2) + xlim(-1000,1000) + facet_wrap(~regulation) #, scales = "free_y")

CLIP_x_table %>% dplyr::group_by(distance, regulation) %>% dplyr::summarise(score = sum(score)/dplyr::n()) %>% 
ggplot( aes(x = distance, y= score)) + geom_smooth(method = "loess", span=0.2) + xlim(-1000,1000) + facet_wrap(~regulation) #, scales = "free_y")

CLIP_x_table %>% dplyr::group_by(regulation) %>% dplyr::summarise(n=dplyr::n())

## distance from miRNA seed seed
tmp = CLIP_x_table %>% dplyr::group_by(seedName, regulation) %>% dplyr::summarise(absolute.distance = min(absolute.distance)) 

summaryCLIP_x_table = left_join(tmp, CLIP_x_table, by = c("seedName" = "seedName", "regulation" = "regulation", "absolute.distance" = "absolute.distance"))

summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") ) %>%
ggplot( aes(x = distance)) + geom_histogram(binwidth = 5) + xlim(-1000,1000) + facet_wrap(~regulation, scales = "free_y")

summaryCLIP_x_table %>% dplyr::group_by(distance, regulation) %>% dplyr::summarise(score = sum(score)/dplyr::n()) %>%
ggplot( aes(x = distance, y= score)) + geom_line() + xlim(-500,500) + facet_wrap(~regulation, scales = "free_y")

summaryCLIP_x_table %>% dplyr::group_by(distance, regulation) %>% dplyr::summarise(score = sum(score)/dplyr::n()) %>%
ggplot( aes(x = distance, y= score)) + geom_smooth(method = "loess", span=0.1, se=F) + xlim(-1000,1000) + facet_wrap(~regulation) #, scales = "free_y")

summaryCLIP_x_table %>% dplyr::group_by(regulation) %>% dplyr::summarise(n=dplyr::n())

## The shortest distance from miRNA seed in an Ago2 peak - some peaks may have more than one seed

tmp = CLIP_x_table %>% dplyr::group_by(peakName, regulation) %>% dplyr::summarise(absolute.distance = min(absolute.distance)) 

summaryCLIP_x_table = left_join(tmp, CLIP_x_table, by = c("peakName" = "peakName", "regulation" = "regulation", "absolute.distance" = "absolute.distance"))

summaryCLIP_x_table %>% dplyr::group_by(regulation) %>% dplyr::summarise(n= dplyr::n())


summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") ) %>%
ggplot( aes(x = distance)) + geom_histogram(binwidth = 5) + xlim(-1000,1000) + facet_wrap(~regulation) #, scales = "free_y")
#
summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") ) %>% 
ggplot( aes(x = regulation, y = distance)) + geom_violin() #+ ylim(-500,500)

summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") ) %>% dplyr::group_by(regulation) %>% dplyr::slice_sample(n=207) %>%
ggplot( aes(x = regulation, y = distance)) + geom_violin() + labs(x = "AGO2 Regulation", y="Distance Between Closest AGO2/Msi-1 Pair") + theme_minimal()
  #+ ylim(-500,500)
ggsave("0_stoilov_microRNA_seed_mapping_with_genomic_locations/varianceLower.svg")


summaryCLIP_x_table %>% dplyr::group_by(distance, regulation) %>% dplyr::summarise(score = sum(score)/dplyr::n()) %>%
ggplot( aes(x = distance, y= score)) + geom_smooth(method = "loess", span=0.1, se=F) + xlim(-1000,1000) + facet_wrap(~regulation) #, scales = "free_y")

#
summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") ) %>% dplyr::group_by(regulation) %>% dplyr::slice_sample(n=207) %>%
ggplot( aes(x = regulation, y = score)) + geom_boxplot() #+ ylim(-500,500)



summaryCLIP_x_table %>% dplyr::group_by(regulation) %>% dplyr::summarise(n=dplyr::n())


summaryCLIP_x_table %>% ungroup() %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>%  rstatix::t_test(score ~ regulation)

summaryCLIP_x_table %>% ungroup() %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>%  rstatix::t_test(distance ~ regulation)

summaryCLIP_x_table %>% ungroup() %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>% aov(distance ~ regulation, data =.) %>% summary()

summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>% ungroup()  %>% var.test(distance ~ regulation, data =.)
summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>% ungroup() %>% group_by(regulation) %>% summarize(var = var(distance))
qf(0.975, 203, 1516)
qf(0.025, 203, 1516)


summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>% dplyr::group_by(regulation) %>% dplyr::slice_sample(n=207) %>% ungroup()  %>% var.test(distance ~ regulation, data =.)


summaryCLIP_x_table %>% dplyr::group_by(regulation) %>% dplyr::summarise(dist = mean(distance))
