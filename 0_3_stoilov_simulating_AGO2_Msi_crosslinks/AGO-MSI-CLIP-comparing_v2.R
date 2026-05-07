library("dplyr")
library("tidyr")
library("stringr")
library("plyranges") 
library("AnnotationHub")
library("ggplot2")
# library("openxlsx")
#library("mirbase.db")
library("ggforce")

library("scanMiR")
library("scanMiRData")
library("BSgenome.Mmusculus.UCSC.mm10")
library("Biostrings")
library("BiocParallel")
#library("miRBaseConverter")
library("rstatix")
library("forcats")
library("ggpubr")
library("ggtext")

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

themeCustom = function (){
  theme_light()+
    theme(axis.line = element_line(color = 'black', linewidth = 1),
          axis.ticks = element_line(color = 'black', linewidth = 1),
          text = element_text(family = "Arial"),
          axis.text.x = element_markdown(halign= 0, size = 11, face = "bold"),
          axis.text.y = element_markdown(size = 11, face= "bold"),
          axis.title = element_markdown(size =13, face = "bold"),
          legend.title = element_markdown(size =13, face = "bold.italic"),
          legend.text = element_markdown(size=11),
          strip.text.x = element_markdown(size = 11, color="black", face="bold"),
          strip.text.y = element_markdown(size = 13, color="black", face="bold.italic"),
          strip.background = element_rect(fill=NA))
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
                     )

## gather anotation and other external data we need
data(mmu, package = "scanMiRData")
ah<-AnnotationHub()
qr<-query(ah, c("EnsDb", "GRCm38"))

ensDB<-qr[["AH89211"]]
seqlevelsStyle(ensDB)<-"UCSC"


proteinGenes<-genes(ensDB, columns=c("gene_id","symbol","description","entrezid"), filter= ~ gene_biotype == "protein_coding") %>% 
  dplyr::filter(!is.na(entrezid), !is.na(gene_id))
selectSeqLevels<-seqlevels(proteinGenes)[str_detect(seqlevels(proteinGenes),regex("^chr.*[0-9XY]"))]


proteinGenes<-proteinGenes %>% dplyr::filter(seqnames %in% selectSeqLevels)
seqlevels(proteinGenes)<-selectSeqLevels


threeUTRs<-threeUTRsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>% 
            unlist() %>%
            reduce_ranges_directed() %>%
            mutate(segmentName = paste(seqnames, ":", start, "-", end, strand, sep =""))

threeUTRs<-threeUTRs %>% dplyr::filter(seqnames %in% selectSeqLevels)
seqlevels(threeUTRs)<-selectSeqLevels

allGenes = genes(ensDB, columns=c("gene_id","symbol","description","entrezid"))
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


fiveUTRs<-fiveUTRsByTranscript(ensDB,filter= ~ tx_biotype == "protein_coding") %>% 
            unlist() %>%
            reduce_ranges_directed()%>%
            mutate(segmentName = paste(seqnames, ":", start, "-", end, strand, sep =""))


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
seedMatches = bplapply( split(filteredGR), findSeed2,BPPARAM = param)

seedMatches = GRangesList(seedMatches) %>% unlist()
saveRDS(seedMatches, "seedMatches.rds")
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


## save a table for the 3' UTRs

peaksList[["UTR3_peaks"]] %>% as.data.frame() %>% dplyr::select(-entrezid) %>% write.csv(., "UTR3_mapped_miRNA_seeds.csv", row.names = F)


peaksList %>% saveRDS("mapped_miRNA_seeds.rds")


################################################################################################################################################
############################################################################################################################################
## recover from saved RDS
seedsList = readRDS("0_2_stoilov_microRNA_seed_mapping_with_genomic_locations/mapped_miRNA_seeds.rds")
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

############################## UAGs infull 3' UTRs
genome(threeUTRs) = "mm10"

threeUTRseq = threeUTRs %>% mutate(sequence = getSeq(BSgenome.Mmusculus.UCSC.mm10, .),
                                    names = paste(seqnames,start,end,strand, sep = ":"))
names(threeUTRseq$sequence) = threeUTRseq$names



############################
## Compare to Msi1 cross-linking sites


MSI1_sites = read.table("0_2_stoilov_microRNA_seed_mapping_with_genomic_locations/MSI1-with_input.sites.ucsc.bed", header = F)
colnames(MSI1_sites ) = c("seqnames","start","end","state","score","strand")
MSI1_sites = makeGRangesFromDataFrame(MSI1_sites, keep.extra.columns = T)

## annotate the x-link sites with the UTR data
UTR3_MSI = join_overlap_inner_directed(MSI1_sites,threeUTRs) %>% 
            as.data.frame()  %>%
            rename(all_of(c( start = "MSI.start", end = "MSI.end", width = "MSI.width")))
            
## combine the miRNA seeds with the MSI1 x-link sites
CLIP_x_table = inner_join(UTR3_seeds, UTR3_MSI, by = c("segmentName" = "segmentName", "strand" = "strand", "seqnames" = "seqnames"), relationship = "many-to-many") %>%
              dplyr::mutate(distance = floor(MSI.start - (seed.start+seed.end)/2) *( (strand == "+")*2 - 1), absolute.distance = abs(distance))

## annotate the regulation
CLIP_x_table = CLIP_x_table %>% dplyr::mutate(regulation = ifelse(padj<0.05 & log2FoldChange<0 , "downregulated", "undetermined"), 
                               regulation = ifelse(padj<0.05 & log2FoldChange>0 , "upregulated", regulation), 
                               regulation = ifelse(padj>0.95 & abs(log2FoldChange)<0.1 , "notregulated", regulation),
                               distance.bin = cut_interval(distance, length=5))

## pointrange plot of score in binned data - seems to be the most readable
## 
plot_data = CLIP_x_table %>% 
    dplyr::filter(distance >-200, distance <= 200, regulation %in% c("downregulated","notregulated")) %>% 
    dplyr::group_by(distance.bin, regulation) %>% 
    dplyr::summarise(sd = sd(score), se = sd/sqrt(dplyr::n()-1), n=dplyr::n(), score = mean(score))  %>% 
    dplyr::mutate(bin = forcats::fct_relabel(distance.bin, ~str_replace_all(.,"[\\[\\(\\]\\) ]","")),
                  regulation = ifelse(regulation == "downregulated", "Downregulated<br>peaks", regulation ),
                  regulation = ifelse(regulation == "notregulated", "Not regulated<br>peaks", regulation ))
plot_data
    ## calculate T-statistics
t_statistics = CLIP_x_table %>% 
    dplyr::filter(distance >-200, distance <= 200, regulation %in% c("downregulated","notregulated")) %>% 
    dplyr::mutate(regulation = ifelse(regulation == "downregulated", "Downregulated<br>peaks", regulation ),
                  regulation = ifelse(regulation == "notregulated", "Not regulated<br>peaks", regulation )) %>%
    dplyr::group_by(distance.bin) %>%
    rstatix::t_test(.,score ~ regulation) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj") 
  #  add_y_position()
    sig_data = plot_data %>% mutate(y.position = score + se + 10) %>% dplyr::group_by(bin, distance.bin) %>% dplyr::summarise(y.position = max(y.position))
    sig_data = left_join(sig_data, t_statistics, by = "distance.bin") %>% ungroup()
     
p = plot_data%>%
        ggplot( aes(x = bin, y= score, color= regulation)) + 
        #geom_col(position = position_dodge()) +
        geom_pointrange(aes(y= score, ymin = score -se, ymax = score + se), size =0.04)+ 
        #scale_x_discrete(labels = plot_data$bin)+
        # themeCustom()+
          scale_color_viridis_d(option="turbo",  begin = .25, end =0.75)+
          scale_x_discrete(breaks = function(x){x[c(T,F)]}, guide = guide_axis(minor.ticks = TRUE))+
        theme(axis.text.x = element_markdown(size = 6, angle = 60, vjust = 1, hjust=1, face = "bold"),
          axis.text.y = element_markdown(size = 7, face= "bold"),
          axis.title = element_markdown(size =8, face = "bold"),
          #legend.title = element_markdown(size =7, face = "bold.italic"),
          legend.text = element_markdown(size=6, face = "bold"),
          strip.text.x = element_markdown(size = 6, color="black", face="bold"),
          strip.text.y = element_markdown(size = 6, color="black", face="bold.italic"),
          legend.title = element_blank(),
          #legend.position = "top",
          legend.box.spacing = unit(0, "pt") ) +
        labs(x = "Distance from the miRNA seed",
             y = "MSI1 cross-link site score") 
        
p = p+  stat_pvalue_manual(sig_data, x="bin" , label = "p.adj.signif", hide.ns = T, remove.bracket = T, size=2)
ggsave("0_3_stoilov_simulating_AGO2_Msi_crosslinks/MSI1-xLinks_relative_to_miRNA-seed_sigLabel.png", plot = p, width = 1600, height = 800, units = "px")

plot_stat_data = left_join(plot_data, t_statistics, by = "distance.bin") %>% ungroup() %>% mutate(`-Log<sub>10</sub>(p.adj)` = ifelse(regulation == "Downregulated<br>peaks", -log10(p.adj),10))
p = plot_stat_data%>%
        ggplot( aes(x = bin, y= score, color= regulation, alpha = `-Log<sub>10</sub>(p.adj)`)) + 
        #geom_col(position = position_dodge()) +
        geom_pointrange(aes(y= score, ymin = score -se, ymax = score + se), size =0.05)+ 
        #scale_x_discrete(labels = plot_data$bin)+
        themeCustom()+
        scale_color_viridis_d(option="turbo",  begin = .25, end =0.75)+
        scale_x_discrete(breaks = function(x){x[c(T,F)]})+
        theme(axis.text.x = element_markdown(size = 6, angle = 60, vjust = 1, hjust=1, face = "bold"),
          axis.text.y = element_markdown(size = 7, face= "bold"),
          axis.title = element_markdown(size =8, face = "bold"),
          legend.title = element_markdown(size =7, face = "bold.italic"),
          legend.text = element_markdown(size=6, face = "bold"),
          strip.text.x = element_markdown(size = 6, color="black", face="bold"),
          strip.text.y = element_markdown(size = 6, color="black", face="bold.italic"),
#           legend.title = element_blank(),
#           legend.position = "top",
          legend.box.spacing = unit(0, "pt") ) +
        labs(x = "Distance from the miRNA seed",
             y = "MSI1 cross-link site score") 
        
# p = p+  stat_pvalue_manual(sig_data, x="bin" , label = "p.adj.signif", hide.ns = T, remove.bracket = T, size=2)
ggsave("MSI1-xLinks_relative_to_miRNA-seed_sigAlpha.png", plot = p, width = 2000, height = 800, units = "px")

themeCustom = function (){
  theme_light()+
    theme(axis.line = element_line(color = 'black', linewidth = 0.5),
          axis.ticks = element_line(color = 'black', linewidth = 0.3),
          text = element_text(family = "Arial"),
          axis.text.x = element_markdown(size = 6, angle = 60, vjust = 1, hjust=1, face = "bold"),
          axis.text.y = element_markdown(size = 7, face= "bold"),
          axis.title = element_markdown(size =8, face = "bold"),
          legend.title = element_markdown(size =7, face = "bold.italic"),
          legend.text = element_markdown(size=6, face = "bold"),
          strip.text.x = element_markdown(size = 6, color="black", face="bold"),
          strip.text.y = element_markdown(size = 6, color="black", face="bold.italic"),
          strip.background = element_rect(fill=NA))
}


CLIP_x_table %>% dplyr::group_by(regulation) %>% dplyr::summarise(n=dplyr::n())


## The distance from miRNA seed in the Ago2 Peak to the nearest Msi1 crosslink - some Ago2 peaks may have more than one seed

tmp = CLIP_x_table %>% dplyr::group_by(peakName, regulation) %>% dplyr::summarise(absolute.distance = min(absolute.distance)) 

summaryCLIP_x_table = left_join(tmp, CLIP_x_table, by = c("peakName" = "peakName", "regulation" = "regulation", "absolute.distance" = "absolute.distance"))

summaryCLIP_x_table %>% dplyr::group_by(regulation) %>% dplyr::summarise(n= dplyr::n()) 
summaryCLIP_x_table = summaryCLIP_x_table  %>%  dplyr::mutate(regulation = ifelse(regulation == "downregulated", "Downregulated <br>peaks", regulation ),
                  regulation = ifelse(regulation == "notregulated", "Not regulated <br>peaks", regulation ))


### Msi1 crosslink score
p = summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("Downregulated <br>peaks","Not regulated <br>peaks") ) %>%# dplyr::group_by(regulation) %>% dplyr::slice_sample(n=207) %>%
    ggplot( aes(x = regulation, y = score)) + geom_boxplot(outliers=F) + #+ ylim(-500,500)
      themeCustom() +
      theme(axis.title.x = element_blank(),
      axis.text.x = element_markdown(size = 8, angle = 60, vjust = 1, hjust=1, face = "bold")) +
      labs(y = "Closest MSI1 cross-link <br> site score") 
    

p = p + stat_compare_means(method = "t.test", aes(label = ..p.signif..), 
                        label.x = 1.25, label.y = 50, size =5)
ggsave("MSI1-closest_peak-score.png", plot = p, width = 500, height = 800, units = "px")

## Closest Msi1 crosslink distance
p = summaryCLIP_x_table %>% 
        dplyr::filter(regulation %in% c("Downregulated <br>peaks","Not regulated <br>peaks") ) %>%# dplyr::group_by(regulation) %>% dplyr::slice_sample(n=207) %>%
        ggplot( aes(x = regulation, y = absolute.distance)) + 
          geom_boxplot(outliers=F) + #+ ylim(-500,500)+
          themeCustom() +
          theme(axis.title.x = element_blank(),
          axis.text.x = element_markdown(size = 8, angle = 60, vjust = 1, hjust=1, face = "bold")) +
          labs(y = "Distance to the closes<br>MSI1 cross-link site") 

p = p + stat_compare_means(method = "t.test", aes(label = ..p.signif..), 
                        label.x = 1.25, label.y = 285, size = 5)
ggsave("MSI1-closest_peak-distance.png", plot = p, width = 500, height = 800, units = "px")



## calculate statistics for various metrics and store them in a table.

p_table = data.frame(mean.score.difference = numeric(), 
                     score.p = numeric(), 
                     mean.absolute.distance.difference = numeric(), 
                     absolute.distance.p = numeric(), 
                     mean.distance.difference = numeric(), 
                     distance.p = numeric(), 
                     anova.score.Df = numeric(),
                     anova.score.F = numeric(),
                     anova.score.p = numeric(),
                     anova.absolute.distance.Df = numeric(),
                     anova.absolute.distance.F = numeric(),
                     anova.absolute.distance.p = numeric(),
                     anova.distance.Df = numeric(),
                     anova.distance.F = numeric(),
                     anova.distance.p = numeric())



p_table[1,c("mean.score.difference","score.p")] = summaryCLIP_x_table %>% 
                                                          ungroup() %>% 
                                                          dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>%  
                                                          rstatix::t_test(score ~ regulation, detailed = T) %>% 
                                                          dplyr::select(estimate,p) 

p_table[1,c("mean.absolute.distance.difference","absolute.distance.p")] = summaryCLIP_x_table %>% 
                                                          ungroup() %>% 
                                                          dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>%  
                                                          rstatix::t_test(absolute.distance ~ regulation, detailed = T) %>% 
                                                          dplyr::select(estimate,p) 

p_table[1,c("mean.distance.difference","distance.p")] = summaryCLIP_x_table %>% 
                                                          ungroup() %>% 
                                                          dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>%  
                                                          rstatix::t_test(distance ~ regulation, detailed = T) %>% 
                                                          dplyr::select(estimate,p) 
                                                          
p_table[1,c("anova.score.Df","anova.score.F","anova.score.p")] = summaryCLIP_x_table %>% 
                                                              ungroup() %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>% 
                                                              aov(score ~ regulation, data =.) %>% 
                                                              summary() %>% 
                                                              unlist() %>% .[c("Df1","F value1","Pr(>F)1")]
                                                              
                                                              
p_table[1,c("anova.absolute.distance.Df","anova.absolute.distance.F","anova.absolute.distance.p")] = summaryCLIP_x_table %>% 
                                                              ungroup() %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>% 
                                                              aov(absolute.distance ~ regulation, data =.) %>% 
                                                              summary() %>% 
                                                              unlist() %>% .[c("Df1","F value1","Pr(>F)1")]
                                                              
                                                              
p_table[1,c("anova.distance.Df","anova.distance.F","anova.distance.p")] = summaryCLIP_x_table %>% 
                                                              ungroup() %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>% 
                                                              aov(distance ~ regulation, data =.) %>% 
                                                              summary() %>% 
                                                              unlist() %>% .[c("Df1","F value1","Pr(>F)1")]
                                                              
p_table = p_table %>% pivot_longer(cols = everything(),names_to = "metric", values_to = "value") 




### Use permutations of the data to calculate the probability 
  ## calculate stats for the two groups and compare them
  
  
## puls two random sets and performs statistical test
## used in repeated testing to determine the significance of the results
## iter parameter is a placeholder so that we can call it with bplapplay()
randomSetTest = function(iter = 1, inputSet, setSize1, setSize2){
    #print(head(inputSet))
  ## generate two random sets and merge them in a table
  seedNames = inputSet %>% distinct(seedName)
  tmp1 = seedNames %>% dplyr::slice_sample(n=setSize1) %>% dplyr::mutate(regulation = "set1")
  tmp2 = seedNames %>% dplyr::filter(!(seedName %in% tmp1$seedName)) %>%  dplyr::slice_sample(n=setSize1)%>% dplyr::mutate(regulation = "set2") 
  tmp = bind_rows(tmp1, tmp2)
  inputSet = inputSet %>% dplyr::select(-regulation)
  rndSet = left_join(tmp, inputSet, by = c("seedName" = "seedName")) 
  tmp = rndSet %>% dplyr::group_by(seedName, regulation) %>% dplyr::summarise(absolute.distance = min(absolute.distance), .groups = "drop") %>% ungroup()
  summary_rndSet_table = left_join(tmp, rndSet, by = c("seedName" = "seedName", "regulation" = "regulation",  "absolute.distance" = "absolute.distance"))
  
  stat_table = data.frame(mean.score.difference = numeric(), 
                     score.p = numeric(), 
                     mean.absolute.distance.difference = numeric(), 
                     absolute.distance.p = numeric(), 
                     mean.distance.difference = numeric(), 
                     distance.p = numeric(), 
                     anova.score.Df = numeric(),
                     anova.score.F = numeric(),
                     anova.score.p = numeric(),
                     anova.absolute.distance.Df = numeric(),
                     anova.absolute.distance.F = numeric(),
                     anova.absolute.distance.p = numeric(),
                     anova.distance.Df = numeric(),
                     anova.distance.F = numeric(),
                     anova.distance.p = numeric())
  stat_table[1,c("mean.score.difference","score.p")] = summary_rndSet_table %>% 
                                                          rstatix::t_test(score ~ regulation, detailed = T) %>% 
                                                          dplyr::select(estimate,p) 

  stat_table[1,c("mean.absolute.distance.difference","absolute.distance.p")] = summary_rndSet_table %>% 
                                                          rstatix::t_test(absolute.distance ~ regulation, detailed = T) %>% 
                                                          dplyr::select(estimate,p) 

  stat_table[1,c("mean.distance.difference","distance.p")] = summary_rndSet_table %>% 
                                                          rstatix::t_test(distance ~ regulation, detailed = T) %>% 
                                                          dplyr::select(estimate,p) 
                                                          
  stat_table[1,c("anova.score.Df","anova.score.F","anova.score.p")] = summary_rndSet_table %>% 
                                                              aov(score ~ regulation, data =.) %>% 
                                                              summary() %>% 
                                                              unlist() %>% .[c("Df1","F value1","Pr(>F)1")]
                                                              
                                                              
  stat_table[1,c("anova.absolute.distance.Df","anova.absolute.distance.F","anova.absolute.distance.p")] = summary_rndSet_table %>% 
                                                              aov(absolute.distance ~ regulation, data =.) %>% 
                                                              summary() %>% 
                                                              unlist() %>% .[c("Df1","F value1","Pr(>F)1")]
                                                              
                                                              
  stat_table[1,c("anova.distance.Df","anova.distance.F","anova.distance.p")] = summary_rndSet_table %>% 
                                                              aov(distance ~ regulation, data =.) %>% 
                                                              summary() %>% 
                                                              unlist() %>% .[c("Df1","F value1","Pr(>F)1")]
  return(stat_table)
}

## number of permutations
iterNum = 100000
iterations = c(1:iterNum)

# set the group sizes to the sizes from the analysis
setSizes = summaryCLIP_x_table %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>% dplyr::group_by(regulation) %>% dplyr::summarise(n=dplyr::n()) %>% pull(n)

# calculate the statistics for the selected number of permutations
iterList =bplapply(iterations, randomSetTest, inputSet = CLIP_x_table, setSize1 = setSizes[1], setSize2 = setSizes[2], BPPARAM=param)
iterTable  = bind_rows(iterList)

# combine the permutations table with the analysis data and calculate the probability that the permuted data will produce the obtained result
iterTable  = iterTable %>% pivot_longer(cols = everything(),names_to = "metric", values_to = "iter_val") %>% left_join(.,p_table, by = "metric")

permutations_summary = iterTable %>% dplyr::group_by(metric, value) %>% dplyr::reframe(p = ifelse(str_detect(metric,".p$|distance.difference$"), sum(iter_val <= value)/dplyr::n(), NA),
                                                                p = ifelse(str_detect(metric,"mean.score.difference$"), sum(iter_val >= value)/dplyr::n(), p)) %>% distinct()
                                                                  
write.csv(permutations_summary, file = "permutations_statistics_summary.csv", row.names = F)
                                                                  
permutations_summary_p = iterTable %>% dplyr::group_by(metric, value) %>% dplyr::reframe(p = ifelse(str_detect(metric,".p$|distance.difference$"), max(1/iterNum,sum(iter_val <= value)/dplyr::n()), NA),
                                                                p = ifelse(str_detect(metric,"mean.score.difference$"), max(1/iterNum,sum(iter_val >= value)/dplyr::n()), p)) %>% distinct()
write.csv(permutations_summary_p, file = "permutations_statistics_summary_p.csv", row.names = F)

##### check for sequence composition bias and triplet frequencies

tmp = CLIP_x_table %>% dplyr::select(seqnames, seed.start, seed.end, strand, regulation) 
colnames(tmp) = c("seqnames", "start", "end", "strand", "regulation")

tmp = makeGRangesFromDataFrame(tmp, keep.extra.columns = T)
tmp = tmp %>% reduce_ranges_directed(regulation = min(regulation))

tmp = tmp %>%  plyranges::stretch(200)
tmp = tmp %>% mutate(sequence = getSeq(BSgenome.Mmusculus.UCSC.mm10, .), 
         names=paste(seqnames,":",start,"-",end, strand, sep=""))

## base frequencies show slightly A/T richer sequence
baseFreqs = alphabetFrequency(tmp$sequence, as.prob = T)
reg = tmp %>% as.data.frame() %>% dplyr::select(regulation)
baseFreqs = bind_cols(reg,baseFreqs)
baseFreqs = baseFreqs %>% group_by(regulation) %>% summarise(n= dplyr::n(),across(!starts_with("r"), ~ mean(.x,na.rm = T)))
baseFreqs 




## TAG is the most enriched triplet in the vicinity of the seed

tripletFreqs = trinucleotideFrequency(tmp$sequence, as.prob = T)
tripletFreqs = bind_cols(reg,tripletFreqs)


tripletFreqs = tripletFreqs %>% 
                  pivot_longer(matches("^[ACGT]"), names_to = "triplet", values_to = "frequency")
                  


tripletFreqsP = tripletFreqs %>% 
                    dplyr::filter(regulation %in% c("downregulated","notregulated") ) %>% 
                    dplyr::group_by(triplet) %>% 
                    rstatix::t_test(frequency ~ regulation) %>% arrange(p)
tripletFreqsP$padj = p.adjust(tripletFreqsP$p)

## plot triplet frequencies
tripletFreqStats = tripletFreqs %>% 
                  group_by(regulation, triplet) %>% 
                  summarise(n = dplyr::n(), mean = mean(frequency), sd = sd(frequency), se = sd/sqrt(n-1)) 
                  

tripletFreqStas = tripletFreqStats %>% 
                  separate(triplet, into = c("triplet","metric"), sep ="_") %>%
                  pivot_wider(names_from = "metric", values_from = "frequency")
                  
                  
                  
tripletFreqStats %>% dplyr::filter(regulation %in% c("downregulated","notregulated") )  %>%
  ggplot(aes(x=triplet, y = mean, fill = regulation)) + 
      geom_col(position = position_dodge()) +
      geom_errorbar(aes(ymax = mean+se, ymin = mean-se), position = position_dodge(), width = 0.9)
      
      
pentaFreq=oligonucleotideFrequency(tmp$sequence, as.prob = T, width=6)
         
pentaFreq = bind_cols(reg,pentaFreq)

pentaFreq = pentaFreq %>% 
                  pivot_longer(matches("^[ACGT]"), names_to = "word", values_to = "frequency")
     
pentaFreqP = pentaFreq %>% 
                    dplyr::filter(regulation %in% c("downregulated","notregulated") ) %>% 
                    dplyr::group_by(word) %>% 
                    rstatix::t_test(frequency ~ regulation) %>% arrange(p)
pentaFreqP$padj = p.adjust(pentaFreqP$p)
pentaFreqP %>% dplyr::filter(statistic>0) %>%print(n=40)



                               
                               
                               
