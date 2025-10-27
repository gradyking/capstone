library("Rsubread")
library("dplyr")
library("tidyr")
library("stringr")
# library("edgeR")
# library("TxDb.Mmusculus.UCSC.mm10.ensGene")
# library("org.Mm.eg.db")
library("plyranges") 
# library("Mus.musculus")
library("DESeq2")
library("AnnotationHub")
library("ggplot2")
library("Rsamtools")    

setwd("/projects/Retina/CLIP/AGO2/difClip")
bamPath = "../files/"
# TxDb(Mus.musculus) = TxDb.Mmusculus.UCSC.mm10.ensGene

#chimeric peaks
#peaksBed = read.table(file="chimericPeaks.bed",header=F,stringsAsFactors=F,sep="\t")
peaksBed = read.table(file="../files/chimeric_peaks.bed",header=F,stringsAsFactors=F,sep="\t")
colnames(peaksBed) = c("Chr","Start","End","miR","Score","Strand")
peaksBed = peaksBed %>% tidyr::unite(GeneID, c(Chr,Start,End,Strand), remove=F) 


peaksGR = makeGRangesFromDataFrame(peaksBed,seqnames.field="Chr", start.field="Start", end.field="End", strand.field="Strand", keep.extra.columns=T)    
peaksGR = reduce_ranges_directed(peaksGR, miR=paste(miR,collapse=";"),Score=paste(Score, collapse=";"))
mcols(peaksGR)$GeneID = paste(seqnames(peaksGR),":",start(peaksGR),"-",end(peaksGR)," ",strand(peaksGR), sep="")
         

#### (ensDB, columns=c("gene_id","symbol","description","entrezid"))
ah = AnnotationHub()
qr = query(ah, c("EnsDb", "GRCm38"))
ensDB = qr[["AH89211"]]
seqlevelsStyle(ensDB) = "UCSC"

exon = exons(ensDB, columns=c("gene_id"))
exon = exon %>% group_by(gene_id) %>% reduce_ranges_directed()
exonSAF =  as.data.frame(exon) %>% dplyr::rename(Chr=seqnames, Start=start,End=end, Strand=strand, GeneID=gene_id) 
genes = genes(ensDB, columns=c("gene_id","symbol","description","entrezid"))
peaksGRAnnotated = join_overlap_left_directed(peaksGR,genes) %>% 
    group_by(GeneID, miR, Score, gene_id,symbol,description) %>% 
    reduce_ranges_directed()


####

# Ugly hack to solve a problem where one peak can be asigned to two
# gene ids. For example a peak in micro-RNA can be asiagned to the mature 
# miR and to the host gene
peaksGRAnnotated = peaksGRAnnotated[!duplicated(peaksGRAnnotated)]
names(peaksGRAnnotated) = mcols(peaksGRAnnotated)$GeneID

###

peaksSAF = as.data.frame(peaksGRAnnotated) %>% dplyr::rename(Chr=seqnames, Start=start,End=end, Strand=strand)


targetsGene = read.table('targetsGene_input.txt', header=T,sep='\t')
targetsPeak = read.table('targetsPeak_IP.txt', header=T,sep='\t')
targets = bind_rows(targetsGene,targetsPeak)


countsGene = featureCounts(files=targetsGene$File, annot.ext=exonSAF, isGTFAnnotationFile=F, countMultiMappingReads=F, isPairedEnd=F, strandSpecific=1,useMetaFeatures = TRUE, nthreads=16)
colnames(countsGene$counts) = targetsGene$Target

countsIP = featureCounts(files=targetsPeak$File, annot.ext=peaksSAF, isGTFAnnotationFile=F, countMultiMappingReads=F, isPairedEnd=F, strandSpecific=1,useMetaFeatures = TRUE, nthreads=16, allowMultiOverlap=T)
colnames(countsIP$counts) = targetsPeak$Target


peakToGene = peaksSAF %>% dplyr::select(GeneID, gene_id)
ipCounts = as.data.frame(countsIP$counts) %>% dplyr::mutate(GeneID=rownames(.))
ipCounts = dplyr::left_join(ipCounts, peakToGene, by=c("GeneID"="GeneID"))

geneCounts = as.data.frame(countsGene$counts) %>% dplyr::mutate(gene_id=rownames(.))
ipCounts = dplyr::left_join(ipCounts, geneCounts, by=c("gene_id"="gene_id")) %>% dplyr::filter(!is.na(gene_id))
rownames(ipCounts) = ipCounts$GeneID
ipCounts = ipCounts %>% dplyr::select(-c(GeneID,gene_id)) %>% as.matrix()


############# using agregate gene input

deSet = DESeqDataSetFromMatrix(countData = ipCounts,
                              colData = targets,
                              rowRanges=peaksGRAnnotated[rownames(ipCounts)],
#                               rowData = ipAnnotation,
                              design= ~ Type + Group + Type:Group)

deSet = DESeq(deSet, test="LRT", reduced= ~ Type + Group)
resultsNames(deSet)  



resGR = results(deSet, format="GRangesList", saveCols=c("symbol","miR","Score","description","gene_id"))

threeUTRs = threeUTRsByTranscript(ensDB) %>% unlist() %>% reduce_ranges_directed()
utr3 = filter_by_overlaps(resGR, threeUTRs) %>% names(.)

cdsTX = cdsBy(ensDB,by="tx") %>% unlist() %>% reduce_ranges_directed()
cds = filter_by_overlaps(resGR, cdsTX) %>% names(.)

fiveUTRs = fiveUTRsByTranscript(ensDB) %>% unlist() %>% reduce_ranges_directed()
utr5 = filter_by_overlaps(resGR, fiveUTRs) %>% names(.)

# reduce-ranges is rather generous as it will include exons in the intronic ranges
# we deal with this during the annotation of the "feature" column:
# start with NA, then annotate introns, then UTRs (will overwrite some of the introns,
# then CDS (will overwrite exons annotated as introns and some of the UTRs)
intronsTX = intronsByTranscript(ensDB) %>% unlist() %>% reduce_ranges_directed()
intron = filter_by_overlaps(resGR, intronsTX) %>% names(.)

resGR = resGR %>% mutate(Feature=NA)
resGR = resGR %>% mutate(Feature=ifelse((names(.) %in% intron),"Intron",Feature))
resGR = resGR %>% mutate(Feature=ifelse((names(.) %in% utr5),"5'-UTR",Feature))
resGR = resGR %>% mutate(Feature=ifelse((names(.) %in% utr3),"3'-UTR",Feature))
resGR = resGR %>% mutate(Feature=ifelse((names(.) %in% cds),"CDS",Feature))

#

# resGR = resGR %>% mutate(threeUTR=ifelse((names(.) %in% utr3),T,F))

res = as.data.frame(resGR) %>% dplyr::mutate(peak=rownames(.)) %>% dplyr::arrange(padj) 
res %>% write.table(., file="diff_chimeric_clip.tsv", row.names=F, sep="\t")



## plot fold change for eCLIP peaks, color by location, transparency by adjusted p-value
res %>% dplyr::filter(padj<0.05, !is.na(Feature)) %>% # removes non-significant changes and non-coding transcripts)
        dplyr::select(log2FoldChange,padj,peak,Feature) %>%
        ggplot(aes(y=log2FoldChange, x=reorder(peak,log2FoldChange),fill=Feature, alpha=padj)) +
            geom_bar(stat="identity")+
            scale_alpha(range = c(1, 0.2))+
            ggtitle("Differentially enriched chimeric-Ago2 eCLIP peaks (Musashi knockout/WT)")+
            ylab("log2 Fold Change KO/WT")+
            theme_classic()+
            theme(
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.title.x=element_blank()
                
            )



res %>% dplyr::filter(padj<0.05, !is.na(Feature)) %>% # removes non-significant changes and non-coding transcripts)
        dplyr::select(log2FoldChange,padj,peak,Feature) %>%
        ggplot(aes(y=log2FoldChange, x=reorder(peak,log2FoldChange),fill=Feature)) +
            geom_bar(stat="identity")+
            scale_alpha(range = c(1, 0.2))+
            ggtitle("Differentially enriched chimeric-Ago2 eCLIP peaks (Musashi knockout/WT)")+
            ylab("log2 Fold Change KO/WT")+
            theme_classic()+
            theme(
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.title.x=element_blank()
                
            )
ggsave('dif_chim-clip_ranked.pdf')

###############################################
## Extract UTR sequences
dna = ah[["AH88475"]]   

utr3SigUp = join_overlap_inner(threeUTRs,resGR) %>% 
            plyranges::filter(Feature=="3'-UTR", log2FoldChange>=1, padj<=0.05) %>%
            reduce()
seqlevelsStyle(utr3SigUp) = "Ensembl"


utr3Contr = join_overlap_inner(threeUTRs,resGR) %>% 
            plyranges::filter(Feature=="3'-UTR", log2FoldChange>=-0.2, log2FoldChange<=0.2, padj>0.8, baseMean>=50) %>%
            reduce() %>% sample(200)
seqlevelsStyle(utr3Contr) = "Ensembl"        

       


getSeq(dna, utr3SigUp) %>% writeXStringSet(.,'utr3SigUp.fa',format='fasta')

getSeq(dna, utr3Contr) %>% writeXStringSet(.,'utr3Contr.fa',format='fasta')
            
##   
utr3SigUp = mutate(anchor_center(resGR), width=50) %>% join_overlap_inner(.,threeUTRs)%>% 
            plyranges::filter(Feature=="3'-UTR", log2FoldChange>=1, padj<=0.05) %>%
            reduce()
seqlevelsStyle(utr3SigUp) = "Ensembl"
mcols(utr3SigUp) = utr3SigUp %>% mutate(name=paste("mm10_",seqnames,":",start,"-",end,":",strand, sep=""))
names(utr3SigUp) = mcols(utr3SigUp)$name
utr3Contr = mutate(anchor_center(resGR), width=50) %>% join_overlap_inner(.,threeUTRs)%>% 
            plyranges::filter(Feature=="3'-UTR", log2FoldChange>=-0.2, log2FoldChange<=0.2, padj>0.8, baseMean>=50) %>%
            reduce() %>% sample(200) 
seqlevelsStyle(utr3Contr) = "Ensembl" 
utr3Contr = utr3Contr %>% mutate(name=paste("mm10_",seqnames,"_",start,"_",end,"_",strand, sep=""))
names(utr3Contr) = mcols(utr3Contr)$name

getSeq(dna, utr3SigUp) %>% writeXStringSet(.,'utr3SigUp_50nt.fa',format='fasta')
getSeq(dna, utr3Contr) %>% writeXStringSet(.,'utr3Contr_50nt.fa',format='fasta')
            
###############################################
## compare to protein expression (just exploring)
            
# read MS protein expression change

msDataFile = "/projects/Retina/MSI/proteomics/out//10ppm_noPhos/MSI_WT-KO_msconvert.txt"

protExpr = read.table(msDataFile, sep="\t",header=T)


## Photoreceptor genes
geneClusters = read.table('/projects/Retina/MSI/proteomics/flat_retina_clusters_Macosko.txt',header=T, sep="\t")

geneClusters = geneClusters %>% separate(Clusters, into=c("ClusterA","ClusterB"),sep="\\+",convert=T)

prClusters = geneClusters %>% dplyr::filter((ClusterA %in% c(24,25)) | (ClusterB %in% c(24,25)) & !(ClusterA %in% c(24,25)) & (ClusterB %in% c(24,25)))
gcA = prClusters %>% dplyr::filter(ClusterA %in% c(24,25))
gcB = prClusters %>% dplyr::filter(ClusterB %in% c(24,25))
gcB = gcB %>% mutate(myDif=(-myDif)) 
colnames(gcB)[c(4,5)] = colnames(gcB)[c(5,4)]

prClusters = bind_rows(gcA,gcB)
prGenes = prClusters %>% dplyr::group_by(gene) %>% summarise(`Photoreceptor protein`=all(myDif>0)) %>% dplyr::select(gene, `Photoreceptor protein`)
protExpr = left_join(protExpr, prGenes, by = c("symbol" = "gene"))
            
#########################

sigClip = res %>% dplyr::filter( Feature == "3'-UTR", !(symbol %in% c("Msi1","Msi2","Esr")), !is.na(padj)) 
protExpr = left_join(protExpr, sigClip[,c("gene_id","log2FoldChange", "padj")], by = c("gene" = "gene_id")) %>% 
    dplyr::mutate(`Ago2 Binding` = ifelse(log2FoldChange>0 & padj < 0.05, "Increased",ifelse(log2FoldChange<0 & padj <0.05 , "Reduced", "Unchanged"))) %>%
    dplyr::group_by(gene, symbol, logFC, sca.adj.pval) %>%
    summarize(padj = min(padj), `Ago2 Binding` = ifelse(all(`Ago2 Binding` %in% c("Reduced", "Unchanged")) & any(`Ago2 Binding` == "Reduced"), "Reduced", 
                                                    ifelse( all(`Ago2 Binding` %in% c("Increased","Unchanged")) & any(`Ago2 Binding` == "Increased"), "Increased",
                                                    ifelse( all(`Ago2 Binding` == "Unchanged"), "Unchanged","Mixed"))))


#     
# sigClip = res %>% dplyr::filter( Feature == "3'-UTR", !(symbol %in% c("Msi1","Msi2","Esr")), !is.na(padj), padj <0.05) 
# protExpr = left_join(protExpr, sigClip[,c("gene_id","log2FoldChange", "padj")], by = c("gene" = "gene_id")) %>% 
#     dplyr::mutate(`Ago2 Binding` = ifelse(log2FoldChange>0 & padj < 0.05, "Increased",ifelse(log2FoldChange<0 & padj <0.05 , "Reduced", "Unchanged"))) %>%
#     dplyr::distinct()
# 
# protExpr %>% dplyr::filter(!is.na(`Ago2 Binding`), `Ago2 Binding` != "Unchanged") %>%
# ggplot(aes(logFC, -log10(adj.P.Val), color = `Ago2 Binding`)) +
#     geom_point()

protExpr %>% dplyr::filter(!is.na(`Ago2 Binding`), `Ago2 Binding` != "Unchanged") %>%
ggplot(aes(logFC, -log10(sca.adj.pval))) +
    geom_point(aes(, color = `Ago2 Binding`)) +
  geom_text_repel( data = subset(protExpr, sca.adj.pval<=0.05 & abs(logFC)>0.5 & !is.na(`Ago2 Binding`) & `Ago2 Binding` != "Unchanged"),
                  aes(label = symbol)) +
  #guides(Direction = "legend", alpha = "none") +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(title = "Protein expression from mRNA with differential Ago2 binding to 3'-UTRs", x ="Log2 Fold change in Msi1, Msi2 double knockout", y = "Log10 adjusted p-value")

ggsave("Protein_expr_from_diff_Ago2-CLIP.png")

protExpr %>% dplyr::filter(!is.na(`Ago2 Binding`), !is.na(`Photoreceptor protein`)) %>%
ggplot(aes(`Photoreceptor protein`, logFC, color = `Ago2 Binding`)) +
    geom_sina()





#########################
res = res %>% mutate(miR= str_split(miR,";"), Score= str_split(Score,";")) %>% 
    unnest(c(miR,Score)) %>% 
    mutate(Score=as.numeric(Score))%>%
    group_by(gene_id,Feature) %>% 
    dplyr::filter(Score==max(Score)) %>%
    dplyr::rename(logFC.eClip=log2FoldChange)

### for each transcript select the peak with maximum log2FC for eCLIP    
resSummary = res %>% #dplyr::filter(Feature %in% c("CDS", "3'-UTR", "5'-UTR" )) %>% 
    dplyr::group_by(gene_id, Feature) %>% 
    dplyr::filter(Score==max(Score)) %>% 
    dplyr::distinct() %>% #head() %>% as.data.frame()
    dplyr::group_by(gene_id) %>%
    summarise(logFC.eClip= logFC.eClip[which.max(abs(logFC.eClip))], 
        Feature= Feature[which.max(abs(logFC.eClip))], 
        padj= padj[which.max(abs(logFC.eClip))],
        miR=miR[which.max(abs(logFC.eClip))]
        ) 

protClip = left_join(protExpr,resSummary, by=c("gene"="gene_id")) %>% dplyr::filter(!is.na(logFC.eClip))

protClip %>%
    ggplot(aes(x=logFC.eClip, y=-logFC, color=miR))+
    geom_point() +
    facet_wrap(~Feature)



    
####### Do not select data
protClip = left_join(protExpr,res, by=c("gene"="gene_id")) %>% dplyr::filter(!is.na(logFC.eClip))

protClip %>% #dplyr::filter(Feature %in% c("3'-UTR" )) %>%
    ggplot(aes(x=logFC.eClip, y=-logFC, color=miR))+
    geom_point() +
    facet_wrap(~Feature)



  
######### using peak input (meh)
peaksBed = read.table(file="eCLIP_peaksall.merged.bed",header=F,stringsAsFactors=F,sep="\t")
colnames(peaksBed) = c("Chr","Start","End","IDR.feat","IDR.glob","Strand")
peaksBed = peaksBed %>% tidyr::unite(GeneID, c(Chr,Start,End,Strand), remove=F) 

peaksSAF = peaksBed %>% 
            tidyr::unite(GeneID, c(Chr,Start,End,Strand), remove=F) %>%
            dplyr::select(GeneID, Chr, Start, End, Strand)

files = list.files(path=bamPath, pattern=".*rmDup\\.sorted\\.bam$")

targets = str_split(files, pattern="_", n=4, simplify=T) %>% 
            as.data.frame() %>% 
            select(1,2,3) %>%
            dplyr::rename("Target"=V1,"Group"=V2,"Replicate"=V3) %>%
            mutate(Type=str_extract(Target,"Input|IP"))
            
            
files = paste(bamPath,files,sep="")

counts = featureCounts(files=files, annot.ext=peaksSAF, isGTFAnnotationFile=F, countMultiMappingReads=F, isPairedEnd=F, strandSpecific=1,useMetaFeatures = TRUE, nthreads=16, allowMultiOverlap=T)

colnames(counts$counts) = targets$Target

filteredCounts = counts$counts %>% as.data.frame() %>% filter_all(all_vars(.>10))  %>% as.matrix()

deSet = DESeqDataSetFromMatrix(countData = filteredCounts,
                              colData = targets,
                              design= ~ Type + Group + Type:Group)

deSet = DESeq(deSet, test="LRT", reduced= ~ Type + Group)
resultsNames(deSet)

results(deSet,name= "Group_FF_vs_CRE") %>% as.data.frame() %>% arrange(padj) %>% head(40)


 
