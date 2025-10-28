# trying to extract gene letters from the ranges
library(tidyverse)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)

clip_data <- read_tsv("diff_chimeric_clip.zip")
clip_data
clip_granges <- clip_data %>% as_granges()

clip_seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, clip_granges)
clip_granges$sequence <- clip_seqs
