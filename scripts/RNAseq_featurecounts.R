library(tidyverse)
# BiocManager::install('Rsubread')
library(Rsubread)
library(readr)
library(dplyr)
dir.create("counts")

aligned3 <-
  featureCounts(
    files=list.files("alignments3", pattern="bam$", full.names=T),
    annot.ext="GCF_000024265.1_ASM2426v1_genomic.gff",
    isGTFAnnotationFile=TRUE,
    GTF.featureType=c("CDS","riboswitch", "RNase_P_RNA", "rRNA", "SRP_RNA", "tmRNA", "tRNA"),
    GTF.attrType="ID", #open GFF file and decide which attribute to use for name
    minMQS=1,
    strandSpecific=0,
    isPairedEnd=TRUE,
    nthreads=4,
    verbose=TRUE
  )

save(aligned3, file = "counts/FeatureCountsResultHisat2.rda")