# RaCHseq Probe Design
## Files required for probe design
- A list of target genes
- gencode.v40.primary_assembly.annotation.gtf (Download from Gencode)
- GRCh38.primary_assembly.genome.fa (Download from Gencode)

## Extract exons of target gene
You can extract the exon information from gencode.v40.primary_assembly.annotation.gtf. 
\
For example, in **bash command line**: 
```
cat gencode.v40.primary_assembly.annotation.gtf | grep 'target_gene_name' > exon_all.tsv
# exon information are now stored in exon_all.tsv
```

All following steps can be done with **R**

## Extract sequences of target gene
Read in **exon_all.tsv** and prepare the information data.frame
```
pathway <- "pathway_to_your_exon.tsv"
targets=read.table(pathway,sep = '\t')

#extract useful information
colnames(targets)=c("chr","database","feature","start","end","score","strand","frame","isoform")
targets=targets[targets$feature=="exon",]
targets$exon_num <- targets$isoform
targets$isoform <- gsub("^.*transcript_name ","",targets$isoform)
targets$isoform <- gsub("; exon_number.*$","",targets$isoform)
targets$exon_num <- gsub("^.*exon_number ","",targets$exon_num)
targets$exon_num <- gsub("; exon_id.*$","",targets$exon_num)
targets$exon_num <- as.numeric(targets$exon_num)

#fix the chromsome name and make sure they are consistent with chr names in the fasta file
targets$chr <- gsub("chr","",targets$chr)
targets$chr[targets$chr=="X"]="23"
targets$chr[targets$chr=="Y"]="24"
targets$chr=as.numeric(targets$chr)

#take 1 probe / per kb, calculate how many probes we need to cover each exon
targets$length <- abs(targets$start - targets$end)
targets$probe_num <- (targets$length %/% 1000) + 1

#find the position of exons. 
#the first exon in a transcripts is called "start"
#the last exon in a transcripts is called "end"
#others are in the "middle"
targets$pos <- "middle"
for (i in unique(targets$isoform)) {
  x <- targets[targets$isoform == i,]
  x[1,]$pos <- "start"
  x[nrow(x),]$pos <- "end"
  targets[targets$isoform == i,] <- x
}
```
Extract sequences from GRCh38.primary_assembly.genome.fa. 
```{r}
extractSeq <- function(x) {
  assign('chr',x['chr'])
  chr <- as.numeric(chr)
  if (x['strand'] == "-") {
    reverseComplement(genome[[chr]][x[,'start']:x[,'end']])
  } else {
    genome[[chr]][x[,'start']:x[,'end']]
  }
  # seq <- DNAStringSet(seq)
  # names(seq) <- paste(x['isoform'], " exon_num", x['exon_num'])
  # writeXStringSet(seq, 'target_exons.fa', append = T)
}

for (i in 1:nrow(targets)) {
  x <- targets[i,]
  seq <- extractSeq(x)
  seq <- DNAStringSet(seq)
  names(seq) <- paste(x$isoform, "_exon_num_", x$exon_num, sep = "")
  writeXStringSet(seq,"all_exon.fa", append = T)
}
```



















