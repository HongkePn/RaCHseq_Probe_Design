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
#code you need to modify =======================================
pathway <- "pathway_to_your_exon.tsv"
targets=read.table(pathway,sep = '\t')

#code you can copy straight forward ============================
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

```
#code you need to modify =======================================
fasta_pathway <- "pathway_to_your_GRCh38.primary_assembly.genome.fa"
save_exon_fasta_pathway <- "pathway_to_where_you_want_to_save_extracted_exon_sequences"

#code you can copy straight forward ============================
library(Biostrings)
#function to extract sequences
extractSeq <- function(x) {
  assign('chr',x['chr'])
  chr <- as.numeric(chr)
  if (x['strand'] == "-") {
    reverseComplement(genome[[chr]][x[,'start']:x[,'end']])
  } else {
    genome[[chr]][x[,'start']:x[,'end']]
  }
}

genome <- readDNAStringSet(fasta_pathway)
for (i in 1:nrow(targets)) {
  x <- targets[i,]
  seq <- extractSeq(x)
  seq <- DNAStringSet(seq)
  names(seq) <- paste(x$isoform, "_exon_num_", x$exon_num, sep = "")
  writeXStringSet(seq, save_exon_fasta_pathway, append = T)
}
```

# Prepare the exon fasta files for probe design
Probes are **120bp** in length. 
- For exon >= 120bp, use the original sequence
- For exon < 120bp, concatenate the short exon to last exon and next exon. Then design probes to fully cover the short exon. 
```
#code you need to modify =======================================
origin_concatenation_pathway <- "pathway_to_where_you_want_to_save_original_exon_sequences"
up_concatenation_pathway <- "pathway_to_where_you_want_to_save_up_concatenations"
down_concatenation_pathway <- "pathway_to_where_you_want_to_save_down_concatenations"

#code you can copy straight forward ============================
#make the origin fasta (containing all exon >= 120bp) and up-concatenations (short exon concatenated to last exon)
for (i in unique(targets$isoform)) {
  x <- targets[targets$isoform == i,]
  x <- x[order(x$exon_num, decreasing = F),]
  seq.last01 <- DNAString("")
  seq.last02 <- DNAString("")
  for (j in 1:nrow(x)) {
    y <- x[j,]
    if(y$length < 120){
      seq.current <- extractSeq(y)
      seq.concatenation <- c(seq.last01, seq.current)
      
      if (length(seq.concatenation) < 120){ # if the up-concatenated exons is still < 120bp, 
        seq.concatenation <- c(seq.last02, seq.concatenation)
        seq.write <- DNAStringSet(seq.concatenation)
        names(seq.write) <- paste(y$isoform, "_exon_num", y$exon_num, "_POS", y$pos, "_L", y$length, "_PN", y$probe_num, "_up2", sep = "")
        writeXStringSet(seq.write, up_concatenation_pathway, append = T)
      } else {
        seq.write <- DNAStringSet(seq.concatenation)
        names(seq.write) <- paste(y$isoform, "_exon_num", y$exon_num, "_POS", y$pos, "_L", y$length, "_PN", y$probe_num, "_up1", sep = "")
        writeXStringSet(seq.write, up_concatenation_pathway, append = T)
      }
      seq.last02 <- seq.last01
      seq.last01 <- seq.current
    } else {
      seq.o <- extractSeq(y)
      seq.write <- DNAStringSet(seq.o)
      names(seq.write) <- paste(y$isoform, "_exon_num", y$exon_num, "_POS", y$pos, "_L", y$length, "_PN", y$probe_num, "_origin", sep = "")
      writeXStringSet(seq.write, origin_concatenation_pathway, append = T)
      seq.last02 <- seq.last01
      seq.last01 <- seq.o
    }
  }
}

#make down-concatenations (short exon concatenated to next exon)
for (i in unique(targets$isoform)) {
  x <- targets[targets$isoform == i,]
  x <- x[order(x$exon_num, decreasing = T),]
  seq.last01 <- DNAString("")
  seq.last02 <- DNAString("")
  for (j in 1:nrow(x)) {
    y <- x[j,]
    if(y$length < 120) {
      seq.current <- extractSeq(y)
      seq.concatenation <- c(seq.current,seq.last01)
      if(length(seq.concatenation) < 120){ # if the down-concatenated exons is still < 120bp,
        seq.concatenation <- c(seq.concatenation,seq.last02)
        seq.write <- DNAStringSet(seq.concatenation)
        names(seq.write) <- paste(y$isoform, "_exon_num", y$exon_num, "_POS", y$pos, "_L", y$length, "_PN", y$probe_num, "_down2", sep = "")
        writeXStringSet(seq.write, down_concatenation_pathway, append = T)
      } else {
        seq.write <- DNAStringSet(seq.concatenation)
        names(seq.write) <- paste(y$isoform, "_exon_num", y$exon_num, "_POS", y$pos, "_L", y$length, "_PN", y$probe_num, "_down1", sep = "")
        writeXStringSet(seq.write, down_concatenation_pathway, append = T)
      }
      seq.last02 <- seq.last01
      seq.last01 <-seq.current
    } else {
      seq.current <- extractSeq(y)
      seq.last02 <- seq.last01
      seq.last01 <- seq.current
    }
  }
}
```

# Make probes through IDT web tool
We can use IDT *Custom Hybridization Capture Panels* to make probes against our fasta files. 
Here is the [link](https://sg.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-hybridization-capture/custom-hyb-panels) to the web tool. 

Please follow the instruction there and submit your job. The unfiltered probe list will be sent to your email box soon. 

# Select the probes

<div align=center>
<img src = "https://github.com/HongkePn/RaCHseq_Probe_Design/blob/main/probe_selection.png">
</div>

\
Aims of this step: 
- select 1 probe / per kb
- for up-concatenations, select a probe from the right side.
- for down-concatenations, select a probe from the left side.
- select probes with similar GC%

```
#code you need to modify =======================================
IDT_unfilter_probe_list_pathway <- "pathway_to_IDT_unfilter_probe_list"
probe.seq <- read.csv(IDT_unfilter_probe_list_pathway)

#code you can copy straight forward ============================
#prepare the data.frame for following probe selection

#exon names
probe.seq$exon_name <- gsub("_POS.*$","",probe.seq$Chromosome)
#probe num, how many probes we need to cover this exon
probe.seq$probe_num <- gsub("^.*_PN","",probe.seq$Chromosome)
probe.seq$probe_num <- gsub("_down.*$","",probe.seq$probe_num)
probe.seq$probe_num <- gsub("_up.*$","",probe.seq$probe_num)
probe.seq$probe_num <- gsub("_origin.*$","",probe.seq$probe_num)
probe.seq$probe_num <- as.numeric(probe.seq$probe_num)
#concatenation
probe.seq$concatenation <- gsub("^.*down.*$","down",probe.seq$Chromosome)
probe.seq$concatenation <- gsub("^.*up.*$","up",probe.seq$concatenation)
probe.seq$concatenation <- gsub("^.*origin.*$","origin",probe.seq$concatenation)
#GC% similarity
probe.seq$drift <- abs(probe.seq$GC - mean(probe.seq$GC))
#length of exon
probe.seq$length <- gsub("^.*_L","",probe.seq$Chromosome)
probe.seq$length <- gsub("_PN.*$","",probe.seq$length)
probe.seq$length <- as.numeric(probe.seq$length)

#probe selection
probe.select <- data.frame()
probe.up <- dplyr::filter(probe.seq, concatenation == "up")
for (i in unique(probe.up$exon_name)) {
  x <- probe.up[probe.up$exon_name == i,]
  x <- x[order(x$Start),]
  probe.select <- rbind(probe.select,x[nrow(x),])
}

probe.down <- dplyr::filter(probe.seq, concatenation == "down")
for (i  in unique(probe.down$exon_name)) {
  x <- probe.down[probe.down$exon_name == i,]
  x <- x[order(x$Start),]
  probe.select <- rbind(probe.select,x[1,])
}

probe.origin <- dplyr::filter(probe.seq, concatenation == "origin")
for (i in unique(probe.origin$exon_name)) {
  x <- probe.origin[probe.origin$exon_name == i,]
  pn <- x[1,]$probe_num
  if(pn == 1){
    x <- x[order(x$drift),]
    probe.select <- rbind(probe.select, x[1,])
  } else if (pn > 1) {
    x <- x[order(x$Start),]
    len <- x$length[1]
    n <- len %/% 1000
    for (j in 0:n) {
      g1 <- j*1000
      g2 <- (j+1)*1000
      y <- x[x$Start>g1 & x$Stop<g2,]
      y <- y[order(y$drift),]
      probe.select <- rbind(probe.select, y[1,])
    }
  }
}

#remove duplicates
probe.select=probe.select[!duplicated(probe.select$Seq),]
#remove NA
probe.select=probe.select[!is.na(probe.select$Seq),]
```

# save selected probe list

```
#code you need to modify =======================================
selected_probe_pathway <- "pathway_to_where_you_want_to_save_selected_probes"

#code you can copy straight forward ============================
write.csv(probe.select,file = selected_probe_pathway,row.names=F)
```
the selected probe list are now ready to go. 
You can send them to companies, sush as IDT, to get your probe libraries. Good Luck :D










