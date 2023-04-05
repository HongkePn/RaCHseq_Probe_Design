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
```{R}
targets=read.table('exon_all.tsv',sep = '\t')
colnames(targets)=c("chr","database","feature","start","end","score","strand","frame","isoform")
targets=targets[targets$feature=="exon",]
targets$exon_num <- targets$isoform

targets$isoform <- gsub("^.*transcript_name ","",targets$isoform)
targets$isoform <- gsub("; exon_number.*$","",targets$isoform)

targets$exon_num <- gsub("^.*exon_number ","",targets$exon_num)
targets$exon_num <- gsub("; exon_id.*$","",targets$exon_num)

targets$length <- abs(targets$start - targets$end)
```
