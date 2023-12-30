Amplicon Analysis

Ben Mellin 12/19/2023

Amplicon Analysis

- I will be analyzing an amplicon dataset generated from DNA found in crushed basalts. The samples were collected in an underwater mountain in the pacific ocean. I will be using a subset of this [data:https://www.frontiersin.org/articles/10.3389/fmicb.2015.01470/full.](data:https://www.frontiersin.org/articles/10.3389/fmicb.2015.01470/full)
- [I will be closely following the analysis performed here: https://astrobiomike.github.io/amplicon/dada2_workflow_ex#extracting-the- standard-goods-from-dada2](https://astrobiomike.github.io/amplicon/dada2_workflow_ex#extracting-the-standard-goods-from-dada2)
- I hope to gain a basic understanding of amplicon analysis in R.

Setting up our conda environment

The following code was executed in terminal.

conda install -y -c conda-forge mamba![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.001.png)

mamba create -y -n hb-dada2-ex-wf -c conda-forge -c bioconda -c defaults \

`             `cutadapt=2.3 r-base=3.6.3 rstudio=1.1.456 r-tidyverse=1.3.0 \

`             `r-vegan=2.5 r-dendextend=1.14.0 r-viridis=0.6 \

`             `bioconductor-phyloseq=1.30 bioconductor-deseq2=1.26 bioconductor-dada2=1.14 \              bioconductor-decipher=2.14 bioconductor-decontam=1.6 r-biocmanager=1.30 \

`             `r-matrix=1.3\_2 libopenblas=0.3.7

conda activate hb-dada2-ex-wf

Now we get the data and make a file with the sample names. The following was exectued in terminal.

cd ~![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.002.png)

curl -L -o dada2\_amplicon\_ex\_workflow.tar.gz https://ndownloader.figshare.com/files/28773936 tar -xzvf dada2\_amplicon\_ex\_workflow.tar.gz

rm dada2\_amplicon\_ex\_workflow.tar.gz

cd dada2\_amplicon\_ex\_workflow/

ls \*\_R1.fq | cut -f1 -d "\_" > samples 

Now we must remove all the primers from our samples. Our primers can be found in the primers.fa file.

**for** sample **in** $(cat samples)![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.003.png)

**do**

echo "On sample: $sample"

`    `cutadapt -a ^GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC \              -A ^GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC \

`             `-m 215 -M 285 --discard-untrimmed \

`             `-o ${sample}\_sub\_R1\_trimmed.fq.gz -p ${sample}\_sub\_R2\_trimmed.fq.gz \ ${sample}\_sub\_R1.fq ${sample}\_sub\_R2.fq \

`             `>> cutadapt\_primer\_trimming\_stats.txt 2>&1

**done**

Check and see what fraction of reads and base pairs were retained in each sample.

paste samples <(grep "passing" cutadapt\_primer\_trimming\_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" ![ref1]cutadapt\_primer\_trimming\_stats.txt | cut -f3 -d "(" | tr -d ")")

Data Processing in R

Next we move to R to begin processing our data.

**library**(dada2)![ref2]

\## Loading required package: Rcpp![ref3]

\## Warning: replacing previous import 'lifecycle::last\_warnings' by ## 'rlang::last\_warnings' when loading 'tibble'![ref4]

\## Warning: replacing previous import 'lifecycle::last\_warnings' by ## 'rlang::last\_warnings' when loading 'pillar'![ref4]

packageVersion("dada2") # import packages ![ref2]## [1] '1.14.0'![ref3]

#setwd("~/dada2\_amplicon\_ex\_workflow") # set working directory![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.008.png)

\## first we're setting a few variables we're going to use ##

- one with all sample names, by scanning our "samples" file we made earlier samples <- scan("samples", what="character")
- we create one holding the file names of all the forward reads forward\_reads <- paste0(samples, "\_sub\_R1\_trimmed.fq.gz")
- and one with the reverse

  reverse\_reads <- paste0(samples, "\_sub\_R2\_trimmed.fq.gz")

- and variables holding file names for the forward and reverse
- filtered reads we're going to generate below filtered\_forward\_reads <- paste0(samples, "\_sub\_R1\_filtered.fq.gz") filtered\_reverse\_reads <- paste0(samples, "\_sub\_R2\_filtered.fq.gz")

Now we’ll take a look at the quality profile of our untrimmed reads

\## Scale for 'y' is already present. Adding another scale for 'y', which will ## replace the existing scale.![ref4]

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.009.jpeg)

\## Scale for 'y' is already present. Adding another scale for 'y', which will ## replace the existing scale.![ref4]

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.010.jpeg)

\## Scale for 'y' is already present. Adding another scale for 'y', which will ## replace the existing scale.![ref4]

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.011.jpeg)

There is a dropoff in quality of the reads starting at around length 200. We will therefore cut the forward reads to length 250 and reverse reads to length 200. Our goal is to keep the quality scores above 30.

filtered\_out <- filterAndTrim(forward\_reads, filtered\_forward\_reads,![ref5]

`                              `reverse\_reads, filtered\_reverse\_reads, maxEE=c(2,2),                               rm.phix=TRUE, minLen=175, truncLen=c(250,200))

Now lets take a look at the quality profile of our trimmed sequences.

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.013.jpeg)

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.014.jpeg)

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.015.jpeg)

Now lets take a look at an error model of our dataset.

err\_forward\_reads <- learnErrors(filtered\_forward\_reads, multithread=TRUE) ![ref6]

\## 45623000 total bases in 182492 reads from 20 samples will be used for learning the error rates. err\_reverse\_reads <- learnErrors(filtered\_reverse\_reads, multithread=![ref7]TRUE) ![ref6]

\## 36498400 total bases in 182492 reads from 20 samples will be used for learning the error rates. ![ref7]Now we take a look at how well our estimated error matches up with our observed error.

plotErrors(err\_forward\_reads, nominalQ=TRUE)![ref6]

\## Warning: Transformation introduced infinite values in continuous y-axis ## Transformation introduced infinite values in continuous y-axis![ref8]

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.019.jpeg)

plotErrors(err\_reverse\_reads, nominalQ=TRUE)![ref6]

\## Warning: Transformation introduced infinite values in continuous y-axis ## Transformation introduced infinite values in continuous y-axis![ref8]

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.020.jpeg)

Now we perform dereplication to remove identical sequences.

derep\_forward <- derepFastq(filtered\_forward\_reads, verbose=TRUE)![ref6]

\## Dereplicating sequence entries in Fastq file: B1\_sub\_R1\_filtered.fq.gz ## Encountered 552 unique sequences from 1498 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: B2\_sub\_R1\_filtered.fq.gz ## Encountered 224 unique sequences from 529 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: B3\_sub\_R1\_filtered.fq.gz ## Encountered 186 unique sequences from 457 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: B4\_sub\_R1\_filtered.fq.gz![ref7]

\## Encountered 203 unique sequences from 475 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: BW1\_sub\_R1\_filtered.fq.gz ## Encountered 723 unique sequences from 2109 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: BW2\_sub\_R1\_filtered.fq.gz ## Encountered 2760 unique sequences from 5527 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: R10\_sub\_R1\_filtered.fq.gz ## Encountered 5501 unique sequences from 10354 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: R11BF\_sub\_R1\_filtered.fq.gz ## Encountered 3452 unique sequences from 8028 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: R11\_sub\_R1\_filtered.fq.gz ## Encountered 4846 unique sequences from 8138 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: R12\_sub\_R1\_filtered.fq.gz ## Encountered 9747 unique sequences from 14423 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: R1A\_sub\_R1\_filtered.fq.gz ## Encountered 6894 unique sequences from 10906 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: R1B\_sub\_R1\_filtered.fq.gz ## Encountered 9380 unique sequences from 14672 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: R2\_sub\_R1\_filtered.fq.gz![ref7]

\## Encountered 9310 unique sequences from 15660 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: R3\_sub\_R1\_filtered.fq.gz![ref7]

\## Encountered 10225 unique sequences from 15950 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: R4\_sub\_R1\_filtered.fq.gz![ref7]

\## Encountered 10407 unique sequences from 17324 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: R5\_sub\_R1\_filtered.fq.gz![ref7]

\## Encountered 10981 unique sequences from 16728 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: R6\_sub\_R1\_filtered.fq.gz![ref7]

\## Encountered 8520 unique sequences from 13338 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: R7\_sub\_R1\_filtered.fq.gz![ref7]

\## Encountered 4432 unique sequences from 7331 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: R8\_sub\_R1\_filtered.fq.gz![ref7]

\## Encountered 6171 unique sequences from 11192 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: R9\_sub\_R1\_filtered.fq.gz![ref7]

\## Encountered 4234 unique sequences from 7853 total sequences read.![ref7]

names(derep\_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.021.png)

derep\_reverse <- derepFastq(filtered\_reverse\_reads, verbose=TRUE)

\## Dereplicating sequence entries in Fastq file: B1\_sub\_R2\_filtered.fq.gz ## Encountered 514 unique sequences from 1498 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: B2\_sub\_R2\_filtered.fq.gz ## Encountered 203 unique sequences from 529 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: B3\_sub\_R2\_filtered.fq.gz ## Encountered 171 unique sequences from 457 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: B4\_sub\_R2\_filtered.fq.gz![ref7]

\## Encountered 189 unique sequences from 475 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: BW1\_sub\_R2\_filtered.fq.gz ## Encountered 660 unique sequences from 2109 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: BW2\_sub\_R2\_filtered.fq.gz ## Encountered 2506 unique sequences from 5527 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: R10\_sub\_R2\_filtered.fq.gz ## Encountered 5054 unique sequences from 10354 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: R11BF\_sub\_R2\_filtered.fq.gz ## Encountered 3113 unique sequences from 8028 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: R11\_sub\_R2\_filtered.fq.gz ## Encountered 4568 unique sequences from 8138 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: R12\_sub\_R2\_filtered.fq.gz ## Encountered 9288 unique sequences from 14423 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: R1A\_sub\_R2\_filtered.fq.gz ## Encountered 6445 unique sequences from 10906 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: R1B\_sub\_R2\_filtered.fq.gz ## Encountered 8799 unique sequences from 14672 total sequences read.![ref7]![ref7]

\## Dereplicating sequence entries in Fastq file: R2\_sub\_R2\_filtered.fq.gz![ref7]

\## Encountered 8762 unique sequences from 15660 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: R3\_sub\_R2\_filtered.fq.gz![ref7]

\## Encountered 9710 unique sequences from 15950 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: R4\_sub\_R2\_filtered.fq.gz![ref7]

\## Encountered 9613 unique sequences from 17324 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: R5\_sub\_R2\_filtered.fq.gz![ref7]

\## Encountered 10432 unique sequences from 16728 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: R6\_sub\_R2\_filtered.fq.gz![ref7]

\## Encountered 8047 unique sequences from 13338 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: R7\_sub\_R2\_filtered.fq.gz![ref7]

\## Encountered 4145 unique sequences from 7331 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: R8\_sub\_R2\_filtered.fq.gz![ref7]

\## Encountered 5696 unique sequences from 11192 total sequences read.![ref7]

\## Dereplicating sequence entries in Fastq file: R9\_sub\_R2\_filtered.fq.gz![ref7]

\## Encountered 3970 unique sequences from 7853 total sequences read. names(derep\_reverse) <- samples![ref7]![ref6]

Lets Infer ASVs

We start by pseudo-pooling our data and merging our forward and reverse reads to generate a complete amplicon. We have set the minimum BP overlap to be 170. We then generate a table count and examine it.

dada\_forward <- dada(derep\_forward, err=err\_forward\_reads, pool="pseudo")![ref6]

\## Sample 1 - 1498 reads in 552 unique sequences.![ref9]

\## Sample 2 - 529 reads in 224 unique sequences.

\## Sample 3 - 457 reads in 186 unique sequences.

\## Sample 4 - 475 reads in 203 unique sequences.

\## Sample 5 - 2109 reads in 723 unique sequences.

\## Sample 6 - 5527 reads in 2760 unique sequences. ## Sample 7 - 10354 reads in 5501 unique sequences. ## Sample 8 - 8028 reads in 3452 unique sequences. ## Sample 9 - 8138 reads in 4846 unique sequences. ## Sample 10 - 14423 reads in 9747 unique sequences. ## Sample 11 - 10906 reads in 6894 unique sequences. ## Sample 12 - 14672 reads in 9380 unique sequences.

\## Sample 13 - 15660 reads in 9310 unique sequences. ## Sample 14 - 15950 reads in 10225 unique sequences. ## Sample 15 - 17324 reads in 10407 unique sequences. ## Sample 16 - 16728 reads in 10981 unique sequences. ## Sample 17 - 13338 reads in 8520 unique sequences. ## Sample 18 - 7331 reads in 4432 unique sequences. ## Sample 19 - 11192 reads in 6171 unique sequences. ## Sample 20 - 7853 reads in 4234 unique sequences. ## 

\##    selfConsist step 2

dada\_reverse <- dada(derep\_reverse, err=err\_reverse\_reads, pool="pseudo")![ref6]

\## Sample 1 - 1498 reads in 514 unique sequences.![ref9]

\## Sample 2 - 529 reads in 203 unique sequences.

\## Sample 3 - 457 reads in 171 unique sequences.

\## Sample 4 - 475 reads in 189 unique sequences.

\## Sample 5 - 2109 reads in 660 unique sequences.

\## Sample 6 - 5527 reads in 2506 unique sequences.

\## Sample 7 - 10354 reads in 5054 unique sequences. ## Sample 8 - 8028 reads in 3113 unique sequences.

\## Sample 9 - 8138 reads in 4568 unique sequences.

\## Sample 10 - 14423 reads in 9288 unique sequences. ## Sample 11 - 10906 reads in 6445 unique sequences. ## Sample 12 - 14672 reads in 8799 unique sequences. ## Sample 13 - 15660 reads in 8762 unique sequences. ## Sample 14 - 15950 reads in 9710 unique sequences. ## Sample 15 - 17324 reads in 9613 unique sequences. ## Sample 16 - 16728 reads in 10432 unique sequences. ## Sample 17 - 13338 reads in 8047 unique sequences. ## Sample 18 - 7331 reads in 4145 unique sequences. ## Sample 19 - 11192 reads in 5696 unique sequences. ## Sample 20 - 7853 reads in 3970 unique sequences. ## 

\##    selfConsist step 2

merged\_amplicons <- mergePairs(dada\_forward, derep\_forward, dada\_reverse,![ref10]

`                               `derep\_reverse, trimOverhang=TRUE, minOverlap=170) seqtab <- makeSequenceTable(merged\_amplicons)

class(seqtab) 

\## [1] "matrix" dim(seqtab)![ref7]![ref6]

\## [1]   20 2521![ref7]

Now we identify any chimera sequences in our data.

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T)![ref6]

\## Identified 17 bimeras out of 2521 input sequences. sum(seqtab.nochim)/sum(seqtab) ![ref7]#Check to see how much we lost in terms of abundance![ref6]

\## [1] 0.9931372![ref7]

Now we check to see how many reads have been dropped throughout, and we save the table.

getN <- **function**(x) sum(getUniques(x))![ref11]

summary\_tab <- data.frame(row.names=samples, dada2\_input=filtered\_out[,1],

`                          `filtered=filtered\_out[,2], dada\_f=sapply(dada\_forward, getN),

`                          `dada\_r=sapply(dada\_reverse, getN), merged=sapply(merged\_amplicons, getN),

`                          `nonchim=rowSums(seqtab.nochim),

`                          `final\_perc\_reads\_retained=round(rowSums(seqtab.nochim)/filtered\_out[,1]\*100, 1)) summary\_tab

\##       dada2\_input filtered dada\_f dada\_r merged nonchim ## B1           1613     1498   1458   1466   1457    1457 ## B2            591      529    523    524    523     523 ## B3            503      457    450    451    450     450 ## B4            507      475    440    447    439     439 ## BW1          2294     2109   2066   2082   2054    2054 ## BW2          6017     5527   5134   5229   4716    4716 ## R10         11258    10354   9658   9819   9009    8847 ## R11BF        8627     8028   7544   7640   7150    6960 ## R11          8927     8138   7279   7511   6694    6577 ## R12         15681    14423  12420  12932  10714   10649 ## R1A         12108    10906   9584   9897   8559    8535 ## R1B         16091    14672  12937  13389  11202   11158 ## R2          17196    15660  14039  14498  12494   12436 ## R3          17494    15950  14210  14662  12503   12444 ## R4          18967    17324  16241  16501  14816   14750 ## R5          18209    16728  14800  15332  12905   12818 ## R6          14600    13338  11934  12311  10459   10448 ## R7           8003     7331   6515   6726   5630    5618 ## R8          12211    11192  10286  10513   9530    9454 ## R9           8600     7853   7215   7390   6740    6695 ##       final\_perc\_reads\_retained![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.025.png)

\## B1                         90.3

\## B2                         88.5

\## B3                         89.5

\## B4                         86.6

\## BW1                        89.5

\## BW2                        78.4

\## R10                        78.6

\## R11BF                      80.7

\## R11                        73.7

\## R12                        67.9

\## R1A                        70.5

\## R1B                        69.3

\## R2                         72.3

\## R3                         71.1

\## R4                         77.8

\## R5                         70.4

\## R6                         71.6

\## R7                         70.2

\## R8                         77.4

\## R9                         77.8

write.table(summary\_tab, "read-count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA) ![ref6]The following code block was skipped because it would be too time consuming to run.

\## downloading DECIPHER-formatted SILVA v138 reference ![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.026.png)download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA\_SSU\_r138\_2019.RData", destfile="S ILVA\_SSU\_r138\_2019.RData")

\## loading reference taxonomy object load("SILVA\_SSU\_r138\_2019.RData")

\## loading DECIPHER

**library**(DECIPHER) packageVersion("DECIPHER") 

\## creating DNAStringSet object of our ASVs

dna <- DNAStringSet(getSequences(seqtab.nochim))

\## and classifying

tax\_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)

Instead, we load the output taxonomy object that would have been created.

load("tax-info.RData") ![ref6]Now we generate a fasta file, count table, and taxonomy table from our DADA2 objects.

- giving our seq headers more manageable names (ASV\_1, ASV\_2...) ![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.027.png)asv\_seqs <- colnames(seqtab.nochim)

  asv\_headers <- vector(dim(seqtab.nochim)[2], mode="character")

**for** (i **in** 1:dim(seqtab.nochim)[2]) {

`  `asv\_headers[i] <- paste(">ASV", i, sep="\_") }

- making and writing out a fasta of our final ASV seqs: asv\_fasta <- c(rbind(asv\_headers, asv\_seqs)) write(asv\_fasta, "ASVs.fa")
- count table:

  asv\_tab <- t(seqtab.nochim)

row.names(asv\_tab) <- sub(">", "", asv\_headers)

write.table(asv\_tab, "ASVs\_counts.tsv", sep="\t", quote=F, col.names=NA)

- tax table:
- creating table of taxonomy and setting any that are unclassified as "NA" ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") asv\_tax <- t(sapply(tax\_info, **function**(x) {

`  `m <- match(ranks, x$rank)

`  `taxa <- x$taxon[m]

`  `taxa[startsWith(taxa, "unclassified\_")] <- NA

`  `taxa

}))

colnames(asv\_tax) <- ranks

rownames(asv\_tax) <- gsub(pattern=">", replacement="", x=asv\_headers) write.table(asv\_tax, "ASVs\_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

Remove contaminant sequences

**library**(decontam) packageVersion(![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.028.png)"decontam") 

packageVersion("decontam") ## [1] '1.6.0'![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.029.png)![ref3]

- identify contaminants ![ref1]

colnames(asv\_tab) # our blanks are the first 4 of 20 samples in this case

\##  [1] "B1"    "B2"    "B3"    "B4"    "BW1"   "BW2"   "R10"   "R11BF" "R11"  ## [10] "R12"   "R1A"   "R1B"   "R2"    "R3"    "R4"    "R5"    "R6"    "R7"   ## [19] "R8"    "R9"![ref12]

vector\_for\_decontam <- c(rep(TRUE, 4), rep(FALSE, 16)) contam\_df <- isContaminant(t(asv\_tab), neg=vector\_for\_decontam) table(contam\_df$contaminant) ![ref5]

\## ![ref12]

\## FALSE  TRUE ##  2498     6

contam\_asvs <- row.names(contam\_df[contam\_df$contaminant == TRUE, ])![ref5]

- check out the contaminants

asv\_tax[row.names(asv\_tax) %**in**% contam\_asvs, ]

\##         domain     phylum             class                 order              ## ASV\_104 "Bacteria" "Proteobacteria"   "Gammaproteobacteria" NA                 ## ASV\_219 "Bacteria" "Proteobacteria"   "Gammaproteobacteria" "Enterobacterales" ## ASV\_230 "Bacteria" "Proteobacteria"   "Gammaproteobacteria" NA                 ## ASV\_274 "Bacteria" "Proteobacteria"   "Gammaproteobacteria" "Pseudomonadales"  ## ASV\_285 "Bacteria" "Proteobacteria"   "Gammaproteobacteria" "Burkholderiales"  ## ASV\_623 "Bacteria" "Actinobacteriota" "Actinobacteria"      "Corynebacteriales" ##         family               genus             species![ref13]

\## ASV\_104 NA                   NA                NA     

\## ASV\_219 NA                   NA                NA     

\## ASV\_230 NA                   NA                NA     ## ASV\_274 "Pseudomonadaceae"   "Pseudomonas"     NA     ## ASV\_285 "Comamonadaceae"     "Tepidimonas"     NA     ## ASV\_623 "Corynebacteriaceae" "Corynebacterium" NA

Now lets remove contaminants from outputs and create some new files

- making new fasta file![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.032.png)

contam\_indices <- which(asv\_fasta %**in**% paste0(">", contam\_asvs)) dont\_want <- sort(c(contam\_indices, contam\_indices + 1)) asv\_fasta\_no\_contam <- asv\_fasta[- dont\_want]

- making new count table

asv\_tab\_no\_contam <- asv\_tab[!row.names(asv\_tab) %**in**% contam\_asvs, ]

- making new taxonomy table

asv\_tax\_no\_contam <- asv\_tax[!row.names(asv\_tax) %**in**% contam\_asvs, ]

\## and now writing them out to files write(asv\_fasta\_no\_contam, "ASVs-no-contam.fa") write.table(asv\_tab\_no\_contam, "ASVs\_counts-no-contam.tsv",

`            `sep="\t", quote=F, col.names=NA) write.table(asv\_tax\_no\_contam, "ASVs\_taxonomy-no-contam.tsv",             sep="\t", quote=F, col.names=NA)

Analysis

**library**(tidyverse) ; packageVersion("tidyverse")![ref2]

\## Warning: replacing previous import 'lifecycle::last\_warnings' by ## 'rlang::last\_warnings' when loading 'hms'![ref4]

\## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──![ref3]

\## ✓ ggplot2 3.3.3     ✓ purrr   0.3.4 ## ✓ tibble  3.1.2     ✓ dplyr   1.0.6 ## ✓ tidyr   1.1.3     ✓ stringr 1.4.0 ## ✓ readr   1.4.0     ✓ forcats 0.5.1![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.033.png)

\## ── Conflicts ────────────────────────────────────────── tidyverse\_conflicts() ── ## x dplyr::filter() masks stats::filter()![ref12]

\## x dplyr::lag()    masks stats::lag()

**library**(phyloseq) ; packageVersion("phyloseq") ![ref1]

**library**(vegan) ; packageVersion("vegan") 

\## Loading required package: permute![ref3]

\## Loading required package: lattice![ref3]

\## This is vegan 2.5-6![ref3]

**library**(DESeq2) ; packageVersion("DESeq2") ## Loading required package: S4Vectors![ref2]![ref3]

\## Loading required package: stats4![ref3]

\## Loading required package: BiocGenerics ## Loading required package: parallel![ref3]![ref3]

\## ![ref4]

\## Attaching package: 'BiocGenerics'

\## The following objects are masked from 'package:parallel':![ref14]

\## 

\##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ, ##     clusterExport, clusterMap, parApply, parCapply, parLapply, ##     parLapplyLB, parRapply, parSapply, parSapplyLB

\## The following objects are masked from 'package:dplyr': ## ![ref12]

\##     combine, intersect, setdiff, union

\## The following objects are masked from 'package:stats': ## ![ref12]

\##     IQR, mad, sd, var, xtabs

\## The following objects are masked from 'package:base':![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.035.png)

\## 

\##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,

\##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep, ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,

\##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,

\##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,

\##     union, unique, unsplit, which, which.max, which.min

\## ![ref4]

\## Attaching package: 'S4Vectors'

\## The following objects are masked from 'package:dplyr': ## ![ref12]

\##     first, rename

\## The following object is masked from 'package:tidyr': ## ![ref12]

\##     expand

\## The following object is masked from 'package:base': ## ![ref12]

\##     expand.grid

\## Loading required package: IRanges![ref3]

\## ![ref4]

\## Attaching package: 'IRanges'

\## The following object is masked from 'package:phyloseq': ## ![ref12]

\##     distance

\## The following objects are masked from 'package:dplyr': ## ![ref12]

\##     collapse, desc, slice

\## The following object is masked from 'package:purrr': ## ![ref12]

\##     reduce

\## Loading required package: GenomicRanges![ref3]

\## Loading required package: GenomeInfoDb![ref3]

\## Loading required package: SummarizedExperiment ## Loading required package: Biobase![ref3]![ref3]

\## Welcome to Bioconductor![ref14]

\## 

\##     Vignettes contain introductory material; view with

\##     'browseVignettes()'. To cite Bioconductor, see

\##     'citation("Biobase")', and for packages 'citation("pkgname")'.

\## ![ref4]

\## Attaching package: 'Biobase'

\## The following object is masked from 'package:phyloseq': ## ![ref12]

\##     sampleNames

\## Loading required package: DelayedArray ## Loading required package: matrixStats![ref3]![ref3]

\## ![ref4]

\## Attaching package: 'matrixStats'

\## The following objects are masked from 'package:Biobase': ## ![ref12]

\##     anyMissing, rowMedians

\## The following object is masked from 'package:dplyr': ## ![ref12]

\##     count

\## Loading required package: BiocParallel![ref3]

\## ![ref4]

\## Attaching package: 'DelayedArray'

\## The following objects are masked from 'package:matrixStats': ## ![ref12]

\##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

\## The following object is masked from 'package:purrr': ## ![ref12]

\##     simplify

\## The following objects are masked from 'package:base': ## ![ref12]

\##     aperm, apply, rowsum

**library**(dendextend) ; packageVersion("dendextend") ![ref6]

\## Registered S3 method overwritten by 'dendextend': ##   method     from ![ref12]

\##   rev.hclust vegan

\## ![ref13]

\## ---------------------

\## Welcome to dendextend version 1.17.1

\## Type citation('dendextend') for how to cite the package.

\## 

\## Type browseVignettes(package = 'dendextend') for the package vignette.

\## The github page is: https://github.com/talgalili/dendextend/

\## 

\## Suggestions and bug-reports can be submitted at: https://github.com/talgalili/dendextend/issues ## You may ask questions at stackoverflow, use the r and dendextend tags: 

\##   https://stackoverflow.com/questions/tagged/dendextend

\## 

\##  To suppress this message use:  suppressPackageStartupMessages(library(dendextend))

\## ---------------------

\## ![ref8]

\## Attaching package: 'dendextend'

\## The following object is masked from 'package:permute': ## ![ref12]

\##     shuffle

\## The following object is masked from 'package:stats': ## ![ref12]

\##     cutree

**library**(viridis) ; packageVersion("viridis") ![ref6]

\## Loading required package: viridisLite ![ref7]Clearing blank space from table, then taking a look

count\_tab <- read.table("ASVs\_counts-no-contam.tsv", header=T, row.names=1,                         check.names=![ref15]F, sep="\t")[ , -c(1:4)]

- Taking a look at our tables

tax\_tab <- as.matrix(read.table("ASVs\_taxonomy-no-contam.tsv", header=T,                                 row.names=1, check.names=F, sep="\t"))

sample\_info\_tab <- read.table("sample\_info.tsv", header=T, row.names=1,                               check.names=F, sep="\t")

- and setting the color column to be of type "character", which helps later sample\_info\_tab$color <- as.character(sample\_info\_tab$color) sample\_info\_tab # to take a peek

  ##       temp    type      char      color ## BW1    2.0   water     water       blue ## BW2    2.0   water     water       blue ## R10   13.7    rock    glassy      black ## R11BF  7.3 biofilm   biofilm  darkgreen ## R11    7.3    rock    glassy      black ## R12     NA    rock   altered chocolate4 ## R1A    8.6    rock   altered chocolate4 ## R1B    8.6    rock   altered chocolate4 ## R2     8.6    rock   altered chocolate4 ## R3    12.7    rock   altered chocolate4 ## R4    12.7    rock   altered chocolate4 ## R5    12.7    rock   altered chocolate4 ## R6    12.7    rock   altered chocolate4 ## R7      NA    rock carbonate  darkkhaki ## R8    13.5    rock    glassy      black ## R9    13.7    rock    glassy      black![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.037.png)

Take a look at Beta Diversity

I start by normalizing for sample depth using a “variance stabilizing transformation.”

- first we need to make a DESeq2 object![ref11]

deseq\_counts <- DESeqDataSetFromMatrix(count\_tab, colData = sample\_info\_tab, design = ~type) deseq\_counts\_vst <- varianceStabilizingTransformation(deseq\_counts)

- and here is pulling out our transformed table

vst\_trans\_count\_tab <- assay(deseq\_counts\_vst)

- and calculating our Euclidean distance matrix

euc\_dist <- dist(t(vst\_trans\_count\_tab))

Let’s make a hierarchical cluster plot of our data using the euclidean distance matrix we just created.

euc\_clust <- hclust(euc\_dist, method="ward.D2")![ref16]

plot(euc\_clust) 

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.039.jpeg)

euc\_dend <- as.dendrogram(euc\_clust, hang=0.1)![ref17]

dend\_cols <- as.character(sample\_info\_tab$color[order.dendrogram(euc\_dend)]) labels\_colors(euc\_dend) <- dend\_cols

plot(euc\_dend, ylab="VST Euc. dist.")

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.041.jpeg)

From this plot, it is clear that the samples from the basalt rocks (black and brown labels) are distinct from the rest of the samples. The basalt rocks also have distinct subclusters, with the brown labels coming from rocks that had highly altered, bumpy outsides, whereas the black labels came from basalt that was smoother and glassier.

Ordination

Here, I am going to try to do some principle coordinates analysis.

vst\_count\_phy <- otu\_table(vst\_trans\_count\_tab, taxa\_are\_rows=T) sample\_info\_tab\_phy <- sample\_data(sample\_info\_tab)![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.042.png)

vst\_physeq <- phyloseq(vst\_count\_phy, sample\_info\_tab\_phy)

- Now we visualize the PCoA with phyloseq

vst\_pcoa <- ordinate(vst\_physeq, method="MDS", distance="euclidean")

eigen\_vals <- vst\_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separatin g apart the samples

plot\_ordination(vst\_physeq, vst\_pcoa, color="char") + 

`  `geom\_point(size=1) + labs(col="type") + 

`  `geom\_text(aes(label=rownames(sample\_info\_tab), hjust=0.3, vjust=-0.4)) + 

`  `coord\_fixed(sqrt(eigen\_vals[2]/eigen\_vals[1])) + ggtitle("PCoA") + 

`  `scale\_color\_manual(values=unique(sample\_info\_tab$color[order(sample\_info\_tab$char)])) +   theme\_bw() + theme(legend.position="none")

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.043.jpeg)

From this, it appears as though the rock samples (black and brown labels) are most similar to eachother, with exterior alteration of the rock correlating with their microbial communities.

Alpha Diversity

First, I will generate rarefaction curves with the data.

rarecurve(t(count\_tab), step=100, col=sample\_info\_tab$color, lwd=2, ylab="ASVs", label=F)![ref10]

- and adding a vertical line at the fewest seqs in any sample

abline(v=(min(rowSums(t(count\_tab)))))

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.044.jpeg)

This suggests that rock samples (brown and black lines) may host greater richness than the other samples (water and biofilm). The samples from the highly altered basalt (brown) may host more microbial communities than the smooth basalt black lines)

Next, I will plot Chao1 richness estimates and shannon diversity values.

count\_tab\_phy <- otu\_table(count\_tab, taxa\_are\_rows=T) tax\_tab\_phy <- tax\_table(tax\_tab)![ref18]

ASV\_physeq <- phyloseq(count\_tab\_phy, tax\_tab\_phy, sample\_info\_tab\_phy)

- and now we can call the plot\_richness() function on our phyloseq object

plot\_richness(ASV\_physeq, color="char", measures=c("Chao1", "Shannon")) + 

`    `scale\_color\_manual(values=unique(sample\_info\_tab$color[order(sample\_info\_tab$char)])) +

`    `theme\_bw() + theme(legend.title = element\_blank(), axis.text.x = element\_text(angle = 90, vjust = 0.5, hjust 

- 1))

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.046.jpeg)

Next, I will plot Chao1 richness estimates and shannon diversity values, while grouping by sample types

plot\_richness(ASV\_physeq, x="type", color="char", measures=c("Chao1", "Shannon")) + ![ref10]

`    `scale\_color\_manual(values=unique(sample\_info\_tab$color[order(sample\_info\_tab$char)])) +

`    `theme\_bw() + theme(legend.title = element\_blank(), axis.text.x = element\_text(angle = 90, vjust = 0.5, hjust 

- 1))

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.047.jpeg)

Again, it appears as though the basalt hosts the most rich and diverse communities. Taxonomic summaries

Here, I will take a look at what microbial species are present in each sample.

- using phyloseq to make a count table that has summed all ASVs![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.048.png)
  - that were in the same phylum

phyla\_counts\_tab <- otu\_table(tax\_glom(ASV\_physeq, taxrank="phylum")) 

- making a vector of phyla names to set as row names

phyla\_tax\_vec <- as.vector(tax\_table(tax\_glom(ASV\_physeq, taxrank="phylum"))[,"phylum"]) rownames(phyla\_counts\_tab) <- as.vector(phyla\_tax\_vec)

- we also have to account for sequences that weren't assigned any
- taxonomy even at the phylum level 
- these came into R as 'NAs' in the taxonomy table, but their counts are
- still in the count table
- so we can get that value for each sample by subtracting the column sums
- of this new table (that has everything that had a phylum assigned to it)
- from the column sums of the starting count table (that has all
- representative sequences)

unclassified\_tax\_counts <- colSums(count\_tab) - colSums(phyla\_counts\_tab)

- and we'll add this row to our phylum count table:

phyla\_and\_unidentified\_counts\_tab <- rbind(phyla\_counts\_tab, "Unclassified"=unclassified\_tax\_counts)

- now we'll remove the Proteobacteria, so we can next add them back in
- broken down by class

temp\_major\_taxa\_counts\_tab <- phyla\_and\_unidentified\_counts\_tab[!row.names(phyla\_and\_unidentified\_counts\_tab) %**i n**% "Proteobacteria", ]

- making count table broken down by class (contains classes beyond the
- Proteobacteria too at this point)

class\_counts\_tab <- otu\_table(tax\_glom(ASV\_physeq, taxrank="class")) 

- making a table that holds the phylum and class level info

class\_tax\_phy\_tab <- tax\_table(tax\_glom(ASV\_physeq, taxrank="class")) 

phy\_tmp\_vec <- class\_tax\_phy\_tab[,2]

class\_tmp\_vec <- class\_tax\_phy\_tab[,3]

rows\_tmp <- row.names(class\_tax\_phy\_tab)

class\_tax\_tab <- data.frame("phylum"=phy\_tmp\_vec, "class"=class\_tmp\_vec, row.names = rows\_tmp)

- making a vector of just the Proteobacteria classes

proteo\_classes\_vec <- as.vector(class\_tax\_tab[class\_tax\_tab$phylum == "Proteobacteria", "class"])

- changing the row names like above so that they correspond to the taxonomy,
- rather than an ASV identifier

rownames(class\_counts\_tab) <- as.vector(class\_tax\_tab$class) 

- making a table of the counts of the Proteobacterial classes

proteo\_class\_counts\_tab <- class\_counts\_tab[row.names(class\_counts\_tab) %**in**% proteo\_classes\_vec, ] 

- there are also possibly some some sequences that were resolved to the level
- of Proteobacteria, but not any further, and therefore would be missing from
- our class table
- we can find the sum of them by subtracting the proteo class count table
- from just the Proteobacteria row from the original phylum-level count table

proteo\_no\_class\_annotated\_counts <- phyla\_and\_unidentified\_counts\_tab[row.names(phyla\_and\_unidentified\_counts\_ta b) %**in**% "Proteobacteria", ] - colSums(proteo\_class\_counts\_tab)

- now combining the tables:

major\_taxa\_counts\_tab <- rbind(temp\_major\_taxa\_counts\_tab, proteo\_class\_counts\_tab, "Unresolved\_Proteobacteria"=p roteo\_no\_class\_annotated\_counts)

- and to check we didn't miss any other sequences, we can compare the column
- sums to see if they are the same
- if "TRUE", we know nothing fell through the cracks

identical(colSums(major\_taxa\_counts\_tab), colSums(count\_tab)) 

\## [1] TRUE![ref7]

- now we'll generate a proportions table for summarizing:![ref17]

major\_taxa\_proportions\_tab <- apply(major\_taxa\_counts\_tab, 2, **function**(x) x/sum(x)\*100)

- if we check the dimensions of this table at this point

dim(major\_taxa\_proportions\_tab)

\## [1] 42 16![ref7]

- we see there are currently 42 rows, which might be a little busy for a![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.049.png)
- summary figure
- many of these taxa make up a very small percentage, so we're going to
- filter some out
- this is a completely arbitrary decision solely to ease visualization and
- intepretation, entirely up to your data and you
- here, we'll only keep rows (taxa) that make up greater than 5% in any
- individual sample

temp\_filt\_major\_taxa\_proportions\_tab <- data.frame(major\_taxa\_proportions\_tab[apply(major\_taxa\_proportions\_tab, 1, max) > 5, ])

- checking how many we have that were above this threshold

dim(temp\_filt\_major\_taxa\_proportions\_tab) 

\## [1] 12 16![ref7]

- now we have 12, much more manageable for an overview figure![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.050.png)
- though each of the filtered taxa made up less than 5% alone, together they
- may add up and should still be included in the overall summary
- so we're going to add a row called "Other" that keeps track of how much we
- filtered out (which will also keep our totals at 100%)

filtered\_proportions <- colSums(major\_taxa\_proportions\_tab) - colSums(temp\_filt\_major\_taxa\_proportions\_tab) filt\_major\_taxa\_proportions\_tab <- rbind(temp\_filt\_major\_taxa\_proportions\_tab, "Other"=filtered\_proportions)

Now, I will make some figures with the summary tables I created above.

- first let's make a copy of our table that's safe for manipulating![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.051.png)

filt\_major\_taxa\_proportions\_tab\_for\_plot <- filt\_major\_taxa\_proportions\_tab

- and add a column of the taxa names so that it is within the table, rather
- than just as row names (this makes working with ggplot easier)

filt\_major\_taxa\_proportions\_tab\_for\_plot$Major\_Taxa <- row.names(filt\_major\_taxa\_proportions\_tab\_for\_plot)

- now we'll transform the table into narrow, or long, format (also makes
- plotting easier)

filt\_major\_taxa\_proportions\_tab\_for\_plot.g <- pivot\_longer(filt\_major\_taxa\_proportions\_tab\_for\_plot, !Major\_Taxa, names\_to = "Sample", values\_to = "Proportion") %>% data.frame()

- take a look at the new table and compare it with the old one

head(filt\_major\_taxa\_proportions\_tab\_for\_plot.g)

\##     Major\_Taxa Sample Proportion ## 1 Nitrospirota    BW1  0.0000000 ## 2 Nitrospirota    BW2  0.0000000 ## 3 Nitrospirota    R10  9.5340421 ## 4 Nitrospirota  R11BF  0.1724138 ## 5 Nitrospirota    R11  6.3098677 ## 6 Nitrospirota    R12  3.7576327![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.052.png)

head(filt\_major\_taxa\_proportions\_tab\_for\_plot)![ref6]

\##                        BW1       BW2      R10      R11BF      R11         R12 ## Nitrospirota      0.000000 0.0000000 9.534042 0.17241379 6.309868  3.75763269 ## Crenarchaeota     5.851344 9.4479830 4.603031 6.06321839 4.926258 15.92296853 ## Bacteroidota     13.125988 0.7430998 1.696449 1.49425287 0.790634  1.67214655 ## Acidobacteriota   6.958355 0.8492569 3.890523 0.05747126 2.113426  6.20009394 ## Desulfobacterota 11.386400 0.0000000 0.000000 0.00000000 0.091227  0.02818225 ## Thermoplasmatota  8.961518 5.7749469 0.000000 0.00000000 0.000000  0.00000000 ##                        R1A       R1B        R2         R3          R4![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.053.png)

\## Nitrospirota      3.479789  3.100914  2.500804  2.8534684  2.75969623

\## Crenarchaeota    14.715876 13.676286 14.144419 15.1515152 12.05587198

\## Bacteroidota      1.909783  1.953755  1.141846  1.1253115  1.15269867

\## Acidobacteriota   4.674868  3.925435  4.575426  4.5333976  4.84811500

\## Desulfobacterota  0.000000  0.000000  0.000000  0.1848726  0.02712232

\## Thermoplasmatota  0.000000  0.000000  0.000000  0.0000000  0.00000000

\##                          R5        R6          R7         R8        R9

\## Nitrospirota     2.66032142  2.192227  1.21212121 3.59674178  3.584765

\## Crenarchaeota    8.66749883 12.387517 13.81461676 4.55939913 13.398058

\## Bacteroidota     1.34186301  1.129619  0.81996435 4.56997778  1.179985

\## Acidobacteriota  4.91496333  6.480950  4.52762923 5.93462393  2.793129

\## Desulfobacterota 0.08581682  0.000000  0.03565062 0.04231461  0.000000

\## Thermoplasmatota 0.00000000  0.000000  0.00000000 0.00000000  0.000000

\##                        Major\_Taxa

\## Nitrospirota         Nitrospirota

\## Crenarchaeota       Crenarchaeota

\## Bacteroidota         Bacteroidota

\## Acidobacteriota   Acidobacteriota

\## Desulfobacterota Desulfobacterota

\## Thermoplasmatota Thermoplasmatota

- now we want a table with "color" and "characteristics" of each sample to![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.054.png)
- merge into our plotting table so we can use that more easily in our plotting
- function
- here we're making a new table by pulling what we want from the sample
- information table

sample\_info\_for\_merge<-data.frame("Sample"=row.names(sample\_info\_tab), "char"=sample\_info\_tab$char, "color"=sampl e\_info\_tab$color, stringsAsFactors=F)

- and here we are merging this table with the plotting table we just made
- (this is an awesome function!)

filt\_major\_taxa\_proportions\_tab\_for\_plot.g2 <- merge(filt\_major\_taxa\_proportions\_tab\_for\_plot.g, sample\_info\_for\_ merge)

- and now we're ready to make some summary figures with our wonderfully
- constructed table

ggplot(filt\_major\_taxa\_proportions\_tab\_for\_plot.g2, aes(x=Sample, y=Proportion, fill=Major\_Taxa)) +     geom\_bar(width=0.6, stat="identity") +

`    `theme\_bw() +

`    `theme(axis.text.x=element\_text(angle=90, vjust=0.4, hjust=1), legend.title=element\_blank()) +

`    `labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples")

` `Now I will try to use boxplots to![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.055.jpeg)

visualize the taxa.

ggplot(filt\_major\_taxa\_proportions\_tab\_for\_plot.g2, aes(Major\_Taxa, Proportion)) +![ref11]![ref19]

`    `geom\_jitter(aes(color=factor(char), shape=factor(char)), size=2, width=0.15, height=0) +

`    `scale\_color\_manual(values=unique(filt\_major\_taxa\_proportions\_tab\_for\_plot.g2$color[order(filt\_major\_taxa\_prop ortions\_tab\_for\_plot.g2$char)])) +

`    `geom\_boxplot(fill=NA, outlier.color=NA) + theme\_bw() +

`    `theme(axis.text.x=element\_text(angle=45, hjust=1), legend.title=element\_blank()) +

`    `labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="All samples")

Now I will make similar plots but seperate water and rock samples, to try to get a cleaner look at things.![ref19]

- let's set some helpful variables first:![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.057.png)

bw\_sample\_IDs <- row.names(sample\_info\_tab)[sample\_info\_tab$type == "water"] rock\_sample\_IDs <- row.names(sample\_info\_tab)[sample\_info\_tab$type == "rock"]

- first we need to subset our plotting table to include just the rock samples to plot

filt\_major\_taxa\_proportions\_rocks\_only\_tab\_for\_plot.g <- filt\_major\_taxa\_proportions\_tab\_for\_plot.g2[filt\_major\_t axa\_proportions\_tab\_for\_plot.g2$Sample %**in**% rock\_sample\_IDs, ]

- and then just the water samples

filt\_major\_taxa\_proportions\_water\_samples\_only\_tab\_for\_plot.g <- filt\_major\_taxa\_proportions\_tab\_for\_plot.g2[filt \_major\_taxa\_proportions\_tab\_for\_plot.g2$Sample %**in**% bw\_sample\_IDs, ]

- and now we can use the same code as above just with whatever minor alterations we want
- rock samples

ggplot(filt\_major\_taxa\_proportions\_rocks\_only\_tab\_for\_plot.g, aes(Major\_Taxa, Proportion)) +

`    `scale\_y\_continuous(limits=c(0,50)) + # adding a setting for the y axis range so the rock and water plots are on the same scale

`    `geom\_jitter(aes(color=factor(char), shape=factor(char)), size=2, width=0.15, height=0) +

`    `scale\_color\_manual(values=unique(filt\_major\_taxa\_proportions\_rocks\_only\_tab\_for\_plot.g$color[order(filt\_major \_taxa\_proportions\_rocks\_only\_tab\_for\_plot.g$char)])) +

`    `geom\_boxplot(fill=NA, outlier.color=NA) + theme\_bw() +

`    `theme(axis.text.x=element\_text(angle=45, hjust=1), legend.position="top", legend.title=element\_blank()) + # m oved legend to top 

`    `labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Rock samples only")

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.058.jpeg)

- water samples![ref18]

ggplot(filt\_major\_taxa\_proportions\_water\_samples\_only\_tab\_for\_plot.g, aes(Major\_Taxa, Proportion)) +

`    `scale\_y\_continuous(limits=c(0,50)) + # adding a setting for the y axis range so the rock and water plots are on the same scale

`    `geom\_jitter(aes(color=factor(char)), size=2, width=0.15, height=0) +

`    `scale\_color\_manual(values=unique(filt\_major\_taxa\_proportions\_water\_samples\_only\_tab\_for\_plot.g$color[order(fi lt\_major\_taxa\_proportions\_water\_samples\_only\_tab\_for\_plot.g$char)])) +

`    `geom\_boxplot(fill=NA, outlier.color=NA) + theme\_bw() +

`    `theme(axis.text.x=element\_text(angle=45, hjust=1), legend.position="none") +

`    `labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Bottom-water samples only")

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.059.jpeg)

Now I will make some pie charts for our final taxonomic summary.

rock\_sample\_major\_taxa\_proportion\_tab <- filt\_major\_taxa\_proportions\_rocks\_only\_tab\_for\_plot.g[, c(1:3)] %>% pivo t\_wider(names\_from = Major\_Taxa, values\_from = Proportion) %>% column\_to\_rownames(![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.060.png)"Sample") %>% t() %>% data.fram e()

water\_sample\_major\_taxa\_proportion\_tab <- filt\_major\_taxa\_proportions\_water\_samples\_only\_tab\_for\_plot.g[, c(1:3)] %>% pivot\_wider(names\_from = Major\_Taxa, values\_from = Proportion) %>% column\_to\_rownames("Sample") %>% t() %>% d ata.frame()

- summing each taxa across all samples for both groups 

rock\_sample\_summed\_major\_taxa\_proportions\_vec <- rowSums(rock\_sample\_major\_taxa\_proportion\_tab) water\_sample\_summed\_major\_taxa\_proportions\_vec <- rowSums(water\_sample\_major\_taxa\_proportion\_tab)

rock\_sample\_major\_taxa\_summary\_tab <- data.frame("Major\_Taxa"=names(rock\_sample\_summed\_major\_taxa\_proportions\_ve c), "Proportion"=rock\_sample\_summed\_major\_taxa\_proportions\_vec, row.names=NULL) water\_sample\_major\_taxa\_summary\_tab <- data.frame("Major\_Taxa"=names(water\_sample\_summed\_major\_taxa\_proportions\_v ec), "Proportion"=water\_sample\_summed\_major\_taxa\_proportions\_vec, row.names=NULL)

- plotting just rocks

ggplot(data.frame(rock\_sample\_major\_taxa\_summary\_tab), aes(x="Rock samples", y=Proportion, fill=Major\_Taxa)) +     geom\_bar(width=1, stat="identity") +

`    `coord\_polar("y") +

`    `scale\_fill\_viridis(discrete=TRUE) +

`    `ggtitle("Rock samples only") +

`    `theme\_void() +

`    `theme(plot.title = element\_text(hjust=0.5), legend.title=element\_blank())

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.061.jpeg)

- and plotting just water samples![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.062.png)

ggplot(data.frame(water\_sample\_major\_taxa\_summary\_tab), aes(x="Bottom water samples", y=Proportion, fill=Major\_Ta xa)) + 

`    `geom\_bar(width=1, stat="identity") +

`    `coord\_polar("y") +

`    `scale\_fill\_viridis(discrete=TRUE) +

`    `ggtitle("Water samples only") +

`    `theme\_void() +

`    `theme(plot.title = element\_text(hjust=0.5), legend.title=element\_blank())

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.063.jpeg)

Statistical analysis

Do the microbial communities significantly differ between the different types of basalt? First: I use the betadisper test to check if there is enough homogeneity of dispersion between the two groups to perform an adonis test.

basalt\_sample\_IDs <- rock\_sample\_IDs[!rock\_sample\_IDs %**in**% "R7"]![ref18]

- new distance matrix of only basalts

basalt\_euc\_dist <- dist(t(vst\_trans\_count\_tab[ , colnames(vst\_trans\_count\_tab) %**in**% basalt\_sample\_IDs]))

- and now making a sample info table with just the basalts

basalt\_sample\_info\_tab <- sample\_info\_tab[row.names(sample\_info\_tab) %**in**% basalt\_sample\_IDs, ]

- running betadisper on just these based on level of alteration as shown in the images above:

anova(betadisper(basalt\_euc\_dist, basalt\_sample\_info\_tab$char)) 

\## Analysis of Variance Table![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.064.png)

\## 

\## Response: Distances

\##           Df Sum Sq Mean Sq F value Pr(>F) ## Groups     1   54.3   54.31  0.1244 0.7316 ## Residuals 10 4364.1  436.41

Now we can perform our adonis test, given that our value was .7.

adonis(basalt\_euc\_dist~basalt\_sample\_info\_tab$char)![ref6]

\## ![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.065.png)

\## Call:

\## adonis(formula = basalt\_euc\_dist ~ basalt\_sample\_info\_tab$char) 

\## 

\## Permutation: free

\## Number of permutations: 999

\## 

\## Terms added sequentially (first to last)

\## 

\##                             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   ## basalt\_sample\_info\_tab$char  1     23678 23677.8  2.4817 0.19883  0.003 \*\* ## Residuals                   10     95408  9540.8         0.80117          ## Total                       11    119086                 1.00000          ## ---

\## Signif. codes:  0 '\*\*\*' 0.001 '\*\*' 0.01 '\*' 0.05 '.' 0.1 ' ' 1

We get a significance level around .003, which gives us evidence that there is a signficant difference between the microbial communities between basalt types.

Now I make a new PCoA of the basalts, and include the signifcance level.

basalt\_vst\_count\_phy <- otu\_table(vst\_trans\_count\_tab[, colnames(vst\_trans\_count\_tab) %**in**% basalt\_sample\_IDs], ta xa\_are\_rows=![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.066.png)T)

basalt\_sample\_info\_tab\_phy <- sample\_data(basalt\_sample\_info\_tab)

basalt\_vst\_physeq <- phyloseq(basalt\_vst\_count\_phy, basalt\_sample\_info\_tab\_phy)

- generating and visualizing the PCoA with phyloseq

basalt\_vst\_pcoa <- ordinate(basalt\_vst\_physeq, method="MDS", distance="euclidean")

basalt\_eigen\_vals <- basalt\_vst\_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitud e of separating apart the samples

- and making our new ordination of just basalts with our adonis statistic

plot\_ordination(basalt\_vst\_physeq, basalt\_vst\_pcoa, color="char") + 

`    `labs(col="type") + geom\_point(size=1) + 

`    `geom\_text(aes(label=rownames(basalt\_sample\_info\_tab), hjust=0.3, vjust=-0.4)) + 

`    `annotate("text", x=25, y=68, label="Highly altered vs glassy") +

`    `annotate("text", x=25, y=62, label="Permutational ANOVA = 0.003") + 

`    `coord\_fixed(sqrt(basalt\_eigen\_vals[2]/basalt\_eigen\_vals[1])) + ggtitle("PCoA - basalts only") + 

`    `scale\_color\_manual(values=unique(basalt\_sample\_info\_tab$color[order(basalt\_sample\_info\_tab$char)])) +     theme\_bw() + theme(legend.position="none")

![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.067.jpeg)

DESeq2

I will use DEseq2 for differential abundance testing to test which representative sequences have different numbers of copy-counts between samples. This is going to help figure out which ASVs are contributing to that difference.

- first making a basalt-only phyloseq object of non-transformed values (as that is what DESeq2 operates on![ref20]

basalt\_count\_phy <- otu\_table(count\_tab[, colnames(count\_tab) %**in**% basalt\_sample\_IDs], taxa\_are\_rows=T) basalt\_count\_physeq <- phyloseq(basalt\_count\_phy, basalt\_sample\_info\_tab\_phy)

- now converting our phyloseq object to a deseq object

basalt\_deseq <- phyloseq\_to\_deseq2(basalt\_count\_physeq, ~char)

\## converting counts to integer mode![ref7]

- and running deseq standard analysis:![ref16]

basalt\_deseq <- DESeq(basalt\_deseq)

\## estimating size factors![ref7]

\## estimating dispersions![ref7]

\## gene-wise dispersion estimates ## mean-dispersion relationship ## final dispersion estimates![ref7]![ref7]![ref7]

\## fitting model and testing![ref7]

\## -- replacing outliers and refitting for 298 genes ## -- DESeq argument 'minReplicatesForReplace' = 7 ## -- original counts are preserved in counts(dds)![ref12]

\## estimating dispersions![ref7]

\## fitting model and testing ![ref7]Now, lets look at the results

- pulling out our results table, we specify the object, the p-value we are going to use to filter our results, ![ref20]

and what contrast we want to consider by first naming the column, then the two groups we care about deseq\_res\_altered\_vs\_glassy <- results(basalt\_deseq, alpha=0.01, contrast=c("char", "altered", "glassy"))

- we can get a glimpse at what this table currently holds with the summary command

summary(deseq\_res\_altered\_vs\_glassy) 

\## ![ref21]

\## out of 1762 with nonzero total read count

\## adjusted p-value < 0.01

\## LFC > 0 (up)       : 7, 0.4%

\## LFC < 0 (down)     : 8, 0.45%

\## outliers [1]       : 91, 5.2%

\## low counts [2]     : 1533, 87%

\## (mean count < 6)

\## [1] see 'cooksCutoff' argument of ?results

\## [2] see 'independentFiltering' argument of ?results

- this tells us out of ~1,800 ASVs, with adj-p < 0.01, there are 7 increased when comparing altered basalts t![ref15]
- glassy basalts, and about 6 decreased
  - "decreased" in this case means at a lower count abundance in the altered basalts than in the glassy basalt

s, and "increased" means greater proportion in altered than in glassy

- remember, this is done with a drastically reduced dataset, which is hindering the capabilities here quite a 

bit i'm sure

- let's subset this table to only include these that pass our specified significance level

sigtab\_res\_deseq\_altered\_vs\_glassy <- deseq\_res\_altered\_vs\_glassy[which(deseq\_res\_altered\_vs\_glassy$padj < 0.01), ]

- now we can see this table only contains those we consider significantly differentially abundant

summary(sigtab\_res\_deseq\_altered\_vs\_glassy) 

\## ![ref21]

\## out of 15 with nonzero total read count

\## adjusted p-value < 0.01

\## LFC > 0 (up)       : 7, 47%

\## LFC < 0 (down)     : 8, 53%

\## outliers [1]       : 0, 0%

\## low counts [2]     : 0, 0%

\## (mean count < 6)

\## [1] see 'cooksCutoff' argument of ?results

\## [2] see 'independentFiltering' argument of ?results

- next let's stitch that together with these ASV's taxonomic annotations for a quick look at both together![ref20]

sigtab\_deseq\_altered\_vs\_glassy\_with\_tax <- cbind(as(sigtab\_res\_deseq\_altered\_vs\_glassy, "data.frame"), as(tax\_tab le(ASV\_physeq)[row.names(sigtab\_res\_deseq\_altered\_vs\_glassy), ], "matrix"))

- and now let's sort that table by the baseMean column

sigtab\_deseq\_altered\_vs\_glassy\_with\_tax[order(sigtab\_deseq\_altered\_vs\_glassy\_with\_tax$baseMean, decreasing=T), ]

\##          baseMean log2FoldChange    lfcSE      stat       pvalue         padj ## ASV\_89  40.473687      -9.793325 1.932354 -5.068080 4.018490e-07 0.0001322083 ## ASV\_94  35.083004      -6.897441 1.853239 -3.721831 1.977837e-04 0.0059685900 ## ASV\_97  34.313175      -5.324613 1.469338 -3.623817 2.902866e-04 0.0079586903 ## ASV\_78  26.343239       7.248110 1.762674  4.111996 3.922527e-05 0.0032262783 ## ASV\_79  26.068816       7.233264 1.723258  4.197436 2.699538e-05 0.0029604930 ## ASV\_81  23.898472       6.381790 1.623360  3.931223 8.451476e-05 0.0046009326 ## ASV\_161 23.660885      -7.132292 1.846040 -3.863564 1.117447e-04 0.0046009326 ## ASV\_88  21.413396       6.218638 1.671868  3.719576 1.995577e-04 0.0059685900 ## ASV\_115 17.811784       6.684276 1.780090  3.755022 1.733259e-04 0.0059685900 ## ASV\_127 15.250724       6.461372 1.672511  3.863275 1.118768e-04 0.0046009326 ## ASV\_256 14.317145      -8.294524 1.965113 -4.220889 2.433411e-05 0.0029604930 ## ASV\_159 13.404936       6.274641 1.752095  3.581222 3.419915e-04 0.0086550145 ## ASV\_323 11.750176      -8.009515 2.274062 -3.522119 4.281120e-04 0.0099042348 ## ASV\_358 10.664628      -7.869840 1.967962 -3.998980 6.361598e-05 0.0041859317 ## ASV\_383  9.509639      -5.435301 1.549420 -3.507958 4.515609e-04 0.0099042348 ##           domain         phylum               class              order![](Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.070.png)

\## ASV\_89      <NA>           <NA>                <NA>               <NA>

\## ASV\_94  Bacteria Proteobacteria Alphaproteobacteria        Rhizobiales

\## ASV\_97  Bacteria Proteobacteria Alphaproteobacteria        Rhizobiales

\## ASV\_78  Bacteria Proteobacteria Gammaproteobacteria               <NA>

\## ASV\_79   Archaea  Crenarchaeota     Nitrososphaeria   Nitrosopumilales

\## ASV\_81  Bacteria   Nitrospirota         Nitrospiria      Nitrospirales

\## ASV\_161 Bacteria Proteobacteria Alphaproteobacteria               <NA>

\## ASV\_88  Bacteria Proteobacteria Gammaproteobacteria               <NA>

\## ASV\_115 Bacteria Proteobacteria Alphaproteobacteria               <NA>

\## ASV\_127 Bacteria   Dadabacteria       Dadabacteriia    Dadabacteriales

\## ASV\_256 Bacteria Proteobacteria Gammaproteobacteria              BD7-8

\## ASV\_159 Bacteria   Dadabacteria       Dadabacteriia    Dadabacteriales

\## ASV\_323 Bacteria Proteobacteria Alphaproteobacteria Puniceispirillales

\## ASV\_358 Bacteria Proteobacteria Alphaproteobacteria    Caulobacterales

\## ASV\_383 Bacteria Proteobacteria Alphaproteobacteria        Rhizobiales

\##                    family      genus species

\## ASV\_89               <NA>       <NA>    <NA>

\## ASV\_94       Rhizobiaceae       <NA>    <NA>

\## ASV\_97       Rhizobiaceae       <NA>    <NA>

\## ASV\_78               <NA>       <NA>    <NA>

\## ASV\_79  Nitrosopumilaceae       <NA>    <NA>

\## ASV\_81     Nitrospiraceae Nitrospira    <NA>

\## ASV\_161              <NA>       <NA>    <NA>

\## ASV\_88               <NA>       <NA>    <NA>

\## ASV\_115              <NA>       <NA>    <NA>

\## ASV\_127              <NA>       <NA>    <NA>

\## ASV\_256              <NA>       <NA>    <NA>

\## ASV\_159              <NA>       <NA>    <NA>

\## ASV\_323       EF100-94H03       <NA>    <NA>

\## ASV\_358  Parvularculaceae       <NA>    <NA>

\## ASV\_383      Rhizobiaceae       <NA>    <NA>

- this puts a sequence derived from a Rhizobiales at the second to highest (first is unclassified) that was d![ref16]

etected in ~7 log2fold greater abundance in the glassy basalts than in the highly altered basalts

Several of these have the same designation, this could be because organisms may have multiple copies of the 16s rRNA gene, which might not be identical and could be altering our results.

[ref1]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.004.png
[ref2]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.005.png
[ref3]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.006.png
[ref4]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.007.png
[ref5]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.012.png
[ref6]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.016.png
[ref7]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.017.png
[ref8]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.018.png
[ref9]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.022.png
[ref10]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.023.png
[ref11]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.024.png
[ref12]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.030.png
[ref13]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.031.png
[ref14]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.034.png
[ref15]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.036.png
[ref16]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.038.png
[ref17]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.040.png
[ref18]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.045.png
[ref19]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.056.jpeg
[ref20]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.068.png
[ref21]: Aspose.Words.c9141775-1976-49d0-9350-ba035ab5f116.069.png
