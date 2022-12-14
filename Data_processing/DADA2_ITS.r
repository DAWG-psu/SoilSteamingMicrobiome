library(dada2)
path_ITS <- "~/storage/work/tuc289/Soil_steaming/ITS"
list.files(path_ITS)

##2. Sort the forward and reverse reads
fnFs <- sort(list.files(path_ITS, pattern="_1.trimmedP.fastq.gz", full.name=TRUE))
fnRs <- sort(list.files(path_ITS, pattern="_2.trimmedP.fastq.gz", full.name=TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)
sample.names

FWD <- "CTTGGTCATTTAGAGGAAGTAA"
REV <- "CTGCGTTCTTCATCGAT"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path_ITS, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path_ITS, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

library(ShortRead)
library(Biostrings)
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#cutadapt <- "/storage/work/tuc289/miniconda3/bin/cutadapt"
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path_ITS, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_1.trimmedP.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.trimmedP.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

##3. Inspect the quality profiles
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

filtFs <- file.path(path.cut, "filtered", paste0(sample.names,
                                                 "_F_filt.fastq.gz"))
filtRs <- file.path(path.cut, "filtered", paste0(sample.names,
                                                 "_R_filt.fastq.gz"))


out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assigning taxonomy 

taxa2 <- assignTaxonomy(seqtab.nochim, "~/scratch/DAWG2022/3_database/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
taxa2 <- addSpecies(taxa, "~/scratch/DAWG2022/3_database/silva_species_assignment_v138.1.fa.gz")

taxa.print2 <- taxa2 # Removing sequence rownames for display only
rownames(taxa.print2) <- NULL
head(taxa.print2)

## Adding metadata and create phyloseq object
metadata <- read.table("~/scratch/DAWG2022/4_metadata/meta_ITS.csv", row.names=1, header=T, sep=",")
library(phyloseq)
ps_ITS <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(metadata), 
                   tax_table(taxa2))
dna <- Biostrings::DNAStringSet(taxa_names(ps_ITS))
names(dna) <- taxa_names(ps_ITS)
ps_ITS <- merge_phyloseq(ps_ITS, dna)
taxa_names(ps_ITS) <- paste0("ASV", seq(ntaxa(ps_ITS)))
ps_ITS
saveRDS(ps_ITS, "~/scratch/DAWG2022/5_phyloseq_results/ps_ITS.rds")
## Pipeline completed ##
