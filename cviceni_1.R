rm(list=ls())

#Path
setwd('V:/MPCPRG/exercise_02')

#instalation
if(!require("seqinr", quietly = TRUE)){
  install.packages("seqinr")
}
library("seqinr")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::available()

BiocManager::install("Biostrings")

library(Biostrings)
library(seqinr)

#sequences_seqinr <- read.fasta(" fishes.fna.gz")

# Open the .fna.gz file
file_connection <- gzfile("fishes.fna.gz", "rt")

# Read the file line by line
file_content <- readLines(file_connection)

file_content

#N
seq <- readDNAStringSet("fishes.fna.gz")
seq

# TASK 2
length(seq)
width(seq[1])
names(seq)
seq1 <- seq[1]
seq1_sequence <- seq[[1]]
seq1_string <- toString(seq[1])
#help(XStringSet)

# TASK 3
seq1 <- AAStringSet(as.character(seq[1]))
seq2 <- AAStringSet(as.character(seq[2]))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pwalign")

library(pwalign)

# Perform global alignment
alignment <- pairwiseAlignment(seq1, seq2, substitutionMatrix = "BLOSUM62", gapOpening = -1, gapExtension = 1, type = "global")



# TASK 4
names_list <- c("anna", "jana", "kamil", "norbert", "pavel",
                "petr", "stanislav", "zuzana")
grep("jana", names_list, perl = TRUE) #vrati index shody
grep("n+", names_list, perl = TRUE) #vsechny jmena obsahujici alespon jedno n
grep("n{2}", names_list, perl = TRUE) # obsahuje usek 'nn'
grep("^n", names_list, perl = TRUE) # zacinajici na n
grep("Anna|Jana", names_list, perl = TRUE) # nebo
grep("^z.*a$", names_list, perl = TRUE)

grep("^.ba.", "tba1aaa")


#TASK 5
mid <- DNAString("ACGAGTGCGT")
rev_mid <- as.character(reverseComplement(mid))
mid <- as.character(mid)

pattern <- paste0("(?=.*", mid, ")(?=.*", rev_mid, ")")
#grep("pattern", seq, perl = TRUE)

matched_seqs <- DNAStringSet()

for (i in seq_along(seq)) {
  current_seq <- as.character(seq[[i]])
  if (grepl(pattern, current_seq, perl = TRUE)) {
    matched_seqs <- append(matched_seqs, seq[i])
  }
}

# 4351
length(matched_seqs)

# TASK 6
mids <- read.csv('fishes_MIDs.csv', sep = ";")
mids


