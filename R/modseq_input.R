###############################################################################
###                             LOAD PACKAGES                               ###
###############################################################################
library(seqinr)
library(parallel)
library(Biostrings)
library(ShortRead)
library(ggplot2)
# NOTES: Install libraries grid and gridExtra - Required and loaded in function 
#        IlluminaStat.
#        Install library scales - Required and loaded in run 3 and 4.
#        Install library Rcpp   - Required and loaded in run 5.

###############################################################################
###                         USER-DEFINED OPTIONS                            ###
###############################################################################
run.info <- paste("ModSeq | ", format(Sys.time(),"%Y%m%d"), sep = "")

### RUNNING MODE
## Which option to run (0=OFF, 1=ON)
run <- rep(0,5)
run[1] <- 1  # Pre-processing (quality trimming) and plots.
run[2] <- 1  # Run paired-end read assembly.
run[3] <- 1  # Run pattern search / Read mapping.
run[4] <- 0  # Run analysis on the search (library composition, modular variants
             # abundacies)
run[5] <- 0  # Run variant calling (mismatches and short indels)

####################### SEQUENCING MODE AND INPUT FILES #######################
# seq.mode:         sequencing mode (Options: SE, for single-read sequencing, 
#                   and PE, for paired-end sequencing. Default: PE).
# in.seqDir:        input directory (path to input sequencing data, by default 
#                   ./data/raw/).
# in.filename:      name of the file containing sequencing reads. Use when 
#                   option seq.mode is set to "SE".
# forward.filename: name of the file containing the set of forward reads without
#                   file-extension (format: fastq). 
# reverse.filename: name of the file containing the set of reverse reads without
#                   file-extension (format: fastq). 
seq.mode         <- "PE"  # Options: "PE" - paired-end reads, "SE" - single-end 
                          #          reads.
in.seqDir        <- "/Users/susanap/Documents/ModSeq/modseq/data/raw/"
in.filename      <- ""  # fastq format (for seq.mode = "SE").
forward.filename <- "4_GTGAA_L001_R1_001"  # fastq format (for seq.mode = "PE").
reverse.filename <- "4_GTGAA_L001_R2_001"  # fastq format (for seq.mode = "PE").

# in.modDir:    input directory (path to module-table file, by default 
#               ./data/modules/).
# mod.filename: name of module-table file without explicitly specifying 
#               the file extension (format: csv).
#               Table-header: names of modules, 
#               first row: name of variants or origin.
in.modDir    <- "/Users/susanap/Documents/ModSeq/modseq/data/modules/" 
mod.filename <- "Modules"  #.csv extension 

################################ OUTPUT FILES #################################
# out.dir:      output directory (path to output files, by default ./output).
# out.filename: prefix for output files (Default: "out").
# out.ssplot:   boolean variable indicating whether or not the summary- 
#               statistics plots should be output (Default: FALSE).
# out.varFiles: boolean variable indicating whether or not intermediary files
#               in the variant discovery module should be output (Default:FALSE)
out.dir      <- "/Users/susanap/Documents/ModSeq/modseq/output"
out.filename <- "4_GTGAA_L001" 
out.ssplot   <- FALSE
out.varFiles <- FALSE

########################## QUALITY TRIMMING options ###########################
# qtrim.thold: Phred value at or below which a nucleotide is removed.
#              It must be greater than 0. Default value: 10. 
#              The range of quality scores depend on the sequencing platform and
#              the base caller. In the case of current Illumina platforms, it 
#              ranges from 0 up to 41.
# qtrim.3end:  Trimming at the 3'-end only (0=OFF, 1=ON), by default 1. 
# qtrim.flag:  Indicated whether reads have been previously trimmed. Flag is 
#              set to 1 after execution of the quality-trimming step (Default
#              value: 1)
qtrim.thold <- 10
qtrim.3end  <- 1  # (0 = trimming at both ends, 1 = trimming at 3'-end only)
qtrim.flag  <- 1  # (0 = untrimmed reads, 1 = trimmed reads)

###################### PAIRED-END READ ASSEMBLY options #######################
# paired.flag:    Indicates which set of reads should be taken for the read 
#                 mapping step: paired-end (1) or single reads (0). By default, 
#                 it is set to 1, unless the sequencing mode is set to "SE". 
#                 Flag is changed to 1 after execution of the paired-end 
#                 assembly step. 
# paired.file:    If paired.flag is set to 0, select set of forward (f) or 
#                 reverse (r) reads. By default, it is set to "f". 
#                 Alternative - select "SE" sequencing mode.
# pandaseq.path:  Path to the pandaseq binary (if not in the PATH).
paired.flag   <- 1  # (0 = unpaired reads, 1 = paired reads)
paired.file   <- "f"  # Options: "f" - forward, "r" - reverse.
pandaseq.path <- character(0)

############################ READ MAPPING options #############################
# map.mode: Options bwa, gls - Grep-like search, and grPA - greedy search/
#           Pairwise alignement. By default bwa.
map.mode <- "bwa"

## Grep-like search options
# gls.ambiguity: Options - TRUE or FALSE. If TRUE (default), an ambiguous
#                character can be matched to any character on the alphabet 
#                associated to the code, according to the object being analysed 
#                (e.g. DNA, RNA, protein). Ambiguous matches enabled. If FALSE,
#                an ambiguous character can only be matched to the corresponding
#                ambiguous character, ambiguous matches disabled (woN in the 
#                naming conventions).
# gls.direction: Search direction ("f" - forward, "r" - reverse). By default "r".
# gls.mma:       Maximal mismatch allowance (*TODO*) 
#
# DNA case, an ambiguous character is any character which differs from those 
# contained in the DNA alphabet = {A,C,G,T}. 
# Example: alphabet DNAStringSet objects: "A" "C" "G" "T" "M" "R" "W" "S" "Y" 
# "K" "V" "H" "D" "B" "N" "-" "+" "."
# 
gls.ambiguity <- TRUE
gls.direction <- "r"
gls.mma       <- 0

## bwa mem
# bwa.path:       path to the bwa binaries (if not in the PATH).
# bwa.cVal:       discard a MEM if it has more than bwa.cVal occurrences in  
#                 the genome.
# bwa.dupl:       mark read duplicates.
# gatk.path:      path to jar file (if not in the PATH).
# picard.path:    path to jar file (if not in the PATH).
# samtools.path:  path to binary (if not in the PATH).
bwa.path      <- character(0)  #"/home/susanap/Documents/bwa/"
bwa.cVal      <- 20000
bwa.dupl      <- TRUE
gatk.path     <- "/Users/susanap/Documents/ModSeq/modseq/bin/GenomeAnalysisTK-3.6"
picard.path   <- "/Users/susanap/Documents/ModSeq/modseq/bin/picard-tools-2.4.1"
samtools.path <- character(0)

## Filtering (0,Inf)
# Filtering options should be greater than 0.
# coverage.left:  default 0.
# coverage.right: default 0.
# mapQ.thold:     default 8.
# editDist.thold: default 8.
coverage.left  <- 0  # leftclipped <= coverage.left
coverage.right <- 0  # rightclipped <= coverage.right
mapQ.thold     <- 8  # mapQ >= maQ.thold
editDist.thold <- 8  # editDistance <= editDist.thold

## Variants discovery
# mismatch.filter: Options - TRUE or FALSE. If TRUE (default) mismatches with 
#                  associated quality score < mismatch.qthold won't be taken 
#                  into account.
# mismatch.qthold: Phred value at or above which a mismatch is taken into 
#                  account. It must be greater than 0. Default value: 15.
# min.coverage:    Mininum number of reads per module combination required for 
#                  the variant calling (default 20).
# seq.error:       Expected sequencing error-rate (specific for the seq. 
#                  platform)
mismatch.filter <- TRUE  
mismatch.qthold <- 15  # quality >= mismatch.qthold
min.coverage    <- 20  
seq.error       <- 0.00304  

## MEMORY CONSUMPTION
## Tracing (FALSE=OFF, TRUE=ON), by default FALSE
mem.trace <- FALSE
