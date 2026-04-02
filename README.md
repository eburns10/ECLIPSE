# ECLIPSE: Enhanced CLonal Inference via Prediction of Single-cell Expression
This repo contains the code for the R package ECLIPSE used for single-cell TCR sequencing (scTCR-seq) analysis.

### Issues/Comments
Please reach out to Ethan Burns (ethan.burns [@] yale.edu) and David Braun (david.braun [@] yale.edu).

### Manuscript
ECLIPSE and its partner, VDJdive, are described in detail in our pre-print on bioRxiv: "VDJdive and ECLIPSE enhance single-cell TCR sequencing analysis through the probabilistic resolution of ambiguous clonotypes" available [here](https://doi.org/10.64898/2026.02.18.706444).  





# How to use:

## Basics:
The ECLIPSE package contains 3 functions:
1. ECLIPSE: Used for running scTCR-seq analysis.
2. findThreeChainClones: A helper function that finds and notes T cell clones that possess 2 alpha or 2 beta chains.
3. tcrDoubletDetect: Based on the ouput of ECLIPSE, this function then can be used to remove cells that possess more chains than desired.


## Setup
1. Install [VDJdive](https://bioconductor.org/packages/release/bioc/html/VDJdive.html).
2. Install [scRepertoire](https://www.bioconductor.org/packages/release/bioc/html/scRepertoire.html).
3. Move “ECLIPSE.R” file onto your computer or HPC
4. Load in the required packages and code
```
library(tidyverse)
library(Seurat)
library(scRepertoire)
library(VDJdive)
source("/path/on/computer/to/file/ECLIPSE.R")
```

## Prepping Input for ECLIPSE
1.	Have your scRNA-seq data as a Seurat object
2.	Make a folder on your computer or the HPC that can be used to store a temporary file. Don’t use the same folder that has all the raw TCR files
3.	Make the following columns in your Seurat metadata. They can be named whatever you want
    a.	Original barcode: the barcode present in the TCR files and raw RNA files. Needs to start with the 16 nucleotide sequence and have -1 at the end with nothing before or behind either
    b.	Index or sample number: tells which TCR file in the vector of folders/files provided to ECLIPSE has the TCR information for each cell. This needs to be a number, i.e. 1 but not “1”, and the         order of files in the vector should match this number
    c.	donor: name that tells which donor mouse or human the cell is from. All TCR clones are called within donors, regardless of which treatment
    d.	group: tells how you want the final statistics on clone size and frequency to be calculated
        i.	Note this column isn’t always necessary. Often it makes sense to not contain this column and have the statistics also calculated on the donor column


### Arguments
    - file_paths = vector of TCR file directories in quotes. The order should match the numbers provided in batch. These should end with “all_contig_annotations.csv”.
    - file_type = "manual"
    - seurat_object is self explanatory
    - batch = name of column in meta data that you made. Should be in quotes, with the column containing numbers.
    - donor = name of column in meta data that you made. Should be in quotes.
    - group = name of column in meta data that you made. Should be in quotes, and you probably want to provide the same column as you did for donor.
    - original_barcode = name of column in meta data that you made. Should be in quotes
    - write_folder = Directory where you want the temporary file written. Should be in quotes














