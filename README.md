# ECLIPSE: Enhanced CLonal Inference via Prediction of Single-cell Expression
This repo contains the code for the R package ECLIPSE used for single-cell TCR sequencing (scTCR-seq) analysis.

### Issues/Comments
Please reach out to Ethan Burns (ethan.burns [@] yale.edu) and David Braun (david.braun [@] yale.edu).

### Manuscript
ECLIPSE and its partner, VDJdive, are described in detail in our pre-print on bioRxiv: "VDJdive and ECLIPSE enhance single-cell TCR sequencing analysis through the probabilistic resolution of ambiguous clonotypes" available [here](https://doi.org/10.64898/2026.02.18.706444).  





# How to use:

## Basics:
The ECLIPSE package contains 3 functions:
1. **ECLIPSE:** Used for running scTCR-seq analysis.
2. **findThreeChainClones:** A helper function that finds and notes T cell clones that possess 2 alpha or 2 beta chains.
3. **tcrDoubletDetect:** Based on the ouput of ECLIPSE, this function then can be used to remove cells that possess more chains than desired.


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
    - **A. Original barcode:** the barcode present in the TCR files and raw RNA files. Needs to start with the 16 nucleotide sequence and have -1 at the end with nothing before or behind either
    - **B. Index or sample number:** tells which TCR file in the vector of folders/files provided to ECLIPSE has the TCR information for each cell. This needs to be a number, i.e. 1 but not “1”,             and the order of files in the vector should match this number
    - **C. donor:** name that tells which donor mouse or human the cell is from. All TCR clones are called within donors, regardless of which treatment
    - **D. group:** tells how you want the final statistics on clone size and frequency to be calculated
        - Note this column isn’t always necessary. Often it makes sense to not contain this column and have the statistics also calculated on the donor column
4. Prepare a vector of either locations of directories that contain cellranger all or filtered contig_annotations.csv files or the locations of the files themselves. There are two options here:
    - **A. Vector of file locations:**
        - Good if all TCR files from multiple sequencing runs are stored in 1 folder, and/or if the files are not called "all_contig_annotations.csv" or "filtered_contig_annotations.csv"
        - An easy way to generate this is the function `list.files()`
        - Usage: feed the vector of file locations into file_paths, don't have anything for folders, and list `file_paths == "manual"`
    - **B. Vector of directory locations:**
        - Good if you have a bunch of folders each with 1 set of TCR files. Each directory must contain a file named "all_contig_annotations.csv" or                                                         "filtered_contig_annotations.csv". There cannot be more than 1 of the file that you want to work with. This is also compatabile with standard cellrange vdj outputs.
        - An easy way to generate this is the function `list.dirs()`
        - Usage: feed the vector of directory locations into folders, don't provide anything for file_paths, and list file_type as "all" or "filtered" depending on which file you would like.
            - We strongly recommend `file_type == "all"`, as ECLIPSE has customize filtering of TCR chains that retains biologically meaningful chains that are lost if file_type == "filtered"`                 is selected. This will also include other TCR chains that are not productive (i.e. they will not contribute to antigen recognition), but the inclusion of these chains greatly                     enhances the sensitivity of clonal tracking to allow for much better clonal prediction.
            - If you are planning on cloning TCR chains or are only interested in TCR chains that are productive, use `file_type == "filtered"`


### Arguments
- folders: 
- file_paths: vector of TCR file directories in quotes. The order should match the numbers provided in batch. These should end with “all_contig_annotations.csv”.
- file_type: "manual"
- seurat_object is self explanatory
- batch = name of column in Seurat meta.data. Should be in quotes, with the column containing numbers.
- donor: name of column in Seurat meta.data. Should be in quotes.
- group: name of column in Seurat meta.data. Should be in quotes, and you probably want to provide the same column as you did for donor (unless you want to the clone sizes/frequency calculated     by some other condition
- original_barcode: name of column in Seurat meta.data. Should be in quotes
- write_folder: Location of a directory where you want the temporary file written. Should be in quotes and not the same directory that has the original TCR files.
- format: Ignore this unsure you are using older TCR contig files that store blank data as "None" instead of "". If this is the case, list `format = "None"`.














