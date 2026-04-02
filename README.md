# ECLIPSE: Enhanced CLonal Inference via Prediction of Single-cell Expression
This repo contains the code for the R package ECLIPSE used for single-cell TCR sequencing (scTCR-seq) analysis.

### Issues/Comments
Please reach out to Ethan Burns (ethan.burns [@] yale.edu) and David Braun (david.braun [@] yale.edu)

### Manuscript
ECLIPSE and its partner, VDJdive, are desribed in detail in our pre-print on bioRxiv: "VDJdive and ECLIPSE enhance single-cell TCR sequencing analysis through the probabilistic resolution of ambiguous clonotypes" available [here](https://doi.org/10.64898/2026.02.18.706444).





# How to use:

### Basics:
ECLIPSE contains 3 functions:
1. ECLIPSE: Used for running scTCR-seq analysis.
2. findThreeChainClones: A helper function that finds and notes T cell clones that posses 2 alpha or 2 beta chains.
3. tcrDoubletDetect: Based on the ouput of ECLIPSE, this function then can be used to remove clones that possess more chains than desired.








