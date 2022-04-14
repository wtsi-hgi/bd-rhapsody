# BD Genomics Rhapsody Analysis Pipeline CWL

This repository contains the CWL and YML files necessary to run the BD Genomics Rhapsody Analysis Pipeline locally.

To obtain the files, click on _Downloads_ in the left navigation panel, and then _Download repository_.

For more details on how to use these files, see the BD Genomics Analysis Setup User Guide (Doc ID: 47383).

# Release notes

## v1.10.1 - April 14, 2022

- Fixed issue with cell label detection on reads from TCR/BCR, when TCR/BCR libraries were combined with other library types (WTA, Targeted, AbSeq) in a single sequencing index.
- Fixed issue with processing FASTQ files whose filenames end in fq.gz


## v1.10 - January 24, 2022

BD Rhapsody™ Targeted Analysis Pipeline and BD Rhapsody™ WTA Analysis Pipeline:

- Updated VDJ pipeline with improved performance, new assembly algorithm, new metrics and new output files containing all available contig sequences
- Added support for Rhapsody Enhanced Beads, with automatic bead version detection
- Added option to call putative cells based on AbSeq read counts (for troubleshooting only)
- Added “Alignment Categories” section to metric summary which provides a breakdown of alignments for read pairs with a valid cell label and UMI
- Added separate metric summary files for each sample tag for experiments using BD Single-Cell Multiplexing kits
- Renamed various metrics in outputs to reflect multiomics nature of data (Target Type → Bioproduct_Type, Gene/Target → Bioproduct)
- Added Pct_Read_Pair_Overlap and Median Reads Per cell metric to metric summary
- Improved support for larger runs on SBG
- Updated workflow on SBG to improve editing of resource requirements
- Optimized pipeline metadata handling
- Improved checking of reference files

## v1.9.1 - October 6, 2020

BD Rhapsody™ WTA Analysis Pipeline:

- Improved putative cell calling algorithm to reduce overcalling of putative cells in high cell input experiments
- Updated alignment settings to improve AbSeq mapping when R2 read length is greater than 75 bases


## v1.9 - July 29, 2020

BD Rhapsody™ Targeted Analysis Pipeline and BD Rhapsody™ WTA Analysis Pipeline:

- Improved FASTQ file pairing - filenames are flexible and pairing is now based on read sequence identifier
- Optimized pipeline in various steps for memory and storage usage
- Fixed bugs related to Sample Multiplexing Kit noise and DBEC mean molecule metric

BD Rhapsody™ Targeted Analysis Pipeline:

- Support for BD Rhapsody™ VDJ CDR3 protocol
- Read and molecule counts for targets from same gene symbol are combined in the output tables
- Updated Bowtie2 alignment parameters for improved sensitivity

BD Rhapsody™ WTA Analysis Pipeline:

- Updated Pct_Cellular Metrics calculations to match Bioinformatics handbook descriptions
- Added support for supplemental reference fasta files, which allow alignment to transgenes, like viral RNA or GFP
- Updated STAR alignment parameters for improved sensitivity


## v1.8 - Oct 4, 2019

BD Rhapsody™ Targeted Analysis Pipeline and BD Rhapsody™ WTA Analysis Pipeline:

- Added Sample_Tag_ReadsPerCell.csv to Multiplex Output
- Optimized pipeline in various steps for memory usage
- Fixed bug in status determination for UMI_Adjusted_Stats.csv file

BD Rhapsody™ Targeted Analysis Pipeline:

- Updated Targets section in Metrics_Summary.csv to calculate metrics based on targets detected in putative cells only
- Removed Clustering Analysis and outputs 

BD Rhapsody™ WTA Analysis Pipeline:

- Added support for BD™ AbSeq libraries
- Removed Targets section in Metrics_Summary.csv for WTA only libraries
- Removed Pct_Error_Reads and Error_Depth in UMI_Adjusted_Stats.csv, which are not applicable to WTA only libraries


## v1.7.1 - August 7, 2019
- Added BD Rhapsody™ WTA Analysis Pipeline
- Fixed bug that can cause stalling when zero putative cells were identified
- Fixed bug that affected runs using Disable Refined Putative Cell Calling option

## v1.6.1 - July 2, 2019
- Increased memory limits for GetDataTable and Metrics
- Fixed bug associated with "No Multiplex" option on SBG
- Uses fewer resources in AddToSam step.
 
## v1.6 - June 10, 2019
- Added new options for putative cell determination:
  - Exact Cell Count: Set a specific number of cells as putative, based on those with the highest error-corrected read count
  - Disable Refined Putative Cell Calling: Determine putative cells using only the basic algorithm
- Updated to Python 3
- Updated alignment defaults (minor molecule count changes expected)
- Local install only - CWL files are bundled into one file

## v1.5 - March 14, 2019
- Added support for BD Single-cell multiplexing kit – Mouse Immune
- Updated various filtering thresholds to support sequencing runs with shorter read length
- Deprecated pipeline input – BAM input
- Fixed bug in Quality Filter (minor metrics changes expected)
- Optimized pipeline (computationally faster, more scalable to support larger input data size, and better logging)

## v1.3 - July 31, 2018
- Added support for BD™ AbSeq assay
- Added support for BD™ single-cell multiplexing kit - Mouse Immune
- New pipeline input - AbSeq Reference
- New pipeline outputs - Unfiltered cell-gene data tables
- Updated Metrics_Summary.csv to support metrics from multiple sequencing libraries
- Updated Recursive Substitution Error Correction (RSEC) algorithm (minor molecule count changes expected)
- Optimized pipeline to run faster

## v1.02  - Nov 27, 2017
- Added support for BD Single-cell multiplexing kit - Human
- Improved pipeline speed by deleting large temp files
- Removed network requirement when running locally
- bug fix for the wrong docker image name - Dec 13, 2017
