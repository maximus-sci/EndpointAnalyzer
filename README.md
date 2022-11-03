# Introduction
Code to analyze and visualize the regions of in vitro generated nucleosome arrays

This code was used in the manuscript:
A pioneer factor locally opens compacted chromatin to enable targeted ATP-dependent nucleosome remodeling

#Preprocessing
Nuclease-digested nucleosomes are sequenced on a Nanopore flow cell to produce fast5 files. Once fast5 files are basecalled, demultiplexed, aligned and converted from bam to bed, they can be used as the input to the "EndpointAnalyzer" R script.

# Endpoint_Analyzer
This is an R script and is easiest to run within R Studio. Within the script, you must specify the location of the reference sequence for the nucleosome array, the location of the input bed files and, optionally, the output location for the filtered and normalized bedgraph files. The script will also generate a shell script for the conversion of the bedgraph files into bigWigs. Options within the script allow you to filter out reads based on length, quality, and location.

#Bedgraph_normalizer_relToCntrl
This python script should be placed in the directory containing normalized and filtered bedgraph files. The script identifies all bedgraphs in the folder and then allows you to specify which one is the control sample. The difference between the experimental and control samples is then output as a new set of bedgraph files.

#Bedgraph_smoother
This takes bedgraph files as input. It contains a series of ggplot graphing functions that can help demonstrate histone positions on the nucleosome array. Note that this script requires some knowledge of R and dplyr coding to change settings and correctly specify input files. 

