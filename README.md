# **ResMobiLys: A nextflow pipeline for Resistome and associated Mobilome anaLysis of field and Clinical Eubacteria isolates**

# INTRODUCTION
Mobile Genetic Elements (MGEs) are DNA segments capable of moving within or between genomes, facilitating the transfer of genetic material across
microbial populations. MGEs play a key role in microbial ecology by spreading genes linked to antimicrobial resistance (AMR), which poses significant
public health challenges. Advances in high-throughput sequencing technologies allow large-scale studies of diverse bacterial genomes from various origins, generating vast
data volumes that require efficient and integrated analysis.

# PURPOSE
ResMobiLys, a Nextflow pipeline designed for comprehensive mobilome and resistome analysis. ResMobiLys provides an end-to-end
solution, from de novo assembly to precise identification of MGEs like plasmids, integrons, prophages, and transposable elements. It also detects Antibiotic Resistance
Genes (ARGs) and virulence factors, whether associated with MGEs or not, enabling streamlined and detailed genomic feature analysis.

# Installation
Git : installation and git clone
Nextflow: install nextflow
Singularity: install singularity

# USAGE
index_file


test with conda : nextflow run waterisk -profile conda

test data: use profile -perso
the test data includes 9 samples: //describe

Tang X, Shang J, Ji Y, Sun Y. PLASMe: a tool to identify PLASMid contigs from short-read assemblies using transformer. Nucleic Acids Res. 2023 Aug 25;51(15):e83. doi: 10.1093/nar/gkad578. PMID: 37427782; PMCID: PMC10450166.

George Bouras, Ghais Houtak, Ryan R Wick, Vijini Mallawaarachchi, Michael J. Roach, Bhavya Papudeshi, Louise M Judd, Anna E Sheppard, Robert A Edwards, Sarah Vreugde - Hybracter: Enabling Scalable, Automated, Complete and Accurate Bacterial Genome Assemblies. (2024) Microbial Genomics doi: https://doi.org/10.1099/mgen.0.001244.
