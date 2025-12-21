---
title: ‘ResMobiLys: a nextflow pipeline for Resistome and associated Mobilome anaLysis of field and clinical Eubacteria isolates’

tags:
    - Nextflow
    - Genomics
    - Antimicrobial resistance
    - Mobile genetic elements

authors:
    - name: Vahiniaina Andriamanga
      orcid: 0000-0001-6165-4858
      corresponding: true
      affiliation: 1,2
    - name: Patricia LICZNAR-FAJARDO
      orcid: 0000-0001-8725-2937
      affiliation: 3
      name: Estelle JUMAS-BILAK
      orcid: 0009-0004-9373-064X
      affiliation: 3
    - name: Stéphanie Bedhomme
      corresponding: true
      orcid: 0000-0001-8075-0968
      affiliation: 2  

affiliations:
    - name: HSM, Univ Montpellier, CNRS, IRD, Montpellier, France
      index: 1
    - name: CEFE, Univ Montpellier, CNRS, EPHE, IRD, Montpellier, France
      index: 2
    - name: HSM, Univ Montpellier, CNRS, IRD, CHU Montpellier, Montpellier, France
      index: 3
date: 19 December 2025
bibliography: paper.bib
---

# Summary  
Mobile Genetic Elements (MGEs) are segments of DNA that can move within or between genomes, enabling the transfer of genetic material between cells of the same species, or between different species. MGEs play a crucial role in microbial ecology because they can transfer genes with diverse functions across microbial populations [@ranking_what_2011]. Among them, AntiMicrobial Resistance (AMR) genes are responsible for one of the largest public health issues. A recent global study, compiling data from more than 200 countries, estimated that in 2019, approximately 5 million deaths were associated with AMR [@murray_global_2022]. Projections suggest that this could rise to 10 million deaths by 2050 [@thompson_staggering_2022]. Besides AMR, MGEs also carry virulence factors (VFs) and heavy metal resistance genes, which enhance bacterial pathogenicity and environmental adaptation [@morales_role_2023].
Advances and availability in high-throughput sequencing now enable large-scale sequencing of diverse bacterial isolates from human, animal, and environment sources, generating large datasets that require automated and integrated analyses to produce interpretable outputs.  
For this purpose, we developed ResMobiLys, an end-to-end workflow that performs de novo assembly and identififies MGEs, including plasmids, integrons, prophages, Integrative and Conjugative Elements (ICEs), and transposable elements. ResMobilYs integrates tools designated to operate accross the wide diversity of bacterial taxa and simultaneously detects antibiotic resistance genes (ARGs), virulence factor, and heavy metal resistance genes, providing an integrated analysis of the mobilome and associated resistome in field and Eubacteria analysis.

# Statement of need  
Antimicrobial resistance (AMR) remains a major public health threat. Mobile genetic elements (MGEs) play a critical role in recruiting and disseminating resistance genes to pathogenic bacteria [@pradier_ecology_2023]. Consequently, identifying and characterizing MGEs is essential for understanding AMR epidemiology. Althought advances of massive parallel sequencing has made genomic analysis widely accessible, user-friendly tools for comprehensively identifying the global mobilome and its association with AMR reamin limited. This lack of accessibile and integrated solutions represent a bottleneck preventing more researchers from conducting analyses to understand the role of mobilome in the propagation of resistance.
Conducting these analyses requires the integration of multiple specialized bioinformatics tools into a unified framework. To address this need, we developed ResMobiLys, an automated, all-in-one workflow for large-scale identification of MGEs and ARGs. ResMobiLys offers scalable, automated, portable, and reproducible solution for processing large genomics datasets and enables consistent comparisons across diverse bacterial samples.  

Few bacterial genomics pipelines currently support comprehensive resistome or mobilome analysis. Examples include MobileElementFinder [@johansson_detection_2020], MGEfinder [@durrant_bioinformatic_2020], and Baargin [@ayer_baargin_2023]. MobileElementFinder relies on homology-based detection using curated MGE databases, limiting its applicability to less-characterized species. MGEfinder is less dependent on genome annotations but is tailored for short-read sequencing and requires a suitable reference genome. Neither tool integrates ARGs detection nor supports comparative genome analysis. Baargin includes AMR detection and plasmid identification but is mainly oriented toward strain-level analyses rather than broader inter-species comparisons.   
ResMobiLys fills these gaps by enabling plasmids recovery from long-read assemblies and identifying MGE while minimizing reliance on reference genomes.Its core components include (i) detection and annotation of AMR genes from long-read assembly,(ii) identification of MGEs combined with plasmid clustering and comparative analysis, and (iii) integration of MGE–AMR associations to characterize resistance dissemination. This versatility allows application across diverse bacterial taxa, including isolates from unidentified lineages and provides a comprehensive framework for mobilome and resistome profiling.

# Materials and methods

## Features  
The pipeline is implemented with Nextflow [@di_tommaso_nextflow_2017] (https://github.com/NjivaAndriamanga/resmobilys), which organizes and automates the different steps of the workflow, ensuring that analyses are carried out in a consistent and reproducible way. It can run on a wide range of computing systems, including high-performance computing (HPC) clusters, and processes samples in parallel to efficiently handle large datasets. To guarantee reliable and consistent execution, the pipeline uses containers to bundle all software dependencies. This approach makes ResMobilYs easy to deploy, portable across different environments, and well-suited for both small-scale and large-scale analyses.

## The workflow  
**Input Data**  
ResMobiLys accepts either long reads alone or a combination of long and short reads as input. Long-read sequencing enables more complete and less fragmented assemblies, significantly improving plasmid resolution.  
**Genome assembly**  
Genome assembly is conducted with Hybracter [@bouras_hybracter_2024], which provides highly accurate bacterial genome assemblies and enables efficient plasmids recovery from de novo assemblies. Taxonomic classification is performed using Kraken2 [@wood_improved_2019]. Assembly quality is assessed using BUSCO [@manni_busco_2021].
**Antimicrobial Resistance genes, Virulence Factors, and heavy metal resistance**  
ARGs are detected using RGI with CARD database [@alcock_card_2023], which is particularly well suited for environmental AMR surveillance [@papp_review_2022]. Resistance determinants for metals and biocides are identified with ABRicate using the MEGARes database [@bonin_megares_2023], and virulence factors are detected using VFDB [@dong_expanded_2024].
**Mobile Genetic Elements annotation**  
An estimated chromosome size can be provided or inferred automatically with Hybracter. Small circular contigs below this threshold are classified as plasmids [@bouras_hybracter_2024], while small linear contigs are further analyzed using PLASMe [@tang_plasme_2023] or PLATON [@schwengers_platon_2020]. These tools use hybrid approaches combining sequence similarity with additional features or models, enabling accurate discrimination between plasmidic and chromosomal contigs and detection of highly divergent plasmids without relying on classical markers such as replicons or relaxases.  
Plasmid typing and clustering are performed with MOB-suite [@robertson_mob-suite_2018], with MOB-typer for plasmid typing and MOB-cluster to identify similar plasmids across isolates based on complete sequence comparisons. 
Integrons are detected using IntegronFinder  [@neron_integronfinder_2022], which detects integrase using HMM profiles, predicts attC sites with sequence and structural motifs. ICEs are annotated using the CONJScan module of MacSyFinder [@cury_integrative_2017], and insertion sequences (ISs) are identified with ISEScan [@xie_isescan_2017]. Prophages are predicted with DBSCAN-SWA [@gan_dbscan-swa_2022]. By prioritizing approaches that operate beyond direct pairwaise homology, the workflow enhances the detection of novel MGEs and expand mobilome exploration to uncharacterized elements.
**Output**  
ResMobiLys produces both isolate-specific and global results, summarizing antimicrobial resistance and mobility features.
At the isolate level, the pipeline identifies ARGs, their genomic locations (chromosomal or plasmid), and associated MGEs while retaining individual outputs from each integrated tool for further downstream analyses.
At the global level, it generates:
- A presence–absence matrix of ARGs across isolates, including genomic context and plasmid cluster information.
- A list of ARG–MGE associations, highlighting genes located within mobile elements.
- A plasmid cluster summary, grouping plasmids by sequence similarity and their associated ARGs.
