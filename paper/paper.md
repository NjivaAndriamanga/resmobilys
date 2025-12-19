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
      affiliation: “1,2”
    - name: Patricia LICZNAR-FAJARDO
      orcid: 0000-0001-8725-2937
      affiliation: “3”
      name: Estelle JUMAS-BILAK
      orcid: 0009-0004-9373-064X
      affiliation: “3”
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

#Summary
Mobile Genetic Elements (MGEs) are segments of DNA that can move within or between genomes, enabling the transfer of genetic material within genomes, between cells of the same species, or between different species. MGEs play a crucial role in microbial ecology  because they can transfer genes with diverse functions across microbial populations [@ranking_what_2011]. Among them, AntiMicrobial Resistance (AMR) genes are responsible for one of the largest public health issues. A recent global study, compiling data from more than 200 countries, estimated that in 2019, approximately 5 million deaths were associated with AMR [@murray_global_2022]. Projections suggest that this could rise to 10 million deaths by 2050 [@thompson_staggering_2022]. Besides AMR, MGEs also carry virulence factors (VFs) and heavy metal resistance genes, which enhance bacterial pathogenicity and environmental adaptation [@morales_role_2023].
The evolution and availability of high-throughput sequencing technologies enable large-scale studies and the rapid sequencing of multiple bacterial isolates from human, animal, and environmental origins, likely representing a wide taxonomic diversity. High-throughput sequencing technologies produce large amounts of data, which require automatic, efficient, and integrated analysis to extract meaningful, synthetic, and interpretable information.
For this purpose, we developed a Nextflow [@di_tommaso_nextflow_2017] workflow called ResMobiLys. ResMobiLys performs an end-to-end analysis, starting with de novo assembly and progressing to the precise identification of MGEs, including plasmids, integrons, prophages, Integrative and Conjugative Elements (ICEs), and transposable elements. To do so, it integrates tools allowing analysis of the wide diversity of environmental bacterial taxa. The workflow also detects antibiotic resistance genes (ARGs), virulence factor genes, and heavy metal resistance genes. By covering the entire process—from initial genomic analysis to detailed MGE genomic features— ResMobiLys offers a robust, streamlined solution for comprehensive mobilome and associated resistome analysis of field and clinical Eubacteria isolates.

#Statement of need
Antimicrobial resistance (AMR) continues to pose a major public health threat. Mobile genetic elements (MGEs) play a critical role in recruiting and disseminating resistance genes to pathogenic bacteria [@pradier_ecology_2023]. Consequently, identifying and characterizing MGEs is essential for understanding AMR epidemiology. The rapid development of massively parallel sequencing has made genomic analysis widely accessible. However, there remains a scarcity of user-friendly tools for comprehensively identifying the global mobilome and its association with AMR. This lack of accessibility may be a bottleneck preventing more researchers from conducting such integrated analyses. Yet, understanding the association of the resistome with the mobilome is vital for elucidating the role of the mobilome in the propagation of resistance.
Conducting these analyses requires the integration of specific and diverse bioinformatics tools into a unified, cohesive pipeline. To address this need, we developed ResMobiLys —an automated, all-in-one workflow tailored for large-scale, high-throughput MGE identification and AMR analysis. ResMobiLys offers scalable, automated, portable, and reproducible workflows, ideally suited for processing extensive datasets and enabling consistent comparisons across diverse samples. 
Few bacterial genomics pipelines currently support comprehensive resistome or mobilome analysis. Examples include MobileElementFinder [@johansson_detection_2020], MGEfinder [@durrant_bioinformatic_2020], and Baargin [@ayer_baargin_2023]. MobileElementFinder employs a homology-based approach, relying on well-annotated genomes and databases of known MGEs, which limits its applicability to less-characterized species. MGEfinder is less dependent on genome annotations but is tailored for short-read sequencing and requires a carefully selected reference genome. Notably, neither tool integrates antimicrobial resistance (AMR) detection and enables comparative genome analysis. Baargin, while providing AMR detection and plasmid identification, is mainly oriented toward strain-level analyses rather than large-scale comparisons across species. 
ResMobiLys fills these gaps by enabling the recovery of plasmids from long-read assemblies and supporting the identification of diverse MGEs, while minimizing the reliance on reference genomes.The central component of the workflow are: (i) detection and annotation of AMR genes from long-read assembly,(ii) the identification of MGEs together with plasmid clustering and comparative analysis, and (iii) the integration of MGE–AMR associations, providing insights into both mobility and resistance dissemination. This versatility allows its application to a broad range of bacterial taxa, including isolates from unidentified lineages. By addressing these limitations, ResMobilYs provides a thorough analysis of MGE diversity, mobility, and resistance potential, enabling comprehensive profiling of the mobilome and resistome in diverse bacterial samples.
#Materials and methods

#Features
The pipeline is implemented with Nextflow (https://github.com/NjivaAndriamanga/resmobilys), which organizes and automates the different steps of the workflow, ensuring that analyses are carried out in a consistent and reproducible way. It can run on a wide range of computing systems, including high-performance computing (HPC) clusters, and processes samples in parallel to efficiently handle large datasets. To guarantee reliable and consistent execution, the pipeline uses containers to bundle all software dependencies. This approach makes ResMobilYs easy to deploy, portable across different environments, and well-suited for both small-scale and large-scale analyses.

#The workflow
**Input Data and Read Processing**
ResMobiLys accepts either long reads alone or a combination of long and short reads as input. Long-read sequencing enables more complete and less fragmented assemblies, significantly improving plasmid resolution.
**Reads quality control and genome assembly**
Genome assembly is conducted with Hybracter [@bouras_hybracter_2024], which provides highly accurate bacterial genome assemblies. Hybracter combines the strengths of long-read and hybrid assembly, making it one of the fastest and most reliable tools for recovering plasmids from whole-genome de novo assemblies. Taxonomic sequence classification is performed using Kraken2 [@wood_improved_2019]. Assembly quality is assessed using BUSCO [@manni_busco_2021], which provides metrics on completeness and accuracy.
**Antimicrobial Resistance genes, Virulence Factors, and heavy metal resistance**
ARGs are detected using RGI with CARD database [@alcock_card_2023], with CARD being the most comprehensive database and particularly well suited for environmental AMR surveillance [@papp_review_2022]. ABRicate (https://github.com/tseemann/abricate), with the MEGARes database [@bonin_megares_2023], is used to identify resistance determinants for metals and biocides. Virulence factors are identified against the VFDB database [@dong_expanded_2024].
**Mobile Genetic Elements annotation**
An estimated chromosome size can be provided or inferred automatically with Hybracter. Small circular contigs, below the estimated chromosome size, are inferred as a plasmid [@bouras_hybracter_2024]. For small linear contigs, plasmid identification is performed with PLASMe [@tang_plasme_2023] or PLATON [@schwengers_platon_2020], which use a hybrid approach combining sequence similarity (homology) with additional features or models to identify plasmid-derived sequences. This ensures precise separation of plasmidic and chromosomal contigs. This approach and these tools detect highly diverged plasmids without relying on traditional genetic markers, such as replicon or relaxase genes, often limited to well-studied bacteria like Enterobacteriaceae.
Plasmid typing and clustering were performed with MOB-suite [@robertson_mob-suite_2018], which includes MOB-typer and MOB-cluster modules. MOB-typer predicts replicon families, relaxase types, mate-pair formation types, and transferability, while MOB-cluster identifies similar plasmids across isolates based on complete sequence comparisons.
Integrons are identified using IntegronFinder [@neron_integronfinder_2022], which detects integrase using HMM profiles, predicts attC sites with sequence and structural motifs, and classifies elements based on their completeness and genomic context. ICEs are annotated using the CONJScan module of MacSyFinder [@cury_integrative_2017], and insertion sequences (ISs) are identified with ISEScan [@xie_isescan_2017]. Prophages are predicted with DBSCAN-SWA [@gan_dbscan-swa_2022], a rapid and accurate approach for identifying prophage regions. By prioritizing approaches that operate beyond direct pairwaise homology, the workflow enhances the detection of novel MGEs and  expand mobilome exploration to uncharacterized elements.
**Output**
ResMobiLys produces both isolate-specific and global results, summarizing antimicrobial resistance and mobility features across all analyzed genomes.
At the isolate level, the pipeline identifies antimicrobial resistance genes (ARGs), their genomic locations (chromosomal or plasmid), and associated mobile genetic elements (MGEs). The individual outputs from each tool integrated into the pipeline are retained and made available, enabling further downstream analyses.
At the global level, it generates:
- A presence–absence matrix of ARGs across isolates, including genomic context and plasmid cluster information.
- A list of ARG–MGE associations, highlighting genes located within mobile elements.
- A plasmid cluster summary, grouping plasmids from different isolates by sequence similarity and their associated ARGs.