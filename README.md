# **ResMobiLys: A Nextflow Pipeline for Resistome and Mobilome Analysis of Field and Clinical Eubacteria Isolates**

## **Introduction**
Mobile Genetic Elements (MGEs) are DNA segments capable of moving within or between genomes, facilitating the transfer of genetic material across microbial populations. MGEs play a key role in microbial ecology by spreading genes linked to antimicrobial resistance (AMR), which poses significant public health challenges. Advances in high-throughput sequencing technologies allow large-scale studies of diverse bacterial genomes from various origins, generating vast data volumes that require efficient and integrated analysis.

## **Purpose**
ResMobiLys is a Nextflow pipeline designed for comprehensive mobilome and resistome analysis. It provides an end-to-end solution, from de novo assembly to precise identification of MGEs such as plasmids, integrons, prophages, and transposable elements. It also detects Antibiotic Resistance Genes (ARGs) and virulence factors, whether associated with MGEs or not, enabling streamlined and detailed genomic feature analysis.

## **Features**
- **Input:**
  - Mandatory: Long-read FASTQ files
  - Optional: Short-read FASTQ files
- **Workflow Overview:**
  - De novo assembly with Hybacter
  - Taxonomy assignment with Kraken2
  - Assembly assessment with busco
  - Identification of antimicrobial resistance genes (ARGs), biocide and heavy metal resistance genes, and virulence factors
  - Identification of MGEs (plasmids, transposable elements, integrons, prophages)
- **Output:**
  - Table associating species with MGEs and resistance genes
  - Table listing resistance genes and their MGE associations
  - Detailed list of plasmids, integrons, and transposable elements with notes on their potential interactions

## **Installation**
Ensure the following dependencies are installed before running ResMobiLys:

To install nextflow: 
Via conda: here
Manually: here

To install singularity:
Follow the instructions at the singularity website

```bash

# Clone the project repository
git clone https://github.com/yourusername/ResMobiLys.git
cd ResMobiLys
```

## **Usage**
Prepare an `index_file.csv` containing metadata for your samples. Then, run the pipeline as follows:

```bash
nextflow run main.nf --input index_file.csv -profile singularity
```

### **Configuration Options**
Modify `nextflow.config` to customize execution parameters, including computing resources and software dependencies.

## **Outputs**
Upon successful execution, ResMobiLys generates:
- Annotated genomic assemblies
- Summary tables linking species, MGEs, and ARGs
- Reports on plasmids, integrons, and transposable elements
- Taxonomic classification of input sequences

## **Contributing**
We welcome contributions! Feel free to open issues or submit pull requests.

## **License**
This project is licensed under the MIT License. See the `LICENSE` file for details.

## **Citation**
If you use ResMobiLys in your research, please cite it as follows:
```
Author(s), ResMobiLys: A Nextflow Pipeline for Resistome and Mobilome Analysis, Year.
```

## **Contact**
For questions or support, please open an issue or contact the developers at [your.email@example.com](mailto:your.email@example.com).
