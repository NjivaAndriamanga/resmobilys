
# **ResMobiLys pipeline is still in development**


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

## 🚀 Installation

Ensure the following dependencies are installed before running **ResMobiLys**:

### ✅ Nextflow

Install Nextflow by following the instructions [here](https://www.nextflow.io/docs/latest/install.html).

### ✅ Apptainer/Singularity or Docker

- To install **Apptainer/Singularity**, follow the guide [here](https://apptainer.org/docs/admin/main/installation.html).
- To install **Docker**, follow the guide [here](https://www.docker.com/get-started/).

### ✅ Clone the ResMobiLys Repository

```bash
# Clone the project repository
git clone https://github.com/yourusername/ResMobiLys.git
cd ResMobiLys
git submodule update --init --recursive
```

## 📦 Database

Most databases are already provided, but the **PLASMe** database must be downloaded manually.

1. Download the database from Zenodo:  
   👉 [Download PLASMe database](https://zenodo.org/record/8046934/files/DB.zip?download=1)

2. Move and unzip the database inside the `ResMobiLys` directory

## **Usage**
Prepare an `index_file.csv` containing metadata for your samples. Then, run the pipeline as follows:

```bash
nextflow run waterisk -profile slurm/local,singularity/apptainer -resume -c waterisk/personal.config
```

### **Configuration Options**
Modify `personal.config` to customize execution parameters, including computing resources and software dependencies.

## **Outputs**
Upon successful execution, ResMobiLys generates:
- Annotated genomic assemblies
- Summary tables linking species, MGEs, and ARGs
- Reports on plasmids, integrons, and transposable elements
- Taxonomic classification of input sequences

## **Citation**
If you use ResMobiLys in your research, please cite it as follows:
```
Publication in process...
```

## **Contact**
If you have any questions or support, please open an issue.
