
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
  - Identification of antimicrobial resistance genes (ARGs) with RGI/CARD, biocide and heavy metal resistance genes, and virulence factors
  - Identification of MGEs (plasmids, transposable elements, integrons, prophages)
## 🧬 Output

- **Table associating species with MGEs and resistance genes**  
  Presence–absence matrix showing antimicrobial resistance genes (ARGs) across samples.  
  Each gene is marked according to its genomic location (chromosome or plasmid) and includes plasmid cluster information when available.

- **Table listing resistance genes and their MGE associations**  
  Detailed list of all identified ARGs with genomic positions and associations with mobile genetic elements (MGEs) such as integrons, prophages, and ICEs.

- **Detailed list of plasmids and clusters**
  Summary table of plasmids with their cluster IDs and associated ARGs, providing an overview of plasmid diversity and possible roles in resistance gene dissemination.

---

### 🧾 Example Entries

#### Presence–Absence Table

| Sample  | soxR | FosA2 | TEM-1 |
|----------|------|-------|-------|
| sample1  | 0 | C | P1(AA019) |
| sample2  | C | 0 | P2(AA019) |

#### Plasmid–ARG Summary

| Sample  | Plasmid      | Cluster | ARGs         |
|----------|--------------|---------|--------------|
| sample1  | plasmid00001 | AA015   | soxR, FosA2  |
| sample1  | plasmid00002 | AA019   | TEM-1        |
| sample2  | plasmid00001 | AA003   |              |
| sample2  | plasmid00002 | AA019   | TEM-1        |


## 🚀 Installation

Ensure the following dependencies are installed before running **ResMobiLys**:

### ✅ Nextflow

Install Nextflow by following the instructions [here](https://www.nextflow.io/docs/latest/install.html)

### ✅ Apptainer/Singularity or Docker

- To install **Apptainer/Singularity**, follow the guide [here](https://apptainer.org/docs/admin/main/installation.html)
- To install **Docker**, follow the guide [here](https://www.docker.com/get-started/)

### ✅ Clone the ResMobiLys Repository

```bash
# Clone the project repository
git clone https://github.com/NjivaAndriamanga/resmobilys.git
cd resmobilys
git submodule update --init --recursive
```

## 📦 Database and test data

Most databases are already provided, but the **PLASMe** database must be downloaded manually.

1. Download the database from Zenodo:  
   👉 [Download PLASMe database](https://zenodo.org/record/8046934/files/DB.zip?download=1)
```bash
#Download and unzip plasme database
wget https://zenodo.org/record/8046934/files/DB.zip
xport UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE
unzip DB.zip
```

2. Move and unzip the database inside the `ResMobiLys` directory

3. Test data
Test datasets can be downloaded from [here]

## 🚀 **Usage**

Prepare an `index_file.csv` containing metadata for your samples.  
For your **first run**, it is recommended to use the provided **test dataset**.  
During this initial execution, all required databases and environments (tools and dependencies) will be downloaded automatically.  
This step may take several minutes depending on your internet connection.

Once the setup is complete, you can run the pipeline on your own dataset using the `-resume` option to avoid re-downloading components.

```bash
# First run with test dataset
nextflow run waterisk -profile test,singularity/apptainer -resume

# Run on your dataset
nextflow run waterisk -profile slurm/local,singularity/apptainer -resume -c waterisk/personal.config -resume

### **Configuration Options**
Modify `personal.config` to customize execution parameters, including computing resources and software parameters.

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
 ResMobilYs integrates the output of the following tools :
- Hybracter (Bouras et al. 2024)
- Kraken2 (Wood, Lu, et Langmead 2019)
- Busco (Manni et al. 2021)
- Abricate (https://github.com/tseemann/abricate)
- MEGARes (Bonin et al. 2023)
- VFDB (Dong et al. 2024)
- PLASMe (Tang et al. 2023)
- IntegronFinder (Néron et al. 2022)
- DBSCAN-SWA (Gan et al. 2022)
- TnFinder (Ross et al. 2021)
- Mob-suite (Robertson et Nash 2018)
## **Contact**
If you have any questions or support, please open an issue.
