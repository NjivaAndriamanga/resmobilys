
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
  - Detection of similar plasmids
    
## 🧬 Output

- **Table associating species with MGEs and resistance genes**  
  Presence–absence matrix showing antimicrobial resistance genes (ARGs) across samples.  
  Each gene is marked according to its genomic location (chromosome or plasmid) and includes plasmid cluster information when available.

| Sample  | soxR | FosA2 | TEM-1 |
|----------|------|-------|-------|
| sample1  | 0 | C | P1(AA019) |
| sample2  | C | 0 | P2(AA019) |


- **Table listing resistance genes and their MGE associations**  
  Detailed list of all identified ARGs with genomic positions and associations with mobile genetic elements (MGEs) such as integrons, prophages, and ICEs.

| Sample   | Localisation | Other informations | ARG | MGEs |
|----------|--------------|---------|--------------|--------|
| sample1 | Chromosome    | ....    | FosA2 | Integron1 |
| sample1 | Plasmid1      | ....    | TEM-1 | Integron2 |
| sample2 | Chromosome    | ....    | soxR | prophage |
| sample1 | Plasmid2      | ....    | TEM-1 | Integron1 |


- **Detailed list of plasmids and clusters**
Summary table of plasmids with their cluster IDs and associated ARGs, providing an overview of plasmid diversity and possible roles in resistance gene dissemination.

| Sample  | Plasmid      | Cluster | ARGs         |
|----------|--------------|---------|--------------|
| sample1  | plasmid00001 | AA015   |              |
| sample1  | plasmid00002 | AA019   | TEM-1        |
| sample2  | plasmid00001 | AA003   |              |
| sample2  | plasmid00002 | AA019   | TEM-1        |


## 🚀 Installation

To run ResMobiLys, ensure the following are available on your system:
- Nextflow
- Apptainer/Singularity or Docker
- Git
- PLASMe database

### ✅ Nextflow

Install Nextflow by following the instructions [here](https://www.nextflow.io/docs/latest/install.html)

### ✅ Apptainer/Singularity or Docker

- To install **Apptainer/Singularity**, follow the guide [here](https://apptainer.org/docs/admin/main/installation.html)
- To install **Docker**, follow the guide [here](https://www.docker.com/get-started/)

### ✅ Clone the ResMobiLys Repository and submodule

To install **git** : [here](https://git-scm.com/install)

```bash
# Clone the project repository
git clone https://github.com/NjivaAndriamanga/resmobilys.git
cd resmobilys
git submodule update --init --recursive
```

## 📦 Database

Most databases are already provided, but the **PLASMe** database must be downloaded manually.

1. Download the database from Zenodo:  
   👉 [Download PLASMe database](https://zenodo.org/record/8046934/files/DB.zip?download=1)
```bash
#Download and unzip plasme database
wget https://zenodo.org/record/8046934/files/DB.zip
export UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE
unzip DB.zip
```

2. Move and unzip the database inside the `ResMobiLys` directory


## 🚀 **Usage**

Prepare an `index_file.csv` containing metadata for your samples.  
For your **first run**, it is recommended to use the provided **test dataset**.  
During this initial execution, all required databases and environments (tools and dependencies) will be downloaded automatically.  
This step may take several minutes depending on your internet connection.

Once the setup is complete, you can run the pipeline on your own dataset using the `-resume` option to avoid re-downloading components.

```bash
# First run with test dataset
nextflow run resmobilys -profile test,singularity/apptainer -resume

# Run on your dataset
nextflow run resmobilys -profile slurm/local,singularity/apptainer -resume -c waterisk/personal.config

'''
```
### **Configuration Options**
Modify `personal.config` to customize execution parameters, including computing resources and software parameters.

## Outputs

Upon successful execution, **ResMobiLys** generates the following in the `resmobylis_output` directory:

### 1. Final output files
- `args_mges.tsv`  
- `merged_plasmid_table.tsv`  
- `presence_absence_with_clusters.tsv`  

### 2. Intermediate/tool-specific outputs
These files can be used for further analyses, allowing exploration of results from each individual tool.

## **Citation**
If you use ResMobiLys in your research, please cite it as follows:

Publication in process...

 ResMobilYs integrates the following tools and databases :
- Hybracter (Bouras et al. 2024)
- Kraken2 (Wood, Lu, et Langmead 2019)
- Busco (Manni et al. 2021)
- Abricate (https://github.com/tseemann/abricate)
- MEGARes (Bonin et al. 2023)
- VFDB (Dong et al. 2024)
- PLASMe (Tang et al. 2023)
- IntegronFinder (Néron et al. 2022)
- DBSCAN-SWA (Gan et al. 2022)
- Mob-suite (Robertson et Nash 2018)

## 🤝 Contributing and Support

We welcome community involvement to improve and extend the **ResMobiLys pipeline**.  
Please follow the guidelines below if you wish to contribute, report issues, or request support.

---

### 🧩 Contributing

Contributions are encouraged and greatly appreciated!  
To contribute:

1. **Fork** the repository on GitHub.  
2. **Create a new branch** for your feature or bugfix:  
   ```bash
   git checkout -b feature/new-feature-name
   ```
3. Commit your changes with clear messages
4. Push to your fork and submit a Pull request describing your contribution

### 🐛 Reporting issues
If you encounter a bug or unexpected behavior:
1. Check the Issues section to see if it has already been reported.
2. If not, open a **new issue** and include:
   - A clear description of the problem
   - The command or configuration you used
   - Relevant log or error messages (if available)
This helps us reproduce and fix the issue efficiently.

### 💬 Seeking support
For questions about usage, workflow customization, or interpretation of results:
  - Open an issue
  - Or contact the maintainers directly
We will do our best to provide assistance and improve the documentation based on user feedback
  
