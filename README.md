
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
Then, clone ResMobiLys repository:
```bash
# Clone the project repository
git clone https://github.com/NjivaAndriamanga/resmobilys.git
cd resmobilys
git submodule update --init --recursive
```

## 📦 Database

Most databases are downloaded automatically; however, the **PLASMe** database must be downloaded manually.

1. Download the database from Zenodo:  
   👉 [Download PLASMe database](https://zenodo.org/record/8046934/files/DB.zip?download=1)
```bash
#Download and unzip plasme database
wget https://zenodo.org/record/8046934/files/DB.zip
```

2. Move and unzip the database inside the `ResMobiLys` directory
```bash
export UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE
unzip DB.zip
```


## 🚀 **Usage**

### 1️⃣ Run ResMobilYs on the test dataset  

For your **first execution**, we strongly recommend running ResMobilYs on the provided **test dataset**.  
This allows you to verify that the pipeline is correctly installed and ensures that all required databases, tools, and environments are downloaded automatically.  
  
⚠️ The test run takes approximately 1 hour on a machine with 20 CPUs and 32 GB of RAM, with the most demanding processes using up to 8 CPUs and 16 GB of memory.  

```bash
# First run with test dataset
nextflow run resmobilys -profile test,singularity/apptainer/docker -resume
```

### 2️⃣ Run ResMobilYs on your own dataset  

**2.1 Prepare the index_file.csv**  
To analyze your own data, you must prepare an index_file.csv containing the input data and metadata for each sample.  
Each row corresponds to one sample.  

An example file is available in the test/ directory.  

The file must include the following columns:  
- **LR_fastq**: Path to the long-read FASTQ file
- **genome_size**: Estimated chromosome size in base pairs. If unknown, use 0.
- **SR1** (optional): Path to the first short-read pair (R1). Leave empty if no short reads are available.
- **SR2** (optional): Path to the second short-read pair (R2). Leave empty if no short reads are available.

Example:  
```
LR_fastq,chrm_size,SR1_fastq,SR2_fastq    
resmobilys/test/17_01_bar09.fastq.gz,0,,  
resmobilys/test/31_03_bar52.fastq.gz,0,,  
resmobilys/test/NB10_LR.fastq.gz,1000000,resmobilys/test/NB10_1.fastq.gz,resmobilys/test/NB10_2.fastq.gz
```
**2.2 Configure and run the pipeline**  
Once the index_file.csv is ready, **edit the personal.config file** by adding the path to your index file in the **index_file** field.  
To run ResMobilYs, only one parameter is mandatory: the path to the index_file.csv, provided either via the --index_file option or defined directly in the personal.config file. For other parameters, a default value will be assigned.  
All other parameters—such as computing resources, input/output paths, and software options—can be customized in personal.config (see Configuration Options for details).  

#### Execution profiles  
You must select **one execution environment** and **one container engine**:
- **Execution environment (choose one):**
    . **local**: run the pipeline on a local workstation or server
    . **slurm/PBS**: run the pipeline on an HPC cluster managed by SLURM/PBS
- **Container engine**
    . **docker/singularity/apptainer**: singularity and apptainer are recommended for HPC environments.

```bash

# Run on your dataset
nextflow run resmobilys -profile slurm/local,singularity/apptainer -resume -c personal.config

```

### **Configuration Options**

Decritption of all parameters: 
```
--help                           [boolean, string] Show the help message for all top level parameters. When a parameter is given to `--help`, the full help message of that parameter will be printed. 
--helpFull                       [boolean]         Show the help message for all non-hidden parameters. 
--showHidden                     [boolean]         Show all hidden parameters in the help message. This needs to be used in combination with `--help` or `--helpFull`. 

##Required inputs
  --output_dir                   [string] Output directory for the pipeline [default: resmobilys_output] 
  --index_file                   [string] The absolute path of the index file that contains estimated chomosome length for each isolate. If set to null, it will be estimated automatically and long read assembly will be performed. 

##Quality control assessment
  --trim_end_size                [integer] Number of bases to remove at each end of the read [default: 0] 
  --quality_trim                 [integer] Parameter can be used to trim low-quality ends from reads. [default: 20] 
  --lineage_db                   [string]  Lineage database for busco. See https://busco.ezlab.org/list_of_lineages.html for a list of available lineage databases. [default: bacteria_odb10] 

##Hybracter parameters
  --read_min_length              [integer] Minimum read length [default: 1000] 
  --medaka                       [boolean] Run medaka [default: true] 
  --flyeModel                    [string]  Flye assembly model  (accepted: --nano-hq, --nano-corr, --nano-raw, --pacbio-raw, --pacbio-corr, --pacbio-hifi) [default: --nano-hq] 

##Tools parameters
  --plasmid                      [string]  Plasmid prediction tool  (accepted: plasme, platon) [default: plasme] 
  --plasme_db                    [string]  Plasme database directory [default: ${projectDir}/plasme_db] 
  --platon_db                    [string]  platon database directory [default: ${projectDir}/platon_db] 
  --kraken_index                 [string]  kraken index directory, by default it will download the smaller index capped at 8 gb [default: https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz] 
  --kraken_db                    [string]  kraken database directory 
  --kraken_taxonomy              [boolean] Enable kraken taxonomy [default: true] 
  --vf_db                        [string]  VF database directory [default: ${projectDir}/vf_db] 
  --rgi_include_nudge            [boolean] include hits nudged from loose to strict hits in the rgi output [default: true] 
  --ice_avg_size                 [integer] ICE average size [default: 0] 
  --compositeIS_size             [integer] CompositeIS average size [default: 0] 

##blastn parameters
  --evalue_vf                    [number]  E-value cutoff for VF detection [default: 1E-10] 
  --pident_vf                    [integer] Percent identity cutoff for VF detection [default: 80] 

```

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
- CONJScan (Cury et al. 2017)

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
  
