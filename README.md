# Pipeline_Gasperini_2019
IGVF Gasperin processing using IGVF single cell like pipeline

## Downloading and unpacking a small fraction of the Gasperini 2019 pilot dataset   
This download and unpack can take a few hours and consume some resourcers.
This some script suggestion of how to download and unpack the dataset.  
change the var __dir_creation__ with the full path of your server/local enviroment.  

```shell
#!/bin/bash
#SBATCH -n 20                              # 10 core
#SBATCH -t 0-07:00                         # in D-HH:MM format
#SBATCH -p short                           # Run in short partition
#SBATCH --mem=80G
#SBATCH -o download.out
#SBATCH -e download.err
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=youremail@gmail.com   # Email to which notifications will be sent

NTHREADS=20; 
#change the dir creation directory  or remove from slurm format to adapt to your server
dir_creation='gasperini_2019_pipeline_crispr'
mkdir $dir_creation
cd $dir_creation
wget https://sra-pub-src-1.s3.amazonaws.com/SRR7967482/pilot_highmoi_screen.1_SI_GA_G1.bam.1;mv pilot_highmoi_screen.1_SI_GA_G1.bam.1 pilot_highmoi_screen.1_SI_GA_G1.bam
./bamtofastq_linux --nthreads="$NTHREADS" pilot_highmoi_screen.1_SI_GA_G1.bam bam_pilot_scrna_1
wget https://github.com/10XGenomics/bamtofastq/releases/download/v1.4.1/bamtofastq_linux; chmod +x bamtofastq_linux
wget https://sra-pub-src-1.s3.amazonaws.com/SRR7967488/pilot_highmoi_screen.1_CGTTACCG.grna.bam.1;mv pilot_highmoi_screen.1_CGTTACCG.grna.bam.1 pilot_highmoi_screen.1_CGTTACCG.grna.bam
./bamtofastq_linux --nthreads="$NTHREADS" pilot_highmoi_screen.1_CGTTACCG.grna.bam bam_pilot_guide_1
```
---
## Installing the PERTURB-SEQ single-cell single cell Environment
```shell
conda create -n perturbseq_like_pipeline_2023 python=3.8  
conda activate perturbseq_like_pipeline_2023  
conda install -c conda-forge mamba -y   
pip install nextflow  
pip uninstall GTFProcessing -y  
pip install git+https://github.com/LucasSilvaFerreira/GTFProcessing.git  
pip install gtfparse==1.3.0  
pip install git+https://github.com/LucasSilvaFerreira/Perturb_Loader.git  
pip install --quiet kb-python  
conda install -c conda-forge mamba -y  
mamba install -c bioconda nextflow -y  
mamba install -c bioconda kallisto -y  
mamba install -c anaconda openpyxl -y  
mamba install -c conda-forge r-base -y  
mamba install -c conda-forge r-gert -y  
mamba install -c conda-forge r-ragg -y  
mamba install -c conda-forge r-ggplot2 -Y  
mamba install -c conda-forge r-biocmanager -Y   
Rscript -e 'BiocManager::install("Rhdf5lib")'  
Rscript -e 'BiocManager::install("rhdf5")'  
Rscript -e 'BiocManager::install("ShortRead")'  
Rscript -e 'install.packages("doParallel", repos = "http://cran.us.r-project.org")'  
mamba conda install -c conda-forge r-devtools -y 
Rscript -e 'devtools::install_github("katsevich-lab/sceptre")'  
Rscript -e 'devtools::install_github("chris-mcginnis-ucsf/MULTI-seq")'  
pip install muon  
mamba install -c bioconda scrublet 
pip install pybiomart  
mamba create -n pygenomictracks pygenometracks==3.8 -c bioconda -y  
mamba install -c conda-forge matplotlib-inline -y
mamba install -c anaconda ipython -y  
```


## Downloading additional required data and config the run
```shell
#Download white list
#Provide the path where the whitelist was saved to the pipeline config file
wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt
```

## Capturing kallisto path bin
```shell
#Add the result of this comand  to the kallisto path in the pipeline config
source activate perturbseq_like_pipeline_2023
which kallisto

```

## Cloning the Pipeline Repository

```shell

git clone https://github.com/LucasSilvaFerreira/pipeline_perturbseq_like.git

```



## Creating the config file

Save this config as gasperini_sample.config  

__OBS:__ This file is assuming our your data is inside the directory __/n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/__  
(modify it to your desired path). 

-  df_from_gasperini_tss.xlsx is provided in this repository. Add this to your path. This file only contains a few guides targeting promoters for different genes.


```shell
params.GTF_GZ_LINK = 'http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz'
params.TRANSCRIPTOME_REFERENCE = "human"
params.KALLISTO_BIN = '/home/lf114/miniconda3/envs/perturbseq_pipeline/envs/perturbseq_like_pipeline_2023/bin/kallisto' //PATHWAY KALLISTO INSTALL check your install path (previous command)
params.GENOME =       'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'

params.GUIDE_FEATURES = '/n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/pipeline_perturbseq_like/df_from_gasperini_tss.xlsx'

params.CHEMISTRY = "10XV2"
params.THREADS = 15
params.DISTANCE_NEIGHBORS = 1000000
params.IN_TRANS = "FALSE"
params.FASTQ_FILES_TRANSCRIPTS = ['/n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L001_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L001_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L001_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L001_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L001_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L001_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HWVT7BGX3/bamtofastq_S1_L001_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HWVT7BGX3/bamtofastq_S1_L001_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L002_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L002_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L002_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L002_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L002_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L002_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HWVT7BGX3/bamtofastq_S1_L002_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HWVT7BGX3/bamtofastq_S1_L002_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L003_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L003_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L003_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L003_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L003_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L003_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HWVT7BGX3/bamtofastq_S1_L003_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HWVT7BGX3/bamtofastq_S1_L003_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L004_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L004_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L004_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L004_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L004_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L004_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HWVT7BGX3/bamtofastq_S1_L004_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_scrna_1/K1000_CRISPRi_cells_r1_SI-GA-G1_MissingLibrary_1_HWVT7BGX3/bamtofastq_S1_L004_R2_001.fastq.gz'] 
params.FASTQ_NAMES_TRANSCRIPTS = ['S1_L1'] 

params.CUSTOM_REFERENCE = false
params.CUSTOM_REFERENCE_IDX = ''
params.CUSTOM_REFERENCE_T2T = ''
params.CUSTOM_GTF_PATH = ''

params.FASTQ_FILES_GUIDES = ['/n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L001_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L001_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L001_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L001_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L001_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L001_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L002_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L002_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L002_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L002_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L002_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L002_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L003_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L003_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L003_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L003_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L003_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L003_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L004_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMKLBGX3/bamtofastq_S1_L004_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L004_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HMMMLBGX3/bamtofastq_S1_L004_R2_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L004_R1_001.fastq.gz  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/bam_pilot_guide_1/K1000_CRISPRi_gRNA_r1_CGTTACCG_MissingLibrary_1_HW5MKBGX3/bamtofastq_S1_L004_R2_001.fastq.gz'] 
params.FASTQ_NAMES_GUIDES = ['S1_L1'] 

params.CREATE_REF = false
params.ADDGENENAMES = ''
params.DIRECTION = 'both'

params.WHITELIST= '/n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/737K-august-2016.txt'

params.EXPECTED_CELL_NUMBER = 10000
params.MITO_SPECIE = 'hsapiens'
params.MITO_EXPECTED_PERCENTAGE = 0.2
params.PERCENTAGE_OF_CELLS_INCLUDING_TRANSCRIPTS=0.01
params.TRANSCRIPTS_UMI_TRHESHOLD = 200
params.GUIDE_UMI_LIMIT = 3

params.MERGE = false


//MULTISEQ

params.RUN_MULTISEQ = false
params.R1_MULTI = 'NOT APPLICABLE'
params.R2_MULTI = 'NOT APPLICABLE'
params.BARCODES_MULTIBAR_LIST_MULTI = 'NOT APPLICABLE'
params.BAR_MULTI= [1,16]
params.UMI_MULTI= [17,28]
params.R2_MULTI_TAG = [1,8]
```


---

## Running the Piepline
-Change the paths to match the paths in  your local machine/server

- main.nf: This will be inside your PERTURBSEQ git cloned repository 
- TOWER_ACCESS_TOKEN: Get your own token in the nextower website (This will keeping updating the pipeline status  and other information over the internet), You can see different runs directly on your webrowser
- These are recommended resourcers but it should run fast with a half of this (10 threads/50G memory)

```shell

#!/bin/bash
#SBATCH -n 20                              # 10 core
#SBATCH -t 0-07:00                         # in D-HH:MM format
#SBATCH -p short                           # Run in short partition
#SBATCH --mem=100G
#SBATCH -o g_debug_transcriptome_01_both.out
#SBATCH -e g_debug_transcriptome_01_both.err
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=your@email.com   # Email to which notifications will be sent
#activating the enviroment
source activate perturbseq_like_pipeline_2023
export TOWER_ACCESS_TOKEN=eyJ0aWQiOiA3MzEyfS42OGZhZWNiZDVlMTE3ODhhMDAxMTgwMGRjOGE5MTZkZDQzNTU3OTU5
export NXF_VER=22.10.6
cd /scratch3/users/l/lf114/gasperini_2019_pipeline_crispr
nextflow run main.nf -c /scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/gasperini_sample.config   -with-tower -w  /n/scratch3/users/l/lf114/gasperini_2019_pipeline_crispr/desired_output_dir  
```