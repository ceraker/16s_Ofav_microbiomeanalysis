# __*Orbicella faveolata* 16S microbiome analysis__
#### Author: Cassie Raker
#### Last updated: March 14, 2024

#### Sequencing performed by Genohub
#### Data uploaded and analyzed on Andromeda
Analysis pipeline adapted from [Dr. Emma Strand](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2021-06-21-16s-Analysis-Pipeline.md)

### Andromeda for URI

Information for creating an account on [Andromeda](https://web.uri.edu/hpc-research-computing/using-andromeda/).

###### Accessing Andromeda account
```
ssh -i ~/andromeda_keys ceraker@ssh3.hac.uri.edu
```

## General Workflow:  
1. [Prepare directory](#Directory)         
5. [QIIME2 Import](#Import)  
3. [Denoising](#Denoise)  
6. [Clustering](#Cluster)  
9. [Align and Classify](#Classify)       
12. [Diversity Analysis](#Diversity)   


### <a name="Directory"></a> **1. Prepare Directory**
###### Prepare a working directory with data
navigate to data folder
```
cd /data/pradalab/ceraker/
```
Make a directory
```
mkdir 16s
cd 16s
```

### <a name="Import"></a> **2. QIIME2 Import**
*Information, documentation, and tutorials for the QIIME2 analysis pipeline can be found [here](https://docs.qiime2.org/2023.5/)*

Directions for formatting and examples of manifest and metadata files compatible with qiime2 can be found [here](https://docs.qiime2.org/2023.5/tutorials/importing/)

###### Create manifest
First moved raw reads into their own subdirectory, then used this directory to create a txt file with all the file names
```
cd og_reads
ls > rawfiles.txt
```
Use this txt file to create `genohub_manifest_raw.txt`

###### Created the metadata in RStudio and Excel
(this is very messy and inefficient but it worked for me)

```
#IN RSTUDIO
#read back in with correct column names
samplekey <- read.csv("Raker_2019_PR_SAMPLEKEY.csv", check.names = FALSE)
View(samplekey)

#merge metadata with lab ID info
meta_lab <- left_join(samplekey, meta2, by = join_by("ATTRIBUTE_sampleID"=="ATTRIBUTE_sampleID"))
View(meta_lab)

#export meta_lab to fill in blanks
write.csv(meta_lab, "meta_lab.csv", row.names = TRUE)

#aaaaaand read it back in
meta_lab <- read.csv("meta_lab.csv", check.names = FALSE)
View(meta_lab)

#read in list of IDs from genohub
labID_only <- read.csv("labID_only.csv", check.names = FALSE)
View(labID_only)

#merge meta_lab with labID_only
geno_meta <- left_join(labID_only, meta_lab, by = join_by(labID==labID))
View(geno_meta)

#get rid of duplicate rows using dplyr
geno_meta_v2 <- geno_meta %>% distinct()

#export geno_meta to use in metadata file for QIIME
write.csv(geno_meta_v2, "geno_meta_v2.csv", row.names = FALSE)
```
Then opened file in Excel to format the way QIIME2 wants it (headers, etc)

- Included a **#q2:types** row to ensure that acclimation is treated as a categorical variable

- Validate that the metadata file follows QIIME2 formatting rules using the **Keemei** add-on in **Google Sheets**

- Saved directly from Google Sheets as **tsv**: _**geno_meta_v2.tsv**_

###### Import raw sequences to QIIME2
```
nano qiime_import_raw.sh
```
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ceraker@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/ceraker/16s/og_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load QIIME2/2023.5

# Metadata path
METADATA="/data/pradalab/ceraker/16s/og_reads/geno_meta_v2.txt"

# Sample manifest path
MANIFEST="/data/pradalab/ceraker/16s/og_reads/genohub_manifest_raw.txt"

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $MANIFEST \
  --output-path 16s_genohub_raw.qza \
  --input-format PairedEndFastqManifestPhred33
```
```
sbatch qiime_import_raw.sh
```

Create a summary of sequences:
```
nano seq_sum.sh
```
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ceraker@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/ceraker/16s/og_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh # this is used so that them function can be found in zsh and bash

module load QIIME2/2023.5

# Metadata path
METADATA="/data/pradalab/ceraker/16s/og_reads/geno_meta_v2.txt"

# Sample manifest path
MANIFEST="/data/pradalab/ceraker/16s/og_reads/genohub_manifest_raw.txt"

qiime demux summarize --i-data /data/pradalab/ceraker/16s/og_reads/16s_genohub_raw.qza --o-visualization /data/pradalab/ceraker/16s/og_reads/16s_genohub_raw.qzv
```
```
sbatch seq_sum.sh
```

### <a name="Denoise"></a> **3. Denoising**
This is a filtering step, so it's important to compare different parameters. The basic script was as follows:
```
nano denoise_BLUE_long.sh
```
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ceraker@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/ceraker/16s/og_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh # this is used so that them function can be found in zsh and bash

module load QIIME2/2023.5

# Metadata path
METADATA="/data/pradalab/ceraker/16s/og_reads/genohub_meta_v2.tsv"

# Sample manifest path
MANIFEST="/data/pradalab/ceraker/16s/og_reads/genohub_manifest_raw.txt"

qiime dada2 denoise-paired \
--i-demultiplexed-seqs /data/pradalab/ceraker/16s/og_reads/16s_genohub_raw.qza \
--p-trim-left-f 54 \
--p-trim-left-r 52 \
--p-trunc-len-f 250 \
--p-trunc-len-r 250 \
--p-n-threads 36 \
--o-table table_blue_l.qza \
--o-representative-sequences /data/pradalab/ceraker/16s/og_reads/rep-seqs_blue_l.qza \
--o-denoising-stats /data/pradalab/ceraker/16s/og_reads/denoising-stats_blue_l.qza

qiime metadata tabulate \
  --m-input-file denoising-stats_blue_l.qza \
  --o-visualization denoising-stats_blue_l.qzv
```
```
sbatch denoise_BLUE_long.sh
```

DADA2 scripts were all variations of the above script, with modifications as follows:

Script Name: **denoise_RED_long.sh**:  
- `--p-trunc-len` choice: 270 reverse and 270 forward. This was based on the quality scores of the reads.  
- `--p-trim-left` choice: 52 reverse and 54 forward. This was based on the primer lengths: 515F = 52 bp long; 806RB = 54 bp long. This include adapter overhang.

Script Name: **denoise_RED_short.sh**:  
- `--p-trunc-len` choice: 270 reverse and 270 forward.  
- `--p-trim-left` choice: 16 reverse and 16 forward

Script Name: **denoise_YELLOW_long.sh**:  
- `--p-trunc-len` choice: 260 for both forward and reverse.  
- `--p-trim-left` choice: 52 reverse and 54 forward.

Script Name: **denoise_YELLOW_short.sh**:  
- `--p-trunc-len` choice: 260 for both forward and reverse.  
- `--p-trim-left` choice: 16 reverse and 16 forward.


#### Summary of the number of samples that passed the denoising stage.
In general, the different denoising parameters ended with the same amount of samples passed and failed. The "yellow" and "long" parameters resulted in a higher percentage saved, a higher percentage merged, and a higher percentage non-chimeric.

**denoise_BLUE_long.sh** was chosen as the best denoising option.

- **Total treatment groups:** 47
- **Treatment groups PASSED:** 42 have three or more replicates
- **Treatment groups FAILED:** 5 have only two replicates
  - Group One DW: three months in cage
  - Group Four BY: two weeks in cage
  - Group Four IR: two weeks in cage
  - Group Five DW: never in cage
  - Group Five GT: never in cage

### <a name="Cluster"></a> **4. Clustering**
##### Clustering with chosen denoising parameters
Cluster raw files that got denoised, using BLUE_long output files
```
nano cluster_blue.sh
```
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ceraker@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/ceraker/16s/og_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load QIIME2/2023.5

# Metadata path
METADATA="/data/pradalab/ceraker/16s/og_reads/genohub_meta_v2.tsv"

# Sample manifest path
MANIFEST="/data/pradalab/ceraker/16s/og_reads/genohub_manifest_raw.txt"

#### CLUSTERING

# Summarize feature table and sequences
qiime metadata tabulate \
  --m-input-file denoising-stats_blue_l.qza \
  --o-visualization denoising-stats_blue_l.qzv
qiime feature-table summarize \
  --i-table table_blue_l.qza \
  --o-visualization table_blue_l.qzv \
  --m-sample-metadata-file $METADATA
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_blue_l.qza \
  --o-visualization rep-seqs_blue_l.qzv
```
```
sbatch cluster_blue.sh
```

### <a name="Align"></a> **5. Align and Classify**

###### Download classifier from QIIME2
```
wget https://data.qiime2.org/2023.5/common/silva-138-99-515-806-nb-classifier.qza
```
###### Classifying samples denoised with blue_long.

```
nano taxonomy_blue.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ceraker@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/ceraker/16s/og_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load QIIME2/2023.5

# Metadata path
METADATA="/data/pradalab/ceraker/16s/og_reads/genohub_meta_v2.tsv"

# Sample manifest path
MANIFEST="/data/pradalab/ceraker/16s/og_reads/genohub_manifest_raw.txt"

#### TAXONOMY CLASSIFICATION

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs_blue_l.qza \
  --o-classification taxonomy_blue.qza

qiime taxa filter-table \
     --i-table table_blue_l.qza \
     --i-taxonomy taxonomy_blue.qza \
     --p-mode contains \
     --p-exclude "Unassigned","Chloroplast","Eukaryota" \
     --o-filtered-table table-filtered_blue.qza

qiime metadata tabulate \
    --m-input-file taxonomy_blue.qza \
    --o-visualization taxonomy_blue.qzv
qiime taxa barplot \
    --i-table table-filtered_blue.qza \
    --i-taxonomy taxonomy_blue.qza \
    --m-metadata-file $METADATA \
    --o-visualization taxa-bar-plots-filtered_blue.qzv
qiime metadata tabulate \
    --m-input-file rep-seqs_blue_l.qza \
    --m-input-file taxonomy_blue.qza \
    --o-visualization tabulated-feature-metadata_blue.qzv

qiime feature-table summarize \
    --i-table table-filtered_blue.qza \
    --o-visualization table-filtered_blue.qzv \
    --m-sample-metadata-file $METADATA

#### CREATES PHYLOGENETIC TREES

# align and mask sequences
qiime alignment mafft \
  --i-sequences rep-seqs_blue_l.qza \
  --o-alignment aligned-rep-seqs_blue_l.qza
qiime alignment mask \
  --i-alignment aligned-rep-seqs_blue_l.qza \
  --o-masked-alignment masked-aligned-rep-seqs_blue_l.qza

# calculate tree
qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs_blue_l.qza \
  --o-tree unrooted-tree_blue.qza
qiime phylogeny midpoint-root \
  --i-tree unrooted-tree_blue.qza \
  --o-rooted-tree rooted-tree_blue.qza
```
```
sbatch taxonomy_blue.sh
```

Looking at the feature table summaries, some samples that passed denoising still have very low reads, and may need to be dropped later. It will depend on if they seem like outliers, and how many replicates we have for each treatment group.

When all samples with <10000 reads are filtered out, the totals are:

- **Total treatment groups:** 47
- **Treatment groups PASSED:** 39 have three or more replicates
- **Treatment groups FAILED:**
  - 6 have only two replicates
    - Group Two FU: one month in cage
    - Group Four BY: two weeks in cage
    - Group Four IR: two weeks in cage
    - Group Five DW: never in cage
    - Group Five FU: never in cage
    - Group Five IR: never in cage
  - 2 have only one replicate
    - Group One DW: three months in cage  
    - Group Five GT: never in cage

When all samples with <1000 reads are filtered out, the totals are:

- **Total treatment groups:** 47
- **Treatment groups PASSED:** 42 have three or more replicates
- **Treatment groups FAILED:**
  - 5 have only two replicates
    - Group One DW: three months in cage
    - Group Four BY: two weeks in cage
    - Group Four IR: two weeks in cage
    - Group Five DW: never in cage
    - Group Five GT: never in cage

### <a name="Diversity"></a> **6. Diversity Analysis**
###### Perform initial diversity analyses in QIIME2

The diversity script creates the output directory core-metrics-results, do not create it ahead of time.
```
nano diversity.sh
```
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ceraker@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/ceraker/16s/og_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load QIIME2/2023.5

# Metadata path
METADATA="/data/pradalab/ceraker/16s/og_reads/genohub_meta_v2.tsv"

# Sample manifest path
MANIFEST="/data/pradalab/ceraker/16s/og_reads/genohub_manifest_raw.txt"

#########################

#### CALCULATES OVERALL DIVERSITY
## change sub sampling depth values

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree_blue.qza \
  --i-table table-filtered_blue.qza \
  --p-sampling-depth 1000 \
  --m-metadata-file $METADATA \
  --output-dir core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $METADATA \
  --m-metadata-column Timepoint \
  --o-visualization core-metrics-results/unweighted-unifrac-station-significance.qzv \
  --p-pairwise
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $METADATA  \
  --m-metadata-column acclimation \
  --o-visualization core-metrics-results/unweighted-unifrac-group-significance.qzv \
  --p-pairwise

# This script calculates the rarefaction curve for the data
  qiime diversity alpha-rarefaction \
    --i-table table-filtered_blue.qza \
    --i-phylogeny rooted-tree_blue.qza \
    --p-max-depth 40000 \
    --m-metadata-file $METADATA \
    --o-visualization alpha-rarefaction.qzv
```
```
sbatch diversity.sh
```

Calculate the beta diversity by genotype
```
nano diversity_gen.sh
```
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ceraker@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/ceraker/16s/og_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load QIIME2/2023.5

# Metadata path
METADATA="/data/pradalab/ceraker/16s/og_reads/genohub_meta_v2.tsv"

# Sample manifest path400
MANIFEST="/data/pradalab/ceraker/16s/og_reads/genohub_manifest_raw.txt"

#########################

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $METADATA  \
  --m-metadata-column genotype \
  --o-visualization core-metrics-results/unweighted-unifrac-group-significance_gen.qzv \
  --p-pairwise
```
```
sbatch diversity_gen.sh
```

Rarefaction curves look to plateau between 20000 and 25000 reads, with some samples cutting off before that (see earlier note).
- Except for the few samples that have fewer than 1000 reads, samples with lower reads look to plateau before 10000 reads and I feel reasonably confident that we are capturing all the diversity.

I removed the five samples with fewer than 1000 reads: 2GT40, 2FU45, 8CX23, 2GT9, 2HS6. An additional outlier, 8GT22, was identified when performing PCAs (see script pca.rmd), and was also removed. Removing these outliers did not result in any additional treatment groups having fewer than three replicates.

I created a new metadata file without these samples, and used that to filter the feature table.

```
nano filter.sh
```
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ceraker@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/ceraker/16s/og_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load QIIME2/2023.5

# Metadata path
METADATA="/data/pradalab/ceraker/16s/og_reads/geno_meta_filt.tsv"

# Sample manifest path
MANIFEST="/data/pradalab/ceraker/16s/og_reads/genohub_manifest_raw.txt"

#########################

qiime feature-table filter-samples \
    --i-table table-filtered_blue.qza \
    --o-filtered-table filtered-table.qza \
    --m-metadata-file $METADATA

qiime feature-table summarize \
    --i-table filtered-table.qza \
    --o-visualization filtered-table.qzv \
    --m-sample-metadata-file $METADATA
```
```
sbatch filter.sh
```

Rerun the diversity analysis with outlier samples removed:
```
nano diversity_filt.sh
```
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ceraker@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/ceraker/16s/og_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load QIIME2/2023.5

# Metadata path
METADATA="/data/pradalab/ceraker/16s/og_reads/genohub_meta_v2.tsv"

# Sample manifest path
MANIFEST="/data/pradalab/ceraker/16s/og_reads/genohub_manifest_raw.txt"

#########################

#### CALCULATES OVERALL DIVERSITY
## change sub sampling depth values

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree_blue.qza \
  --i-table filtered-table.qza \
  --p-sampling-depth 1000 \
  --m-metadata-file $METADATA \
  --output-dir core-metrics-filt

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-filt/faith_pd_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization core-metrics-filt/faith-pd-group-significance.qzv
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-gilt/evenness_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization core-metrics-filt/evenness-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-filt/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $METADATA \
  --m-metadata-column Timepoint \
  --o-visualization core-metrics-filt/unweighted-unifrac-station-significance.qzv \
  --p-pairwise
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-filt/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $METADATA  \
  --m-metadata-column acclimation \
  --o-visualization core-metrics-filt/unweighted-unifrac-group-significance.qzv \
  --p-pairwise

# This script calculates the rarefaction curve for the data
  qiime diversity alpha-rarefaction \
    --i-table filtered-table.qza \
    --i-phylogeny rooted-tree_blue.qza \
    --p-max-depth 40000 \
    --m-metadata-file $METADATA \
    --o-visualization alpha-rarefaction.qzv
```
```
sbatch diversity_filt.sh
```

Re-calculate the beta diversity by genotype

```
nano diversity_gen_filt.sh
```
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ceraker@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/ceraker/16s/og_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load QIIME2/2023.5

# Metadata path
METADATA="/data/pradalab/ceraker/16s/og_reads/genohub_meta_v2.tsv"

# Sample manifest path
MANIFEST="/data/pradalab/ceraker/16s/og_reads/genohub_manifest_raw.txt"

#########################

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-filt/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $METADATA  \
  --m-metadata-column genotype \
  --o-visualization core-metrics-filt/unweighted-unifrac-group-significance_gen.qzv \
  --p-pairwise
```
```
sbatch diversity_gen_filt.sh
```
*Copy all .qzv files to location outside Andromeda to visualize and export stats*
