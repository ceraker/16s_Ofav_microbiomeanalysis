# __*Orbicella faveolata* Omic Analysis__
## __Comparing microbial and metabolomic data__
#### Author: Cassie Raker
#### Last updated: September 2, 2024

#### 16S sequencing performed by Genohub
#### Metabolite extraction performed at UT Austin ([protocol](https://github.com/ceraker/CatherineERaker_Notebook/blob/master/protocols/Raker_chemicalextraction.md))
#### Data uploaded and analyzed on the Unity cluster

#### Unity for URI

Information for creating an account on [Unity](https://its.uri.edu/research-computing/using-unity/).

###### Accessing Unity account
```
ssh -l ceraker_uri_edu unity.uri.edu
```

### Prepare a working directory with data
Navigate to Prada PI group
```
cd /work/pi_prada_uri_edu/
```
Make a directory
```
mkdir ceraker
cd ceraker
```
Copy over necessary microbe, metabolite, and metadata text files using Cyberduck
```
metabolites2.txt
metabolites.txt
met_obs_meta.txt
meta_samp_met.txt
mic_obs_meta.txt
microbes_b.txt #use this one
microbes.txt
mic_tax_meta.txt
taxonomy.tsv
```

### Installing mmvec alone

##### Create environment with necessary modules
First, request an interactive node
```
salloc -c 1 -n 1 -p uri-cpu
```
Load miniconda and create a conda environment
```
module load miniconda/22.11.1-1 #have to do this every time
conda create -n omics
conda activate omics
```

All mmvec install information can be found on the [Github](https://github.com/biocore/mmvec/blob/master/README.md)
```
conda create -n mmvec_env mamba python=3.7 -c conda-forge
conda activate mmvec_env
mamba install mmvec -c conda-forge
```

Now that programs are installed and environment is created, just need to perform these commands when starting a new session:
```
salloc -c 1 -n 1 -p uri-cpu
module load miniconda/22.11.1-1
conda activate #whichever other environment I am using
```

### Installing mmvec as a QIIME2 plugin
*Note: mmvec is only supported by QIIME2 2020.6 or earlier*

Install [QIIME2 2020.6](https://docs.qiime2.org/2020.6/install/native/):
```
salloc --mem=30G --partition=cpu-preempt
conda install wget
wget https://data.qiime2.org/distro/core/qiime2-2020.6-py36-linux-conda.yml
```

Modify .yml file per instructions from HPC (i.e. use nano and go in and change the code manually). Diff below:
```
- scipy=1.4.1  #don't change: for orientation
-  - seaborn=0.10.1  #remove this
-  - seaborn-base=0.10.1  #remove this
+  - seaborn  #add this
- send2trash=1.5.0  #don't change: for orientation
```
Resume building environment
```
conda env create -n qiime2-2020.6 --file qiime2-2020.6-py36-linux-conda.yml

#optional cleanup
rm qiime2-2020.6-py36-linux-conda.yml  #I didn't do this because I wanted to save the file in case I needed to modify it again

#check install
conda activate qiime2-2020.6
qiime --help
```

Pip install mmvec:
```
pip install mmvec
qiime dev refresh-cache

#check install
mmvec --help
```

### Convert files to mmvec compatible formats
Convert metabolite data to .biom
```
# convert metabolite aboundance to .biom
biom convert -i metabolites.txt -o metabolites_hdf5.biom --table-type="OTU table" --to-hdf5

# incorporate sample metadata (same for both)
biom add-metadata -i metabolites_hdf5.biom -o metabolites_meta.biom -m meta_samp_met.txt

# incorporate observational metadata
biom add-metadata -i metabolites_meta.biom -o metabolites_done.biom -m met_obs_meta.txt
```

Manipulated microbiome data in R so that column headers are coral IDs. May still be an issue with metabolite data having more samples than microbe data, but we'll see.

Upload new microbiome data 'mic_fin.txt' using Cyberduck

Create new .biom file for microbes
```
# convert microbiome aboundance to .biom
biom convert -i mic_fin.txt -o microbiome_hdf5.biom --table-type="OTU table" --to-hdf5

# incorporate sample metadata (same for both)
biom add-metadata -i microbiome_hdf5.biom -o microbiome_meta.biom --sample-metadata-fp meta_samp_met.txt

# incorporate observational metadata
biom add-metadata -i microbiome_meta.biom -o microbiome_done.biom --observation-metadata-fp mic_obs_meta.txt --sc-separated taxonomy
```

### Run mmvec solo
```
mmvec paired-omics \
        --microbe-file microbiome_done.biom \
        --metabolite-file metabolites_done.biom \
        --summary-dir summary
```

Output:
- latent_dim_3_input_prior_1.00_output_prior_1.00_beta1_0.90_beta2_0.95
  - checkpoint
  - events.out.tfevents.1723752746.uri-cpu014
  - model.ckpt-0.data-00000-of-00001  model.ckpt-0.meta
  - model.ckpt-0.index
- latent_dim_3_input_prior_1.00_output_prior_1.00_beta1_0.90_beta2_0.95_embedding.txt
- latent_dim_3_input_prior_1.00_output_prior_1.00_beta1_0.90_beta2_0.95_ordination.txt
- latent_dim_3_input_prior_1.00_output_prior_1.00_beta1_0.90_beta2_0.95_ranks.txt

### Run mmvec as a QIIME2 plugin
##### Import data
Import metabolite data into QIIME2:
```
qiime tools import \
        --input-path metabolites_done.biom \
        --output-path metabolites.qza \
        --type FeatureTable[Frequency]
```

Import microbiome data into QIIME2:
```
qiime tools import \
        --input-path microbiome_done.biom \
        --output-path microbiome_done.qza \
        --type FeatureTable[Frequency]
```
*Note: I found it worth it to export the microbiome data from qiime, convert it to .biom format and then import it back into qiime. This guaranteed that all the metadata was incorporated correctly, and in a format consistent with use in mmvec.*

##### Run analysis in qiime
```
qiime mmvec paired-omics \
        --i-microbes microbiome_done.qza \
        --i-metabolites metabolites.qza \
        --p-summary-interval 1 \
        --output-dir model_summary1
```

Output in model_summary1:
- conditional_biplot.qza
- conditionals.qza
- model_stats.qza

*I first ran the model without any additional parameters, just to ensure everything was formatted correctly and mmvec was generating all the proper output files.*

##### Train model
Examine model stats file:
```
qiime mmvec summarize-single \
        --i-model-stats model_stats.qza \
        --o-visualization model-summary.qzv
```

Generate null model with only biases:
```
qiime mmvec paired-omics \
        --i-microbes microbiome_done.qza \
        --i-metabolites metabolites.qza \
        --p-latent-dim 0 \
        --p-summary-interval 1 \
        --output-dir null_summary
```

Compare model output to null and generate pseudo-Q2 value:
```
qiime mmvec summarize-paired \
        --i-model-stats model_summary1/model_stats.qza \
        --i-baseline-stats null_summary/model_stats.qza \
        --o-visualization paired-summary.qzv
```
Pseudo Q-squared: -0.359562

This initial model was overfit, to I experimented with different parameters, as chronicles on this [spreadsheet](https://github.com/ceraker/16s_Ofav_microbiomeanalysis/blob/main/omics/mmvec_stats.xlsx).

I ended up running mmvec on individual experimental groups (acclimation time or original colony), as this yielded much more powerful models.


### Treatment specific models
##### Acclimation Groups
Subset microbiome data by acclimation time. Should only need to do this with one omics dataset: mmvec will automatically discard any samples that aren't present in both.

Created new OTU tables in R, then loaded into Unity using Cyberduck

Create new microbe files for mmvec (repeat for each acclimation time)
```
# convert microbiome aboundance to .biom
biom convert -i mic_one.txt -o microbiome_hdf5_1.biom --table-type="OTU table" --to-hdf5

# incorporate sample metadata (same for both)
biom add-metadata -i microbiome_hdf5_1.biom -o microbiome_meta_1.biom --sample-metadata-fp meta_samp_met.txt

# incorporate observational metadata
biom add-metadata -i microbiome_meta_1.biom -o microbiome_done_acc1.biom --observation-metadata-fp mic_obs_meta.txt --sc-separated taxonomy
```

Import into qiime
```
qiime tools import \
        --input-path microbiome_done_acc1.biom \
        --output-path microbiome_done_acc1.qza \
        --type FeatureTable[Frequency]
```

Run null model:
```
qiime mmvec paired-omics \
        --i-microbes microbiome_done_acc1.qza \
        --i-metabolites metabolites.qza \
        --p-latent-dim 0 \
        --p-summary-interval 1 \
        --output-dir null_summary_acc1
```

Run mmvec with best parameters (see spreadsheet). These parameters resulted in the model needing a much longer amount of time, so I submitted these runs as bash scripts:
```
nano mmvec_acc1_2.sh
```
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -c 4
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ceraker@uri.edu
#SBATCH --account=pi_prada_uri_edu
#SBATCH -D /work/pi_prada_uri_edu/ceraker
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load miniconda/22.11.1-1

conda activate qiime2-2020.6

qiime mmvec paired-omics \
        --i-microbes acclimation/microbiome_done_acc1.qza \
        --i-metabolites metabolites.qza \
        --p-training-column train \
        --p-learning-rate 1e-3 \
        --p-summary-interval 1 \
        --p-input-prior 1 \
        --p-output-prior 1 \
        --p-epochs 200 \
        --p-batch-size 5 \
        --p-latent-dim 2 \
        --p-min-feature-count 5 \
        --output-dir model_summary_acc1_2
```
```
sbatch mmvec_acc1_2.sh
```

Compare to null:
```
qiime mmvec summarize-paired \
        --i-model-stats model_summary_acc1_2/model_stats.qza \
        --i-baseline-stats null_summary_acc1/model_stats.qza \
        --o-visualization paired-summary_acc1_2.qzv
```

Pseudo Q-squared values:
- Acc 1: 0.412999
- Acc 2: 0.635704
- Acc 3: -0.492223
- Acc 4: 0.399474
- Acc 5: 0.800511

##### Original Colony ("Genotype") Groups
Subset microbiome data by original colony. Should only need to do this with one omics dataset: mmvec will automatically discard any samples that aren't present in both.

Created new OTU tables in R, then loaded into Unity using Cyberduck

Create new microbe files for mmvec (repeat for each original colony)
```
# convert microbiome aboundance to .biom
biom convert -i mic_BY.txt -o microbiome_hdf5_BY.biom --table-type="OTU table" --to-hdf5

# incorporate sample metadata (same for both)
biom add-metadata -i microbiome_hdf5_BY.biom -o microbiome_meta_BY.biom --sample-metadata-fp meta_samp_met.txt

# incorporate observational metadata
biom add-metadata -i microbiome_meta_BY.biom -o microbiome_done_BY.biom --observation-metadata-fp mic_obs_meta.txt --sc-separated taxonomy
```

Import into qiime
```
qiime tools import \
        --input-path microbiome_done_BY.biom \
        --output-path microbiome_done_BY.qza \
        --type FeatureTable[Frequency]
```

Run null model:
```
qiime mmvec paired-omics \
        --i-microbes microbiome_done_BY.qza \
        --i-metabolites metabolites.qza \
        --p-latent-dim 0 \
        --p-summary-interval 1 \
        --output-dir null_summary_BY
```

*Note: there were too few samples from colony AZ to run mmvec on that group alone*

Ran these groups in interactive mode since they didn't need to run as long:
```
qiime mmvec paired-omics \
        --i-microbes microbiome_done_BY.qza \
        --i-metabolites metabolites.qza \
        --p-training-column train \
        --p-learning-rate 1e-3 \
        --p-summary-interval 1 \
        --p-input-prior 0.1 \
        --p-output-prior 0.1 \
        --p-epochs 200 \
        --p-batch-size 5 \
        --p-latent-dim 2 \
        --p-min-feature-count 5 \
        --output-dir model_summary_BY_3
```
Compare to null:
```
qiime mmvec summarize-paired \
        --i-model-stats model_summary_BY_3/model_stats.qza \
        --i-baseline-stats null_summary_BY/model_stats.qza \
        --o-visualization paired-summary_BY_3.qzv
```

Pseudo-Q2:
- BY: 0.713771
- CX: 0.750141
- DW: 0.760828
- EV: 0.706305
- FU: 0.589150
- GT: 0.716759
- HS: 0.756788
- IR: 0.813571
- JQ: 0.514692


### Data visualization
*Note: before this stage, I reorganized files into different subdirectories, so some filepaths have changed.*
Tabulate to view log conditional probabilities:
```
qiime metadata tabulate \
        --m-input-file model_summary_acc5_2/conditionals.qza \
        --o-visualization model_summary_acc5_2/conditionals-viz.qzv
```

Basic heatmap, make one for each acclimation group:
```
qiime mmvec heatmap \
        --i-ranks model_summary_acc1_2/conditionals.qza \
        --m-microbe-metadata-file mic_obs_meta.txt \
        --m-microbe-metadata-column taxonomy \
        --m-metabolite-metadata-file met_obs_meta.txt \
        --m-metabolite-metadata-column subclass \
        --p-level 5 \
        --o-visualization model_summary_acc1_2/ranks-heatmap_acc1.qzv
```

Paired heatmaps:
```
qiime mmvec paired-heatmap \
  --i-ranks model_summary_acc1_2/conditionals.qza \
  --i-microbes-table ../microbiome_done.qza \
  --i-metabolites-table ../metabolites.qza \
  --m-microbe-metadata-file mic_obs_meta.txt \
  --m-microbe-metadata-column taxonomy \
  --p-top-k-microbes 5 \
  --p-normalize rel_row \
  --p-top-k-metabolites 100 \
  --p-level 5 \
  --o-visualization model_summary_acc1_2/paired-heatmap_acc1_top3.qzv
```
- Acc 1: top 5
- Acc 2: top 5
- Acc 3: top 5
- Acc 4: top 3 (bug for top 5)
- Acc 5: top 3 (bug for top 5)

Can view these heatmaps in the [QIIME2 viewer](https://view.qiime2.org/), and export them as pdfs as well as export the raw data.

***Further analysis and data visualization in R***
