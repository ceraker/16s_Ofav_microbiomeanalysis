# __*Orbicella faveolata* 16S microbiome analysis__
## __Comparing microbial and metabolomic data__
#### Author: Cassie Raker
#### Last updated: January 26, 2024

#### 16S sequencing performed by Genohub
#### Metabolite extraction performed at UT Austin ([protocol](https://github.com/ceraker/CatherineERaker_Notebook/blob/master/protocols/Raker_chemicalextraction.md))
#### Data uploaded and analyzed on Andromeda

#### Andromeda for URI

Information for creating an account on [Andromeda](https://web.uri.edu/hpc-research-computing/using-andromeda/).

###### Accessing Andromeda account
```
ssh -i ~/andromeda_keys ceraker@ssh3.hac.uri.edu
```

### Prepare a working directory with data
Navigate to data folder
```
cd /data/pradalab/ceraker/
```
Make a directory
```
mkdir omics
cd omics
```
Copy over necessary files from previous 16S analysis directory
```
cp /data/pradalab/ceraker/16s/og_reads/geno_meta_filt.tsv .
```

### Import files
First, load desired OTU tables into Andromeda.

Then convert OTU tables and metadata to .biom files
```
nano met_convert.sh
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
#SBATCH -D /data/pradalab/ceraker/omics
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load mmvec/1.0.5

# convert metabolite aboundance to .biom
biom convert -i metabolites.txt -o metabolites_hdf5.biom --table-type="OTU table" --to-hdf5

# incorporate sample metadata (same for both)
biom add-metadata -i metabolites_hdf5.biom -o metabolites_meta.biom -m meta_samp_met.txt

# incorporate observational metadata
biom add-metadata -i metabolites_meta.biom -o metabolites_done.biom -m met_obs_meta.txt
```
```
sbatch met_convert.sh
```
Initially had trouble adding the metadata: error messages stating that there was no header line in the file. Make sure that the first column is labeled the way BIOM expects (don't forget the #!):
- for sample metadata, #SampleID
- for observation metadata, #OTUID

Additionally, make sure that the first column has no repeating values.

Convert 16s data
```
nano mic_convert.sh
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
#SBATCH -D /data/pradalab/ceraker/omics
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load mmvec/1.0.5

# convert microbiome aboundance to .biom
biom convert -i metabolites.txt -o microbiome_hdf5.biom --table-type="OTU table" --to-hdf5

# incorporate sample metadata (same for both)
biom add-metadata -i microbiome_hdf5.biom -o microbiome_meta.biom --sample-metadata-fp meta_samp_met.txt

# incorporate observational metadata
biom add-metadata -i microbiome_meta.biom -o microbiome_done.biom --observation-metadata-fp mic_obs_meta.txt --sc-separated taxonomy
```
```
sbatch mic_convert.sh
```

Finally, import .biom formatted data to QIIME2
- Note: I realize it seems inefficient to convert the 16s data from QIIME data formats to .biom and back again, but since I did a fair bit of datacrunching in R post QIIME, this made the most sense to me.

```
nano import_biom.sh
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
#SBATCH -D /data/pradalab/ceraker/omics
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load QIIME2/2023.5

# Import microbiome data

qiime tools import \
        --input-path microbiome_done.biom \
        --output-path microbiome.qza \
        --type FeatureTable[Frequency]


# Import metabolomics data

qiime tools import \
        --input-path metabolites_done.biom \
        --output-path metabolites.qza \
        --type FeatureTable[Frequency]
```
```
sbatch import_biom.sh
```

### Run mmvec
##### Running mmvec on its own:

```
nano mmvec2.sh
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
#SBATCH -D /data/pradalab/ceraker/omics
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load mmvec/1.0.5

mmvec paired-omics \
        --microbe-file microbiome_done.biom \
        --metabolite-file metabolites_done.biom \
        --summary-dir summary
```
```
sbatch mmvec2.sh
```
This ran, but did produce some error messages. Output files were:
- latent_dim_3_input_prior_1.00_output_prior_1.00_beta1_0.90_beta2_0.95/
  - checkpoint                                       
  - model.ckpt-0.data-00000-of-00001  
  - model.ckpt-0.meta
  - events.out.tfevents.1702845687.n074.cluster.com  
  - model.ckpt-0.index
- latent_dim_3_input_prior_1.00_output_prior_1.00_beta1_0.90_beta2_0.95_embedding.txt
- latent_dim_3_input_prior_1.00_output_prior_1.00_beta1_0.90_beta2_0.95_ordination.txt
- latent_dim_3_input_prior_1.00_output_prior_1.00_beta1_0.90_beta2_0.95_ranks.txt

##### Running mmvec as QIIME2 plugin

Just loading the module in a bash script:
```
nano mmvec.sh
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
#SBATCH -D /data/pradalab/ceraker/omics
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load QIIME2/2023.5
module load mmvec/1.0.5
qiime dev refresh-cache  #need this command to enable mmvec plugin

qiime mmvec paired-omics \
        --i-microbes microbiome.qza \
        --i-metabolites metabolites.qza \
        --p-summary-interval 1 \
        --output-dir model_summary

qiime metadata tabulate \
        --m-input-file results/conditionals.qza \
        --o-visualization conditionals-viz.qzv
```
```
sbatch mmvec.sh
```

QIIME2 isn't recognizing the mmvec plugin.

_Troubleshooting_
- use an older QIIME2 version, 2019.7: failed
- conda install mmvec -c conda-forge: install process killed
- from URI HPC: It might be best to create your own conda environment with Miniconda3/22.11.1-1. Do this in an interactive session on a node (use the interactive command to get a shell on a compute node), since the login environment is too constrained to run conda.

```
interactive
curl -sL \
  "https://data.qiime2.org/distro/core/qiime2-2020.8-py36-osx-conda.yml" > \
  "qiime2.yml"
conda env create -n omics --file qiime2.yml

```
Install mmvec in the existing qiime2 environment:
```
conda install mmvec -c conda-forge
```
This timed out at the ```conda env create``` step, got alternate code from URI HPC:
```
interactive
conda create -n mmvec_env mamba python=3.7 -c conda-forge #worked
conda activate mmvec_env #worked (duh)
mamba install mmvec -c conda-forge
```
Got this error message, but mmvec seems to have installed?
```
Download error (28) Timeout was reached [https://conda.anaconda.org/conda-forge/noarch/keras-applications-1.0.8-py_1.tar.bz2]
SSL connection timeout
 ????U?4???Uyhd8ed1ab_0.tar.? extraction failed
error    libmamba Error when extracting package: [json.exception.type_error.316] invalid UTF-8 byte at index 2: 0xC2
terminate called after throwing an instance of 'std::bad_alloc'
  what():  std::bad_alloc
Aborted
```
Try to run mmvec: did not recognize command

##### Run mmvec in a container

Download `qiime_mmvec.def` file from URI HPC (thank you, Cecile!), and upload to Andromeda.

Here are the different steps to create the container:
1. Start an interactive session on andromeda: `srun --mem=50G -c 5 --time=02:00:00 --pty /bin/bash`
2. Build the container: `apptainer build qiime2_mmvec.sif qiime2_mmvec.def`
3. You can verify that qiime and mmvec are correctly installed by running: `apptainer exec qiime2_mmvec.sif qiime mmvec --help`
4. When running any subsequent commands, you must begin with `apptainer exec -B /data/pradalab/ceraker/omics/ qiime2_mmvec.sif` to operate in teh container and have access to files

THIS WORKED!!!

###### Re-run earlier file conversions within the container, so everything is using the same version of QIIME2 and mmvec

**Microbiome**

Convert microbiome aboundance to .biom
```
apptainer exec -B /data/pradalab/ceraker/omics qiime2_mmvec.sif qiime biom convert -i metabolites.txt -o microbiome_hdf5.biom --table-type="OTU table" --to-hdf5
```

Incorporate sample metadata (same for both)
```
apptainer exec -B /data/pradalab/ceraker/omics qiime2_mmvec.sif biom add-metadata -i microbiome_hdf5.biom -o microbiome_meta.biom --sample-metadata-fp meta_samp_met.txt
```

Incorporate observational metadata
```
apptainer exec -B /data/pradalab/ceraker/omics qiime2_mmvec.sif biom add-metadata -i microbiome_meta.biom -o microbiome_done.biom --observation-metadata-fp mic_obs_meta.txt --sc-separated taxonomy
```

**Metabolome**

Convert metabolite aboundance to .biom
```
apptainer exec -B /data/pradalab/ceraker/omics qiime2_mmvec.sif biom convert -i metabolites.txt -o metabolites_hdf5.biom --table-type="OTU table" --to-hdf5
```

Incorporate sample metadata (same for both)
```
apptainer exec -B /data/pradalab/ceraker/omics qiime2_mmvec.sif biom add-metadata -i metabolites_hdf5.biom -o metabolites_meta.biom -m meta_samp_met.txt
```

Incorporate observational metadata
```
apptainer exec -B /data/pradalab/ceraker/omics qiime2_mmvec.sif biom add-metadata -i metabolites_meta.biom -o metabolites_done.biom -m met_obs_meta.txt
```


Re-import data as qiime artifacts using QIIME2 2020.6, since this is the version loaded in the container.

```
apptainer exec -B /data/pradalab/ceraker/omics qiime2_mmvec.sif qiime tools import --input-path microbiome_done.biom --output-path microbiome2.qza --type FeatureTable[Frequency]

apptainer exec -B /data/pradalab/ceraker/omics qiime2_mmvec.sif qiime tools import --input-path metabolites_done.biom --output-path metabolites2.qza --type FeatureTable[Frequency]
```

###### Basic comparison
```
apptainer exec -B /data/pradalab/ceraker/omics/ qiime2_mmvec.sif qiime mmvec paired-omics \
          --i-microbes microbiome2.qza \
          --i-metabolites metabolites2.qza \
          --m-metadata-file meta_samp_met.txt \
          --p-summary-interval 60 \
          --output-dir model_summary2
```

```
apptainer exec -B /data/pradalab/ceraker/omics/ qiime2_mmvec.sif qiime mmvec paired-omics \
          --i-microbes microbiome2.qza \
          --i-metabolites metabolites2.qza \
          --m-metadata-file meta_samp_met.txt \
          --p-min-feature-count 10 \
          --p-num-testing-examples 10 \
          --p-learning-rate 1e-5 \
          --p-input-prior 1 \
          --p-output-prior 1 \
          --p-latent-dim 3 \
          --output-dir model_summary3 \
          --p-batch-size 1000 \
          --p-epochs 1000 \
          --p-summary-interval 10
```

Tabulate metadata
```
apptainer exec -B /data/pradalab/ceraker/omics/ qiime2_mmvec.sif qiime metadata tabulate \
        --m-input-file model_summary2/conditionals2.qza \
        --o-visualization model_summary2/conditionals-viz2.qzv
```


##### Further analysis
**QIIME2 emperor biplot**
```
apptainer exec -B /data/pradalab/ceraker/omics/ qiime2_mmvec.sif qiime emperor biplot \
        --i-biplot model_summary2/conditional_biplot.qza \
        --m-sample-metadata-file met_obs_meta.txt \
        --m-feature-metadata-file mic_tax_meta.txt \
        --p-ignore-missing-samples \
        --o-visualization emperor.qzv
```
*Plugin error from emperor:*

  *None of the feature identifiers match between the metadata and the coordinates. Verify that you are using metadata and coordinates corresponding to the same dataset.*

*Debug info has been saved to /tmp/qiime2-q2cli-err-4gb_c408.log*

**Heatmap**
```
apptainer exec -B /data/pradalab/ceraker/omics/ qiime2_mmvec.sif qiime mmvec heatmap \
        --i-ranks model_summary2/conditionals.qza \
        --m-microbe-metadata-file mic_tax_meta.txt \
        --m-microbe-metadata-column family \
        --m-metabolite-metadata-file met_obs_meta.txt \
        --m-metabolite-metadata-column subclass \
        --p-level 2 \
        --o-visualization mmvec-heatmap.qzv
```
*Plugin error from mmvec:*

  *index 0 is out of bounds for axis 0 with size 0*

*Debug info has been saved to /tmp/qiime2-q2cli-err-t9gt1pfs.log*


Export conditionals
```
apptainer exec -B /data/pradalab/ceraker/omics qiime2_mmvec.sif qiime tools export --input-path model_summary2/conditionals.qza --output-path model_summary2
```

Okay the conditionals table doesn't look correct: it looks like it's only looking at the metabolome, not the microbiome.

Try importing files from solo mmvec run (MUST run this in container: `FeatureData[Conditional]` does not exist in QIIME2/2023.5):

```
apptainer exec -B /data/pradalab/ceraker/omics/ qiime2_mmvec.sif qiime tools import --input-path summary/latent_dim_3_input_prior_1.00_output_prior_1.00_beta1_0.90_beta2_0.95_ranks.txt --output-path conditionals2.qza --type FeatureData[Conditional]
```
```
apptainer exec -B /data/pradalab/ceraker/omics/ qiime2_mmvec.sif qiime tools import --input-path summary/latent_dim_3_input_prior_1.00_output_prior_1.00_beta1_0.90_beta2_0.95_ordination.txt --output-path ordination2.qza --type 'PCoAResults % Properties("biplot")'
```

```
apptainer exec -B /data/pradalab/ceraker/omics/ qiime2_mmvec.sif qiime mmvec heatmap \
        --i-ranks conditionals2.qza \
        --m-microbe-metadata-file mic_tax_meta.txt \
        --m-microbe-metadata-column family \
        --m-metabolite-metadata-file met_obs_meta.txt \
        --m-metabolite-metadata-column subclass \
        --p-level 5 \
        --o-visualization mmvec-heatmap.qzv \
        --verbose
```

**Model training troubleshooting**

Submitted a queryto the [QIIME2 forum](https://forum.qiime2.org/t/no-model-stats-qza-file/29044) about the lack of `model_stats` file in model output. Developer suggested it was because the model did not train for long enough.
```
apptainer exec -B /data/pradalab/ceraker/omics/ qiime2_mmvec.sif qiime mmvec paired-omics \
          --i-microbes microbiome2.qza \
          --i-metabolites metabolites2.qza \
          --m-metadata-file meta_samp_met.txt \
          --p-learning-rate 1e-1 \
          --p-latent-dim 1 \
          --p-epochs 100 \
          --p-summary-interval 1 \
          --output-dir model_summary4
```
This ran longer, but still did not generad a `model_stats` file.

Up `--p-epochs` to 1000:
```
apptainer exec -B /data/pradalab/ceraker/omics/ qiime2_mmvec.sif qiime mmvec paired-omics \
          --i-microbes microbiome2.qza \
          --i-metabolites metabolites2.qza \
          --m-metadata-file meta_samp_met.txt \
          --p-learning-rate 1e-1 \
          --p-latent-dim 1 \
          --p-epochs 1000 \
          --p-summary-interval 1 \
          --output-dir model_summary5
```
```
apptainer exec -B /data/pradalab/ceraker/omics/ qiime2_mmvec.sif qiime emperor biplot \
        --i-biplot model_summary5/conditional_biplot5.qza \
        --m-sample-metadata-file met_obs_meta.txt \
        --m-feature-metadata-file mic_tax_meta.txt \
        --o-visualization emperor5.qzv
```































.
