# __*Orbicella faveolata* 16S microbiome analysis__
#### Author: Cassie Raker
#### Last updated: December 11, 2023

#### Sequencing performed by Genohub
#### Data uploaded and analyzed on Andromeda

### Andromeda for URI

Information for creating an account on [Andromeda](https://web.uri.edu/hpc-research-computing/using-andromeda/).

###### Accessing Andromeda account
```
ssh -i ~/andromeda_keys ceraker@ssh3.hac.uri.edu
```

## General Workflow:  
1. [Prepare directory](#Directory)      
2. [QC sequences](#QC)    
5. [QIIME2 Import](#Import)    
6. [Aligning](#Align)    
7. [Preclustering](#Precluster)    
8. [Identify chimeras](#Chimera)  
9. [Classify sequences](#Classify)     
10. [Cluster OTU's](#Cluster)    
11. [Subsampling](#Subsample)  
12. [Calculate ecological statistics](#Statistics)  
13. [Output data for R analysis](#Output)   


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
Load modules
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ceraker@uri.edu
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load Miniconda3/4.9.2
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2
module load QIIME2/2023.5
```
Other packages that gave me errors:
```
module load fastp/0.23.2-GCC-11.2.0
module load cutadapt/2.10-GCCcore-9.3.0-Python-3.8.2
module load Mothur/1.46.1-foss-2020b
```
Create a conda environment
```
conda create -n 16s
conda activate 16s
```
Transfer fastq files from home computer
```
scp -r -i /Users/cassieraker/Documents/_URI/_Research/Data/2019_PR_ch3_microbiome/Raker_genohub_fall2023/*fastq.gz ceraker@ssh3.hac.uri.edu:/data/pradalab/ceraker/16s
```
##### Quality Analysis using FastQC
FastQC is used to identify any problems with reads, such as low quality or contamination from TruSeq adapters.
```
mkdir fastqc_results
cd fastqc_results
```
Create fastqc script
```
nano fastqc.sh
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
#SBATCH -D /data/pradalab/ceraker/16s                   
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in /data/pradalab/ceraker/16s/*fastq.gz
do
fastqc $file --outdir /data/pradalab/ceraker/16s/fastqc_results         
done

multiqc --interactive fastqc_results  

mv multiqc_report.html 16S_genohub_qc_multiqc_report_CR.html
```
```
sbatch fastqc.sh
```
Move MultiQC results to home computer
```
scp -i ceraker@ssh3.hac.uri.edu:/data/pradalab/ceraker/16s/16S_genohub_qc_multiqc_report_CR.html /Users/cassieraker/Documents/_URI/_Research/Data/2019_PR_ch3_microbiome
```
cannot get it to copy: no such file or directory

Move MultiQC results to github repository
```
cp /data/pradalab/ceraker/16s/16S_genohub_qc_multiqc_report_CR.html /data/pradalab/ceraker/repos/16s_Ofav_microbiomeanalysis/genohub_test
```
##### Quality control using fastp
fastp is used to trim low quality reads, remaining adapter sequences, etc.
```
nano fastp.sh
```
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ceraker@uri.edu #your email to send notifications
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/ceraker/16s/scripts                   
#SBATCH --error="script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="output_script" #once your job is completed, any final job report comments will be put in this file

source /usr/share/Modules/init/sh # this is used so that them function can be found in zsh and bash

module load fastp/0.23.2-GCC-11.2.0

#use fastp to clean files
for filename in /data/pradalab/ceraker/16s/*.fastq.gz
  do
  fastp -i /data/pradalab/ceraker/16s/*.fastq.gz -o ${filename}.clean.fastq.gz \
  --adapter_sequence= TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG --adapter_sequence_r2= GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC
  done
```
run script
```
sbatch fastp.sh
```
Willow's script:
```
for r1_file in ./*R1*.fastq.gz
do r2_file="${r1_file/_R1_/_R2_}"
fastp -i $r1_file -I $r2_file -o ${r1_file/.fastq.gz/}.fq.gz -O ${r2_file/.fastq.gz/}.fq.gz --detect_adapter_for_pe --$
done
```
FastQC on clean files
```
mkdir fastqc_results_clean
```
```
nano fastqc_clean.sh
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
#SBATCH -D /data/pradalab/ceraker/16s                   
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in /data/pradalab/ceraker/16s/*fastq.gz.clean
do
fastqc $file --outdir /data/pradalab/ceraker/16s/fastqc_results_clean        
done

multiqc --interactive fastqc_results_clean

mv multiqc_report_1.html 16S_genohub_qc_multiqc_report_clean_CR.html
```
Move MultiQC results to github repository
```
cp /data/pradalab/ceraker/16s/16S_genohub_qc_multiqc_report_clean_CR.html /data/pradalab/ceraker/repos/16s_Ofav_microbiomeanalysis/genohub_test
```

There was no checksum file, so I created one with clean reads using md5sum
```
# create file
md5sum *.fastq.gz.clean > 16s_genohub_checksum.md5
```
check that individual files went into new md5 file
```
md5sum -c 16s_genohub_checksum.md5
```

### QUIIME2
#### Create metadata
```
mkdir metadata
```
Move clean files into their own directory
```
mkdir fastq_clean
mv *clean.fastq.gz ./clean_reads
cd clean_reads
```
Create a text file with just the file names
```
ls > ../cleanfiles.txt
```
Copy this txt file to wherever allows you to download it to your own computer and edit it into the proper manifest format

```
scp -i ceraker@ssh3.hac.uri.edu:/data/pradalab/ceraker/16s/cleanfiles.txt /Users/cassieraker/Documents/_URI/_Research/Data/2019_PR_ch3_microbiome
```
import cleaned sequences into QIIME2 as QUIIME artifact

samples are Casava 1.8 paired-end demultiplexed fastq format

make sure script is in same location as directory with files
```
nano qiime_import.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/clean_reads                  
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load QIIME2/2023.5

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path ​/data/pradalab/ceraker/16s/clean_reads \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path 16s_genohub.qza
```

```
sbatch qiime_import.sh
```
FAILED: exit code 1
input path does not exist

Attempt with only a subset of ten files:
```
106_S216_R1_001.fastq.gz.clean.fastq.gz
137_S65_R1_001.fastq.gz.clean.fastq.gz
1_S173_R1_001.fastq.gz.clean.fastq.gz
222_S203_R1_001.fastq.gz.clean.fastq.gz
31_S127_R1_001.fastq.gz.clean.fastq.gz
40_S13_R1_001.fastq.gz.clean.fastq.gz
54_S190_R1_001.fastq.gz.clean.fastq.gz
63_S98_R1_001.fastq.gz.clean.fastq.gz
81_S162_R1_001.fastq.gz.clean.fastq.gz
94_S100_R1_001.fastq.gz.clean.fastq.gz
106_S216_R2_001.fastq.gz.clean.fastq.gz
137_S65_R2_001.fastq.gz.clean.fastq.gz
1_S173_R2_001.fastq.gz.clean.fastq.gz
222_S203_R2_001.fastq.gz.clean.fastq.gz
31_S127_R2_001.fastq.gz.clean.fastq.gz
40_S13_R2_001.fastq.gz.clean.fastq.gz
54_S190_R2_001.fastq.gz.clean.fastq.gz
63_S98_R2_001.fastq.gz.clean.fastq.gz
81_S162_R2_001.fastq.gz.clean.fastq.gz
94_S100_R2_001.fastq.gz.clean.fastq.gz
```

Made a manifest on home computer (genohub_manifest_subset.txt)
```
nano qiime_import_sub.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/subset                    
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load QIIME2/2023.5

# Sample manifest path
MANIFEST="/data/pradalab/ceraker/16s/subset/genohub_manifest_subset.tsv"

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $MANIFEST \
  --output-path 16s_genohub_sub.qza \
  --input-format PairedEndFastqManifestPhred33V2
```

```
sbatch qiime_import_sub.sh
```
This actually made a file! Hooray!!!

Create a summary of sequences
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
#SBATCH -D /data/pradalab/ceraker/16s/subset
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh # this is used so that them function can be found in zsh and bash

module load QIIME2/2023.5

qiime demux summarize --i-data /data/pradalab/ceraker/16s/subset/16s_genohub_sub.qza --o-visualization /data/pradalab/ceraker/16s/subset/16s_genohub_sub.qzv
```
```
sbatch seq_sum.sh
```
Nope it looks like nothing

### SWITCHING GEARS
#### BLAST a subset of files

Since I can't get QIIME2 to work, and we just need to see if there are marine bacteria in the sequences, I'm going to just blast the subset of samples that I chose earlier.

Test on one file:
```
module load seqtk/1.3-GCC-9.3.0
seqtk seq -a 106_S216_R1_001.fastq.gz.clean.fastq.gz > 106_S216_R1_001.fastq.gz.clean.fasta
```
appears to have worked!

Script:
```
nano fasta.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/subset
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh # this is used so that them function can be found in zsh and bash

module load seqtk/1.3-GCC-9.3.0

for file in /data/pradalab/ceraker/16s/subset/*fastq.gz
do
seqtk seq -a $file > $file.fasta          
done
```
```
sbatch fasta.sh
```
Also appears to have worked: at least, the files exist, and look as I expect them to.

Concat all fasta files into one:
```
cat *.fasta > genohub_sub.fasta
```

Now blast my files against a database:
```
nano blast_sub.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/subset
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load BLAST+/2.13.0-gompi-2022a

blastn -query genohub_sub.fasta -db nt -remote -outfmt 7 -out genohub_sub.blastn
```

```
sbatch blast_sub.sh
```
This timed out and also the output was full of errors, which makes sense. You really shouldn't blast raw files, I was just trying a shortcut.

### BACK TO QIIME
#### Import again
_Troubleshooting:_

- Make a new manifest as a csv, then manually change it to a txt (per chat with Willow)

```
nano qiime_import.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/clean_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load QIIME2/2023.5

# Sample manifest path
MANIFEST="/data/pradalab/ceraker/16s/clean_reads/genohub_manifest_v2.txt"

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $MANIFEST \
  --output-path 16s_genohub_clean.qza \
  --input-format PairedEndFastqManifestPhred33
```
_Troubleshooting:_

```
An unexpected error has occurred:

  Each sample id can have only one forward read record in a paired-end read manifest, but the following sample ids associated with more than one forward read record: 129AND, 187AND, 213AND, 76AND

See above for debug info.
```

- _Forgot to properly deal with the "AND" files: go fix the manifest!_

```
sbatch qiime_import.sh
```
It worked!
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
#SBATCH -D /data/pradalab/ceraker/16s/clean_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh # this is used so that them function can be found in zsh and bash

module load QIIME2/2023.5

qiime demux summarize --i-data /data/pradalab/ceraker/16s/clean_reads/16s_genohub_clean.qza --o-visualization /data/pradalab/ceraker/16s/clean_reads/16s_genohub_clean.qzv
```
```
sbatch seq_sum.sh
```
Okay this file kind of looks like nothing in the QIIME viewer, but I actually think it's okay and am going to move forward.

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

- Saved directly from Google Sheets as **tsv**: _**genohub_meta_v2.tsv**_

##### DADA2
Experiment with different denoising parameters and compare statistics to choose the best one
```
nano denoise_RED_long.sh
```
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH —nodes=1 -c 36
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ceraker@uri.edu
#SBATCH --account=pradalab
#SBATCH -D /data/pradalab/ceraker/16s/clean_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh # this is used so that them function can be found in zsh and bash

module load QIIME2/2023.5

qiime dada2 denoise-paired \
--i-demultiplexed-seqs /data/pradalab/ceraker/16s/clean_reads/16s_genohub_clean.qza \
--p-trim-left-f 54 \
--p-trim-left-r 52 \
--p-trunc-len-f 270 \
--p-trunc-len-r 270 \
--p-n-threads 36 \
--o-table table_red_l.qza \
--o-representative-sequences /data/pradalab/ceraker/16s/clean_reads/rep-seqs_red_l.qza \
--o-denoising-stats /data/pradalab/ceraker/16s/clean_reads/denoising-stats_red_l.qza
```
```
sbatch denoise_RED_long.sh
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


_Troubleshooting:_

```
Plugin error from dada2:

  An error was encountered while running DADA2 in R (return code 1), please inspect stdout and stderr to learn more.

Debug info has been saved to /tmp/qiime2-q2cli-err-5c3m_1um.log
Plugin error from dada2:

  An error was encountered while running DADA2 in R (return code 1), please inspect stdout and stderr to learn more.

Debug info has been saved to /tmp/qiime2-q2cli-err-6n6xdz1m.log
```
Debug info:
```
5) Remove chimeras (method = consensus)
Error in isBimeraDenovoTable(unqs[[i]], ..., verbose = verbose) :
  Input must be a valid sequence table.
Calls: removeBimeraDenovo -> isBimeraDenovoTable
3: stop("Input must be a valid sequence table.")
2: isBimeraDenovoTable(unqs[[i]], ..., verbose = verbose)
1: removeBimeraDenovo(seqtab, method = chimeraMethod, minFoldParentOverAbundance = minParentFold,
       allowOneOff = allowOneOff, multithread = multithread)
Running external command line application(s). This may print messages to stdout and/or stderr.
The command(s) being run are below. These commands cannot be manually re-run as they will depend on temporary files that no longer exist.

Command: run_dada.R --input_directory /tmp/tmpjiu3ru_7/forward --input_directory_reverse /tmp/tmpjiu3ru_7/reverse --output_path /tmp/tmpjiu3ru_7/output.tsv.biom --output_track /tmp/tmpjiu3ru_7/track.tsv --filtered_directory /tmp/tmpjiu3ru_7/filt_f --filtered_directory_reverse /tmp/tmpjiu3ru_7/filt_r --truncation_length 270 --truncation_length_reverse 270 --trim_left 16 --trim_left_reverse 16 --max_expected_errors 2.0 --max_expected_errors_reverse 2.0 --truncation_quality_score 2 --min_overlap 12 --pooling_method independent --chimera_method consensus --min_parental_fold 1.0 --allow_one_off False --num_threads 1 --learn_min_reads 1000000

Traceback (most recent call last):
  File "/opt/software/QIIME2/2023.5/lib/python3.8/site-packages/q2_dada2/_denoise.py", line 326, in denoise_paired
    run_commands([cmd])
  File "/opt/software/QIIME2/2023.5/lib/python3.8/site-packages/q2_dada2/_denoise.py", line 36, in run_commands
    subprocess.run(cmd, check=True)
  File "/opt/software/QIIME2/2023.5/lib/python3.8/subprocess.py", line 516, in run
    raise CalledProcessError(retcode, process.args,
subprocess.CalledProcessError: Command '['run_dada.R', '--input_directory', '/tmp/tmpjiu3ru_7/forward', '--input_directory_reverse', '/tmp/tmpjiu3ru_7/reverse', '--output_path', '/tmp/tmpjiu3ru_7/output.tsv.biom', '--output_track', '/tmp/tmpjiu3ru_7/track.tsv', '--filtered_directory', '/tmp/tmpjiu3ru_7/filt_f', '--filtered_directory_reverse', '/tmp/tmpjiu3ru_7/filt_r', '--truncation_length', '270', '--truncation_length_reverse', '270', '--trim_left', '16', '--trim_left_reverse', '16', '--max_expected_errors', '2.0', '--max_expected_errors_reverse', '2.0', '--truncation_quality_score', '2', '--min_overlap', '12', '--pooling_method', 'independent', '--chimera_method', 'consensus', '--min_parental_fold', '1.0', '--allow_one_off', 'False', '--num_threads', '1', '--learn_min_reads', '1000000']' returned non-zero exit status 1.

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/opt/software/QIIME2/2023.5/lib/python3.8/site-packages/q2cli/commands.py", line 468, in __call__
    results = action(**arguments)
  File "<decorator-gen-72>", line 2, in denoise_paired
  File "/opt/software/QIIME2/2023.5/lib/python3.8/site-packages/qiime2/sdk/action.py", line 274, in bound_callable
    outputs = self._callable_executor_(
  File "/opt/software/QIIME2/2023.5/lib/python3.8/site-packages/qiime2/sdk/action.py", line 509, in _callable_executor_
    output_views = self._callable(**view_args)
  File "/opt/software/QIIME2/2023.5/lib/python3.8/site-packages/q2_dada2/_denoise.py", line 339, in denoise_paired
    raise Exception("An error was encountered while running DADA2"
Exception: An error was encountered while running DADA2 in R (return code 1), please inspect stdout and stderr to learn more.
```

##### Try denoising with deblur instead of dada2

Deblur tutorial from: https://telatin.github.io/microbiome-bioinformatics/Metabarcoding-deblur/

Deblur can only deal with single end reads, so merge reads using VSEARCH
```
nano vsearch.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/clean_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh # this is used so that them function can be found in zsh and bash

module load QIIME2/2023.5

qiime vsearch merge-pairs \
     --i-demultiplexed-seqs /data/pradalab/ceraker/16s/clean_reads/16s_genohub_clean.qza \
     --o-merged-sequences /data/pradalab/ceraker/16s/clean_reads/joined-reads.qza \
     --p-threads 8
```
```
sbatch vsearch.sh
```
Sequence quality control and feature table construction: filter sequences based on quality scores and the presence of ambiguous base calls.

```
nano qualfilt.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/clean_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh # this is used so that them function can be found in zsh and bash

module load QIIME2/2023.5
qiime quality-filter q-score \
      --i-demux /data/pradalab/ceraker/16s/clean_reads/joined-reads.qza \
      --o-filtered-sequences /data/pradalab/ceraker/16s/clean_reads/joined-filtered.qza \
      --o-filter-stats /data/pradalab/ceraker/16s/clean_reads/joined-filter-stats.qza
```
```
sbatch qualfilt.sh
```

Now denoise with Deblur
```
nano deblur_RED_long.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/clean_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh # this is used so that them function can be found in zsh and bash

module load QIIME2/2023.5

qiime deblur denoise-16S \
--i-demultiplexed-seqs /data/pradalab/ceraker/16s/clean_reads/joined-filtered.qza \
--p-trim-length 270 \
--o-table table_red_l.qza \
--o-representative-sequences /data/pradalab/ceraker/16s/clean_reads/rep-seqs_red_l.qza \
--o-stats /data/pradalab/ceraker/16s/clean_reads/denoising-stats_red_l.qza
```

_Troubleshooting:_

```
Plugin error from deblur:

  No sequences passed the filter. It is possible the trim_length (270) may exceed the longest sequence, that all of $

Debug info has been saved to /tmp/qiime2-q2cli-err-mnmlsyzl.log
```

##### Try with raw files
```
cd og_reads
ls > rawfiles.txt
```
Use this txt file to create genohub_manifest_raw.txt (see earlier strategy)

Can use the same metadata file (genohub_meta.txt), as the samples are the same.

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

"Denoise" cleaned files: no additional filtering, just creating the output files.
```
nano denoise_clean.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/clean_reads
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh # this is used so that them function can be found in zsh and bash

module load QIIME2/2023.5

qiime dada2 denoise-paired \
--i-demultiplexed-seqs /data/pradalab/ceraker/16s/clean_reads/16s_genohub_clean.qza \
--p-n-threads 36 \
--p-trunc-len-f 300 \
--p-trunc-len-r 300 \
--o-table table_red_l.qza \
--o-representative-sequences /data/pradalab/ceraker/16s/clean_reads/rep-seqs_red_l.qza \
--o-denoising-stats /data/pradalab/ceraker/16s/clean_reads/denoising-stats_red_l.qza
```
```
sbatch denoise_clean.sh
```
Job ID: 289186, n074
Plugin error from dada2:

  An error was encountered while running DADA2 in R (return code 1), please inspect stdout and stderr to learn more.

Debug info has been saved to /tmp/qiime2-q2cli-err-a4r4tvve.log

##### DADA2 with raw reads
Experiment with different denoising parameters and compare statistics to choose the best one
```
nano denoise_RED_long.sh
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

qiime dada2 denoise-paired \
--i-demultiplexed-seqs /data/pradalab/ceraker/16s/og_reads/16s_genohub_raw.qza \
--p-trim-left-f 54 \
--p-trim-left-r 52 \
--p-trunc-len-f 280 \
--p-trunc-len-r 280 \
--p-n-threads 36 \
--o-table table_red_l.qza \
--o-representative-sequences /data/pradalab/ceraker/16s/clean_reads/rep-seqs_red_l.qza \
--o-denoising-stats /data/pradalab/ceraker/16s/clean_reads/denoising-stats_red_l.qza
```
```
sbatch denoise_RED_long.sh
```

### Clustering
Clutering with files cleaned in fastp and not denoised, and raw files NOT cleaned in fastp that were later denoised.

Cluster raw files that got denoised, using RED_short output files
```
nano cluster_og.sh
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
METADATA="/data/pradalab/ceraker/16s/og_reads/geno_meta_v2.txt"

# Sample manifest path
MANIFEST="/data/pradalab/ceraker/16s/og_reads/genohub_manifest_raw.txt"

#### CLUSTERING

# Summarize feature table and sequences
qiime metadata tabulate \
  --m-input-file denoising-stats_red_s.qza \
  --o-visualization denoising-stats_red_s.qzv
qiime feature-table summarize \
  --i-table table_red_s.qza \
  --o-visualization table_red_s.qzv \
  --m-sample-metadata-file $METADATA
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_red_s.qza \
  --o-visualization rep-seqs_red_s.qzv
```
```
sbatch cluster_og.sh
```

### Classifying
Attempting classifying with files cleaned in fastp and not denoised, and raw files NOT cleaned in fastp that were later denoised.

##### Download classifier from QIIME2
```
wget https://data.qiime2.org/2023.5/common/silva-138-99-515-806-nb-classifier.qza
```

```
nano taxonomy_og.sh
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

#### TAXONOMY CLASSIFICATION

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs_red_s.qza \
  --o-classification taxonomy_og.qza

qiime taxa filter-table \
     --i-table table_red_s.qza \
     --i-taxonomy taxonomy_og.qza \
     --p-mode contains \
     --p-exclude "Unassigned","Chloroplast","Eukaryota" \
     --o-filtered-table table-filtered_og.qza

qiime metadata tabulate \
    --m-input-file taxonomy_og.qza \
    --o-visualization taxonomy_og.qzv
qiime taxa barplot \
    --i-table table-filtered_og.qza \
    --i-taxonomy taxonomy_og.qza \
    --m-metadata-file $METADATA \
    --o-visualization taxa-bar-plots-filtered_og.qzv
qiime metadata tabulate \
    --m-input-file rep-seqs_red_s.qza \
    --m-input-file taxonomy_og.qza \
    --o-visualization tabulated-feature-metadata_og.qzv

qiime feature-table summarize \
    --i-table table-filtered_og.qza \
    --o-visualization table-filtered_og.qzv \
    --m-sample-metadata-file $METADATA

#### CREATES PHYLOGENETIC TREES

# align and mask sequences
qiime alignment mafft \
  --i-sequences rep-seqs_red_s.qza \
  --o-alignment aligned-rep-seqs_og.qza
qiime alignment mask \
  --i-alignment aligned-rep-seqs_og.qza \
  --o-masked-alignment masked-aligned-rep-seqs_og.qza

# calculate tree
qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs_og.qza \
  --o-tree unrooted-tree_og.qza
qiime phylogeny midpoint-root \
  --i-tree unrooted-tree_og.qza \
  --o-rooted-tree rooted-tree_og.qza
```
```
sbatch taxonomy_og.sh
```

```
nano denoise_vis.sh
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

qiime metadata tabulate \
  --m-input-file denoising-stats_red_l.qza \
  --o-visualization denoising-stats_red_l.qzv
```
```
sbatch denoise_vis.sh
```

OTHER DENOISE PARAMETERS
```
nano denoise_YELLOW_long.sh
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
--p-trunc-len-f 260 \
--p-trunc-len-r 260 \
--p-n-threads 36 \
--o-table table_yellow_l.qza \
--o-representative-sequences /data/pradalab/ceraker/16s/og_reads/rep-seqs_yellow_l.qza \
--o-denoising-stats /data/pradalab/ceraker/16s/og_reads/denoising-stats_yellow_l.qza

qiime metadata tabulate \
  --m-input-file denoising-stats_yellow_l.qza \
  --o-visualization denoising-stats_yellow_l.qzv
```
```
sbatch denoise_YELLOW_long.sh
```

```
nano denoise_YELLOW_short.sh
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
--p-trim-left-f 16 \
--p-trim-left-r 16 \
--p-trunc-len-f 260 \
--p-trunc-len-r 260 \
--p-n-threads 36 \
--o-table table_yellow_s.qza \
--o-representative-sequences /data/pradalab/ceraker/16s/og_reads/rep-seqs_yellow_s.qza \
--o-denoising-stats /data/pradalab/ceraker/16s/og_reads/denoising-stats_yellow_s.qza

qiime metadata tabulate \
  --m-input-file denoising-stats_yellow_s.qza \
  --o-visualization denoising-stats_yellow_s.qzv
```
```
sbatch denoise_YELLOW_short.sh
```

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

### SAMPLE SUMMARY

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


### Clustering with chosen denoising parameters
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

### Classifying
Classifying samples denoised with blue_long.

Classifier already downloaded from QIIME2

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


### Diversity Analyses
##### Perform initial diversity analyses in QIIME2

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
Copy all .qzv files to location outside Andromeda to visualize and export stats
