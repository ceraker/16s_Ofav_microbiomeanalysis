# __*Orbicella faveolata* 16S microbiome analysis__
#### Author: Cassie Raker
#### Last updated: December 29, 2023

#### Sequencing performed by Genohub
#### Data uploaded and analyzed on Andromeda
#### Workflow modified from [Dr Ariana Huffmyer](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/16S/Mothur_bioinformatics.md) and [Dr Emma Strand](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2022-02-17-Holobiont-Integration-16S-Mothur-Pipeline.md#holobiont-integration-16s-mothur-pipeline)

### Andromeda for URI

Information for creating an account on [Andromeda](https://web.uri.edu/hpc-research-computing/using-andromeda/).

###### Accessing Andromeda account
```
ssh -i ~/andromeda_keys ceraker@ssh3.hac.uri.edu
```

General Workflow:  
1. [Prepare directory](#Directory)   
2. [Start mothur](#Start)    
3. [Preparing sequences](#Prepare)   
4. [QC sequences](#QC)    
5. [Unique sequences](#Unique)    
6. [Aligning](#Align)    
7. [Preclustering](#Precluster)    
8. [Identify chimeras](#Chimera)  
9. [Classify sequences](#Classify)     
10. [Cluster OTU's](#Cluster)    
11. [Subsampling](#Subsample)  
12. [Calculate ecological statistics](#Statistics)  
13. [Output data for R analysis](#Output)   

### <a name="Directory"></a> **1. Prepare Directory**
Navigate to data folder
```
cd /data/pradalab/ceraker/16s
```
Make a directory
```
mkdir mothur
cd mothur
mkdir scripts
mkdir processed_data
```
Activate previously created conda environment
```
conda activate 16s
```
Create symbolic links to sequencing files (to save space)
```
ln -s /data/pradalab/ceraker/16s/og_reads/*.gz .
```
Primers used:

16S-V3V4-Bact-For:
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG

16S-V3V4-Bact-Rev:
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC

Create a text file with this sequence information
```
nano oligos.oligos

## copy and paste the following text into that file

primer TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC
```

The current mothur module available is: Mothur/1.46.1-foss-2020b

As of 11 December 2023, this is the current version of Mothur available.  

### <a name="Start"></a> **2. Start Mothur**  

Start in interactive mode to see how we can get mothur to run and test that the module is working.  

```
interactive
module load Mothur/1.46.1-foss-2020b
cd mothur/
mothur
```

This successfully activated mothur. This should see a display of citation, program information, and `mothur >`.  

At any time, use `quit` to leave mothur.  

### <a name="Prepare"></a> **3. Preparing Sequences: make.file, make.contig, and summary.seq**  

These can be run separately, but it's easy to combine them into one script.
```
cd scripts
nano contigs.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/mothur
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load Mothur/1.46.1-foss-2020b

mothur "#make.file(inputdir=., type=gz, prefix=ofav)"

mothur "#make.contigs(inputdir=., outputdir=., file=ofav.files, trimoverlap=T, oligos=oligos.oligos, pdiffs=2, checkorient=t)"

mothur "#summary.seqs(fasta=ofav.trim.contigs.fasta)"
```
```
cd .. #have to be in the mothur directory when running script
sbatch scripts/contigs.sh
```
This failed out, with no error message, and all files seem to be empty. I'm going to run each step individually.

- ```make.sh```: COMPLETED, ~2 minutes, outputs a three column file ```ofav.files```
```
100     100_S46_R1_001.fastq.gz 100_S46_R2_001.fastq.gz
100_S86 100_S86_R1_001.fastq.gz 100_S86_R2_001.fastq.gz
101     101_S28_R1_001.fastq.gz 101_S28_R2_001.fastq.gz
101_S68 101_S68_R1_001.fastq.gz 101_S68_R2_001.fastq.gz
102     102_S46_R1_001.fastq.gz 102_S46_R2_001.fastq.gz
102_S6  102_S6_R1_001.fastq.gz  102_S6_R2_001.fastq.gz
```
- ```contig.sh```: COMPLETED, 1.5 hours, but ```ofav.trim.contigs.fasta``` is empty. All output files:
  - ```ofav.trim.contigs.fasta```
  - ```ofav.scrap.contigs.fasta```
  - ```ofav.contigs.report```
  - ```mothur.1703879368.logfile```

Head of ```ofav.scrap.contigs.fasta``` is:
```
>VH01338_43_AAF7L7LM5_1_1101_42166_1076 | f(f)(bf)	ee=0.670777	fbdiffs=1000(noMatch), rbdiffs=1000(noMatch) fpdiffs=36(noMatch), rpdiffs=1002(noMatch)
CNTACGGGCGGCAGCAGGTCGGAAGAGCGTCGTTTGGGATTAGATACCCGTGTAGTC
>VH01338_43_AAF7L7LM5_1_1101_38511_1057 | f(f)(bf)	ee=0.079656	fbdiffs=1000(noMatch), rbdiffs=1000(noMatch) fpdiffs=72(noMatch), rpdiffs=1002(noMatch)
CCTACGGGGGGCAGCAGCGGCAGCAGCGGCAGCAGCGGCAGCAGCGGTAGCAGCGGTAGCAGCGGTAGTAGCGGTAGTAGCGGTAGTAGCGGTAGTGGCGGTAGTAGCGGTAGTAGCGGTAGTAGCGGTAGTAGCAGTAGGGTTAGGGTTTAGGGATTAGATACCCGGGTAGTC
>VH01338_43_AAF7L7LM5_1_1101_34573_1095 | f(f)(bf)	ee=0.0955442	fbdiffs=1000(noMatch), rbdiffs=1000(noMatch) fpdiffs=89(noMatch), rpdiffs=1002(noMatch)
GGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAAGGGAGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGTGTGTAGGCGGATTATCAAGTTAGGGGTGAAATCCCGGGGCTCAACCTCGGCACTGCCTTTAAGGCTGATAATCTAGAGTATGTGAGGG
>VH01338_43_AAF7L7LM5_1_1101_46256_1076 | f(f)(bf)	ee=0.67077	fbdiffs=1000(noMatch), rbdiffs=1000(noMatch) fpdiffs=66(noMatch), rpdiffs=1002(noMatch)
CNTACGGGAGGCAGCAGCGTCAGATGTGGATTAGATACCCCAGTAGTC
>VH01338_43_AAF7L7LM5_1_1101_40859_1133 | f(f)(bf)	ee=0.680791	fbdiffs=1000(noMatch), rbdiffs=1000(noMatch) fpdiffs=27(noMatch), rpdiffs=1002(noMatch)
GGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTNTCCGGAATAATTGGGCGTAAAGCGTTCGTAGGTGGTTTTGTAAGTCTACTGTTAAAGCGTGTAGCTCAACTACATATAGGCAGTGGAAACTACAAGACTTGAGTGCGTTCGGG
```

_**Troubleshooting:**_

Mothur isn't finding any matches between the primers and the sequences, therefore it is filtering all of them out and the ```ofav.trim.contigs.fasta``` file is empty. Possible that GenoHub already trimmed the primers?

Changes in ```contig.sh``` to try to fix this issue:
- ```trimoverlap=F```: trim file is still empty
- assemble without checking for primers (i.e. no ```oligos``` argument, and change ```trimoverlap``` back to ```T```): 1 hr 40 min, sequences are in ```ofav.trim.contigs.fasta``` file!

Head of ```ofav.trim.contigs.fasta``` is:
```
>VH01338_43_AAF7L7LM5_1_1101_38511_1057	ee=0.079656
CCTACGGGGGGCAGCAGCGGCAGCAGCGGCAGCAGCGGCAGCAGCGGTAGCAGCGGTAGCAGCGGTAGTAGCGGTAGTAGCGGTAGTAGCGGTAGTGGCGGTAGTAGCGGTAGTAGCGGTAGTAGCGGTAGTAGCAGTAGGGTTAGGGTTTAGGGATTAGATACCCGGGTAGTC
>VH01338_43_AAF7L7LM5_1_1101_34573_1095	ee=0.0955442
GGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAAGGGAGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGTGTGTAGGCGGATTATCAAGTTAGGGGTGAAATCCCGGGGCTCAACCTCGGCACTGCCTTTAAGGCTGATAATCTAGAGTATGTGAGGG
>VH01338_43_AAF7L7LM5_1_1101_40859_1133	ee=0.680791
GGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTNTCCGGAATAATTGGGCGTAAAGCGTTCGTAGGTGGTTTTGTAAGTCTACTGTTAAAGCGTGTAGCTCAACTACATATAGGCAGTGGAAACTACAAGACTTGAGTGCGTTCGGG
>VH01338_43_AAF7L7LM5_1_1101_64794_1170	ee=0.00361062
GGCTAACTCCGTGCCAGCAACCGCGGTAAGACGGAGGGTGCAAACGTTGTTCGGAATCACTGGGCATAAAGGGCACGTAGGCGGTCAACTAAGTCAGATGTGAAAGCCCCCGGCTCAACCGGGGAACGGCATTTGATACT
>VH01338_43_AAF7L7LM5_1_1101_63695_1246	ee=0.0134992
GGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAAGGGAGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGTGTAGGCGGATTATCAAGTTAGGGGTGAAATCCCGGGGCTCAACCTCGGCACTGCCTTTAAGACTGATAATCTAGAGTATGTGAGGG
```
Head of ```ofav.scrap.contigs.fasta``` is:
```
>VH01338_43_AAF7L7LM5_1_1101_40443_1057 | l	ee=2.52383
NNNN
>VH01338_43_AAF7L7LM5_1_1101_50327_1095 | l	ee=2.52383
NNNN
>A00572_431_HMCWKDRX2_2_2108_25843_23234 | l	ee=2.52383
NNNN
>VH01338_43_AAF7L7LM5_1_1101_26336_28716 | l	ee=2.52383
NNNN
>A00572_431_HMCWKDRX2_2_2116_23204_26271 | l	ee=2.52383
NNNN
```

**_Now that contigs have been assembled, run summary script_**

```
nano seqsum.sh
sbatch seqsum.sh
```
Took ~20 minutes.

Head of ```ofav.trim.contigs.summary```:
```
VH01338_43_AAF7L7LM5_1_1101_38511_1057	1	174	174	0	6	1
VH01338_43_AAF7L7LM5_1_1101_34573_1095	1	162	162	0	4	1
VH01338_43_AAF7L7LM5_1_1101_40859_1133	1	162	162	1	4	1
VH01338_43_AAF7L7LM5_1_1101_64794_1170	1	140	140	0	5	1
VH01338_43_AAF7L7LM5_1_1101_63695_1246	1	162	162	0	4	1
VH01338_43_AAF7L7LM5_1_1101_41673_1265	1	162	162	0	4	1
```
Count number of "good" and number of "bad" sequences:
```
grep -c "^>" ofav.trim.contigs.fasta #47865613
grep -c "^>" ofav.scrap.contigs.fasta #2601
```
**47865613** sequences were saved and **2601** were removed. These may not be the best parameters, but I'm going to move ahead for now.

Summary table from ```ofav.trim.contigs.summary```
```
  Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	1	1	0	1	1
2.5%-tile:  1 2 2 0 1 1196641
25%-tile:	1	42	42	0	3	11966404
Median:         1	138     138     0	4	23932807
75%-tile:	1	159     159     0	4	35899210
97.5%-tile:     1	216     216     16	8	46668973
Maximum:        1	310     310     126     190     47865613
Mean:   1	110     110     1	3
# of Seqs:	47865613
```
_(fix this formatting later)_

### <a name="QC"></a> **4. QC'ing sequences with screen.seqs**
In the `screen.seqs()` function, we will specify the fasta file of the contigs generated in the previous step and remove any sequence with an ambiguous call ("N"). We will also remove sequences >350 nt. We will also set a minimum size amount (200). These parameters could be adjusted based on specific experiment and variable region.

```
nano screen.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/mothur
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#screen.seqs(inputdir=., outputdir=., fasta=ofav.trim.contigs.fasta, group=ofav.contigs.groups, maxambig=0, maxlength=350, minlength=200)"

mothur "#summary.seqs(fasta=ofav.trim.contigs.good.fasta)"
```
```
sbatch screen.sh
```
Count number of "good" and number of "bad" sequences:
```
grep -c "^>" ofav.trim.contigs.good.fasta #1503388
grep -c "^>" ofav.trim.contigs.bad.accnos #0
```
These parameters did not filter out any sequences.
```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	200     200     0	2	1
2.5%-tile:	1	206     206     0	4	37585
25%-tile:	1	215     215     0	6	375848
Median:         1	288     288     0	8	751695
75%-tile:	1	289     289     0	8	1127542
97.5%-tile:     1	298     298     0	9	1465804
Maximum:        1	301     301     0	189     1503388
Mean:   1	260     260     0	7
# of Seqs:	1503388
```
_(fix this formatting later)_

### <a name="Unique"></a> **5. Determining and counting unique sequences**  

Next, determine the number of unique sequences.
```
nano unique.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/mothur
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#unique.seqs(fasta=ofav.trim.contigs.good.fasta)"

mothur "#count.seqs(name=ofav.trim.contigs.good.names, group=ofav.contigs.good.groups)"

mothur "#summary.seqs(fasta=ofav.trim.contigs.good.unique.fasta, count=ofav.trim.contigs.good.count_table)"
```
```
sbatch unique.sh
```

After determining the unique sequences, we can use these unique sequences to count how many times each sequence (which will later be classified to OTU) shows up in each sample.  

The whole script generate the following outputs:

```
ofav.trim.contigs.good.names
ofav.trim.contigs.good.unique.fasta
ofav.trim.contigs.good.unique.summary
ofav.trim.contigs.good.count_table
```
Summary table:
```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	200     200     0	2	1
2.5%-tile:	1	206     206     0	4	37585
25%-tile:	1	215     215     0	6	375848
Median:         1	288     288     0	8	751695
75%-tile:	1	289     289     0	8	1127542
97.5%-tile:     1       298     298     0       9       1465804
Maximum:        1	301     301     0	189     1503388
Mean:   1       260     260     0 	7
# of unique seqs:	773970
total # of seqs:  	1503388
```
_(fix this formatting later)_

Now we can align just the unique sequences, which will be much faster than aligning the full data set and is an indicator of how polished and clean the data are.  

*From this, we have our unique sequences identified and can proceed with further cleaning and polishing of the data. Next we will look at alignment, error rate, chimeras, classification and further analysis.*  

### <a name="Align"></a> **6. Aligning to reference database**

#### Prepare the reference sequences  

First, download the Silva reference files from the [Mothur Wiki](https://mothur.org/wiki/silva_reference_files/) at the latest release.

The silva reference is used and recommended by the Mothur team. It is a manually curated data base with high diversity and high alignment quality.  

In your mothur directory, run `wget` to download the silva reference and training sets from the mothur website and unzip them and move into the right directory.

```
wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.bacteria.zip

wget https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip

unzip silva.bacteria.zip
cp silva.bacteria.fasta ../silva.bacteria.fasta

unzip trainset9_032012.pds.zip
```
Now we have the reference files that we need in our directory.  

Now we are going to take our Silva database reference alignment and select the V4 region.  

We can do this with the `pcr.seqs` function.

This region includes the primers, even though our sequences dont include the primers. If we did keepdots=T then the first hundreds of columns would be periods, which is not useful for us. These periods are placeholders in the silva database.   

```
nano silva_ref.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/mothur
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F)"

mothur "#summary.seqs(fasta=silva.bacteria.pcr.fasta)"

mothur "#rename.file(input=silva.bacteria.pcr.fasta, new=silva.v4.fasta)"
```
```
sbatch silva_ref.sh
```
The script outputs the following file:

```
silva.bacteria.pcr.fasta
```
Which was then renamed to `silva.v4.fasta` as specified in the script.

We now have a reference to align to.  

#### Align sequences to the reference  

We next align sequences with the `align.seqs` command.  

Write a script to do this, generate a new summary of the output, and run.  

```
nano align.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/mothur
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#align.seqs(fasta=ofav.trim.contigs.good.unique.fasta, reference=silva.v4.fasta)"

mothur "#summary.seqs(fasta=ofav.trim.contigs.good.unique.align)"
```
```
sbatch align.sh
```
The script generates these output files:

```
Output File Names:
ofav.trim.contigs.good.unique.align
ofav.trim.contigs.good.unique.align.report
ofav.trim.contigs.good.unique.flip.accnos
```

View the report file.  

```
head ofav.trim.contigs.good.unique.align.report

#output
QueryName	QueryLength	TemplateName	TemplateLength	SearchMethod	SearchScore	AlignmentMethod	QueryStart	QueryEnd	TemplateStart	TemplateEnd	PairwiseAlignmentLength	GapsInQuery	GapsInTemplate	LongestInsert	SimBtwnQuery&Template
VH01338_43_AAF7L7LM5_1_2104_52959_12170	289	D38625.1	293	kmer	2.12	needleman	1	17	277	293	17	0	0	0	52.94
VH01338_43_AAF7L7LM5_1_1404_47884_50601	214	EF574945.1	293	kmer	2.89	needleman	1	43	251	293	43	0	0	0	51.16
VH01338_43_AAF7L7LM5_1_2102_49626_18569	237	X84241.1	294	kmer	3.47	needleman	217	237	1	21	21	0	0	0	57.14
VH01338_43_AAF7L7LM5_1_1201_75682_25157	289	DQ181669.1	293	kmer	2.48	needleman	1	43	251	293	43	0	0	0	51.16
VH01338_43_AAF7L7LM5_1_1101_53849_27486	214	AF387306.1	293	kmer	2.41	needleman	202	214	1	13	13	0	0	0	61.53
VH01338_43_AAF7L7LM5_1_2403_24404_7437	293	AJ292602.1	293	kmer	3.14	needleman	267	293	1	27	27	0	0	0	55.55
VH01338_43_AAF7L7LM5_1_1202_71573_20935	298	AJ278167.1	293	kmer	33.33	needleman	8	298	1	291	291	0	0	0	82.81
VH01338_43_AAF7L7LM5_1_2304_66479_29738	289	D38625.1	293	kmer	2.12	needleman	1	5	289	293	5	0	0	0	60.00
VH01338_43_AAF7L7LM5_1_1203_12456_48614	215	AF387306.1	293	kmer	2.40	needleman	1	17	277	293	17	0	0	0	52.94

```

View the accnos file.  

```
head ofav.trim.contigs.good.unique.flip.accnos

#output
VH01338_43_AAF7L7LM5_1_2104_52959_12170
VH01338_43_AAF7L7LM5_1_1404_47884_50601
VH01338_43_AAF7L7LM5_1_2102_49626_18569
VH01338_43_AAF7L7LM5_1_1201_75682_25157
VH01338_43_AAF7L7LM5_1_1101_53849_27486
VH01338_43_AAF7L7LM5_1_2403_24404_7437
VH01338_43_AAF7L7LM5_1_2304_66479_29738
VH01338_43_AAF7L7LM5_1_1203_12456_48614
VH01338_43_AAF7L7LM5_1_2103_23723_19175
VH01338_43_AAF7L7LM5_1_1202_25635_3518
```
The summary now looks like this (in `output_script`):
```

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        0	0	0	0	1	1
2.5%-tile:	1	1253    5	0	1	19350
25%-tile:	1	1968    7	0	1	193493
Median:         10657   13425   17	0	3	386986
75%-tile:	13399   13425   21	0	3	580478
97.5%-tile:     13404   13425   291     0	6	754621
Maximum:        13425   13425   296     0	10	773970
Mean:   6962    8862    32	0	2
# of Seqs:	773970
```
_(fix this formatting later)_

From this, we see that 773,970 sequences aligned to the reference, which matches the number of sequences that we had after the unique.sh step (773,970). In the next steps we will filter out any sequences that don't meet alignment settings.  

#### QC sequences according to alignment to the reference  

Our sequences now align at the correct positions on the reference.

Now remove sequences that are outside the alignment window (1968-11550bp). This removes anything that starts after `start` and ends before `end`. Maxhomop=8 argument removes anything that has repeats greater than the threshold - e.g., 8 A's in a row = polymer 8. Here we will removes polymers >8 because we are confident these are likely not high quality sequences (see mothur MiSeq SOP for more information).  

```
nano screen2.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/mothur
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#screen.seqs(fasta=ofav.trim.contigs.good.unique.align, count=ofav.trim.contigs.good.count_table, start=1968, end=11550, maxhomop=8)"

mothur "#summary.seqs(fasta=ofav.trim.contigs.good.unique.good.align, count=ofav.trim.contigs.good.good.count_table)"
```
```
sbatch screen2.sh
```
This will output the following files:
```
Output File Names:
ofav.trim.contigs.good.unique.good.align
ofav.trim.contigs.good.unique.bad.accnos
ofav.trim.contigs.good.good.count_table
```

View the accnos file to see why sequences will be removed and count the number of "bad" sequences.  

```
head ofav.trim.contigs.good.unique.bad.accnos

#output
VH01338_43_AAF7L7LM5_1_1303_24404_26709	start
VH01338_43_AAF7L7LM5_1_1303_29233_24002	end
VH01338_43_AAF7L7LM5_1_2103_66668_23756	start|<length
VH01338_43_AAF7L7LM5_1_1302_35216_25138	end
VH01338_43_AAF7L7LM5_1_2302_33588_20746	start
VH01338_43_AAF7L7LM5_1_1402_65116_53858	start
VH01338_43_AAF7L7LM5_1_1204_31524_20273	end
VH01338_43_AAF7L7LM5_1_1401_24120_3972	end
VH01338_43_AAF7L7LM5_1_1101_38170_23510	start
VH01338_43_AAF7L7LM5_1_2204_25938_18474	end
```
```
grep -c ".*" ofav.trim.contigs.good.unique.bad.accnos
```

729,440 uniques are tagged to be removed due to filtering at this step.

Summary table:
```

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	11550   271     0	3	1
2.5%-tile:	1	13422   290     0	4	1449
25%-tile:	1	13422   291     0	4	14485
Median:         1	13422   291     0	4	28969
75%-tile:	1	13422   291     0	5	43453
97.5%-tile:     1	13422   292     0	7	56489
Maximum:        1968    13425   296     0	8	57937
Mean:   15	13421   290     0	4
# of unique seqs:	44530
total # of seqs:        57937
```

#### Filter sequences  

Now we can filter out sequences that didn't meet our criteria above, which will generate a report and a new summary of our sequences.  

We will run the following code. We align vertically and use trump=. to align the sequences accounting for periods in the reference.   

```
filter.seqs(fasta=mcap.trim.contigs.good.unique.good.align, vertical=T, trump=.)
```

Write and run a script.  

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
#SBATCH -D /data/pradalab/ceraker/16s/mothur
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#filter.seqs(fasta=ofav.trim.contigs.good.unique.good.align, vertical=T, trump=.)"

mothur "#summary.seqs(fasta=ofav.trim.contigs.good.unique.good.filter.fasta, count=ofav.trim.contigs.good.good.count_table)"
```
```
sbatch filter.sh
```
This script outputs the following files:  

```
Output File Names:
ofav.filter
ofav.trim.contigs.good.unique.good.filter.fasta
```

We get a report on the filtering in the `output_script` file that looks like this:

```
Length of filtered alignment: 451
Number of columns removed: 12974
Length of the original alignment: 13425
Number of sequences used to construct filter: 44530
```

We also get a new summary that looks like this:  
```

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	450     241     0	3	1
2.5%-tile:	1	451     252     0	4	1449
25%-tile:	1	451     253     0	4	14485
Median:         1	451     253     0	4	28969
75%-tile:	1	451     253     0	5	43453
97.5%-tile:     1	451     254     0	7	56489
Maximum:        2	451     268     0	8	57937
Mean:   1	450     252     0	4
# of unique seqs:	44530
total # of seqs:        57937
```
From this summary we see that the alignment window spans ~450 bp and the length of our sequences is about 253 nt. We have a maximum polymer of 8 as specified in our settings above.  

### <a name="Precluster"></a> **7. Polish the data with pre clustering**     

Now we need to further polish and cluster the data with pre.cluster. The purpose of this step is to remove noise due to sequencing error. The rational behind this step assumes that the most abundant sequences are the most trustworthy and likely do not have sequencing errors. Pre-clustering then looks at the relationship between abundant and rare sequences - rare sequences that are "close" (e.g., 1 nt difference) to highly abundant sequences are likely due to sequencing error. This step will pool sequences and look at the maximum differences between sequences within this group to form ASV groupings.

In this step, the number of sequences is not reduced, but they are grouped into amplicon sequence variants ASV's which reduces the error rate. V4 region has the lowest likelihood of errors, so the error rate is going to be lower than for other variable regions.  

Other programs that conduct this "denoising" are DADA2, UNOISE, and DEBLUR. However, these programs remove the rare sequences, which can distort the relative abundance of remaining sequences. DADA2 also removes all sigletons (sequences with single representation) which disproportionately affects the sequence relative abundance. Mothur avoids the removal of rare sequences for this reason.

We will first add code to identify unique sequences after the filtering steps above.  

We will then perform the pre-clustering a default of 1 nt difference. Diffs can be changed according to your requirements.  

Finally, we will run another summary.  

Write and run the script.  

```
nano precluster.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/mothur
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#unique.seqs(fasta=ofav.trim.contigs.good.unique.good.filter.fasta, count=ofav.trim.contigs.good.good.count_table)"

mothur "#pre.cluster(fasta=ofav.trim.contigs.good.unique.good.filter.unique.fasta, count=ofav.trim.contigs.good.unique.good.filter.count_table, diffs=1)"

mothur "#summary.seqs(fasta=ofav.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=ofav.trim.contigs.good.unique.good.filter.unique.precluster.count_table)"
```
```
sbatch precluster.sh
```
Mothur output the following files:   

```
Output File Names:
ofav.trim.contigs.good.unique.good.filter.count_table
ofav.trim.contigs.good.unique.good.filter.unique.fasta
```
Then, the pre-clustering step outputs a lot of files, but the two most important are:

```
ofav.trim.contigs.good.unique.good.filter.unique.precluster.fasta
ofav.trim.contigs.good.unique.good.filter.unique.precluster.count_table
```

The other files have text for maps of sequence name, errors, abundance, differences, and the filtered sequence for each sample.   

Finally, we get the output from the summary:   
```
Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	450     241     0	3	1
2.5%-tile:	1	451     252     0	4	1449
25%-tile:	1	451     253     0	4	14485
Median:         1	451     253     0	4	28969
75%-tile:	1	451     253     0	5	43453
97.5%-tile:     1	451     254     0	7	56489
Maximum:        2	451     268     0	8	57937
Mean:   1	450     252     0	4
# of unique seqs:	16622
total # of seqs:        57937
```
Note that the number of unique sequences has decreased from 44,530 to 16,622 as expected since we are clustering sequences that are within 1 nt difference from each other.

### <a name="Chimera"></a> **8. Identify chimeras**  

Now we will remove chimeras using the dereplicate method. In this method, we are again using the assumption that the highest abundance sequences are most trustworthy. Chimeras are sequences that did not extend during PCR and then served as templates for other PCR products, forming sequences that are partially from one PCR product and partially from another. This program looks for chimeras by comparing each sequences to the next highest abundance sequences to determine if a sequence is a chimera of the more abundance sequences.  

We will use the `chimera.vsearch` function to identify chimeras.

We will then remove the identified chimeras with `remove.seqs`.

Finally, we will run a new summary.

This step requires an executable program called "vsearch". This is now available as a module on Andromeda. Since I am working in Andromeda, I will load the module.

Write and run a script to do these steps.  

```
nano chimera.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/mothur
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load Mothur/1.46.1-foss-2020b

module load VSEARCH/2.18.0-GCC-10.2.0

mothur

mothur "#chimera.vsearch(fasta=ofav.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=ofav.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=T)"

mothur "#remove.seqs(fasta=ofav.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=ofav.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)"

mothur "#summary.seqs(fasta=ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=ofav.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)"

mothur "#count.groups(count=ofav.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)"

```
```
sbatch chimera.sh
```

This script outputs the following files:

```
Output File Names:
ofav.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table
ofav.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.chimeras
ofav.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos
ofav.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count.summary
```
The new summary looks like this:
```

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	450     241     0	3	1
2.5%-tile:	1	451     252     0	4	1442
25%-tile:	1	451     253     0	4	14419
Median:         1	451     253     0	4	28838
75%-tile:	1	451     253     0	5	43257
97.5%-tile:     1	451     254     0	7	56234
Maximum:        2	451     258     0	8	57675
Mean:   1	450     252     0	4
# of unique seqs:	16375
total # of seqs:        57675
```
The program identified and removed ~1.5% chimeras.

We can also look at information about sampling depth in mothur.1705444376.logfile (summary stats generated in R).

```      
Min.   :   1.0  
1st Qu.:  19.0  
Median :  63.0  
Mean   : 293.8  
3rd Qu.: 238.5  
Max.   :9245.0  
```

Sampling depth is extremely low; only ten samples have >1000 seqs. I'm going to proceed anyway, mostly just for the practice.

### <a name="Classify"></a> **9. Classifying sequences**  

Now our sequences are clean and ready for classification!  

We will use the training set downloaded above from the silva database through the [Mothur wiki](https://mothur.org/wiki/classify.seqs/).

#### Classify sequences  

We will use the `classify.seqs` command.

The training and taxonomy files were from downloads. The mothur [classify.seq wiki page](https://mothur.org/wiki/classify.seqs/) has the reference and training sets that you can download and put into the directory to use for analyzing data. These sets include adding chlorophyll, mitchondria, etc to identify and remove these.

*Outside andromeda, download these files and then transfer to Andromeda if you need to.*     

We have already downloaded the silva database in previous step from the Mothur wiki page. Download the mothur-formatted version of training set - this is version 9.

```
wget https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip

#unzip and rename if you need to
```

The output file from `classify.seqs()` ending in `.taxonomy` has the name of sequence and the classification with % confidence in parentheses for each level. It will end at the level that is has confidence.  

The tax.summary file has the taxonimc level, the name of the taxonomic group, and the number of sequences in that group for each sample.  

We will also remove sequences that are classified to Chloroplast, Mitochondria, Unknown (not bacteria, archaea, or eukaryotes), Archaea, and Eukaryotes.

Now write and run a `classify` script:
```
nano classify.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/mothur
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#classify.seqs(fasta=ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=ofav.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax)"

mothur "#remove.lineage(fasta=ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=ofav.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Cyanobacteria_Chloroplast-Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"

mothur "#summary.seqs(fasta=ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=ofav.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)"

```
```
sbatch classify.sh
```

Output files:

- The output file `.taxonomy` has name of sequence and the classification with % confidence in parentheses for each level. It will end at the level that is has confidence.

```
head ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy
```

- The `tax.summary` file has the taxonimc level, the name of the taxonomic group, and the number of sequences in that group for each sample.

```
head ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary
```

Updated summary table looks like:
```

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	450     241     0	3	1
2.5%-tile:	1	451     252     0	4	1140
25%-tile:	1	451     253     0	4	11400
Median:         1	451     253     0	4	22800
75%-tile:	1	451     253     0	5	34199
97.5%-tile:     1	451     254     0	7	44459
Maximum:        2	451     258     0	8	45598
Mean:   1	450     252     0	4
# of unique seqs:	13126
total # of seqs:        45598

```

### <a name="Cluster"></a> **10. Cluster for OTUs**  

*In this analysis, we will cluster to OTU level (cutoff=0.03). For ASV clustering, you can move directly to the make.shared step, skipping the dist.seqs and cluster steps because mothur pre-clustering occurs as the ASV level.*  

First, we will calculate the pairwise distances between sequences.  

We will first use the `dist.seqs` command.  

We will then cluster using the `cluster` command.

This will run a line for each iteration of clustering. This is run until the Matthews correlation coefficient (MCC) value is maximized. A high MCC = high confidence in clustering. MCC is optimized by randomly aligning sequences to OTU's and calculating the correlation coefficient. Then sequences are moved between OTU's to see if the MCC is improved. This is repeated many times until the MCC is maximized. This method is fast and RAM efficient. AKA Opticlust.

Next we will make a shared file. This `.shared` file has the label for the OTU, sample name, the number of OTU's and then the number of time each OTU appears in each sample.

This file will be the basis of what we will do to measure richness of communities compared to each other.  

We want to keep the shared file and a consensus taxonomy file.

The classify.otu command will then output a concensus cons.taxonomy file that has the taxonomic information for each OTU. Use label=ASV to specify ASV in taxonomy names if starting from the make.shared step for ASV's.

Then we can rename the files to something more useful.

Finally, view a count of the number of sequences in each sample.

Create a script to accomplish all of these steps:
```
nano cluster.sh
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
#SBATCH -D /data/pradalab/ceraker/16s/mothur
#SBATCH --error="script_error"
#SBATCH --output="output_script"

source /usr/share/Modules/init/sh

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#dist.seqs(fasta=ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta)"

mothur "#cluster(column=ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=ofav.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff=0.03)"

mothur "#make.shared(list=ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=ofav.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)"

mothur "#classify.otu(list=ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=ofav.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy)"

mothur "#rename.file(taxonomy=ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy, shared=ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)"

mothur "#count.groups(shared=ofav.opti_mcc.shared)"

```
```
sbatch cluster.sh
```

`Dist.seqs` outputs the following files:
```
ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist
```

`Cluster` outputs the following files:  
```
ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list
ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.steps
ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.sensspec
```

`Make.shared` outputs the following file:
```
ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared
```

Finally, `classify.otu` outputs the following file:  

```
ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy
ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.tax.summary
```

The `rename` function at the end will rename our files to something more useful. We now have a "taxonomy" file and a "shared" file.  

The files we now care about are:
```
list = ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list
shared = ofav.opti_mcc.shared
taxonomy = ofav.taxonomy

constaxonomy = ofav.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy
count = ofav.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table
```

*Now we have classified taxonomic data and we are ready to move onto subsampling and calculating statistics*  

***_START HERE NEXT TIME_***

### <a name="Subsample"></a> **11. Subsampling for Sequencing Depth**   

**The remaining commands can all be run in interactive mode. In Andromeda, use the following commands:**   

```
interactive
module load Mothur/1.46.1-foss-2020b
mothur
```

This will open the interactive mothur terminal.  In order to use any bash commands like `nano`, you have to first `quit` in mothur. Once you want to return to mothur, use the `mothur` command.  
























































.
