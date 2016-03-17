## Arabidopsis pdx3 mutants

> Arabidopsis mutants in the PMP/PNP oxidase PDX3 show abberant growth and development.RNA sequencing reveals strong induction of stress-related genes in pdx3, particularly those associated with biotic stress. Overall design: Whole rosettes (21- days old) of three biological replicates of wild-type and pdx3-3 (SALK_054167C) and pdx3-4 (GK-260E03) mutants, all ecotype Col-0, were sequenced and a differential gene expression analysis was performed.

### External resources describing the data

- GEO: http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-77428/
- ENA: http://www.ebi.ac.uk/ena/data/view/SRP069130

The metadata from both data repositories is also archived in the `./metadata` ciredtory.

### Basic information

- **species: _Arabidopsis thaliana_**
- **protocol: single-end, stranded**
- **conditions: 3**
- **replicates: 3**
- **total reads: 191,035,737**

### Workshop steps

1. Look at the experiment methods and metadata
  - discuss the protocol choices (PE/SE, reps, depth etc.)
  - discuss the metadata quality
2. Get the reference transcriptome and index it with Salmon
  - (optional) wget the reference
  - run salmon index on the reference
3. Quantify a single rep and condition with Salmon
  - construct the salmon command for one sample
  - (optional) create a script that would run the command for all samples
  - run Salmon for however many samples are being analysed
4. Inspect the Salmon outputs (pre-generated ones are provided if needed)
  - look at the output files at the command line
  - load the data into R
  - (optional) write an R function to load a single sample
  - loop over all samples to load all the data
  - look at the loaded data frames
  - discuss TPM vs expected counts, effective length
5. Run some QC on the data
  - look at expression distributions
  - look at correlation between samples and conditions
  - normalisation
6. Perform differential expression analysis

### Setup

1. Install the required software on each student machine
2. Decide how much of the raw data students will analyse themselves
3. Run the setup script to download the reads, reference and expression data


### Resources needed

- Salmon for one sample: ~2 minutes @ 8 cores, 1 minute @ 20 cores
- Salmon for all samples: ~11 minutes @ 20 cores (digitalocean)
