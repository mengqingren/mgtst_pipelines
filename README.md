
## Directory Structure 
```
├── bin - executable binaries
├── code - commands used to run the pipelines
│   ├── dada 
│   └── mothur
├── data - sequencing data
├── dada2
├── mothur
├── qiime
└── reference - reference files 
```

Pipeline results are stored in the `dada2`, `mothur`, and `qiime` directories. 

## QIIME

The QIIME pipeline is based off of the Illumina overview tutorial (http://nbviewer.jupyter.org/github/biocore/qiime/blob/1.9.1/examples/ipynb/illumina_overview_tutorial.ipynb). 
The pipeline was run using MacQIIME 1.9.1 20150604. 

Run QIIME pipeline from within `qiime` directory.
`bash ../code/qiime_pipeline.sh` 

## Mothur

Mothur pipeline is based off of the mothur MiSeq SOP, https://www.mothur.org/wiki/MiSeq_SOP. Makefile and code based on code in the https://github.com/SchlossLab/Zackular_AbAOMDSS_mSphere_2015 repository. 

Run Mothur pipeline from the repository root directory. 
`make -f code/mothur/Makefile all`

## DADA2 

The DADA2 pipeline is based on the DADA2 tutorial http://benjjneb.github.io/dada2/tutorial.html. 

Run DADA2 pipeline from the repository root directory. 
`R CMD BATCH code/dada2_pipeline.R` 