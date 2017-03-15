#!/usr/bin/sh
## Uses md5deep for recursive MD5 checksum, installed with brew

## Mothur MD5
md5deep -r mothur > mothur_md5.txt

## QIIME MD5 
md5deep -r qiime > qiime_md5.txt 

## DADA2 MD5 
md5deep -r dada2 > dada2_md5.txt

## Data MD5
md5deep -r data > data_md5.txt

## Reference MD5
md5deep -r reference > reference_md5.txt
