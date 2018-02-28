# Summary

Processing pipeline for the PCAWG dataset to be used in conjunction with the CNN arm of Rotation_1 of my DPhil.
This process is very data heavy and requires large amounts of memory. Partitioning of data is necessary for speed and memory utilization.
The process is as follows.

###### Please note that no data is to be stored within this repository. All data storage complies with proper patient data protection.

## Preparing PCAWG Data for Analysis
```bash
cd ./PCAWGData # Not on github
bash ./SplitPCAWGData.sh
```
1. This does the following:
   - Reads from the concatenated, compressed .maf.gz file and parses the mutation data into a 
   subdirectory structure as follows:
     - Cancers
       - Type1
       - Type2
       - TypeN
       
## Parsing data into individual patient maf files for each Cancer type and converts to vcf file.
```bash
# From the github repository directory
python ./DataGrooming/ExtractData.py
```
