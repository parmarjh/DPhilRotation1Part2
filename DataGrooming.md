# Summary

Convert MAF files to VCF files and append CNV information (if needed).

###### Please note that no data is to be stored within this repository. All data storage complies with proper patient data protection.

## Dependencies
1. Python (>=3.6)
2. samtools (==1.4.1)
3. Ensembl GRCh37.75 Reference Build as a single *.fa (indexed)

## Preparing PCAWG Data for Analysis
```bash
cd ./PCAWGData # Not on github with data
python ./SplitPCAWGData.py
```
1. This does the following:
   - Reads from the concatenated, compressed .maf.gz file and parses the mutation data into a 
   subdirectory structure as follows:
     - Cancers
       - Type1
       - Type2
       - TypeN
       
## Parsing data into individual patient maf files for each Cancer type and convert to vcf file.
```bash
# From the github repository directory
python ./DataGrooming/ExtractData.py --ref_genome=/Path/To/Reference/Genome/Reference.fa

# Optionally this can be ran with an individual cancer sub-directory
python ./DataGrooming/ExtractData.py --ref_genome=/Path/To/Reference/Genome/Reference.fa -l --cancer_name=MELA

# If errors occur along the way or to verify the scripts outputs run --repair_struct with -l and --cancer_name
python ./DataGrooming/ExtractData.py --ref_genome=/Path/To/Reference/Genome/Reference.fa --repair_struct -l --cancer_name=MELA
```
1. This script will process each patients set of variants out of the maf files created for each cancer type above.
   - Optionally, if you prefer to only examine one cancer type, specify the -l flag and feed the Cancer Directory name into the scipt using --cancer_name
2. Then it will create a sorted, compressed vcf file.

## Alternatively running maf2vcf conversion alone.
```bash
python ./DataGrooming/maf2vcf.py --help

# Yields the following:
Usage: maf2vcf.py -i <*.maf> -o <directory> -r <ref.fa>

Options:
  -h, --help            show this help message and exit
  -i MAF, --input_maf=MAF
                        .maf file to be converted.
  -o OUTDIR, --output_dir=OUTDIR
                        Output directory for .vcf file
  -r REFGENOME, --ref_genome=REFGENOME
                        Reference genome to be used for maf2vcf conversion.
  -s, --spotCheckMaf    Use this flag to verify reference matching to maf
                        file. Default=False
  -v, --verbose         Use this flag to turn on verbose mode. Default=False
  
# To run (maf file should be decompressed at this point)
# Highly recommend using --spotCheckMaf (will take longer)
python ./DataGrooming/maf2vcf.py --spotCheckMaf --input_maf=/Path/To/Your/PCAWG/maf/file.maf --output_dir=/Your/Output/Dir/ --ref_genome=/Reference/Ref.fa
```
1. The purpose of the maf -> vcf conversion is to allow for analysis of mutational classification that may differ for various transcript isoforms.
2. The output of maf2vcf is a sorted and compressed .vcf file.
   - Two files will be created.
     - A *.ignoredSNVs.maf.gz file containing unprocessed variants (primarily due to insufficient information).
     - A *.sorted.vcf.gz file containing the processed variants.
     
     
## Append copy number information...
# Ignoring for now.
```bash
python ./DataGrooming/AppendCNV2VCF.py --options...
```
1. May come back to this if needed.