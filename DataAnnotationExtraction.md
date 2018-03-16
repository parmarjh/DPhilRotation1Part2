# Summary
Scripts responsible for annotating cancer vcf files with mutation region and getting any additional rs ID's and pulling candidate sites for saturation mutagenesis with the CNN.

### Dependencies
1. Python (>=3.6)
3. Java (>=1.8.0)
2. Ensembl GRCh37.75 Reference Build as a single *.fa (indexed)
3. snpEFF (With GRCh37.75 downloaded)
4. DPhilRotation1 repository cloned and DataPreProcessing step completed for AllENCODEData [repository](https://github.com/rschenck/DPhilRotation1)

### Pre-requisites
1. Ensure that snpEFF is installed and GRCh37.75 is downloaded and ready for use.
```bash
# Download reference using...
java -jar snpEff.jar download -v GRCh37.75
```
2. Process the data used for DPhilRotation1 bed files to be used in snpEFF
   - This means configuring ./DataAnnotationExtraction/usr_paths.ini with the appropriate paths.
     - These paths lead to the processed data and tables from DPhilRotation1

# Part 1
```bash
# Place holder

```
1. Prepares snpEFF to annotate regulatory and other regions.
   - To save time, for now, only the first biological replicates are used in snpEFF if a replicate exists.
2. Appends up to date rs IDs across all sites.
# Part 2
```bash
# Placeholder
```
1. Extracts noncoding mutations to construct candidate sites and proper file for saturation mutagenesis.



# No Longer Doing this for now:
```bash
# Will adjust later to be a different script.
# Step 1, preparing the data for snpEFF
python ./DataAnnotationExtraction/Prepare_snpEFF_information.py

# Step 2, build custom regulatory annotation database
# Unable to complete locally due to initial heap size, Please see any notes about this step here: 
# http://snpeff.sourceforge.net/SnpEff_manual.html#databasesNc
cd /path/to/snpEFF/
java -Xmx40G -jar snpEff.jar build -v -onlyReg GRCh37.75
```