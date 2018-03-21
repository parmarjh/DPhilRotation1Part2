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
python ./DataAnnotationExtraction/AnnotateMutRegion.py
```
1. Annotates variants using snpEFF to allow for isoforms
2. Appends up to date rs IDs across all sites (CURRENTLY DISABLED).
# Part 2
```bash
python ./DataAnnotationExtraction/DataExtractionForMutagenesis.py --help
python ./DataAnnotationExtraction/DataExtractionForMutagenesis.py -l --cancer_name=MELA
```
1. Extracts and formats noncoding mutations as annotated by snpEFF.
   - Different isoforms are allowed. If any fall outside of a coding region they are considered.
2. This script yields a pickle file containing the following:
   - CancerData Class with the following attributes:
     - Cancer type
     - PatientData Class with the following attributes:
       - vcf Files
       - bed files (created within the class)
       - vcf Mutations (not saved due to the requirement of large memory)
         - If you want this, access it by running CancerData.vcfData.BuildMutsStruct()
       - Matrices of mutations for all noncoding sites, priority targets (e.g. TF binding sites), and genomic position.
3. Merge all mutation matrices across all cancer types.
```bash
python ./DataAnnotationExtraction/DataExtractionForMutagenesis.py --build_final
```
   - Note this can only be ran once all cancer types have been processed as outlined above.
   - This will merge all matrices to be used for targeted saturation mutagenesis.

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