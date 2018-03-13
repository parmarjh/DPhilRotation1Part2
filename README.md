# Summary

Processing pipeline for the PCAWG dataset to be used in conjunction with the CNN arm of Rotation_1 of my DPhil.
This process is very data heavy and requires large amounts of memory. Partitioning of data is necessary for speed and memory utilization.
The steps for this process are outlined below.

###### Please note that no data is to be stored within this repository. All data storage complies with proper patient data protection.

## Steps
1. [Convert PCAWG MAF files to VCF files.](DataGrooming.md)
2. [Annotate VCF files for variant regions (allow for isoforms) Part 1](DataAnnotationExtraction.md)
3. [Extract candidate sites for saturation mutagenesis. Part 2](https://github.com/rschenck/DPhilRotation1Part2/blob/master/DataAnnotationExtraction.md#part-2)
