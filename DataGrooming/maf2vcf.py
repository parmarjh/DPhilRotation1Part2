'''
Script to convert a MAF to a vcf4.2 file using python >=3.6.
Created by Ryan Schenck
8 March 2018
'''

import os
import sys
from optparse import OptionParser
import subprocess
from functools import wraps
import datetime
import time
import numpy as np

def OptionParsing():
    usage = 'usage: %prog -i <*.maf> -o <directory> -r <ref.fa>'
    parser = OptionParser(usage)
    parser.add_option('-i', '--input_maf', dest="maf", default=None, help=".maf file to be converted.")
    parser.add_option('-o', '--output_dir', dest="outDir", default=None, help="Output directory for .vcf file")
    parser.add_option('-r', '--ref_genome', dest="refGenome", default="/Users/schencro/Desktop/Bioinformatics_Tools/Ref_Genomes/Ensembl/GRCh37.75/GRCh37.75.fa", help="Reference genome to be used for maf2vcf conversion.")
    parser.add_option('-s', '--spotCheckMaf', dest='spotcheck', default=False, action='store_true', help="Use this flag to verify reference matching to maf file. Default=False")
    parser.add_option('-v', '--verbose', dest='verbose', default=False, action='store_true', help="Use this flag to turn on verbose mode. Default=False")
    (options, args) = parser.parse_args()
    if options.maf is None or options.outDir is None or options.refGenome is None:
        print("ERROR: Please include arguments for maf file, output directory, and reference genome (single fasta file).")
        sys.exit()
    else:
        pass
    return (options, parser)

def fn_timer(function):
    '''
    Use this as a wrapper at the top of any function you want to get run time information about.

    :param function: Function of interest.
    :return: A function to wrap around a function.
    '''
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print ("INFO: Total time running %s: %s minutes" %
               (function.__name__, str(round((t1-t0)/60.,2)))
               )
        return result
    return function_timer

def UpdateProgressGetN(fileName):
    if fileName[len(fileName)-1]=="z":
        cmd = "gzip -cd %s | wc -l" % (fileName)
        pipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
    else:
        cmd = "wc -l %s" % (fileName)
        pipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
    return(int(pipe.read().decode("utf-8").lstrip(" ").split(" ")[0]))

def UpdateProgress(i, n, DisplayText):
    '''
    Prints a progress bar where appropriate.

    :param i: Current Step
    :param n: Total number of steps.
    :param DisplayText: A string that you want to print out that is informative.
    :return: None
    '''
    sys.stdout.write('\r')
    j = (i + 1) / n
    sys.stdout.write("[%-20s] %d%%\t INFO: %s" % ('=' * int(20 * j), 100 * j, DisplayText))
    sys.stdout.flush()

def SamtoolsFaidx(refGenome, genomicPos, ref='', check=True):
    '''
    Obtain reference sequence and perform check if needed.

    :param check: Whether or not to throw error if the provided reference matches
    :param refGenome: Reference Fasta file
    :param genomicPos: Genomic Position of interest.
    :param ref: Reference sequence to compare to fetched sequence.
    :return: Fetched reference sequence.
    '''
    proc = subprocess.Popen(['samtools','faidx',refGenome, genomicPos], stdout=subprocess.PIPE)
    proc.wait()
    outInfo = proc.stdout.readlines()
    refSeq = ''.join([line.decode('utf-8').rstrip('\n') for line in outInfo[1:]])

    if check:
        if refSeq == ref:
            return(True)
        else:
            print('ERROR: May not be proper reference genome')
            print('ERROR: Improper reference. Found %s at %s. Reference genome shows %s' % (ref, genomicPos, refSeq))
            sys.exit()
        return(None)
    else:
        return(refSeq)

def SpotCheckProperReference(mafFile, Options, fileLength):
    '''
    Randomly samples the file to ensure proper reference file is used. Random sampling is employed to ensure proper
    reference is used. Will spot check 2% of a file of more than 200 variants.

    :param mafFile: Input mafFile object (opened)
    :param Options: Parser Options
    :param fileLength: Length of the file being read
    :return: None
    '''
    print("INFO: Verifying maf file.")
    if fileLength > 200:
        n=0.02
    else:
        n=1.
    a = np.arange(fileLength)
    np.random.shuffle(a)
    a = list(a[:int(fileLength*n)])
    i = 0
    count = 0
    for line in mafFile:
        if i != 0:
            checkIt = len([k for k in a if k==i])
            if checkIt==1:
                UpdateProgress(count, len(a), "INFO: Verifying maf file")
                count+=1
                line = line.rstrip('\n').split('\t')
                genomicPos = line[1] + ":" + line[2] + "-" + line[3]
                ref = line[7]
                mutType = line[5]
                variantClass = line[6]
                if variantClass != "INS" and variantClass != "TNP" and variantClass !="ONP":
                    toContinue = SamtoolsFaidx(Options.refGenome, genomicPos, ref)
                if count == len(a):
                    print('')
                    return(toContinue)
        elif i != 0 and line.startswith('Hugo_Symbol	Chromosome	Start_position') == False:
            print("")
            print("ERROR: No header found in maf file.")
        i+=1
    print('')
    return(toContinue)

def processSNP(line, chrom, pos, rsid, mutType, variantType, strand, errorFile, Options):
    ref = line[7]
    tAllele1 = line[8]  # Normal Allele
    tAllele2 = line[9]  # Alt Allele
    QUAL = line[42]

    if QUAL == 'None' or QUAL == 'NA' or QUAL == '':
        QUAL = '.'
    if ref == tAllele1:
        altAllele = tAllele1
        refAllele = tAllele2
    else:
        altAllele = tAllele2
        refAllele = tAllele1

    ref_reads = line[39]
    alt_reads = line[38]
    reportedVAF = line[28]
    # Get phasing information and determine reads for vaf==1
    if ref_reads == 'NA' or alt_reads == 'NA' and reportedVAF == '1':
        GT = "1/1" # Appears to be homozygous for alternative allele (germline unlikely since it is called w.r.t normal?)
        vaf = reportedVAF # Sets VAF equal to 1
        if ref_reads == 'NA':
            ref_reads = '.'
            total_reads = alt_reads
        else:
            alt_reads = '.'
            total_reads = ref_reads

        sampleField = ':'.join([GT, ','.join([ref_reads, alt_reads]), total_reads, vaf])

    # Tossing these very strange mutations within the MAF file.
    elif ref_reads == 'NA' or alt_reads == 'NA' and reportedVAF == 'NA':
        with open(errorFile, 'a') as errerOut:
            errerOut.write('\t'.join(line)+'\n')
        if Options.verbose:
            print("WARNING: %s" % '\t'.join(line))
        return(None)

    # Simple SNV cases
    else:
        total_reads = str(int(ref_reads) + int(alt_reads))
        vaf = repr(round(int(alt_reads) / float(total_reads), 4))

        if vaf != '1.' and strand=="+" or strand=="-":
            GT="0|1"
        else:
            GT="0/1"

        sampleField = ':'.join([GT, ','.join([ref_reads, alt_reads]), total_reads, vaf])

    # Last check for interesting but unresolved MAF line
    if (ref != tAllele1 and ref != tAllele2) or (strand != '+' and strand != '-'):
        with open(errorFile, 'a') as errerOut:
            errerOut.write('\t'.join(line)+'\n')
        if Options.verbose:
            print("WARNING: %s" % '\t'.join(line))
        return(None)

    # Create INFO field
    INFO = "MAF_Hugo_Symbol=" + line[0] + ";MAF_ref_context=" + line[15].upper() + ";MAF_Genome_Change=" + line[14] + ";MAF_Variant_Type=" + variantType + ";MAF_Variant_Classification=" + mutType +";DCC_Project_Code=" + line[44]

    # Normal variant field if anything
    if line[41]=="NA":
        normalGenotype = ".:.,.:.:."
    else:
        normalGenotype = ".:.,.:.:%s"%(line[41])

    # Final vcf line out
    lineOut = [chrom, pos, rsid, refAllele, altAllele, QUAL, '.', INFO, "GT:AD:DP:VF", normalGenotype, sampleField]

    return(lineOut)

def processDEL(line, chrom, pos, rsid, mutType, variantType, strand, errorFile, Options):
    ref = line[7]
    tAllele1 = line[8]  # Normal Allele Typically
    tAllele2 = line[9]  # Alt Allele Typically

    QUAL = line[42]
    if QUAL == 'None' or QUAL == 'NA' or QUAL == '':
        QUAL = '.'

    if ref == tAllele1:
        altAllele = tAllele1
        refAllele = tAllele2
    else:
        altAllele = tAllele2
        refAllele = tAllele1

    # Obtain the reference sequence + 1 preceding base for the DEL
    refAnchorPos = str(int(pos)-1) # Fetch the base that precedes the deletion.
    refSeq = SamtoolsFaidx(Options.refGenome, chrom + ":" + refAnchorPos + "-" + line[3], check=False)

    if refSeq[1:] != altAllele:
        print("ERROR: Deletion alternative allele does not match reference sequence. %s" % ('\t'.join(line)))
        sys.exit()

    # VCF reference is the preceding base plus the reported deletion in the MAF file.
    vcfRef = refSeq

    # VCF has base directly preceding the deletion as the alternative base and the variant pos
    vcfAlt=refSeq[0]
    vcfPos=refAnchorPos

    # Get read information
    iref_reads = line[37]
    ialt_reads = line[36]
    ref_reads = line[39]
    alt_reads = line[38]
    reportedVAF = line[28]
    i_t_vaf = line[43]
    # Get phasing information and determine reads for vaf==1
    if (ref_reads != 'NA' or iref_reads!='NA') and (alt_reads != 'NA' or ialt_reads!='NA'):
        GT="0/1"
        ref_reads = [read for read in [ref_reads, iref_reads] if read != "NA"][0]
        alt_reads = [read for read in [alt_reads, ialt_reads] if read != "NA"][0]
        total_reads = str(int(ref_reads) + int(alt_reads))
        vaf = str(int(alt_reads)/float(total_reads))
    elif i_t_vaf!="" and i_t_vaf!="NA" and ref_reads == 'NA' and iref_reads=='NA' and alt_reads == 'NA' and ialt_reads=='NA':
        vaf=i_t_vaf
        GT="./."
        ref_reads = '.'
        alt_reads = '.'
        total_reads = '.'
    elif (i_t_vaf=="" or i_t_vaf=="NA") and ref_reads == 'NA' and iref_reads=='NA' and alt_reads == 'NA' and ialt_reads=='NA':
        GT='./.'
        ref_reads='.'
        alt_reads='.'
        total_reads='.'
        vaf='.'
    else:
        sys.exit("ERROR: Problem processing DEL %s"%('\t'.join(line)))

    sampleField = ':'.join([GT, ','.join([ref_reads, alt_reads]), total_reads, vaf])

    # Create INFO field
    INFO = "MAF_Hugo_Symbol=" + line[0] + ";MAF_ref_context=" + line[15].upper() + ";MAF_Genome_Change=" + line[
        14] + ";MAF_Variant_Type=" + variantType + ";MAF_Variant_Classification=" + mutType + ";DCC_Project_Code=" + \
           line[44]

    # Normal variant field if anything
    if line[41] == "NA":
        normalGenotype = ".:.,.:.:."
    else:
        normalGenotype = ".:.,.:.:%s" % (line[41])

    lineOut = [chrom, vcfPos, rsid, vcfRef, vcfAlt, QUAL, '.', INFO, "GT:AD:DP:VF", normalGenotype, sampleField]

    return(lineOut)

def processINS(line, chrom, pos, rsid, mutType, variantType, strand, errorFile, Options):
    ref = line[7]
    tAllele1 = line[8]  # Normal Allele Typically
    tAllele2 = line[9]  # Alt Allele Typically

    QUAL = line[42]
    if QUAL == 'None' or QUAL == 'NA' or QUAL == '':
        QUAL = '.'

    if tAllele1 == '-':
        altAllele = tAllele2
    else:
        altAllele = tAllele1

    # Obtain the reference sequence + 1 preceding base for the DEL
    refAnchorPos = str(int(pos) - 1)  # Fetch the base that precedes the deletion.
    refSeq = SamtoolsFaidx(Options.refGenome, chrom + ":" + refAnchorPos + "-" + line[3], check=False)

    # VCF reference is the preceding base in the insertion in MAF
    vcfRef = refSeq[0]

    # VCF has base directly preceding the deletion as the alternative base and the variant pos
    vcfAlt = refSeq[0]+altAllele
    vcfPos = refAnchorPos

    # Get read information
    iref_reads = line[37]
    ialt_reads = line[36]
    ref_reads = line[39]
    alt_reads = line[38]
    reportedVAF = line[28]
    i_t_vaf = line[43]
    # Get phasing information and determine reads for vaf==1
    if (ref_reads != 'NA' or iref_reads != 'NA') and (alt_reads != 'NA' or ialt_reads != 'NA'):
        GT = "0/1"
        ref_reads = [read for read in [ref_reads, iref_reads] if read != "NA"][0]
        alt_reads = [read for read in [alt_reads, ialt_reads] if read != "NA"][0]
        total_reads = str(int(ref_reads) + int(alt_reads))
        vaf = str(int(alt_reads) / float(total_reads))
    elif i_t_vaf != "" and i_t_vaf != "NA" and ref_reads == 'NA' and iref_reads == 'NA' and alt_reads == 'NA' and ialt_reads == 'NA':
        vaf = i_t_vaf
        GT = "./."
        ref_reads = '.'
        alt_reads = '.'
        total_reads = '.'
    elif (
            i_t_vaf == "" or i_t_vaf == "NA") and ref_reads == 'NA' and iref_reads == 'NA' and alt_reads == 'NA' and ialt_reads == 'NA':
        GT = './.'
        ref_reads = '.'
        alt_reads = '.'
        total_reads = '.'
        vaf = '.'
    else:
        sys.exit("ERROR: Problem processing INS %s" % ('\t'.join(line)))

    sampleField = ':'.join([GT, ','.join([ref_reads, alt_reads]), total_reads, vaf])

    # Create INFO field
    INFO = "MAF_Hugo_Symbol=" + line[0] + ";MAF_ref_context=" + line[15].upper() + ";MAF_Genome_Change=" + line[
        14] + ";MAF_Variant_Type=" + variantType + ";MAF_Variant_Classification=" + mutType + ";DCC_Project_Code=" + \
           line[44]

    # Normal variant field if anything
    if line[41] == "NA":
        normalGenotype = ".:.,.:.:."
    else:
        normalGenotype = ".:.,.:.:%s" % (line[41])

    lineOut = [chrom, vcfPos, rsid, vcfRef, vcfAlt, QUAL, '.', INFO, "GT:AD:DP:VF", normalGenotype, sampleField]

    return (lineOut)

def CreateVCFLine(line, errorFile, Options):
    line = line.rstrip('\n').split('\t')

    # Genomic Position
    chrom, pos, id = line[1], line[2], line[10]

    # Get rs ID
    rsid = line[10]
    if rsid == '':
        rsid = '.'
    elif rsid.startswith("rs") == False:
        if Options.verbose:
            print("ERROR: %s"%(line))
        sys.exit("ERROR: Problem in id column")

    # Strand Information
    strand = line[4]

    # Variant Classification/Type (Type is SNP, INS, DEL, etc.)
    mutType = line[5]
    variantType = line[6]

    # Create proper vcf formatted information
    if mutType == '':
        mutType = '.'
    if variantType == '':
        variantType = '.'

    # Determine type of variant to continue processing.
    linetowrite = None
    if variantType=="SNP":
        linetowrite = processSNP(line, chrom, pos, rsid, mutType, variantType, strand, errorFile, Options)
    elif variantType=="DEL":
        linetowrite = processDEL(line, chrom, pos, rsid, mutType, variantType, strand, errorFile, Options)
    elif variantType=="INS":
        linetowrite = processINS(line, chrom, pos, rsid, mutType, variantType, strand, errorFile, Options)
    elif variantType=="TNP" or variantType=="ONP":
        with open(errorFile, 'a') as errerOut:
            errerOut.write('\t'.join(line)+'\n')
    else: # This may seem duplicitious, but I explicityly want to know as much of what I'm choosing to filter out as possible...
        if Options.verbose:
            print("WARNING: Malformed MAF entry. %s"%('\t'.join(line)))
        print('')
        with open(errorFile, 'a') as errerOut:
            errerOut.write('\t'.join(line)+'\n')
        # print(line)
        # sys.exit("ERROR: Malformed MAF entry.")

    return(linetowrite)


def CreateHeader(ioObject, Options, tumorID, normalID):
    now = datetime.datetime.now()
    ioObject.write("##fileformat=VCFv4.2\n")
    ioObject.write("##fileDate=%s\n"%(now.date()))
    ioObject.write("##source=maf2vcf.py\n")
    ioObject.write("##reference=%s\n"%(Options.refGenome))
    ioObject.write("##sampleColumns=Normal.Tumor\n")
    ioObject.write("##INFO=<ID=MAF_Hugo_Symbol,Number=1,Type=String,Description=\"HUGO Symbol in original MAF file.\">\n")
    ioObject.write("##INFO=<ID=MAF_ref_context,Number=1,Type=String,Description=\"Reference context in original MAF file.\">\n")
    ioObject.write("##INFO=<ID=MAF_Genome_Change,Number=1,Type=String,Description=\"Genome change in original MAF file.\">\n")
    ioObject.write("##INFO=<ID=MAF_Variant_Type,Number=1,Type=String,Description=\"Variant type (SNP,INS,DEL) in original MAF file.\">\n")
    ioObject.write("##INFO=<ID=MAF_Variant_Classification,Number=1,Type=String,Description=\"Variant Classification (if SNP) in original MAF file.\">\n")
    ioObject.write("##INFO=<ID=DCC_Project_Code,Number=1,Type=String,Description=\"DCC Project Code in original MAF file.\">\n")
    ioObject.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    ioObject.write("##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"Allelic depths of REF and ALT(s) in the order listed\">\n")
    ioObject.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth across this site\">\n")
    ioObject.write("##FORMAT=<ID=VF,Number=1,Type=Float,Description=\"Variant Allele Frequency.\">\n")
    ioObject.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\t%s\n"%(normalID,tumorID))

@fn_timer
def ProcessFile(Options):
    n = UpdateProgressGetN(Options.maf)

    if Options.spotcheck:
        with open(Options.maf, 'r') as inFile:
            SpotCheckProperReference(inFile, Options, n)

    with open(Options.maf,'r') as inFile:
        i = 0
        for line in inFile:
            if i == 1:
                toPullIDs = line.rstrip('\n').split('\t')
                break
            else:
                header = line
            i+=1
    tumorID = toPullIDs[12]
    normalID = toPullIDs[13]

    count = 0
    i = 0
    with open(Options.maf, 'r') as inFile:
        with open(Options.outDir + Options.maf.split('/')[len(Options.maf.split('/'))-1].replace('.maf','.vcf'), 'w') as outVCF:
            errorFile = Options.outDir + Options.maf.split('/')[len(Options.maf.split('/')) - 1].replace('.maf', '.ignoredSNVs.maf')
            with open(errorFile, 'w') as errorOut:
                errorOut.write(header)

            CreateHeader(outVCF, Options, tumorID, normalID)
            for line in inFile:
                UpdateProgress(i, n, "Processing Maf File")
                if line.startswith('Hugo_Symbol	Chromosome	Start_position'):
                    count+=1
                    i += 1
                else:
                    i += 1
                    linetoWrite = CreateVCFLine(line, errorFile, Options)

                    if linetoWrite is not None:
                        outVCF.write('\t'.join(linetoWrite)+'\n')

    print('')
    print("INFO: Sorting vcf file.")
    vcfFile = Options.outDir + Options.maf.split('/')[len(Options.maf.split('/'))-1].replace('.maf','.vcf')
    vcfFileSorted = Options.outDir + Options.maf.split('/')[len(Options.maf.split('/'))-1].replace('.head.maf','.sorted.vcf.gz')
    os.system("cat %s | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"LC_ALL=C sort -k1,1 -k2,2n\"}' | gzip > %s"%(vcfFile, vcfFileSorted))
    os.system("rm %s"%(vcfFile))
    os.system("gzip %s"%(errorFile))

def main():
    print("INFO: Processing MAF file.")
    FilePath = os.path.dirname(os.path.abspath(__file__))
    (Options, Parser) = OptionParsing()

    ProcessFile(Options)

if __name__=="__main__":
    main()