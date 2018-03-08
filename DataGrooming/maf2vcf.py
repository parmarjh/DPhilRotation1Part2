'''
Custom Script to convert an individual patient MAF to a vcf version 4.2 file using python 3.6.
Requries Samtools!!!
Created by Ryan Schenck
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
    parser.add_option('-s', '--spotCheckMaf', dest='spotcheck', default=False, action='store_true', help="Use this flag if ")
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

def SamtoolsFaidx(refGenome, genomicPos, ref):
    proc = subprocess.Popen(['samtools','faidx',refGenome, genomicPos], stdout=subprocess.PIPE)
    proc.wait()
    outInfo = proc.stdout.readlines()
    refSeq = ''.join([line.decode('utf-8').rstrip('\n') for line in outInfo[1:]])
    if refSeq == ref:
        return(True)
    else:
        print('ERROR: May not be proper reference genome')
        print('ERROR: Improper reference. Found %s at %s. Reference genome shows %s' % (ref, genomicPos, refSeq))
        sys.exit()

def SpotCheckProperReference(mafFile, Options, fileLength):
    '''
    Randomly samples the file to ensure proper reference file is used. Random sampling is employed to ensure proper
    reference is used. Will spot check 10% of a file of more than 200 variants.

    :param mafFile: Input mafFile object (opened)
    :param Options: Parser Options
    :param fileLength: Length of the file being read
    :return: None
    '''
    print("INFO: Verifying maf file.")
    if fileLength > 200:
        n=0.1
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
                if len(line) == 46:
                    del line[5]  # Ensures that the proper format is met by forcefully deleting the strange insterted \t
                genomicPos = line[1] + ":" + line[2] + "-" + line[3]
                ref = line[7]
                mutType = line[5]
                if mutType != "INS":
                    toContinue = SamtoolsFaidx(Options.refGenome, genomicPos, ref)
                if count == len(a):
                    print('')
                    return(toContinue)
        elif i == 0 and line.startswith('Hugo_Symbol	Chromosome	Start_position') == False:
            print("ERROR: No header found in maf file.")
        i+=1
    print('')
    return(toContinue)

def CreateHeader(ioObject, Options, tumorID, normalID):
    now = datetime.datetime.now()
    ioObject.write("##fileformat=VCFv4.2\n")
    ioObject.write("##fileDate=%s\n"%(now.date()))
    ioObject.write("##source=maf2vcf.py\n")
    ioObject.write("##reference=%s\n"%(Options.refGenome))
    ioObject.write("##sampleColumns=Normal.Tumor\n")
    ioObject.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    ioObject.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths of REF and ALT(s) in the order listed\">\n")
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
            i+=1
    tumorID = toPullIDs[12]
    normalID = toPullIDs[13]

    count = 0
    i = 0
    with open(Options.maf, 'r') as inFile:
        with open(Options.outDir + Options.maf.split('/')[len(Options.maf.split('/'))-1].replace('.maf','.vcf'), 'w') as outFile:
            CreateHeader(outFile, Options, tumorID, normalID)
            for line in inFile:
                # UpdateProgress(i, n, "Processing Maf File")
                if line.startswith('Hugo_Symbol	Chromosome	Start_position'):
                    count+=1
                    i += 1
                else:
                    i += 1
                    line = line.rstrip('\n').split('\t')
                    if len(line)==46:
                        del line[11] # Ensures that the proper format is met by forcefully deleting the strange insterted \t
                    chrom = line[1]
                    pos = line[2]
                    id = line[10]
                    if id == '':
                        id = '.'
                    elif id.startswith("rs") == False:
                        print(line)
                        sys.exit("Problem in id column")
                    strand = line[4]
                    mutType = line[5]
                    if mutType == '':
                        mutType='.'
                    variantType = line[6]
                    if variantType == '':
                        variantType = '.'
                    ref = line[7]
                    tAllele1 = line[8] # Normal Allele
                    tAllele2 = line[9] # Alt Allele
                    QUAL = line[41]
                    if QUAL=='None' or QUAL=='NA' or QUAL=='':
                        QUAL='.'
                    if ref==tAllele1:
                        altAllele=tAllele1
                        refAllele=tAllele2
                    else:
                        altAllele=tAllele2
                        refAllele=tAllele1

                    if strand=="+" or strand=="-":
                        GT="0|1"
                    else:
                        GT="0/1"

                    if variantType != 'INS' and variantType != 'DEL':
                        ref_reads = line[38]
                        alt_reads = line[37]
                        if ref_reads=='NA' or alt_reads=='NA' and line[27]=='1':
                            if ref_reads =='NA':
                                ref_reads='.'
                                total_reads=alt_reads
                            else:
                                alt_reads='.'
                                total_reads=ref_reads
                            GT="1|1"
                            vaf = line[27]
                            sampleField = ':'.join([GT, ','.join([ref_reads, alt_reads]), total_reads, vaf])
                        elif ref_reads=='NA' or alt_reads=='NA' and line[27]=='NA':
                            print("WARNING: %s" % '\t'.join(line))
                        else:
                            total_reads = str(int(ref_reads) + int(alt_reads))
                            vaf = repr(round(int(alt_reads)/float(total_reads),4))

                            sampleField = ':'.join([GT, ','.join([ref_reads, alt_reads]), total_reads, vaf])
                    else:
                        if variantType=="INS":
                            print(line)
                        pass
                    if (ref != tAllele1 and ref != tAllele2) or strand != '+':
                        print('\t'.join(line))
                        sys.exit("Differences in Alleles or malformed strand information.")
                    INFO = "Hugo_Symbol=" + line[0] + ";ref_context=" + line[14].upper() + ";Genome_Change=" + line[13] + ";Variant_Type=" + variantType + ";Variant_Classification=" + mutType
                    # if variantType != 'INS' and variantType != 'DEL':
                    lineOut = [chrom, pos, id, refAllele, altAllele, QUAL, '.', INFO, "GT:AD:DP:VF", "0|0:.:.:.",]
                    # print(lineOut)
    print('')


def main():
    print("INFO: Processing MAF file.")
    FilePath = os.path.dirname(os.path.abspath(__file__))
    (Options, Parser) = OptionParsing()

    ProcessFile(Options)

if __name__=="__main__":
    main()