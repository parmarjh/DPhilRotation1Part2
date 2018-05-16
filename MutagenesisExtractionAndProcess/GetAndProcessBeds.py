'''
Created by Ryan Schenck
20 March 2018

This script is designed to work on the new dataset of final variant calls from the PCAWG dataset.
'''
import sys
import os
import glob
from optparse import OptionParser
from collections import OrderedDict
import time
import subprocess
from functools import wraps
import numpy as np
import gzip
import pickle
import subprocess

try:
    import ConfigParser as configparser # for python 2
except:
    import configparser # for python 3

def OptionParsing():
    usage = 'usage: %prog [options]'
    parser = OptionParser(usage)
    parser.add_option('-r', '--ref_genome', dest="refGenome",
                      default="/Users/schencro/Desktop/Bioinformatics_Tools/Ref_Genomes/Ensembl/GRCh37.75/GRCh37.75.fa",
                      help="Reference genome to be used.")
    parser.add_option('-i', '--input', dest="vcfs", help="Directory where the vcf files are stored. This will process all files in a directory.")
    parser.add_option('-e', '--extension', dest='ext', default=600, help="Size of the sequence for the model. Default=600.")
    (options, args) = parser.parse_args()
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

def UpdateProgressGetN(fileName):
    if fileName[len(fileName)-1]=="z":
        cmd = "gzip -cd %s | wc -l" % (fileName)
        pipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
    else:
        cmd = "wc -l %s" % (fileName)
        pipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
    return(int(pipe.read().decode("utf-8").lstrip(" ").split(" ")[0]))

def PullFiles(filepath):
    gzs = glob.glob(filepath+"*.vcf.gz")
    nongzs = glob.glob(filepath+"*.vcf")
    if len(gzs)>0:
        return(gzs)
    else:
        return(nongzs)

class VCFFile:

    def __init__(self, Options, vcfFile):
        self.vcfFile = vcfFile
        self.bedinfo = None
        self.bedtojoin = []

    def GetBedFile(self):
        cmd = ' '.join(['gzip -cd', self.vcfFile, "| sed -e 's/chr//'", "| awk '{OFS=\"\\t\"; if (!/^#/){print $1,$2-1,$2,$4\"/\"$5,\"+\",$8}}'"])
        result = subprocess.check_output(cmd, shell=True)
        result = result.decode('UTF-8').split('\n')
        self.bedinfo=result
        for v in result:
            v = v.split('\t')
            if v!=['']:
                try:
                    info = dict(item.split("=") for item in v[5].split(";") if '=' in item)
                except ValueError:
                    sys.exit("Error building info dictionary.")
                varclass= info['Variant_Classification']

                self.bedtojoin.append('\t'.join([v[0],v[1],v[2],v[3],v[4],varclass]))

class FinalBeds:

    def __init__(self, vcfClass):
        self.vcfClasses = vcfClass

    def SortAndMerge(self):
        if os.path.isfile('./allSites.merged.sorted.bed')==False:
            with open('./allSitesUnmerged.bed','w') as outIt:
                for vcf in self.vcfClasses:
                    for entry in self.vcfClasses[vcf].bedtojoin:
                        outIt.write(entry + '\n')

            os.system("awk '!seen[$0]++' ./allSitesUnmerged.bed | sort -k 1,1 -k2,2n > ./allSites.merged.sorted.bed")
            os.remove('./allSitesUnmerged.bed')

    def BuildWTFasta(self, ext, ref):
        '''
        Constructs a fasta of the proper window length for the model where the mutation falls within the middle of the sequence.
        :return: Sets self.wildtypeFasta as a fasta file list where i and i+1 are the fasta header and sequence, respectively.
        '''
        if os.path.isfile('./SitesToBuildFasta.bed')==False:
            with open('./allSites.merged.sorted.bed','r') as inFile:
                with open('./SitesToBuildFasta.bed','w') as outFile:
                    for pos in inFile:
                        p = pos.split('\t')
                        start = int(p[1])
                        end = int(p[2])
                        chrom = p[0]
                        start, end = self.__extend(start, end, ext)
                        outBedLine = '\t'.join([chrom,str(start),str(end),p[3],p[4],p[3]])
                        outFile.write(outBedLine + '\n')
        else:
            pass

        if os.path.isfile('./WTseqs.fasta')==False:
            cmd = 'bedtools getfasta -fi %s -bed ./SitesToBuildFasta.bed -s -fo /dev/stdout'%(ref)
            result = subprocess.check_output(cmd, shell=True)
            result = result.decode('UTF-8').split('\n')
            with open('./WTseqs.fasta','w') as outFile:
                for i in range(0,len(result)-1,2):
                    outFile.write(result[i]+'\n'+result[i+1]+'\n')
        else:
            pass

    @fn_timer
    def BuildMUTFasta(self):
        '''
        Constructs mutations based on the observed PCAWG Muts.

        :return:
        '''
        # n = UpdateProgressGetN('./WTseqs.fasta')
        i=0
        seq = ''
        with open('./WTseqs.fasta','r') as inFasta:
            with open('./MUTseqs.fasta','w') as outFasta:
                for line in inFasta:
                    # UpdateProgress(i, n, str(i) + "/" + str(n) + "Building MUT Fasta")
                    i += 1
                    if line[0] == '>':
                        if seq:
                            snv = header.split('(')[1].replace(')','')
                            ref = snv.split('/')[0]
                            alt = snv.split('/')[1]
                            if ref==seq[300]:
                                mutableSeq = [base for base in seq]
                                mutableSeq[300]=alt
                                seqout=''.join(mutableSeq)
                                assert seqout[300]==alt, "Mutation did not take place."
                                assert seq!=seqout, "Sequences are identical, but should be different."
                                assert len(seqout)==600, "Sequence is longer than expected."
                                outFasta.write('>'+header.split('(')[0]+'('+ref+'/'+alt+'.Mut)'+'\n')
                                outFasta.write(seqout+'\n')
                            else:
                                sys.exit("Error. Unable to determine appropriate snv.")
                        header = line[1:].rstrip()
                        seq = ''
                    else:
                        seq += line.rstrip()

    def __extend(self, start, end, ext):
        '''
        Extends a sequence
        :param start: Start from Bed
        :param end: End from Bed
        :param ext: Length of extension
        :return: A start, an end, and the location of the mutation.
        '''
        bedstart = int(max(0, start - ext / 2))
        end = int(bedstart + ext)
        mut=start
        assert (bedstart-end)==-(ext), "Length of sequence was unexpected."
        return(bedstart, end)

@fn_timer
def GatherVCFData(Options, file_list):
    n = len(file_list)
    i=1
    files = dict.fromkeys(file_list)
    for vcf in files:
        files[vcf]=VCFFile(Options, vcf)
        files[vcf].GetBedFile()
        UpdateProgress(i,n,"Reading vcf files.")
        i+=1
    return(files)

@fn_timer
def main():
    FilePath = os.path.dirname(os.path.abspath(__file__))
    localpath = os.path.abspath(__file__).rstrip('DataExtractionForMutagenesis.py')  # path to scripts working directory

    (Options, Parser) = OptionParsing()
    Config = configparser.ConfigParser()

    file_list = PullFiles(Options.vcfs)

    if os.path.isfile('./vcfClasses.p')==False:
        vcfClasses = GatherVCFData(Options, file_list)
        pickle.dump(vcfClasses,open('./vcfClasses.p','wb'))
    if os.path.isfile('./WTseqs.fasta')==False:
        vcfClasses = pickle.load(open('./vcfClasses.p','rb'))
        FinalBedClass = FinalBeds(vcfClasses)
        FinalBedClass.SortAndMerge()
    else:
        FinalBedClass = FinalBeds(None)

    FinalBedClass.BuildWTFasta(Options.ext, Options.refGenome)
    FinalBedClass.BuildMUTFasta()


if __name__=='__main__':
    main()