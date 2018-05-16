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
    parser.add_option('-c', '--cancer_name', dest='cancerName', default=None,
                      help="Cancer directory name to be processed. List of names can be found in CancerTypes.txt")
    parser.add_option('-f', '--build_final', dest='buildFinal', default=False, action='store_true', help="Instructions to build the final matrix for extracting sequences for CNN predictions.")
    parser.add_option('-u', '--unit_test', dest='unitTest', default=False, action='store_true', help='Use with --build_final for development on only a subset of the cancer groups.')
    parser.add_option('-s', '--no_stats', dest="noStats", default=True, action="store_false", help="Flag to disable statistics")
    parser.add_option('-z', '--clean', dest="clean", default=False, action="store_true")
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

class FinalBed:

    def __init__(self, vcfClass):
        self.vcfClasses = vcfClass

    def SortAndMerge(self):
        with open('./allSitesUnmerged.bed','w') as outIt:
            for vcf in self.vcfClasses:
                for entry in self.vcfClasses[vcf].bedtojoin:
                    outIt.write(entry + '\n')

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
    else:
        vcfClasses = pickle.load(open('./vcfClasses.p','rb'))

    # print(vcfClasses)
    # FinalBed(vcfClasses).SortAndMerge()


if __name__=='__main__':
    main()