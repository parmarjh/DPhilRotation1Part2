'''
Created by Ryan Schenck
9 March 2018
'''
import sys
import os
import gzip
import glob
import io
from optparse import OptionParser
import datetime
import time
import subprocess
from functools import wraps
import pandas as pd


def OptionParsing():
    usage = 'usage: %prog [options] -f <*.h5>'
    parser = OptionParser(usage)
    parser.add_option('-c', '--cnvDir', dest="cnv_dir", default=None, help="Directory containing the Battenburg CNV calls.")
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

class PCAWGData:
    def __init__(self, FilePath, Options, CancerType):
        self.CancerType = CancerType
        self.patients = []
        self.tumors = []
        self.cnvFile = []
        self.cnvFileExists = []
        self.vcfFiles = self.GetVCFFiles(FilePath, Options)
        self.ProcessCNVandVCF(FilePath)

    def GetVCFFiles(self, FilePath, Options):
        vcfFiles = glob.glob(FilePath.replace("DataGrooming","PCAWGData/Cancers/%s/*.sorted.vcf.gz"%(self.CancerType)))
        for vcf in vcfFiles:
            self.patients.append(vcf.split('/')[len(vcf.split('/'))-1].split('.',1)[0])
            self.tumors.append(vcf.split('/')[len(vcf.split('/'))-1].split('.')[1])
            cnvFileName = Options.cnv_dir + vcf.split('/')[len(vcf.split('/'))-1].split('.')[1] + '.tar.gz'
            if os.path.isfile(cnvFileName):
                self.cnvFileExists.append(True)
                self.cnvFile.append(cnvFileName)
            else:
                self.cnvFileExists.append(False)
                self.cnvFile.append(None)
        return(vcfFiles)

    def ProcessCNVandVCF(self, FilePath):
        '''
        Process vcf file and initiate decompression of tar.gz and orchestrate extraction of CNV information.
        Will decompress, extract, and delete the extracted information.

        :param FilePath: Working directory
        :return: CNV information at vcf specific location.
        '''
        for k, vcf in enumerate(self.vcfFiles):
            cnvList = {}
            if self.cnvFileExists[k]:
                cnvLines = self.extractCopyNumberDir(self.cnvFile[k], FilePath)
                del cnvLines[0]

                for mut in cnvLines:
                    mutPos = mut.split('\t')[0]+':'+mut.split('\t')[1]+"-"+mut.split('\t')[2]
                    info = mut.split('\t')[3:]
                    cnvList.update({mutPos:info})
            else:
                # TODO Figure out what to append based on what gets added...
                pass
            self.AppendCopyNumber(FilePath, self.cnvFileExists[k], cnvList, vcf)

            sys.exit()

    def extractCopyNumberDir(self, cnvdir, FilePath):
        '''
        Extracts lines into memory of the appropriate file by creating a directory locally. This file is subsequently removed
        after reading in the lines of the appropriate file.
        :param cnvdir: The directory containing the CNV information this is a .tar.gz file.
        :param FilePath: The working path of the executable.
        :return: The lines from the CNV info file.
        '''
        os.system("tar -xvzf %s"%cnvdir)
        dirName = FilePath + '/' + cnvdir.split('/')[len(cnvdir.split('/'))-1].split('.')[0] + '/'
        fileWithCNVInfo = cnvdir.split('/')[len(cnvdir.split('/'))-1].replace('.tar.gz','_allDirichletProcessInfo.txt')
        file2Extract = dirName  + fileWithCNVInfo
        with open(file2Extract, 'r') as inFile:
            lines = [line.rstrip('\n') for line in inFile.readlines()]
        os.system("rm -r %s"%(dirName))
        return(lines)

    def AppendCopyNumber(self, FilePath, cnvExists, cnvList, vcfFile):
        if cnvExists:
            inFile = gzip.open(vcfFile, 'rb')
            outFile = gzip.open(vcfFile.replace('.vcf.gz','.cnv.vcf.gz'), 'wb')

            present = 0
            absent = 0
            for line in inFile:
                line = line.decode('utf-8')
                if line.startswith('#'):
                    if line.startswith('##INFO=<ID=DCC_Project_Code'):
                        outFile.write(line.encode('utf-8'))
                        outFile.write('##INFO=<ID=Subclonal.CN,Number=1,Type=Integer,Description=\"CNV Subclonal Copy Number.\">\n'.encode('utf-8'))
                        outFile.write('##INFO=<ID=Subclonal.CN,Number=1,Type=Integer,Description=\"CNV Subclonal Copy Number.\">\n'.encode('utf-8'))
                        outFile.write('##INFO=<ID=nMaj1,Number=1,Type=Integer,Description=\"CNV Major Allele Copy Number clone 1.\">\n'.encode('utf-8'))
                        outFile.write('##INFO=<ID=nMin1,Number=1,Type=Integer,Description=\"CNV Minor Allele Copy Number clone 1.\">\n'.encode('utf-8'))
                        outFile.write('##INFO=<ID=frac1,Number=1,Type=Float,Description=\"CNV Clone 1 fraction.\">\n'.encode('utf-8'))
                        outFile.write('##INFO=<ID=nMaj2,Number=1,Type=Integer,Description=\"CNV Major Allele Copy Number clone 2.\">\n'.encode('utf-8'))
                        outFile.write('##INFO=<ID=nMin2,Number=1,Type=Integer,Description=\"CNV Major Allele Copy Number clone 2.\">\n'.encode('utf-8'))
                        outFile.write('##INFO=<ID=frac2,Number=1,Type=Float,Description=\"CNV Clone 2 fraction.\">\n'.encode('utf-8'))
                    else:
                        outFile.write(line.encode('utf-8'))
                else:
                    line = line.rstrip('\n')
                    chr = line.split('\t')[0]
                    start = line.split('\t')[1]
                    # For SNVs
                    if len(line.split('\t')[3])==1 and len(line.split('\t')[4])==1:
                        end = int(line.split('\t')[1]) + 1
                    elif len(line.split('\t')[3])>1: # DELs
                        end = int(line.split('\t')[1]) + 1
                    elif len(line.split('\t')[4])>1: # INSs
                        end = int(line.split('\t')[1]) + len(line.split('\t')[4])

                    try:
                        region = cnvList['%s:%s-%s'%(chr,start,end)]
                        present+=1
                    except KeyError:
                        absent+=1
            print(present)
            print(absent)
            print(len(cnvList))

            outFile.close()
            inFile.close()

@fn_timer
def PrepareCancerClasses(Options, FilePath):
    with open(FilePath.rstrip("DataGrooming")+"PCAWGData/CancerTypes.txt", 'r') as inFile:
        cancerTypes = [line.rstrip('\n').replace("-","") for line in inFile.readlines()]

    allData = {}
    count = 0
    for cancer in cancerTypes:
        print("INFO: Processing %s"%(cancer))
        allData.update({cancer:PCAWGData(FilePath, Options, cancer)})
        count+=1

        if count == 1:
            sys.exit()

    return(allData)

def main():
    FilePath = os.path.dirname(os.path.abspath(__file__))
    now = datetime.datetime.now()
    (Options, Parser) = OptionParsing()

    PrepareCancerClasses(Options, FilePath)

if __name__=="__main__":
    main()