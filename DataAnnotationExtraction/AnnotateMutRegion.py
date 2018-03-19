'''
Created by Ryan Schenck
8 March 2018
'''
import sys
import os
import glob
from optparse import OptionParser
import time
import subprocess
from functools import wraps
try:
    import ConfigParser as configparser # for python 2
except:
    import configparser # for python 3

def OptionParsing():
    usage = 'usage: %prog [options] -f <*.h5>'
    parser = OptionParser(usage)
    parser.add_option('-r', '--ref_genome', dest="refGenome",
                      default="/Users/schencro/Desktop/Bioinformatics_Tools/Ref_Genomes/Ensembl/GRCh37.75/GRCh37.75.fa",
                      help="Reference genome to be used for maf2vcf conversion.")
    parser.add_option('-l', '--onecancer', dest="oneCancer", default=False, action="store_true",
                      help="Used in conjunction with --cancer_dir to only process one cancer type directory.")
    parser.add_option('-c', '--cancer_name', dest='cancerName', default=None,
                      help="Cancer directory name to be processed. List of names can be found in CancerTypes.txt")
    (options, args) = parser.parse_args()
    return (options, parser)

def ConfigSectionMap(section, Config):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                print("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return(dict1)

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

class CancerData:
    def __init__(self, FilePath, Options, CancerType, vcfFilePath, snpEff):
        self.FilePath = FilePath
        self.cancer = CancerType
        self.vcfParentDir = vcfFilePath
        self.vcfFiles = glob.glob(self.vcfParentDir + "*.sorted.vcf.gz")
        self.AnnotateVCFs(snpEff)
        print('INFO: Annotations Completed.')

    def AnnotateVCFs(self, snpEff):
        count = 0
        for vcf in self.vcfFiles:
            outPath = vcf.split('/')[0:len(vcf.split('/'))-1]
            outputFile = '/'.join(outPath) + '/' + vcf.split('/')[len(vcf.split('/'))-1].replace('.vcf.gz', '.ann.vcf')
            UpdateProgress(count, len(self.vcfFiles), "INFO: Annotating vcf file.")
            self.AnnotateRegions(snpEff, vcf, outputFile)
            count += 1

    def AddRSids(self):
        # TODO Will add this at a later point. Not needed right now.
        pass

    def AnnotateRegions(self, snpEff, inputFile, outputFile):
        os.system('java -Xmx10G -jar %s -t -noStats GRCh37.75 %s > %s'%(snpEff['snpeff'], inputFile, outputFile))

@fn_timer
def ProcessFiles(Options, FilePath, snpEFF):
    print("INFO: Begging the process...")
    with open(FilePath.rstrip("DataAnnotationExtraction") + "PCAWGData/CancerTypes.txt", 'r') as inFile:
        cancerTypes = [line.rstrip('\n') for line in inFile.readlines()]

    allData = {}
    if Options.oneCancer:
        if Options.cancerName not in cancerTypes:
            sys.exit("ERROR: Unrecognized cancer_name argument provided.")
        print("INFO: Processing %s" % (Options.cancerName))
        pathToVCFs = "%sPCAWGData/Cancers/%s/" % (FilePath.rstrip("DataAnnotationExtraction"), Options.cancerName)
        allData.update({Options.cancerName: CancerData(FilePath, Options, Options.cancerName, pathToVCFs, snpEFF)})
    else:
        for cancer in cancerTypes:
            print("INFO: Processing %s"%(cancer))
            pathToVCFs = "%sPCAWGData/Cancers/%s/"%(FilePath.rstrip("DataAnnotationExtraction"), cancer)
            allData.update({cancer:CancerData(FilePath, Options, cancer, pathToVCFs, snpEFF)})


def main():
    FilePath = os.path.dirname(os.path.abspath(__file__))
    localpath = os.path.abspath(__file__).rstrip('AnnotateMutRegion.py')  # path to scripts working directory

    (Options, Parser) = OptionParsing()
    Config = configparser.ConfigParser()
    Config.read(localpath + "usr_paths.ini")
    snpEFF = ConfigSectionMap(Config.sections()[0], Config)  # get snpEFF path

    ProcessFiles(Options, FilePath, snpEFF)

if __name__=="__main__":
    main()