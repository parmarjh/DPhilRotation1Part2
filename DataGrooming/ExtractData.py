import sys
import os
import gzip
import io
from optparse import OptionParser
import datetime
import time
import subprocess
from functools import wraps

def OptionParsing():
    usage = 'usage: %prog [options] -f <*.h5>'
    parser = OptionParser(usage)
    parser.add_option('-i', '--inputFile', dest='inputMaf', default=None, help="Raw maf file.")
    parser.add_option('-e', '--releasenotes', dest='releaseNotes', default=None, help="Release Data corresponding to MAF file.")
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

@fn_timer
def ReadFile(Options):
    n = UpdateProgressGetN(Options.inputMaf)

    inFile = gzip.open(Options.inputMaf, 'r')
    i = 0
    genes = []
    site = []
    for line in inFile:
        if line.decode("UTF-8").startswith("Hugo_Symbol\tChromosome\tStart_position") == False:
            genes.append(line.decode("UTF-8").split('\t')[0])
            site.append(line.decode("UTF-8").split('\t')[5])

            i+=1
        else:
            print(line.decode("UTF-8").split('\t'))

        UpdateProgress(i, n, "Reading maf file...")
    inFile.close()

    print(list(set(genes)))
    print(len(genes))
    print(list(set(site)))
    print(len(site))

class PCAWGData:

    def __init__(self, CancerType, header):
        self.CancerType = CancerType
        self.header = header
        self.AllPatientDataUnparsed = []
        self.PatientDataParsed = [] # Each list entry is a dict for a patient of cancer type X

    def ParseAllPatsIntoInd(self):
        for i in self.AllPatientDataUnparsed:
            self.PatientDataParsed.append(dict(zip(self.header, i)))
        self.AllPatientDataUnparsed = None

def PatientInfo(Options):
    n = UpdateProgressGetN(Options.releaseNotes)

    with open(Options.releaseNotes, 'r') as inFile:
        lines = [line.rstrip('\n').split('\t') for line in inFile.readlines()]
    head = lines[0]
    del lines[0]
    print(head)

    cancerTypes = {}
    for i,line in enumerate(lines):
        try:
            cancerTypes[line[2].split('-')[0]].AllPatientDataUnparsed.append(line)
        except KeyError:
            cancerTypes.update({line[2].split('-')[0]:PCAWGData(line[2].split('-')[0], head)})

    for item in cancerTypes:
        print(item + "-")
        cancerTypes[item].ParseAllPatsIntoInd()

    for item in head:
        print(item)


if __name__=="__main__":
    FilePath = os.path.dirname(os.path.abspath(__file__))
    now = datetime.datetime.now()
    (Options, Parser) = OptionParsing()
    allOutDir = FilePath.rstrip("DataGrooming") + "PCAWGData"

    PatientInfo(Options)
    # ReadFile(Options)
