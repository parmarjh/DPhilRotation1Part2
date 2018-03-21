'''
Created by Ryan Schenck
20 March 2018
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

try:
    import ConfigParser as configparser # for python 2
except:
    import configparser # for python 3

def OptionParsing():
    usage = 'usage: %prog [options]'
    parser = OptionParser(usage)
    parser.add_option('-r', '--ref_genome', dest="refGenome",
                      default="/Users/schencro/Desktop/Bioinformatics_Tools/Ref_Genomes/Ensembl/GRCh37.75/GRCh37.75.fa",
                      help="Reference genome to be used for maf2vcf conversion.")
    parser.add_option('-l', '--onecancer', dest="oneCancer", default=False, action="store_true",
                      help="Used in conjunction with --cancer_dir to only process one cancer type directory.")
    parser.add_option('-c', '--cancer_name', dest='cancerName', default=None,
                      help="Cancer directory name to be processed. List of names can be found in CancerTypes.txt")
    parser.add_option('-f', '--build_final', dest='buildFinal', default=None, help="Instructions to build the final matrix for extracting sequences for CNN predictions.")
    parser.add_option('-s', '--no_stats', dest="noStats", default=True, action="store_false", help="Flag to disable statistics")
    parser.add_option('-z', '--clean', dest="clean", default=False, action="store_true")
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
    def __init__(self, FilePath, Options, CancerType):
        self.FilePath = FilePath
        self.cancer = CancerType
        self.vcfData = PatientData(Options, "%sPCAWGData/Cancers/%s/annotated/" % (FilePath.rstrip("DataAnnotationExtraction"), Options.cancerName))

class PatientData:
    def __init__(self, Options, vcfPath):
        if Options.clean:
            os.system('rm -r %s'%(vcfPath+'tmpBeds'))
            sys.exit("Complete")
        else:
            self.vcfFiles = glob.glob(vcfPath + "*.sorted.ann.vcf.gz")
            self.bedFiles = self.BuildBeds(vcfPath) # Returns a list of all mutation sites across all patients.
            self.allSitesBedFile = self.GetUniqueRegions(vcfPath)
            self.vcfMuts = self.BuildMutsStruct()
            self.binaryMatrix, self.priorityMatrix, self.posMatrix, self.headerMatrix = self.GetMutMatrix(vcfPath)

    def BuildBeds(self, vcfPath):
        '''
        Creates a temporary directory for bed files and then fills the directory with patint specific bed files for mutation positions

        :param vcfPath: Path to the vcf files used to create tmpBeds directory
        :return: A list of the bed files that can be read from.
        '''
        tmpFilePath = vcfPath + 'tmpBeds'
        out = []
        try:
            os.mkdir(tmpFilePath)
        except FileExistsError:
            pass
        if os.path.isfile(tmpFilePath+'/allSites.sorted.bed') == False:
            for i,indFile in enumerate(self.vcfFiles):
                tmpBed = tmpFilePath + '/' + indFile.split('/')[len(indFile.split('/'))-1].replace('.sorted.ann.vcf.gz','.bed')
                UpdateProgress(i, len(self.vcfFiles), 'Building Beds')
                if os.path.isfile(tmpBed) == True:
                    out.append(tmpBed)
                else:
                    cmd = "gzip -cd %s | sed -e 's/chr//' | awk '{OFS=\"\\t\"; if (!/^#/){print $1,$2-1,$2,$4\"/\"$5,\"+\"}}' > %s"%(indFile.replace('.gz',''), tmpBed)
                    os.system(cmd)
                    out.append(tmpBed)
            print('')
            return(out)
        else:
            print("INFO: All sites, sorted bed file found.")
            return(None)

    def GetUniqueRegions(self, vcfPath):
        '''
        Extracts unique sites from bed files and cleans up mutations.

        :return: Filename of a sorted bed file of all unique sites.
        '''
        if os.path.isfile(vcfPath + 'tmpBeds/allSites.sorted.bed') == False:
            allPos = []
            for i, indFile in enumerate(self.bedFiles):
                with open(indFile, 'r') as inFile:
                    lines = [line.replace('\n','') for line in inFile.readlines()]
                for line in lines:
                    allPos.append(line)

            unique = list(set(allPos))
            print("INFO: Total mutated sites: %s"%(len(allPos)))
            print("INFO: Total shared mutations: %s"%(len(allPos)-len(unique)))

            genomicPosOnly = ['\t'.join(pos.split('\t',3)[0:3]) for pos in allPos]
            print("INFO: Total shared genomic sites (ignoring specific mutations): %s"%(len(allPos)-len(list(set(genomicPosOnly)))))

            allBedsOut = '/'.join(self.bedFiles[0].split('/')[0:len(self.bedFiles[0].split('/'))-1]) + '/allSites.bed'
            if os.path.isfile(allBedsOut)==False:
                with open(allBedsOut, 'w') as outFile:
                    for i, line in enumerate(genomicPosOnly):
                        outFile.write(line + "\n")

            if os.path.isfile(allBedsOut)==True:
                cmd = 'bedtools sort -chrThenSizeA -i %s > %s'%(allBedsOut, allBedsOut.replace('.bed','.sorted.bed'))
                os.system(cmd)
            else:
                sys.exit("ERROR: No bed file found.")


            if os.path.isfile(allBedsOut) == True:
                for toDel in self.bedFiles:
                    os.system('rm %s' % (toDel))
                os.system('rm %s'%(allBedsOut))

            return(allBedsOut.replace('.bed','.sorted.bed'))
        else:
            return(vcfPath + 'tmpBeds/allSites.sorted.bed')

    @fn_timer
    def BuildMutsStruct(self):
        '''
        Extracts and builds a dictionary structure for each vcfFile. { Patient : { Pos : MutInfo } }
        :param vcfPath: path to the vcf files.
        :return: A dictionary of patients with genomic positions and mutations.
        '''
        finalDict = OrderedDict()
        for i, patVCF in enumerate(self.vcfFiles):
            patID = patVCF.split('/')[len(patVCF.split('/'))-1]
            UpdateProgress(i, len(self.vcfFiles), "Processing vcf mutations")

            mutDict = OrderedDict()
            f = gzip.open(patVCF, 'rb')
            for line in f:
                line = line.decode('utf-8')
                if line.startswith('#')==False:
                    line = line.replace('\n', '')
                    genPos = line.split('\t')[0:2]
                    key = "%s\t%s\t%s"%(genPos[0], str(int(genPos[1])-1), str(genPos[1]))
                    val = line
                    mutDict.update({key:val})
            f.close()
            finalDict.update({ patID : mutDict })
        print('')

        return(finalDict)

    @fn_timer
    def GetMutMatrix(self, vcfPath):
        '''
        This is the primary workhorse. This function loops through the bed file and pulls out all mutations.
        Each of these mutations is determined to be present or absent in each file and what region.
        It yields a row identifier matrix, a header matrix, a priority matrix (referencing a position/row id),
        and a matrix.
        Within the matrix the following IDs hold:
            SNP-NonCoding = 1.
            INS-NonCoding = 2.
            DEL-NonCoding = 3.
            CodingOnly = 0.

        :return: Row ID matrix, Header ID matrix, Priority Site Matrix, Binary Matrix
        '''
        mapping = {'SNP':1., 'INS':2., 'DEL':3., 'C':0.}

        n = UpdateProgressGetN(self.allSitesBedFile)

        headerIDs = [patID.split('/')[len(patID.split('/'))-1].replace('.sorted.ann.vcf.gz','') for patID in self.vcfFiles]
        binaryMatrix = np.zeros(shape=(n,len(self.vcfFiles)), dtype='float32')
        headerMatrix = np.asarray(headerIDs)
        posMatrix = []
        priorityMatrix = []

        i = 0 # Row (Access positions using matrix[row,column] or matrix[i,j]
        with open(self.allSitesBedFile, 'r') as bedFile:
            for line in bedFile:
                UpdateProgress(i, n, "Extracting mutation information.")
                genomicPos = line.replace('\n','')

                posMatrix.append(genomicPos)

                j=0 # Column
                for pat in self.vcfMuts:
                    try:
                        info = self.vcfMuts[pat][genomicPos]
                    except KeyError:
                        info = None


                    if info is not None:
                        info=info.split('\t')
                        if len(info[3])==1 and len(info[4])==1: # Capture SNP
                            key, priority = self.GetImpact(info, 'SNP')
                            matrixVal = mapping[key]
                        elif len(info[4])>1: # Caputre DEL
                            key, priority = self.GetImpact(info, 'DEL')
                            matrixVal = mapping[key]
                        else: # Caputre INS
                            key, priority = self.GetImpact(info, 'INS')
                            matrixVal = mapping[key]
                        if priority==1:
                            priorityMatrix.append(genomicPos)

                        binaryMatrix[i,j]=matrixVal
                    j += 1
                i+=1
        print('')
        priorityMatrix = list(set(priorityMatrix))
        priorityMatrix.sort()
        priorityMatrix = np.asarray(priorityMatrix)
        posMatrix = np.asarray(posMatrix)

        # Clean up matrix, if any row is all zeros it should be removed from all appropriate lists
        binaryMatrix, posMatrix = self.CleanMatrixOfCoding(binaryMatrix, posMatrix)

        # Clean up env
        os.system('rm -r %s'%(vcfPath + 'tmpBeds'))
        self.vcfMuts = None # Dump for memory management.

        return(binaryMatrix, priorityMatrix, posMatrix, headerMatrix)

    def CleanMatrixOfCoding(self, primaryMatrix, posMatrix):
        '''
        Removes all coding rows from the matrix

        :param primaryMatrix: Matrix with values representing mutations
        :param posMatrix: Genomic position matrix
        :return: a cleaned matrix and a cleaned genomic position matrix.
        '''
        tmp = np.sum(primaryMatrix, axis=1)
        cleanMatrix = primaryMatrix[np.where(tmp != 0)[0],:]
        cleanPos = posMatrix[np.where(tmp!=0)[0]]
        return(cleanMatrix, cleanPos)

    def GetImpact(self, info, vType):
        '''
        Pulls variant information, specifically coding versus non-coding.

        :param info: VCF Info Field
        :param vType: Identified variant type.
        :return: string(SNP | INS | DEL | C)
        '''
        priority = set(['regulatory_region_variant', 'TF_binding_site_variant'])
        # Current annotations from snpEFF do not have regulatory region variants information.
            # Could use to validate???
        info = dict(item.split("=") for item in info[7].split(";"))['ANN'].split(',')
        infoPrim = list(set([item.split('|')[1] for item in info]))
        infoAlt = list(set([item.split('|')[7] for item in info]))
        if 'intergenic_region' in infoPrim and len(infoAlt)==1:
            return(vType, None)
        elif 'protein_coding' in infoAlt and len(infoAlt) == 1:
            return('C', None)
        elif priority.intersection(set(infoPrim)) or priority.intersection(set(infoAlt)):
            return(vType, 1)
        else:
            return(vType, None)

def GatherData(FilePath, Options):
    print("INFO: Begging the process...")
    with open(FilePath.rstrip("DataAnnotationExtraction") + "PCAWGData/CancerTypes.txt", 'r') as inFile:
        cancerTypes = [line.rstrip('\n') for line in inFile.readlines()]

    allData = {}
    if Options.oneCancer:
        if Options.cancerName not in cancerTypes:
            sys.exit("ERROR: Unrecognized cancer_name argument provided.")
        pickleFile = FilePath.replace('DataAnnotationExtraction', 'PCAWGData/Cancers/%s/CancerDataClass.p') % (Options.cancerName)
        if os.path.isfile(pickleFile):
            with open(pickleFile,'rb') as f:
                allData = pickle.load(f)
            print("INFO: Pickled class exists. Nothing to do.")
        else:
            print("INFO: Processing %s" % (Options.cancerName))
            allData.update({Options.cancerName: CancerData(FilePath, Options, Options.cancerName)})
            with open(pickleFile,'wb') as f:
                pickle.dump(allData, f)
    else:
        for cancer in cancerTypes:
            pickleFile = FilePath.replace('DataAnnotationExtraction', 'PCAWGData/Cancers/%s/CancerDataClass.p')%(cancer)
            if os.path.isfile(pickleFile):
                pass
            print("INFO: Processing %s" % (cancer))
            allData.update({cancer: CancerData(FilePath, Options, cancer)})

@fn_timer
def main():
    FilePath = os.path.dirname(os.path.abspath(__file__))
    localpath = os.path.abspath(__file__).rstrip('DataExtractionForMutagenesis.py')  # path to scripts working directory

    (Options, Parser) = OptionParsing()
    Config = configparser.ConfigParser()
    Config.read(localpath + "usr_paths.ini")
    snpEFF = ConfigSectionMap(Config.sections()[0], Config)  # get snpEFF path

    GatherData(FilePath, Options)


if __name__=='__main__':
    main()