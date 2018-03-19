import sys
import os
import gzip
import subprocess

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

def PrimaryParser(FilePath, cancers):
    readFile = '%s/%s'%(FilePath,'final_consensus_12aug_passonly_snv_indel.maf.gz')
    n = UpdateProgressGetN(readFile)

    fileHandleIn = gzip.open(readFile, 'rb')

    i = 0
    currentCode = None
    outFile = None
    for line in fileHandleIn:
        line = line.decode('utf-8').split('\t')
        if line[0]=='Hugo_Symbol':
            pass
        else:
            code = line[44].split('-')[0]

            if code != currentCode:
                currentCode = code
                try:
                    outFile.close()
                except AttributeError:
                    pass
                if os.path.isfile(FilePath+"/Cancers/" + code + '/' + code + '-.snvs.indels.maf'):
                    outFile = open(FilePath+"/Cancers/" + code + '/' + code + '-.snvs.indels.maf', 'a')
                else:
                    outFile = open(FilePath+"/Cancers/" + code + '/' + code + '-.snvs.indels.maf', 'w')

            outFile.write('\t'.join(line))

            UpdateProgress(i, n, "Parsing into cancer specific .maf file.")
            i+=1

    try:
        # Closes the file for the last line
        outFile.close()
    except:
        pass

    fileHandleIn.close()

    for cancer in cancers:
        os.system('gzip %s'%(FilePath+"/Cancers/" + cancer + '/' + cancer + '-.snvs.indels.maf'))

    print("Complete")

def main():
    FilePath = os.path.dirname(os.path.abspath(__file__))
    with open('%s/%s'%(FilePath,'CancerTypes.txt'), 'r') as cancerFile:
        cancers = [line.replace('\n', '') for line in cancerFile.readlines()]

    if os.path.exists(FilePath+"/Cancers"):
        pass
    else:
        os.mkdir(FilePath+"/Cancers") # Throws error if it already exists.
        for type in cancers:
            os.mkdir(FilePath+"/Cancers/" + type)

    PrimaryParser(FilePath, cancers)

if __name__=='__main__':
    main()