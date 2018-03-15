import sys
import os
try:
    import ConfigParser as configparser # for python 2
except:
    import configparser # for python 3

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
    return dict1

def ENCODEtoGet(dataFiles):
    with open("%s/ENCODE_BedFileInfo.txt"%(dataFiles['encodedata']), 'r') as ENCODEfile:
        lines = [line.replace('\n','') for line in ENCODEfile.readlines()]

    with open(dataFiles['databedsfile'], 'r') as inFile:
        filePaths = [line.replace('\n','') for line in inFile.readlines()]
        myList = {line.split('\t')[0]:line.split('\t')[1] for line in filePaths}

    filesToUse = []
    for line in lines:
        if line.split('\t')[4]=='1':
            try:
                filesToUse.append("%s\t%s\t%s"%(line.split('\t')[0],myList[line.split('\t')[0]], line.split('\t')[6]))
            except KeyError:
                pass
    return(filesToUse)

def RoadmapToGet(dataFiles):
    with open(dataFiles['databedsfile'], 'r') as inFile:
        lines = [line.replace('\n','') for line in inFile.readlines()]

    dataToUse = []
    for line in lines:
        if "RoadmapGenomics" in line:
            dataToUse.append(line)

    return(dataToUse)

def MoveBeds(encodeIDs, roadmapIDs, snpEFF):
    # Prepare snpEFF to recieve bed files for database building
    os.system("mkdir %s/data/GRCh37.75/regulation.bed"%(snpEFF['snpeff']))

    for line in encodeIDs:
        line = line.split('\t')
        targetDir = "%s/data/GRCh37.75/regulation.bed/regulation.%s.%s.bed.gz" % (snpEFF['snpeff'], line[2].replace(" ","_").replace(',',''), line[0])
        os.system("cp %s %s"%(line[1], targetDir))
        os.system("gzip -d %s"%(targetDir))

    for line in roadmapIDs:
        line = line.split('\t')
        nam = line[1].split('/')[len(line[1].split('/'))-1].split('.',1)[0].split('-DNase')[0]
        targetDir = "%s/data/GRCh37.75/regulation.bed/regulation.%s.%s.bed.gz" % (snpEFF['snpeff'], line[0].replace(".","_"), nam)
        os.system("cp %s %s"%(line[1], targetDir))
        os.system("gzip -d %s"%(targetDir))

    # Cleanup any failed...
    os.system("rm %s/data/GRCh37.75/regulation.bed/*.bed.gz"%(snpEFF['snpeff']))


def main():
    localpath = os.path.abspath(__file__).rstrip('Prepare_snpEFF_information.py')  # path to scripts working directory
    Config = configparser.ConfigParser()
    Config.read(localpath + "usr_paths.ini")
    snpEFF = ConfigSectionMap(Config.sections()[0], Config)  # get snpEFF path
    dataFiles = ConfigSectionMap(Config.sections()[1], Config)  # get regulatory data

    encodeIDs = ENCODEtoGet(dataFiles)
    roadmapIDs = RoadmapToGet(dataFiles)

    MoveBeds(encodeIDs, roadmapIDs, snpEFF)

    print("WARNING: Please mind any errors that were created during the processing of files.")
    print("INFO: Complete.")


if __name__=="__main__":
    main()