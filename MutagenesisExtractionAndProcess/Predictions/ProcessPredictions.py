import sys
import numpy as np
import glob
import gzip

def CreateDifference():
    mt = '/Users/schencro/Desktop/Oxford/Rotation_1/PCAWGArm/MutagenesisExtractionAndProcess/Predictions/MUTPreds.txt'
    wt = '/Users/schencro/Desktop/Oxford/Rotation_1/PCAWGArm/MutagenesisExtractionAndProcess/Predictions/WTPreds.Match.txt'

    sites = []
    mtpreds = []
    with open(mt, 'r') as inMt:
        for line in inMt:
            line = line.split('\t')
            sites.append(line[0])
            mtpreds.append(float(line[1]))

    wtpreds = []
    with open(wt, 'r') as inWT:
        for line in inWT:
            line = line.split('\t')
            wtpreds.append(float(line[1]))

    mtpreds = np.asarray(mtpreds)
    wtpreds = np.asarray(wtpreds)

    diffs = wtpreds - mtpreds
    diffs = list(diffs)

    with open(
            '/Users/schencro/Desktop/Oxford/Rotation_1/PCAWGArm/MutagenesisExtractionAndProcess/Predictions/Diffs.txt',
            'w') as outFile:
        # Pos WT MUT Diff
        for i in range(0, len(diffs)):
            outFile.write('\t'.join([sites[i], str(wtpreds[i]), str(mtpreds[i]), str(diffs[i])]) + '\n')

def __extend(start, end, ext):
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

def PullDiffs():
    with open('/Users/schencro/Desktop/Oxford/Rotation_1/PCAWGArm/MutagenesisExtractionAndProcess/Predictions/DiffsHighLow', 'r') as inFile:
        lines = inFile.readlines()
        del lines[0]

    lines = [line.replace('"','').replace('\n','').split(',') for line in lines]

    diffDict = {}
    for line in lines:
        if line[2]=='WT':
            wtline = line
            for getMutLine in lines:
                if line[0] == getMutLine[0] and getMutLine[2] == 'SNV':
                    mtline = getMutLine
                    diffDict.update({wtline[0]:[wtline[3],mtline[3],wtline[1]]})

    return(diffDict)

def GetPatientMuts(diffDict):
    files = glob.glob('/Users/schencro/Desktop/Oxford/Rotation_1/PCAWGArm/PCAWGData/Cancers/*.vcf.gz')

    hits = []
    for vcffile in files:
        pat = vcffile.split('/')[len(vcffile.split('/'))-1].split('.')[0]
        with gzip.open(vcffile, 'rb') as inFile:
            for line in inFile:
                line = line.decode('utf-8')
                if line.startswith('#')==False:
                    line = line.split('\t')
                    start, end = __extend(int(line[1])-1,int(line[1]), 600)
                    posMut = line[0]+':'+ str(start) + '-' + str(end) + '(' + line[3] + '/' + line[4] + '.Mut)'
                    infoDict = {}
                    for item in line[7].split(';'):
                        if '=' in item:
                            infoDict.update({item.split('=')[0]: item.split('=')[1]})

                    try:
                        values = diffDict[posMut]
                        hits.append([''.join(posMut), '\t'.join(values), pat, infoDict['VAF'], infoDict['Variant_Classification']])
                    except KeyError:
                        pass
    hits = ['\t'.join(item) for item in hits]
    with open('/Users/schencro/Desktop/Oxford/Rotation_1/PCAWGArm/MutagenesisExtractionAndProcess/Predictions/COAD.matches.txt','w') as outFile:
        outFile.write('\n'.join(hits))

def main():
    # CreateDifference() # Already processed!
    diffDict = PullDiffs()
    GetPatientMuts(diffDict)



if __name__=="__main__":
    main()