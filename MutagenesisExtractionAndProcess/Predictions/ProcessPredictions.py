import sys
import numpy as np

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

def main():
    pass

if __name__=="__main__":
    main()