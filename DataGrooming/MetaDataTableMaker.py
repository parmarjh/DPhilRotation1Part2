import sys
import os
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
from collections import OrderedDict

'''
IMPORTANT! This is only to be used after you run the ExtractData.py script and a grep -v "===" nohup.out for input.
'''

def ParseFile(FilePath, ParserInfo):
    with open(ParserInfo, 'r') as inFile:
        lines = [line.rstrip('\n') for line in inFile]

    data = {}
    count=0
    for line in lines:
        if " Patients: " in line:
            count+=1
            cancer = line.split('INFO: ')[1].replace(" Patients: ",":").split(':')[0]
            pats = line.split('INFO: ')[1].replace(" Patients: ",":").split(':')[1]
        elif " Tumours: " in line:
            count+=1
            tums = line.split('INFO: ')[1].replace(" Tumours: ",":").split(':')[1]
        else:
            count=0
            multiSite = line.split("patient: ")[1]
            if int(multiSite) > 0:
                cancer += "*"

            data.update({cancer:{"Patients":pats,"Tumours":tums}})

    WriteTable(FilePath, data)

def WriteTable(FilePath, data):
    with open(FilePath.rstrip("DataGrooming")+"DataTables/Patient_and_tumor_counts.txt", 'w') as outFile:
        outFile.write("Cancer\tPatient\tTumours\n")
        for item in data:
            outFile.write("%s\t%s\t%s\n"%(item, data[item]['Patients'], data[item]['Tumours']))

def end_of_loop():
     raise StopIteration

def GetMutationsPerMB(FilePath):
    genomeLen = 3137161264

    with open(FilePath.rstrip("DataGrooming")+"PCAWGData/MutationCounts.txt", 'r') as inFile:
        lines = [line.rstrip('\n') for line in inFile.readlines()]

    allTypes = list(set([line.lstrip().split(" ")[1] for line in lines if line.startswith("./Cancers/")==False]))

    codingMuts = ['De_novo_Start_InFrame','Nonsense_Mutation','Nonstop_Mutation','Silent','Splice_Site','Missense_Mutation','De_novo_Start_OutOfFrame','Start_Codon_SNP','Intron']
    indels = ['DEL', 'INS']
    metaInfo = 'Cancer\tPatient\tTotalVariants\tTotalSNVs\tTotalIndels\tTotalCodingRegion\tTotalNonCodingRegion'.split('\t')

    data = [] # Cancer\tPatient\tTotalVariants\tTotalSNVs\tTotalIndels\tTotalCodingRegion\tTotalNonCodingRegion\tInsertions\tDeletions\t...SNV type breakdown.
    for i, line in enumerate(lines):
        if line.startswith('./Cancers/'):
            mutInfo = dict.fromkeys(allTypes, 0)

            cType = line.replace('./Cancers/','').split('/')[0]
            pat = line.replace('./Cancers/','').split('/')[1].replace('.maf.gz','')

            Chunk = list(end_of_loop() if val.startswith('./Cancers/') else k+i for k, val in enumerate(lines[i+1:]))
            Chunk.append(Chunk[len(Chunk)-1]+1)
            for j in Chunk:
                if lines[j].startswith('./Cancers/')==False:
                    count, mType = lines[j].lstrip().split(' ')
                    mutInfo[mType] = count

            mutInfo.update({'Cancer':cType, 'Patient':pat})

            TotalVariants = sum([int(mutInfo[key]) for key in allTypes])
            TotalSNVs = sum([int(mutInfo[key]) for key in allTypes if key not in indels])
            TotalIndels = sum([int(mutInfo[key]) for key in allTypes if key in indels])
            TotalCodingRegion = sum([int(mutInfo[key]) for key in allTypes if key in codingMuts])
            TotalNonCodingRegion = sum([int(mutInfo[key]) for key in allTypes if key not in codingMuts])

            mutInfo.update({'TotalVariants':TotalVariants, 'TotalSNVs':TotalSNVs, 'TotalIndels':TotalIndels, 'TotalCodingRegion':TotalCodingRegion, 'TotalNonCodingRegion':TotalNonCodingRegion})

            data.append(mutInfo)

    OutLineHeader = metaInfo + allTypes
    with open(FilePath.rstrip('DataGrooming') + "DataTables/MutationBreakdown.txt", 'w') as outFile:
        outFile.write('\t'.join(OutLineHeader) + '\n')
        for tumor in data:
            outFile.write('\t'.join([str(tumor[key]) for key in OutLineHeader]) + '\n')

def CreatePlot(FilePath):
    df = pd.read_csv(FilePath.strip("DataGrooming") + "DataTables/Patient_and_tumor_counts.txt", sep="\t", index_col=False, header=0)
    sums = sum(df['Patient'])
    sumt = sum(df['Tumours'])
    df['FracPat'] = df['Patient']/float(sums)
    df['FracTum'] = df['Tumours']/float(sumt)

    dfVis = df.loc[ df['FracPat']>=0.015 ]
    dfOt = df.loc[ df['FracPat']<0.015 ]

    dfOt = pd.DataFrame({"Cancer":['Other'], "Patient":sums-sum(dfVis['Patient']), "Tumours":sumt-sum(dfVis['Tumours']), "FracPat":sum(dfOt['FracPat']), "FracTum":sum(dfOt['FracTum']) })

    dfVis = pd.concat([dfVis, dfOt])

    trace = go.Pie(labels=dfVis['Cancer'], values=dfVis['Patient'], textinfo='label+percent')

    layout = go.Layout(title='Cancer Types', width=900, height=900, showlegend=False)
    fig = go.Figure(data=[trace], layout=layout)
    py.image.save_as(fig, filename=FilePath.strip("DataGrooming") + "/DataTables/CancerTypes.svg")

def main():
    FilePath = os.path.dirname(os.path.abspath(__file__))
    ParserInfo = FilePath.rstrip("DataGrooming") + "PCAWGData/ParserInformation.txt"

    # ParseFile(FilePath, ParserInfo)

    # CreatePlot(FilePath)

    # GetMutationsPerMB(FilePath)

    # os.system("Rscript %s/VisMutationBreakdown.R %s %s" % (FilePath, FilePath.rstrip("DataGrooming") + "DataTables/", 'MutationBreakdown.txt'))

if __name__=="__main__":
    main()