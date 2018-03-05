import sys
import os
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go

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

    ParseFile(FilePath, ParserInfo)

    CreatePlot(FilePath)


if __name__=="__main__":
    main()