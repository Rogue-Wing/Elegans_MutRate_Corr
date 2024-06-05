# Gets the dN/dS values from the .json file of values, as well as the previously calculated pN/pS values, and appends them together along with background information on each gene

import json
import requests

import time
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project

dNdS_FILE = f'{PATH}/Actual dNdS.renamed.json' # dN/dS ratios
OTHER_DATA_FILE = f'{PATH}/Combined NS-MA-TajD Raw Data v2.json' # pN/pS and background gene information
NAME_FILE = f'{PATH}/Gene Name to WBID (Stranded).json' # Look-up table

NUMERAL_TO_NUM = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10, 'XI':11, 'XII':12, 'MtDNA':13}

def progressBar(currentPercent, tick): # Outputs a progress bar to the command line
    symbol = {0:'|', 1:'/', 2:'-', 3:'\\'}[tick]

    tick += 1
    if tick > 3:
        tick = 0

    return f"{'-' * (currentPercent)}{' ' * (100 - currentPercent)} [{currentPercent}%] {symbol}\r", tick

dNdSRatios = []
with open(dNdS_FILE, 'r') as f: # Opens the file of dN/dS ratios
    for line in f:
        dNdSRatios.append(json.loads(line.strip()))

print(len(dNdSRatios))

count = 0
for gene in dNdSRatios:
    if not gene['ratio']:
        count += 1


with open(NAME_FILE, 'r') as f: # Opens the file that acts as a lookup table to convert between gene names and WBIDs
    geneNames = json.loads(f.read())

#? Code for renaming the ratios from their original gene names to WBIDs. Once done, the .renamed.json file is created, and then used for subsequent runs of this code
'''
renamedRatios = []
errs = 0
for gene in dNdSRatios:
    try:
        #- gene['strand'] = geneNames[gene]['strand']
        renamedRatios.append({'geneID':geneNames[gene]['geneID'], 'ratio':dNdSRatios[gene]})
    except KeyError: #! If not in the Lookup table
        found = False
        for entry in geneNames:
            for transcriptID in geneNames[entry]['transcriptID']:
                if gene in transcriptID:
                    found = True
                    #- print(geneName, transcriptID)
                    #- gene['strand'] = geneNames[entry]['strand']
                    renamedRatios.append({'geneID':geneNames[entry]['geneID'], 'ratio':dNdSRatios[gene]})
                    break
            if found:
                break
        if not found:
            print(gene)
            errs += 1
            #- print(geneName)
print(errs)
print(len(renamedRatios))
'''

with open(OTHER_DATA_FILE, 'r') as f: # Open the pN/pS data
    otherData = json.loads(f.read())

#? Old code to calculate the number of genes in the pN/pS dataset using out-dated names
'''
bar, tick = progressBar(0, 0)
print(f'{bar}', end = '')
genesPerPercent = len(otherData) / 100

count = 0
for geneNum, gene in enumerate(otherData):
    merged = requests.get(f"http://api.wormbase.org//rest/field/gene/{gene}/merged_into").json()['merged_into']['data']
    if merged:
        count += 1
    
    bar, tick = progressBar(int((geneNum + 1) / genesPerPercent), tick)
    print(f'{bar}', end = '')
print(f'Done! {" " * 95}')
print(count)
'''

errs = 0
count = 0
dNdSGenes = {}
for gene in dNdSRatios: # For each gene with a dN/dS ratio, get it's other data, and combined all this together into a single entry
    try:
        geneData = otherData[gene['geneID']]
        if gene['geneID'] in dNdSGenes:
            count += 1
        dNdSGenes[gene['geneID']] = {'length':geneData['length'], 'numberOfHits':geneData['numberOfHits'], 'mutsPerBase':geneData['numberOfHits']/geneData['length'], 'dRatio':float(gene['ratio']) if gene['ratio'] != None else None, 'pRatio':geneData['ratio'], 'elegansLocus':geneData['locus'], 'ma':geneData['MAGene']}

    except KeyError:
        errs += 1

print(errs, count)
print(len(dNdSGenes))

with open(f'{PATH}/Actual dNdS v2.renamed.json', 'w') as f: # Write all these dN/dS-pN/pS entries to a new file
    f.write(json.dumps(dNdSGenes))











'''
with open(f'{PATH}/Gene Name to WBID (Stranded).json', 'r') as f:
    geneNames = json.loads(f.read())


ahringerGenes = []
with open(f'{PATH}/Ahringer Coding Promoter ATAC-seq Data.json', 'r') as f:
    for line in f:
        ahringerGenes.append(json.loads(line.strip()))

print(len(ahringerGenes))

renamedAhringerGenes = []
errs = 0
for gene in ahringerGenes:
    geneName = gene['geneID']
    try:
        gene['geneID'] = geneNames[geneName]['geneID']
        gene['strand'] = geneNames[geneName]['strand']
        renamedAhringerGenes.append(gene)
    except KeyError:
        found = False
        for entry in geneNames:
            for transcriptID in geneNames[entry]['transcriptID']:
                if geneName in transcriptID:
                    found = True
                    #- print(geneName, transcriptID)
                    gene['geneID'] = geneNames[entry]['geneID']
                    gene['strand'] = geneNames[entry]['strand']
                    renamedAhringerGenes.append(gene)
                    break
            if found:
                break
        if not found:
            print(geneName)
            errs += 1
            #- print(geneName)
print(errs)
print(len(renamedAhringerGenes))

with open(f'{PATH}/Ahringer Coding Promoter ATAC-seq Data.renamed.json', 'w') as f:
    for gene in renamedAhringerGenes:
        f.write(json.dumps(gene) + '\n')


requiredLines = []
with open(f'{PATH}/tissue-specific.ATAC-seq.dataset (Ahringer).txt', 'r') as f:
    for lineNum, line in enumerate(f):
        if lineNum != 0:
            data = line.strip().split('\t')
            if data[4] == 'coding_promoter':
                if ',' in data[3]:
                    for geneID in data[3].split(','):
                        requiredLines.append({'chromosome':data[0], 'start':data[1], 'stop':data[2], 'geneID':geneID, 'regulatoryClass':data[4], 'tissue':data[5]})
                else:
                    requiredLines.append({'chromosome':data[0], 'start':data[1], 'stop':data[2], 'geneID':data[3], 'regulatoryClass':data[4], 'tissue':data[5]})

with open(f'{PATH}/Ahringer Coding Promoter ATAC-seq Data.json', 'w') as f:
    for line in requiredLines:
        f.write(json.dumps(line) + '\n')
'''


        
    