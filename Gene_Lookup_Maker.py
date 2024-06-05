# This takes a input of an annotated genome (Annotated Ce11 Genome.gff), or a BioMart .txt and creates a JSON file of GeneName:WBGene ID pairs

import json
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project
GENOME = f'{PATH}/Annotated Ce11 Genome.gff'
BIOMART_FILE = f'{PATH}/Raw Gene ID Strand Data.txt'

def saveJSONFile(fileData, fileName): #* Saves the look-up table to a JSON file
    with open(f'{PATH}/{fileName}.json', 'w') as FoI:
        if isinstance(fileData, dict):
            FoI.write(json.dumps(fileData))
        elif isinstance(fileData, list):
            for element in fileData:
                FoI.write(json.dumps(element) + '\n')


def getNamesFromBiomart(): #* Gets the gene names and WBIDs from the Biomart file, and creates two dicts to convert between the two of them
    nameToID = {}
    IDToName = {}
    with open(BIOMART_FILE, 'r') as FoI:
        for lineNum, line in enumerate(FoI):
            if lineNum != 0:
                data = line.strip().split('\t')
                geneID = data[1]
                transcriptID = data[2]
                geneName = data[3]
                strand = data[4]

                if geneID in IDToName:
                    IDToName[geneID]['transcriptID'].append(transcriptID) # Only the transcriptID can vary - the geneName will always be the same for any WBID, same as with the strand
                else:
                    IDToName[geneID] = {'transcriptID':[transcriptID], 'name':geneName, 'strand':strand}
                
                if geneName in nameToID:
                    nameToID[geneName]['transcriptID'].append(transcriptID) # Only the transcriptID can vary - the geneName will always be the same for any WBID, same as with the strand
                else:
                    nameToID[geneName] = {'transcriptID':[transcriptID], 'geneID':geneID, 'strand':strand}

    return nameToID, IDToName             


nameToID, IDToName = getNamesFromBiomart()
saveJSONFile(nameToID, 'Gene Name to WBID (Stranded)')
saveJSONFile(IDToName, 'WBID to Gene Name (Stranded)')







                        

