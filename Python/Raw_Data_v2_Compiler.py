# This takes Combined NS-MA-TajD Raw Data v2.json and Actual dNdS v2.renamed.json and combines them into one universal file, because I'm INSANELY dumb and didn't just do this when we first got the dNdS data...

import json
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project

pNpS_FILE = f'{PATH}/Combined NS-MA-TajD Raw Data v3.json' # The pN/pS, TajD, pi and background data for each gene in the dataset (both AG and MA)
dNdS_FILE = f'{PATH}/Actual dNdS v2.renamed.json' # The dN/dS data (along with background information - created with Actual_dNdS_Getter.py)

OUTPUT_FILE = f'{PATH}/Combined Selection Measures v3.json' # The name for the output file

with open(pNpS_FILE, 'r') as FoI: # Opens and reads in the pN/pS file
    pNpSData = json.loads(FoI.read())

with open(dNdS_FILE, 'r') as FoI: # Opens and reads in the dN/dS file
    dNdSData = json.loads(FoI.read())

combinedData = {}
for gene in pNpSData: # For each gene in the dataset
    pGeneData = pNpSData[gene]
    dGeneData = None
    if gene in dNdSData: # If these gene has a dN/dS value, grab its data
        dGeneData = dNdSData[gene]
    if dGeneData: # Create a combined data entry including the dN/dS and McKt values
        combinedGeneData = {'length':pGeneData['length'], 'numberOfHits':pGeneData['numberOfHits'], 'mutsPerBase':pGeneData['numberOfHits'] / pGeneData['length'], 
                            'chromosome':pGeneData['chromosome'], 'locus':pGeneData['locus'], 'ma':pGeneData['MAGene'], 'pRatio':pGeneData['ratio'], 
                            'tajD':pGeneData['tajD'], 'pi':pGeneData['piNorm'], 'dRatio':dGeneData['dRatio'], 'McKt':dGeneData['dRatio'] / pGeneData['ratio'] if dGeneData['dRatio'] != None else None}
    else: # There's no associated dN/dS value - dN/dS and McKt = None
        combinedGeneData = {'length':pGeneData['length'], 'numberOfHits':pGeneData['numberOfHits'], 'mutsPerBase':pGeneData['numberOfHits'] / pGeneData['length'], 
                            'chromosome':pGeneData['chromosome'], 'locus':pGeneData['locus'], 'ma':pGeneData['MAGene'], 'pRatio':pGeneData['ratio'], 
                            'tajD':pGeneData['tajD'], 'pi':pGeneData['piNorm'], 'dRatio':None, 'McKt':None}
    combinedData[gene] = combinedGeneData

print(len(pNpSData), len(dNdSData), len(combinedData))

with open(OUTPUT_FILE, 'w') as FoI: # Save the combined data
    FoI.write(json.dumps(combinedData))
