# Determines the similarity between the genes of a dataset after binning by pN/pS, pi and/or dN/dS

import json
import math
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project

MA_FILE = f'{PATH}/Combined MA Data v3.json' # The file of all the MA mutants (and thus genes)
PROMOTER_FILE = f'{PATH}/Ahringer Coding Promoter ATAC-seq Data.renamed.json'

pNpS_RATIO_FILE = f'{PATH}/dNdS_Ratios v2.renamed.json' # The file of all the pN/pS ratios
PI_FILE = f'{PATH}/Combined NS-MA-TajD Raw Data v3.json' # The file that contains the pi values
dNdS_RATIO_FILE = f'{PATH}/Actual dNdS v2.renamed.json' # The file of dN/dS ratios for every gene that has one

OUTPUT_FILE = f'{PATH}/Binned MA Genes Similarity Output.json'

def getLabelledMAGenes():
    maMutants = []
    with open(MA_FILE, 'r') as f: # Load in the MA mutants
        for line in f:
            if line[0] != '#':
                maMutants.append(json.loads(line.strip()))

    with open(pNpS_RATIO_FILE, 'r') as f: # Load in the pNpS ratios
        ratios = json.loads(f.read())
    for mutant in maMutants: # Assign the pN/pS ratios to every MA mutant that has a WBID with associated pN/pS data
        try:
            mutant['pRatio'] = ratios[mutant['gene']]['ratio']
        except: # If there is no associated pN/pS value, just mark it as None
            mutant['pRatio'] = None

    with open(PI_FILE, 'r') as f: # As above, but with pi values
        data = json.loads(f.read())
    for mutant in maMutants:
        try:
            mutant['pi'] = data[mutant['gene']]['piNorm']
        except:
            mutant['pi'] = None

    with open(dNdS_RATIO_FILE, 'r') as f: # Also as above, but with dN/dS values instead
        ratios = json.loads(f.read())
    for mutant in maMutants:
        try:
            mutant['dRatio'] = ratios[mutant['gene']]['dRatio']
        except:
            mutant['dRatio'] = None

    maGenes = []
    geneNames = set()
    for mutant in maMutants:
        if mutant['gene'] not in geneNames:
            geneNames.add(mutant['gene'])
            maGenes.append({'gene':mutant['gene'], 'pRatio':mutant['pRatio'], 'pi':mutant['pi'], 'dRatio':mutant['dRatio']})

    return maGenes

def getLabelledPromoters(): #* Reads in the promoters and labels them with their appropriate pN/pS, pi and dN/dS bins, before removing any duplicates
    promoters = []
    with open(PROMOTER_FILE, 'r') as f: # Load in the promoters
        for line in f:
            data = json.loads(line.strip())
            promoters.append({'gene':data['geneID'], 'midPos':(int(data['stop']) + int(data['start'])) / 2})

    with open(pNpS_RATIO_FILE, 'r') as f: # Load in the pNpS ratios
        ratios = json.loads(f.read())
    for promoter in promoters: # Assign the pN/pS ratios to every promoter that has a WBID with associated pN/pS data
        try:
            promoter['pRatio'] = ratios[promoter['gene']]['ratio']
        except: # If there is no associated pN/pS value, just mark it as None
            promoter['pRatio'] = None

    with open(PI_FILE, 'r') as f: # As above, but with pi values
        data = json.loads(f.read())
    for promoter in promoters:
        try:
            promoter['pi'] = data[promoter['gene']]['piNorm']
        except:
            promoter['pi'] = None

    with open(dNdS_RATIO_FILE, 'r') as f: # Also as above, but with dN/dS values instead
        ratios = json.loads(f.read())
    for promoter in promoters:
        try:
            promoter['dRatio'] = ratios[promoter['gene']]['dRatio']
        except:
            promoter['dRatio'] = None

    uniquePromoters = []
    promoterNames = set()
    for promoter in promoters:
        if promoter['gene'] not in promoterNames:
            promoterNames.add(promoter['gene'])
            uniquePromoters.append(promoter)

    return uniquePromoters

def binGenes(genes): #* Bins all the genes into low and high pN/pS, pi and dN/dS (or neither, as it's just a TRUE/FALSE system)
    lowpNpS, highpNpS, lowpi, highpi, lowdNdS, highdNdS = [], [], [], [], [], []
    for gene in genes:
        gene['lowpNpS'], gene['highpNpS'], gene['lowpi'], gene['highpi'], gene['lowdNdS'], gene['highdNdS'] = False, False, False, False, False, False # Set all the 'is it in this bin' flags to false
        if gene['pRatio']: # If the gene has a pN/pS ratio associated with it
            if 0.05 <= gene['pRatio'] < 0.5: #? If the pN/pS ratio is between 0.05 and 0.5, the gene is considered strongly purifying 
                lowpNpS.append(gene)
                gene['lowpNpS'] = True
            elif 1.2 <= gene['pRatio'] < 2.5: #? If the pN/pS ratio is between 1.2 and 2.5, the gene is considered diversifying
                highpNpS.append(gene)
                gene['highpNpS'] = True
        if gene['pi']:
            if math.log10(gene['pi']) < -3.75: #? If the log of the pi value is less than -3.75, the gene is considered extremely purifying 
                lowpi.append(gene)
                gene['lowpi'] = True
            elif -2.75 <= math.log10(gene['pi']): #? If the log of the pi value is greater than -2.75, the gene is considered weakly purifying 
                highpi.append(gene)
                gene['highpi'] = True
        if gene['dRatio']:
            if 0 <= gene['dRatio'] < 0.03: #? If the dN/dS ratio is less than 0.03, the gene is considered extremely purifying 
                lowdNdS.append(gene)
                gene['lowdNdS'] = True
            elif 0.15 <= gene['dRatio'] < 0.25: #? If the dN/dS ratio is between 0.15 and 0.25, the gene is considered weakly purifying
                highdNdS.append(gene)
                gene['highdNdS'] = True
    
    return genes, lowpNpS, highpNpS, lowpi, highpi, lowdNdS, highdNdS

def combineBins(pNpS, pi, dNdS, promoterMidPos = False): #* Combines the three bins (either low or high) together. If promoters is true, then both the WBID and midPos (locus) must match. Else, only the WBID must match
    combinedBin = []
    geneNames = set()
    for gene in (pNpS + dNdS + pi):
        if promoterMidPos:
            if (gene['gene'], gene['midPos']) not in geneNames:
                geneNames.add((gene['gene'], gene['midPos']))
                combinedBin.append(gene)
        else:
            if gene['gene'] not in geneNames:
                geneNames.add(gene['gene'])
                combinedBin.append(gene)
    
    return combinedBin

maGenes = getLabelledMAGenes()
promoters = getLabelledPromoters()

maGenes, lowpNpS, highpNpS, lowpi, highpi, lowdNdS, highdNdS = binGenes(maGenes)
#- promoters, lowpNpS, highpNpS, lowpi, highpi, lowdNdS, highdNdS = binGenes(promoters)

#- lowBin = combineLowBins([promoter for promoter in promoters if promoter['lowpNpS'] == True], [promoter for promoter in promoters if promoter['lowpi'] == True], [promoter for promoter in promoters if promoter['lowdNdS'] == True], True)
lowBin = combineBins(lowpNpS, lowpi, lowdNdS)
#- highBin = combineHighBins([promoter for promoter in promoters if promoter['highpNpS'] == True], [promoter for promoter in promoters if promoter['highpi'] == True], [promoter for promoter in promoters if promoter['highdNdS'] == True], True)
highBin = combineBins(highpNpS, highpi, highdNdS)
#- print(len(lowBin), len(highBin))

def getOverlap(genes1, genes2, completeMatch = False): #* Calculates the number of genes that are the same between two groups. If completeMatch is true, then the entire dictionary entry must match. Else, only the WBIDs must match
    overlap = 0
    for gene1 in genes1:
        for gene2 in genes2:
            if completeMatch:
                if gene1 == gene2:
                    overlap += 1
            else:
                if gene1['gene'] == gene2['gene']:
                    overlap += 1
    
    return overlap

def outputOverlap(genes1, genes2): #* Prints the overlap in the genes between two groups
    overlap = getOverlap(genes1, genes2)
    print(f'Overlap of {overlap} ({round(overlap / (len(genes1) + len(genes2)) * 100, 2)}%)')

# For each pair of interest, calculate and output the overlap
pairs = [(lowpNpS, lowdNdS), (lowpNpS, highdNdS), (lowpNpS, lowpi), (lowpNpS, highpi), (highpNpS, lowdNdS), (highpNpS, highdNdS), 
         (highpNpS, lowpi), (highpNpS, highpi), (lowpi, lowdNdS), (lowpi, highdNdS), (highpi, lowdNdS), (highpi, highdNdS)]
for pair in pairs:
    outputOverlap(pair[0], pair[1])

with open(OUTPUT_FILE, 'w') as FoI: # Saves the data on the binned genes
    for gene in maGenes:
        FoI.write(json.dumps(gene) + '\n')