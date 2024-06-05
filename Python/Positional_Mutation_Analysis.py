# Calculates the relative rate of mutation around the TTS and TSS

import json
import requests

import math
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project

MA_FILE = f'{PATH}/Combined MA Data v3.json' # MA mutants file
CAENDR_FILE = f'{PATH}/Collated CaeNDR Mutants v2.json' # CaeNDR (AG) SNPs

SELECTION_MEASURES_FILE = f'{PATH}/Combined Selection Measures v3.json' # All five selection measures for each gene
RATIO_FILE = f'{PATH}/dNdS_Ratios v2.renamed.json' # pN/pS ratios
PI_FILE = f'{PATH}/Combined NS-MA-TajD Raw Data v3.json' # pi values
ACTUAL_RATIO_FILE = f'{PATH}/Actual dNdS v2.renamed.json' # dN/dS ratios

LENGTHS_FILE = f'{PATH}/Ce11 Protein Coding Genes.txt' #! This gets the lengths that have been used in other analyses. While the length could be calculated by summing the nucleotide counts together, this often gets different total lengths, if only slightly
CODING_PROMOTERS = f'{PATH}/Ahringer Coding Promoter ATAC-seq Data.renamed.json' # A TSV file of all the genes' TTSs and TSSs, along with the WBIDs
TSS_POSITIONS = f'{PATH}/TTS and TSS Positions (All Genes).txt'

#- OUTPUT_FILE = f'{PATH}/MA Relative Extremely Purifying Coding Promoter Offsets (Nearest, dNdS).json' 
OUTPUT_FILE = f'{PATH}/MA Relative High Coding Promoter Offsets (Nearest).json' 

OFFSET = 10000 #? How many bp either side of the midpoint should be considered for a SNP to be counted
MIN_LENGTH = 3000
SILENT = False

NUMERAL_TO_CHROM = {'I':'chr1', 'II':'chr2', 'III':'chr3', 'IV':'chr4', 'V':'chr5', 'X':'chr10', 'MtDNA':'MtDNA'} # Converts numeral chromosome names to shorthand variants

def progressBar(currentPercent, tick): # Outputs a progress bar as a string
    symbol = {0:'|', 1:'/', 2:'-', 3:'\\'}[tick]

    tick += 1
    if tick > 3:
        tick = 0

    return f"{'-' * (currentPercent)}{' ' * (100 - currentPercent)} [{currentPercent}%] {symbol}\r", tick

def getLengths(): # Get all the lengths of the genes from the lengths file (gene start and end positions), along with the chromosome number
    lengths = {}
    with open(LENGTHS_FILE, 'r') as FoI:
        for lineNum, line in enumerate(FoI):
            if lineNum != 0:
                data = line.strip().split('\t')
                geneID = data[1]
                start, end = int(data[3]), int(data[4])
                length = end - start
                chromosome = NUMERAL_TO_CHROM[data[5]]

                lengths[geneID] = {'length':length, 'start':start, 'end':end, 'chrom':chromosome}
    
    with open(PI_FILE, 'r') as f: # For every gene in the database, and add the length to the list of lengths if it doesn't exist in the original file (i.e. if it was missing from LENGTHS_FILE)
        data = json.loads(f.read())
    for gene in data:
        if gene not in lengths:
            lengths[gene] = {'length':data[gene]['length'], 'start':None, 'end':None}
    
    return lengths

def getRatios(): #* Gets the pN/pS ratios to apply to the data (provided the gene ID can be found in the ratio file)
    with open(RATIO_FILE, 'r') as FoI:
        data = json.loads(FoI.read())
    
    ratios = {}
    for gene in data:
        ratios[gene] = data[gene]['ratio'] 
    return ratios

def getActualRatios(): #* Gets the dN/dS ratios to apply to the data (provided the gene ID can be found in the ratio file)
    with open(ACTUAL_RATIO_FILE, 'r') as FoI:
        data = json.loads(FoI.read())
    
    actualRatios = {}
    for gene in data:
        actualRatios[gene] = data[gene]['dRatio']
    return actualRatios

def getPis(): #* Gets the Pi values to apply to the data (provided the gene ID can be found in the ratio file)
    with open(PI_FILE, 'r') as FoI:
        data = json.loads(FoI.read())
    
    pis = {}
    for gene in data:
        pis[gene] = data[gene]['piNorm']
    return pis

def getSelectionMeasures(): #* Gets the pN/pS, pi and dN/dS ratios (IN ONE FILE) to apply to the data (provided the gene ID can be found in the file)
    with open(SELECTION_MEASURES_FILE, 'r') as FoI:
        data = json.loads(FoI.read())
    
    measures = {}
    for gene in data:
        measures[gene] = {'pRatio':data[gene]['pRatio'], 'pi':data[gene]['pi'], 'dRatio':data[gene]['dRatio']}
    return measures

def getAllCodingPromoters(): #* Gets all the coding promoters from the Ahringer file, specifically by grabbing the gene ID that the promoter is for, as well as the midpoint of the promoter and the strand its on 
    genePositions = []
    with open(CODING_PROMOTERS, 'r') as FoI:
        for lineNum, line in enumerate(FoI):
            data = json.loads(line.strip()) # Each line is a dictionary, so load it via JSON parsing
            geneID = data['geneID'] # The gene that the promoter is associated with
            pos = (int(data['stop']) + int(data['start'])) // 2 # Find the midpoint between the start and end of the coding promoter, rounded down
            strand = data['strand'] # The strand that the promoter is on
            #- if strand == '1': #! All the other data is forward-strand exclusive, so this should be too
            genePositions.append({'promoterID':geneID, 'midPos':pos, 'strand':strand}) # Add the promoter to the list of all promoter location-strand entries
    
    return genePositions

def getAllTSSs(): #* Gets all the TSSs from the Biomart file, grabbing the geneID, position and strand data 
    genePositions = []
    with open(TSS_POSITIONS, 'r') as FoI:
        for lineNum, line in enumerate(FoI):
            if lineNum != 0:
                data = line.strip().split('\t') # Each line consists of tab-separated data, so split it into a list of elements [source, TSS, TTS, geneStart, geneEnd, strand]
                geneID = data[1]
                pos = int(data[2])
                strand = int(data[6])
                genePositions.append({'promoterID':geneID, 'midPos':pos, 'strand':strand})
    
    return genePositions

def getNearbyGenes(genePositions, pos, geneID): #* For a single position, gets all coding promoter regions/TSSs that contain the position within the offset boundaries. A negative value is upstream, a positive one is downstream (on the positive [1] strand - this is inverted on the -1 strand)
    nearbyGenes = []
    for gene in genePositions: # For each gene in the list of gene positions (usually either the TSS or coding promoter regions)
        if (gene['midPos'] - OFFSET) <= pos <= (gene['midPos'] + OFFSET): # If the position of the mutation falls within the offset region around the gene
            nearbyGenes.append({'geneID':geneID, 'distance':pos - gene['midPos'], 'strand':gene['strand'], 'locus':pos}) # Add it to the list
    
    return nearbyGenes

def getAllNearbyGenes(genePositions, genes): #* For each mutation/SNP in the dataset (MA or CaeNDR), take its position and find any coding promoter regions/TSSs that it is close (within the specified offset) too  
    allNearbyGenes = []
    if not SILENT:
        bar, tick = progressBar(0, 0)
        print(bar, end = '')
        mutsPerPercent = len(genes) // 100

    missingGenes = 0
    for mutNum, mutant in enumerate(genes): # For each gene in the dataset
        nearbyGenes = getNearbyGenes(genePositions, int(mutant['orderingData'][1]), mutant['gene']) # Get all the coding promoter regions/TSSs that are nearby to the position of the gene
        if nearbyGenes == []: # If there are no genes within the specified offset
            missingGenes += 1
        allNearbyGenes += nearbyGenes # Add the locations for that gene to the list of all locations

        if not SILENT:
            bar, tick = progressBar(int(mutNum / mutsPerPercent), tick)
            print(bar, end = '')
    
    if not SILENT:
        print(f"Done! {95 * ' '}")

    return allNearbyGenes, missingGenes 

def getClosestGene(genePositions, pos, geneID): #* For a single position, gets the CLOSEST coding promoter region/TSS to it, provided that it falls within the offset region
    closestGene = None
    for gene in genePositions: # For each gene in the list of gene positions (usually either the TSS or coding promoter regions) 
        if (gene['midPos'] - OFFSET) <= pos <= (gene['midPos'] + OFFSET): # If the position falls within the region (+- the offset)
            distance = pos - gene['midPos'] # Calculate the distance between the position and the centre of the region/the TSS. Positive values mean that the position is downstream, negative values mean upstream
            if closestGene: # If there's already a gene recorded
                if abs(distance) < abs(closestGene['distance']): #! If the newly calculated difference is less than the currently stored distance, replace the stored distance with the new (lower) difference
                    closestGene = {'geneID':geneID, 'promoterGeneID':gene['promoterID'], 'distance':distance, 'strand':gene['strand'], 'locus':pos}
            else: #? There's no gene already stored, so just store this one (some distance is less than an infinite distance)
                closestGene = {'geneID':geneID, 'promoterGeneID':gene['promoterID'], 'distance':distance, 'strand':gene['strand'], 'locus':pos}
    
    return closestGene

def getAllClosestGenes(genePositions, genes): #* For each mutation/SNP in the dataset (MA or CaeNDR), take its position and find the closest coding promoter regions/TSSs to it (provided it's within the offset region) 
    allClosestGenes = []
    if not SILENT:
        bar, tick = progressBar(0, 0)
        print(bar, end = '')
        mutsPerPercent = len(genes) // 100

    missingGenes = 0
    for mutNum, mutant in enumerate(genes): # For each gene in the dataset
        closestGene = getClosestGene(genePositions, int(mutant['orderingData'][1]), mutant['gene']) # Get the closest coding promoter region/TSS to the gene
        if not closestGene: # If there are no genes within the specified offset
            missingGenes += 1
        else:
            allClosestGenes.append(closestGene) # Add that location for that gene to the list of all locations

        if not SILENT:
            bar, tick = progressBar(int(mutNum / mutsPerPercent), tick)
            print(bar, end = '')
    
    if not SILENT:
        print(f"Done! {95 * ' '}")

    return allClosestGenes, missingGenes 

def getMAMutants(): #* Gets all the mutations from the MA dataset
    maMutants = []
    with open(MA_FILE, 'r') as FoI:
        for line in FoI:
            if line[0] != '#':
                newMutant = json.loads(line.strip())
                present = False
                for mutant in maMutants:
                    if mutant['locus'] == newMutant['locus'] and mutant['substitution'] == newMutant['substitution']:
                        present = True
                if not present:
                    maMutants.append(newMutant)
    
    return maMutants

def getCaeNDRMutants(): #* Gets all the SNPs from the CaeNDR dataset
    CaeNDRMutants = []
    with open(CAENDR_FILE, 'r') as FoI:
        for line in FoI:
            if line[0] != '#':
                CaeNDRMutants.append(json.loads(line.strip())) #! No need to do a check for duplicates, as these data have already been checked for that (any duplicates were summed into the 'strainNum' value)

    return CaeNDRMutants

def addRatios(ratios, mutants): #* For each mutant, the dN/dS ratio is added to the dictionary if its gene ID can be found in the list of dN/dS ratios
    for mutant in mutants:
        try:
            mutant['ratio'] = ratios[mutant['gene']]
        except KeyError:
            mutant['ratio'] = None
    
    return mutants

def addLength(lengths, mutants): #* For each mutant, add the length to the dictionary (if its gene ID can be found in the list of lengths)
    for mutant in mutants:
        try:
            mutant['length'] = lengths[mutant['gene']]['length']
        except KeyError:
            mutant['length'] = None
    
    return mutants

def getMidpointMutants(lengths, mutants): # Calculate the midpoint of the genes
    for gene in lengths:
        if (lengths[gene]['start'] != None) and (lengths[gene]['end'] != None):
            midpoint = (lengths[gene]['start'] + lengths[gene]['end']) / 2
            mutants.append({'gene':gene, 'orderingData':[lengths[gene]['chrom'], midpoint]})
    
    return mutants

def filterByLength(mutants): #* Remove any mutants below the length specified by MIN_LENGTH
    filteredMutants = []
    for mutant in mutants:
        if mutant['length']: #? If the length is not None
            if mutant['length'] >= MIN_LENGTH:
                filteredMutants.append(mutant)
    
    return filteredMutants

def binPromotersByRatio(ratios, promoters, binLowRange, binHighRange, silent = False): #* Takes the list of promoters, with their geneIDs, and bins them based on the pN/pS value of that gene
    errs = 0
    for promoter in promoters:
        try: # Try to get the pN/pS ratio for that WBID. If it doesn't exist, then consider it to be None
            promoter['ratio'] = ratios[promoter['promoterID']]
        except KeyError:
            errs += 1
            promoter['ratio'] = None
    
    if not silent:
        print(f'Of {len(promoters)} promoter(s), {len(promoters) - errs} have had pN/pS ratios assigned ({round((len(promoters) - errs) / len(promoters) * 100, 3)}% coverage)')

    binnedPromoters = []
    for promoter in promoters:
        if promoter['ratio']: # If the promoter actually has a ratio
            if binLowRange <= promoter['ratio'] < binHighRange: # If the ratio falls within the desired range for the bin
                binnedPromoters.append(promoter)
    
    if not silent:
        print(f'Of the {len(promoters) - errs} promoter(s) with pN/pS ratios, {len(binnedPromoters)} have ratios that fall within the range of {binLowRange} <= pN/pS < {binHighRange} ({round(len(binnedPromoters) / (len(promoters) - errs) * 100, 3)}%)')

    return binnedPromoters

def binPromotersByActualRatio(actualRatios, promoters, binLowRange, binHighRange, silent = False): #* Takes the list of promoters, with their geneIDs, and bins them based on the dN/dS value of that gene
    errs = 0
    for promoter in promoters:
        try: # Try to get the dN/dS ratio for that WBID. If it doesn't exist, then consider it to be None
            promoter['ratio'] = actualRatios[promoter['promoterID']]
        except KeyError:
            errs += 1
            promoter['ratio'] = None
    
    if not silent:
        print(f'Of {len(promoters)} promoter(s), {len(promoters) - errs} have had dN/dS ratios assigned ({round((len(promoters) - errs) / len(promoters) * 100, 3)}% coverage)')

    binnedPromoters = []
    for promoter in promoters:
        if promoter['ratio']: # If the promoter actually has a ratio
            if binLowRange <= promoter['ratio'] < binHighRange: # If the ratio falls within the desired range for the bin
                binnedPromoters.append(promoter)
    
    if not silent:
        print(f'Of the {len(promoters) - errs} promoter(s) with dN/dS ratios, {len(binnedPromoters)} have ratios that fall within the range of {binLowRange} <= dN/dS < {binHighRange} ({round(len(binnedPromoters) / (len(promoters) - errs) * 100, 3)}%)')

    return binnedPromoters

def binPromotersByPi(pis, promoters, binLowRange, binHighRange, silent = False): #* Takes the list of promoters, with their geneIDs, and bins them based on the pi value of that gene
    errs = 0
    for promoter in promoters:
        try: # Try to get the pi value for that WBID. If it doesn't exist, then consider it to be None
            promoter['pi'] = pis[promoter['promoterID']]
        except KeyError:
            errs += 1
            promoter['pi'] = None
    
    if not silent:
        print(f'Of {len(promoters)} promoter(s), {len(promoters) - errs} have had pi values assigned ({round((len(promoters) - errs) / len(promoters) * 100, 3)}% coverage)')

    binnedPromoters = []
    for promoter in promoters:
        if promoter['pi']: # If the promoter actually has a ratio
            if binLowRange <= math.log10(promoter['pi']) < binHighRange: # If the ratio falls within the desired range for the bin (note that the pi values get logged)
                binnedPromoters.append(promoter)
    
    if not silent:
        print(f'Of the {len(promoters) - errs} promoter(s) with pi values, {len(binnedPromoters)} have ratios that fall within the range of {binLowRange} <= pi < {binHighRange} ({round(len(binnedPromoters) / (len(promoters) - errs) * 100, 3)}%)')

    return binnedPromoters

def combinePromoters(ratios, pis, actualRatios, genePromoters): #* Creates two lists of promoters, low and high, which are all the non-overlapping promoters from the pN/pS, pi AND dN/dS bins
    lowRatioBin = binPromotersByRatio(ratios, genePromoters, 0.05, 0.5, True) # Gets the pN/pS low and high bins
    highRatioBin = binPromotersByRatio(ratios, genePromoters, 1.2, 2.5, True)

    lowPiBin = binPromotersByPi(pis, genePromoters, float('-inf'), -3.75, True)  # Gets the pi low and high bins
    highPiBin = binPromotersByPi(pis, genePromoters, -2.75, float('inf'), True)

    lowActualRatioBin = binPromotersByActualRatio(actualRatios, genePromoters, 0, 0.03, True)  # Gets the dN/dS low and high bins
    highActualRatioBin = binPromotersByActualRatio(actualRatios, genePromoters, 0.15, 0.25, True)

    uniquePromoters = set()
    lowBin = []
    for newPromoter in lowRatioBin + lowPiBin + lowActualRatioBin: # Gets all the unique promoters in the low bins
        if (newPromoter['promoterID'], newPromoter['midPos']) not in uniquePromoters:
            uniquePromoters.add((newPromoter['promoterID'], newPromoter['midPos']))
            lowBin.append(newPromoter)
    
    print(f"There are {len(lowBin)} promoters in the low bin, out of a total of {len(lowRatioBin + lowPiBin + lowActualRatioBin)} possible promoters.")

    uniquePromoters = set()
    highBin = []
    for newPromoter in highRatioBin + highPiBin + highActualRatioBin: # Gets all the unique promoters in the high bins
        if (newPromoter['promoterID'], newPromoter['midPos']) not in uniquePromoters:
            uniquePromoters.add((newPromoter['promoterID'], newPromoter['midPos']))
            highBin.append(newPromoter)
    
    print(f"There are {len(highBin)} promoters in the high bin, out of a total of {len(highRatioBin + highPiBin + highActualRatioBin)} possible promoters.")

    uniqueLowPromoterIDs = set([promoter['promoterID'] for promoter in lowBin])
    uniqueHighPromoterIDs = set([promoter['promoterID'] for promoter in highBin])
    
    # Calculates the promoters that exist in both the low-combined and high-combined bins, and removes them from both to leave only the unique promoters for each bin
    promotersToRemove = []
    overlaps = 0
    for promoter in uniqueLowPromoterIDs:
        if promoter in uniqueHighPromoterIDs:
            promotersToRemove.append(promoter)
            overlaps += 1
    
    for promoter in promotersToRemove:
        uniqueLowPromoterIDs.remove(promoter)
        uniqueHighPromoterIDs.remove(promoter)
    
    uniqueLowBin, uniqueHighBin = [], []
    for promoter in lowBin:
        if promoter['promoterID'] in uniqueLowPromoterIDs:
            uniqueLowBin.append(promoter)
    for promoter in highBin:
        if promoter['promoterID'] in uniqueHighPromoterIDs:
            uniqueHighBin.append(promoter)
    
    print(f'There is an overlap of {overlaps} promoters, with {len(uniqueLowBin)} low promoters (across {len(uniqueLowPromoterIDs)} genes), and {len(uniqueHighBin)} high promoters (across {len(uniqueHighPromoterIDs)} genes).')

    return uniqueLowBin, uniqueHighBin

ratios = getRatios()
pis = getPis()
actualRatios = getActualRatios()
measures = getSelectionMeasures()
lengths = getLengths()

mutants = getMAMutants()
#- mutants = getCaeNDRMutants()
genePromoters = getAllCodingPromoters()
geneTSSs = getAllTSSs()

mutants = addLength(lengths, mutants)
#- mutants = filterByLength(mutants)

BINNED = False

if BINNED: # If trying to get the binned files, run through each of the measures independently, and then the combined one. In each case, calculate the promoters for the bin, then get the closest promoter for each SNP from that list. Save the result
    runs = [('pN/pS', 'Strongly Purifying', 'Diversifying'), ('pi', 'Extremely Purifying', 'Weakly Purifying'), ('dN/dS', 'Extremely Purifying', 'Weakly Purifying'), ('Combined', 'Low', 'High')]
    for run in runs:
        if run[0] == 'pN/pS':
            lowBinPromoters = binPromotersByRatio(ratios, genePromoters, 0.05, 0.5)
            highBinPromoters = binPromotersByRatio(ratios, genePromoters, 1.2, 2.5)
        elif run[0] == 'pi':
            lowBinPromoters = binPromotersByPi(pis, genePromoters, float('-inf'), -3.75)
            highBinPromoters = binPromotersByPi(pis, genePromoters, -2.75, float('inf'))
        elif run[0] == 'dN/dS':
            lowBinPromoters = binPromotersByActualRatio(actualRatios, genePromoters, 0, 0.03)
            highBinPromoters = binPromotersByActualRatio(actualRatios, genePromoters, 0.15, 0.25)
        elif run[0] == 'Combined':
            lowBinPromoters, highBinPromoters = combinePromoters(ratios, pis, actualRatios, genePromoters)
        else:
            print('ERR: Selection measure unknown!')     
        
        allNearestGenes, missingGenes = getAllClosestGenes(lowBinPromoters, mutants)
        print(f'{run[1]}: {len(allNearestGenes)} data point(s) acquired! Out of {len(mutants)} SNPs, {missingGenes} were missed - coverage of {round((len(mutants) - missingGenes) / len(mutants) * 100, 3)}%')

        with open(f"{PATH}/MA Relative {run[1]} Coding Promoter Offsets (Nearest, {run[0].replace('/', '')}).json", 'w') as FoI:
            for position in allNearestGenes:
                FoI.write(json.dumps(position) + '\n')
        
        allNearestGenes, missingGenes = getAllClosestGenes(highBinPromoters, mutants)
        print(f'{run[2]}: {len(allNearestGenes)} data point(s) acquired! Out of {len(mutants)} SNPs, {missingGenes} were missed - coverage of {round((len(mutants) - missingGenes) / len(mutants) * 100, 3)}%')

        with open(f"{PATH}/MA Relative {run[2]} Coding Promoter Offsets (Nearest, {run[0].replace('/', '')}).json", 'w') as FoI:
            for position in allNearestGenes:
                FoI.write(json.dumps(position) + '\n')
        
        print('')
else: # Else, you're just calculating from the entire promoter dataset (no binning required). Do the same thing against this larger list of promoters
    allNearestGenes, missingGenes = getAllClosestGenes(genePromoters, mutants)
    print(f'ALL: {len(allNearestGenes)} data point(s) acquired! Out of {len(mutants)} SNPs, {missingGenes} were missed - coverage of {round((len(mutants) - missingGenes) / len(mutants) * 100, 3)}%')
    with open(OUTPUT_FILE, 'w') as FoI:
        for position in allNearestGenes:
            FoI.write(json.dumps(position) + '\n')
