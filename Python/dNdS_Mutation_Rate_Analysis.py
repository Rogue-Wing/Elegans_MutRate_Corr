# Takes a series of dN/dS ratios covering a genome, and mutation data for the same genome, and highlights whether, if any, relationship exists between the two

import json
import requests
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project
MA_FILE_NAME = 'Combined MA Data v3.json'

SILENT = False

def progressBar(currentPercent, tick): #* Creates a progress bar as a string
    symbol = {0:'|', 1:'/', 2:'-', 3:'\\'}[tick]

    tick += 1
    if tick > 3:
        tick = 0

    return f"{'-' * (currentPercent)}{' ' * (100 - currentPercent)} [{currentPercent}%] {symbol}\r", tick

def getGeneLengths(): #* Gets the lengths of all the genes from the Ce11 Genome (an annotated genome)
    geneNames = {}
    with open(f'{PATH}/Annotated Ce11 Genome.gff', 'r') as FoI:
        for line in FoI:
            proteinCoding = False
            if line[0] != '#':
                data = line.split('\t')
                if data[2] == 'gene':
                    features = data[8].split(';')
                    for feature in features:
                        if 'protein_coding' in feature:
                            proteinCoding = True
                    
                    if proteinCoding:
                        if 'WormBase' in features[1]:
                            wormBaseID = features[1].split(',')[0].split(':')[1] # The WormBase ID, if present, always comes first
                            startPos, endPos = int(data[3]), int(data[4])
                            length = endPos - startPos
                            geneNames[wormBaseID] = length
    
    return geneNames

def ratioMAGenes(NSRatios, genes): #* Gets all the MA mutants, and assigns them a length and pN/pS ratio 
    MAGenes = {} # Dict of dicts ({<WBID>:{'length', 'numberOfMutHits', 'ratio', 'mutationType', 'chrom', 'locus'}})
    if not SILENT:
        with open(f'{PATH}/{MA_FILE_NAME}', 'r') as FoI:
            lineCount = 0
            for line in FoI:
                lineCount += 1
        
        linesPerPercent = lineCount // 100

        print(f'There are {lineCount} line(s) in the file! Now reading...')

    with open(f'{PATH}/{MA_FILE_NAME}', 'r') as FoI:
        errScore = 0

        if not SILENT:
            bar, tick = progressBar(0, 0) 
            print(bar, end='')

        for lineNum, line in enumerate(FoI): # For each mutant
            if line[0] != '#': # If a mutant has been commented (#'d) out, it's because it's a dead gene
                mutation = json.loads(line.strip()) # Reads in the mutant
                geneName = mutation['gene']
                if not '-' in geneName: # If there's a -, then the gene name is two genes (and the MA mutant comes from an intergenic region between the two). As the CaeNDR data has been santitised of intergenic sequences, this data is redundant
                    mutType = mutation['mutationType'].split('/')
                    if mutType[0] not in ['Modifier', 'NC', 'IG']: # Gets the type of mutation (in a more easily read format with less variation)
                        if '->' in mutType[1]:
                            if mutType[1][0] == mutType[1][-1]:
                                mutType[1] = 'Synonymous'
                            if '*' in mutType[1]:
                                if mutType[1][0] == '*':
                                    mutType[1] = 'Stop_Loss'
                                elif mutType[1][-1] == '*':
                                    mutType[1] = 'Nonsense'
                            elif mutType[1][0] == 'M':
                                mutType[1] = 'Start_Loss'
                            else:
                                mutType[1] = 'Missense'
                        mutType = f'{mutType[0]}/{mutType[1]}'
                        try: # Try to get the length of the gene from the annotated genome
                            try:
                                length = genes[geneName]
                            except KeyError: # If the length can't be found in the database, API request to wormbase to attempt to find a gene matching that WBID, and pull the start and end position data
                                location = requests.get(f"http://api.wormbase.org//rest/field/gene/{geneName}/genomic_position").json()['genomic_position']['data'][0]
                                startPos, endPos = location['start'], location['stop']
                                length = endPos - startPos
                            numberOfHits = 0
                            if geneName in MAGenes: # If the gene already exists in the list of mutants, then we're updating the entry, rather than creating a new one, so grab the data to make sure we don't override it
                                numberOfHits = MAGenes[geneName]['numberOfHits'] # Update the number of hits to what's currently stored
                                MAMutType = MAGenes[geneName]['mutationType']
                                if MAMutType != mutType: # If there's more than one mutation type for this gene, consider it to experience 'mixed' mutation
                                    mutType = 'Mixed'
                            numberOfHits += 1 # Add one hit, representing the number of hits (mutations) for that gene in the MA dataset 
                            ratio = NSRatios[geneName]['ratio'] # Add the pN/pS ratio for the gene
                            chrom, locus = mutation['locus'].split(':')
                            MAGenes[geneName] = {'length':length, 'numberOfHits':numberOfHits, 'dN':NSRatios[geneName]['count'][0], 'dS':NSRatios[geneName]['count'][1], 'ratio':ratio, 'mutationType':mutType, 'chromosome':chrom, 'locus':locus}
                        except KeyError: # The gene name does not exist in the dN/dS database
                            errScore += 1
                        except TypeError:
                            print(f'{geneName} has no data')
        
            if not SILENT:
                bar, tick = progressBar(lineNum // linesPerPercent, tick)
                print(bar, end='')
    
    return MAGenes

def ratioNSGenes(NSRatios, MAGenes, geneLengths): #* Assigns the pN/pS ratio to the genes of the AG dataset
    errs = 0
    if not SILENT:
        print(f"There are {len(NSRatios)} genes in the file! Now reading... {' ' * 50}")
        bar, tick = progressBar(0, 0)
        print(bar, end='')

        genesPerPercent = len(NSRatios) // 100

    for geneNum, geneName in enumerate(NSRatios): # For each gene
        if geneName in MAGenes: # If it's also in the MA dataset, add one mutation onto the experimental data's hits
            MAGenes[geneName]['numberOfHits'] += 1
            MAGenes[geneName]['MAGene'] = True
        else: # It's not, so we need to grab the same data as we did in ratioMAGenes (length, primarily)
            try:
                try: # Attempts to get the length data from the annotated genome. If not, make an API call to Wormbase
                    length = geneLengths[geneName]
                except KeyError:
                    location = requests.get(f"http://api.wormbase.org//rest/field/gene/{geneName}/genomic_position").json()['genomic_position']['data'][0]
                    startPos, endPos = location['start'], location['stop']
                    length = endPos - startPos
                MAGenes[geneName] = {'length':length, 'numberOfHits':1, 'dN':NSRatios[geneName]['count'][0], 'dS':NSRatios[geneName]['count'][1], 'ratio':NSRatios[geneName]['ratio'], 'mutationType':'Mixed', 'chromosome':NSRatios[geneName]['ratio'], 'locus':-1, 'MAGene':False} # Hits are set to one for all AG genes that aren't MA genes too
            except TypeError: # If the gene doesn't exist under that specific WBID, it maybe an out-dated ID. Try to get the history of the WBID and find an updated ID instead
                merged = False
                history = requests.get(f"http://api.wormbase.org//rest/field/gene/{geneName}/history").json()['history']['data']
                for timepoint in history:
                    if timepoint['action'].lower() == 'merged_into':
                        merged = True
                        try:
                            location = requests.get(f"http://api.wormbase.org//rest/field/gene/{timepoint['gene']['id']}/genomic_position").json()['genomic_position']['data'][0]
                            startPos, endPos = location['start'], location['stop']
                            length = endPos - startPos
                            MAGenes[geneName] = {'length':length, 'numberOfHits':1, 'dN':NSRatios[geneName]['count'][0], 'dS':NSRatios[geneName]['count'][1], 'ratio':NSRatios[geneName]['ratio'], 'mutationType':'Mixed', 'chromosome':NSRatios[geneName]['ratio'], 'locus':-1, 'MAGene':False}
                        except TypeError:
                            errs += 1
                            #- print(f'{geneName} has no data')
                if not merged:
                    errs += 1
                    #- print(f'{geneName} has no data')

        
        if not SILENT:
            bar, tick = progressBar(geneNum // genesPerPercent, tick)
            print(bar, end = '')
    
    if not SILENT:
        print(f"Done! {' ' * 110}")

    print(errs)

    return MAGenes 

with open(f'{PATH}/dNdS_Ratios v2.json', 'r') as FoI:
    NSRatios = json.loads(FoI.read())

geneLengths = getGeneLengths()

MAGenes = ratioMAGenes(NSRatios, geneLengths)
MAGenes = ratioNSGenes(NSRatios, MAGenes, geneLengths)

def saveToJSON(MAGenes): #* Saves the data to an output file
    with open(f'{PATH}/Combined NS-MA Raw Data v3.json', 'w') as FoI:
        FoI.write(json.dumps(MAGenes))

saveToJSON(MAGenes)
