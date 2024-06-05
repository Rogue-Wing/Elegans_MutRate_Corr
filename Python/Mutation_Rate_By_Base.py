# Takes the coding sequence of every gene in the C. elegans genome, the mutants/polymorphisms and pN/pS, pi and dN/dS values to calculate the mutation rate of each nucleotide
#! Note that it adds one to every counted mutation (i.e. a gene with a single A > G ([1, 0, 0, 0]) mutation would have [2, 1, 1, 1] for [A, C, G, T] mutations)

import json
import requests
from datetime import datetime
from datetime import timedelta
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project

CaeNDR_MUTANTS_FILE = f'{PATH}/Collated CaeNDR Mutants v2.json' # The file for SNPs
MA_MUTANTS_FILE = f'{PATH}/Combined MA Data v3.json' # The file for MA mutants
RATIO_FILE = f'{PATH}/dNdS_Ratios v2.renamed.json' # The file of pN/pS ratios
SELECTION_MEASURE_FILE = f'{PATH}/Combined Selection Measures v3.json' # The file of all selection measure data
POSITIONS_FILE = f'{PATH}/Ce11 Protein Coding Genes.txt' #! This gets the lengths that have been used in other analyses. While the length could be calculated by summing the nucleotide counts together, this often gets different total lengths, if only slightly
GENOME = f'{PATH}/WBcel235 Genome.json'

NUMERAL_TO_CHROM = {'I':'chr1', 'II':'chr2', 'III':'chr3', 'IV':'chr4', 'V':'chr5', 'X':'chr10', 'MtDNA':'MtDNA'} # Converts numeral chromosome names to shorthand variants
COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a'} # Converts nucleotides owing to the complementary base pairing rule

OUTPUT_FILE = f'{PATH}/NEW Individual CaeNDR Nucleotide Mutation Rates (+1).json'
ALL_GENES = True # Whether or not to consider the entire AG dataset, or just the MA genes

def progressBar(currentPercent, tick, withTime = False, currentTime = None, previousPercent = None, timePerPercent = None): #* Outputs a progress bar, and ETA, as a string
    symbol = {0:'|', 1:'/', 2:'-', 3:'\\'}[tick]

    tick += 1
    if tick > 3:
        tick = 0
    
    if withTime:
        if currentPercent > previousPercent:
            timePerPercent = ((timePerPercent * previousPercent) + (datetime.now() - currentTime)) / currentPercent
            currentTime = datetime.now()
        timeLeft = timePerPercent * (100 - currentPercent)
        return f"{'-' * (currentPercent)}{' ' * (100 - currentPercent)} [{currentPercent}%] {symbol} <ETA: {timeLeft}>\r", tick, currentTime, timePerPercent 
            
    return f"{'-' * (currentPercent)}{' ' * (100 - currentPercent)} [{currentPercent}%] {symbol}\r", tick

def getMAMutants(): #* Reads in the MA mutants
    mutants = []
    with open(MA_MUTANTS_FILE, 'r') as FoI:
        for line in FoI:
            if line[0] != '#': #? If the entry hasn't been commented out (likely because it's a dead gene)
                mutants.append(json.loads(line.strip())) # Add the mutant data to the list
    
    return mutants

def getCaeNDRMutants(): #* Reads in the SNPs from the CaeNDR database
    mutants = []
    with open(CaeNDR_MUTANTS_FILE, 'r') as FoI:
        for line in FoI:
            if line[0] != '#': #? If the entry hasn't been commented out (likely because it's a dead gene)
                mutants.append(json.loads(line.strip())) # Add the mutant data to the list
    
    return mutants

def getRatios(): #* Gets the pN/pS ratios to apply to the data (provided the gene ID can be found in the ratio file)
    with open(RATIO_FILE, 'r') as FoI:
        ratios = json.loads(FoI.read())
    return ratios

def getMeasures(): #* Gets the pN/pS, pi and dN/dS values to apply to the data (provided the gene ID can be found in the ratio file)
    with open(SELECTION_MEASURE_FILE, 'r') as FoI:
        measures = json.loads(FoI.read())
    return measures

def getGenome(): #* Reads in the N2 reference genome
    with open(GENOME, 'r') as FoI:
        return json.loads(FoI.read())

def assignSelectionMeasures(selectionMeasures, mutant, geneID): #* Assigns the pN/pS, pi and dN/dS values to a gene, if it has them
    try:
        mutant['pRatio'] = selectionMeasures[geneID]['pRatio']
    except KeyError:
        mutant['pRatio'] = None
    
    try:
        mutant['pi'] = selectionMeasures[geneID]['pi']
    except KeyError:
        mutant['pi'] = None
    
    try:
        mutant['dRatio'] = selectionMeasures[geneID]['dRatio']
    except KeyError:
        mutant['dRatio'] = None
    
    return mutant

def assignMutantRatio(ratios, maMutant, geneID): #* Assigns the pN/pS ratio to a gene, if it has one
    try:
        maMutant['ratio'] = ratios[geneID]['ratio']
    except KeyError:
        maMutant['ratio'] = None
    
    return maMutant

def getPositions(): #* Reads in the start and end positions, and length data, for each gene
    positions = {}
    with open(POSITIONS_FILE, 'r') as FoI:
        for lineNum, line in enumerate(FoI):
            if lineNum != 0:
                data = line.strip().split('\t')
                geneID = data[1]
                strand = data[2]
                start, end = int(data[3]), int(data[4])
                length = end - start
                chromosome = NUMERAL_TO_CHROM[data[5]]

                positions[geneID] = {'start':start, 'end':end, 'length':length, 'strand':strand, 'chromosome':chromosome}
    
    return positions

def assignMutantLength(lengths, geneID): #* Gets the length for a gene
    try:
        length = lengths[geneID]
    except KeyError:
        length = None
    
    return length

def assignSequenceWithoutMutations(genome, positions, mutants): #* Grabs the sequence for every gene, discarding any information pertaining to individual mutations/SNPs 
    bar, tick = progressBar(0, 0)
    print(f'{bar}', end = '')
    mutsPerPercent = len(mutants) // 100
    currentTime = datetime.now()
    timePerPercent = timedelta(seconds = 0)
    
    errs = 0
    genes = {}
    for mutNum, mutant in enumerate(mutants):
        if '-' not in mutant['gene']: #! If it's not an intergenic region
            try:
                geneData = positions[mutant['gene']] # Grab the data on the start, end, length, strand and chromosome for that gene
                seq = genome[geneData['chromosome']][geneData['start'] - 1:geneData['end']]
                if geneData['strand'] == '-1': #! If the strand is negative, then the complement of the genome (which is of the forward strand) needs to be taken
                    compSeq = ''
                    for nuc in seq:
                        compSeq += COMPLEMENT[nuc]
                    seq = compSeq
                if mutant['gene'] not in list(genes.keys()):
                    genes[mutant['gene']] = {'mutations':[], 'sequence':seq}
            except KeyError: # If the WBID cannot be found in the set of coding genes, then either it's not a coding gene, or the name used by the MA data is out-of-date
                errs += 1

        bar, tick, currentTime, timePerPercent = progressBar(int(mutNum / mutsPerPercent), tick, True, currentTime, int((mutNum - 1) / mutsPerPercent), timePerPercent)
        print(f'{bar}', end = '')
        
    print(f'Done! {" " * 95}')       

    return genes, errs

def assignSequenceToMutants(genome, positions, maMutants): #* Grabs the sequence for every gene, and compresses all mutations of the same gene into a single data entry (the 'mutations' list)
    bar, tick = progressBar(0, 0)
    print(f'{bar}', end = '')
    mutsPerPercent = len(maMutants) // 100

    errs = 0
    mutants = {}
    for mutNum, mutant in enumerate(maMutants):
        if '-' not in mutant['gene']: #! If it's not an intergenic region
            try:
                geneData = positions[mutant['gene']] # Grab the data on the start, end, length, strand and chromosome for that gene
                seq = genome[geneData['chromosome']][geneData['start'] - 1:geneData['end']]
                if geneData['strand'] == '-1': #! If the strand is negative, then the complement of the genome (which is of the forward strand) needs to be taken
                    compSeq = ''
                    for nuc in seq:
                        compSeq += COMPLEMENT[nuc]
                    seq = compSeq
                if mutant['gene'] not in list(mutants.keys()):
                    mutants[mutant['gene']] = {'mutations':[{'locus':mutant['locus'], 'substitution':mutant['substitution'], 'mutationType':mutant['mutationType']}], 'sequence':seq}
                else: # The WBID is already being used as a key - this is another mutation present within the same gene
                    duplicate = False
                    for storedMutant in mutants[mutant['gene']]['mutations']:
                        if (storedMutant['locus'] == mutant['locus']) and (storedMutant['substitution'] == mutant['substitution']):
                            duplicate = True
                    if not duplicate:
                        mutants[mutant['gene']]['mutations'].append({'locus':mutant['locus'], 'substitution':mutant['substitution'], 'mutationType':mutant['mutationType']})
            except KeyError: # If the WBID cannot be found in the set of coding genes, then either it's not a coding gene, or the name used by the MA data is out-of-date
                errs += 1
        
        bar, tick = progressBar(int(mutNum / mutsPerPercent), tick)
        print(f'{bar}', end = '')
    print(f'Done! {" " * 95}')
    
    return mutants, errs

def countNucleotides(sequence, desiredNuc): #* Counts the number of a particular nucleotide (A, C, G or T) present in a sequence
    count = 0
    for nuc in sequence:
        if nuc.lower() == desiredNuc.lower():
            count += 1
    
    return count

def calculateNucleotideMutationRates(mutant, geneID, count, rate = True, plusOne = False, fullLength = False, lengths = None): #* Calculates the mutation rate (if rate is true) for each nucleotide in that particular gene. If fullLength is true, then the rate is calculated depending on the total number of nucleotides, not just that specific type of nucleotide (i.e. A + C + G + T, rather than just A). The lengths parameter provides a dict for getting the length of the gene - if empty, the length is calculated by summing the nucleotide counts
    if plusOne:
        nucleotideMutations = {'A':[1, countNucleotides(mutant['sequence'], 'A')], 'C':[1, countNucleotides(mutant['sequence'], 'C')], 'G':[1, countNucleotides(mutant['sequence'], 'G')], 'T':[1, countNucleotides(mutant['sequence'], 'T')]} #! Initialises the data structure to store the mutation counts. Starts with 1 in all cases - this is like doing +1 to the MA hits
    else:
        nucleotideMutations = {'A':[0, countNucleotides(mutant['sequence'], 'A')], 'C':[0, countNucleotides(mutant['sequence'], 'C')], 'G':[0, countNucleotides(mutant['sequence'], 'G')], 'T':[0, countNucleotides(mutant['sequence'], 'T')]} #! Initialises the data structure to store the mutation counts. Starts with 0 in all cases - this is like testing *just* the mutant hits
    for mutation in mutant['mutations']: # For each mutation of the gene
        nuc = mutation['substitution'][0] # Get the nucleotide that it originally was 
        count += 1 # Add one to the total number of mutations/polymorphisms counted
        nucleotideMutations[nuc][0] += 1 #? Add one to the number of recorded mutations of that nucleotide in that gene
    if fullLength:
        if lengths: #! If a lengths dict has been provided (in the form {geneID:length})
            length = assignMutantLength(lengths, geneID)
            if not length: #! If a length could not be found in the database, simply get it from adding the nucleotide sequences
                length = nucleotideMutations['A'][1] + nucleotideMutations['C'][1] + nucleotideMutations['G'][1] + nucleotideMutations['T'][1]
        else: #! If no length dict is provided, simply calculate the length by adding the nucleotide sequences 
            length = nucleotideMutations['A'][1] + nucleotideMutations['C'][1] + nucleotideMutations['G'][1] + nucleotideMutations['T'][1]
    for nucleotide in nucleotideMutations:
        if nucleotideMutations[nucleotide]: # If it's not None (i.e. it has a mutation associated with it)
            if rate:
                if fullLength:
                    nucleotideMutations[nucleotide] = nucleotideMutations[nucleotide][0] / length # Calculate the mutation rate compared to the length of the ENTIRE gene
                else:
                    nucleotideMutations[nucleotide] = nucleotideMutations[nucleotide][0] / nucleotideMutations[nucleotide][1] # Calculate the mutation rate as though the gene only had one nucleotide type
            else:
                nucleotideMutations[nucleotide] = [nucleotideMutations[nucleotide][0], nucleotideMutations[nucleotide][1]] # Gets the number of mutations, and the total number of that nucleotide in the sequence, and stores it

    return nucleotideMutations, count

def calculateIndividualMutationRates(mutants, ratios, plusOne = False, rate = True): #* Takes a list of mutants and runs calculateNucleotideMutationRates on each of them, collating the data into a single dict entry before assigning the selection measures
    nucleotideMutations = {}
    count = 0

    bar, tick = progressBar(0, 0)
    print(f'{bar}', end = '')
    mutsPerPercent = len(mutants) // 100

    for mutNum, mutant in enumerate(mutants): # For each mutant/SNP, calculate the mutation rate and assign the pN/pS, pi and dN/dS values
        nucleotideMutations[mutant], count = calculateNucleotideMutationRates(mutants[mutant], mutant, count, rate, plusOne)
        nucleotideMutations[mutant] = assignSelectionMeasures(ratios, nucleotideMutations[mutant], mutant)

        bar, tick = progressBar(int(mutNum / mutsPerPercent), tick)
        print(f'{bar}', end = '')
    print(f'Done! {95 * " "}')
    
    print(count)
    
    return nucleotideMutations

selectionMeasures = getMeasures()
positions = getPositions()
genome = getGenome()

maMutants = getMAMutants()
allGenesCount = len(maMutants)
mutants, missingGenes = assignSequenceToMutants(genome, positions, maMutants) # Get the sequences for all MA mutants

if ALL_GENES: # If the CaeNDR SNPs are also being considered
    AGPolyMorphs = getCaeNDRMutants()
    caeNDRGenes, CaeNDRMissingGenes = assignSequenceWithoutMutations(genome, positions, AGPolyMorphs) # Assign the sequence to the SNPs
    mutants = {**caeNDRGenes, **mutants} # Merge the two dicts of sequences
    missingGenes += CaeNDRMissingGenes
    allGenesCount = len(mutants)

nucleotideMutations = calculateIndividualMutationRates(mutants, selectionMeasures, True) # Calculate the nucleotide-specific mutation rates

with open(OUTPUT_FILE, 'w') as FoI: #* Save the nucleotide-specific mutation rates to an output .json file
    FoI.write(json.dumps(nucleotideMutations))
