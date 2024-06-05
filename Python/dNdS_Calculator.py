# This takes a .json mutation data file and calculates dN/dS ratios from it

import json
import re
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project

AA_THREE_TO_ONE = {'ala':'A', 'arg':'R', 'asn':'N', 'asp':'D', 'cys':'C', 'glu':'E', 'gln':'Q', 'gly':'G', 'his':'H', 'ile':'I', 
                  'leu':'L', 'lys':'K', 'met':'M', 'phe':'F', 'pro':'P', 'ser':'S', 'thr':'T', 'trp':'W', 'tyr':'Y', 'val':'V', 'ter':'*', 'unk':'?'}

def progressBar(currentPercent, tick): #* Outputs a progress bar as a string
    symbol = {0:'|', 1:'/', 2:'-', 3:'\\'}[tick]

    tick += 1
    if tick > 3:
        tick = 0

    return f"{'-' * (currentPercent)}{' ' * (100 - currentPercent)} [{currentPercent}%] {symbol}\r", tick

def getRatios(allMutants, cumulativeMutations = {}, silent = False): #* Takes all the SNPs, extracts any that are present within protein-coding stretches of the DNA, then sums (for each gene) the number of synonymous and non-synonymous SNPs to calculate a pN/pS ratio
    if not silent:
        linesPerPercent = len(allMutants) // 100

        bar, tick = progressBar(0, 0)
        print(bar, end = '')

    mutantsCheck = set() # Set of all SNPs. As its a set, they must be unique (i.e. cannot be the same gene name, locus *and* mutationType (synon or non-synon))
    uniqueMutants = {}
    for mutantNum, mutant in enumerate(allMutants):
        # Get the SNP (JSON-ised). This takes the form of {locus, gene, substitution, mutationType, orderingData}. Important here are the 'gene' (the name of the gene <GENE_NAME (ALLELE_NAME)>, although the allele is not always present) and the 'mutationType' (usually either Exon/Synonymous or Exon/<AA_SUB>)
        chromosome = mutant['locus'].split(':')[0] # Split the locus to only have the chromosome number (as a roman numeral/mtDNA)
        locus = mutant['locus'].split(':')[1] # Split the locus to only have the genomic position 
        geneName = mutant['gene'].split(' ')[0] # Split the gene to only have the gene name (remove the variant/allele name)
        
        mutType = mutant['mutationType'] # Grab the specific type of mutation it is
        aminos = re.findall(r'([A-Z][a-z][a-z])|\*', mutType)

        if len(aminos) == 0: #! There's no amino acids nor stop codons predicted - we really don't care about it, as we can't tell whether it's synonymous or non-synonymous, so we just get rid of it...
            aminos = ['unk', 'unk']
        if len(aminos) == 1: #! If this is only one, then there's no second amino acid, nor stop codon (i.e. there's no predicted change). Usually, this means that it's unknown, '?'.
            if '?' in mutType:
                aminos.append('unk')
        for i in range(2): # If there's more than two, it'll be because of extension (HGVS) notation or something, which we don't particularly care about
            if aminos[i] in ['', '*']: # If it's a stop codon (denoted with a *, because sometime's they're 'ter' instead. Also, re seems to not read the star because it's a special character, so outputs '' to the list)
                aminos[i] = 'ter'
            elif aminos[i] == '?': # A question mark means that the outcome is unknown (unk)
                aminos[i] = 'unk' 
            aminos[i] = AA_THREE_TO_ONE[aminos[i].lower()]
        wtaa, mutaa = aminos[0], aminos[1] #? re.findall outputs in order, so the wildtype will always come before the mutant in the list

        mutType = f'{wtaa} -> {mutaa}'
        if (wtaa == '?') and (mutaa == '?'): # If both elements are unknown, we just scrap the SNP
            pass
        else:
            if wtaa == mutaa: # If they're the same, then the SNP is synonymous (TRUE). Else, it's non-synonymous (FALSE)
                synon = True
            else:
                synon = False
    
            previousLength = len(mutantsCheck) # Checks to see if a new (unique) SNP has been added to the set of unique SNPs
            mutantsCheck.add((geneName, mutant['locus'], synon))
            if len(mutantsCheck) != previousLength:
                uniqueMutants[f"{geneName}/{locus}"] = {'chromosome':chromosome, 'synon':synon} # If the SNP IS unique, then add it to the dict of unique SNPs
        
            if not silent:
                try:
                    bar, tick = progressBar((mutantNum // linesPerPercent), tick)
                except ZeroDivisionError:
                    bar, tick = progressBar((mutantNum), tick)
                print(f'Analysing mutants: {bar}', end='')
    
    mutsPerPercent = len(uniqueMutants) // 100
    for mutNum, mutant in enumerate(uniqueMutants): # For each gene, calculate a pN and pS count
        geneName = mutant.split('/')[0]
        currentRatio = [0, 0] # Defaults to the gene being new (i.e. no synonymous or non-synonmous variants have yet been found for it)
        if geneName in cumulativeMutations.keys(): # If the gene has already been recorded (i.e. it's been encountered already), update the ratio to whatever is currently on file
            currentRatio = cumulativeMutations[geneName]['count']
        if uniqueMutants[mutant]['synon']:
            currentRatio[1] += 1 # The ratio is in the form [non-synonymous, synonymous]
        else:
            currentRatio[0] += 1
        cumulativeMutations[geneName] = {'count':currentRatio, 'chromosome':uniqueMutants[mutant]['chromosome']} # Either way, (re-)assign the ratio back to its dictionary entry

        if not silent:
            try:
                bar, tick = progressBar((mutNum // mutsPerPercent), tick)
            except ZeroDivisionError:
                bar, tick = progressBar(mutNum, tick)
            print(f'Calculating ratios: {bar}', end='')
    
    if not silent:
        print(f" Done! {' ' * 120}")

    return cumulativeMutations # Return all the unique variants in that file


silent = False

def removeUncounted(mutations): #* Removes all SNPs that have 0 for either the synonymous or non-synonymous count (as these cannot produce a ratio)
    viableMutations = {}
    nonViableMutations = {}
    for mutant in mutations:
        if mutations[mutant]['count'][0] != 0 and mutations[mutant]['count'][1] != 0:
            viableMutations[mutant] = mutations[mutant]
            viableMutations[mutant]['ratio'] = mutations[mutant]['count'][0] / mutations[mutant]['count'][1]
        else:
            nonViableMutations[mutant] = mutations[mutant]
            nonViableMutations[mutant]['ratio'] = 'ERR'
    
    return viableMutations, nonViableMutations

cumulativeMutations = {}

mutants = []
with open(f'{PATH}/Collated CaeNDR Mutants.json', 'r') as FoI: #* Read in the list of all SNPs from the CaeNDR database
    for line in FoI:
        mutants.append(json.loads(line.strip()))

cumulativeMutations = getRatios(mutants, cumulativeMutations) # Calculate the pN/pS ratios

viableMutations, nonViableMutations = removeUncounted(cumulativeMutations) # Remove all the genes for which a pN/pS ratio cannot be calculated
if not silent:
    print(f'{len(viableMutations)} different unique genes found. Writing to output file (/Resources/dNdS_Ratios)')

def saveToJSON(mutations, fileName): #* Save the counts to a JSON file
    with open(f'{PATH}/{fileName}', 'w') as FoI:
        FoI.write(json.dumps(mutations))

def saveToCSV(mutations): #* Save the counts to a csv file
    with open(f'{PATH}/Resources/NEW dNdS_Ratios.csv', 'w') as FoI:
        FoI.write("Gene Name, Non-Synonymous Mutations, Synonymous Mutations, dN/dS Ratio, Locus" + '\n')
        for mutation in mutations:
            FoI.write(f"{mutation}, {mutations[mutation]['count'][0]}, {mutations[mutation]['count'][1]}, {mutations[mutation]['ratio']}, {mutations[mutation]['chromosome']}" + '\n')

#- saveToCSV(viableMutations)
saveToJSON(viableMutations, 'dNdS_Ratios v2.json')
