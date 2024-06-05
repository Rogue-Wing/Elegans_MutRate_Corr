# Takes the C. elegans genome sequence, gets the start and end positions of all the genes (along with their WBID) and returns the sequence of that gene

import json
import math
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project
POSITIONS_FILE = f'{PATH}/Ce11 Protein Coding Genes.txt'
GENES_FILE = f'{PATH}/Combined Selection Measures v3.json'
GENOME_FILE = f'{PATH}/WBcel235 Genome.json'

OUTPUT_FILE = f'{PATH}/WS235 TEST.fa'

NUMERAL_TO_CHROM = {'I':'chr1', 'II':'chr2', 'III':'chr3', 'IV':'chr4', 'V':'chr5', 'X':'chr10', 'MtDNA':'MtDNA'} # Converts numeral chromosome names to shorthand variants
COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a'} # Converts nucleotides owing to the complementary base pairing rule

def progressBar(currentPercent, tick): #* Outputs a string for a progress bar
    symbol = {0:'|', 1:'/', 2:'-', 3:'\\'}[tick]

    tick += 1
    if tick > 3:
        tick = 0

    return f"{'-' * (currentPercent)}{' ' * (100 - currentPercent)} [{currentPercent}%] {symbol}\r", tick

def getPositions(): #* Gets the (start and end) positions of the genes, as well as their lengths
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

def getGenome(): #* Reads in the N2 genome
    with open(GENOME_FILE, 'r') as FoI:
        return json.loads(FoI.read())

#* Bins by the pN/pS ratio into specifically assigned intervals
def binBypNpS(genes):
    for gene in genes:
        if genes[gene]['pRatio'] < 0.05:
            genes[gene]['bin'] = 'Extremely Purifying'
        elif 0.05 <= genes[gene]['pRatio'] < 0.5:
            genes[gene]['bin'] = 'Strongly Purifying'
        elif 0.5 <= genes[gene]['pRatio'] < 0.8:
            genes[gene]['bin'] = 'Weakly Purifying'
        elif 0.8 <= genes[gene]['pRatio'] < 1.2:
            genes[gene]['bin'] = 'Neutral'
        elif 1.2 <= genes[gene]['pRatio'] < 2.5:
            genes[gene]['bin'] = 'Diversifying'
        else:
            genes[gene]['bin'] = 'Extemely Diversifying'

    return genes

#* Bins by the pi value into specifically assigned intervals
def binBypi(genes):
    for gene in genes:
        genes[gene]['bin'] = None
        if genes[gene]['pi'] != None:
            if genes[gene]['pi'] > 0:
                if math.log10(genes[gene]['pi']) < -3.75:
                    genes[gene]['bin'] = 'Extremely Purifying'
                elif -3.75 <= math.log10(genes[gene]['pi']) < -3.25:
                    genes[gene]['bin'] = 'Strongly Purifying'
                elif -3.25 <= math.log10(genes[gene]['pi']) < -2.75:
                    genes[gene]['bin'] = 'Purifying'
                else:
                    genes[gene]['bin'] = 'Weakly Purifying'

    return genes

#* Bins by the dN/dS ratio into specifically assigned intervals
def binBydNdS(genes):
    for gene in genes:
        genes[gene]['bin'] = None
        if genes[gene]['dRatio'] != None:
            if genes[gene]['dRatio'] < 0.03:
                genes[gene]['bin'] = 'Extremely Purifying'
            elif 0.03 <= genes[gene]['dRatio'] < 0.08:
                genes[gene]['bin'] = 'Strongly Purifying'
            elif 0.08 <= genes[gene]['dRatio'] < 0.15:
                genes[gene]['bin'] = 'Purifying'
            elif 0.15 <= genes[gene]['dRatio'] < 0.25:
                genes[gene]['bin'] = 'Weakly Purifying'
            else:
                genes[gene]['bin'] = 'Diversifying'

    return genes

def binGenes(genes, selectionMeasure, bin): #* A generic binning frame, which calls to the specific functions above
    binnedGenes = {}
    if selectionMeasure == 'pN/pS':
        genes = binBypNpS(genes)
        selectionMeasure = 'pRatio'
    elif selectionMeasure == 'pi':
        genes = binBypi(genes)
        selectionMeasure = 'pi'
    elif selectionMeasure == 'dN/dS':
        genes = binBydNdS(genes)
        selectionMeasure = 'dRatio'
    else:
        print('ERR: Unknown selection measure')
        return ValueError

    for gene in genes:
        if genes[gene]['bin'] == bin:
            binnedGenes[gene] = genes[gene]

    return binnedGenes

def assignSequenceToGenes(genome, positions, geneNames): #* Assigns the sequence (from the genome file) to all provided genes by reading in the stretch of nucleotides between the genes' start and end positions
    bar, tick = progressBar(0, 0)
    print(f'{bar}', end = '')
    genesPerPercent = len(geneNames) / 100

    errs = 0
    seqGenes = {}
    for geneNum, geneName in enumerate(geneNames):
        if '-' not in geneName: #! If it's not an intergenic region (genes need binning, so we can't use intergenic regions)
            try:
                geneData = positions[geneName] # Grab the data on the start, end, length, strand and chromosome for that gene
                seq = genome[geneData['chromosome']][geneData['start'] - 1:geneData['end']]
                if geneData['strand'] == '-1': #! If the strand is negative, then the complement of the genome (which is of the forward strand) needs to be taken
                    compSeq = ''
                    for nuc in seq:
                        compSeq += COMPLEMENT[nuc]
                    seq = compSeq
                seqGenes[geneName] = seq
            except KeyError: # If the WBID cannot be found in the set of coding genes, then either it's not a coding gene, or the name used by the MA data is out-of-date
                errs += 1
        
        bar, tick = progressBar(int((geneNum + 1) / genesPerPercent), tick)
        print(f'{bar}', end = '')
    print(f'Done! {" " * 95}')
    
    return seqGenes, errs

def writeBinnedSNPs(seqGenes, binName, selectionMeasure, dataset): #* Outputs the binned sequences to a FASTA file
    with open(f'{PATH}/WS235 {binName} ({selectionMeasure}) Transcripts ({dataset}).fa', 'w') as FoI:
        bar, tick = progressBar(0, 0)
        print(f'{bar}', end = '')
        genesPerPercent = len(seqGenes) / 100

        for geneNum, gene in enumerate(seqGenes):
            FoI.write(f'>{gene}\n') # Write the header of the FASTA file
            for nucNum, nuc in enumerate(seqGenes[gene]): # Write the sequence of the gene in 50bp lines
                FoI.write(nuc)
                if ((nucNum + 1) % 50 == 0) and (nucNum > 0):
                    FoI.write('\n')
            FoI.write('\n')

            bar, tick = progressBar(int((geneNum + 1) / genesPerPercent), tick)
            print(f'{bar}', end = '')
        print(f'Done! {" " * 95}')


with open(GENES_FILE, 'r') as FoI: #* Reads in the genes of the AG dataset
    allGenes = json.loads(FoI.read())

genome = getGenome()
genePositions = getPositions()

maGenes = {}
caendrGenes = {}
for gene in allGenes:
    if allGenes[gene]['ma']: # Separates out the MA mutants into another separate list
        maGenes[gene] = allGenes[gene]
    caendrGenes[gene] = allGenes[gene] #? The AG dataset should include *all* genes
print(f"Number of MA genes - {len(maGenes)}, number of AG genes - {len(caendrGenes)}")


# Iterates through the different bins (low and high for pN/pS, pi and dN/dS) for the two datasets, saving each to its own FASTA file
runs = [('CaeNDR', 'pN/pS', 'Strongly Purifying'), ('CaeNDR', 'pN/pS', 'Diversifying'), ('CaeNDR', 'pi', 'Extremely Purifying'), ('CaeNDR', 'pi', 'Weakly Purifying'), 
        ('CaeNDR', 'dN/dS', 'Extremely Purifying'), ('CaeNDR', 'dN/dS', 'Weakly Purifying'), ('MA', 'pN/pS', 'Strongly Purifying'), ('MA', 'pN/pS', 'Diversifying'), 
        ('MA', 'pi', 'Extremely Purifying'), ('MA', 'pi', 'Weakly Purifying'), ('MA', 'dN/dS', 'Extremely Purifying'), ('MA', 'dN/dS', 'Weakly Purifying')]
for run in runs:
    if run[0] == 'CaeNDR':
        genes = caendrGenes
    elif run[0] == 'MA':
        genes = maGenes
    binnedGenes = binGenes(genes, run[1], run[2])
    seqGenes, errs = assignSequenceToGenes(genome, genePositions, list(binnedGenes.keys())) 

    print(f'There are {len(seqGenes)} gene(s). Sequence data could not be found for {errs} gene(s).')
    writeBinnedSNPs(seqGenes, run[2], run[1].replace('/', ''), run[0])
