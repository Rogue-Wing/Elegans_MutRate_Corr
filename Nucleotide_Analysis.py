# Produces summarised nucleotide substitution data of the six different substitutions (UPDATED 18/04/2024 to incorporate the dN/dS and McKt stats)

import json
import re

import time
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project

COMBINED_MUTANTS_FILE = f'{PATH}/Combined Selection Measures v3.json' # Get the file that has all the selection measures in it
CaeNDR_MUTANTS_FILE = f'{PATH}/Collated CaeNDR Mutants v2.json' # Gets the file of all CaeNDR substitutions
MA_MUTANTS_FILE = f'{PATH}/Combined MA Data v3.json' # Gets the file of all MA substitutions

AA_THREE_TO_ONE = {'ala':'A', 'arg':'R', 'asn':'N', 'asp':'D', 'cys':'C', 'glu':'E', 'gln':'Q', 'gly':'G', 'his':'H', 'ile':'I', 
                  'leu':'L', 'lys':'K', 'met':'M', 'phe':'F', 'pro':'P', 'ser':'S', 'thr':'T', 'trp':'W', 'tyr':'Y', 'val':'V', 'ter':'*', 'unk':'?'}

NUC_CONV = {'A':'T', 'C':'G', 'G':'C', 'T':'A'} # Convert between complementary base pairs

def progressBar(currentPercent, tick): # Creates a string to represent a progress bar
    symbol = {0:'|', 1:'/', 2:'-', 3:'\\'}[tick]

    tick += 1
    if tick > 3:
        tick = 0

    return f"{'-' * (currentPercent)}{' ' * (100 - currentPercent)} [{currentPercent}%] {symbol}\r", tick

def getMutType(mutationType): # Gets the type of mutation/SNP (i.e. whether it falls within a protein-coding region or not)
        aminos = re.findall(r'([A-Z][a-z][a-z])|\*', mutationType) # Grabs the amino acids

        if len(aminos) == 0: #! There's no amino acids nor stop codons predicted - we really don't care about it, as we can't tell whether it's synonymous or non-synonymous, so we just get rid of it...
            aminos = ['unk', 'unk']
        if len(aminos) == 1: #! If this is only one, then there's no second amino acid, nor stop codon (i.e. there's no predicted change). Usually, this means that it's unknown, '?'.
            if '?' in mutationType:
                aminos.append('unk')
        for i in range(2): # If there's more than two, it'll be because of extension (HGVS) notation or something, which we don't particularly care about
            if aminos[i] in ['', '*']: # If it's a stop codon (denoted with a *, because sometime's they're 'ter' instead. Also, re seems to not read the star because it's a special character, so outputs '' to the list)
                aminos[i] = 'ter'
            elif aminos[i] == '?': # A question mark means that the outcome is unknown (unk)
                aminos[i] = 'unk' 
            aminos[i] = AA_THREE_TO_ONE[aminos[i].lower()]
        wtaa, mutaa = aminos[0], aminos[1] #? re.findall outputs in order, so the wildtype will always come before the mutant in the list

        if wtaa == mutaa: # If the mutant and reference sequences are the same, then there's no change - it's synonymous
            return 'Exon/Synonymous'
        elif wtaa == '*': #? We already know that it's not synoynmous, so this must be stop-loss
            return 'Exon/Stop_Loss'
        elif mutaa == '*': #? This must be nonsense
            return 'Exon/Nonsense'
        elif wtaa == 'M': #! It should be noted that this technically isn't correct, as it considers *all* Ms, not just the first one in the sequence. For the sake of the analysis performed, however, this does not matter. Note that this should be fixed if anyone ever wants to consider specific types of mutations, however!!
            return 'Exon/Start_Loss'
        else:
            return 'Exon/Missense'

def writeAllSNPs(): # Collates all the SNPs and mutants together, marks whether they come from the MA (mutant) or CaeNDR (SNP) data, and writes them to a file
    with open(COMBINED_MUTANTS_FILE, 'r') as FoI: # Reads in the selection measures
        combinedMutants = json.loads(FoI.read())

    caeNDRMutants = []
    with open(CaeNDR_MUTANTS_FILE, 'r') as FoI: # Reads in the SNPs
        for line in FoI:
            if line[0] != '#':
                caeNDRMutants.append(json.loads(line.strip()))

    MAMutants = []
    with open(MA_MUTANTS_FILE, 'r') as FoI: # Reads in the mutants
        for line in FoI:
            if line[0] != '#': # Ignores the line if it's #'d out
                MAMutants.append(json.loads(line.strip()))

    allMutants = {}

    mutsPerPercent = len(combinedMutants) // 100
    print(len(combinedMutants))
    bar, tick = progressBar(0, 0)
    print(bar, end = '')

    for mutNum, mutant in enumerate(combinedMutants): # For every gene in the combined selection measures (i.e. the AG dataset)
        length = combinedMutants[mutant]['length'] # Get all the background data for that gene
        hits = combinedMutants[mutant]['numberOfHits']
        pRatio = combinedMutants[mutant]['pRatio']
        tajD = combinedMutants[mutant]['tajD']
        pi = combinedMutants[mutant]['pi']
        dRatio = combinedMutants[mutant]['dRatio']
        McKt = combinedMutants[mutant]['McKt']
        allMutants[mutant] = {'length':length, 'hits':hits, 'ratio':pRatio, 'tajD':tajD, 'pi':pi, 'dRatio':dRatio, 'McKt':McKt, 'snps':[]}
        for entry in caeNDRMutants: # If the gene is found within the list of SNPs, add the SNP data
            if entry['gene'] == mutant:
                geneData = {'locus':entry['locus'], 'substitution':entry['substitution'], 'mutationType':getMutType(entry['mutationType']), 'numStrains':entry['strainNum'], 'ma':False}
                allMutants[mutant]['snps'].append(geneData)
        for entry in MAMutants: # If the gene is found within the list of mutants, add the mutant data
            if entry['gene'] == mutant:
                geneData = {'locus':entry['locus'], 'substitution':entry['substitution'], 'mutationType':entry['mutationType'],'ma':True}
                allMutants[mutant]['snps'].append(geneData)
        
        bar, tick = progressBar(mutNum // mutsPerPercent, tick)
        print(bar, end = '')

    with open(f'{PATH}/All SNPs v3.json', 'w') as FoI: # Write all SNPs/mutants to a single file
        FoI.write(json.dumps(allMutants))

def writeAllSNPWithTypes(): #! UNUSED
    subs = []
    with open(f'{PATH}/All SNPs v3.json', 'r') as FoI:
        genes = json.loads(FoI.read())

    mutTypeSNPs = []
    for geneNum, gene in enumerate(genes):
        for SNP in genes[gene]['snps']:
            mutType = SNP['mutationType']
            if mutType.split('/')[0] == 'Exon':
                if mutType.split('/')[1] != 'Synonymous':
                    mutType = 'Exon/Non-synonymous'

            sub = f"{SNP['substitution'][0]}>{SNP['substitution'][-1]}"
            compSub = f"{NUC_CONV[SNP['substitution'][0]]}>{NUC_CONV[SNP['substitution'][-1]]}"
            if compSub in subs:
                substitution = compSub
            else:
                substitution = sub
                if sub not in subs:
                    subs.append(sub)
            mutTypeSNPs.append({'geneName':gene, 'ratio':genes[gene]['ratio'], 'tajD':genes[gene]['tajD'], 'piNorm':genes[gene]['piNorm'], 'substitution':substitution, 'mutationType':mutType, 'ma':SNP['ma']})

    with open(f'{PATH}/Simplified SNPs v3.json', 'w') as FoI:
        for SNP in mutTypeSNPs:
            FoI.write(json.dumps(SNP) + '\n')

def summariseSNPs(): # Takes the file of all SNPs/mutants and merges them back down into to one entry per gene
    with open(f'{PATH}/All SNPs v3.json', 'r') as FoI:
        genes = json.loads(FoI.read())

    summarisedGenes = {}
    for geneNum, gene in enumerate(genes): # For each SNP/mutant
        pRatio = float(genes[gene]['ratio']) if genes[gene]['ratio'] != None else None # Get the measures of selection, if they exist
        tajD = float(genes[gene]['tajD']) if genes[gene]['tajD'] != None else None 
        pi = float(genes[gene]['pi']) if genes[gene]['pi'] != None else None 
        dRatio = float(genes[gene]['dRatio']) if genes[gene]['dRatio'] != None else None 
        McKt = float(genes[gene]['McKt']) if genes[gene]['McKt'] != None else None 
        geneData = {'length':genes[gene]['length'], 'hits':genes[gene]['hits'], 'pRatio':pRatio, 'tajD':tajD, 'pi':pi, 'dRatio':dRatio, 'McKt':McKt, 'NSSubs':{}, 'MASubs':{}} # Write in the selection measures and the background information
        for snp in genes[gene]['snps']: # For each SNP/mutant, get the subsitition it correlates to (in the form WT_Nuc>Mut_Muc), and flip it by complementary base-pairing rules if it's not one of the listed six
            sub = f"{snp['substitution'][0]}>{snp['substitution'][-1]}"
            compSub = f"{NUC_CONV[snp['substitution'][0]]}>{NUC_CONV[snp['substitution'][-1]]}"
            if sub not in ['A>G', 'G>C','C>T' ,'C>A' ,'T>A' ,'A>C']:
                sub = compSub

            if snp['ma']: # If it's a mutant
                if sub in geneData['MASubs']: # Add it to the list of mutants for that gene
                    geneData['MASubs'][sub] += 1
                else:
                    geneData['MASubs'][sub] = 1
            else: # It's a SNP
                if sub in geneData['NSSubs']: # Add it to the list of SNPs for that gene
                    geneData['NSSubs'][sub] += 1
                else:
                    geneData['NSSubs'][sub] = 1
        if geneData['NSSubs'] == {}:
            geneData['NSSubs'] = None
        if geneData['MASubs'] == {}:
            geneData['MASubs'] = None
        summarisedGenes[gene] = geneData
    
    return summarisedGenes

def saveToJSON(summarisedSNPs, fileName = 'Summarised SNPs v3'): # Save all the data to a combined JSON file
    with open(f'{PATH}/{fileName}.json', 'w') as FoI:
        FoI.write(json.dumps(summarisedSNPs))

def saveNSToCSV(summarisedSNPs): # Save the CaeNDR SNPs to a .csv file (for R analysis)
    with open(f'{PATH}/Summarised NS SNPs v3.csv', 'w') as FoI:
        FoI.write('geneName,length,hits,pRatio,tajD,piNorm,dRatio,McKt,A>G,G>C,C>T,C>A,T>A,A>C\n')
        for gene in summarisedSNPs:
            if summarisedSNPs[gene]['NSSubs']: # If there's at least one substitution denoted in NSSubs (there should always be)
                subs = f"{summarisedSNPs[gene]['NSSubs']['A>G'] if 'A>G' in summarisedSNPs[gene]['NSSubs'] else ''},{summarisedSNPs[gene]['NSSubs']['G>C'] if 'G>C' in summarisedSNPs[gene]['NSSubs'] else ''},{summarisedSNPs[gene]['NSSubs']['C>T'] if 'C>T' in summarisedSNPs[gene]['NSSubs'] else ''},{summarisedSNPs[gene]['NSSubs']['C>A'] if 'C>A' in summarisedSNPs[gene]['NSSubs'] else ''},{summarisedSNPs[gene]['NSSubs']['T>A'] if 'T>A' in summarisedSNPs[gene]['NSSubs'] else ''},{summarisedSNPs[gene]['NSSubs']['A>C'] if 'A>C' in summarisedSNPs[gene]['NSSubs'] else ''}" 
            else:
                subs = f"{''},{''},{''},{''},{''},{''}"
            FoI.write(f"{gene},{summarisedSNPs[gene]['length']},{summarisedSNPs[gene]['hits']},{summarisedSNPs[gene]['pRatio']},{summarisedSNPs[gene]['tajD']},{summarisedSNPs[gene]['pi']},{summarisedSNPs[gene]['dRatio']},{summarisedSNPs[gene]['McKt']},{subs}" + '\n')

def saveMAToCSV(summarisedSNPs): # Save the MA mutants to a .csv file (for R analysis)
    with open(f'{PATH}/Summarised MA SNPs v3.csv', 'w') as FoI:
        FoI.write('geneName,length,hits,pRatio,tajD,piNorm,dRatio,McKt,A>G,G>C,C>T,C>A,T>A,A>C\n')
        for gene in summarisedSNPs:
            if summarisedSNPs[gene]['MASubs']: # If there's at least one substitution denoted in MASubs
                subs = f"{summarisedSNPs[gene]['MASubs']['A>G'] if 'A>G' in summarisedSNPs[gene]['MASubs'] else ''},{summarisedSNPs[gene]['MASubs']['G>C'] if 'G>C' in summarisedSNPs[gene]['MASubs'] else ''},{summarisedSNPs[gene]['MASubs']['C>T'] if 'C>T' in summarisedSNPs[gene]['MASubs'] else ''},{summarisedSNPs[gene]['MASubs']['C>A'] if 'C>A' in summarisedSNPs[gene]['MASubs'] else ''},{summarisedSNPs[gene]['MASubs']['T>A'] if 'T>A' in summarisedSNPs[gene]['MASubs'] else ''},{summarisedSNPs[gene]['MASubs']['A>C'] if 'A>C' in summarisedSNPs[gene]['MASubs'] else ''}"
                FoI.write(f"{gene},{summarisedSNPs[gene]['length']},{summarisedSNPs[gene]['hits']},{summarisedSNPs[gene]['pRatio']},{summarisedSNPs[gene]['tajD']},{summarisedSNPs[gene]['pi']},{summarisedSNPs[gene]['dRatio']},{summarisedSNPs[gene]['McKt']},{subs}" + '\n')


writeAllSNPs()
summarisedGenes = summariseSNPs()
saveToJSON(summarisedGenes)
saveNSToCSV(summarisedGenes)
saveMAToCSV(summarisedGenes)