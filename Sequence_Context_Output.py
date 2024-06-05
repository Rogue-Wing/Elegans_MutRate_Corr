# For each mutant/SNP, outputs it and its surrounding nucleotide context (one basepair either side of it) into a txt file

import json
import requests
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project
MA_DATA = f'{PATH}/Combined MA Data.named.json'
CAENDR_DATA = f'{PATH}/Collated CaeNDR Mutants v2.json'
#- RAW_GENOME_DATA = f'{PATH}/WBcel235 Genome.fna'
GENOME_DATA = f'{PATH}/WBcel235 Genome.json'
MINE_DATA = f'{PATH}/WS285 Gene Start-Ends.txt'

def progressBar(currentPercent, tick):
    symbol = {0:'|', 1:'/', 2:'-', 3:'\\'}[tick]

    tick += 1
    if tick > 3:
        tick = 0

    return f"{'-' * (currentPercent)}{' ' * (100 - currentPercent)} [{currentPercent}%] {symbol}\r", tick

'''
lineCount = 0
with open(RAW_GENOME_DATA, 'r') as f:
    for line in f:
        lineCount += 1
linesPerPercent = lineCount // 100

genome = {}
with open(RAW_GENOME_DATA, 'r') as f:
    bar, tick = progressBar(0, 0)
    print(bar, end = '')
    currentChromosome = None
    for lineNum, line in enumerate(f):
        if line[0] != '#':
            if line[0] == '>':
                if currentChromosome:
                    genome[currentChromosome] = nucleotideSequence
                nucleotideSequence = ''
                currentChromosome = line.strip()[1:]
            else:
                nucleotideSequence += line.strip()
    
        bar, tick = progressBar(int(lineNum / linesPerPercent), tick)
        print(bar, end = '')
    
    genome[currentChromosome] = nucleotideSequence # Need to ensure the final chromosome (the final lines) are actually put into the dictionary

print('Done!')

with open(f'{PATH}/WBcel235 Genome.json', 'w') as f:
    f.write(json.dumps(genome))
'''

refGenes = {}
with open(MINE_DATA, 'r') as f:
    for lineNum, line in enumerate(f):
        if lineNum != 0: 
            data = line.strip().split('\t')
            start, end = int(data[4]), int(data[3])
            length = end - start
            refGenes[data[1]] = {'start':start, 'end':end, 'length':length, 'strand':data[2]}

with open(GENOME_DATA, 'r') as f:
    genome = json.loads(f.read())

'''
testGene = requests.get(f"http://api.wormbase.org//rest/widget/gene/{'WBGene00001214'}/location").json()['fields']
print(f"{testGene['genomic_position']}" + '\n')
print(f"{testGene['genomic_image']}" + '\n')
print(f"{testGene['name']}" + '\n')
print(f"{testGene['jbrowse_tracks']}" + '\n')
print(f"{testGene['genetic_position']}" + '\n')
'''

maGenes = []
with open(MA_DATA, 'r') as f:
    for line in f:
        if line[0] != '#':
            maGenes.append(json.loads(line.strip()))
            maGenes[-1]['ma'] = True

caeNDRGenes = []
with open(CAENDR_DATA, 'r') as f:
    for line in f:
        if line[0] != '#':
            caeNDRGenes.append(json.loads(line.strip()))
            caeNDRGenes[-1]['ma'] = False

allGenes = maGenes + caeNDRGenes

errs = 0
mergedGenes = {}
genesPerPercent = len(allGenes) // 100
bar, tick = progressBar(0, 0)
print(bar, end = '')
with open(f'{PATH}/All SNPs Context Data v2.txt', 'w') as f:
    f.write('gene\tsource\tposition\tstart\tend\tstrand\tref\tmut\tcontext\n')
    for geneNum, gene in enumerate(allGenes):
        skip = False
        geneName = gene['gene']
        chromosome = f"chr{gene['orderingData'][0]}"
        if chromosome == 'chr13':
            chromosome = 'chrMT'
        position = int(gene['orderingData'][1])
        refNuc = gene['substitution'][0]
        mutNuc = gene['substitution'][-1]
        if '-' in geneName: #? This means the region is intergenic
            geneNames = geneName.split('-')
            datas = []
            for tempGeneName in geneNames:
                try:
                    datas.append(refGenes[tempGeneName])
                except KeyError:
                    if tempGeneName in mergedGenes:
                        tempGeneName = mergedGenes[tempGeneName]
                    else:
                        mergedGene = requests.get(f"http://api.wormbase.org//rest/field/gene/{tempGeneName}/merged_into").json()['merged_into']['data']
                        if mergedGene:
                            mergedGenes[tempGeneName] = mergedGene['id']
                            tempGeneName = mergedGene['id']
                    try:
                        datas.append(refGenes[tempGeneName])
                    except KeyError:
                        errs += 1
            if len(datas) == 2:
                geneStart, geneEnd = datas[0]['end'], datas[1]['start']
                firstStrand, secondStrand = datas[0]['strand'], datas[1]['strand']
                if firstStrand == secondStrand:
                    strand = firstStrand
                else:
                    strand = '='
            else:
                skip = True

        else:
            try:
                data = refGenes[geneName]
            except KeyError:
                if geneName in mergedGenes:
                    geneName = mergedGenes[geneName]
                else:
                    mergedGene = requests.get(f"http://api.wormbase.org//rest/field/gene/{geneName}/merged_into").json()['merged_into']['data']
                    if mergedGene:
                        mergedGenes[geneName] = mergedGene['id']
                        geneName = mergedGene['id']
                try:
                    data = refGenes[geneName]
                except KeyError:
                    skip = True
                    errs += 1
        
        if not skip:
            geneStart, geneEnd = data['start'], data['end']
            length = data['length']
            strand = data['strand']

        seqContext = f'{genome[chromosome][position - 2]}{genome[chromosome][position - 1]}{genome[chromosome][position]}'.lower() # Index starts from 0 in Python/the genome file, but these data are +1 indexed. Hence, 1 must be removed from the position to translate it from the +1 indexed SNP data to the genome data  
        if seqContext[1].lower() != refNuc.lower():
            print(f'Oh darn - {geneName} (Recorded as {seqContext[1].lower()}, but should be {refNuc.lower()})')


        source = 'MA' if gene['ma'] else 'Polymorphism'

        strand = '+' if strand == '1' else '-' if strand == '-1' else strand

        # 'gene, source, position, start, end, strand, ref, mut, context'
        if skip:
            f.write(f'{geneName}\t{source}\t{chromosome}:{position}\tN/A\tN/A\tN/A\t{refNuc}\t{mutNuc}\t{seqContext}\n')
        else:
            f.write(f'{geneName}\t{source}\t{chromosome}:{position}\t{geneStart}\t{geneEnd}\t{strand}\t{refNuc}\t{mutNuc}\t{seqContext}\n')
        
        bar, tick = progressBar(int(geneNum / genesPerPercent), tick)
        print(bar, end = '')

        #- if geneNum > 50000:
            #- break

print(errs, len(allGenes))

# Several genes are somewhat outdated, for one reason or another, or are intergenic (and thus lack a start/end). These genes are marked with 'N/A' in their start, end and strand fields
# The start and end positions are annotated with WS285, NOT WS235 (i.e. the genome version for the MA/polymorphism data). Hence, the position of the SNP appears to be outside of the start-end range.