# Takes a .vcf file and converts it into a .json (it won't carry over annotations)

import json

PATH = 'C:/Users/ginge/Documents/Super Cool Poop/University/Part II Project/Coding/Resources'
#- FILE_NAME = 'Konrad Mutants (Zeroed)'
FILE_NAME = 'Saxena LiftOver.annotated'

CHROM_TO_NUMERAL = {'chr1':'I', 'chr2':'II', 'chr3':'III', 'chr4':'IV', 'chr5':'V', 'chr10':'X', 'chrMT':'MtDNA'}
NUMERAL_TO_NUM = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'X':10, 'MtDNA':13}

vcfLines = []
with open(f'{PATH}/{FILE_NAME}.vcf', 'r') as FoI:
    for line in FoI:
        if line[0] != '#':
            vcfLines.append(line.strip().split('\t'))

with open(f'{PATH}/{FILE_NAME}.json', 'w') as FoI:
    for line in vcfLines:
        locus = f'{line[0]}:{line[1]}'
        gene = None
        substitution = f'{line[3]} -> {line[4]}'
        mutationType = None
        orderingData = [NUMERAL_TO_NUM[line[0]], line[1]]
        FoI.write(json.dumps({"locus":locus, "gene":gene, "substitution":substitution, "mutationType":mutationType, "orderingData":orderingData}) + '\n')


'''
I	92951	.	G	A
{"locus":"Chrom:Pos", "gene":None, "substitution":"X -> X", "mutationType": None, "orderingData": [Chrom, Pos]}
'''
