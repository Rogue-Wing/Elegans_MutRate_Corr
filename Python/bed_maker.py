# Makes a .bed file from a standard JSON mutant file

import json
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project

CHROM_NAME_CONVERT = {'i':'chr1', 'ii':'chr2', 'iii':'chr3', 'iv':'chr4', 'v':'chr5', 'x':'chr10', 'mtdna':'chrMT'}  # Converts numeral chromosome names to shorthand variants
COMP_CONVERT = {'a':'t', 'c':'g', 'g':'c', 't':'a'}
INDEX_OFFSET = 1

#- DATA = 'Saxena A. S. et al CE Data Formatted.json'
DATA = 'Konrad Unfiltered Data.json'
DATA_NAME = f"{DATA.split(' ')[0]} {DATA.split(' ')[0]}"

mutants = []
with open(f'{PATH}/{DATA}', 'r') as FoI:
    for line in FoI:
        mutants.append(json.loads(line.strip()))

def writeMutantsToBed(mutants): #* Writes all the individual mutants/SNPs to a bed file (TSV txt file)
    with open(f'{PATH}/{DATA_NAME} Mutants [{INDEX_OFFSET}].bed', 'w') as FoI:
        FoI.write('#Chrom' + '\t' + 'chromStart' + '\t' + 'chromEnd' + '\t' + 'name' + '\n') # Writes the header of the bed file
        for mutant in mutants: # For each mutant, write its data as a line in the file
            chrom = mutant['locus'].split(':')[0]
            chromStart = str(int(mutant['locus'].split(':')[1]) - INDEX_OFFSET)
            chromEnd = str(int(chromStart) + 1)
            name = mutant['substitution']
            FoI.write(CHROM_NAME_CONVERT[chrom.lower()] + '\t' + chromStart + '\t' + chromEnd + '\t' + name + '\n')

def testBed(mutants, fileName = None): #* Tests the bed file to check that the nucleotides it has provided for the position data are the same as those that are wt in the mutant data (i.e. it checks whether the data is 0 indexed or not!)
    if not fileName:
        fileName = f'{PATH}/{DATA_NAME} (234) [{INDEX_OFFSET}].out' 
    with open(fileName, 'r') as FoI: # Read in the file
        fullyCorrect, complementaryCorrect, incorrect = 0, 0, 0
        for lineNum, line in enumerate(FoI): # For each line, check if the data in the bed file matches to the .json data
            line = line.strip()
            if line[0] == '>':
                currentGenomePosition = line[1:]
            else:
                chromStart = currentGenomePosition.split(':')[1].split('-')[0]
                if str(int(mutants[lineNum // 2]['locus'].split(':')[1]) - INDEX_OFFSET) == chromStart:
                    if line.lower() == mutants[lineNum // 2]['substitution'][0].lower():
                        fullyCorrect += 1
                    elif COMP_CONVERT[line.lower()] == mutants[lineNum // 2]['substitution'][0].lower():
                        complementaryCorrect += 1 
                    else:
                        incorrect += 1
                        print(f"Problem with {mutants[lineNum // 2]['gene']}. Bed file ({line.lower()}) compared to JSON data {mutants[lineNum // 2]['substitution'][0].lower()}")
                else:
                    print(f"ERR: Incorrect position! {mutants[lineNum // 2]['locus'].split(':')[1]} to {chromStart}")
    
    print(f'Correct: {fullyCorrect}, Complementary: {complementaryCorrect}, Incorrect: {incorrect}')

#- writeMutantsToBed(mutants)
testBed(mutants, f'{PATH}/Konrad Unfiltered [1].out')



