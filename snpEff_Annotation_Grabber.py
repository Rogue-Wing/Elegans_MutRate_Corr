# Takes an annotated (by snpEff) .vcf file as an input, and pulls out all the synonymous and non-synonymous exon mutants from the file
#? Unlike the Annotated_Analyser, the Grabber doesn't care about the type of mutation - that is meant to be dealt with locally later on. This was done to decrease Genoa processing times, and because it was difficult to account for all mutation types without being able to see all the data

import json
import os
import sys

PATH = f'{os.path.dirname(__file__)}'.replace('\\', '/') # The path to the working folder of the project

NUMERAL_TO_NUM = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10, 'XI':11, 'XII':12, 'MtDNA':13}
AA_THREE_TO_ONE = {'ala':'A', 'arg':'R', 'asn':'N', 'asp':'D', 'cys':'C', 'glu':'E', 'gln':'Q', 'gly':'G', 'his':'H', 'ile':'I', 
                  'leu':'L', 'lys':'K', 'met':'M', 'phe':'F', 'pro':'P', 'ser':'S', 'thr':'T', 'trp':'W', 'tyr':'Y', 'val':'V'}

def parseData(line, oldAlleles): #* Takes in a line from the .vcf file, and parses it to return a standardised mutation/allele JSON-ised structure of any SNPs
        line = line.strip().split('\t') # Split the line by tabs. This produces eight elements (CHROM, POS, ID, REF, ALT, QUAL, FILTER/ANN, INFO)
        locus = f'{line[0]}:{line[1]}' #? Gets the chromosome number and position (within the chromosome) - i.e. the locus
        if not (len(line[3]) == 1 and len(line[4]) == 1): #! If either nucleotide string is longer (or shorter) than one, it's not an substitution SNP  - scrap it
            return None, oldAlleles
        substitution = f'{line[3]} -> {line[4]}' #? Gets the SNP
        filterData = line[7].split(';') # Split the FILTER/ANN element into filter parameters
        for filter in filterData:
            if filter[0:4] == 'ANN=': # If the filter element starts with 'ANN=', then it's the annotations!
                annotationData = filter.split(',') # Grab the annotationData and split it into individual annotations
        alleleData = []
        for annotation in annotationData:
            annotation = annotation.split('|')
            if annotation[10] != '': #? If it's protein coding
                #- wtaa, mutaa = AA_THREE_TO_ONE[annotation[10][2:5].lower()] if annotation[10][2] != '*' else '*', AA_THREE_TO_ONE[annotation[10][-3:].lower()] if annotation[10][-1] != '*' else '*' # Converts the amino acids into their single characteer representations
                alleleData.append([annotation[4], annotation[10]]) # Append the WBID and the amino acid change

        newAlleles = []
        for allele in alleleData: # For all the alleles that have been collected from this line
            if allele not in oldAlleles: # If the allele is new (i.e. a different locus and/or different mutation) 
                oldAlleles.append(allele) # Add the allele to the oldAlleles for future checking 
                newAlleles.append({'locus':locus, 'gene':allele[0], 'substitution':substitution, 'mutationType':allele[1], 'orderingData':(NUMERAL_TO_NUM[line[0]], int(line[1]))}) # Add the information to the values list
        
        return newAlleles, oldAlleles

def getSNPs(strain:str, silent:bool = False): #* Get the SNPs from an annotated .vcf file
    if not os.path.exists(f'{PATH}/Annotated_VCFs'): # If the input folder doesn't exist, make it. Note that if this mkdir() has to occur, the program will fail (elegantly)
        os.mkdir(f'{PATH}/Annotated_VCFs')
    
    if not os.path.exists(f'{PATH}/Summarised_Mutants'): # If the output folder doesn't exist, make it
        os.mkdir(f'{PATH}/Summarised_Mutants')
    
    fileName = f'{PATH}/Annotated_VCFs/{strain}.annotated.vcf' # Open the file of the specified strain
    try:
        with open(fileName, 'r') as FoI: #? Opens the file and gets the number of lines (for the purposes of producing a progress bar)
            totalLineCount = 0 # The total number of lines (including commented ones)
            dataLineCount = 0 # The number of lines that contain variant data
            for line in FoI:
                totalLineCount += 1
                if line[0] != '#':
                    dataLineCount += 1
    except FileNotFoundError:
        print(f"ERR: File '{strain}.annotated.vcf' does not exist")
        return FileNotFoundError
    
    if not silent: #? If output to the console is desired
        print(f'{strain}.annotated.vcf consists of {dataLineCount} line(s) of data')
        linesPerPercent = totalLineCount // 100 # Calculate how many lines read constitutes 1% of the file
        print('Progress: ', end='')

    numberOfAlleles = 0
    if os.path.isfile(f'{PATH}/Summarised_Mutants/{strain} Exon Mutants.json'):
        os.remove(f'{PATH}/Summarised_Mutants/{strain} Exon Mutants.json')
    with open(fileName, 'r') as inputFoI:
        addedAlleleData = [] #? Holds a easily-searched (tuple) list of the unique data for all the alleles in the file
        for lineNum, line in enumerate(inputFoI):
            if line[0] != '#': # If the line is not a comment (i.e. it's a line of *actual* data)
                newAlleles, addedAlleleData = parseData(line, addedAlleleData) #* Parse the data
                if newAlleles: # If any new alleles were found
                    with open(f'{PATH}/Summarised_Mutants/{strain} Exon Mutants.json', 'a') as outputFoI: # Write the alleles (JSON-ified) to the output file, labelled as coming from this strain
                        for allele in newAlleles:
                            numberOfAlleles += 1
                            outputFoI.write(json.dumps(allele) + '\n') #! Note this writes them immediately, rather than storing them in a list first, in the hope to optimise memory
            
            if not silent:
                if lineNum % (linesPerPercent * 2) == 0:
                    print('-', end='') # Print a progress bar

    if not silent: # Output some useful information (if desired) once the data has been read
        print(' Done!')
        print(f'{numberOfAlleles} exon variants identified and written to output')

if __name__ == '__main__':
    getSNPs(sys.argv[1])