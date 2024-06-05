# This takes the count .json files produced from setting rate to FALSE in Mutation_Rate_By_Base.py, and formats them into a human-readable tab-separated .txt file

import json
import os

PATH = f'{os.path.dirname(__file__)}'.replace('\\', '/') # The path to the working folder of the project

INPUT_FILE = f'{PATH}/OLD Individual CaeNDR Nucleotide Mutation Counts (+1).json'
OUTPUT_FILE = f'{PATH}/OLD Individual AG Nucleotide Mutation Counts.txt'

with open(INPUT_FILE, 'r') as FoI:
    genes = json.loads(FoI.read())

with open(OUTPUT_FILE, 'w') as FoI:
    FoI.write('Gene ID\tTotal Mutations\tdN/dS\tpi\tpN/pS\tLength\tA Mutations\tA Length\tC Mutations\tC Length\tG Mutations\tG Length\tT Mutations\tT Length\n')
    for gene in genes:
        data = genes[gene]
        aMuts, aLength = data['A']
        cMuts, cLength = data['C']
        gMuts, gLength = data['G']
        tMuts, tLength = data['T']
        totalLength = aLength + cLength + gLength + tLength
        totalMuts = aMuts + cMuts + gMuts + tMuts
        FoI.write(f"{gene}\t{totalMuts}\t{data['dRatio'] if data['dRatio'] != None else '-'}\t{data['pi'] if data['pi'] != None else '-'}\t{data['pRatio'] if data['pRatio'] != None else '-'}\t{totalLength}\t{aMuts}\t{aLength}\t{cMuts}\t{cLength}\t{gMuts}\t{gLength}\t{tMuts}\t{tLength}\n")
