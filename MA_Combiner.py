# This file takes the named Saxena and Konrad files and merges them into a single combined list of mutations, removing any duplicates in the process

import json
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project

SAXENA_FILE = f'{PATH}/Saxena LiftOver Data.named.json' # The Saxena data (post-LiftOver and with WBIDs assigned)
KONRAD_FILE = f'{PATH}/Konrad Data.named.json' # The Konrad data (with WBIDs assigned)

OUTPUT_FILE = f'{PATH}/Combined MA Data v3.json' # The name for the output file

konradMutants = []
with open(KONRAD_FILE, 'r') as f: # Open the Konrad mutants
    for line in f:
        if line[0] != '#': # Ignore the mutant if it's #'d out
            konradMutants.append(json.loads(line.strip()))

saxenaMutants = []
with open(SAXENA_FILE, 'r') as f: # Open the Saxena mutants
    for line in f:
        if line[0] != '#': # Ignore the mutant if it's #'d out
            saxenaMutants.append(json.loads(line.strip()))

print(len(konradMutants) + len(saxenaMutants))

def reformat(mutant): #* Takes a mutant and ensures all the types for the values in the dict are correct
    locus = str(mutant['locus'])
    gene = str(mutant['gene'])
    substitution = str(mutant['substitution'])
    mutationType = str(mutant['mutationType'])
    orderingData = [int(mutant['orderingData'][0]), int(mutant['orderingData'][1])]
    return {'locus': locus, 'gene': gene, 'substitution': substitution, 'mutationType': mutationType, 'orderingData': orderingData}

maMutants = []
for mutant in konradMutants:
    formattedMutant = reformat(mutant)
    if formattedMutant not in maMutants: # Add all unique mutants
        maMutants.append(formattedMutant)

for mutant in saxenaMutants:
    formattedMutant = reformat(mutant)
    if formattedMutant not in maMutants: # Add all unique mutants
        maMutants.append(formattedMutant)

with open(OUTPUT_FILE, 'w') as f: # Save the file
    for mutant in maMutants:
        f.write(json.dumps(mutant) + '\n')