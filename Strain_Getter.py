# Combined all the strain names from two different files into one, removing duplicates

import json
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project

S1 = f"{PATH}/Strains_1.txt"
S2 = f"{PATH}/Strains_2.txt"

allStrains = []
strains1 = []
strains2 = []
with open(S1, 'r') as f:
    for line in f:
        strains1.append(line.strip())

allStrains += strains1

with open(S2, 'r') as f:
    for line in f:
        line = line.split('\t')
        if line[1] == 'C. elegans':
            strains2.append(line[0].split(' ')[0].strip())

for strain in strains2:
    duplicate = False
    for confirmedStrain in allStrains:
        if confirmedStrain == strain:
            duplicate = True
    
    if not duplicate:
        allStrains.append(strain)

allStrains.sort()

with open(f'{PATH}/C. elegans NV Strains.json', 'w') as f:
    f.write(json.dumps({'strains':allStrains}))