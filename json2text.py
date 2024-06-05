#

import json
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project
FILE_NAME = f'{PATH}/MA Relative Coding Promoter Offsets (Nearest, 10000, Labelled)'

jsonData = []
with open(f'{FILE_NAME}.json', 'r') as FoI:
    for line in FoI:
        jsonData.append(json.loads(line.strip()))

colHeads = list(jsonData[0].keys())
with open(f'{FILE_NAME}.txt', 'w') as FoI:
    FoI.write('\t'.join(colHeads) + '\n')
    for line in jsonData:
        FoI.write('\t'.join([str(line[element]) for element in line]) + '\n')
