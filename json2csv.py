# Takes a JSON file (either a dict, a dict of dicts or a list of dicts), and converts it into a csv file

import os
import json

#- PATH = os.path.dirname(__file__) # Gets the path of the working directory (locally, this is [...]/Part II Project/Coding)
PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project

FILE_NAME = '/MA Relative High TSS Offsets (Nearest, Combined)'
isDict = False

def saveDictDictCSV(data, colHeads):
    with open(f'{PATH}{FILE_NAME}.csv', 'w') as FoI:
        FoI.write(f",{','.join(colHeads)}" + '\n')
        for entry in data:
            FoI.write(f"{entry},{','.join([(str(data[entry][key]) if data[entry][key] != None else '') for key in data[entry].keys()])}" + '\n') 

def saveListDictCSV(data, colHeads):
    with open(f'{PATH}{FILE_NAME}.csv', 'w') as FoI:
        FoI.write(f"{','.join(colHeads)}" + '\n')
        for entry in data:
            FoI.write(f"{','.join([(str(entry[key]) if entry[key] != None else '') for key in entry.keys()])}" + '\n')

with open(f'{PATH}{FILE_NAME}.json', 'r') as FoI:
    if isDict:
        data = FoI.read()
        jsonData = json.loads(data)
    else:
        jsonData = []
        for line in FoI:
            if line[0] != '#':
                jsonData.append(json.loads(line.strip()))
        
if isinstance(jsonData, dict):
    if isinstance(jsonData[list(jsonData.keys())[0]], dict): # If true, then it's a dict of dicts. Hence, the outer keys are our 'row headers', and the inner keys are the column headers
        colHeads = jsonData[list(jsonData.keys())[0]].keys() #! Note this assumes that all the entries (the inner dicts) have the same keys (as the first element)
        saveDictDictCSV(jsonData, colHeads)
    else: # Else it's just a dict
        pass

elif isinstance(jsonData, list):
    if isinstance(jsonData[0], dict): # If true, then it's a list of dicts
        colHeads = list(jsonData[0].keys())
        saveListDictCSV(jsonData, colHeads)
else:
    print('ERR: Currently unconvertable format')
    exit