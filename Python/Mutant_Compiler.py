# Takes all the individual 'Summarised_Mutant' .json files of strain-specific SNPs and compiles them into a single .json file

import os
import json

PATH = f'{os.path.dirname(__file__)}'.replace('\\', '/') # The path to the working folder of the project

def progressBar(currentPercent, tick): # Outputs a progress bar to the command line
    symbol = {0:'|', 1:'/', 2:'-', 3:'\\'}[tick]

    tick += 1
    if tick > 3:
        tick = 0

    return f"{'-' * (currentPercent)}{' ' * (100 - currentPercent)} [{currentPercent}%] {symbol}\r", tick

def getMutants(FoI, uniqueMutants): # Gets all the SNPs  in the .json file that are unique compared to every other previously read file
    for line in FoI:
        data = json.loads(line.strip())
        hashKey = f"{data['locus']}{data['substitution']}" # Creates a 'key' of the locus and substitution, as this is the minimum amount of information that determines a unique SNP
        if hashKey not in uniqueMutants:
            data['strainNum'] = 1
            uniqueMutants[hashKey] = data
        else:
            uniqueMutants[hashKey]['strainNum'] += 1
    
    return uniqueMutants

def getFiles(): # Gets all the files in the files in the directory as a list of filenames
    fileNames = []
    for fileName in os.listdir(f'{PATH}/Summarised_Mutants'): # For each item in the path
        if os.path.isfile(f'{PATH}/Summarised_Mutants/{fileName}'): # If the file is actually a file (and not a folder). Note that this assumes all files at the PATH are .json, and are summarised mutant files
            fileNames.append(fileName)
    
    return fileNames

uniqueMutants = {}
fileNames = getFiles()
fileCount = len(fileNames)
filesPerPercent = fileCount / 100

bar, tick = progressBar(0, 0)
print(f"(0/{fileCount}) - {bar}", end = '')

for fileNum, fileName in enumerate(fileNames): # Iterate through each file
    bar, tick = progressBar(int(fileNum // filesPerPercent), tick)
    print(f"({fileNum}/{fileCount}) - {bar}", end = '')
    with open(f'{PATH}/Summarised_Mutants/{fileName}', 'r') as FoI:
        uniqueMutants = getMutants(FoI, uniqueMutants)

with open(f'{PATH}/Resources/CaeNDR Test.json', 'w') as FoI: # Save all the unique SNPs to a single file
    for mutant in uniqueMutants:
        FoI.write(json.dumps(uniqueMutants[mutant]) + '\n')
