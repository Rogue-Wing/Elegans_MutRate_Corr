# V2 of the pipeline, with slightly improved annotation grabbing (Annotation_Grabber, rather than Annotated_Analyser). Also has quadrupled the number of bots (the genoa cluster has 48 cores), although this is set to 16 so it doesn't kill the local PC

import threading
import os
import subprocess
import snpEff_Annotation_Grabber as sAG
import sys
from datetime import datetime
import time

import random

PATH = os.path.dirname(__file__) 
MAX_BOTS = 16 # The maximum number of bots. Each both takes a thread on the CPU, so this shouldn't exceed 2x the number of available cores (assuming an Intel processor)

class Bot(): #* The class for bots, which run a process, in this case the analytical pipeline for a single strain of the CaeNDR data

    botList = []

    def __init__(self, sample:str, id:int):
        self.sample = sample
        self.id = id

        self.thread = threading.Thread(target = self.run)
        self.startTime = datetime.now() # Doesn't have to be massively precise, so we'll use datetime

        Bot.botList.append(self)
        self.thread.start()

    def run(self, silent = False): #* The process for the bot
        if not silent:
            timeNow = datetime.now()
            print(f'[{timeNow.day:02}/{timeNow.month:02}/{timeNow.year}] @ {timeNow.hour:02}:{timeNow.minute:02}:{timeNow.second:02} - {self.sample} started split...')
        subprocess.run(f'bcftools view -c1 -Oz -s {self.sample} -o ./Split_VCFs/{self.sample}.vcf.gz ./CaeNDR.220216.HF-Full.vcf.gz'.split(' ')) # Runs a command-line command that runs bcftools to split out a single strain out of the overall vcf file
        if not silent:
            timeNow = datetime.now()
            print(f'[{timeNow.day:02}/{timeNow.month:02}/{timeNow.year}] @ {timeNow.hour:02}:{timeNow.minute:02}:{timeNow.second:02} - {self.sample} split complete! Started annotation...')
        with open(f'./Annotated_VCFs/{self.sample}.annotated.vcf', 'w') as output:
            subprocess.run(f'java -jar ./snpEff/snpEff.jar -noStats WBcel235.75 ./Split_VCFs/{self.sample}.vcf.gz'.split(' '), stdout = output) # Takes that split off .vcf file and annotates it using SnpEff
        os.remove(f'./Split_VCFs/{self.sample}.vcf.gz') # Removes the split .vcf to save space (it's a different file to the annotated version)
        if not silent:
            timeNow = datetime.now()
            print(f'[{timeNow.day:02}/{timeNow.month:02}/{timeNow.year}] @ {timeNow.hour:02}:{timeNow.minute:02}:{timeNow.second:02} - {self.sample} annotation complete! Started grabbing...')
        sAG.getSNPs(self.sample, True) # Uses snpEff_Annotation_Grabber.py to get the SNPs present within protein-coding sequences and saves them to a .json file
        os.remove(f'./Annotated_VCFs/{self.sample}.annotated.vcf') # Removes the annotated .vcf to save space
        if not silent:
            timeNow = datetime.now()
            runTime = timeNow - self.startTime
            print(f'[{timeNow.day:02}/{timeNow.month:02}/{timeNow.year}] @ {timeNow.hour:02}:{timeNow.minute:02}:{timeNow.second:02} - {self.sample} annotation complete! Total run time = {runTime.seconds}s')

        self.kill() # Kills the bot, which triggers the creation of the next one (to deal with the next strain)
    
    def kill(self): #* Kill the bot, ending its processes and removing it from the list of bots
        Bot.botList.remove(self)
        del(self)

def startProcesses(sampleList:list): #! THE MAIN THREAD
    for sample in sampleList:
        while len(Bot.botList) >= MAX_BOTS: # If the bot count is at max, wait. If not, add a new bot. Repeat for every strain in the list of strains
            time.sleep(1)
        Bot(sample, len(Bot.botList))

def program(fileName:str):
    sampleList = []
    with open(fileName, 'r') as FoI:
        for line in FoI:
            sampleList.append(line.strip()) # Take a list of strains (line-separated, .txt file) and put them into a list

    startProcesses(sampleList)

if __name__ == '__main__':
    program(sys.argv[1]) #? As this file is run from the command line, the 0th argument is the file name (VCF_SA_Parallel.py), while the 1st argument is the file name of the list of strains