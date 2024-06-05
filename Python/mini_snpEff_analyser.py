# Takes a .vcf file as an input, annotates it, and pulls out all (MA) mutants (primarily, this is for getting consistent WB gene IDs from the MA line files). This version is unbiased towards the original annotation results (i.e. from the papers), and prioritises mutants that will be within genes
#? This is primarily used for the MA data - a similar (non-mini) version was deployed onto Genoa for the CaeNDR database analysis

import json
import os
import subprocess
import re

PATH = f'{os.path.dirname(__file__)}'.replace('\\', '/') # The path to the working folder of the project

VCF_FILE_NAME = 'Konrad Mutants' # The input annotated VCF file
JSON_FILE_NAME = 'Konrad Data' # The name of the JSON output file
CONVERT_TO_SNPEFF = True

NUMERAL_TO_NUM = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10, 'XI':11, 'XII':12, 'MtDNA':13}
SNPEFF_MUTANT_TO_JSON_MUTANT = {'upstream_gene_variant':'Modifier/-', 'intron_variant':'Intron/-', 'downstream_gene_variant':'Modifier/-', 'intergenic_region':'IG/-',
                                '3_prime_UTR_variant':'UTR/-', 'splice_donor_variant&intron_variant':'Intron/-', 'splice_region_variant&intron_variant':'Intron/-',
                                '5_prime_UTR_variant':'UTR/-', 'splice_acceptor_variant&intron_variant':'Intron/-', '5_prime_UTR_premature_start_codon_gain_variant':'UTR/-',
                                'splice_region_variant&non_coding_transcript_exon_variant':'NC/-', 'splice_region_variant':'Intron/-', 'non_coding_transcript_exon_variant':'NC/-',
                                'splice_region_variant&non_coding_exon_variant':'NC/-'}
AA_THREE_TO_ONE = {'ala':'A', 'arg':'R', 'asn':'N', 'asp':'D', 'cys':'C', 'glu':'E', 'gln':'Q', 'gly':'G', 'his':'H', 'ile':'I', 
                  'leu':'L', 'lys':'K', 'met':'M', 'phe':'F', 'pro':'P', 'ser':'S', 'thr':'T', 'trp':'W', 'tyr':'Y', 'val':'V', 'ter':'*', 'unk':'?'}

with open(f'{PATH}/Coding/Resources/Gene Lookup (JSON [Transcript ID~Gene ID]).json', 'r') as FoI:
    NAME_TO_WBID = json.loads(FoI.read())

def loadJSONMutants(): #* Reads in the mutants from the JSON file
    mutants = [] #? The list of all the mutants from the MA data, in the standard JSON format (locus, gene [name], subsitution, mutationType)
    with open(f'{PATH}/Coding/Resources/{JSON_FILE_NAME}.json', 'r') as FoI: # Open the file and grab all the MA mutants
        for line in FoI:
            mutants.append(json.loads(line.strip()))
    return mutants

def saveJSONMutants(mutants): #* Saves the json file of mutants, now properly named after the annotation
    with open(f'{PATH}/Coding/Resources/{JSON_FILE_NAME}.named.json', 'w') as FoI:
        for mutant in mutants:
            FoI.write(json.dumps(mutant) + '\n')

def annotateVCF(): #* Annotates the VCF file using SnpEff
    with open(f'{PATH}/Coding/Resources/{VCF_FILE_NAME}.annotated.vcf', 'w') as output:
        subprocess.run(f'java+-jar+{PATH}/snpEff/snpEff.jar+-noStats+WBcel235.75+{PATH}/Coding/Resources/{VCF_FILE_NAME}.vcf'.split('+'), stdout = output)

scoring = [0, 0, 0, 0, 0]
errs = [0]   
def getWBID(data, mutant): #* For a mutant, reads all the annotations and determines the best one (i.e. the one most likely to be associated with a protein-coding region)
    locus = f'{data[0]}:{data[1]}' #? Gets the chromosome number and position (within the chromosome) - i.e. the locus
    if locus == mutant['locus']: #! If the locus derived from the VCF file does NOT match that of the JSON data, then something's gone wrong - error
        substitution = f'{data[3]} -> {data[4]}' #? Gets the SNP
        annotationData = data[7].split(',') #? Gets the annotation data from the VCF file (a list of annotations, which are comma-separated)
        featureScore = 6
        for annotation in annotationData:
            annotation = annotation.split('|')
            features = annotation[1].split('&')
            for feature in features:
                if feature in ['missense_variant', 'synonymous_variant', 'stop_gained', 'stop_lost', 'stop_retained_variant', 'start_lost']:
                    if 1 < featureScore:
                        topAnnotation = annotation
                        featureScore = 1
                elif feature in ['splice_donor_variant', 'splice_acceptor_variant', 'start gain', '5_prime_UTR_premature_start_codon_gain_variant']:
                    if 2 < featureScore:
                        topAnnotation = annotation
                        featureScore = 2
                elif feature in ['intron_variant', 'splice_region_variant']:
                    if 3 < featureScore:
                        topAnnotation = annotation
                        featureScore = 3
                elif feature in ['3_prime_UTR_variant', '5_prime_UTR_variant']:
                    if 4 < featureScore:
                        topAnnotation = annotation
                        featureScore = 4
                elif feature in ['intergenic_region', 'upstream_gene_variant', 'downstream_gene_variant', 'non_coding_exon_variant', 'non_coding_transcript_exon_variant']:
                    if 5 < featureScore:
                        topAnnotation = annotation
                        featureScore = 5
                else:
                    print('ERR: Feature not listed!')
        
        scoring[featureScore - 1] += 1

        mutationType = topAnnotation[1]
        geneID = topAnnotation[4]
        aaSub = topAnnotation[10]

        if aaSub != '':
            aminos = re.findall(r'([A-Z][a-z][a-z])|\*', aaSub)
            if len(aminos) == 1: #! If this is only one, then there's no second amino acid, nor stop codon (i.e. there's no predicted change). Usually, this means that it's unknown, '?'.
                if '?' in aaSub:
                    aminos.append('unk')
            for i in range(2): # If there's more than two, it'll be because of extension (HGVS) notation or something, which we don't particularly care about
                if aminos[i] in ['', '*']: # If it's a stop codon (denoted with a *, because sometime's they're 'ter' instead. Also, re seems to not read the star because it's a special character, so outputs '' to the list)
                    aminos[i] = 'ter'
                elif aminos[i] == '?': # A question mark means that the outcome is unknown (unk)
                    aminos[i] = 'unk' 
                aminos[i] = AA_THREE_TO_ONE[aminos[i].lower()]
            wtaa, mutaa = aminos[0], aminos[1] #? re.findall outputs in order, so the wildtype will always come before the mutant in the list
            mutType = (f'Exon/{wtaa} -> {mutaa}')
            if wtaa == mutaa: # If it's a synonymous mutation
                mutType = 'Exon/Synonymous'
        else:
            mutType = SNPEFF_MUTANT_TO_JSON_MUTANT[mutationType]
        mutant['mutationType'] = mutType

        if geneID[:6] == 'WBGene':
            mutant['gene'] = geneID
        else:
            try:
                mutant['gene'] = NAME_TO_WBID[geneID]
            except KeyError:
                print(f'ERR: Name ({topAnnotation[3]}) is not known.')
                errs[0] += 1
        
        return mutant
    else: #! If the loci don't match
        print(locus, mutant)
        print('ERR: JSON and VCF data do not match. Are they both properly ordered?')
        return None


#- annotateVCF()
mutants = loadJSONMutants()
with open(f'{PATH}/Coding/Resources/{VCF_FILE_NAME}.annotated.vcf', 'r') as FoI: # Open the VCF file that has the same data as the JSON file (i.e. this is a VCF that has been converted from a bed file, which itself was made from the JSON data)
    mutantNum = 0 # Keep track of the list index for the 'mutants' list - every time a line of annotated data is read, the next mutant is being considered
    errors = 0
    for line in FoI:
        if line[0] != '#': # If the line actually has MA data on it (and isn't just a comment in the VCF file)
            mutant = mutants[mutantNum] # Get the associated JSON data (namely for the mutationType)
            data = line.strip().split('\t') # Split the line by tabs. This produces eight elements (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO/ANN)
            mutant = getWBID(data, mutant)
            if mutant: # If a valid (updated) mutant was returned
                mutants[mutantNum] = mutant # Update the mutant in the list of mutants
            else:
                print('ERR: Mutant has not been updated.')
            mutantNum += 1

print(errs)

saveJSONMutants(mutants)

'''
[1] - missense_variant, synonymous_variant, stop_gained, stop_lost, stop_retained_variant, start_lost
[2] - splice_donor_variant, splice_acceptor_variant, start gain, 5_prime_UTR_premature_start_codon_gain_variant
[3] - intron_variant, splice_region_variant
[4] - 3_prime_UTR_variant, 5_prime_UTR_variant 
[5] - intergenic_region, upstream_gene_variant, downstream_gene_variant, non_coding_exon_variant, non_coding_transcript_exon_variant 
'''

