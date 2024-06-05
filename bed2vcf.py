# Takes a .bed (or .json) file and converts it into a .vcf

from datetime import datetime
import os

PATH = f'{os.path.dirname(__file__)}/Resources'.replace('\\', '/') # The path to the 'Resources' folder of the project
#- FILE_NAME = 'Konrad Mutants (Zeroed)'
FILE_NAME = 'Saxena LiftOver'

CHROM_TO_NUMERAL = {'chr1':'I', 'chr2':'II', 'chr3':'III', 'chr4':'IV', 'chr5':'V', 'chr10':'X', 'chrMT':'MtDNA'}

bedLines = []
with open(f'{PATH}/{FILE_NAME}.bed', 'r') as FoI:
    for line in FoI:
        if line[0] != '#':
            bedLines.append(line.strip())

with open(f'{PATH}/{FILE_NAME}.vcf', 'w') as FoI:
    currentDateTime = datetime.now()
    FoI.write('##fileformat=JAK_VCFv1\n')
    FoI.write('##reference="C. elegans N2 WBcel235"\n')
    FoI.write(f'##filedate={currentDateTime.day:02}/{currentDateTime.month:02}/{currentDateTime.year}@{currentDateTime.hour:02}:{currentDateTime.minute:02}:{currentDateTime.second:02}\n')
    FoI.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    for line in bedLines:
        data = line.split('\t')
        chrom = CHROM_TO_NUMERAL[data[0]]
        pos = data[2]
        id = '.'
        ref = data[3][0]
        alt = data[3][-1]
        qual = '.'
        filter = 'PASS'
        info = '.'
        FoI.write(f'{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\n')




'''
1	CHROM	The name of the sequence (typically a chromosome) on which the variation is being called. This sequence is usually known as 'the reference sequence', i.e. the sequence against which the given sample varies.
2	POS	The 1-based position of the variation on the given sequence.
3	ID	The identifier of the variation, e.g. a dbSNP rs identifier, or if unknown a ".". Multiple identifiers should be separated by semi-colons without white-space.
4	REF	The reference base (or bases in the case of an indel) at the given position on the given reference sequence.
5	ALT	The list of alternative alleles at this position.
6	QUAL	A quality score associated with the inference of the given alleles.
7	FILTER	A flag indicating which of a given set of filters the variation has failed or PASS if all the filters were passed successfully.
8	INFO    	An extensible list of key-value pairs (fields) describing the variation. See below for some common fields. Multiple fields are separated by semicolons with optional values in the format: <key>=<data>[,data].

java -jar ./snpEff/snpEff.jar -v -noStats WBcel235.75 "./Coding/Resources/Konrad Unfiltered Mutants.vcf" > "./Coding/Resources/Konrad Unfiltered Mutants.annotated.vcf"
"./Coding/Resources/Saxena/Saxena 234 Mutants.vcf" > "./Coding/Resources/Saxena/Saxena 234 Mutants.annotated.vcf"
'''