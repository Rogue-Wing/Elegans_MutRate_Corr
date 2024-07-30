### The Mutation Rate of _Caenorhabditis elegans_ Genes is Influenced by Selection
## GitHub File Dump
# Introduction
This is the collection of code and files that I produced during my Masters' project. The primary purpose of all the code is to sanitise and manipulate large volumes of _C. elegans_ genomics data in order to calculate the mutation rates, and amount of purifying selection, of _C. elegans_ genes, from which analyses of the relationship between these two variables can be drawn. This GitHub contains the Python and R code, as well as example files for these program's inputs/outputs.

While the Python and R code can be found in full, the resources files are truncated, in almost all cases, to the first 10 entries (as many are .json or .csv lists of some description or another). Equally, those produced during the processing of the AG dataset have not been included here, as the pipelines that convert the AG dataset from it's initial collated .vcf.gz file into strain-separated .json files automatically delete the intermediate files to conserve space.

# Running the code
The R and Python code can be kept within the same directory, and both require the 'Resources' file to sit as a child directory within this. In the majority of cases, the code will take any input files from this Resources section, and output any files back to it.

Please note that running the AG (CaeNDR) pipeline requires the use of bcftools and SnpEff, as well as the hard-filtered variants file from the CaeNDR database. This is best run on a dedicated cluster, owing to the compute time and power required to properly analyse these data.
