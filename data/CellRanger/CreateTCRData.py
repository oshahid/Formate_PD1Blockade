"""
This Python3 script process the TCR data for a given sample.


Inputs
------
Both inputs here are produced by the 'cellranger vdj' command.
- clonotypes.csv
- filtered_contig_annotations.csv

The third input is produced by the 'cellranger count' command (found in the filtered_feature_bc_matrix folder).
This file should remain gzipped.
- barcodes.tsv.gz

Outputs
-------
- TCR_metadata.csv
This single file has four columns: Cell barcode, TCR sequence, clone size, and AB status.
AB status describes whether or not a TCR has both an Alpha and a Beta chain, and this column
can take one of four values: 'notcr', 'alpha_chain_only', 'beta_chain_only', or 'AB_cell'.

More Info on AB Cells
---------------------
All TCRs are written in the form

 alpha-chain|beta-chain

AB cells have at least one alpha chain AND one beta chain.
For example, the TCR sequence CALSDRGGSNAKLTF|CASTNSAETLYF belongs
to an AB cell, but the TCR sequence |CASSDGGAYAEQFF is missing an
alpha chain, and would therefore not be considered an AB cell.
Instead, this TCR would be classified as 'beta_chain_only'.
Additionally, the clone size of all non-AB cells are defined to be 0.

Usage
____
# Call this script from the command line (single line):
python3 CreateTCRData.py clonotypes.csv filtered_contig_annotations.csv filtered_feature_bc_matrix/barcodes.tsv.gz TCR_metadata.csv

# Call this script from command line (multiple lines):
python3 CreateTCRData.py \
clonotypes.csv \
filtered_contig_annotations.csv \
filtered_feature_bc_matrix/barcodes.tsv.gz \
TCR_metadata.csv

"""

# Import packages
import sys                          # To read inputs and write output
import pandas as pd                 # To use data frames
from collections import Counter     # To count clone sizes

# Load data from command line
clonotypes = pd.read_csv(sys.argv[1], index_col = 0)
contigs = pd.read_csv(sys.argv[2], index_col = 0)
barcodes = pd.read_csv(sys.argv[3], compression = 'gzip', index_col = 0, header = None)

# Initialize TCR metadata
TCR_metadata = pd.DataFrame(index = barcodes.index,
                            columns = ['TCR', 'CloneSize_pre_filter', 'AB_status'])
TCR_metadata.index.name = 'Barcode'


# Create a function to get the alpha and beta chains given a clonotype
def Clonotype_to_TCR(clonotype):
    
    # If the clonotype is 'None', write the TCR as 'notcr'
    if str(clonotype) == 'None' or str(clonotype) == 'nan':
        TCR = 'notcr'
        return TCR
    
    # Get the Amino Acid chains from the clonotypes file
    cdr3s_aa = clonotypes.loc[clonotype, 'cdr3s_aa']
    
    # All chains are separated by semicolons
    all_chains = cdr3s_aa.split(';')
    
    # all_chains is a list with each entry beginning with either 'TRA:' or 'TRB:'
    
    # First get all the alpha and beta chains using list comprehensions
    # Use the [4:] notation to remove the first four characters in the chain
    # which correspond to the 'TRA:' and 'TRB:'
    alphas = [chain[4:] for chain in all_chains if 'TRA:' in chain]
    betas = [chain[4:] for chain in all_chains if 'TRB:' in chain]
    
    # Now that we have the alpha and beta chains separated,
    # start working with the alpha chains
    
    # If there are no alpha chains,
    if len(alphas) == 0:
        
        # Simply add a vertical bar
        TCR = '|'
        
    # Otherwise, use the first alpha chain
    else:
        TCR = alphas[0]
    
        # If there are more alpha chains, add them to the TCR, separated by a hypen
        if len(alphas) > 1:
            for chain in alphas[1:]:
                TCR += '-' + chain
    
        # Add a vertical bar to separate the alpha from the beta chains
        TCR += '|'
    
    # If there are no beta chains,
    if len(betas) == 0:
        
        # Return the TCR
        return TCR
    
    # Otherwise, add the first beta chain
    TCR += betas[0]
    
    # If there are more beta chains, add them to the TCR, separated by a hypen
    if len(betas) > 1:
        for chain in betas[1:]:
            TCR += '-' + chain
    
    return TCR


# Remove duplicates from filtered contig annotations
# This is because all the duplicated cell barcodes have the same clonotype
contigs_no_duplicates = contigs[~contigs.index.duplicated()]

# Loop over each cell, add the TCR to the metadata
for i,cell in enumerate(TCR_metadata.index):
    
    # If the cell is not in the contigs file, there is no TCR
    # Set the clone size to 0
    if cell not in contigs_no_duplicates.index:
        TCR_metadata.loc[cell, ['TCR', 'AB_status', 'CloneSize_pre_filter']] = ['notcr', 'notcr', 0]
    
    else:
        # Get clonotype
        clonotype = contigs_no_duplicates.loc[cell, 'raw_clonotype_id' ]
    
        # Add TCR to metadata
        TCR = Clonotype_to_TCR(clonotype)
        TCR_metadata.loc[cell, 'TCR'] = TCR

        # Check to see if the TCR exists
        if TCR == 'notcr':
            TCR_metadata.loc[cell, ['TCR', 'AB_status', 'CloneSize_pre_filter']] = ['notcr', 'notcr', 0] 
 
        # Check to see if the TCR has only an alpha chain
        # If it does, set the clone size to 0
        elif TCR[-1] == '|':         # Check last letter of TCR
            TCR_metadata.loc[cell, ['AB_status', 'CloneSize_pre_filter']] = ['alpha_chain_only', 0]

        # Check to see if the TCR has only a beta chain
        # If it does, set the clone size to 0
        elif TCR[0] == '|':          # Check first letter of TCR
            TCR_metadata.loc[cell, ['AB_status', 'CloneSize_pre_filter']] = ['beta_chain_only', 0]

        # Otherwise, it must be an AB cell
        else:
            TCR_metadata.loc[cell, 'AB_status'] = 'AB_cell'
 
# Use the Counter function to get the TCR frequencies
Frequencies = Counter(TCR_metadata['TCR'])

# Loop over each cell, add the frequency to each TCR if it belong to an AB cell
for i,cell in enumerate(TCR_metadata.index):
    
    # Check the status, see if it's an AB cell
    if TCR_metadata.loc[cell, 'AB_status'] == 'AB_cell':
        
        # Get the TCR
        TCR = TCR_metadata.loc[cell, 'TCR']
   
        # Add the frequency
        TCR_metadata.loc[cell, 'CloneSize_pre_filter'] = Frequencies[TCR]

# Write output
TCR_metadata.to_csv(sys.argv[4])
