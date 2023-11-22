import pandas as pd
import numpy as np
import scanpy as sc
import os


###################################
#                                 #
#     Function: Sort DE Genes     #
#                                 #
###################################

def SortDEgenes(df, pvals = 'Adj pvals', logFC = 'log2FC', GSEA_style = False):
    """
    Given a dataframe (df) where each row is a gene from a DE gene analysis, this function will
    sort the dataframe by doing one of the two:
    
    Regular Sorting:
        - Rank by adjusted p-values (lowest adjust p-value at the top of the list)
        - Rank by log-fold change (highest log-fold change at top of list)
        - Take the average rank of the two
        
    GSEA-Style Sorting:
        - Rank each gene by the product of (1) the sign of the log-fold change and (2) -log10(adjusted p-value)
    
    Inputs
    ------
    df: A dataframe that holds the results of a DE gene analysis.
    pvals (str): The column name that holds the adjusted p-values.
    logFC (str): The column name that holds the log-fold change.
    GSEA_style (bool): Whether or not to use GSEA-style sorting. Default is False.
    
    Output
    ------
    df: The sorted dataframe. If using GSEA-style, returns only the rank column.
    """
    
    df = df.copy()
    
    # Make sure the inputs are in the dataframe
    assert pvals in df.columns, "Column {} not found in the dataframe.".format(pvals)
    assert logFC in df.columns, "Column {} not found in the dataframe.".format(logFC)
    
    if GSEA_style == True:
        
        # Set all 0 pvals to minimum non-zero pval
        df.loc[df[pvals] == 0, pvals] = df.loc[df[pvals] > 0, pvals].min() 

        # Take sign(Log2FC) * -log10(Adj pval) as the rank of each gene
        df['Rank'] = np.sign(df[logFC]) * -np.log10(df[pvals])
        
        # Sort the ranks and keep only the rank column
        df = df.sort_values(by = 'Rank', ascending = False)
        df = df['Rank']
        
        
    else:

        # Rank by adjusted p-value (lowest p-value at top)
        df['Rank1'] = df[pvals].rank(ascending = False)

         # Rank by adjusted log-fold change (highest log-fold change at top)
        df['Rank2'] = df[logFC].rank(ascending = True)

        # Take the averange rank, sort by average rank
        df['Avg_Rank'] = (df['Rank1'] + df['Rank2'])/2
        
        # Sort the rank and drop three columns
        df = df.sort_values(by = 'Avg_Rank', ascending = False)
        df = df.drop(columns = ['Rank1', 'Rank2', 'Avg_Rank'])
    
    return df

##################################################
#                                                #
#     Function: Find DE Genes (One vs. Rest)     #
#                                                #
##################################################

def RankGenesGroupsOneVsRest(Adata, groupby, method = 'wilcoxon', sig = True, positive = True,
                            GSEA = False, save_name = None):
    """
    This function uses the scanpy 'rank_genes_groups' function to find DE genes in a 'one vs. rest' fashion
    for each group (Cluster 0 vs. the rest, Cluster 1 vs. the rest, etc.)
    By default, it uses the wilcoxon rank-sum mehod for comparison, and also has the option to save
    the gene list in such a way that it can be loaded directly into GSEA for further analysis.
    
    Inputs
    ------
    Adata: AnnData object.
    groupby (str): Which grouping to compare (can be 'seurat_clusters', 'samples', etc).
    method (str): Which statistical test to use.
    sig (bool): Whether or not to keep the significant genes (genes with adjusted p-values < 0.05).
    positive (bool): Whether or not to keep genes with a positive log-fold change.
    GSEA (bool): Whether or not to save data in a way that is directly compatible with GSEA input.
    save_name (str): If given, then save the DE gene list to an .xlsx file.
    
    Outputs
    -------
    DEgenes (dict): A Python dictionary that holds upregulated genes for each group.    
    """
    
    # Define dictionary keys, to be used later
    old_keys = ['names', 'pvals_adj', 'pvals', 'logfoldchanges', 'scores']
    new_keys = ['Gene Names', 'Adj pvals', 'Pvals', 'Log2FC', 'Scores']
    
    # Create an emtpy dictionary, DEgenes, to store the results
    DEgenes = {}

    # Find the DE genes, put them in the Adata object
    sc.tl.rank_genes_groups(Adata, groupby = groupby, method = method, n_genes = len(Adata.var_names),
                            use_raw=False)

    # Plot the DE genes (if not using GSEA pipeline)
    if GSEA == False:
        sc.pl.rank_genes_groups(Adata, n_genes = 20)

    # Extract the DE genes and statistics from the Adata object
    results = Adata.uns['rank_genes_groups']

    # groups are the names of the groups in which DE genes are upregulated
    groups = results['names'].dtype.names
    
    if groupby in ['leiden', 'seurat_clusters']:
        groups = sorted(groups, key = int)
    
    else:
        groups = sorted(groups)
        
    # If GSEA is true, we must be using all genes, not only significant genes
    # Force positive and sig = False
    if GSEA == True:
        positive = False
        sig = False
    
    
    for i,group in enumerate(groups):
        
        # Create an emtpy pandas dataframe
        df = pd.DataFrame(columns = new_keys)

        # Take results for each group and add them to a dataframe
        for old_key, new_key in zip(old_keys,new_keys):
            df[new_key] = results[old_key][group]

        # Keep genes that are significant
        if sig == True:
            df = df[ df['Adj pvals'] < 0.05 ]

        # Keep only positive log-fold changes
        if positive == True:
            df = df[ df['Log2FC'] > 0 ]

        # If there are any genes left,
        if len(df) > 0:
            
            # Set gene names as index
            df = df.set_index('Gene Names')
                
            # If saving for the GSEA pipeline, score each gene by
            # sign of Log2FC * -log10(Adjusted pvalue)
            # Since log(pval) is used here, set any adjusted pvalues that are currently 0 to a slightly
            # higher value
            
            if GSEA == True:
                
                # Make all the gene names upper case
                df.index = list(map(str.upper, df.index))
                
                # Sort values
                df = SortDEgenes(df, pvals = 'Adj pvals', logFC = 'Log2FC', GSEA_style = True)

                
            # Otherwise, take the average rank of Adj pvals and log fold change
            else:
                
                # Sort
                df = SortDEgenes(df, pvals = 'Adj pvals', logFC = 'Log2FC', GSEA_style = False)
                    
        else:
            print('Group {} has no significant genes.'.format(group))

        # Save the dataframe into the DEgenes object
        DEgenes[group] = df
        
    # If save_name is given, save results to an Excel file
    if save_name is not None:

        with pd.ExcelWriter(save_name) as writer:

            # For each key,
            for key in DEgenes.keys():

                # Extract the pandas dataframe, save it to an Excel sheet
                DEgenes[key].to_excel(writer, sheet_name = key)
                
    return DEgenes


##############################################
#                                            #
#     Function: Find DE Genes (Pairwise)     #
#                                            #
##############################################

def RankGenesGroupsPairwise(Adata, groupby, method = 'wilcoxon', sig = True, positive = True, GSEA = False,
                            save_name = None):

    """
    This function uses the scanpy 'rank_genes_groups' function to find DE genes in a pairwise fashion
    for each group (Group A vs. Group B, Group B vs. Group A, etc.)
    By default, it uses the wilcoxon rank-sum mehod for comparison, and also has the option to save
    the gene list in such a way that it can be loaded directly into GSEA for further analysis.
    
    Inputs
    ------
    Adata: AnnData object.
    groupby (str): Which grouping to compare (default is 'leiden' or other grouping method.)
    method (str): Which statistical test to use.
    sig (bool): Whether or not to keep the significant genes (genes with adjusted p-values < 0.05).
    positive (bool): Whether or not to keep genes with a positive log-fold change.
    GSEA (bool): Whether or not to save data in a way that is directly compatible with GSEA input.
    save_name (str): If given, then save the DE gene list to an .xlsx file.
    
    Outputs
    -------
    DEgenes (dict): A Python dictionary that holds upregulated genes for each group.    
    """

    # Define dictionary keys, to be used later
    old_keys = ['names', 'pvals_adj', 'pvals', 'logfoldchanges', 'scores']
    new_keys = ['Gene Names', 'Adj pvals', 'Pvals', 'Log2FC', 'Scores']
    
    # Create an emtpy dictionary, DEgenes, to store the results
    DEgenes = {}

    # For each category in the 'groupby' data, make sure each group has more than 10 cells
    group_counts = Adata.obs[groupby].value_counts()
    groups = list(group_counts[group_counts > 10].index)
    
    # If GSEA is true, we must be using all genes, not only significant genes
    # Force positive and sig = False
    if GSEA == True:
        positive = False
        sig = False
    
    # For each group in groups, make the group a reference
    for ref in groups:
        print('The reference group is {}\n'.format(ref))
    
        # Find the DE genes, put them in the Adata object
        sc.tl.rank_genes_groups(Adata, groupby = groupby, groups = groups, reference = ref, method = method,
                                n_genes = len(Adata.var_names), use_raw = False )

        # Plot the DE genes (if not using GSEA pipeline)
        if GSEA == False:
            sc.pl.rank_genes_groups(Adata, n_genes = 20)
    
        # Extract the DE genes and statistics from the Adata object
        results = Adata.uns['rank_genes_groups']
        
        # s_groups (sample groups) are the names of the groups in which DE genes are upregulated
        s_groups = results['names'].dtype.names
        
        
        # For each sample in the comparison,
        for sample in s_groups:
            
            # Create an emtpy pandas dataframe
            df = pd.DataFrame(columns = new_keys)

            # Load the data from the results object to the dataframe
            for new_key,old_key in zip(new_keys,old_keys):
                df[new_key] = results[old_key][sample]
                
            # Set gene names as index
            df = df.set_index('Gene Names')
                
            # If saving for the GSEA pipeline, score each gene by
            # sign of Log2FC * -log10(Adjusted pvalue)
            # Since log(pval) is used here, set any adjusted pvalues that are currently 0 to a slightly
            # higher value
            
            if GSEA == True:
                
                # Make all the gene names upper case
                df.index = list(map(str.upper, df.index))
                
                # Sort values
                df = SortDEgenes(df, pvals = 'Adj pvals', logFC = 'Log2FC', GSEA_style = True)

            # Otherwise, take the average rank of Adj pvals and log fold change
            else:
                
                # Keep genes that are significant
                if sig == True:
                    df = df[ df['Adj pvals'] < 0.05 ]

                # Keep only the positive log fold changes
                if positive == True:
                    df = df[ df['Log2FC'] > 0 ]
                
                # Sort
                df = SortDEgenes(df, pvals = 'Adj pvals', logFC = 'Log2FC', GSEA_style = False)
            
            # Save the dataframe into the DEgenes object
            DEgenes[sample + ' vs. ' + ref] = df
            
    # If save_name is given, save results to an Excel file
    if save_name is not None:
        
        with pd.ExcelWriter(save_name) as writer:

            # For each key,
            for key in DEgenes.keys():

                # Save the pandas dataframe to an Excel sheet
                DEgenes[key].to_excel(writer, sheet_name = key)
        
        
    return DEgenes

###########################################
#                                         #
#     Function: Prepare Volcano Plots     #
#                                         #
###########################################


def PrepareForVolcanoPlots(Adata, groupby, key, save_name = None):
    """
    Given an Adata object, run a pair-wise comparison between two groups and save the results in a way
    that it can loaded into R to create a Volcano plot.
    
    Inputs
    ------
    Adata: AnnData object.
    groupby (str): Name of metadata column to run comparison between.
    key (str): Name of comparison to use (format: 'GroupA vs. GroupB')
    save_name (str): Name of file when saving to disk.
    
    Outputs
    -------
    df: Pandas Dataframe that holds the results of the DE gene analysis.
    """
    
    # Make sure 'groupby' is in Adata object
    assert groupby in list(Adata.obs.columns), "Variable {} not found in Adata.obs.".format(groupby)
    
    # Find DE genes
    DEgenes = RankGenesGroupsPairwise(Adata, groupby = groupby, sig = True, positive = False)
    
    # Make sure the key is a valid comparison
    assert key in list(DEgenes.keys()), "Key {} not found among the comparisons.".format(key)
    
    # Keep only the specified comparison
    df = DEgenes[key]
    
    # If there are any adjusted p-values == 0, set them to the minimum non-zero adjusted p-value
    df.loc[df['Adj pvals'] == 0, 'Adj pvals'] = df.loc[df['Adj pvals'] > 0, 'Adj pvals'].min() 
    
    # Add a column to calculate the -Log10(Adj pvals)
    df['-Log10(Adj pvals)'] = -np.log10(df['Adj pvals'])
    
    # Save file
    if save_name is not None:
        df.to_csv(save_name)
    return df

###########################################
#                                         #
#     Function: GSEA Helper Functions     #
#                                         #
###########################################

def ReversedKey(key):
    """
    ReversedKey will take a key from a given DEgenes dictionary and return the reversed version of the key. For example, given
    the key 'Group1 vs. Group2', the function call ReversedKey(key) will return 'Group2 vs. Group1'.

    Input
    -----
    key (str): A string that names the comparison between 2 groups, with the phrase ' vs. ' inside of it (Ex: 'A vs. B').

    Output
    ------
    reversedkey (str) The reversed key (Ex: 'B vs. A').
    """
    return key.split()[2] + ' vs. ' + key.split()[0]

def GetUniqueComparisons(DEgenes, sort = True):
    """
    GetUniqueComparisons will take all the keys of the DEgenes dictionary and will return only the unique keys
    (where we treat the key and its reverse as the same comparison). For example, 'A vs. B' and 'B vs. A' will
    be considered as the same comparison, and only 'A vs. B' will be returned.

    Inputs
    ------
    DEgenes (dictionary): Holds the results of any DE gene analyses.

    Outputs
    -------
    UniqueComparisons: A Python list that holds the names of all the unique comparisons in the DEgenes object.
    """

    # Get all the keys/comparisons from the DEgenes dictionary
    all_keys = list(DEgenes.keys())
    
    # Initialize empty list
    UniqueComparisons = []
    
    # Go through each key, check to see if the key (or the key's reverse) is in UniqueComparisons
    for key in all_keys:

        # If it is, then skip
        if key in UniqueComparisons or ReversedKey(key) in UniqueComparisons:
            continue
        
        # Otherwise, add it to the list
        else:
            UniqueComparisons.append(key)

    # Sort, if specified
    if sort == True:
        UniqueComparisons = sorted(UniqueComparisons)
            
    return UniqueComparisons

########################################
#                                      #
#     Function: GSEA CLI Interface     #
#                                      #
########################################

def GSEA_CLI(name_of_analysis, rank_file, output, gene_set_name):
    """
    GSEA_CLI is a Python interface for running the GSEA CLI.

    Inputs
    ------
    name_of_analysis (str): The name of the analysis, usually in the form of 'key_geneset'. For example,
    this could be 'GroupA_vs_GroupB_hallmark'.
    rank_file (str): The file path to the rank file.
    output (str): The output directory path.
    gene_set_name (str): The name of the gene set to run GSEA against.

    Outputs
    -------
    None, since the outputs will saved to disk.
    """

    # Set up the possible gene set names, and file paths to each gene set
    gene_set_names = ['hallmark', 'C2', 'C5', 'C7', 'Topics']
    
    gene_sets = {'hallmark':'ftp.broadinstitute.org://pub/gsea/gene_sets/h.all.v7.1.symbols.gmt',
                 'C2':'ftp.broadinstitute.org://pub/gsea/gene_sets/c2.all.v7.1.symbols.gmt',
                 'C5':'ftp.broadinstitute.org://pub/gsea/gene_sets/c5.all.v7.1.symbols.gmt',
                 'C7':'ftp.broadinstitute.org://pub/gsea/gene_sets/c7.all.v7.1.symbols.gmt',
                 'Topics':'data/GeneLists/TopicSigsUpper.gmx'}

    # Make sure the input 'gene_set_name' is valid
    assert gene_set_name in gene_set_names, "gene_set_name must be one of: hallmark, C2, C5, C7, Topics"

    # Define command-line variables
    os.environ['TITLE'] = name_of_analysis
    os.environ['RANK_FILE'] = rank_file
    os.environ['OUT'] = output
    os.environ['GENE_SET'] = gene_sets[gene_set_name]
    os.environ['GSEA'] = 'results/scRNA/OldResults/GSEA/GSEA_4.0.3/gsea-cli.sh'

    # This format can be used with the Python 'subprocess.call' module
    GSEA_commands = ["bash", "$GSEA GSEAPreranked", "-gmx $GENE_SET", "-collapse No_Collapse", "-mode Max_probe",
                     "-norm meandiv", "-nperm 1000", "-rnk $RANK_FILE", "-scoring_scheme classic", "-rpt_label $TITLE",
                     "-create_svgs true", "-include_only_symbols true", "-make_sets true", "-plot_top_x 20", "rnd_seed 1",
                     "-set_max 500", "-set_min 15", "-zip_report false", "-out $OUT"]

    # This format can be used with the Python 'os.system' module
    GSEA_commands = ' '.join(GSEA_commands)

    # Run GSEA
    os.system(GSEA_commands)
    
    print('\tFinished analysis on gene set {}'.format(gene_set_name))
    return

###################################
#                                 #
#     Function: GSEA Pipeline     #
#                                 #
###################################

def GSEA_Pipeline(DEgenes, output_dir, keys = None, gene_set_names = None):
    """
    GSEA_Pipeline is a Python wrapper designed to run GSEA for each comparison in the DEgenes object.

    Inputs
    ------
    DEgenes (dict): A Python dictionary that comes as a result of the DE gene analyses.
    output_dir (str): The directory path and name where the rank file will be created and the GSEA results stored.
    keys (list): A list of the DEgenes keys to use. By default, all unique keys will be used.
    gene_set_names (list): A list of the gene sets to use for GSEA. By default,
    all gene sets (hallmark, C2, C5, C7, and topic-modeling signatures) will be used.

    Outputs
    -------
    None, since the outputs will saved to disk.
    """

    # Make sure the output directory exists
    if not os.path.isdir(output_dir):
        print('Output directory does not exist - making one now: {}'.format(output_dir))
        os.mkdir(output_dir)
    else:
        print('Output directory is {}'.format(output_dir))

    # Set up the pipeline
    # Default: All gene sets are used
    if gene_set_names is None:
        gene_set_names = ['hallmark', 'C2', 'C5', 'C7', 'Topics']
    
    # Otherwise, make sure the gene sets are valid
    else:
        for gene_set in gene_set_names:
            assert gene_set in ['hallmark', 'C2', 'C5', 'C7', 'Topics'], "Gene sets must be one of: hallmark, C2, C5, C7, Topics"
    
    gene_sets = {'hallmark':'ftp.broadinstitute.org://pub/gsea/gene_sets/h.all.v7.1.symbols.gmt',
                 'C2':'ftp.broadinstitute.org://pub/gsea/gene_sets/c2.all.v7.1.symbols.gmt',
                 'C5':'ftp.broadinstitute.org://pub/gsea/gene_sets/c5.all.v7.1.symbols.gmt',
                 'C7':'ftp.broadinstitute.org://pub/gsea/gene_sets/c7.all.v7.1.symbols.gmt',
                 'Topics':'data/GeneLists/TopicSigsUpper.gmx'}

    # Check the keys. By default, use all the unique comparisons
    if keys is None:
        keys = GetUniqueComparisons(DEgenes)

    # If keys are given, make sure they exist in the DEgenes object
    else:
        for key in keys:
            assert key in list(DEgenes.keys()), "All keys must be in the DEgenes.keys() slot."

    # Start running the pipeline
    for key in keys:
        print('The comparison being made is: {}'.format(key))

        # Extract the data frame
        df = DEgenes[key]

        # Parse the key to create the name of the rank file by replacing spaces with underscores
        # Example: 'GroupA vs. GroupB' will become 'GroupA_vs_GroupB'
        key = key.replace(' vs. ', '_vs_')

        # Create the rank file and save it to disk
        rank_file = os.path.join(output_dir, key+'.rnk')
        df.to_csv(rank_file, header = False, sep = '\t')

        # Print output for user
        print('The rank file is: {}'.format(rank_file))

        # Run GSEA for the DEgenes key/comparison for each gene set (hallmark, C2, C5, C7, topics)
        for gene_set_name in gene_set_names:

            # Print output for user
            print('\n\tRunning on gene set: {}\n'.format(gene_set_name))

            # The name of the analysis for this key/comparison should be 'key_geneset'
            # Example: 'GroupA_vs_GroupB_hallmark'
            name_of_analysis = key + '_' + gene_set_name

            # Run GSEA CLI
            GSEA_CLI(name_of_analysis, rank_file, output_dir, gene_set_name)

    print('Finished GSEA on all gene sets.')
    return