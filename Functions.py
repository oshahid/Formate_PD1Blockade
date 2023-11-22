import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
import scanpy as sc
import anndata



################################################
#                                              #
#     Function: Make Non-Normalized Counts     #
#                                              #
################################################

def RawCounts(Adata_, target_sum = 1e4):
    """
    This function returns the Adata object with the original counts by using the 'norm_factor' column
    in the metadata. This should be very rarely used, as most of the analyses use normalized data.
    """

    Adata = Adata_.copy()

    
    # Get the total number of counts per cell
    counts_per_cell = np.array(Adata.obs['norm_factor'])
    
    # Get the count matrix X, un-logarithmize it, and subtract 1
    X = Adata.X
    X = np.power(np.e, X)    # Since natural log was used for normalization
    X = X - 1                # Since log(1+X) was used for normalization
    
    # Divide by the given target_sum and multiply each row by the normalization factor
    X = X / target_sum * counts_per_cell[:, np.newaxis]
    
    # Round each element, since there was likely rounding error, and add to Adata object
    Adata.X = np.round(X)
    Adata.X = Adata.X.astype(int)
    
    return Adata

############################################
#                                          #
#     Function: Initialize Empty Plots     #
#                                          #
############################################

def MakeFigAxes(ncols, n_axes, figsize=(10,8), single = False, constrained_layout = False):
    """
    Helper function to create matplotlib subplots.
    
    Inputs
    ------
    ncols (int): Number of columns to make.
    n_axes (int): Number of axes (subplots) that are needed.
    figsize (tuple): Width and height of total figure.
    single (bool): Whether or not to return a single ax object instead of an array of axes.
    constrained_layout (bool): Whether or not to use the constrained_layout option for matplotlib.
    
    Outputs
    -------
    fig: Matplotlib figure.
    axes: Matplotlib axes as an array.
    """
    
    
    # Determine how many rows are needed
    if n_axes%ncols == 0:
        nrows = int(n_axes/ncols)
    else:
        nrows = int(1+ int(n_axes/ncols))
    
    # Special case if only using one row
    if nrows == 1:
        fig, axes = plt.subplots(nrows = nrows, ncols = n_axes, figsize = figsize,
                                 constrained_layout = constrained_layout)
        axes = np.array(axes)
        axes.shape = (n_axes,)
    
    # For more than 1 row
    else:
        fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize,
                                 constrained_layout = constrained_layout)
    
        # Delete unnecessary axes
        if n_axes != nrows*ncols:
            n_axes_to_delete = (nrows*ncols - n_axes)
            for i in range(1,n_axes_to_delete+1):
                ax = axes.flatten()[-i]
                fig.delaxes(ax)
                
    # Return single axis, if specified
    if single == True:
        axes = axes[0]
    return fig, axes

###################################################
#                                                 #
#     Function: Plot Categorical Data on UMAP     #
#                                                 #
###################################################

def SplitPlot(Adata, groupby, ncols = 3, palette = None, figsize = (8,6), save_name = None):
    
    # Get groups
    groups = sorted(Adata.obs[groupby].unique())
    
    # Make figure and axes
    # Number of axes/subplots is equal to number of gene lists
    # Number of columns is default to 3, unless there are fewer than 3 gene lists
    n_axes = len(groups)
    ncols = min(ncols, n_axes)
    fig, axes = MakeFigAxes(ncols = ncols, n_axes = n_axes, figsize = figsize)

    if palette is not None:
        palette = palette[:len(groups)]
    
    # Loop over each group
    for i,group in enumerate(groups):
        
        ax = axes.flatten()[i]
        
        sc.pl.umap(Adata, color = groupby, groups = [group], ax = ax, show = False, palette = palette,
                   title = str(group), legend_loc = None, return_fig=False)
        
    plt.tight_layout()
    if save_name is not None:
        plt.savefig(save_name)
    plt.show()
    return

################################################
#                                              #
#     Function: Find Intersection of Genes     #
#                                              #
################################################

def IntersectionOfGenes(object1, object2, sort = False, upper = False):
    """
    This function returns the genes in the intersection of two objects. If the object is an Adata object,
    then the gene names are used (Adata.var_names). Otherwise, the function expects a set or a list of genes.
    
    This code is NOT optimized for speed - rather, the goal is to return the intersection of genes:
    - in the order of the first list
    - regardless of whether or not either list is upper or lower case
    
    Inputs
    ------
    object1: An Adata object, a list of genes, or a set of genes.
    object2: An Adata object, a list of genes, or a set of genes.
    sort (Bool): Sorts the final output. Default is False (genes will be return in the order of object1)
    upper (Bool): Upper-cases the final output. Default is False (genes will be returned in the original case)
    """
    
    ### Parse the inputs and see if they're lists or Adata objects ###
    
    # Initialize an empty object to store the lists
    Lists = [0, 0]
    
    # Loop over each input object
    for i,obj in enumerate([object1, object2]):
        
        # If the object is an Adata object, extract the genes as a list
        if type(obj) == anndata._core.anndata.AnnData:
            Lists[i] = list(obj.var_names)
        
        else:
            # Otherwise,  object is assumed to be a list or set - therefore, turn it into a list.
            Lists[i] = list(obj)
            
            
    # Now we have two gene lists
    list1 = Lists[0]
    list2 = Lists[1]
    
    ### Find the intersection of both lists ###
    
    # Initialize a list of False values
    TrueFalse = [False] * len(list1)
    
    # Make both genes upper case to find the intersection
    List1 = list(map(str.upper, list1))
    List2 = list(map(str.upper, list2))
    
    # For each gene in the first list, if it exists in the second list,
    # then set its truth value to True
    for i,gene in enumerate(List1):
        if gene in List2:
            TrueFalse[i] = True
            
    # Create an object to be returned by doing the following:
    ListToReturn = np.array(list1)             # Turn the original list into a NumPy array
    ListToReturn = ListToReturn[TrueFalse]     # Keep only the genes that are in the intersection
    ListToReturn = list(ListToReturn)          # Return it to a list to to use list methods
          
    # Sort the genes, if specified
    if sort == True:
        ListToReturn.sort()
        
    # Upper case all the genes, if specified
    if upper == True:
        ListToReturn = list(map(str.upper, ListToReturn))
    
    # Return the list
    return ListToReturn

#####################################
#                                   #
#     Function: Cap Gene Scores     #
#                                   #
#####################################

def CapScores(score_list, thresh_low = 3, thresh_high = 3):
    """
    This function takes a list of gene scores (one score per cell) and fits it to a normal distribution,
    estimating a mean (mu) and standard deviation (SD). Then it caps the scores that are above a certain
    threshold and caps the scores that are below another threshold.
    
    This action should only be done for visualization purposes, and not for doing any statistical tests.
    
    Inputs
    ------
    score_list: A list of gene scores.
    thresh_low:  The number of standard deviations away from the mean to cap the lower scores.
    thresh_high: The number of standard deviations away from the mean to cap the upper scores.
    
    For both thresholds, the default behavior is that any scores that are 3 SDs above or
    below the mean will be capped.
    
    Output
    ------
    score_list: The capped list of gene scores.
    """
    
    # Fit the gene scores to a normal distribution
    mu, SD = scipy.stats.norm.fit(score_list)
    
    # Cap the scores
    thresh_high = mu + thresh_high * SD
    thresh_low  = mu - thresh_low  * SD
    
    for e,ele in enumerate(score_list):
        
        if ele >= thresh_high:
            score_list[e] = thresh_high
            
        if ele <= thresh_low:
            score_list[e] = thresh_low
            
    return score_list

#####################################
#                                   #
#     Function: Load Gene Lists     #
#                                   #
#####################################

def LoadGeneSig(file_name, sheet_name = None):
    """
    Function to load gene lists from either .txt or Excel files
    """
    
    # If last 4 characters of file name are .txt,
    if file_name[-4:] == '.txt':
        with open(file_name, 'r') as f:
            genes = f.read().splitlines()
    
    # Otherwise, assume it's an Excel file
    else:
        genes = pd.read_excel(file_name, sheet_name = sheet_name, index_col = 0, header = None, engine = 'openpyxl')
        genes = list(genes.index.str.upper())
        
    return genes

###################################
#                                 #
#     Function: Keep AB Cells     #
#                                 #
###################################

def KeepABcells(Adata_):
    """
    This function returns an Adata object of only AB cells, which are cells that have at least
    one alpha chain and one beta chain.
    """
    
    # Make a copy of the data, as to not change the original object
    # Keep only AB cells
    Adata = Adata_[ Adata_.obs['AB_status'] == 'AB_cell' ].copy()
    
    # Remve unused categories
    Adata.obs['AB_status'].cat.remove_unused_categories(inplace=True)
    
    return Adata

######################################
#                                    #
#     Function: Check Gene Lists     #
#                                    #
######################################

def CheckGeneLists(gene_lists = None, gene_lists_names = None, gene_dict = None):
    """
    This helper function ensures that gene lists and their names are of the Python list type.
    If given a Python dictionary, this function converts the keys and values into two lists.
    """
    
    # If using gene_lists and gene_lists_names (not using gene_dict)
    if gene_dict is None:
    
        # Check to see if gene_lists is a list of gene signatures. If it isn't (meaning it's a single
        # gene signature), turn it into a list of gene signature for convenience.
        if not any([isinstance(element, list) for element in gene_lists]):
            gene_lists = [gene_lists]

        # Check to see if gene_lists_names is a list
        if not isinstance(gene_lists_names, list):
            gene_lists_names = [gene_lists_names]
            
    # Turn the dictionary into Python lists
    else:
        gene_lists_names, gene_lists = zip(*gene_dict.items())
        gene_lists_names = list(gene_lists_names)
        gene_lists = list(gene_lists)
        
    return gene_lists, gene_lists_names

###########################################
#                                         #
#     Function: Calculate Gene Scores     #
#                                         #
###########################################

def CalculateGeneScores(Adata_, gene_lists = None, gene_lists_names = None, gene_dict = None, copy = True, scale = False):
    """
    Given a set of gene lists and their names, CalculateGeneScores will calculate a score for each cell
    for each gene list.
    """
    
    # Make a copy of the Adata object
    if copy == True:
        Adata = Adata_.copy()
    else:
        Adata = Adata_
    
    # Get the gene lists and names
    gene_lists, gene_lists_names = CheckGeneLists(gene_lists = gene_lists, gene_lists_names = gene_lists_names,
                                                  gene_dict = gene_dict)
    
    # Scale if needed
    if scale == True:
        sc.pp.scale(Adata)
        
    # Loop over gene lists and gene names simultaneously
    for gene_list, gene_list_name in zip(gene_lists,gene_lists_names):
    
        # Find intersection of genes between Adata object and gene list
        gene_list = IntersectionOfGenes(Adata, gene_list)

        # Run the score_genes module
        sc.tl.score_genes(Adata, gene_list, score_name = gene_list_name, use_raw = False)
        
    return Adata

######################################
#                                    #
#     Function: Get Axes Min/Max     #
#                                    #
######################################

def GetMinMax(ax, Adata, padding = 1):
    """
    Given an Adata object, GetMinMax will adjust the given matplotlib axis so that the x-limits and y-limits
    fit to the Adata UMAP coordinates. The padding provides a little extra space around the points near the 
    edge of the plot.
    """

    # Get the x-limits and y-limits while adding padding
    xmin = min(Adata.obs['UMAP_1']) - padding
    xmax = max(Adata.obs['UMAP_1']) + padding
    ymin = min(Adata.obs['UMAP_2']) - padding
    ymax = max(Adata.obs['UMAP_2']) + padding

    # Adjust the matplotlib
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    return ax

#####################################
#                                   #
#     Function: Formalize Names     #
#                                   #
#####################################

def Formalize(name):
    """
    For computational analyses, sample names are written out in full,
    with no spaces (for example, 'CreNeg_PD1Pos'). We'd rather see this as 'Cre-PD1+' in publications.
    This function formalizes the name of the sample for the sake of plotting.
    """

    # Force the object to be a string
    name = str(name)

    # Sample names currently use the old syntax
    # We will replace these with a new syntax
    old_syntax = ['CreNeg_PD1Neg', 'CreNeg_PD1Pos', 'CrePos_PD1Neg', 'CrePos_PD1Pos', 'CrePos', 'CreNeg', 'PD1Pos', 'PD1Neg',
    'seurat_clusters', 'leiden', 'combined', 'C_I', 'C_P', 'F_I', 'F_P']
    new_syntax = ['Cre- PD1-', 'Cre- PD1+', 'Cre+ PD1-', 'Cre+ PD1+', 'Cre+', 'Cre-', 'PD1+', 'PD1-',
    'Seurat Clusters', 'Leiden', 'Combined', 'C.I.', 'C.P.', 'F.I.', 'F.P.']

    # Replace the names
    for old,new in zip(old_syntax,new_syntax):
        name = name.replace(old, new)

    return name

###################################
#                                 #
#     Function: Two-Way Stats     #
#                                 #
###################################

def TwoWayStats(df, groupby, col, groups = None, paired = False):
    """
    Given a dataframe (df), this will run a two-way statistical comparison between them and return
    a test statistic and p-value.
    
    Inputs
    ------
    df: A pandas DataFrame.
    groupby (str): A column name that labels each entry with one of two groups.
    col (str): A column name by which to run the statistical test on.
    groups (list): A list that has the two groups by which to run the test.
    paired (bool): If true, run a paired test (Wilcoxon signed-rank test).
                   Otherwise (default), run unpaired test (Wilcoxon rank-sum test).
                   
    Outputs
    -------
    stat: Test statistic.
    pval: P-value.
    diff_mean: Difference in means between two groups.
    diff_median: Difference in medians between two groups.
    """
    
    # Get the groups
    if groups is None:

        # If no groups are specified, find them
        groups = list(df[groupby].unique())

    # Make sure there are only two groups
    assert len(groups) == 2, 'There are not two groups: {}'.format(', '.join(groups))
    
    # Split dataframe
    # Each of these is a NumPy array with the column values
    group1 = df.loc[ df[groupby] == groups[0], col ].values
    group2 = df.loc[ df[groupby] == groups[1], col ].values
    
    # If running paired test,
    if paired == True:
        
        # Calculate the (paired) Wilcoxon signed-rank test
        stat, pval = scipy.stats.wilcoxon(group1, group2)
        
    else:
        
        # Calculate the (non-paired) Wilcoxon rank-sum test
        stat, pval = scipy.stats.mannwhitneyu(group1, group2)
    
    # Get the difference in means and medians
    diff_mean = group1.mean() - group2.mean()
    diff_median = np.median(group1) - np.median(group2)
    
    return stat, pval, diff_mean, diff_median


######################################################
#                                                    #
#     Function: Helper Functions for Comparisons     #
#                                                    #
######################################################

def ReversedComp(comp):
    """
    A helper function that will take a string of the format 'A vs. B' and reverse the comparison,
    giving 'B vs. A' instead.
    """
    
    # Make sure the comparison is a string value.
    assert isinstance(comp, str), 'Comparison is not a string.'
    
    # Make the reverse comparison
    reverse_comp = comp.split(' vs. ')[1] + ' vs. ' + comp.split(' vs. ')[0]
    
    return reverse_comp

def GetUniqueComps(all_groups):
    """
    Given a list (all_groups), create a new list of all possible pairwise comparisons. Order doesn't matter,
    so if there are n groups, this function will return n*(n-1)/2 comparisons.
    
    Example: If all_groups = ['A', 'B', 'C'], this function will return
    ['A vs. B', 'A vs. C', 'B vs. C']
    as output.
    
    Input
    -----
    all_groups (list): A list of all the groups.
    
    Output
    ------
    UniqueComparisons (list): A list of all possible pairwise comparisons.
    """
    
    # Initialize empty list
    UniqueComparisons = []
    
    # Loop over each group
    for i,group1 in enumerate(all_groups):
        
        # Loop over each group again
        for j,group2 in enumerate(all_groups):
            
            # Make sure that we don't include self-comparisons
            if i == j:
                continue
            
            # Make a comparison
            comp = '{} vs. {}'.format(group1, group2)
        
            # If the comparison (or opposite comparison) is already in UniqueComparisons, then don't do anything
            if comp in UniqueComparisons or ReversedComp(comp) in UniqueComparisons:
                continue
            
            # Otherwise, add the comparison
            else:
                UniqueComparisons.append(comp)
            
    return sorted(UniqueComparisons)


#####################################
#                                   #
#     Function: Mutli-Way Stats     #
#                                   #
#####################################

def MultiWayStats(df_, groupby, col, groups = None):
    """
    Given a dataframe (df), this will run all possible statistical comparisons between the groups
    (along with FDR correction) and return the p-values.
    
    Inputs
    ------
    df: A pandas DataFrame.
    groupby (str): A column name that labels each entry with one of four groups.
    col (str): A column name by which to run the statistical test on.
    groups (list): A list that has the groups by which to run the test.
    
    Outputs
    -------
    results: A pandas DataFrame that holds the p-values and q-values for each pairwise comparison.
    """
    
    from statsmodels.stats.multitest import multipletests
    
    # Make a copy of the dataframe
    df = df_.copy()
    
    # Make sure the column is string-formatted
    df[groupby] = df[groupby].astype(str)

    # Formalize names
    df[groupby] = df[groupby].apply(Formalize)
    
    # Get the groups
    if groups is None:

        # If no groups are specified, find them
        groups = list(df[groupby].unique())
    
    print('There are {} groups: {}'.format(len(groups), ', '.join(groups)))
    
    # Make dataframe to store information
    results = pd.DataFrame(index = GetUniqueComps(groups), columns = ['pval', 'qval', 'Higher In'], data = 0)
    
    # Loop over groups
    for i,group1_name in enumerate(groups):
        for j,group2_name in enumerate(groups):
            
            if j > i:
            
                # Subset data
                group1 = df.loc[ df[groupby] == group1_name, col ].values
                group2 = df.loc[ df[groupby] == group2_name, col ].values
            
                # Run statistical test
                stat, pval = scipy.stats.mannwhitneyu(group1, group2)
                
                # Add result
                res_name = '{} vs. {}'.format(group1_name, group2_name)
                results.loc[res_name, 'pval'] = pval
                
                if group1.mean() > group2.mean():
                    results.loc[res_name, 'Higher In'] = group1_name
                else:
                    results.loc[res_name, 'Higher In'] = group2_name

    # Run FDR correction
    _, qvals, _, _ = multipletests(results['pval'].values, method = 'fdr_bh')
    results['qval'] = qvals

    # Add information stating whether or not the result is significant
    results['Sig'] = ['Yes' if q < 0.05 else 'No' for q in results['qval']]

    # For non-significant results, change the 'Higher In' column to '--'
    results.loc[ results['qval'] > 0.05, 'Higher In'] = '--'
    
    return results[['pval', 'qval', 'Sig', 'Higher In']]


###########################################
#                                         #
#     Function: Make Heatmap of Genes     #
#                                         #
###########################################

def MakeHeatmapGenes(Adata, genes, groupby = 'seurat_clusters', row_cluster = False, col_cluster = False, linewidth = 0.1,
                     figsize = (10,8), vmin = None, vmax = None, cmap = 'RdBu_r', zscore = True, rotate = False,
                     order = None, title = None,
                     save_name = None):
    """
    Given an Adata object and a list of genes, make a heatmap showing the expression of each gene
    in a given grouping.
    
    Inputs
    ------
    Adata: AnnData object.
    genes: A Python list of genes.
    groupby (str): The name of the column in the Adata.obs metadata to group the data by.
    row_cluster (bool): Whether or not to put a dendrogram along the rows (genes).
    col_cluster (bool): Whether or not to put a dendrogram along the columns (groups).
    figsize (tuple): The size of the figure in (width, height).
    vmin/vmax (float): Cap the mininmum/maximum values of the colormap.
    cmap: The colormap to use
    zscore (bool): Whether or not to take Z-scores for each gene.
    rotate (bool): Whether or not to rotate the heatmap.
    order (list): The names of the categories across the x-axis.
    title (str): The title of the figure.
    save_name (str): If given, save the heatmap to this file.
    
    Outputs
    -------
    None
    """
    
    # Get all the genes that are in the Adata object
    genes = IntersectionOfGenes(genes, Adata)

    # Make sure 'seurat_clusters' is in the metadata - otherwise, use Leiden clustering
    if groupby == 'seurat_clusters':
        if groupby not in Adata.obs.columns:
            groupby = 'leiden'
        
    # Create a dataframe for the heatmap and add the metadata
    df = pd.DataFrame(index = Adata.obs_names, columns = genes, data = Adata[:, genes].X)
    df[groupby] = list(Adata.obs[groupby])
    
    # Take the average gene expression for each group, and then the Z-score for visualization
    df = df.groupby(groupby).mean().T
    if zscore == True:
        df = pd.DataFrame(index = df.index, columns = df.columns, data = scipy.stats.zscore(df, axis = 1))
        center = 0
    else:
        center = None
        
    # Generate a title
    if title is None:
        title = 'Average Gene Expression Across {}'.format(Formalize(groupby))
        if zscore == True:
            title += ' (Z-Scores)'
        else:
            title += ' (Avg Counts)'
            
    # Rotate the heatmap
    if rotate == True:
        df = df.T
        
    # Re-order the heatmap
    if order is not None:
        df = df[order]
    
    # Make the heatmap, set the title, axes labels, and rotate the gene labels
    g = sns.clustermap(df, row_cluster = row_cluster, col_cluster = col_cluster, center = center, cmap = cmap,
                       linewidth = linewidth, linecolor = 'gray', vmin = vmin, vmax = vmax, figsize = figsize,
                       dendrogram_ratio = 0.15)
    g.fig.suptitle(title)
    if rotate == True:
        g.ax_heatmap.set_ylabel(Formalize(groupby))
        g.ax_heatmap.set_xlabel('Genes')
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation = 0)
    else:
        g.ax_heatmap.set_xlabel(Formalize(groupby))
        g.ax_heatmap.set_ylabel('Genes')
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation = 0)
    
    # Save the figure
    if save_name is not None:
        g.savefig(save_name)
    plt.show()
    return

########################################
#                                      #
#     Function: Plot Genes on UMAP     #
#                                      #
########################################

def PlotGenesUMAP(Adata, genes, ncols = 3, cmap = None, vmin = None, vmax = None, save_name = None):
    """
    PlotGenesUMAP will take a given gene list (genes) and plot each on a UMAP. This is essentially a wrapper around
    Scanpy's plotting function.

    Inputs
    -----
    Adata: AnnData object.
    genes: A Python list that holds the genes.
    ncols (int): The number of columns to plot in the figure.
    cmap: The colormap to use for the heatmap.
    vmin/vmax (float): Cap the mininmum/maximum values of the colormap.
    save_name: A file name of where to save the plot, if specified.

    Outputs
    -------
    None
    """
    
    # Make sure 'genes' is a list
    if type(genes) != list:
            genes = [genes]    

    # Get all the genes in the gene list that are also in the Adata object
    Genes = IntersectionOfGenes(genes, Adata)

    # If there are no genes, check the metadata
    if len(Genes) == 0:
            Genes = IntersectionOfGenes(genes, Adata.obs.columns)

    # Make sure there are genes in the list
    assert len(Genes) > 0, "None of the genes are in the Adata object."    
 
    # Make sure that the UMAP coordinates are in the Adata.obsm slot.
    err_msg = "UMAP coordinates are not found. Make sure they are placed as a NumPy array in Adata.obsm['X_umap']."
    assert 'X_umap' in list(Adata.obsm.keys()), err_msg
    
    # Plot the genes
    fig = sc.pl.umap(Adata, color = Genes, ncols = ncols, cmap = cmap, vmin = vmin, vmax = vmax,
                     show = False, return_fig = True)
    
    # Save the figure
    if save_name is not None:
        plt.savefig(save_name)
    plt.show()
    return

#####################################################
#                                                   #
#     Function: Make Heatmap of Gene Signatures     #
#                                                   #
#####################################################

def MakeHeatmapSignatures(Adata_, gene_lists = None, gene_lists_names = None, gene_dict = None, scale = False,
                          groupby = 'seurat_clusters', row_cluster = False, col_cluster = False,
                          figsize = (10,8), cmap = 'RdBu_r', vmin = None, vmax = None, order = None,
                          title = None, save_name = None):
    """
    MakeHeatmapSignatures works on a list of gene signatures. For each gene signature/gene list,
    it calculates the gene scores for all the cells in a given AnnData object, and
    displays them on a heatmap.
    This function will accept one of either:
        1.) A gene_lists input and an associated gene_lists_names input
        2.) A gene_dict object
    
    Inputs
    ------
    Adata: AnnData object.
    gene_lists: A list of gene signatures to calculate gene scores.
    gene_lists_names: A list of the names of the gene signatures (ex: Activation, Naive, etc).
    gene_dict: A Python dictionary that has keys as gene list names and values as the gene lists.
    scale: Whether or not to scale the data so that each gene has unit variance and zero mean (default False).
    groupby (str): The name of the column in the Adata.obs metadata to group the data by.
    row_cluster (bool): Whether or not to put a dendrogram along the rows (genes).
    col_cluster (bool): Whether or not to put a dendrogram along the columns (groups).
    figsize (tuple): The size of the figure in (width, height).
    cmap: The colormap to use for the heatmap.
    vmin/vmax (float): Cap the mininmum/maximum values of the colormap.
    title (str): The title of the figure.
    save_name: A file name of where to save the plot, if specified.
    
    Outputs
    -------
    None
    """
    
    # Suppress warning
    pd.options.mode.chained_assignment = None
    
    # Copy the Adata object, as to not change the global object
    Adata = Adata_.copy()

    # Make sure 'seurat_clusters' is in the metadata - otherwise, use Leiden clustering
    if groupby == 'seurat_clusters':
        if groupby not in Adata.obs.columns:
            groupby = 'leiden'
    
    # Generate a title
    if title is None:
        title = 'Average Score for Gene Signatures by {} (Z-Scores)'.format(Formalize(groupby))
    
    # Check the gene lists
    gene_lists, gene_lists_names = CheckGeneLists(gene_lists = gene_lists, gene_lists_names = gene_lists_names,
                                                  gene_dict = gene_dict)
    
    # Calculate the gene scores for each gene list
    Adata = CalculateGeneScores(Adata, gene_lists = gene_lists, gene_lists_names = gene_lists_names, copy = False, scale = scale)
        
    # Make a heatmap
    # Get the cell scores for each gene list as a dataframe, and add the grouping 
    df = Adata.obs[gene_lists_names]
    df.loc[:, groupby] = Adata.obs[groupby]
    
    # Find the average score for gene list for each group, and then the Z-score for visualization
    df = df.groupby(groupby).mean().T
    df = pd.DataFrame(index = df.index, columns = df.columns, data = scipy.stats.zscore(df, axis=1))
    
    # Re-order the dataframe
    if order is not None:
        df = df[order]
    
    # Make the heatmap, set the title and axes labels, and rotate the signature labels
    g = sns.clustermap(df, row_cluster = row_cluster, col_cluster = col_cluster, center = 0, cmap = cmap,
                       linewidth = 1, linecolor = 'gray', vmin = vmin, vmax = vmax,
                       figsize = figsize, dendrogram_ratio = 0.15)
    g.fig.suptitle(title)
    g.ax_heatmap.set_xlabel(Formalize(groupby))
    g.ax_heatmap.set_ylabel('Gene Signatures')
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation = 0)

    # Save the figure
    if save_name is not None:
        g.savefig(save_name)
    plt.show()
    return

#################################################
#                                               #
#     Function: Plot Gene Signatures on UMAP    #
#                                               #
#################################################

def GeneSignatureUMAP(Adata_, gene_lists = None, gene_lists_names = None, gene_dict = None, scale = False,
                      cmap = 'RdBu_r', vmin = None, vmax = None, ncols = 3, save_name = None):
    """
    GeneSignatureUMAP works on a list of gene signatures. For each gene signature/gene list,
    it calculates the gene scores for all the cells in a given AnnData object, and
    displays them on a UMAP.
    This function will accept one of either:
        1.) A gene_lists input and an associated gene_lists_names input
        2.) A gene_dict object
    
    Inputs
    ------
    Adata: AnnData object.
    gene_lists: A list of gene signatures to calculate gene scores.
    gene_lists_names: A list of the names of the gene signatures (ex: Activation, Naive, etc)
    gene_dict: A Python dictionary that has keys as gene list names and values as the gene lists.
    scale: Whether or not to scale the data so that each gene has unit variance and zero mean (default False).
    cmap: The colormap to use for the heatmap.
    vmin/vmax (float): Cap the mininmum/maximum values of the colormap.
    ncols (int): The number of columns to plot in the figure.
    save_name: A file name of where to save the plot, if specified.
    
    Outputs
    -------
    None
    """
    
    # Copy the Adata object, as to not change the global object
    Adata = Adata_.copy()
    
    # Check the inputs
    gene_lists, gene_lists_names = CheckGeneLists(gene_lists = gene_lists, gene_lists_names = gene_lists_names,
                                                  gene_dict = gene_dict)
    
    # Calculate the gene scores for each gene list
    Adata = CalculateGeneScores(Adata, gene_lists = gene_lists, gene_lists_names = gene_lists_names, copy = False, scale = scale)
    
    # Plot the results
    fig = sc.pl.umap(Adata, color = gene_lists_names, cmap = cmap, ncols = ncols, vmin = vmin, vmax = vmax,
                     show = False, return_fig = False)
    
    # Save the figure
    if save_name is not None:
        plt.savefig(save_name)
    plt.show()
    return

#########################################################
#                                                       #
#     Function: Plot Gene Signatures on Violin Plot     #
#                                                       #
#########################################################

def MakeViolinSignatures(Adata_, gene_lists = None, gene_lists_names = None, gene_dict = None, scale = False,
                         groupby = 'sample', Cap = False, order = None, palette = None, ncols = 3, Title = None,
                         figsize=(10,8), save_name = None):
    """
    ViolinPlotGeneLists calculates the gene scores for all the cells in a given AnnData object, and
    displays them on a violin plot, separated by the 'groupby' variable.
    This function will accept one of either:
        1.) A gene_lists input and an associated gene_lists_names input
        2.) A gene_dict object
    
    Inputs
    ------
    Adata: AnnData object.
    gene_lists: A list of gene signatures to calculate gene scores.
    gene_lists_names: A list of the names of the gene signatures (ex: Activation, Naive, etc)
    gene_dict: A Python dictionary that has keys as gene list names and values as the gene lists.
    scale: Whether or not to scale the data so that each gene has unit variance and zero mean (default False).
    groupby (str): The name of the column in the Adata.obs metadata to group the data by.
    Cap (bool): Whether or not to cap the signature scores before plotting.
    order (list): A list of groups in the order of displaying the violin.
    palette (list): A list of colors for each violin.
    ncols (int): The number of columns to make in the figure.
    figsize (tuple): The size of the figure in (width, height).
    save_name: A file name of where to save the plot, if specified.
    
    Outputs
    -------
    None
    """
    
    # Copy the Adata object, as to not change the global object
    Adata = Adata_.copy()
    
    # Make sure 'seurat_clusters' is in the metadata - otherwise, use Leiden clustering
    if groupby == 'seurat_clusters':
        if groupby not in Adata.obs.columns:
            groupby = 'leiden'

    # Create a sample palette if none is specified
    if palette is None:
        from matplotlib import colors
        palette = sns.xkcd_palette(['scarlet', 'black', 'ultramarine'])
        palette.append(colors.to_rgb('green'))
            
    # Check the gene lists
    gene_lists, gene_lists_names = CheckGeneLists(gene_lists = gene_lists, gene_lists_names = gene_lists_names,
                                                  gene_dict = gene_dict)
            
    # Make figure and axes
    # Number of axes/subplots is equal to number of gene lists
    # Number of columns is default to 3, unless there are fewer than 3 gene lists
    n_axes = len(gene_lists)
    ncols = min(ncols, n_axes)
    fig, axes = MakeFigAxes(ncols = ncols, n_axes = n_axes, figsize = figsize)
    
    # Calculate the gene scores for each gene list
    Adata = CalculateGeneScores(Adata, gene_lists = gene_lists, gene_lists_names = gene_lists_names, copy = False, scale = scale)
    
    # Formalize names
    Adata.obs[groupby] = Adata.obs[groupby].astype(str).apply(Formalize)
        
    # Loop over gene lists and gene names simultaneously
    for i,(gene_list,gene_list_name) in enumerate(zip(gene_lists,gene_lists_names)):
    
        # Get the ith axis
        ax = axes.flatten()[i]

        # Cap gene scores
        if Cap == True:
            Adata.obs[gene_list_name] = CapScores(Adata.obs[gene_list_name])

        # Make a violin plot
        sns.violinplot(data = Adata.obs, x = groupby, y = gene_list_name, order = order, scale = 'width',
                       ax = ax, palette = palette)

        # Set the y-label and initialize the title
        ax.set_ylabel('Gene Score')
        gene_list = IntersectionOfGenes(gene_list, Adata)
        if Title is None:
            title = '{} ({} cells, {} genes)'.format(gene_list_name, len(Adata), len(gene_list))
        else:
            title = Title

        # Add two-way statistics if there are only 2 groups
        groups = list(Adata.obs[groupby].unique())
        if len(groups) == 2:

            # Calculate the statistics
            _, pval, diff_mean, _ = TwoWayStats(Adata.obs, groupby = groupby, col = gene_list_name, groups = [groups[0], groups[1]])

            # Add the p-value to the title of the axies
            title += '\np-val = {:.3g}'.format(pval)
            
            # For significant differences, add the corresponding info
            if pval < 0.05:
                if diff_mean > 0:
                    title += '\nHigher in {} by {:.2f}'.format(groups[0], diff_mean)
                else:
                    title += '\nHigher in {} by {:.2f}'.format(groups[1], -1*diff_mean)
                    
            # Otherwise, label as 'No difference'
            else:
                title += '\nNo difference'

        # Set the title
        ax.set_title(title)
    
    # Save the figure
    plt.tight_layout()
    if save_name is not None:
        plt.savefig(save_name)
    plt.show()
    return

############################################
#                                          #
#     Function: Make Stacked Bar Plots     #
#                                          #
############################################

def StackedBarPlot(Adata, palette = None, xaxis = 'seurat_clusters', yaxis = 'combined', order = None,
                   orientation = 'horizontal', title = None, save_name = None):
    """
    Using the metadata stored in the Adata.obs object, create a stacked bar plot.
    
    Inputs
    ------
    Adata: AnnData object to use.
    palette (list): List of colors to use.
    xaxis (str): Name of metadata column used to make bars across the x-axis.
    yaxis (str): Name of metadata column used to make bars across the y-axis.
    order (list): Order to plot the groups in (x-axis if orientation='vertical', y=axis if orientation='horizontal').
    orientation (str): Whether the bars should be horizontal or vertical.
    title (str): Title of the plot.
    save_name (str): The name of the file to save on disk.
    
    Outputs
    -------
    None
    """
    
    # Make sure the orientation is either horizontal or vertical.
    assert orientation in ['horizontal', 'vertical'], "Orientation must be either 'horizontal' or 'vertical.'"
    
    # Assign a color palette
    if palette is None:
        palette = ['black','tan','green','purple','orange','red','magenta','blue','pink','gray','brown',
                 'yellow','lightblue','lightgreen','maroon','aqua','lime','yellow']
    
    # Copy the metadata
    df = Adata.obs.copy()
    
    # Create a cross-tabulation dataframe
    df = pd.crosstab(df[yaxis], df[xaxis])
    
    # Assume horizontal bars
    if orientation == 'horizontal':
        bar = 'barh'
        if title is None:
            title = 'Percentage of {} per {}'.format(Formalize(xaxis), Formalize(yaxis))
    
    # Vertical bars
    if orientation == 'vertical':
        
        # Transpose the dataframe
        df = df.T
        bar = 'bar'
        if title is None:
            title = 'Percentage of {} per {}'.format(Formalize(yaxis), Formalize(xaxis))
    
    # Normalize each row so that it sums to 1
    df = df.div(df.sum(axis=1), axis=0)

    # Re-order dataframe
    if order is not None:
        if orientation == 'vertical':
            df.index = list(df.index)
            df = df.loc[order]
        else:
            df.columns = list(df.columns)
            df = df[order]
            
    # Formalize the names of the groups
    df.columns = [Formalize(x) for x in df.columns]
    df.index = [Formalize(x) for x in df.index] 
    
    # Make a stacked bar plot
    ax = df.plot(kind = bar, stacked = True, color = palette)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))

    # Make the title and label the axes
    ax.set_title(title)
    
    # Assume horizontal bars
    if orientation == 'horizontal':
        ax.set_ylabel(Formalize(yaxis))
        ax.set_xlabel('Percent of Cells')
        
    # Vertical bars
    if orientation == 'vertical':
        ax.set_xlabel(Formalize(xaxis))
        ax.set_ylabel('Percent of Cells')
        ax.set_xticklabels(ax.get_xticklabels(), rotation = 0)
        
    if save_name is not None:
        plt.savefig(save_name, bbox_extra_artists = (lgd,), bbox_inches = 'tight')
    plt.show()
    return

