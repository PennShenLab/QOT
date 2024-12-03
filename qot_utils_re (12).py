import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture
import ot
import importlib
import scanpy as sc
import importlib
import scanpy as sc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve, auc
import phate
import warnings
warnings.filterwarnings("ignore")
import scanpy as sc
import time
import importlib
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve, auc
import phate
import warnings
warnings.filterwarnings("ignore")
import scanpy as sc
import time
import scipy.stats as sps
import scipy.linalg as spl
import seaborn as sns
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import precision_recall_curve, auc
import numpy as np
import phate as ph

def cell_type_diff_two_sub_patient_groups(proportions: pd.DataFrame = None,
                                          cell_types: list = None,
                                          labels:str = 'Predicted_Labels',
                                          group1: str = 'Tumor 1',
                                          group2: str = 'Tumor 2',
                                          pval_thr: float = 0.05,
                                          figsize: tuple = (15, 4),
                                          file_path: str = None,fontsize:int = 19):
    """


    Parameters
    ----------
    proportions : pd.DataFrame, optional
        cell types proportion in each samples. The default is None.
    cell_types : list, optional
        number of cell types to be considered. The default is None.
    labels : str, optional
        name of the column containing the labels of patient sub-groups. The default is 'Predicted_Labels'.
    group1 : str, optional
        Name of the first patients sub-group to check the differentiations. The default is 'Tumor 1'.
    group2 : str, optional
        Name of the second patients sub-group to check the differentiations. The default is 'Tumor 2'.
    pval_thr : float, optional
        P-value threshold. The default is 0.05.
    figsize : tuple, optional
        Determine the figure size. The default is (15, 4).
    file_path : str, optional
        Determine the path to store the figure. The default is None.
    fontsize:int, the size of the font
    Returns
    -------
    None.
    Plot the statistical scores showing the how significant each cell type
    differentiated in each patients sub-group (filter based on p-value threshold)
    Save table of statisical information of all cell types,
    sorted by their p-value and statistic score

    """

    group1_proportions = proportions[proportions['Predicted_Labels'] == group1]
    group2_proportions = proportions[proportions['Predicted_Labels'] == group2]

    stats_group12 = []
    scores_group12 = []
    # using Welch’s t-test to detect the statistically significant difference
    ## between group 1 and group 2
    for cell_type in cell_types:
        ttest_result = ttest_ind(group1_proportions[cell_type],
                                 group2_proportions[cell_type], equal_var = False)
        scores_group12.append(ttest_result[0])
        stats_group12.append(ttest_result[1])

    # adjust the p-values
    stats_group12 = multipletests(list(stats_group12), method='fdr_bh')[1]

    stats_bar_group12 = pd.DataFrame()
    stats_bar_group12['cell_type'] = cell_types
    stats_bar_group12['adjPval'] = stats_group12
    stats_bar_group12['-logPval'] = -np.log(stats_group12)
    stats_bar_group12['score'] = scores_group12

    # sort the data based on p-value and statistic score
    stats_bar_group12 = stats_bar_group12.sort_values(by = ['score', '-logPval'],
                                                      ascending = [False, False])

    # save the table of statistical differentiation of two groups
    ## as positive score shows cell types statistically significant on the first group
    ## and negative scrore shows cell types statistically significant on the second group



    # filter data based on a p-value threshold
    stats_bar_group12 = stats_bar_group12[stats_bar_group12['adjPval'] < pval_thr]


    fig, ax = plt.subplots(figsize=figsize)

# Call the modified plotting function
    plot_hor_vs_vert(ax, stats_bar_group12, x='score', y='cell_type', c='type',
                 xlabel='statistic score', ylabel=None,
                 rotation=None, tick_bottom=True, tick_left=False,
                 title=  group1 + " vs " + group2, group1=group1, group2=group2, fontsize=fontsize)

    fig.tight_layout()
  
def plot_hor_vs_vert(ax, data, x, y, c, xlabel, ylabel, rotation, tick_bottom, tick_left, title, group1, group2, color1='tab:blue', color2='tab:red', fontsize=24):
    '''
    Plot horizontal bar charts using Seaborn with legends for score significance.

    Parameters:
        ax : matplotlib.axes.Axes
            The axes on which to plot.
        data : DataFrame
            The data for plotting.
        x : str
            The column name for x values.
        y : str
            The column name for y values.
        c : str
            The column name for color coding based on conditions.
        xlabel : str
            The label for the x-axis.
        ylabel : str
            The label for the y-axis.
        rotation : int
            The rotation angle of the x-tick labels.
        tick_bottom : bool
            Whether to display bottom ticks.
        tick_left : bool
            Whether to display left ticks.
        title : str
            The title of the plot.
        group1 : str
            Name of the first group (used for legend and color mapping).
        group2 : str
            Name of the second group (used for legend and color mapping).
        color1 : str, optional
            Color associated with group1 (default is 'tab:blue').
        color2 : str, optional
            Color associated with group2 (default is 'tab:red').
        fontsize : int
            The font size for labels and titles.
    '''
    # Determine colors based on the score values
    cols = [color1 if val >= 0 else color2 for val in data[x]]
    bar_plot = sns.barplot(x=x, y=y, data=data, ci=None, palette=cols, ax=ax)

    # Set titles and labels with specified font settings
    ax.set_title(title, fontsize=fontsize, fontweight='bold')
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=fontsize, rotation=rotation)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=fontsize)

    # Customize seaborn style
    sns.despine(ax=ax, bottom=False, left=True)
    ax.grid(False)
    ax.tick_params(bottom=tick_bottom, left=tick_left)

    # Create a custom legend for the color coding
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color1, edgecolor=color1, label=group1),
                       Patch(facecolor=color2, edgecolor=color2, label=group2)]
    ax.legend(handles=legend_elements, fontsize=fontsize)

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


from sklearn import metrics
import os
from sklearn.neighbors import NearestCentroid
from sklearn.metrics.cluster import rand_score
def select_best_sil(adata,resolutions=[],marker='o',figsize=(6,6),facecolor="white",metric='cosine',path=None,start=0.2,step=0.1,end=2):
    """
    Parameters
    ----------
    adata : adata,
    EMD : W distance,
    path_to_results:str, path to save the plot
    resolutions: list,a list of your desire resulotions
    marker : str, optional (default='o')
            Marker style for data points in the plot.
    figsize : tuple, optional
        Figure size (width, height) in inches. Default is (10, 10).
    facecolor : str, optional
        Background color of the saved plot image. Default is 'white'.

    metric: str, metric for leiden clustering and calculating sil. 

    Returns
    -------
    None,
    plot Silhouette Score vs. Resolution to figure out the best sil.
    plot Silhouette Score vs. Number of clusters.
    """
    
    resolutions = [start + step * i for i in range(int((end - start) / step) + 1)]
    # Create a list to store the Silhouette Scores
    best_res=start
    best_sil=0

    sil_scores = []
    number_of_clusters=[]
    EMD=adata.uns['QOT_Distance']
    # Define the range of resolutions/number of clusters to test
    # Loop through the resolutions and calculate the Silhouette Score for each
    for resolution in resolutions:
        # Cluster the data using Louvain clustering at the specified resolution
        adata_emd = sc.AnnData(EMD)
        sc.pp.neighbors(adata_emd, metric=metric)
        sc.tl.leiden(adata_emd, resolution = resolution)
        # Calculate the Silhouette Score
        predicted_labels = np.array(adata_emd.obs.leiden)
        
        Silhouette = metrics.silhouette_score(EMD, predicted_labels, metric =metric)
        if float(Silhouette) > best_sil:
            best_sil=Silhouette
            best_res=resolution
            

        # Append the Silhouette Score to the list
        sil_scores.append(Silhouette)
        number_of_clusters.append(len(np.unique(predicted_labels)))
    adata.uns['best_res']=best_res
    # Plot the Silhouette Scores against the resolutions
    plt.figure(figsize = figsize,facecolor=facecolor)
    plt.figure(figsize = figsize,facecolor=facecolor)
    plt.plot(resolutions, sil_scores, marker=marker)
    plt.xlabel('Resolution')
    plt.ylabel('Silhouette Score')
    plt.title('Silhouette Score vs. Resolution')
    plt.show()
    
    
    # Plot the Silhouette Scores against the resolutions
    plt.figure(figsize = figsize,facecolor=facecolor)
    plt.plot(resolutions, number_of_clusters, marker=marker)
    plt.xlabel('Resolution')
    plt.ylabel('Number of Clusters')
    plt.title('Number of Clusters vs. Resolution')
    plt.show()
def Sil_computing(EMD, real_labels, metric='cosine'):
    """
    Compute the Silhouette score based on Wasserstein distances.

    Parameters
    ----------
    EMD : numpy.ndarray
        Wasserstein distances matrix.
    real_labels : list or numpy.ndarray
        True labels or ground truth.
    metric : str, optional
        Metric for calculating pairwise distances, by default 'cosine'.

    Returns
    -------
    float
        Silhouette score indicating cluster quality.
    """
    Silhouette = silhouette_score(EMD, real_labels, metric =metric)
    #print("Silhouette score: ", Silhouette)
    return Silhouette
def clustering_emd(adata,res=0.3,metric='cosine',groupby_col='Leiden',swap_axes=False,cmap="Blues_r",dendrogram=True,show_gene_labels=True,var_group_rotation=45,figsize=[6,6],save=False,sorter_leiden=None):
    
    """
    Perform clustering and visualization of EMD (Earth Mover's Distance) data in AnnData object.

    Parameters:
    adata (AnnData): Input AnnData object containing EMD data.
    res (float): Resolution parameter for Leiden clustering. Default is 0.3.
    metric (str): Distance metric for clustering. Default is 'cosine'.
    groupby_col (str): Grouping variable for plotting. 'Leiden' groups by predicted clusters, 'status' groups by real labels. Default is 'Leiden'.
    swap_axes (bool): Swap the axes in the heatmap. Default is False.
    cmap (str): Colormap for the heatmap. Default is "Blues_r".
    dendrogram (bool): Display dendrograms in the heatmap. Default is True.
    show_gene_labels (bool): Show gene labels in the heatmap. Default is True.
    var_group_rotation (int): Rotation angle for gene labels. Default is 45 degrees.
    figsize (list): Size of the heatmap figure. Default is [12, 12].
    save (bool): Save the heatmap figure. Default is False.
    sorter_leiden (list or None): Custom order for Leiden clusters. If not provided, the default order is used.

    Returns:
    proportion_df (DataFrame): DataFrame containing proportions of sub-clusters in each sample.
    """
    
    EMD=adata.uns['QOT_Distance']
    proportions=adata.uns['proportions']
    df = adata.uns['Datafame_for_use']
    sources = df['sampleID'].unique()
    annot=adata.uns['annot']
    real_labels=adata.uns['real_labels']
    
    EMD_df=pd.DataFrame(EMD,columns=sources)
    EMD_df['sampleID']=sources
    EMD_df['status']=list(real_labels)
    
    
    adata_emd = sc.AnnData(EMD)
    sc.pp.neighbors(adata_emd, metric=metric)
    sc.tl.leiden(adata_emd, resolution = res)
    predicted_labels = np.array(adata_emd.obs.leiden)
    Silhouette = Sil_computing(EMD/EMD.max(), predicted_labels,metric=metric)
    
    proportion_df=pd.DataFrame(proportions)
    proportion_df=proportion_df.T
    proportion_df.columns=annot.cell_type.unique()
    
    # proportion_df['Predicted_Labels']=predicted_labels
    # proportion_df['sampleID']=list(sources)
    
    EMD_df['Leiden']=predicted_labels

    # print(EMD_df)
    proportion_df['sampleID'] = proportion_df.index
    proportion_df = pd.merge(EMD_df[['sampleID', 'Leiden']], proportion_df, on='sampleID', how='inner')

    if groupby_col=='status':
        sorter=np.unique(real_labels)
        EMD_df['status'] = EMD_df.status.astype("category")
        EMD_df['status'] = EMD_df['status'].cat.set_categories(sorter)
        EMD_df=EMD_df.sort_values(["status"])
    elif groupby_col=='Leiden':
        if sorter_leiden==None:
            sorter=EMD_df.Leiden.unique()
        else:
            sorter=sorter_leiden
        EMD_df['Leiden'] = EMD_df.Leiden.astype("category")
        EMD_df['Leiden'] = EMD_df['Leiden'].cat.set_categories(sorter)
        EMD_df=EMD_df.sort_values(["Leiden"])
    obs = pd.DataFrame()
    obs['sampleID']=EMD_df.sampleID.astype(str)
    obs['status']=EMD_df.status.astype(str)
    obs['Leiden']=EMD_df.Leiden.astype(str)
    proportion_df = proportion_df.rename(columns={'Leiden': 'Predicted_Labels'})
    df_genes = pd.DataFrame(index = EMD_df.columns[0:EMD_df.shape[0]])
    adata_emd = sc.AnnData(X = EMD_df[ EMD_df.columns[0:EMD_df.shape[0]]].values, var =df_genes, obs = obs )
    sc.pl.heatmap(adata_emd,adata_emd.obs.sampleID,groupby=[groupby_col],swap_axes=swap_axes,cmap=cmap,dendrogram=dendrogram,show_gene_labels=show_gene_labels,var_group_rotation=var_group_rotation,figsize=figsize,save=save)
    return proportion_df
def Extract_Info(adata, gene_matrix='X_PCA', type_cell='cell_types',id='sampleID',progession='status',dataset_type = 'pathomic'):
    """
    Create a combined DataFrame from PCA results and annotations in adata_T.

    Parameters:
    adata_T (Anndata object): An object containing PCA results and annotations.

    Returns:
    DataFrame: A combined DataFrame with PCA results and annotations.
    """

    # Convert PCA results to a DataFrame and reset index
    if dataset_type =='pathomic':
        pca_results_df = pd.DataFrame(adata.X).reset_index(drop=True)
    else:
        pca_results_df = pd.DataFrame(adata.obsm[gene_matrix]).reset_index(drop=True)
    # Reset the index for 'sampleID' and 'cell_subtype'
    sample_ids = adata.obs[id].reset_index(drop=True)
    cell_subtypes = adata.obs[type_cell].reset_index(drop=True)
    status = adata.obs[progession].reset_index(drop=True)

    # Concatenate the PCA results with 'sampleID', 'cell_subtype', and 'status'
    combined_pca_df = pd.concat([pca_results_df, sample_ids, cell_subtypes, status], axis=1)

    # Rename columns
    combined_pca_df = combined_pca_df.rename(columns={type_cell: 'Cell_type',
                                                      id: 'sampleID',
                                                      progession: 'status'})

    adata.uns['Datafame_for_use'] = combined_pca_df
    return adata



def Gaussian_Mixture_Representation(adata, num_components=5, random_state=2, min_samples_gmm=0):
    """
    Process the given DataFrame using Gaussian Mixture Models (GMM).

    Parameters:
    df (DataFrame): The input data frame.
    num_components (int): The number of components for GMM.
    random_state (int): The seed for the random number generator.
    min_samples_for_gmm (int): Minimum number of samples required for fitting GMM.

    Returns:
    dict: A dictionary containing GMM parameters for each (source, label) pair.
    """
    df = adata.uns['Datafame_for_use']
    # Uncomment if you need to filter or drop columns initially
    df = df[df['Cell_type'] != 'Unknown']
    df = df.drop(['status'], axis=1)

    grouped = df.groupby(['sampleID', 'Cell_type'])
    gmm_count_per_source = {}
    params = {}

    for (source, label), group in grouped:
        data = group.drop(['sampleID', 'Cell_type'], axis=1)
        num_samples = data.shape[0]

        if num_samples >= num_components + min_samples_gmm:
            gmm = GaussianMixture(n_components=num_components, random_state=random_state).fit(data)

            key = (source, label)
            params[key] = {
                'means': gmm.means_,
                'covariances': gmm.covariances_,
                'weights': gmm.weights_,
                'proportion': len(group) / len(df[df['sampleID'] == source])
            }

        elif num_samples == 1:  # Exactly 1 samples, use itself
            mean_sample = data.iloc[0].values.reshape(1, -1)
            key = (source, label)
            params[key] = {
            'means': mean_sample,
            'covariances': np.zeros((1, data.shape[1], data.shape[1])),
            'weights': np.array([1]),
            'proportion': len(group) / len(df[df['sampleID'] == source])
            }
        elif num_samples == 2:  # Exactly 2 samples, use mean
            mean_sample = data.mean(axis=0).values.reshape(1, -1)
            key = (source, label)
            params[key] = {
            'means': mean_sample,
            'covariances': np.zeros((1, data.shape[1], data.shape[1])),
            'weights': np.array([1]),
            'proportion': len(group) / len(df[df['sampleID'] == source])
            }

    adata.uns['GMM_Representation'] = params
    return adata


def euclidean_distance(vec1, vec2):
    """Calculate the Euclidean distance between two vectors."""
    return np.linalg.norm(vec1 - vec2)

def cosine_similarity(vec1, vec2):
    """Calculate the cosine similarity between two vectors."""
    return np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))


def GaussianW2(m0,m1,Sigma0,Sigma1):
    # compute the quadratic Wasserstein distance between two Gaussians with means m0 and m1 and covariances Sigma0 and Sigma1
    Sigma00  = spl.sqrtm(Sigma0)
    Sigma010 = spl.sqrtm(Sigma00@Sigma1@Sigma00)
    d        = np.linalg.norm(m0-m1)**2+np.trace(Sigma0+Sigma1-2*Sigma010)
    return d



def calculate_qot(adata, method = "alternative",normalized = False):
    """
    Calculate the GW2 distance between sources.

    Parameters:
    sources (list): List of source names.
    params (dict): Dictionary containing GMM parameters.
    df (pandas.DataFrame): DataFrame containing cell types.
    GW2 (function): Function to calculate GW2 distance.

    Returns:
    np.ndarray: Distance matrix containing GW2 distances between sources.
    """
    params = adata.uns['GMM_Representation']
    df = adata.uns['Datafame_for_use']
    sources = df['sampleID'].unique()

    num_sources = len(sources)
    total_calculations = (num_sources * (num_sources - 1)) // 2
    distance_matrix = np.zeros((num_sources, num_sources))
    calculation_count = 0
    for i, source1 in enumerate(sources):
        for j, source2 in enumerate(sources):
            if i < j:
                calculation_count += 1
                # print(f"Calculating GW2 distance {calculation_count}/{total_calculations}: {source1} to {source2}")

                # print("celltype_lens  ", len(df['Cell_type'].unique()))
                arrays_to_concat1 = [params[(source1, cell_type)]['means']
                                     for cell_type in df['Cell_type'].unique()
                                     if (source1, cell_type) in params and params[(source1, cell_type)]['means'].size > 0]

                arrays_to_concat2 = [params[(source2, cell_type)]['means']
                                     for cell_type in df['Cell_type'].unique()
                                     if (source2, cell_type) in params and params[(source2, cell_type)]['means'].size > 0]

                # Check if either list is empty after filtering
                if not arrays_to_concat1 or not arrays_to_concat2:
                    print(f"Skipping distance calculation for {source1} and {source2} due to missing data.")
                    continue


                means1 = np.concatenate(arrays_to_concat1) if arrays_to_concat1 else np.array([])
                means2 = np.concatenate(arrays_to_concat2) if arrays_to_concat2 else np.array([])



                covs1 = np.concatenate([params[(source1, cell_type)]['covariances']
                                        for cell_type in df['Cell_type'].unique()
                                        if (source1, cell_type) in params])

                covs2 = np.concatenate([params[(source2, cell_type)]['covariances']
                                        for cell_type in df['Cell_type'].unique()
                                        if (source2, cell_type) in params])
                weights1 = np.concatenate([params[(source1, cell_type)]['weights'] * params[(source1, cell_type)]['proportion']
                           for cell_type in df['Cell_type'].unique()
                           if (source1, cell_type) in params and params[(source1, cell_type)]['weights'].size > 0])

                weights2 = np.concatenate([params[(source2, cell_type)]['weights'] * params[(source2, cell_type)]['proportion']
                           for cell_type in df['Cell_type'].unique()
                           if (source2, cell_type) in params and params[(source2, cell_type)]['weights'].size > 0])


                weights1 = weights1 / np.sum(weights1)
                weights2 = weights2 / np.sum(weights2)



                if len(means1) > 0 and len(means2) > 0 and len(weights1) > 0 and len(weights2) > 0:
                    _, distGW2 = GW2(weights1, weights2, means1, means2, covs1, covs2,method ,normalized)
                    distance_matrix[i, j] = distGW2
                    distance_matrix[j, i] = distGW2
    adata.uns['QOT_Distance'] = distance_matrix
    return adata



def normalize_matrix(mat):
    min_val = np.min(mat)
    max_val = np.max(mat)
    range_val = max_val - min_val
    if range_val == 0:
        return mat - min_val
    else:
        return (mat - min_val) / range_val

def GW2(pi0,pi1,mu0,mu1,S0,S1,method,normalized):
    alpha = 1

    # return the GW2 discrete map and the GW2 distance between two GMM
    K0 = mu0.shape[0]
    K1 = mu1.shape[0]
    d  = mu0.shape[1]
    M  = np.zeros((K0,K1))
    trace_cov_sim_matrix= np.zeros((K0,K1))
    angle_sim_matrix= np.zeros((K0,K1))
    # First we compute the distance matrix between all Gaussians pairwise
    for k in range(K0):
        for l in range(K1):

            if method== "cosine":
                angle_similarity = cosine_similarity(mu0[k, :], mu1[l, :])
                M[k, l] = (1 - angle_similarity)

            elif method =="euclidean":


                M[k, l]  = euclidean_distance(mu0[k, :], mu1[l, :])

            elif method== "exact":

                M[k, l] = GaussianW2(mu0[k,:], mu1[l,:], S0[k,:,:], S1[l,:,:])

    if normalized == True:
      M = normalize_matrix(M)


    wstar     = ot.emd(pi0,pi1,M)         # discrete transport plan

    distGW2   = np.sum(wstar*M)
    return wstar,distGW2


import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score

def visualize_and_evaluate(adata, cmap='Blues_r', metric='precomputed'):
    """
    Normalizes the matrix, creates a heatmap, and calculates the silhouette score.

    Parameters:
    qot_labeled (array-like): The data to be normalized and visualized.
    labels (array-like): The labels for silhouette score calculation.
    cmap (str): The colormap for the heatmap.
    metric (str): The metric to use for silhouette score calculation.

    Returns:
    float: The silhouette score.
    """
    qot = adata.uns['QOT_Distance']
    df =  adata.uns['Datafame_for_use']
    # Normalize the matrix (assuming normalize_matrix is a defined function)
    normalized_matrix = normalize_matrix(qot)
    # Get unique sampleIDs in the order they appear from end to beginning
    unique_sampleIDs = df['sampleID'].unique()

# Initialize an empty list to hold the status for each unique sampleID
    true_labels = []

# Iterate over the unique sampleIDs
    for sampleID in unique_sampleIDs:
    # Find the status for the current sampleID
        status = df[df['sampleID'] == sampleID]['status'].iloc[0]
    # Append the found status to the true_labels list
        true_labels.append(status)

# Convert true_labels to a numpy array
    true_labels_array = np.array(true_labels)
    # Create heatmap
    sources = df['sampleID'].unique()
    qot_labeled_df = pd.DataFrame(normalized_matrix, index=sources, columns=sources)
    sns.clustermap(qot_labeled_df, cmap='Blues_r', figsize=(6, 6))
    plt.suptitle("QOT", y=1.05)
    plt.show()

    # Calculate and print the silhouette score
    score_1 = silhouette_score(normalized_matrix, true_labels_array, metric=metric)
    print("Silhouette Score:", score_1)

    score_2 = silhouette_score(normalized_matrix, true_labels_array, metric="cosine")
    print("Silhouette Score (PILOT Version):", score_2)
    return score_1,score_2
def custom_encode(labels,dataset_name):
    if dataset_name == "MYIO":
        mapping = {'IZ': 1, 'control': 0}
    elif dataset_name == "KID_T":
        mapping = {'<30': 2, '30-60': 1, '>60': 0}
    elif dataset_name == "KID_G":
        mapping = {'<30': 2, '30-60': 1, '>60': 0}
    elif dataset_name == "PDAC":
        mapping = {'T': 1, 'N': 0}
    elif dataset_name == "Kidney_AC_G":
        mapping = {'Case': 1, 'Normal': 0}
    elif dataset_name == "Kidney_AC_T":
        mapping = {'Case': 1, 'control': 0}
    elif dataset_name == "Covid":
        mapping = {'severe/critical': 2, 'mild/moderate': 1, 'control': 0}
    elif dataset_name == "Lupus":
        mapping = {'Case': 1, 'Healthy': 0}
    elif dataset_name == "FOLL":
        mapping = {'follicular lymphoma': 1, 'normal': 0}
    elif dataset_name == "Kidney_RNA":
        mapping ={'Normal': 0, 'chronic kidney disease': 1, 'acute kidney failure': 2}
    elif dataset_name == "diabete":
        mapping ={'normal': 0, 'type 2 diabetes mellitus': 2, 'endocrine pancreas disorder': 1, 'type 1 diabetes mellitus': 2}
    return np.array([mapping[label] for label in labels])

def trajectory_analysis(adata,knn_number=5,dataset_name = 'PDAC',flip_condition = False ):
    """
    Perform PHATE analysis and calculate AUCPR.

    Parameters:
    data (DataFrame or ndarray): The high-dimensional data.
    labels (Series or list): The labels for the data.
    custom_mapping (dict): Mapping for custom encoding of labels.
    true_label_array (list): Array of true labels for AUCPR calculation.

    Returns:
    float: The average AUCPR score.
    """
    qot = adata.uns['QOT_Distance']
    df =  adata.uns['Datafame_for_use']
    # Initialize PHATE
    true_labels = []
    unique_sampleIDs = df['sampleID'].unique()

# Initialize an empty list to hold the status for each unique sampleID
    true_labels = []
    # Iterate over the unique sampleIDs
    for sampleID in unique_sampleIDs:
    # Find the status for the current sampleID
        status = df[df['sampleID'] == sampleID]['status'].iloc[0]
    # Append the found status to the true_labels list
        true_labels.append(status)

    # Convert true_labels to a numpy array
    true_labels_array = np.array(true_labels)

    phate_op = ph.PHATE( n_components=2, knn_dist='precomputed_distance', knn = knn_number,random_state=42)
    # Normalize the matrix (assuming normalize_matrix is a defined function
    qot = normalize_matrix(qot)
    # Fit the model and obtain the two-dimensional embedding
    data_phate = phate_op.fit_transform(qot)

    # Create DataFrame for embedding
    embedding_df = pd.DataFrame({
        'Dim1': data_phate[:, 0],
        'Dim2': data_phate[:, 1],
        'Label': true_labels
    })


    # Sort embedding DataFrame by 'Dim1'
    sorted_embedding_df = embedding_df.sort_values(by='Dim1').reset_index(drop=True)
    sorted_labels = sorted_embedding_df['Label']



    if flip_condition:
      sorted_embedding_df = sorted_embedding_df.iloc[::-1].reset_index(drop=True)

# Add a rank column
    sorted_embedding_df['Rank'] = sorted_embedding_df.index + 1


    plot_embedding(sorted_embedding_df)
    # print(sorted_embedding_df)
    # Encode labels

    sorted_labels = sorted_embedding_df['Label']
    numeric_labels = custom_encode(sorted_labels,dataset_name)
    numeric_labels_sp = custom_encode(sorted_labels,dataset_name)

    # Calculate AUCPR
    if dataset_name == "MYIO":
        array =[0] * 13 +[1] * 7
    elif dataset_name == "KID_T":
        array =[0] * 400+ [1] * 177+ [2] * 57
        numeric_labels = np.array(numeric_labels)
        numeric_labels[numeric_labels > 0] = 1
    elif dataset_name == "KID_G":
        array =[0] * 400+ [1] * 177+ [2] * 57
        numeric_labels = np.array(numeric_labels)
        numeric_labels[numeric_labels > 0] = 1
    elif dataset_name == "PDAC":
        array =[0] * 11 +[1] * 24
    elif dataset_name == "Kidney_AC_G":
        array =[0] * 17 +[1] * 40
    elif dataset_name == "Kidney_AC_T":
        array =[0] * 17 +[1] * 40
    elif dataset_name == "Covid":
        array =[0] * 20 +[1] * 68+[2]*84
        numeric_labels = np.array(numeric_labels)
        numeric_labels[numeric_labels > 0] = 1
    elif dataset_name == "Lupus":
        array =[0] * 99 +[1] * 162
    elif dataset_name == "FOLL":
        array =[0] * 3 +[1] * 20
    elif dataset_name == "Kidney_RNA":
        array =[0] * 18 +[1] * 13 +[2] * 5
        numeric_labels = np.array(numeric_labels)
        numeric_labels[numeric_labels > 0] = 1
    elif dataset_name == "diabete":
        array =[0] * 26 +[1] * 6 +[2] * 24
        numeric_labels[numeric_labels > 0] = 1

    aucpr_scores = []
    precision, recall, _ = precision_recall_curve(numeric_labels, array)
    aucpr = auc(recall, precision)
    aucpr_scores.append(aucpr)


    average_aucpr = np.mean(aucpr_scores)
    print("AUCPR:", average_aucpr)

    # print(numeric_labels_sp)
    corr, p_value = spearmanr( sorted_embedding_df['Rank'] , numeric_labels_sp)

# Print the results
    print("Spearman's correlation coefficient:", corr)
    print("P-value:", p_value)

# Define a function to plot the embedding and save it as a PDF with smaller axis numerical values
def plot_embedding(embedding_df, figsize=(6, 6), font_size=12, legend_font_size=10, axis_num_size=10):
    plt.figure(figsize=figsize)
    for label, group in embedding_df.groupby('Label'):
        plt.scatter(group['Dim1'], group['Dim2'], label=label, alpha=0.7)
    plt.xlabel('Dimension 1', fontsize=font_size)
    plt.ylabel('Dimension 2', fontsize=font_size)
    plt.title('QOT', fontsize=font_size)
    plt.xticks(fontsize=axis_num_size)
    plt.yticks(fontsize=axis_num_size)
    plt.legend(fontsize=legend_font_size)
    plt.show()
def Run_QOT_shape(adata, gene_matrix='X_pca', type_cell='cell_types',
                id_col='sampleID', progession='status', dataset_type='rna',
                num_components_list=[1], random_state=2, min_samples_for_gmm=0,
                qot_method="cosine"):
    """
    Processes the PDAC dataset and evaluates Gaussian Mixture Models.

    Args:
        file_path (str): Path to the PDAC .h5ad file.
        gene_matrix (str): Key for the gene matrix in adata.
        type_cell (str): Key for cell types in adata.
        id_col (str): Key for sample IDs in adata.
        progession (str): Key for progression status in adata.
        dataset_type (str): Type of dataset (e.g., 'rna').
        num_components_list (list): List of num_components to iterate over.
        random_state (int): Random state for reproducibility.
        min_samples_for_gmm (int): Minimum samples required for GMM.
        qot_method (str): Method to calculate QOT (e.g., "cosine").

    Returns:
        dict: A dictionary containing timing and evaluation scores for each num_components.
    """
    # Read the dataset


    # Extract information from the dataset
    adata = Extract_Info(adata, gene_matrix=gene_matrix, type_cell=type_cell,
                id=id_col, progession=progession, dataset_type=dataset_type)
    # print("---Information extracted from the dataset---")
    # print()
    # Initialize results dictionary
    results = {}

    for num_components in num_components_list:
        # print(f"\nProcessing for num_components = {num_components}")
        
        total_start_time = time.time()  # Start time for the total process

        # Apply Gaussian Mixture Representation
        gmm_start_time = time.time()
        adata = Gaussian_Mixture_Representation(
            adata,
            num_components=num_components,
            random_state=random_state,
            min_samples_gmm=min_samples_for_gmm
        )
        gmm_end_time = time.time()
        gmm_time = gmm_end_time - gmm_start_time
        # print(f"---GMM completed in {gmm_time:.2f} seconds---")
        # print()
        # Calculate QOT (including WD)
        wd_start_time = time.time()
        adata = calculate_qot(adata, method=qot_method)
        wd_end_time = time.time()
        wd_time = wd_end_time - wd_start_time
        # print(f"---QOT calculation completed in {wd_time:.2f} seconds---")


        total_end_time = time.time()
        total_time = total_end_time - total_start_time
        # print(f"Total processing time: {total_time:.2f} seconds.")


    return adata
def Run_QOT(adata, gene_matrix='X_pca', type_cell='cell_types',
                id_col='sampleID', progession='status', dataset_type='rna',
                num_components_list=[1], random_state=2, min_samples_for_gmm=0,
                qot_method="cosine",normalized_set = False):
    """
    Processes the PDAC dataset and evaluates Gaussian Mixture Models.

    Args:
        file_path (str): Path to the PDAC .h5ad file.
        gene_matrix (str): Key for the gene matrix in adata.
        type_cell (str): Key for cell types in adata.
        id_col (str): Key for sample IDs in adata.
        progession (str): Key for progression status in adata.
        dataset_type (str): Type of dataset (e.g., 'rna').
        num_components_list (list): List of num_components to iterate over.
        random_state (int): Random state for reproducibility.
        min_samples_for_gmm (int): Minimum samples required for GMM.
        qot_method (str): Method to calculate QOT (e.g., "cosine").

    Returns:
        dict: A dictionary containing timing and evaluation scores for each num_components.
    """
    # Read the dataset


    # Extract information from the dataset
    adata = Extract_Info(adata, gene_matrix=gene_matrix, type_cell=type_cell,
                id=id_col, progession=progession, dataset_type=dataset_type)
    print("---Information extracted from the dataset---")
    print()
    # Initialize results dictionary
    results = {}

    for num_components in num_components_list:
        # print(f"\nProcessing for num_components = {num_components}")
        
        total_start_time = time.time()  # Start time for the total process

        # Apply Gaussian Mixture Representation
        gmm_start_time = time.time()
        adata = Gaussian_Mixture_Representation(
            adata,
            num_components=num_components,
            random_state=random_state,
            min_samples_gmm=min_samples_for_gmm
        )
        gmm_end_time = time.time()
        gmm_time = gmm_end_time - gmm_start_time
        print(f"---GMM completed in {gmm_time:.2f} seconds---")
        print()
        # Calculate QOT (including WD)
        wd_start_time = time.time()
        adata = calculate_qot(adata, method=qot_method,normalized=normalized_set)
        wd_end_time = time.time()
        wd_time = wd_end_time - wd_start_time
        print(f"---QOT calculation completed in {wd_time:.2f} seconds---")


        total_end_time = time.time()
        total_time = total_end_time - total_start_time
        # print(f"Total processing time: {total_time:.2f} seconds.")


    return adata



   
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text
from matplotlib.lines import Line2D
from gprofiler import GProfiler
import pandas as pd
pd.DataFrame.iteritems = pd.DataFrame.items

def extract_cells_from_gene_expression_for_clustering(adata,sample_col,col_cell,cell_list,path_results=None,normalization=True,n_top_genes=2000,highly_variable_genes_=False):


    """
    Extract and save gene expression data for specific cells for clustering analysis.

    Parameters:
        adata : AnnData object
            An Annotated Data (AnnData) object containing gene expression data.
        sample_col : str
            The column name in the adata.obs DataFrame containing sample IDs.
        col_cell : str
            The column name in the adata.obs DataFrame containing cell type labels.
        cell_list : list
            A list of cell types for which gene expression data will be extracted and saved.
        path_results : str or None, optional (default=None)
            The path to the directory where the extracted data will be saved as CSV files.
            If None, the default path 'Results_PILOT/cells/' will be used.
        normalization : bool, optional (default=True)
            Whether to normalize the gene expression data by total count and apply log1p transformation.
        n_top_genes : int, optional (default=2000)
            The number of top highly variable genes to select. Only applicable if highly_variable_genes_ is True.
        highly_variable_genes_ : bool, optional (default=False)
            Whether to select highly variable genes for analysis.

    Returns:
        Gene exp. dataframe
    """
    for cell in cell_list:
        adata_new = adata[adata.obs[col_cell].isin([cell]),:]
        if normalization:
            sc.pp.normalize_total(adata_new, target_sum=1e4)
            sc.pp.log1p(adata_new)

        if highly_variable_genes_:

            sc.pp.highly_variable_genes(adata_new, n_top_genes=n_top_genes)
                # Access the list of highly variable genes
            highly_variable_genes = adata_new.var['highly_variable']
            df=adata_new[:,highly_variable_genes].X
            df=pd.DataFrame(df.toarray())
            highly_variable_gene_names = adata_new.var_names[np.array(adata_new.var['highly_variable'])]
            df.columns=list(highly_variable_gene_names)
        else:

            df=adata_new[:,adata_new.var_names].X
            df=pd.DataFrame(df.toarray())
            df.columns=adata_new.var_names


        df['sampleID']=list(adata_new.obs[sample_col])

 



        return df
def volcano_plot(scores, foldchanges, p_values, cell_type, feature1, feature2, fc_thr = 1, pv_thr = 1,
                 figsize = (20,20), output_path = None,n_p=5,n_n=5,font_size=18, marker='o',
                             color='w',
                             markersize=8,
                             font_weight_legend='normal',
                             size_legend=12,dpi=100


                             ):

    """
    Generate a volcano plot to visualize gene expression significance.

    Parameters:
        scores : pandas Series
            A pandas Series containing the expression scores for genes.
        foldchanges : array-like
            An array-like containing the fold changes for genes.
        p_values : pandas Series
            A pandas Series containing the p-values for genes.
        cell_type : str
            The name of the cell type being analyzed.
        feature1 : str
            The name of the first feature being compared.
        feature2 : str
            The name of the second feature being compared.
        fc_thr : float, optional (default=1)
            The threshold for log2FoldChange to determine significance.
        pv_thr : float, optional (default=1)
            The threshold for negative log10 of p-value to determine significance.
        figsize : tuple, optional (default=(15, 15))
            The size of the plot figure.
        output_path : str, optional (default=None)
            The path to save the output plot. If None, the plot will be displayed.
        n_p : int, optional (default=5)
            The number of labels that the user wants to show over the plot for positive threshold.
        n_n : int, optional (default=5)
            The number of labels that the user wants to show over the plot for negative threshold.
        font_size : int, optional (default=18)
            Font size for the plot.
        marker : str, optional (default='o')
            Marker style for data points in the plot.
        color : str, optional (default='w')
            Color for data points in the plot.
        markersize : int, optional (default=8)
            Marker size for data points in the plot.
        font_weight_legend : str, optional (default='normal')
            Font weight for legend text.
        size_legend : int, optional (default=12)
            Font size for legend text.
        dpi : int, optional
            Dots per inch for the saved plot image. Default is 100.

    Returns:
        None
    """



    df = pd.DataFrame(columns=['log2FoldChange', 'nlog10', 'symbol'])
    df['log2FoldChange'] = foldchanges
    df['nlog10'] = -np.log10(p_values.values)
    df['symbol'] = scores.index.values

    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.dropna(subset=["nlog10"], how="all", inplace=True)


    selected_labels = df.loc[ (np.abs(df.log2FoldChange) >= fc_thr) & (df['nlog10'] >= pv_thr)]['symbol'].values
    group1_selected_labels = df.loc[ (df.log2FoldChange <= -fc_thr) & (df['nlog10'] >= pv_thr)]['symbol'].values
    # pd.DataFrame(group1_selected_labels).to_csv(output_path + "/significant_genes_" + str(cell_type) + "_" + str(feature1) + ".csv")

    group2_selected_labels = df.loc[ (df.log2FoldChange >= fc_thr) & (df['nlog10'] >= pv_thr)]['symbol'].values
    # pd.DataFrame(group2_selected_labels).to_csv(output_path + "/significant_genes_" + str(cell_type) + "_" + str(feature2) + ".csv")

    def map_shape(symbol):
        if symbol in selected_labels:
            return 'important'
        return 'not'

    df['color'] = df[['log2FoldChange', 'symbol', 'nlog10']].apply(map_color, fc_thrr = fc_thr, pv_thrr = pv_thr, axis = 1)
    df['shape'] = df.symbol.map(map_shape)
    df['baseMean'] = df.nlog10*10


    plt.figure(figsize = figsize, frameon=False, dpi=100)
    plt.style.use('default')


    #plt.xlim(-xlim, xlim)
    ax = sns.scatterplot(data = df, x = 'log2FoldChange', y = 'nlog10',
                         hue = 'color', hue_order = ['no', 'very higher','higher', 'mix', 'very lower', 'lower'],
                         palette = ['lightgrey', '#d62a2b', '#D62A2B7A',
                                    'lightgrey', '#1f77b4', '#1F77B47D'],
                         style = 'shape', style_order = ['not', 'important'],
                         markers = ['o', 'o'],
                         size = 'baseMean', sizes = (40, 800)
                        )

    ax.axhline(pv_thr, zorder = 0, c = 'k', lw = 2, ls = '--')
    ax.axvline(fc_thr, zorder = 0, c = 'k', lw = 2, ls = '--')
    ax.axvline(-fc_thr, zorder = 0, c = 'k', lw = 2, ls = '--')

    texts = []
    filtered_df = df.loc[df['nlog10'] >= pv_thr]
    subset_labels_fold_change_pos = filtered_df.loc[filtered_df['log2FoldChange'] >= fc_thr]
    subset_labels_fold_change_pos = subset_labels_fold_change_pos.sort_values(by='nlog10', ascending=False)
    subset_labels_fold_change_pos = subset_labels_fold_change_pos.head(n_p)['symbol'].values

    subset_labels_fold_change_neg = filtered_df.loc[filtered_df['log2FoldChange'] <= -fc_thr]
    subset_labels_fold_change_neg = subset_labels_fold_change_neg.sort_values(by='nlog10', ascending=False)
    subset_labels_fold_change_neg = subset_labels_fold_change_neg.head(n_n)['symbol'].values
    # Combine the subsets of genes
    subset_labels = np.concatenate([subset_labels_fold_change_pos, subset_labels_fold_change_neg])
    for i in range(len(df)):
        if df.iloc[i].symbol in subset_labels:
            if df.iloc[i].nlog10 >= pv_thr and (df.iloc[i].log2FoldChange >= fc_thr):
                texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].symbol,
                                     fontsize = font_size, weight = 'bold', family = 'sans-serif'))
            if df.iloc[i].nlog10 >= pv_thr and ( df.iloc[i].log2FoldChange <= -fc_thr):
                texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].symbol,
                                     fontsize = font_size, weight = 'bold', family = 'sans-serif'))
    adjust_text(texts)
   # for i in range(len(df)):
    #    if df.iloc[i].symbol in subset_labels:
     #       if df.iloc[i].nlog10 >= pv_thr and (df.iloc[i].log2FoldChange >= fc_thr):
      #          texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].symbol,
          #                           fontsize = 16, weight = 'bold', family = 'sans-serif'))
       #     if df.iloc[i].nlog10 >= pv_thr and ( df.iloc[i].log2FoldChange <= -fc_thr):
        #        texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].symbol,
         #                            fontsize = 16, weight = 'bold', family = 'sans-serif'))
    #adjust_text(texts)

    custom_lines = [Line2D([0], [0], marker=marker, color=color, markerfacecolor='#d62a2b', markersize=markersize),
                   Line2D([0], [0], marker=marker, color=color, markerfacecolor='#1f77b4', markersize=markersize)]

    plt.legend(custom_lines, ['Higher expressions in ' + feature2, 'Higher expressions in ' + feature1],loc = 1,
               bbox_to_anchor = (1,1.1), frameon = False, prop = {'weight': font_weight_legend, 'size': size_legend})

    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(2)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(width = 2)

    plt.title("Expression Score \n "+feature1+" - "+feature2, fontsize = font_size)
    plt.xticks(size = font_size, weight = 'bold')
    plt.yticks(size = font_size, weight = 'bold')
    plt.xlabel("$log_{2}$ (Fold Change)", size = font_size)
    plt.ylabel("-$log_{10}$ (P-value)", size = font_size)

#     plt.savefig(filename, dpi = 100, bbox_inches = 'tight', facecolor = 'white')
    
    plt.show()
def compute_diff_expressions(adata,cell_type: str = None,
                             proportions: pd.DataFrame = None,
                             selected_genes: list = None,
                             font_size:int=18,
                             group1: str = 'Tumor 1',
                             group2: str = 'Tumor 2',
                             label_name: str = 'Predicted_Labels',
                             fc_thr: float = 0.5,
                             pval_thr: float = 0.01,
                             sample_col:str='sampleID',
                             col_cell:str ='cell_types',
                             path=None,
                             normalization=False,
                             n_top_genes=2000,
                             highly_variable_genes_=True,
                             number_n=5,
                             number_p=5,
                             marker='o',
                             color='w',
                             markersize=8,
                             font_weight_legend='normal',
                             size_legend=12,
                             figsize=(6,6),dpi=100
                             ):

    """
    Using limma R package, lmFit fits a linear model using weighted least squares for each gene.
    Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models.
    Empirical Bayes smoothing of standard errors (shrinks standard errors
    that are much larger or smaller than those from other genes towards the average standard error).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    cell_type : str, optional
        Specify cell type name to check its differential expression genes. The default is None.
    proportions : pd.DataFrame, optional
        Cell types proportions in each sample. The default is None.
    selected_genes : list, optional
        Specify gene names to be considered for checking their differentiation.
    font_size : int, optional
        Font size for plot labels and legends. The default is 18.
    group1 : str, optional
        Name of the first patient sub-group for comparison. The default is 'Tumor 1'.
    group2 : str, optional
        Name of the second patient sub-group for comparison. The default is 'Tumor 2'.
    label_name : str, optional
        Name of the column containing the labels of patient sub-groups. The default is 'Predicted_Labels'.
    fc_thr : float, optional
        Specify the fold change threshold. The default is 0.5.
    pval_thr : float, optional
        Specify the adjusted p-value threshold. The default is 0.01.
    sample_col : str, optional
        Name of the column containing sample IDs. The default is 'sampleID'.
    col_cell : str, optional
        Name of the column containing cell type annotations. The default is 'cell_types'.
    path : str, optional
        Path to save the results. The default is None.
    normalization : bool, optional
        Perform gene expression normalization. The default is False.
    n_top_genes : int, optional
        Number of top variable genes to consider. The default is 2000.
    highly_variable_genes_ : bool, optional
        Determine highly variable genes. The default is True.
    number_n : int, optional
        The number of labels that the user wants to show over the plot for negative thresholds. The default is 5.
    number_p : int, optional
        The number of labels that the user wants to show over the plot for positive thresholds. The default is 5.
    marker : str, optional
        Marker style for the labels in the volcano plot. The default is 'o'.
    color : str, optional
        Marker color for the labels in the volcano plot. The default is 'w'.
    markersize : int, optional
        Marker size for the labels in the volcano plot. The default is 8.
    font_weight_legend : str, optional
        Font weight for legend labels. The default is 'normal'.
    size_legend : int, optional
        Font size for legend labels. The default is 12.
    figsize: tuple, optional
        Figure size. The default is (15,15).
    dpi : int, optional
        Dots per inch for the saved plot image. Default is 100.

    Returns
    -------
    None

    Generates and displays a volcano plot of fold changes between two interested patient sub-groups.
    Saves a statistical table of each gene.
    Saves significantly differentiated genes in each group.
    """

    
    if cell_type not in adata.uns:
        cells=extract_cells_from_gene_expression_for_clustering(adata,sample_col=sample_col,col_cell=col_cell,cell_list=[cell_type],path_results=path,normalization=normalization,n_top_genes=n_top_genes,highly_variable_genes_=highly_variable_genes_)
        #cells = pd.read_csv(path_to_result + "/cells/" + cell_type + ".csv", index_col = 0)

    elif cell_type in adata.uns :
         cells=adata.uns[cell_type]


    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    pandas2ri.activate()

    # prepare data for R
    proportions.index = proportions['sampleID']

    if selected_genes is None:
        selected_genes = cells.iloc[:,1:-1].columns
    data = cells[selected_genes]
    pred_labels = pd.DataFrame()
    pls = proportions.loc[cells['sampleID']]
    pred_labels['Predicted_Labels'] = pls[label_name]
    pred_labels['sampleID'] = pls['sampleID']

    # load R packages and data
    R=robjects.r
    R('library(limma)')
    R.assign('data',data)
    R.assign('pred_labels', pred_labels)
    R.assign('selected_groups', [group1, group2])
    R('selected_pred_labels <- pred_labels[which(pred_labels$Predicted_Labels %in% selected_groups),]')
    R('subresult <- data[row.names(selected_pred_labels),]')

    # delete for memory
    del data
    del pred_labels

    # run limma
    print('run limma lmFit')
    R('fit <- limma::lmFit(t(subresult), design = unclass(as.factor(selected_pred_labels$Predicted_Labels)))')
    print('run limma eBayes')
    R('fit <-  limma::eBayes(fit)')
    R('res <- limma::topTable(fit, n = 2000)')
    R('res <- res[colnames(data), ]')

    # get results
    res = R('''res''')


    pv_thr = -np.log10(pval_thr)
    volcano_plot(cells[selected_genes].transpose(), res['logFC'], res['adj.P.Val'],
                 cell_type, group1, group2, fc_thr, pv_thr,
                 figsize = figsize, output_path = path_to_result,n_p=number_p,n_n=number_n,font_size=font_size, marker=marker,
                             color=color,
                             markersize=markersize,
                             font_weight_legend=font_weight_legend,
                             size_legend=size_legend,dpi=dpi)

def compute_shapley_values(
    adata,
    n_components=50,
    num_components_list=[1],
    random_state=2,
    min_samples_for_gmm=16,
    qot_method="cosine",
    gene_matrix='X_pca',
    type_cell='cell_types',
    id_col='sampleID',
    progression='status',
    dataset_type='rna'
):
    import numpy as np
    import pandas as pd
    import time
    from sklearn.decomposition import PCA
    from sklearn.metrics import silhouette_score

    # Define helper functions
    def normalize_matrix(mat):
        min_val = np.min(mat)
        max_val = np.max(mat)
        range_val = max_val - min_val
        if range_val == 0:
            return mat - min_val
        else:
            return (mat - min_val) / range_val

    def evaluate(adata, cmap='Blues_r', metric='precomputed'):
        """
        Normalizes the matrix, creates a heatmap, and calculates the silhouette score.

        Parameters:
        adata: AnnData object

        Returns:
        float: The silhouette score.
        """
        qot = adata.uns['QOT_Distance']
        df = adata.uns['Datafame_for_use']
        # Normalize the matrix
        normalized_matrix = normalize_matrix(qot)
        # Get unique sampleIDs in the order they appear from end to beginning
        unique_sampleIDs = df['sampleID'].unique()

        # Initialize an empty list to hold the status for each unique sampleID
        true_labels = []

        # Iterate over the unique sampleIDs
        for sampleID in unique_sampleIDs:
            # Find the status for the current sampleID
            status = df[df['sampleID'] == sampleID]['status'].iloc[0]
            # Append the found status to the true_labels list
            true_labels.append(status)

        # Convert true_labels to a numpy array
        true_labels_array = np.array(true_labels)
        # Calculate and print the silhouette score
        score_1 = silhouette_score(normalized_matrix, true_labels_array, metric=metric)
        # print("Silhouette Score:", score_1)

        score_2 = silhouette_score(normalized_matrix, true_labels_array, metric="cosine")
        # print("Silhouette Score (PILOT Version):", score_2)
        return score_1, score_2

    print("---PCA---")
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(adata.X)
    adata.obsm['X_pca'] = principal_components
    weights = pca.components_
    

    # Call the processing function (ensure Run_QOT is defined)
    adata = Run_QOT_shape(
        adata,
        gene_matrix=gene_matrix,
        type_cell=type_cell,
        id_col=id_col,
        progession=progression,
        dataset_type=dataset_type,
        num_components_list=num_components_list,
        random_state=random_state,
        min_samples_for_gmm=min_samples_for_gmm,
        qot_method=qot_method
    )
    print("---QOT Complete for Base Model---")
    # Evaluate the base model
    score_1_base, score_2_base = evaluate(adata)

    results = {}
    total_start_time = time.time()
    n_components_actual = adata.obsm['X_pca'].shape[1]

    # Iterate over each PCA component
    for i in range(n_components_actual):
        print(f"---Processing leave-out PCA component: {i}---")
        adata_temp = adata.copy()

        # Remove the i-th column from the PCA matrix
        PCA_temp = np.delete(adata_temp.obsm['X_pca'], i, axis=1)
        adata_temp.obsm['X_pca'] = PCA_temp

        # Run QOT again
        adata_temp = Run_QOT_shape(
            adata_temp,
            gene_matrix=gene_matrix,
            type_cell=type_cell,
            id_col=id_col,
            progession=progression,
            dataset_type=dataset_type,
            num_components_list=num_components_list,
            random_state=random_state,
            min_samples_for_gmm=min_samples_for_gmm,
            qot_method=qot_method
        )

        # Evaluate the model without the i-th PCA component
        score_1, score_2 = evaluate(adata_temp)
        results[i] = {'score_1': score_1, 'score_2': score_2}

    # Compute adjusted Shapley values
    
    v_N = score_1_base
    n = len(results)
    shapley_values = {}

    # print("Adjusted Shapley Values for each PCA component:")
    for col_index, data in results.items():
        v_N_minus_i = data['score_1']
        phi_i = (1 / n) * (v_N - v_N_minus_i)
        shapley_values[col_index] = phi_i
        print(f"PCA Component {col_index}: Shapley Value = {phi_i:.6f}")

    # Compute gene contributions
    phi_components = np.array([shapley_values[i] for i in range(n_components_actual)])
    weights_T = weights.T  # Shape: (n_genes, n_components)
    phi_genes = np.dot(weights_T, phi_components)  # Shape: (n_genes,)

    # Normalize gene contributions
    total_phi_components = np.sum(phi_components)
    total_phi_genes = np.sum(phi_genes)

    if total_phi_genes != 0:
        phi_genes_normalized = phi_genes * (total_phi_components / total_phi_genes)
    else:
        phi_genes_normalized = phi_genes

    # Create a DataFrame for gene contributions
    gene_names = adata.var_names  # Shape: (n_genes,)
    df_gene_contributions = pd.DataFrame({
        'Gene': gene_names,
        'Contribution': phi_genes_normalized
    })

    # Sort genes by their contributions
    df_gene_contributions_sorted = df_gene_contributions.sort_values(by='Contribution', ascending=False)

    # Display the top contributing genes
    # print("Top 20 Genes Contributing to Model Performance:")
    # print(df_gene_contributions_sorted.head(20))

    # Optionally, save the results to a CSV file

    return df_gene_contributions_sorted




def plot_top_gene_contributions(df_gene_contributions_sorted, top_n=20, output_file='rank_gene.pdf'):
    """
    Plots the top N genes contributing to model performance.

    Parameters:
    - df_gene_contributions_sorted (DataFrame): A DataFrame sorted by 'Contribution' containing gene names and their contributions.
    - top_n (int): The number of top genes to display in the plot.
    - output_file (str): The filename for saving the plot.

    Returns:
    - None: The function displays and saves the plot.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np

    # Select the top_n genes
    df = df_gene_contributions_sorted.head(top_n)

    # Set the Seaborn style without gridlines
    sns.set(style="white", context="talk")

    # Number of colors needed
    n_colors = len(df)

    # Define the start and stop range for the colormap (values between 0 and 1)
    start = 0.3  # Start at 30% into the colormap
    stop = 0.7   # Stop at 70% into the colormap

    # Get the colormap
    cmap = plt.get_cmap('Blues')

    # Generate the colors within the specified range
    colors = cmap(np.linspace(start, stop, n_colors))

    # Reverse the colors array to reverse the color mapping
    colors = colors[::-1]

    # Initialize the figure
    plt.figure(figsize=(6, 6))

    # Create the barplot with reversed colors
    ax = sns.barplot(
        x='Contribution',
        y='Gene',
        data=df,
        palette=colors
    )

    # Remove gridlines
    ax.grid(False)

    # Remove top and right spines for a cleaner look
    sns.despine(ax=ax)

    # Set background color to white (optional)
    ax.set_facecolor('white')

    # Add title and labels
    plt.title(f'Top {top_n} Genes Contributing to Model Performance', fontsize=20, pad=20)
    plt.xlabel('Contribution', fontsize=16)
    plt.ylabel('Gene', fontsize=16)

    # Adjust layout
    plt.tight_layout()

    # Save the plot

    # Display the plot
    plt.show()

def map_color(a, fc_thrr, pv_thrr):
    """
    Map colors based on specified thresholds for Fold Change and p-value.

    Parameters:
        a : tuple
            A tuple containing log2FoldChange, symbol, and negative log10 of p-value.
        fc_thrr : float
            The threshold for log2FoldChange to determine different color mappings.
        pv_thrr : float
            The threshold for negative log10 of p-value to determine different color mappings.

    Returns:
        str
            A string indicating the color mapping based on the provided thresholds and input values.
    """
    log2FoldChange, symbol, nlog10 = a
    if log2FoldChange >= fc_thrr and nlog10 >= pv_thrr:
        return 'very higher'
    elif log2FoldChange <= -fc_thrr and nlog10 >= pv_thrr:
        return 'very lower'
    elif log2FoldChange >= fc_thrr and nlog10 < pv_thrr:
        return 'higher'
    elif log2FoldChange <= -fc_thrr and nlog10 < pv_thrr:
        return 'lower'
    elif abs(log2FoldChange) < fc_thrr and nlog10 >= pv_thrr:
        return 'mix'
    else:
        return 'no'
    
def compute_diff_expressions(adata,cell_type: str = None,
                             proportions: pd.DataFrame = None,
                             selected_genes: list = None,
                             font_size:int=18,
                             group1: str = 'Tumor 1',
                             group2: str = 'Tumor 2',
                             label_name: str = 'Predicted_Labels',
                             fc_thr: float = 0.5,
                             pval_thr: float = 0.01,
                             sample_col:str='sampleID',
                             col_cell:str ='cell_types',
                             path=None,
                             normalization=False,
                             n_top_genes=2000,
                             highly_variable_genes_=True,
                             number_n=5,
                             number_p=5,
                             marker='o',
                             color='w',
                             markersize=8,
                             font_weight_legend='normal',
                             size_legend=12,
                             figsize=(15,15),dpi=100
                             ):

    """
    Using limma R package, lmFit fits a linear model using weighted least squares for each gene.
    Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models.
    Empirical Bayes smoothing of standard errors (shrinks standard errors
    that are much larger or smaller than those from other genes towards the average standard error).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    cell_type : str, optional
        Specify cell type name to check its differential expression genes. The default is None.
    proportions : pd.DataFrame, optional
        Cell types proportions in each sample. The default is None.
    selected_genes : list, optional
        Specify gene names to be considered for checking their differentiation.
    font_size : int, optional
        Font size for plot labels and legends. The default is 18.
    group1 : str, optional
        Name of the first patient sub-group for comparison. The default is 'Tumor 1'.
    group2 : str, optional
        Name of the second patient sub-group for comparison. The default is 'Tumor 2'.
    label_name : str, optional
        Name of the column containing the labels of patient sub-groups. The default is 'Predicted_Labels'.
    fc_thr : float, optional
        Specify the fold change threshold. The default is 0.5.
    pval_thr : float, optional
        Specify the adjusted p-value threshold. The default is 0.01.
    sample_col : str, optional
        Name of the column containing sample IDs. The default is 'sampleID'.
    col_cell : str, optional
        Name of the column containing cell type annotations. The default is 'cell_types'.
    path : str, optional
        Path to save the results. The default is None.
    normalization : bool, optional
        Perform gene expression normalization. The default is False.
    n_top_genes : int, optional
        Number of top variable genes to consider. The default is 2000.
    highly_variable_genes_ : bool, optional
        Determine highly variable genes. The default is True.
    number_n : int, optional
        The number of labels that the user wants to show over the plot for negative thresholds. The default is 5.
    number_p : int, optional
        The number of labels that the user wants to show over the plot for positive thresholds. The default is 5.
    marker : str, optional
        Marker style for the labels in the volcano plot. The default is 'o'.
    color : str, optional
        Marker color for the labels in the volcano plot. The default is 'w'.
    markersize : int, optional
        Marker size for the labels in the volcano plot. The default is 8.
    font_weight_legend : str, optional
        Font weight for legend labels. The default is 'normal'.
    size_legend : int, optional
        Font size for legend labels. The default is 12.
    figsize: tuple, optional
        Figure size. The default is (15,15).
    dpi : int, optional
        Dots per inch for the saved plot image. Default is 100.

    Returns
    -------
    None

    Generates and displays a volcano plot of fold changes between two interested patient sub-groups.
    Saves a statistical table of each gene.
    Saves significantly differentiated genes in each group.
    """

    path_to_result='Results_PILOT'

    if os.path.exists(path_to_result + "/cells/" + cell_type + ".csv"):

        cells = pd.read_csv(path_to_result + "/cells/" + cell_type + ".csv", index_col = 0)

    elif cell_type not in adata.uns:
        cells=extract_cells_from_gene_expression_for_clustering(adata,sample_col=sample_col,col_cell=col_cell,cell_list=[cell_type],path_results=path,normalization=normalization,n_top_genes=n_top_genes,highly_variable_genes_=highly_variable_genes_)
        #cells = pd.read_csv(path_to_result + "/cells/" + cell_type + ".csv", index_col = 0)

    elif cell_type in adata.uns :
         cells=adata.uns[cell_type]


    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    pandas2ri.activate()

    # prepare data for R
    proportions.index = proportions['sampleID']

    if selected_genes is None:
        selected_genes = cells.iloc[:,1:-1].columns
    data = cells[selected_genes]
    pred_labels = pd.DataFrame()
    # print(proportions)
    pls = proportions.loc[cells['sampleID']]
    pred_labels['Predicted_Labels'] = pls[label_name]
    pred_labels['sampleID'] = pls['sampleID']

    # load R packages and data
    R=robjects.r
    R('library(limma)')
    R.assign('data',data)
    R.assign('pred_labels', pred_labels)
    R.assign('selected_groups', [group1, group2])
    R('selected_pred_labels <- pred_labels[which(pred_labels$Predicted_Labels %in% selected_groups),]')
    R('subresult <- data[row.names(selected_pred_labels),]')

    # delete for memory
    del data
    del pred_labels

    # run limma
    print('run limma lmFit')
    R('fit <- limma::lmFit(t(subresult), design = unclass(as.factor(selected_pred_labels$Predicted_Labels)))')
    print('run limma eBayes')
    R('fit <-  limma::eBayes(fit)')
    R('res <- limma::topTable(fit, n = 2000)')
    R('res <- res[colnames(data), ]')

    # get results
    res = R('''res''')

    if not os.path.exists(path_to_result+'/Diff_Expressions_Results/'+cell_type):
            os.makedirs(path_to_result+'/Diff_Expressions_Results/'+cell_type)
    path_to_result=path_to_result+'/Diff_Expressions_Results/'+cell_type+'/'
    res.to_csv(path_to_result + "/diff_expressions_stats_" + cell_type + ".csv")

    pv_thr = -np.log10(pval_thr)
    volcano_plot(cells[selected_genes].transpose(), res['logFC'], res['adj.P.Val'],
                 cell_type, group1, group2, fc_thr, pv_thr,
                 figsize = figsize, output_path = path_to_result,n_p=number_p,n_n=number_n,font_size=font_size, marker=marker,
                             color=color,
                             markersize=markersize,
                             font_weight_legend=font_weight_legend,
                             size_legend=size_legend,dpi=dpi)