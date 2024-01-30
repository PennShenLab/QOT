import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture
import ot 
import scipy.stats as sps
import scipy.linalg as spl
import seaborn as sns
import pandas as pd
from sklearn.metrics import precision_recall_curve, auc
import numpy as np
import phate

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
    


def Gaussian_Mixture_Representation(adata, num_components=5, random_state=2, min_samples_for_gmm=0):
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

        if num_samples >= num_components + min_samples_for_gmm:
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
        elif num_samples == 2:  # Exactly 1 samples, use itself
            mean_sample = data.mean(axis=0).values.reshape(1, -1)
            key = (source, label)
            params[key] = {
            'means': mean_sample,
            'covariances': np.zeros((1, data.shape[1], data.shape[1])),
            'weights': np.array([1]),
            'proportion': len(group) / len(df[df['sampleID'] == source])
            }
        # elif num_samples == 3:  # Exactly 1 samples, use itself
        #     mean_sample = data.mean(axis=0).values.reshape(1, -1)
        #     key = (source, label)
        #     params[key] = {
        #     'means': mean_sample,
        #     'covariances': np.zeros((1, data.shape[1], data.shape[1])),
        #     'weights': np.array([1]),
        #     'proportion': len(group) / len(df[df['sampleID'] == source])
        #     }
        # else:
        #     print(f"Skipped source {source}, label {label}, number of samples {num_samples} - not enough data for {num_components} components")

    
    adata.uns['GMM_Representation'] = params


def euclidean_distance(vec1, vec2):
    """Calculate the Euclidean distance between two vectors."""
    return np.linalg.norm(vec1 - vec2)

def cosine_similarity(vec1, vec2):
    """Calculate the cosine similarity between two vectors."""
    return np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))

def trace_covariance_similarity(cov1, cov2):
    """Calculate a similarity measure based on the trace of the product of covariance matrices."""
    prod = np.dot(cov1, cov2)
    return np.trace(prod)



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
                print(f"Calculating GW2 distance {calculation_count}/{total_calculations}: {source1} to {source2}")

    
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
            
            if method== "alternative":
            
            
            
                angle_similarity = cosine_similarity(mu0[k, :], mu1[l, :])
                #trace_cov_sim_matrix[k, l] = trace_covariance_similarity(S0[k, :, :], S1[l, :, :])
                # M[k, l] = euclidean_distance(mu0[k, :], mu1[l, :])


                M[k, l] = (1 - angle_similarity) 
             
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
    score = silhouette_score(normalized_matrix, true_labels_array, metric=metric)
    print("Silhouette Score:", score)

    score = silhouette_score(normalized_matrix, true_labels_array, metric="cosine")
    print("Silhouette Score (PILOT Version):", score)

def custom_encode(labels,dataset_name):
    if dataset_name == "MYIO":
        mapping = {'IZ': 1, 'control': 0}
    elif dataset_name == "KID_T":
        mapping = {'<30': 2, '30-60': 1, '>60': 0}
    elif dataset_name == "KID_G":
        mapping = {'<30': 2, '30-60': 1, '>60': 0}
    elif dataset_name == "PDAC":
        mapping = {'T': 0, 'N': 1}
    # Mapping: '<30' -> 2, '30-60' -> 1, '>60' -> 0
    
    return np.array([mapping[label] for label in labels])

def perform_analysis_phate_and_calculate_aucpr(adata,knn_number=5,dataset_name = 'PDAC'):
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

    phate_op = phate.PHATE(n_components=2, knn_dist='precomputed_distance', knn = knn_number,random_state=2)
    qot = normalize_matrix(qot)
    # Fit the model and obtain the two-dimensional embedding
    data_phate = phate_op.fit_transform(qot)

    # Create DataFrame for embedding
    embedding_df = pd.DataFrame({
        'Dim1': data_phate[:, 0],
        'Dim2': data_phate[:, 1],
        'Label': true_labels
    })
    plot_embedding(embedding_df)
    # Sort embedding DataFrame by 'Dim1'
    sorted_embedding_df = embedding_df.sort_values(by='Dim1')

    sorted_labels = sorted_embedding_df['Label']
    sorted_labels
    # Encode labels

    numeric_labels = custom_encode(sorted_labels,dataset_name)

    # Calculate AUCPR
    if dataset_name == "MYIO":
        array =[0] * 13 +[1] * 7
    elif dataset_name == "KID_T":
        array =[0] * 400+ [1] * 177+ [2] * 57
    elif dataset_name == "KID_G":
        array =[0] * 400+ [1] * 177+ [2] * 57
    elif dataset_name == "PDAC":
        array =[0] * 24 +[1] * 11

    aucpr_scores = []

    for class_label in np.unique(array):
    
        true_binary = (array == class_label).astype(int)
        
        pred_binary = (numeric_labels == class_label).astype(int)
        
        precision, recall, _ = precision_recall_curve(true_binary, pred_binary)
        aucpr = auc(recall, precision)
        aucpr_scores.append(aucpr)

    # Calculate the average AUCPR
    average_aucpr = np.mean(aucpr_scores)
    print("AUCPR:", average_aucpr)

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

