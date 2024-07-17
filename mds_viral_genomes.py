#%%
#import libraries
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from sklearn.manifold import MDS
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from mpl_toolkits.mplot3d import Axes3D
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

#%%
# upload data
data = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/genomes_to_pcs.csv').dropna()
lifestyle_labels = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/network_labels.csv')

#%%

pcs_presence_absence_df = pd.DataFrame(index = data['contig_id'].unique(), columns = data['pc_id'].unique())
#creates dataframe where the rows are the contig (virus) names and the columns are the protein cluster ids

#%%
#only contains ind, non, and pre samples
pcs_presence_absence_df = pcs_presence_absence_df[pcs_presence_absence_df.index.str.contains('pal-|ind-|pre-')]

#sort in order for easy finding
pcs_presence_absence_cols = pcs_presence_absence_df.sort_index(axis=1)
pcs_presence_absence_df = pcs_presence_absence_cols.sort_index()

#%%
#adds 0 where there is no protein cluster for that genome and adds 1 where there is that protein cluster present
pcs_presence_absence_df = pcs_presence_absence_df.fillna(0)
for _, row in data.iterrows():
    genome = row['contig_id']
    protein_cluster = row['pc_id']
    if genome in pcs_presence_absence_df.index and protein_cluster in pcs_presence_absence_df.columns:
        pcs_presence_absence_df.at[genome, protein_cluster] = 1

#%%
#only with contig_id/genome and lifestyle
metadata = lifestyle_labels.iloc[:,[0,-1]]
metadata['Time'] = 'Industrial'
metadata.loc[lifestyle_labels['Genome'].str.contains('pal-', na=False), 'Time'] = 'Pre-modern'
metadata.loc[lifestyle_labels['Genome'].str.contains('pre-', na=False), 'Time'] = 'Non-Industrial'

#%%
#put genomes in order
metadata = metadata[metadata['Genome'].str.contains('pal-|ind-|pre-')]
metadata = metadata.sort_values(by='Genome')

#%%
# Apply MDS
mds = MDS(n_components=2, metric = True, random_state=42)
embedding = mds.fit_transform(pcs_presence_absence_df)

#%%
embedding_df_lifestyle = pd.DataFrame(embedding, columns=['MDS1', 'MDS2'])
embedding_df_lifestyle['metadata'] = metadata['Lifestyle']

#%%
embedding_df_time = pd.DataFrame(embedding, columns=['MDS1', 'MDS2'])
embedding_df_time['metadata'] = metadata['Time']


#%%
#plot MDS projection colored by lifestyle
plt.figure(figsize=(15,15))
sns.scatterplot(data=embedding_df_lifestyle, x='MDS1', y='MDS2', hue='metadata', palette='icefire', alpha = 0.7)
plt.title('MDS Projection of Viral Genomes')
plt.xlabel('MDS1')
plt.ylabel('MDS2')
plt.legend(title='Lifestyle')
plot_filename = 'Viral_Genomics/outputs/mds_viral_genomes_by_lifestyle.png'
plt.savefig(plot_filename)

#%%
#plot MDS projection colored by Time period
embedding_df_time['metadata'] = metadata['Time']
plt.figure(figsize=(15,15))
sns.scatterplot(data=embedding_df_time, x='MDS1', y='MDS2', hue='metadata', palette='icefire', alpha = 0.7)
plt.title('MDS Projection of Viral Genomes')
plt.xlabel('MDS1')
plt.ylabel('MDS2')
plt.legend(title='Time Period')
plot_filename = 'Viral_Genomics/outputs/mds_viral_genomes_by_time.png'
plt.savefig(plot_filename)

#%%

# Extract the MDS coordinates
mds_coords = embedding_df_time[['MDS1', 'MDS2']]

# Determine the optimal number of clusters using the elbow method and silhouette score
def determine_optimal_clusters(data, max_clusters=10):
    inertia = []
    silhouette_scores = []
    cluster_range = range(2, max_clusters + 1)

    for n_clusters in cluster_range:
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        cluster_labels = kmeans.fit_predict(data)

        inertia.append(kmeans.inertia_)
        silhouette_scores.append(silhouette_score(data, cluster_labels))

    return cluster_range, inertia, silhouette_scores


# Plot the elbow method and silhouette scores
def plot_metrics(cluster_range, inertia, silhouette_scores):
    fig, ax1 = plt.subplots()

    ax1.set_xlabel('Number of clusters')
    ax1.set_ylabel('Inertia', color='tab:blue')
    ax1.plot(cluster_range, inertia, 'o-', color='tab:blue', label='Inertia')
    ax1.tick_params(axis='y', labelcolor='tab:blue')

    ax2 = ax1.twinx()
    ax2.set_ylabel('Silhouette Score', color='tab:red')
    ax2.plot(cluster_range, silhouette_scores, 's-', color='tab:red', label='Silhouette Score')
    ax2.tick_params(axis='y', labelcolor='tab:red')

    fig.tight_layout()
    plt.title('Elbow Method and Silhouette Scores')
    plt.savefig('Viral_Genomics/outputs/mds_viral_genomes_metrics.png')


# Determine optimal clusters and plot metrics
cluster_range, inertia, silhouette_scores = determine_optimal_clusters(mds_coords)
plot_metrics(cluster_range, inertia, silhouette_scores)

#%%
# Optimal number of clusters from metrics plot
optimal_clusters = 4 #

# Perform k-means clustering with the optimal number of clusters
kmeans = KMeans(n_clusters=optimal_clusters, random_state=42)
mds_coords['Cluster'] = kmeans.fit_predict(mds_coords)

# Visualize the clusters in the MDS plot
plt.figure(figsize=(10, 8))
plt.scatter(mds_coords['MDS1'], mds_coords['MDS2'], c=mds_coords['Cluster'], cmap='viridis', s=50)
plt.xlabel('MDS1')
plt.ylabel('MDS2')
plt.title('K-means Clustering of Viral Genomes')
plt.colorbar(label='Cluster')
plt.savefig('Viral_Genomics/outputs/mds_viral_genomes_clusters.png')

#%%
embedding_df_time['Cluster'] = mds_coords['Cluster']

#%%
cluster_dist = embedding_df_time.groupby(['Cluster', 'metadata']).size().unstack(fill_value = 0)
colors = {'Industrial' : 'blue', 'Non-Industrial':'red', 'Pre-modern' : 'green'}

cluster_dist.plot(kind='bar', stacked=True, figsize=(10,7),color=[colors[col] for col in cluster_dist.columns])
plt.xlabel('Cluster')
plt.ylabel('Count')
plt.title('Distribution of Time Periods within Clusters')
plt.legend(title='Time Period')
plt.savefig('Viral_Genomics/outputs/mds_viral_genomes_distribution_in_clusters.png')