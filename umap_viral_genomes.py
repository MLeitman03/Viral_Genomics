#%%
#import libraries
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
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
# Apply UMAP
reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
embedding = reducer.fit_transform(pcs_presence_absence_df)

#%%
embedding_df_lifestyle = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'])
embedding_df_lifestyle['metadata'] = metadata['Lifestyle']

#%%
embedding_df_time = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'])
embedding_df_time['metadata'] = metadata['Time']


#%%
#plot UMAP projection colored by lifestyle
plt.figure(figsize=(15,15))
sns.scatterplot(data=embedding_df_lifestyle, x='UMAP1', y='UMAP2', hue='metadata', palette='viridis', alpha = 0.7)
plt.title('UMAP Projection of Viral Genomes')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.legend(title='Lifestyle')
plot_filename = 'Viral_Genomics/outputs/umap_viral_genomes_by_lifestyle.png'
plt.savefig(plot_filename)

#%%
#plot UMAP projection colored by Time period
embedding_df_time['metadata'] = metadata['Time']
plt.figure(figsize=(15,15))
sns.scatterplot(data=embedding_df_time, x='UMAP1', y='UMAP2', hue='metadata', palette='viridis', alpha = 0.7)
plt.title('UMAP Projection of Viral Genomes')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.legend(title='Time Period')
plot_filename = 'Viral_Genomics/outputs/umap_viral_genomes_by_time.png'
plt.savefig(plot_filename)

#%%
#Separate metadatas
metadata_lytic = metadata.loc[metadata['Lifestyle']=='Lytic',]
metadata_temp = metadata.loc[metadata['Lifestyle']=='Temperate',]
metadata_ind = metadata.loc[metadata['Time']=='Industrial',]
metadata_non = metadata.loc[metadata['Time']=='Non-Industrial',]
metadata_pre = metadata.loc[metadata['Time']=='Pre-modern',]


#%%
df_lytic = pcs_presence_absence_df[pcs_presence_absence_df.index.isin(metadata_lytic['Genome'])]
df_temp = pcs_presence_absence_df[pcs_presence_absence_df.index.isin(metadata_temp['Genome'])]
df_ind = pcs_presence_absence_df[pcs_presence_absence_df.index.isin(metadata_ind['Genome'])]
df_non = pcs_presence_absence_df[pcs_presence_absence_df.index.isin(metadata_non['Genome'])]
df_pre = pcs_presence_absence_df[pcs_presence_absence_df.index.isin(metadata_pre['Genome'])]

#%%
# Apply UMAP for lytic viruses separated by time period
reducer_lytic = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
embedding_lytic = reducer_lytic.fit_transform(df_lytic)

embedding_df_lytic = pd.DataFrame(embedding_lytic, columns=['UMAP1', 'UMAP2'])
embedding_df_lytic['metadata'] = metadata_lytic['Time'].tolist()

#plot UMAP projection colored by time
plt.figure(figsize=(15,15))
sns.scatterplot(data=embedding_df_lytic, x='UMAP1', y='UMAP2', hue='metadata', palette='viridis', alpha = 0.7)
plt.title('UMAP Projection of Lytic Viral Genomes')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.legend(title='Time')
plot_filename = 'Viral_Genomics/outputs/umap_viral_genomes_by_time_for_lytic_viruses.png'
plt.savefig(plot_filename)

#%%
# Apply UMAP for temperate viruses separated by time period
reducer_temp = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
embedding_temp = reducer_temp.fit_transform(df_temp)

embedding_df_temp = pd.DataFrame(embedding_temp, columns=['UMAP1', 'UMAP2'])
embedding_df_temp['metadata'] = metadata_temp['Time'].tolist()

#plot UMAP projection colored by time
plt.figure(figsize=(15,15))
sns.scatterplot(data=embedding_df_temp, x='UMAP1', y='UMAP2', hue='metadata', palette='viridis', alpha = 0.7)
plt.title('UMAP Projection of Temperate Viral Genomes')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.legend(title='Time')
plot_filename = 'Viral_Genomics/outputs/umap_viral_genomes_by_time_for_temp_viruses.png'
plt.savefig(plot_filename)

#%%
# Apply UMAP for industrial viruses separated by lifestyle
reducer_ind = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
embedding_ind = reducer_ind.fit_transform(df_ind)

embedding_df_ind = pd.DataFrame(embedding_ind, columns=['UMAP1', 'UMAP2'])
embedding_df_ind['metadata'] = metadata_ind['Lifestyle'].tolist()

#plot UMAP projection colored by lifestyle
plt.figure(figsize=(15,15))
sns.scatterplot(data=embedding_df_ind, x='UMAP1', y='UMAP2', hue='metadata', palette='viridis', alpha = 0.7)
plt.title('UMAP Projection of Industrial Viral Genomes')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.legend(title='Lifestyle')
plot_filename = 'Viral_Genomics/outputs/umap_viral_genomes_by_time_for_ind_viruses.png'
plt.savefig(plot_filename)

#%%
# Apply UMAP for non-industrial viruses separated by lifestyle
reducer_non = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
embedding_non = reducer_non.fit_transform(df_non)

embedding_df_non = pd.DataFrame(embedding_non, columns=['UMAP1', 'UMAP2'])
embedding_df_non['metadata'] = metadata_non['Lifestyle'].tolist()

#plot UMAP projection colored by lifestyle
plt.figure(figsize=(15,15))
sns.scatterplot(data=embedding_df_non, x='UMAP1', y='UMAP2', hue='metadata', palette='viridis', alpha = 0.7)
plt.title('UMAP Projection of Non-Industrial Viral Genomes')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.legend(title='Lifestyle')
plot_filename = 'Viral_Genomics/outputs/umap_viral_genomes_by_time_for_non_viruses.png'
plt.savefig(plot_filename)

#%%
# Apply UMAP for pre-modern viruses separated by lifestyle
reducer_pre = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
embedding_pre = reducer_pre.fit_transform(df_pre)

embedding_df_pre = pd.DataFrame(embedding_pre, columns=['UMAP1', 'UMAP2'])
embedding_df_pre['metadata'] = metadata_pre['Lifestyle'].tolist()

#plot UMAP projection colored by lifestyle
plt.figure(figsize=(15,15))
sns.scatterplot(data=embedding_df_pre, x='UMAP1', y='UMAP2', hue='metadata', palette='viridis', alpha = 0.7)
plt.title('UMAP Projection of Pre-modern Viral Genomes')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.legend(title='Lifestyle')
plot_filename = 'Viral_Genomics/outputs/umap_viral_genomes_by_time_for_pre_viruses.png'
plt.savefig(plot_filename)



#%%
#make dendogram and hiearchical cluster
def dendrogram_HC(embedding_df, output_dir, file_prefix):
    # dendrogram
    Z = linkage(embedding_df[['UMAP1', 'UMAP2']], method='ward')
    #ward minimizes the within-cluster variance

    # plot deondrogram
    plt.figure(figsize=(40, 15))
    dendrogram(Z, labels=embedding_df['Genome'].to_list(), leaf_rotation=90, leaf_font_size=8)
    plt.title('Hieararchical Clustering Dendrogram of Viral Genomes')
    plt.xlabel('Viral Genomes')
    plt.ylabel('Distance')
    plt.savefig(f"{output_dir}{file_prefix}_dendrogram.png")


    clusters = fcluster(Z, 7.0, criterion='distance')
    embedding_df['Hieararchical_Cluster'] = clusters

    #plot HC
    plt.figure(figsize=(15, 15))
    sns.scatterplot(data=embedding_df, x='UMAP1', y='UMAP2', hue='Hieararchical_Cluster', palette='Spectral', alpha = 0.7)
    plt.title('Hieararchical Clustering of Viral Genomes')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.legend(title='Hieararchical Cluster')
    plt.savefig(f"{output_dir}{file_prefix}_hiearchical.png")

#%%
#save to 'Viral_Genomics/outputs/'

embedding_df_lifestyle['Genome'] = ' '
embedding_df_lifestyle['Genome'] = metadata['Genome'].tolist()

embedding_df_lytic['Genome'] = ' '
embedding_df_lytic['Genome'] = metadata_lytic['Genome'].tolist()

embedding_df_temp['Genome'] = ' '
embedding_df_temp['Genome'] = metadata_temp['Genome'].tolist()

embedding_df_ind['Genome'] = ' '
embedding_df_ind['Genome'] = metadata_ind['Genome'].tolist()

embedding_df_non['Genome'] = ' '
embedding_df_non['Genome'] = metadata_non['Genome'].tolist()

embedding_df_pre['Genome'] = ' '
embedding_df_pre['Genome'] = metadata_pre['Genome'].tolist()

embeddings = {'total' :embedding_df_lifestyle,
              'industrial' :embedding_df_ind,
              'non-industrial' : embedding_df_non,
              'pre-modern':embedding_df_pre,
              'lytic':embedding_df_lytic,
              'temperate': embedding_df_temp}

for file_prefix, embedding in embeddings.items() :
    dendrogram_HC(embedding, 'Viral_Genomics/outputs/', file_prefix)

#%%
embedding_df_lifestyle.to_csv('Viral_Genomics/embedding_df_lifestyle.csv')
embedding_df_ind.to_csv('Viral_Genomics/embedding_df_ind.csv')
embedding_df_non.to_csv('Viral_Genomics/embedding_df_non.csv')
embedding_df_pre.to_csv('Viral_Genomics/embedding_df_pre.csv')
embedding_df_lytic.to_csv('Viral_Genomics/embedding_df_lytic.csv')
embedding_df_temp.to_csv('Viral_Genomics/embedding_df_temp.csv')
