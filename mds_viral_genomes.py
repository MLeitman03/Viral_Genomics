#%%
#import libraries
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from sklearn.manifold import MDS
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
mds = MDS(n_components=2, random_state=42)
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
