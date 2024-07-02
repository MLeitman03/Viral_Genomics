#%%
# Import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#%%
# Upload data
edge_df = pd.read_table(filepath_or_buffer='/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/edge_weights.ntw', sep=' ', header=None)
network_labels_df = pd.read_csv(filepath_or_buffer='/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/network_labels.csv')


#%%
# Merge dataframes on Genome label
edge_df_genome_col = edge_df[[0]]
network_labels_df_genome_col = network_labels_df[['Genome']]
df = pd.merge(edge_df, network_labels_df, left_on=0, right_on='Genome', how='left')
df.rename(columns={0:'Source', 1:'Target', 2:'edge_weight', 'Lifestyle':'Source Lifestyle'}, inplace=True)

#%%
#assign lifestyle to target genomes
target_genome = df[['Target']]
genome_with_lifestyle = network_labels_df.loc[:, ['Genome', 'Lifestyle']]

# Merge dataframes
target_genome_lifestyle = pd.merge(target_genome, genome_with_lifestyle, left_on='Target', right_on='Genome', how='left')

# Drop duplicates and create a new DataFrame
target_genome_lifestyle = target_genome_lifestyle.drop_duplicates()

# Rename the column
target_genome_lifestyle.rename(columns={'Lifestyle': 'Target Lifestyle'}, inplace=True)
target_genome_lifestyle = target_genome_lifestyle.drop(columns=['Genome'])

df = pd.merge(df, target_genome_lifestyle, on = 'Target', how='left')

#%%
# Add time period column for source and target genomes

# Initialize the new columns with default values
df['Source Time Label'] = 'Database Virus'
df['Target Time Label'] = 'Database Virus'

# Update 'Source Time Label' based on conditions
df.loc[df['Source'].str.contains('ind', na=False), 'Source Time Label'] = 'Industrial'
df.loc[df['Source'].str.contains('pal', na=False), 'Source Time Label'] = 'Pre-modern'
df.loc[df['Source'].str.contains('pre', na=False), 'Source Time Label'] = 'Non-Industrial'

# If you also need to update 'Target Time Label', you can do the same
df.loc[df['Target'].str.contains('ind', na=False), 'Target Time Label'] = 'Industrial'
df.loc[df['Target'].str.contains('pal', na=False), 'Target Time Label'] = 'Pre-modern'
df.loc[df['Target'].str.contains('pre', na=False), 'Target Time Label'] = 'Non-Industrial'



#%%
# New dataframe with count
ind_source_df = df.loc[df['Source Time Label'] == 'Industrial'].loc[:, ['Source', 'Target', 'Source Time Label', 'Target Time Label', 'edge_weight', 'Source Lifestyle', 'Target Lifestyle']]
non_source_df = df.loc[df['Source Time Label'] == 'Non-Industrial'].loc[:, ['Source', 'Target', 'Source Time Label', 'Target Time Label', 'edge_weight', 'Source Lifestyle','Target Lifestyle']]
pre_source_df = df.loc[df['Source Time Label'] == 'Pre-modern'].loc[:, ['Source', 'Target', 'Source Time Label', 'Target Time Label', 'edge_weight', 'Source Lifestyle','Target Lifestyle']]

#does ind-DNK_MH0192_k99_103327 source go to pal-BMS_UT30.3_k141_1024356 target?

testing = ind_source_df.loc[(ind_source_df['Source'] == 'ind-DNK_MH0192_k99_103327') & (ind_source_df['Target'].str.contains('pal')), 'Target']

#yes, meaning each source is a target and vice versa

#%%
#target time label counts for industiral source
ind_target_counts = ind_source_df['Target Time Label'].value_counts()
#normalized
ind_target_counts_normalized = ind_target_counts / ind_target_counts.sum()

#target time label counts for non-industrial source
non_target_counts = non_source_df['Target Time Label'].value_counts()
#normalized
non_target_counts_normalized = non_target_counts / non_target_counts.sum()

#target time label counts for pre-modern source
pre_target_counts = pre_source_df['Target Time Label'].value_counts()
#normalized
pre_target_counts_normalized = pre_target_counts / pre_target_counts.sum()

#maybe see what type of database viruses they are most similar to?

#%%
#Combine normalized target counts

target_counts_df = pd.DataFrame(index = ['Industrial','Non-Industrial', 'Pre-modern'], columns = ['Industrial','Non-Industrial', 'Pre-modern', 'Database Virus'])
target_counts_df.loc['Industrial'] = ind_target_counts_normalized
target_counts_df.loc['Non-Industrial'] = non_target_counts_normalized
target_counts_df.loc['Pre-modern'] = pre_target_counts_normalized

#%%
#Make heatmap
target_counts_df = target_counts_df.astype(float)
plt.figure()
sns.heatmap(target_counts_df)
plt.title('Normalized Target Time Label Counts by Source Time Label')
plt.ylabel('Source Time Label')
plt.xlabel('Target Time Label')
plt.show()

import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg', 'Agg', 'TkAgg', etc.

