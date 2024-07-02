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
# Plot edge weights

plt.figure(figsize=(12, 6))

# Create a violin plot
sns.violinplot(x='Source Time Label', y='edge_weight', hue='Target Time Label', data=df, split=True, inner="quart", palette="muted")

# Add titles and labels
plt.title('Violin Plot of Source Genomes by Target Genomes')
plt.xlabel('Source Time Label')
plt.ylabel('Edge Weight')

# Show plot
plt.show()