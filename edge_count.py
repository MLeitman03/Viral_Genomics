#%%
# Import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#%%
# Upload data (vcontact2 output)
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
df['Source Time Label'] = df['Source']
df['Target Time Label'] = df['Target']

# Update 'Source Time Label' based on conditions
df.loc[df['Source'].str.contains('ind-', na=False), 'Source Time Label'] = 'Industrial'
df.loc[df['Source'].str.contains('pal-', na=False), 'Source Time Label'] = 'Pre-modern'
df.loc[df['Source'].str.contains('pre-', na=False), 'Source Time Label'] = 'Non-Industrial'

# If you also need to update 'Target Time Label', you can do the same
df.loc[df['Target'].str.contains('ind-', na=False), 'Target Time Label'] = 'Industrial'
df.loc[df['Target'].str.contains('pal-', na=False), 'Target Time Label'] = 'Pre-modern'
df.loc[df['Target'].str.contains('pre-', na=False), 'Target Time Label'] = 'Non-Industrial'



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

target_counts_df = pd.DataFrame(index = ['Source- Industrial','Source- Non-Industrial', 'Source- Pre-modern'], columns = ['Industrial','Non-Industrial', 'Pre-modern'])
target_counts_df.loc['Source- Industrial'] = ind_target_counts_normalized
target_counts_df.loc['Source- Non-Industrial'] = non_target_counts_normalized
target_counts_df.loc['Source- Pre-modern'] = pre_target_counts_normalized

#%%
#Make heatmap
target_counts_df = target_counts_df.astype(float)
plt.figure()
sns.set(style="whitegrid")
sns.heatmap(target_counts_df, cmap='Blues')
plt.title('Normalized Target Time Label Counts by Source Time Label')
plt.ylabel('Source Time Label')
plt.xlabel('Target Time Label')
plot_filename = 'Viral_Genomics/outputs/heatmap_edge_counts_by_time.png'
plt.savefig(plot_filename)

#%%
def create_dict(df):
    '''
    Create a dictionary where each key is a unique genome and the values are an array of target lifestyles associated with that genome.
    Parameters: df: the input dataframe (ie. ind_source_df)
    '''

    sorted_df = df.sort_values(by='Source', ascending=True, inplace=False)
    target_count_dict = {}

    for index, row in sorted_df.iterrows():
        key = row['Source']
        value = row['Target Time Label']
        if key in target_count_dict:
            target_count_dict[key].append(value)
        else:
            target_count_dict[key] = [value]

    return target_count_dict

#%%
#create dictionaries for each time period source genome
ind_dict = create_dict(ind_source_df)
non_dict = create_dict(non_source_df)
pre_dict = create_dict(pre_source_df)

#%%
unique_ind_genome_count = len(ind_dict.keys())
unique_non_genome_count = len(non_dict.keys())
unique_pre_genome_count = len(pre_dict.keys())

#%%

#initialize lists for dataframe columns
ind_genomes = []
ind_ind_counts = []
ind_non_ind_counts = []
ind_pre_counts = []
ind_database_counts = []
ind_lengths = []

for key, values in ind_dict.items():
    ind_genomes.append(key)
    ind_lengths.append(len(values))
    ind_ind_counts.append(values.count('Industrial'))
    ind_non_ind_counts.append(values.count('Non-Industrial'))
    ind_pre_counts.append(values.count('Pre-modern'))
    ind_database_counts.append(len(values) - (values.count('Industrial') + values.count('Non-Industrial') + values.count('Pre-modern')))

#dataframe where each genome has a count of industrial targets, non-industrial, and pre-modern
ind_target_per_genome = pd.DataFrame({'Genome': ind_genomes,
                                      'Industrial Targets': ind_ind_counts,
                                      'Non-Industrial Targets': ind_non_ind_counts,
                                      'Pre-modern Targets': ind_pre_counts,
                                      'Database Targets': ind_database_counts,
                                      'Number of connections': ind_lengths})

#%%

#normalized by total number of edges/connections
normalized_ind_target_per_genome = ind_target_per_genome
for column in ['Industrial Targets', 'Non-Industrial Targets', 'Pre-modern Targets', 'Database Targets']:
    normalized_ind_target_per_genome[column] = ind_target_per_genome[column] / ind_target_per_genome['Number of connections']

#%%

#initialize lists for dataframe columns
non_genomes = []
non_ind_counts = []
non_non_ind_counts = []
non_pre_counts = []
non_database_counts = []
non_lengths = []

for key, values in non_dict.items():
    non_genomes.append(key)
    non_lengths.append(len(values))
    non_ind_counts.append(values.count('Industrial'))
    non_non_ind_counts.append(values.count('Non-Industrial'))
    non_pre_counts.append(values.count('Pre-modern'))
    non_database_counts.append(len(values) - (values.count('Industrial') + values.count('Non-Industrial') + values.count('Pre-modern')))

#dataframe where each genome has a count of industrial targets, non-industrial, and pre-modern
non_target_per_genome = pd.DataFrame({'Genome': non_genomes,
                                      'Industrial Targets': non_ind_counts,
                                      'Non-Industrial Targets': non_non_ind_counts,
                                      'Pre-modern Targets': non_pre_counts,
                                      'Database Targets': non_database_counts,
                                      'Number of connections': non_lengths})

#%%

#normalized by total number of edges/connections
normalized_non_target_per_genome = non_target_per_genome
for column in ['Industrial Targets', 'Non-Industrial Targets', 'Pre-modern Targets', 'Database Targets']:
    normalized_non_target_per_genome[column] = non_target_per_genome[column] / non_target_per_genome['Number of connections']

#%%
#initialize lists for dataframe columns
pre_genomes = []
pre_ind_counts = []
pre_non_ind_counts = []
pre_pre_counts = []
pre_database_counts = []
pre_lengths = []

for key, values in pre_dict.items():
    pre_genomes.append(key)
    pre_lengths.append(len(values))
    pre_ind_counts.append(values.count('Industrial'))
    pre_non_ind_counts.append(values.count('Non-Industrial'))
    pre_pre_counts.append(values.count('Pre-modern'))
    pre_database_counts.append(len(values) - (values.count('Industrial') + values.count('Non-Industrial') + values.count('Pre-modern')))

#dataframe where each genome has a count of industrial targets, non-industrial, and pre-modern
pre_target_per_genome = pd.DataFrame({'Genome': pre_genomes,
                                      'Industrial Targets': pre_ind_counts,
                                      'Non-Industrial Targets': pre_non_ind_counts,
                                      'Pre-modern Targets': pre_pre_counts,
                                      'Database Targets': pre_database_counts,
                                      'Number of connections': pre_lengths})

#%%

#normalized by total number of edges/connections
normalized_pre_target_per_genome = pre_target_per_genome
for column in ['Industrial Targets', 'Non-Industrial Targets', 'Pre-modern Targets', 'Database Targets']:
    normalized_pre_target_per_genome[column] = pre_target_per_genome[column] / pre_target_per_genome['Number of connections']


#%%
normalized_combined = pd.concat([
    normalized_ind_target_per_genome.assign(Source='Industrial'),
    normalized_non_target_per_genome.assign(Source='Non-Industrial'),
    normalized_pre_target_per_genome.assign(Source='Pre-modern')
])

# Melt the combined DataFrame to long format
df_melted_combined = pd.melt(
    normalized_combined,
    id_vars=['Genome', 'Source'],
    value_vars=['Industrial Targets', 'Non-Industrial Targets', 'Pre-modern Targets'],
    var_name='Target Type',
    value_name='Proportion'
)

# Create the violin plot
plt.figure(figsize=(14, 10))
sns.violinplot(data=df_melted_combined, x='Source', y='Proportion', hue='Target Type', palette='rocket')
plt.title('Proportion of Targets by Genome for Each Source Type (Edge Count by Source)')
plt.ylabel('Proportion')
plt.xlabel('Source Type')
plt.legend(title='Target Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plot_filename = 'Viral_Genomics/outputs/violinplot_edge_counts_by_time.png'
plt.savefig(plot_filename)

#Add confidence intervals, look for statistical significance

#%%
#export source files
ind_source_df.to_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/industrial_source_genomes.csv', index=False)
non_source_df.to_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/nonindustrial_source_genomes.csv', index=False)
pre_source_df.to_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/pre_modern_source_genomes.csv', index=False)