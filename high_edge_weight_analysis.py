#%%
# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#%%
# Import data
ind = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/industrial_source_genomes.csv')
non = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/nonindustrial_source_genomes.csv')
pre = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/pre_modern_source_genomes.csv')
network_labels_df = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/network_labels.csv')


#%%
ind_high_edge_weights = ind.loc[ind['edge_weight'] == 300]

non_high_edge_weights = non.loc[non['edge_weight'] == 300]

pre_high_edge_weights = pre.loc[pre['edge_weight'] == 300]

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
ind_dict = create_dict(ind_high_edge_weights)
non_dict = create_dict(non_high_edge_weights)
pre_dict = create_dict(pre_high_edge_weights)

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
normalized_ind_target_per_genome = ind_target_per_genome.copy()
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
normalized_non_target_per_genome = non_target_per_genome.copy()
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

plot_filename = 'Viral_Genomics/outputs/violinplot_edge_counts_by_time_of_high_edge_weight_300.png'
plt.savefig(plot_filename)

#%%
network_labels_df = network_labels_df.iloc[:,[0,5,6,7,8,9,10,11,12]]

#%%
ind_high_edge_weights_with_cluster_label = pd.merge(ind_high_edge_weights, network_labels_df, left_on='Source', right_on='Genome', how='left')
non_high_edge_weights_with_cluster_label = pd.merge(non_high_edge_weights, network_labels_df, left_on='Source', right_on='Genome', how='left')
pre_high_edge_weights_with_cluster_label = pd.merge(pre_high_edge_weights, network_labels_df, left_on='Source', right_on='Genome', how='left')

#%%
ind_target_counts = ind_high_edge_weights_with_cluster_label['Target Time Label'].value_counts()
#normalized
ind_target_counts_normalized = ind_target_counts / ind_target_counts.sum()

#target time label counts for non-industrial source
non_target_counts = non_high_edge_weights_with_cluster_label['Target Time Label'].value_counts()
#normalized
non_target_counts_normalized = non_target_counts / non_target_counts.sum()


#%%
#Combine normalized target counts

target_counts_df = pd.DataFrame(index = ['Source- Industrial','Source- Non-Industrial'], columns = ['Industrial','Non-Industrial'])
target_counts_df.loc['Source- Industrial'] = ind_target_counts_normalized
target_counts_df.loc['Source- Non-Industrial'] = non_target_counts_normalized


#%%
#Make heatmap
targets_counts_df = target_counts_df.T
target_counts_df = target_counts_df.astype(float)
plt.figure()
sns.set(style="whitegrid")
sns.heatmap(target_counts_df, cmap='Blues')
plt.title('Normalized Target Time Label Counts by Source Time Label')
plt.ylabel('Source Time Label')
plt.xlabel('Target Time Label')
plot_filename = 'Viral_Genomics/outputs/heatmap_edge_counts_by_time_for_edge_weight_300.png'
plt.savefig(plot_filename)

#%%

genomes_with_cluster_label = pd.concat([ind_high_edge_weights_with_cluster_label, non_high_edge_weights_with_cluster_label, pre_high_edge_weights_with_cluster_label])