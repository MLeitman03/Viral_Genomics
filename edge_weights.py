#%%
# Import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#%%
# Import data
ind = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/industrial_source_genomes.csv')
non = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/nonindustrial_source_genomes.csv')
pre = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/pre_modern_source_genomes.csv')

#%%
combined = pd.concat([ind, non, pre])

#%%

desired_labels = ['Industrial', 'Non-Industrial', 'Pre-modern']
filtered_combined = combined[combined['Target Time Label'].isin(desired_labels)]

# Melt the combined DataFrame to long format
df_melted_combined = pd.melt(
    filtered_combined,
    id_vars=['Source', 'Source Time Label', 'Target Time Label'],
    value_vars=['edge_weight'],
    var_name='Metric',
    value_name='Value'
)

# Create the violin plot
plt.figure(figsize=(14, 10))
sns.violinplot(data=df_melted_combined, x='Source Time Label', y='Value', hue='Target Time Label', order = ['Industrial', 'Non-Industrial', 'Pre-modern'], hue_order = ['Industrial', 'Non-Industrial', 'Pre-modern'], palette='rocket')
plt.title('Edge Weights by Source and Target Time Labels')
plt.ylabel('Edge Weight')
plt.xlabel('Source Time Label')
plt.legend(title='Target Time Label', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plot_filename = 'Viral_Genomics/outputs/violinplot_edge_weights_by_time.png'
plt.savefig(plot_filename)

#%%

desired_labels_lifestyle = ['Lytic', 'Temperate']
filtered_combined_lifestyle = combined[combined['Target Lifestyle'].isin(desired_labels_lifestyle)]

# Melt the combined DataFrame to long format
df_melted_combined_lifestyle = pd.melt(
    filtered_combined_lifestyle,
    id_vars=['Source', 'Source Lifestyle', 'Target Lifestyle'],
    value_vars=['edge_weight'],
    var_name='Metric',
    value_name='Value'
)

# Create the violin plot
plt.figure(figsize=(14, 10))
sns.violinplot(data=df_melted_combined_lifestyle, x='Source Lifestyle', y='Value', hue='Target Lifestyle', palette='rocket')
plt.title('Edge Weights by Source and Target Lifestyle')
plt.ylabel('Edge Weight')
plt.xlabel('Source Lifestyle')
plt.legend(title='Target Lifestyle', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plot_filename = 'Viral_Genomics/outputs/violinplot_edge_weights_by_lifestyle.png'
plt.savefig(plot_filename)

#%%
# distribution of edge weights
plt.figure(figsize=(10, 6))
plt.hist(combined['edge_weight'], bins=30, edgecolor='k', color='blue')
plt.title('Distribution of Edge Weights')
plt.xlabel('Edge Weight')
plt.ylabel('Frequency')
plt.savefig('Viral_Genomics/outputs/histogram_edge_weights.png')

#%%
plt.figure(figsize=(10, 6))
plt.hist(combined['edge_weight'], bins=30, edgecolor='k', color='blue')
plt.yscale('log')
plt.title('Distribution of Edge Weights')
plt.xlabel('Edge Weight')
plt.ylabel('Frequency (log-scaled)')
plt.savefig('Viral_Genomics/outputs/histogram_edge_weights_log_scale.png')

#%%
mean = np.mean(combined['edge_weight'])
std = np.std(combined['edge_weight'])
plt.figure(figsize=(10, 6))
plt.hist(combined['edge_weight'], bins=30, edgecolor='k', color='lightsteelblue')
plt.yscale('log')

plt.axvline(mean, color='green', linestyle='dashed', linewidth=2, label='Mean')
plt.axvline(mean + std, color='orange', linestyle='dashed', linewidth=2, label='Mean + 1 Std')
plt.axvline(mean + 2 * std, color='yellow', linestyle='dashed', linewidth=2, label='Mean + 2 Std')
plt.axvline(mean + 3 * std, color='purple', linestyle='dashed', linewidth=2, label='Mean + 3 Std')

plt.legend()

plt.title('Distribution of Pre-modern Edge Weights')
plt.xlabel('Edge Weight')
plt.ylabel('Frequency (log-scaled)')
plt.savefig('Viral_Genomics/outputs/histogram_edge_weights_all_log_scale_with_cutoffs.png')

#%%
pre_only_combined = combined.loc[combined['Source Time Label'] == 'Pre-modern']

#%%
# distribution of pre-modern source edge weights

plt.figure(figsize=(10, 6))
plt.hist(pre_only_combined['edge_weight'], bins=30, edgecolor='k', color='blue')
plt.title('Distribution of Pre-modern Edge Weights')
plt.xlabel('Edge Weight')
plt.ylabel('Frequency')
plt.savefig('Viral_Genomics/outputs/histogram_edge_weights_pre_only.png')

#%%
# distribution of pre-modern source edge weights now with log-scaled y-axis

plt.figure(figsize=(10, 6))
plt.hist(pre_only_combined['edge_weight'], bins=30, edgecolor='k', color='blue')
plt.yscale('log')
plt.title('Distribution of Pre-modern Edge Weights')
plt.xlabel('Edge Weight')
plt.ylabel('Frequency (log-scaled)')
plt.savefig('Viral_Genomics/outputs/histogram_edge_weights_pre_only_log_scale.png')

#%%

excluded_labels = ['Industrial', 'Non-Industrial', 'Pre-modern']
only_database_target_combined = combined[~combined['Target Time Label'].isin(excluded_labels)]

#%%
# Melt the combined DataFrame to long format
df_melted_combined_database = pd.melt(
    only_database_target_combined,
    id_vars=['Source', 'Source Time Label'],
    value_vars=['edge_weight'],
    var_name='Metric',
    value_name='Value'
)

# Create the violin plot
plt.figure(figsize=(14, 10))
sns.violinplot(data=df_melted_combined_database, x='Source Time Label', y='Value',order = ['Industrial', 'Non-Industrial', 'Pre-modern'])
plt.title('Edge Weights by Source to Database Viruses')
plt.ylabel('Edge Weight')
plt.xlabel('Source Time Label')
plt.tight_layout()

plot_filename = 'Viral_Genomics/outputs/violinplot_edge_weights_to_database_by_time.png'
plt.savefig(plot_filename)

#%%
mean_weight_pre = np.mean(pre_only_combined['edge_weight'])
std_weight_pre = np.std(pre_only_combined['edge_weight'])

#%%
one = mean_weight_pre + std_weight_pre
two = mean_weight_pre + 2 * std_weight_pre
three = mean_weight_pre + 3 * std_weight_pre
four = mean_weight_pre + 4 * std_weight_pre

#%%
plt.figure(figsize=(10, 6))
plt.hist(pre_only_combined['edge_weight'], bins=30, edgecolor='k', color='lightsteelblue')
plt.yscale('log')

plt.axvline(mean_weight_pre, color='green', linestyle='dashed', linewidth=2, label='Mean')
plt.axvline(one, color='orange', linestyle='dashed', linewidth=2, label='Mean + 1 Std')
plt.axvline(two, color='yellow', linestyle='dashed', linewidth=2, label='Mean + 2 Std')
plt.axvline(three, color='purple', linestyle='dashed', linewidth=2, label='Mean + 3 Std')
plt.axvline(four, color='red', linestyle='dashed', linewidth=2, label='Mean + 4 Std')

plt.legend()

plt.title('Distribution of Pre-modern Edge Weights')
plt.xlabel('Edge Weight')
plt.ylabel('Frequency (log-scaled)')
plt.savefig('Viral_Genomics/outputs/histogram_edge_weights_pre_only_log_scale_with_cutoffs.png')