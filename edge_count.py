#%%
# Import libraries
import pandas as pd
import numpy as np

#%%
# Upload data
edge_df = pd.read_table(filepath_or_buffer='/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/edge_weights.ntw', sep=' ', header=None)
network_labels_df = pd.read_csv(filepath_or_buffer='/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/network_labels.csv')


#%%
# Merge dataframe on Genome label
edge_df_genome_col = edge_df[[0]]
network_labels_df_genome_col = network_labels_df[['Genome']]
df = pd.merge(edge_df, network_labels_df, left_on=0, right_on='Genome', how='left')
df.rename(columns={0:'Source', 1:'Target', 2:'edge_weight'}, inplace=True)



