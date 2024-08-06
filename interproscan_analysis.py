#%%
import pandas as pd
import numpy as np

#%%
df = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/pcs_interproscan_again.tsv', sep = '\t', header = None)

my_dict = {}

# Iterate through the DataFrame rows
for index, row in df.iterrows():
    key = row[0]
    value = row[12]

    # If the key is not in the dictionary, add it with an empty list
    if key not in my_dict:
        my_dict[key] = []

    # Append the value to the list for this key
    my_dict[key].append(value)

#%%
keys_with_dash = [key for key, values in my_dict.items() if '-' in values]

#%%
df[df[12] == '-']

#%%
num_of_proteins = []
for key, values in my_dict.items():
    num_of_proteins.append(len(values))