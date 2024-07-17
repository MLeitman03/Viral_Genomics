from sys import argv
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
from collections import defaultdict
import csv

# Read the CSV files
data = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/genomes_to_pcs.csv') # Lists which protein clusters are located within each genome
file = argv[1] # This is the input file that the user inputs, it says the genome names within each cluster

vf_df = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/blastp_to_VFDB.csv')
lifestyle_labels = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/network_labels.csv')

# Turn the input file into a list of contigs
cluster_contigs = pd.read_table(file, header=None)[0].tolist()

# Filter the data for relevant contigs
only_cluster_contigs_with_pcid = data[data['contig_id'].isin(cluster_contigs)]
only_cluster_contigs_with_pcid = only_cluster_contigs_with_pcid.dropna()

# Create a DataFrame for the heatmap
heatmap_df = pd.DataFrame(index=cluster_contigs, columns=only_cluster_contigs_with_pcid['pc_id'].unique())
heatmap_df = heatmap_df[heatmap_df.index.str.contains('pal-|ind-|pre-')]
heatmap_df = heatmap_df.fillna(0)

# Populate the heatmap DataFrame
PC_choice = []
for _, row in only_cluster_contigs_with_pcid.iterrows():
    genome = row['contig_id']
    protein_cluster = row['pc_id']
    PC_choice.append(protein_cluster)
    if genome in heatmap_df.index and protein_cluster in heatmap_df.columns:
        heatmap_df.at[genome, protein_cluster] = 1

# Create a dictionary for lifestyle labels
genome_names = heatmap_df.index
lifestyle = []
for index, row in lifestyle_labels.iterrows():
    if row.iloc[0] in genome_names:
        lifestyle.append(row.iloc[16])
dict_of_lifestyle_to_genome = dict(zip(genome_names, lifestyle))

# Update row names of the heatmap DataFrame to include lifestyle
new_row_names = [f"{key}({value})" for key, value in dict_of_lifestyle_to_genome.items()]
heatmap_df.index = new_row_names

# Sort heatmap DataFrame
lytic_rows = heatmap_df[heatmap_df.index.str.contains('Lytic')]
temperate_rows = heatmap_df[heatmap_df.index.str.contains('Temperate')]
sorted_heatmap_df = pd.concat([lytic_rows, temperate_rows])

# Plot and save the heatmap
plt.figure(figsize=(30, 10))
heatmap = sns.heatmap(sorted_heatmap_df)
target_directory = os.path.dirname(file)
file_path_heatmap = os.path.join(target_directory, "heatmap_of_genomes_pc.png")
plt.savefig(file_path_heatmap)

# Add lifestyle and time period columns
sorted_heatmap_df = sorted_heatmap_df.assign(Lifestyle='', Time_Period='')
sorted_heatmap_df.loc[sorted_heatmap_df.index.str.contains('Lytic'), 'Lifestyle'] = 'Lytic'
sorted_heatmap_df.loc[sorted_heatmap_df.index.str.contains('Temperate'), 'Lifestyle'] = 'Temperate'
sorted_heatmap_df.loc[sorted_heatmap_df.index.str.contains('pre'), 'Time Period'] = 'Pre-Industrial'
sorted_heatmap_df.loc[sorted_heatmap_df.index.str.contains('ind'), 'Time Period'] = 'Industrial'
sorted_heatmap_df.loc[sorted_heatmap_df.index.str.contains('pal'), 'Time Period'] = 'Paleo'


# Save the protein cluster of interest
mydict = defaultdict(list)
for index, row in vf_df.iterrows():
    if row['query_id'][-9:] in PC_choice:
        contig = row['query_id'][:-10]
        description = row['description']
        seq = row['sequence']
        e_value = row['e_value']
        bit_score = row['bit_score']
        pc_id = row['query_id'][-9:]
        mydict[pc_id].append({
            'contig_id': contig,
            'description': description,
            'sequence': seq,
            'e_value': e_value,
            'bit_score': bit_score
        })

# Convert the dictionary to a list of dictionaries for CSV writing
csv_data = []
for k, v in mydict.items():
    for item in v:
        row = {'pc_id': k}
        row.update(item)
        csv_data.append(row)

# Define  CSV file path
file_path_pc_interest = os.path.join(target_directory, "pc_interest.csv")

# Write the data to a CSV file
with open(file_path_pc_interest, 'w', newline='') as csvfile:
    fieldnames = ['pc_id', 'contig_id', 'description', 'sequence', 'e_value', 'bit_score']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(csv_data)

# Save the sorted heatmap DataFrame
file_path_heatmap_df = os.path.join(target_directory, "heatmap_dataframe.csv")
sorted_heatmap_df.to_csv(file_path_heatmap_df, index=True)

print('Saved files to ' + target_directory)
