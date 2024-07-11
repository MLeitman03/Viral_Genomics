
from sys import argv
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os

data = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/genomes_to_pcs.csv') #lists which protein clusters are located within each genome
file = argv[1] #this is the input file that the user inputs, it says the genome names within each cluster

pc_id_to_function = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/vConTACT_pcs.csv') #this is each protein cluster and its function

cluster_contigs = pd.read_table(file, header = None)[0].tolist() #turns the input file into a list of contigs

only_cluster_contigs_with_pcid = data[data['contig_id'].isin(cluster_contigs)] 

only_cluster_contigs_with_pcid = only_cluster_contigs_with_pcid.dropna()

heatmap_df = pd.DataFrame(index = cluster_contigs, columns = only_cluster_contigs_with_pcid['pc_id'].unique())
#creates dataframe where the rows are the contig (virus) names and the columns are the protein cluster ids in this cluster

heatmap_df = heatmap_df[heatmap_df.index.str.contains('pal|ind|pre')]
heatmap_df = heatmap_df.fillna(0)

#this puts a 1 in the cell where that contig contains that protein cluster
for _, row in only_cluster_contigs_with_pcid.iterrows(): 
    genome = row['contig_id']
    protein_cluster = row['pc_id']
    if genome in heatmap_df.index and protein_cluster in heatmap_df.columns:
        heatmap_df.at[genome, protein_cluster] = 1

#these grab the paleo contigs
target_indices = heatmap_df.index[heatmap_df.index.str.contains("pal")].tolist()
#indices that are paleo


PC_choice = []
for index, row in heatmap_df.iterrows():
    for col in range(len(heatmap_df.columns)):
        if index in target_indices and row[heatmap_df.columns[col]] == 0:
            PC_choice.append(heatmap_df.columns[col])

#this file has all contig names and the lifestyle that it is
lifestyle_labels = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/network_labels.csv')

genome_names = heatmap_df.index
lifestyle = []

#this loops through the huge file of contigs and their lifestyle and creates a dictionary of the contigs in this cluster with its lifestyle
for index, row in lifestyle_labels.iterrows():
        if row.iloc[0] in genome_names:
                lifestyle.append(row.iloc[16])
dict_of_lifestyle_to_genome = dict(zip(genome_names, lifestyle))

#changes the row names of the heatmap df to contain the lifestyle of the virus
new_row_names = [f"{key}({value})" for key, value in dict_of_lifestyle_to_genome.items()]

heatmap_df.index = new_row_names

#extract lytic and temperate rows separately
lytic_rows = heatmap_df[heatmap_df.index.str.contains('Lytic')]
temperate_rows = heatmap_df[heatmap_df.index.str.contains('Temperate')]
sorted_heatmap_df = pd.concat([lytic_rows, temperate_rows])
#moved lytic rows to the top

#plot heatmap
plt.figure(figsize = (30,10))
heatmap = sns.heatmap(sorted_heatmap_df)


target_directory = argv[1].strip().split('/')[0]
file_path_heatmap = os.path.join(target_directory, "heatmap_of_genomes_pc.png")

plt.savefig(file_path_heatmap)

sorted_heatmap_df = sorted_heatmap_df.assign(Lifestyle='', Time_Period='')

sorted_heatmap_df.loc[sorted_heatmap_df.index.str.contains('Lytic'), 'Lifestyle'] = 'Lytic'
sorted_heatmap_df.loc[sorted_heatmap_df.index.str.contains('Temperate'), 'Lifestyle'] = 'Temperate'

sorted_heatmap_df.loc[sorted_heatmap_df.index.str.contains('pre'), 'Time Period'] = 'Pre-Industrial'
sorted_heatmap_df.loc[sorted_heatmap_df.index.str.contains('ind'), 'Time Period'] = 'Industrial'
sorted_heatmap_df.loc[sorted_heatmap_df.index.str.contains('pal'), 'Time Period'] = 'Paleo'

#dictionary of function of protein clusters of interest (aka the ones that are paleo and not the others)
values = []
for index, row in pc_id_to_function.iterrows():
    if row.iloc[0] in PC_choice:
        values.append(row.iloc[3])
mydict = dict(zip(PC_choice, values))

content = '\n'.join(f"{k}: {v}" for k, v in mydict.items())
file_path_pc_interest = os.path.join(target_directory, "pc_interest.txt")
with open(file_path_pc_interest, 'w') as f:
	f.write(content)

file_path_heatmap_df = os.path.join(target_directory, "heatmap_dataframe.csv")
sorted_heatmap_df.to_csv(file_path_heatmap_df, index = True)      
print('Saved heatmap to ' + target_directory)
