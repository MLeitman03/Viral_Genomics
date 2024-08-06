#%%
#Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import MDS
from Bio import SeqIO
from scipy.stats import f_oneway
import numpy as np
from scipy.stats import stats

#%%
lifestyle_time = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/cleaned_genome_to_time_lifestyle.csv')

column_names = ['query_id', 'virulence_factor', 'percent_identity', "alignment_length", "mismatches",
           "gap_opens", "q_start", "q_end", "s_start", "s_end", "e_value", "bit_score"]
blast = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/blast_results_predicted_proteins.txt', sep ='\t', header=None, names=column_names)

def read_fasta_to_dataframe(fasta_file):
    records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        records.append({
            'sequence_id': record.id,
            'description': record.description,
            'sequence': str(record.seq)
        })
    df = pd.DataFrame(records)
    return df

vfdb_df = read_fasta_to_dataframe('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/VFDB_setB_pro.fas')


# Replace 'predicted_genes.fna' with the path to your .fna file
file_path = '/Users/madelaineleitman/Downloads/KnowlesLab/prodigal_output/predicted_genes.fna'

# Read the .fna file
sequences = list(SeqIO.parse(file_path, 'fasta'))

# Create lists to store IDs and sequences
ids = []
seqs = []

# Iterate over each sequence record and store the ID and sequence
for seq_record in sequences:
    ids.append(seq_record.id)
    seqs.append(str(seq_record.seq))

# Create a DataFrame
genes = pd.DataFrame({'ID': ids, 'Sequence': seqs})

#%%
#merge blast output with virulence factor database
df = pd.merge(blast, vfdb_df, left_on = 'virulence_factor', right_on = 'sequence_id', how = 'left')

#%%

df['Genome'] = df['query_id'].apply(lambda x : x.split('_')[:-1])
df['Genome'] = df['Genome'].apply(lambda x : '_'.join(x))


#%%
merged = pd.merge(df, lifestyle_time, on = 'Genome', how = 'left')


#%%
unique_counts = merged.groupby('Genome')['query_id'].nunique().reset_index()

unique_counts.columns = ['Genome', 'Unique_VF']

#%%
#no_db = lifestyle_time[lifestyle_time['Genome'].str.contains('pal-|ind-|pre-')]
no_db = lifestyle_time

#%%
all_genomes = no_db.loc[:,['Genome', 'Time Period', 'Lifestyle']]

#%%
all_genomes = pd.merge(all_genomes, unique_counts, on = 'Genome', how = 'left').fillna(0)


#%%
genes['Genome'] = genes['ID'].apply(lambda x : x.split('_')[:-1])
genes['Genome'] = genes['Genome'].apply(lambda x : '_'.join(x))


#%%
genes_per_genome = genes.groupby('Genome')['ID'].nunique().reset_index()

genes_per_genome.columns = ['Genome', 'Number of Genes']

#%%
all_genomes = pd.merge(all_genomes, genes_per_genome, on = 'Genome', how = 'left')

#%%
all_genomes['VF per gene'] = (all_genomes['Unique_VF'] / all_genomes['Number of Genes'])
#%%
# Drop rows where 'Time Period' is 'Database'
all_genomes = all_genomes[all_genomes['Time Period'] != 'Database']

#%%
plt.figure(figsize=(12,8))
sns.violinplot(x='Lifestyle', y = 'VF per gene', data = all_genomes, edgecolor = 'black', hue = 'Time Period')
plt.xlabel('Lifestyle')
plt.ylabel('Number of Virulence Factors per Gene in Genome')
plt.title('Number of Virulence Factors per Gene in Genome by Lifestyle')
plt.grid(True)
plt.show()


#%%

plt.figure(figsize=(12,8))
sns.violinplot(x='Time Period', y = 'VF per gene', data = all_genomes, edgecolor = 'black', cut=0)
plt.xlabel('Time Period')
plt.ylabel('Number of Virulence Factors per ORF in Genome')
plt.title('Number of Virulence Factors per ORF in Genome by Time Period')
plt.grid(True)
plt.show()

#%%
all_genomes = all_genomes[all_genomes['Lifestyle'] != 0]
all_genomes['VF per 100 ORFs'] = all_genomes['VF per gene'] * 100
plt.figure(figsize=(10, 6))
palette = sns.color_palette("Accent", 3)
sns.violinplot(x='Time Period', y='VF per 100 ORFs', hue='Lifestyle', data=all_genomes, split=True, palette=palette[1:3], order = ['Industrial', 'Non-Industrial', 'Pre-Modern'])
plt.title('Percentage of ORFs that are VFs by Time Period')
plt.ylabel('% of ORFs that are VFs')

#plt.show()
#plt.savefig('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/outputs/VF_per_100_ORF_violin.png')

#%%
industrial = all_genomes[all_genomes['Time Period'] == 'Industrial']
ind_lytic = industrial.loc[industrial['Lifestyle'] == 'Lytic', 'VF per 100 ORFs']
ind_temp = industrial.loc[industrial['Lifestyle'] == 'Temperate', 'VF per 100 ORFs']

non_ind = all_genomes[all_genomes['Time Period'] == 'Non-Industrial']
non_lytic = non_ind.loc[non_ind['Lifestyle'] == 'Lytic', 'VF per 100 ORFs']
non_temp = non_ind.loc[non_ind['Lifestyle'] == 'Temperate', 'VF per 100 ORFs']

pre = all_genomes[all_genomes['Time Period'] == 'Pre-Modern']
pre_lytic = pre.loc[pre['Lifestyle'] == 'Lytic', 'VF per 100 ORFs']
pre_temp = pre.loc[pre['Lifestyle'] == 'Temperate', 'VF per 100 ORFs']

#%%
from scipy.stats import ttest_ind

ttest_result = ttest_ind(ind_lytic, ind_temp)
print("t-Test results:", ttest_result)

#%%
from scipy.stats import mannwhitneyu

mannwhitney_result = mannwhitneyu(pre_lytic, pre_temp)
print("Mann-Whitney U Test results:", mannwhitney_result)

#%%
from scipy.stats import ks_2samp

ks_result = ks_2samp(non_lytic, non_temp)
print("Kolmogorov-Smirnov Test results:", ks_result)

#%%

vf_presence_absence_df = pd.DataFrame(index = df['query_id'].unique(), columns = df['description'].unique())

#creates dataframe where the rows are the contig (virus) names and the columns are the virulence factors

#%%
#only contains ind, non, and pre samples
vf_presence_absence_df = vf_presence_absence_df[vf_presence_absence_df.index.str.contains('pal-|ind-|pre-')]

#%%
#sort in order for easy finding
vf_presence_absence_df = vf_presence_absence_df.sort_index(axis=1)
vf_presence_absence_df = vf_presence_absence_df.sort_index()

#%%
#adds 0 where there is no protein cluster for that genome and adds 1 where there is that protein cluster present
vf_presence_absence_df = vf_presence_absence_df.fillna(0)
for _, row in df.iterrows():
    genome = row['query_id']
    vf = row['description']
    if genome in vf_presence_absence_df.index and vf in vf_presence_absence_df.columns:
        vf_presence_absence_df.at[genome, vf] = 1

#%%
#only with contig_id/genome, time and lifestyle
metadata = lifestyle_time.iloc[:,[2,3,4]]

#%%
#put genomes in order
metadata = metadata[metadata['Genome'].str.contains('pal-|ind-|pre-')]
metadata = metadata.sort_values(by='Genome')

#%%
# Apply MDS
mds = MDS(n_components=2, metric = True, random_state=42)
embedding = mds.fit_transform(vf_presence_absence_df)

#%%
embedding_df_lifestyle = pd.DataFrame(embedding, columns=['MDS1', 'MDS2'])
embedding_df_lifestyle['metadata'] = metadata['Lifestyle']

#%%
embedding_df_time = pd.DataFrame(embedding, columns=['MDS1', 'MDS2'])
embedding_df_time['metadata'] = metadata['Time Period']


#%%
#plot MDS projection colored by lifestyle
plt.figure(figsize=(15,15))
sns.scatterplot(data=embedding_df_lifestyle, x='MDS1', y='MDS2', hue='metadata', palette='icefire', alpha = 0.7)
plt.title('MDS Projection of Viral Genomes by VFs')
plt.xlabel('MDS1')
plt.ylabel('MDS2')
plt.legend(title='Lifestyle')
plot_filename = 'Viral_Genomics/outputs/mds_vfs_by_lifestyle.png'
plt.savefig(plot_filename)

#%%
#plot MDS projection colored by Time period

plt.figure(figsize=(15,15))
sns.scatterplot(data=embedding_df_time, x='MDS1', y='MDS2', hue='metadata', palette='icefire', alpha = 0.7)
plt.title('MDS Projection of Viral Genomes by VFs')
plt.xlabel('MDS1')
plt.ylabel('MDS2')
plt.legend(title='Time Period')
plot_filename = 'Viral_Genomics/outputs/mds_vfs_by_time.png'
plt.savefig(plot_filename)

#run vcontact2 only on contig virulence factors to make .ntw file then put into cytoscape

#%%
filtered_df = df.loc[df.groupby('query_id')['bit_score'].idxmax()]
#only keep the highest bit score for each orf

#%%
# Function to write FASTA file
def write_fasta(df, output_file):
    with open(output_file, 'w') as f:
        for index, row in df.iterrows():
            f.write(f">{row['description']}")
            f.write(f"{row['sequence']}\n")

#%%

output_file = '/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/vf_sequences.fasta'
# Write the DataFrame to a FASTA file
write_fasta(filtered_df, output_file)

#%%
# Initialize the gene_to_genome DataFrame
gene_to_genome = pd.DataFrame(columns=['protein_id', 'contig_id', 'keywords'])

# List to collect rows for the new DataFrame
rows = []

# Extract relevant columns from filtered_df and append to the rows list
for index, row in filtered_df.iterrows():
    protein_id = row['description']
    contig_id = row['Genome']
    keywords = row['query_id']  # Replace with actual keyword extraction logic if needed

    # Create a dictionary for the row and append it to the rows list
    rows.append({'protein_id': protein_id, 'contig_id': contig_id, 'keywords': keywords})

# Convert the list of rows to a DataFrame
new_rows_df = pd.DataFrame(rows)

# Concatenate the new rows DataFrame with the existing gene_to_genome DataFrame
gene_to_genome = pd.concat([gene_to_genome, new_rows_df], ignore_index=True)

# Save the gene_to_genome DataFrame to a TSV file
gene_to_genome.to_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/gene_to_genome_mapping.tsv', sep='\t', index=False)

