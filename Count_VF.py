#%%
#Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import numpy as np

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
lifestyle_time['Time Period'] = lifestyle_time['Time Period'].astype(str)

# Correcting the 'Time Period' values in the DataFrame
lifestyle_time.loc[lifestyle_time['Time Period'] == 'Non_Industrial', 'Time Period'] = 'Non-Industrial'
lifestyle_time.loc[lifestyle_time['Time Period'] == 'Pre-modern', 'Time Period'] = 'Pre-Modern'

#%%
merged = pd.merge(df, lifestyle_time, on = 'Genome', how = 'left')


#%%
unique_counts = merged.groupby('Genome')['query_id'].nunique().reset_index()

unique_counts.columns = ['Genome', 'Unique_VF']

#%%
no_db = lifestyle_time[lifestyle_time['Genome'].str.contains('pal-|ind-|pre-')]

#%%
all_genomes = no_db.iloc[:,[1,2,3]]

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
plt.figure(figsize=(12,8))
sns.violinplot(x='Lifestyle', y = 'VF per gene', data = all_genomes, edgecolor = 'black', hue = 'Time Period')
plt.xlabel('Lifestyle')
plt.ylabel('Number of Virulence Factors per Gene in Genome')
plt.title('Number of Virulence Factors per Gene in Genome by Lifestyle')
plt.grid(True)
plt.show()


#%%
plt.figure(figsize=(12,8))
sns.violinplot(x='Time Period', y = 'VF per gene', data = all_genomes, edgecolor = 'black', hue = 'Lifestyle', cut=0, order = ['Industrial', 'Non-Industrial', 'Pre-Modern'])
plt.xlabel('Time Period')
plt.ylabel('Number of Virulence Factors per ORF in Genome')
plt.title('Number of Virulence Factors per ORF in Genome by Time Period')
plt.grid(True)
plt.show()

#%%
import scipy.stats as stats


# Extracting the data for each time period
industrial = all_genomes[all_genomes['Time Period'] == 'Industrial']['VF per gene']
non_industrial = all_genomes[all_genomes['Time Period'] == 'Non-Industrial']['VF per gene']
pre_modern = all_genomes[all_genomes['Time Period'] == 'Pre-Modern']['VF per gene']

# Perform Kruskal-Wallis Test
kruskal_result = stats.kruskal(industrial, non_industrial, pre_modern)
print(kruskal_result)

#%%

def remove_outliers(df, column):
    Q1 = df[column].quantile(0.25)
    Q3 = df[column].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    filtered_df = df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]
    return filtered_df

#%%
# Assuming your DataFrame is named 'all_genomes' and the column to filter is 'Unique_VF'

all_genomes_filtered = remove_outliers(all_genomes, 'VF per gene')

plt.figure(figsize=(10, 6))
palette = sns.color_palette("Accent", 3)
sns.violinplot(x='Time Period', y='VF per gene', hue='Lifestyle', data=all_genomes_filtered, split=True, palette=palette[1:3])
plt.title('Number of Virulence Factors per ORF in Genome by Time Period (Without Outliers)')

plt.title('Number of Virulence Factors per ORF in Genome by Time Period (Without Outliers)')
plt.savefig('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/outputs/VF_per_ORF_violin.png')

#%%

# Kruskal-Wallis Test
industrial_filtered = all_genomes_filtered[all_genomes_filtered['Time Period'] == 'Industrial']['VF per gene']
non_industrial_filtered = all_genomes_filtered[all_genomes_filtered['Time Period'] == 'Non-Industrial']['VF per gene']
pre_modern_filtered = all_genomes_filtered[all_genomes_filtered['Time Period'] == 'Pre-Modern']['VF per gene']

kruskal_result_filtered = stats.kruskal(industrial_filtered, non_industrial_filtered, pre_modern_filtered)
print("Kruskal-Wallis results on filtered data:")
print(kruskal_result_filtered)