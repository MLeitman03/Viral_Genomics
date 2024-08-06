#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#genome length by lifestyle and time period

#%%
lifestyle_time = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/cleaned_genome_to_time_lifestyle.csv')

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

genomes = read_fasta_to_dataframe('/Users/madelaineleitman/Downloads/KnowlesLab/all_genomes_singleline.fa')

#%%
df = pd.merge(genomes, lifestyle_time.iloc[:,2:], left_on = 'sequence_id', right_on = 'Genome')

#%%
df['sequence_length'] = df['sequence'].apply(lambda x: len(x))