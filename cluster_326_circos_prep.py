#%%
import pandas as pd
import numpy as np

#%%
heatmap_df = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/big_dir_cluster_326_0/heatmap_dataframe.csv')
pc_interest = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/big_dir_cluster_326_0/pc_interest.csv')
lifestyle_time = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/cleaned_genome_to_time_lifestyle.csv')

#%%
# Define a custom sorting key
def sorting_key(chr_id):
    if chr_id.startswith('ind-'):
        return 0
    elif chr_id.startswith('pre-'):
        return 1
    elif chr_id.startswith('pal-'):
        return 2
    return 3

# Sort the DataFrame based on the custom key
heatmap_df = heatmap_df.sort_values(by='Unnamed: 0', key=lambda x: x.map(sorting_key))
#%%

for i, row in heatmap_df.iterrows():
    heatmap_df.at[i,'Unnamed: 0'] = row['Unnamed: 0'].split('(')[0]


#%%
heatmap_df.rename(columns={'Unnamed: 0':'chr_id'}, inplace=True)

#%%
karyotype_file = pd.DataFrame(columns=['chr', 'chr_type', 'chr_id', 'label', 'start', 'end', 'color'])
#colored by lifestyle
#%%
karyotype_file['label'] = heatmap_df.iloc[:,0].astype(str)
karyotype_file['chr'] = 'chr'
karyotype_file['chr_type'] = '-'
karyotype_file['chr_id'] =  heatmap_df['chr_id']
karyotype_file['start'] = 0
karyotype_file['end'] = np.shape(heatmap_df)[1] - 4
for i, row in heatmap_df.iterrows():
    karyotype_index = row['chr_id']
    if row['Lifestyle'] == 'Lytic':
        karyotype_file.loc[karyotype_file['chr_id'] == karyotype_index, 'color'] = 'purple'
    elif row['Lifestyle'] == 'Temperate':
        karyotype_file.loc[karyotype_file['chr_id'] == karyotype_index, 'color'] = 'orange'

#%%
#colored by time period
karyotype_file_time = pd.DataFrame(columns=['chr', 'chr_type', 'chr_id', 'label', 'start', 'end', 'color'])
#colored by lifestyle
#%%
karyotype_file_time['label'] = heatmap_df.iloc[:,0].astype(str)
karyotype_file_time['chr'] = 'chr'
karyotype_file_time['chr_type'] = '-'
karyotype_file_time['chr_id'] =  heatmap_df['chr_id']
karyotype_file_time['start'] = 0
karyotype_file_time['end'] = np.shape(heatmap_df)[1] - 4
for i, row in heatmap_df.iterrows():
    karyotype_index = row['chr_id']
    if row['Time Period'] == 'Industrial':
        karyotype_file_time.loc[karyotype_file_time['chr_id'] == karyotype_index,'color'] = 'red'
    elif row['Time Period'] == 'Pre-Industrial':
        karyotype_file_time.loc[karyotype_file_time['chr_id'] == karyotype_index,'color'] = 'blue'
    elif row['Time Period'] == 'Paleo':
        karyotype_file_time.loc[karyotype_file_time['chr_id'] == karyotype_index,'color'] = 'green'

#%%
for pc in range(1, heatmap_df.shape[1]):
    pc_column = heatmap_df.columns[pc]
    if pc_column in pc_interest['pc_id'].values:
        for chr_index in range(len(heatmap_df)):
            if heatmap_df.iloc[chr_index, pc] == 1:
                heatmap_df.at[chr_index, pc_column] = 0.5

#%%
heatmap_df.drop(columns=['Lifestyle', 'Time_Period', 'Time Period'], inplace=True)

heatmap_df_fixed = pd.DataFrame(columns = ['chr', 'start', 'end', 'value', 'options'])
chrs = []
starts = []
ends = []
values = []
for chr in range(0,np.shape(heatmap_df)[0]):
    for pc in range(1,np.shape(heatmap_df)[1]):
        chrs.append(heatmap_df.iloc[chr,0])
        starts.append(pc)
        ends.append(pc+1)
        values.append(heatmap_df.iloc[chr,pc])
heatmap_df_fixed['chr'] = chrs
heatmap_df_fixed['start'] = starts
heatmap_df_fixed['end'] = ends
heatmap_df_fixed['value'] = values

#%%
links = pd.DataFrame(columns=['Source', 'Source_start', 'Source_end', 'Target', 'Target_start', 'Target_end'])
source = []
source_start = []
source_end = []
target = []
target_start = []
target_end = []
for pc in range(1,heatmap_df.shape[1]):
    for i in range(0,heatmap_df.shape[0] - 1):
        for j in range(1,heatmap_df.shape[0]):
            if i != j:
                if (heatmap_df.iloc[i,pc] == 1 and heatmap_df.iloc[j,pc] == 1) or (heatmap_df.iloc[i,pc] == 0.5 and heatmap_df.iloc[j,pc] == 0.5):
                    source.append(heatmap_df.iloc[i,0])
                    source_start.append(pc)
                    source_end.append(pc+1)
                    target.append(heatmap_df.iloc[j,0])
                    target_start.append(pc)
                    target_end.append(pc+1)

#%%
links['Source'] = source
links['Source_start'] = source_start
links['Source_end'] = source_end
links['Target'] = target
links['Target_start'] = target_start
links['Target_end'] = target_end


#%%
vf_links = pd.DataFrame(columns=['Source', 'Source_start', 'Source_end', 'Target', 'Target_start', 'Target_end'])
source = []
source_start = []
source_end = []
target = []
target_start = []
target_end = []

for pc in range(1,heatmap_df.shape[1]):
    for i in range(0,heatmap_df.shape[0] - 1):
        for j in range(1,heatmap_df.shape[0]):
            if i != j:
                if (heatmap_df.iloc[i,pc] == 0.5 and heatmap_df.iloc[j,pc] == 0.5):
                    source.append(heatmap_df.iloc[i,0])
                    source_start.append(pc)
                    source_end.append(pc+1)
                    target.append(heatmap_df.iloc[j,0])
                    target_start.append(pc)
                    target_end.append(pc+1)

vf_links['Source'] = source
vf_links['Source_start'] = source_start
vf_links['Source_end'] = source_end
vf_links['Target'] = target
vf_links['Target_start'] = target_start
vf_links['Target_end'] = target_end


#%%
karyotype_file_time.to_csv('/Users/madelaineleitman/Downloads/KnowlesLab/big_dir_cluster_326_0/karyotype_time.txt', sep='\t', index=False, header=False)
heatmap_df_fixed.to_csv('/Users/madelaineleitman/Downloads/KnowlesLab/big_dir_cluster_326_0/heatmap.txt', sep='\t', index=False, header=False)
links.to_csv('/Users/madelaineleitman/Downloads/KnowlesLab/big_dir_cluster_326_0/links.txt', sep='\t', index=False, header=False)
vf_links.to_csv('/Users/madelaineleitman/Downloads/KnowlesLab/big_dir_cluster_326_0/vf_links.txt', sep='\t', index=False, header=False)