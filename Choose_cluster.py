#%%
import numpy as np
import pandas as pd
import os
import subprocess

#%%
lifestyle_time = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/cleaned_genome_to_time_lifestyle.csv')
pre = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/pre_modern_source_genomes.csv')
network_labels_df = pd.read_csv('/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/data/network_labels.csv')

#%%
# Create a mapping from 'Genome' to 'VC'
vc_to_genomes = network_labels_df.groupby('VC')['Genome'].apply(list).to_dict()
genome_to_lifestyle = network_labels_df.groupby('Genome')['Lifestyle'].apply(list).to_dict()

#%%
# Function to check the conditions
def check_time(genomes):
    #has_ind = any(genome.startswith('ind-') for genome in genomes)
    has_pre_or_ind = any(genome.startswith('pre-') or genome.startswith('ind-') for genome in genomes)
    has_pal = any(genome.startswith('pal-') for genome in genomes)
    return has_pre_or_ind and has_pal

def check_lifestyle_switch(genomes, genome_to_lifestyle):
    lifestyles = set()
    for genome in genomes:
        lifestyles.update(genome_to_lifestyle.get(genome, []))
    return len(lifestyles) > 1  # True if there are at least two different lifestyles

# Base directory to store the files
base_dir = "/Users/madelaineleitman/Downloads/KnowlesLab/Viral_Genomics/"


# Iterate through each VC and check the conditions and count the viruses
for vc, genomes in vc_to_genomes.items(): # clusters that have at least a pre modern and have a switch in lifestyle
    if check_time(genomes) and check_lifestyle_switch(genomes, genome_to_lifestyle):
        num_viruses = len(genomes)
        if 3 <= num_viruses <= 20:
            for genome in genomes:
                print(f"Genome: {genome}, Lifestyles: {genome_to_lifestyle.get(genome, [])}")
            # Create directory for the VC
            vc_dir = os.path.join(base_dir, vc)
            os.makedirs(vc_dir, exist_ok=True)

            # Write the genomes to a text file in the VC directory
            text_file_path = os.path.join(vc_dir, f"{vc}_genomes.txt")
            with open(text_file_path, "w") as file:
                for genome in genomes:
                    file.write(genome + "\n")

            print(f"{vc} contains {num_viruses} viruses. Genomes saved to {text_file_path}")

            # Call the protein_cluster_finder.py script with the text file as an argument
            subprocess.run(["python3", "protein_cluster_finder.py", text_file_path])

            # Paths for the input files for circos_prep.py
            heatmap_dataframe_path = os.path.join(vc_dir, "heatmap_dataframe.csv")
            pc_interest_path = os.path.join(vc_dir, "pc_interest.csv")

            # Ensure these files exist before proceeding
            if os.path.exists(heatmap_dataframe_path) and os.path.exists(pc_interest_path):
                # Call the circos_prep.py script with the required arguments
                subprocess.run(["python3", "Viral_Genomics/circos_prep.py", heatmap_dataframe_path, pc_interest_path, vc_dir])
#120, 149, 177_1, 233_0, 267_0, 289_0, 386_0, 396_0, 453_0, 755_0, 82