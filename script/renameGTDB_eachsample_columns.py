import pandas as pd
import os

input_dir = './result/Emu/gtdb/'
output_dir = './result/Emu/gtdb_cleaned/'

os.makedirs(output_dir, exist_ok=True)

rename_dict = {
    'd__Bacteria': 'superkingdom',
    'p__Pseudomonadota': 'phylum',
    'c__Gammaproteobacteria': 'class',
    'o__Enterobacterales': 'order',
    'f__Enterobacteriaceae': 'family',
    'g__Escherichia': 'genus',
    's__Escherichia': 'species'
}

for sample in os.listdir(input_dir):
    file_path = os.path.join(input_dir, sample, f"{sample}.filtered.qc_rel-abundance.tsv")
    if os.path.exists(file_path):
        df = pd.read_csv(file_path, sep='\t')
        df.rename(columns=rename_dict, inplace=True)
        df.to_csv(os.path.join(output_dir, f"{sample}.filtered.qc_rel-abundance-renamed.tsv"), sep='\t', index=False)

print("✅ GTDB samples renamed and saved.")
