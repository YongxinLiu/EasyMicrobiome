import os
import pandas as pd

input_dir = './result/Emu/gtdb_cleaned/'
output_dir = './result/Emu/gtdb_cleaned_noprefix/'
os.makedirs(output_dir, exist_ok=True)

prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
taxonomy_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def remove_prefix(value):
    if pd.isna(value):
        return value
    for prefix in prefixes:
        if value.startswith(prefix):
            return value[len(prefix):]
    return value

for file in os.listdir(input_dir):
    if file.endswith('.tsv'):
        file_path = os.path.join(input_dir, file)
        df = pd.read_csv(file_path, sep='\t', dtype=str)
        
        for rank in taxonomy_ranks:
            if rank in df.columns:
                df[rank] = df[rank].apply(remove_prefix)
        
        output_path = os.path.join(output_dir, file)
        df.to_csv(output_path, sep='\t', index=False)
        print(f"✅ Prefix removed and saved: {file}")

print("All files processed successfully.")
