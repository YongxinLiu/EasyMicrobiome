import os
import pandas as pd

input_dir = './result/Emu/silva/'
output_dir = './result/Emu/silva_cleaned/'
os.makedirs(output_dir, exist_ok=True)

# Taxonomic ranks
ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
columns_required = ['tax_id', 'abundance'] + ranks

# Helper to split lineage
def split_lineage(lineage):
    if pd.isna(lineage) or lineage.strip() == '':
        return pd.Series([None] * 7)
    parts = lineage.rstrip(';').split(';')
    parts += [None] * (7 - len(parts))  # pad if shorter
    return pd.Series(parts[:7])

for sample in os.listdir(input_dir):
    file_path = os.path.join(input_dir, sample, f"{sample}.filtered.qc_rel-abundance.tsv")
    if os.path.isfile(file_path):
        try:
            df = pd.read_csv(file_path, sep='\t', dtype=str)

            if 'lineage' not in df.columns and df.shape[1] == 3:
                df.columns = ['tax_id', 'abundance', 'lineage']
            elif 'lineage' not in df.columns and 'lineage' not in df.columns:
                print(f"⚠️ Skipped (no lineage column found): {sample}")
                continue

            # Split the lineage into 7 taxonomic ranks
            split_cols = df['lineage'].apply(split_lineage)
            split_cols.columns = ranks

            # Combine cleaned table
            df_cleaned = pd.concat([df[['tax_id', 'abundance']], split_cols], axis=1)

            # Ensure all columns exist
            for col in columns_required:
                if col not in df_cleaned.columns:
                    df_cleaned[col] = ''

            df_cleaned = df_cleaned[columns_required]

            # Write cleaned file
            output_file = os.path.join(output_dir, f"{sample}.filtered.qc_rel-abundance-renamed.tsv")
            df_cleaned.to_csv(output_file, sep='\t', index=False)
            print(f"✅ Cleaned and saved: {output_file}")
        except Exception as e:
            print(f"❌ Error in {sample}: {e}")
    else:
        print(f"❌ Missing file: {file_path}")
