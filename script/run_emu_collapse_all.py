import os
import subprocess
import sys
import time  # Small delay for WSL file sync if needed

# Define the base path (WSL format)
base_path = "/mnt/d/Amplicon/results_emu2"

# Subfolders containing EMU output
db_folders = ["gtdb_cleaned"]

# Taxonomic ranks to collapse to
ranks = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]

print(f"Starting EMU taxonomy collapsing process from: {base_path}")
print(f"Processing folders: {', '.join(db_folders)}")
print(f"Collapsing to ranks: {', '.join(ranks)}")
print("-" * 50)

for db in db_folders:
    folder_path = os.path.join(base_path, db)

    # Check that the folder exists
    if not os.path.isdir(folder_path):
        print(f"⚠️ Warning: Folder not found: {folder_path}. Skipping.")
        continue

    print(f"\n📁 Processing folder: {folder_path}")

    # Find input files
    for file in os.listdir(folder_path):
        if file.endswith("_rel-abundance-renamed.tsv"):
            print(f"\n🔍 Found file: {file}")
            original_basename = os.path.splitext(file)[0]  # e.g., HS_1.filtered...

            for rank in ranks:
                print(f"    📦 Collapsing to: {rank}")

                expected_output_file = f"{original_basename}-{rank}.tsv"
                expected_output_path = os.path.join(folder_path, expected_output_file)

                # Run EMU
                command = [
                    "emu",
                    "collapse-taxonomy",
                    file,
                    rank
                ]

                try:
                    result = subprocess.run(
                        command,
                        cwd=folder_path,
                        capture_output=True,
                        text=True,
                        check=False
                    )

                    time.sleep(0.5)  # WSL sometimes needs a slight delay for file I/O

                    if result.returncode == 0:
                        if os.path.exists(expected_output_path):
                            print(f"    ✅ Done: Output saved as '{expected_output_file}'")
                        else:
                            print(f"    ⚠️ EMU reported success, but output file not found: {expected_output_file}")
                            print(f"       Stdout: {result.stdout.strip()}")
                            print(f"       Stderr: {result.stderr.strip()}")
                    else:
                        print(f"    ❌ EMU failed for '{file}' to '{rank}'")
                        print(f"       Stdout: {result.stdout.strip()}")
                        print(f"       Stderr: {result.stderr.strip()}")

                except FileNotFoundError:
                    print(f"    ❌ Error: 'emu' command not found. Please check your PATH or installation.")
                    sys.exit(1)
                except Exception as e:
                    print(f"    ❌ Unexpected error: {e}")

print("\n" + "-" * 50)
print("🎉 EMU taxonomy collapsing process completed.")
