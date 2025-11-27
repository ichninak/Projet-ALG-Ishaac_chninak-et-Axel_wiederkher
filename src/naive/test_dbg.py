#!/usr/bin/env python3
import os
import subprocess

# --- Dossier temporaire pour les tests ---
test_dir = "test_tmp"
os.makedirs(test_dir, exist_ok=True)

# --- Création de fichiers FASTA fictifs ---
genome1 = ">genome1\nATGCGTACGTTAGC"
genome2 = ">genome2\nCGTACGTTAGCGTA"
query1 = ">query1\nATGCGTAC"
query2 = ">query2\nCGTTAGC"

# --- Fichiers génomes ---
genome_list_file = os.path.join(test_dir, "genomes_list.txt")
genome1_file = os.path.join(test_dir, "genome1.fasta")
genome2_file = os.path.join(test_dir, "genome2.fasta")

with open(genome1_file, "w") as f: f.write(genome1)
with open(genome2_file, "w") as f: f.write(genome2)
with open(genome_list_file, "w") as f:
    f.write(f"{genome1_file}\n{genome2_file}\n")

# --- Fichier de requêtes ---
query_file = os.path.join(test_dir, "queries.fasta")
with open(query_file, "w") as f:
    f.write(query1 + "\n" + query2 + "\n")

# --- Fichiers de sortie ---
index_file = os.path.join(test_dir, "index.pkl")
output_file = os.path.join(test_dir, "results.txt")

# --- Construction de l'index ---
subprocess.run([
    "python", "dbg_indexer.py", "build",
    "-i", genome_list_file,
    "-k", "3",
    "-o", index_file
], check=True)

# --- Interrogation de l'index ---
subprocess.run([
    "python", "dbg_indexer.py", "query",
    "-q", query_file,
    "-i", index_file,
    "-k", "3",
    "-o", output_file
], check=True)

# --- Affichage des résultats ---
print("\n--- Contenu du fichier de résultats ---")
with open(output_file, "r") as f:
    print(f.read())
