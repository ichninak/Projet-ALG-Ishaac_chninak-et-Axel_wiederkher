#!/usr/bin/env python3
import time
import pickle

def read_fasta(path):
    """Lit un FASTA (une ou plusieurs lignes) et renvoie la séquence collée."""
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq_parts.append(line.strip())
    return "".join(seq_parts)

def get_kmers(sequence, k):
    """Générateur simple de k-mers (sliding window)."""
    n = len(sequence)
    for i in range(n - k + 1):
        yield sequence[i:i+k]

def build_index(list_file_path, k, output_file):
    """
    Construit l'index naïf : dict{kmer -> set(color_id)} et le sérialise en pickle.
    Imprime les logs conformes V2 :
      OUT TIME_BUILD: <float>
      OUT TIME_SERIALISATION: <float>
    """
    t0 = time.time()

    graph = {}  # kmer -> set(color_id)

    # lire la liste des chemins (relatifs au répertoire courant d'exécution)
    with open(list_file_path, "r") as fh:
        genome_paths = [ln.strip() for ln in fh if ln.strip()]

    # build
    start_build = time.time()
    for color_id, path in enumerate(genome_paths):
        seq = read_fasta(path)
        for kmer in get_kmers(seq, k):
            if kmer in graph:
                graph[kmer].add(color_id)
            else:
                graph[kmer] = {color_id}
    end_build = time.time()

    # sérialisation
    start_ser = time.time()
    with open(output_file, "wb") as out:
        pickle.dump(graph, out, protocol=pickle.HIGHEST_PROTOCOL)
    end_ser = time.time()

    print(f"OUT TIME_BUILD: {end_build - start_build}")
    print(f"OUT TIME_SERIALISATION: {end_ser - start_ser}")

if __name__ == "__main__":
    print("Ce fichier ne doit pas être exécuté directement. Utilisez dbg_indexer.py")
