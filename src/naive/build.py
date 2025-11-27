#!/usr/bin/env python3
import time
import pickle

def read_fasta(path):
    seq = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)

def get_kmers(sequence, k):
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i+k]

def build_index(list_file_path, k, output_file):
    start_build = time.time()
    graph = {}

    with open(list_file_path, "r") as file_list:
        genomes_paths = [line.strip() for line in file_list if line.strip()]

    for color_id, genome_path in enumerate(genomes_paths):
        seq = read_fasta(genome_path)

        for kmer in get_kmers(seq, k):
            if kmer not in graph:
                graph[kmer] = set()
            graph[kmer].add(color_id)

    end_build = time.time()

    start_ser = time.time()
    with open(output_file, "wb") as out:
        pickle.dump(graph, out)
    end_ser = time.time()

    print(f"OUT TIME_BUILD: {end_build - start_build}")
    print(f"OUT TIME_SERIALISATION: {end_ser - start_ser}")

if __name__ == "__main__":
    print("Ce fichier ne doit pas être appelé directement. Utilisez dbg_indexer.py")
