#!/usr/bin/env python3
import pickle
import time

def read_fasta_queries(path):
    queries = []
    header = None
    seq = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    queries.append((header, "".join(seq)))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            queries.append((header, "".join(seq)))
    return queries

def get_kmers(sequence, k):
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i+k]

def query_index(index_file, fasta_queries, k, output_file):
    graph = pickle.load(open(index_file, "rb"))
    nb_colors = 1 + max(max(colors) for colors in graph.values())

    queries = read_fasta_queries(fasta_queries)

    with open(output_file, "w") as out:
        for header, seq in queries:
            kmers = list(get_kmers(seq, k))
            scores = [0] * nb_colors
            for kmer in kmers:
                if kmer in graph:
                    for c in graph[kmer]:
                        scores[c] += 1
            total = len(kmers)
            ratios = [s / total for s in scores]
            out.write(header + "\t" + "\t".join(map(str, ratios)) + "\n")

if __name__ == "__main__":
    print("Ce fichier ne doit pas être appelé directement. Utilisez dbg_indexer.py")
