#!/usr/bin/env python3
import time
import pickle


def read_fasta_queries(path):
    """
    Lit un fichier FASTA contenant les séquences requêtes
    et retourne une liste de tuples (header, sequence).
    """
    queries = []
    header = None
    seq_parts = []

    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    queries.append((header, "".join(seq_parts)))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip().upper())

        if header is not None:
            queries.append((header, "".join(seq_parts)))

    return queries


def get_kmers(seq, k):
    """
    Générateur de k-mers pour une séquence requête.
    """
    n = len(seq)
    for i in range(n - k + 1):
        yield seq[i:i+k]


def query_index(index_file, query_file, k, output_file):
    """
    Interroge un graphe de De Bruijn coloré compacté.

    Les k-mers des requêtes sont associés aux unitigs,
    puis aux couleurs correspondantes afin de calculer
    les ratios de similarité avec chaque génome.
    """
    t_deser_start = time.time()
    with open(index_file, "rb") as f:
        data = pickle.load(f)
    t_deser_end = time.time()

    print(f"OUT TIME_DESERIALISATION: {t_deser_end - t_deser_start}")

    unitigs = data["unitigs"]
    unitig_colors = data["unitig_colors"]
    kmer_to_unitig = data["kmer_to_unitig"]

    nb_colors = max((max(c) for c in unitig_colors if c), default=-1) + 1

    t_query_start = time.time()
    queries = read_fasta_queries(query_file)

    with open(output_file, "w") as out:
        for header, seq in queries:
            scores = [0] * nb_colors
            kmers = list(get_kmers(seq, k))
            total = len(kmers)

            for kmer in kmers:
                if kmer in kmer_to_unitig:
                    uid = kmer_to_unitig[kmer]
                    for c in unitig_colors[uid]:
                        scores[c] += 1

            if total > 0:
                ratios = [s / total for s in scores]
            else:
                ratios = [0.0] * nb_colors

            out.write(header + "\t" + "\t".join(f"{r:.4f}" for r in ratios) + "\n")

    t_query_end = time.time()
    print(f"OUT TIME_QUERY: {t_query_end - t_query_start}")


if __name__ == "__main__":
    print("Ce fichier ne doit pas être exécuté directement. Utilisez dbg_indexer.py")
