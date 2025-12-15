#!/usr/bin/env python3
import time
import pickle


def read_fasta_queries(path):
    """
    Lit un fichier FASTA contenant plusieurs séquences requêtes.

    Retourne une liste de tuples (header, sequence), où :
    - header correspond à la ligne FASTA sans le caractère '>'
    - sequence est la séquence nucléotidique concaténée
    """
    queries = []
    header = None
    seq_parts = []

    with open(path, "r") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    queries.append((header, "".join(seq_parts)))
                header = line[1:].strip()
                seq_parts = []
            else:
                if line.strip():
                    seq_parts.append(line.strip().upper())

        if header is not None:
            queries.append((header, "".join(seq_parts)))

    return queries


def get_kmers(seq, k):
    """
    Générateur de k-mers pour une séquence donnée.
    """
    n = len(seq)
    if n < k:
        return
    for i in range(n - k + 1):
        yield seq[i:i+k]


def query_index(index_file, fasta_queries, k, output_file):
    """
    Interroge un graphe de De Bruijn coloré naïf.

    Pour chaque séquence requête, on calcule le ratio de k-mers
    partagés avec chacun des génomes de référence.

    Les résultats sont écrits dans un fichier texte :
    header <TAB> ratio_genome1 <TAB> ratio_genome2 ...

    Affiche les temps d'exécution au format imposé :
    - OUT TIME_DESERIALISATION
    - OUT TIME_QUERY
    """
    # Désérialisation du graphe
    t_deser_start = time.time()
    with open(index_file, "rb") as fh:
        graph = pickle.load(fh)
    t_deser_end = time.time()

    print(f"OUT TIME_DESERIALISATION: {t_deser_end - t_deser_start}")

    # Détermination du nombre de génomes (couleurs)
    if graph:
        nb_colors = 1 + max(max(c) for c in graph.values())
    else:
        nb_colors = 0

    queries = read_fasta_queries(fasta_queries)

    # Phase de requête
    t_query_start = time.time()
    with open(output_file, "w") as out:
        for header, seq in queries:
            kmers = list(get_kmers(seq, k))
            scores = [0] * nb_colors  # compteur par génome

            for kmer in kmers:
                if kmer in graph:
                    for cid in graph[kmer]:
                        scores[cid] += 1

            total = len(kmers)
            if total > 0:
                ratios = [s / total for s in scores]
            else:
                ratios = [0.0] * nb_colors

            out.write(header + "\t" + "\t".join(f"{r:.4f}" for r in ratios) + "\n")

    t_query_end = time.time()
    print(f"OUT TIME_QUERY: {t_query_end - t_query_start}")


if __name__ == "__main__":
    print("Ce fichier ne doit pas être exécuté directement. Utilisez dbg_indexer.py")
