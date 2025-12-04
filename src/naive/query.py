#!/usr/bin/env python3
import time
import pickle

def read_fasta_queries(path):
    """
    Lit un FASTA multi-entries et retourne une liste de tuples (header, sequence).
    Le header retourné est la ligne sans '>' (trim autour) — on préserve le texte après le nom.
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
    n = len(seq)
    if n < k:
        return
    for i in range(n - k + 1):
        yield seq[i:i+k]

def query_index(index_file, fasta_queries, k, output_file):
    """
    Charge l'index (pickle), mesure le temps de désérialisation, puis exécute la query.
    Écrit la sortie sous la forme demandée (header TAB ratios(4 décimales) TAB-sep).
    Imprime les logs conformes V2 :
      OUT TIME_DESERIALISATION: <float>
      OUT TIME_QUERY: <float>
    """
    # désérialisation (mesurée)
    t0 = time.time()
    with open(index_file, "rb") as fh:
        graph = pickle.load(fh)
    t_deser = time.time() - t0
    print(f"OUT TIME_DESERIALISATION: {t_deser}")

    # déterminer nb couleurs (s'il y a des k-mers)
    if graph:
        try:
            nb_colors = 1 + max(max(s) for s in graph.values() if s)
        except ValueError:
            nb_colors = 0
    else:
        nb_colors = 0

    queries = read_fasta_queries(fasta_queries)

    t_query_start = time.time()
    with open(output_file, "w") as out:
        for header, seq in queries:
            kmers = list(get_kmers(seq, k))
            if not kmers or nb_colors == 0:
                ratios = [0.0] * nb_colors
            else:
                scores = [0] * nb_colors
                for kmer in kmers:
                    if kmer in graph:
                        for cid in graph[kmer]:
                            scores[cid] += 1
                total = len(kmers)
                ratios = [s / total for s in scores]
            # écriture : header <TAB> r1 <TAB> r2 ...
            out.write(header + "\t" + "\t".join(f"{r:.4f}" for r in ratios) + "\n")
    t_query = time.time() - t_query_start
    print(f"OUT TIME_QUERY: {t_query}")

if __name__ == "__main__":
    print("Ce fichier ne doit pas être exécuté directement. Utilisez dbg_indexer.py")
