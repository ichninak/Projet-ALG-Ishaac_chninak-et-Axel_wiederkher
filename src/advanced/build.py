#!/usr/bin/env python3
import time
import pickle
from collections import defaultdict


def read_fasta(path):
    """
    Lit un fichier FASTA et renvoie la séquence nucléotidique concaténée.
    """
    seq = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


def get_kmers(seq, k):
    """
    Générateur de k-mers à partir d'une séquence.
    """
    n = len(seq)
    for i in range(n - k + 1):
        yield seq[i:i+k]


def build_dbg(genomes, k):
    """
    Construit le graphe de De Bruijn à partir d'une liste de génomes.

    - succ : dictionnaire des successeurs
    - pred : dictionnaire des prédécesseurs
    - colors : k-mer -> ensemble des couleurs (genomes)
    """
    succ = defaultdict(set)
    pred = defaultdict(set)
    colors = defaultdict(set)

    for color_id, seq in enumerate(genomes):
        for kmer in get_kmers(seq, k):
            prefix = kmer[:-1]
            suffix = kmer[1:]
            succ[prefix].add(suffix)
            pred[suffix].add(prefix)
            colors[kmer].add(color_id)

    return succ, pred, colors


def compact_dbg(succ, pred, k):
    """
    Compacte le graphe de De Bruijn en unitigs.

    Un unitig est un chemin maximal sans branchement
    (degré entrant = 1 et degré sortant = 1).

    Retourne :
    - la liste des unitigs (séquences)
    - un dictionnaire k-mer -> identifiant de unitig
    """
    visited = set()
    unitigs = []
    kmer_to_unitig = {}

    nodes = set(succ.keys()) | set(pred.keys())
    starts = []

    # Détection des noeuds de départ (non linéaires)
    for node in nodes:
        indeg = len(pred[node]) if node in pred else 0
        outdeg = len(succ[node]) if node in succ else 0
        if indeg != 1 or outdeg != 1:
            starts.append(node)

    # Construction des unitigs
    for start in starts:
        if start in visited:
            continue

        if len(succ[start]) == 0:
            unitigs.append(start)
            visited.add(start)
            continue

        for nxt in succ[start]:
            seq = start
            cur = nxt

            while len(pred[cur]) == 1 and len(succ[cur]) == 1:
                seq += cur[-1]
                cur = next(iter(succ[cur]))

            seq += cur[-1]

            uid = len(unitigs)
            unitigs.append(seq)

            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                kmer_to_unitig[kmer] = uid

            visited.add(cur)

    return unitigs, kmer_to_unitig


def color_unitigs(unitigs, kmer_to_unitig, colors):
    """
    Attribue à chaque unitig l'ensemble des couleurs
    des k-mers qui le composent.
    """
    unitig_colors = [set() for _ in unitigs]

    for kmer, cset in colors.items():
        if kmer in kmer_to_unitig:
            uid = kmer_to_unitig[kmer]
            unitig_colors[uid].update(cset)

    return unitig_colors


def build_index(list_file, k, output_file):
    """
    Construit un graphe de De Bruijn coloré compacté (version avancée).

    L'utilisation des unitigs permet de factoriser l'information
    de couleur et de réduire la taille de la structure.
    """
    with open(list_file, "r") as f:
        genome_paths = [line.strip() for line in f if line.strip()]

    genomes = [read_fasta(p) for p in genome_paths]

    t_build_start = time.time()
    succ, pred, colors = build_dbg(genomes, k)
    unitigs, kmer_to_unitig = compact_dbg(succ, pred, k)
    unitig_colors = color_unitigs(unitigs, kmer_to_unitig, colors)
    t_build_end = time.time()

    print(f"OUT TIME_BUILD: {t_build_end - t_build_start}")

    t_ser_start = time.time()
    data = {
        "k": k,
        "unitigs": unitigs,
        "unitig_colors": unitig_colors,
        "kmer_to_unitig": kmer_to_unitig
    }

    with open(output_file, "wb") as out:
        pickle.dump(data, out, protocol=pickle.HIGHEST_PROTOCOL)
    t_ser_end = time.time()

    print(f"OUT TIME_SERIALISATION: {t_ser_end - t_ser_start}")


if __name__ == "__main__":
    print("Ce fichier ne doit pas être exécuté directement. Utilisez dbg_indexer.py")
