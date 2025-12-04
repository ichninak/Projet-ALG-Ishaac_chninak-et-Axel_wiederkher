#!/usr/bin/env python3
import time
import pickle
from collections import defaultdict


def read_fasta(path):
    seq = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


def get_kmers(seq, k):
    n = len(seq)
    for i in range(n - k + 1):
        yield seq[i:i+k]


def build_dbg(genomes, k):
    """
    Construit un DBG simple : 
        succ[prefix] = set of suffixes
        pred[suffix] = set of prefixes
        colors[kmer] = set(colors)
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
    Construit les unitigs à partir du DBG.
    Retourne:
        - unitigs (list of strings)
        - kmer_to_unitig (dict kmer -> unitig ID)
    """

    visited = set()
    unitigs = []
    kmer_to_unitig = {}

    # Trouver les noeuds start = indegree != 1 ou outdegree != 1
    nodes = set(succ.keys()) | set(pred.keys())
    starts = []

    for node in nodes:
        indeg = len(pred[node]) if node in pred else 0
        outdeg = len(succ[node]) if node in succ else 0
        if indeg != 1 or outdeg != 1:
            starts.append(node)

    # Construire les unitigs depuis les start nodes
    for start in starts:
        if start in visited:
            continue

        # Si pas de succ, c'est un unitig trivial
        if len(succ[start]) == 0:
            unitigs.append(start)
            visited.add(start)
            continue

        # Pour chaque successeur
        for nxt in succ[start]:

            seq = start
            cur = nxt

            # Avancer tant que noeud linéaire
            while len(pred[cur]) == 1 and len(succ[cur]) == 1:
                seq += cur[-1]   # ajouter dernière base
                next_node = next(iter(succ[cur]))
                cur = next_node

            # Ajouter dernier noeud
            seq += cur[-1]

            # Enregistrer le unitig
            uid = len(unitigs)
            unitigs.append(seq)

            # Marquer tous les k-mers du unitig
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                kmer_to_unitig[kmer] = uid

            visited.add(cur)

    return unitigs, kmer_to_unitig


def color_unitigs(unitigs, kmer_to_unitig, colors):
    """
    Attribue les couleurs aux unitigs.
    """
    unitig_colors = [set() for _ in unitigs]

    for kmer, cset in colors.items():
        if kmer in kmer_to_unitig:
            uid = kmer_to_unitig[kmer]
            unitig_colors[uid].update(cset)

    return unitig_colors


def build_index(list_file, k, output_file):
    t0 = time.time()

    # Charger génomes
    with open(list_file, "r") as f:
        genome_paths = [line.strip() for line in f if line.strip()]
    genomes = [read_fasta(p) for p in genome_paths]

    # BUILD DBG
    t_build0 = time.time()
    succ, pred, colors = build_dbg(genomes, k)
    unitigs, kmer_to_unitig = compact_dbg(succ, pred, k)
    unitig_colors = color_unitigs(unitigs, kmer_to_unitig, colors)
    t_build1 = time.time()

    print(f"OUT TIME_BUILD: {t_build1 - t_build0}")

    # SERIALISATION
    t_ser0 = time.time()
    data = {
        "k": k,
        "unitigs": unitigs,
        "unitig_colors": unitig_colors,
        "kmer_to_unitig": kmer_to_unitig
    }

    with open(output_file, "wb") as out:
        pickle.dump(data, out, protocol=pickle.HIGHEST_PROTOCOL)

    t_ser1 = time.time()
    print(f"OUT TIME_SERIALISATION: {t_ser1 - t_ser0}")


if __name__ == "__main__":
    print("Use dbg_indexer.py")
