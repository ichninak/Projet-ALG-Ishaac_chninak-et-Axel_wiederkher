#!/usr/bin/env python3
import time
import pickle


def read_fasta(path):
    """
    Lit un fichier FASTA (une ou plusieurs lignes par séquence)
    et renvoie la séquence concaténée sous forme de chaîne.
    """
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq_parts.append(line.strip().upper())
    return "".join(seq_parts)


def get_kmers(sequence, k):
    """
    Générateur de k-mers à partir d'une séquence,
    en utilisant une fenêtre glissante de taille k.
    """
    n = len(sequence)
    for i in range(n - k + 1):
        yield sequence[i:i+k]


def build_index(list_file_path, k, output_file):
    """
    Construit un graphe de De Bruijn coloré naïf.

    Chaque k-mer est stocké comme clé d'un dictionnaire,
    et la valeur associée est l'ensemble des identifiants
    des génomes (couleurs) dans lesquels ce k-mer apparaît.

    Le graphe est ensuite sérialisé sur disque à l'aide de pickle.

    Affiche les temps d'exécution au format imposé :
    - OUT TIME_BUILD
    - OUT TIME_SERIALISATION
    """
    # Lecture de la liste des fichiers génomes
    with open(list_file_path, "r") as fh:
        genome_paths = [ln.strip() for ln in fh if ln.strip()]

    graph = {}  # dictionnaire : k-mer -> ensemble des couleurs (genome IDs)

    # Construction du graphe
    t_build_start = time.time()
    for color_id, path in enumerate(genome_paths):
        seq = read_fasta(path)
        for kmer in get_kmers(seq, k):
            if kmer in graph:
                graph[kmer].add(color_id)
            else:
                graph[kmer] = {color_id}
    t_build_end = time.time()

    print(f"OUT TIME_BUILD: {t_build_end - t_build_start}")

    # Sérialisation du graphe
    t_ser_start = time.time()
    with open(output_file, "wb") as out:
        pickle.dump(graph, out, protocol=pickle.HIGHEST_PROTOCOL)
    t_ser_end = time.time()

    print(f"OUT TIME_SERIALISATION: {t_ser_end - t_ser_start}")


if __name__ == "__main__":
    print("Ce fichier ne doit pas être exécuté directement. Utilisez dbg_indexer.py")
