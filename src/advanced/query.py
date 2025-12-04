#!/usr/bin/env python3
import time
import pickle


def read_fasta_queries(path):
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
    n = len(seq)
    for i in range(n - k + 1):
        yield seq[i:i+k]


def query_index(index_file, query_file, k, output_file):
    # DESERIALISATION
    t0 = time.time()
    with open(index_file, "rb") as f:
        data = pickle.load(f)
    t1 = time.time()

    print(f"OUT TIME_DESERIALISATION: {t1 - t0}")

    unitigs = data["unitigs"]
    unitig_colors = data["unitig_colors"]
    kmer_to_unitig = data["kmer_to_unitig"]

    nb_colors = max((max(c) for c in unitig_colors if c), default=-1) + 1

    # QUERY
    t_q0 = time.time()
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

            # FORMAT CORRECT
            out.write(header + "\t" + "\t".join(f"{r:.4f}" for r in ratios) + "\n")

    t_q1 = time.time()
    print(f"OUT TIME_QUERY: {t_q1 - t_q0}")


if __name__ == "__main__":
    print("Use dbg_indexer.py")
