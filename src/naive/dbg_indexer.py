#!/usr/bin/env python3
import argparse
from build import build_index
from query import query_index

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd")

    # ---- BUILD ----
    b = sub.add_parser("build")
    b.add_argument("-i", dest="input_list", required=True, help="Fichier liste de génomes")
    b.add_argument("-k", type=int, required=True, help="Taille des k-mers")
    b.add_argument("-o", dest="output_index", required=True, help="Fichier de sortie (index)")

    # ---- QUERY ----
    q = sub.add_parser("query")
    q.add_argument("-q", dest="query_fasta", required=True, help="Fichier FASTA des requêtes")
    q.add_argument("-i", dest="input_index", required=True, help="Index sérialisé")
    q.add_argument("-o", dest="output_results", required=True, help="Fichier de sortie (similarités)")

    args = parser.parse_args()

    if args.cmd == "build":
        print(">>> Building index...")
        build_index(args.input_list, args.k, args.output_index)

    elif args.cmd == "query":
        print(">>> Querying index...")
        query_index(args.input_index, args.query_fasta, args.k, args.output_results)

    else:
        parser.print_help()

if __name__ == "__main__":
    main()
