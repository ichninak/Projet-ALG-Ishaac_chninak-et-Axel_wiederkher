#!/usr/bin/env python3
"""
Interface en ligne de commande pour la version avancée
du graphe de De Bruijn coloré basé sur les unitigs.
"""
import argparse
from .build import build_index
from .query import query_index


def main():
    parser = argparse.ArgumentParser(
        description="Advanced colored De Bruijn graph indexer"
    )
    sub = parser.add_subparsers(dest="cmd")

    # BUILD
    b = sub.add_parser("build", help="Construire l'index avancé")
    b.add_argument("-i", required=True, help="Fichier listant les génomes")
    b.add_argument("-k", type=int, required=True, help="Taille des k-mers")
    b.add_argument("-o", required=True, help="Fichier de sortie de l'index")

    # QUERY
    q = sub.add_parser("query", help="Interroger l'index avancé")
    q.add_argument("-q", required=True, help="Fichier FASTA des requêtes")
    q.add_argument("-i", required=True, help="Index sérialisé")
    q.add_argument("-k", type=int, required=True, help="Taille des k-mers")
    q.add_argument("-o", required=True, help="Fichier de sortie")

    args = parser.parse_args()

    if args.cmd == "build":
        build_index(args.i, args.k, args.o)
    elif args.cmd == "query":
        query_index(args.i, args.q, args.k, args.o)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
