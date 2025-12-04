#!/usr/bin/env python3
import argparse
from .build import build_index
from .query import query_index


def main():
    parser = argparse.ArgumentParser(description="Advanced cDBG builder & querier")
    sub = parser.add_subparsers(dest="cmd")

    # BUILD
    b = sub.add_parser("build", help="Build advanced cDBG")
    b.add_argument("-i", required=True, help="list of genome paths")
    b.add_argument("-k", type=int, required=True, help="k-mer size")
    b.add_argument("-o", required=True, help="output file")

    # QUERY
    q = sub.add_parser("query", help="Query advanced cDBG")
    q.add_argument("-q", required=True)
    q.add_argument("-i", required=True)
    q.add_argument("-k", type=int, required=True)
    q.add_argument("-o", required=True)

    args = parser.parse_args()

    if args.cmd == "build":
        build_index(args.i, args.k, args.o)
    elif args.cmd == "query":
        query_index(args.i, args.q, args.k, args.o)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
