#!/usr/bin/env bash
set -euo pipefail

# Fichiers d'entrée
LIST="First_set/list_genomes.txt"
QUERY="First_set/query.fa"

# Fichier de sortie
OUT="bench_results.tsv"

# En-tête du tableau
echo -e "k\tversion\tbuild_s\tserial_s\tdeser_s\tquery_s\tindex_bytes" > "$OUT"

# Valeurs de k testées
for k in 15 21 27 31 41; do
  echo "=== k=$k ==="

  # ======================
  # NAIVE
  # ======================
  idx="index_naive_k${k}.bin"
  res="results_naive_k${k}.txt"

  b=$(python -m naive.dbg_indexer build -i "$LIST" -k "$k" -o "$idx")
  tb=$(echo "$b" | awk '/OUT TIME_BUILD:/ {print $3}')
  ts=$(echo "$b" | awk '/OUT TIME_SERIALISATION:/ {print $3}')

  q=$(python -m naive.dbg_indexer query -q "$QUERY" -i "$idx" -k "$k" -o "$res")
  td=$(echo "$q" | awk '/OUT TIME_DESERIALISATION:/ {print $3}')
  tq=$(echo "$q" | awk '/OUT TIME_QUERY:/ {print $3}')

  size=$(wc -c < "$idx" | tr -d ' ')
  echo -e "${k}\tnaive\t${tb}\t${ts}\t${td}\t${tq}\t${size}" >> "$OUT"

  # ======================
  # ADVANCED
  # ======================
  idx="index_advanced_k${k}.bin"
  res="results_advanced_k${k}.txt"

  b=$(python -m advanced.dbg_indexer build -i "$LIST" -k "$k" -o "$idx")
  tb=$(echo "$b" | awk '/OUT TIME_BUILD:/ {print $3}')
  ts=$(echo "$b" | awk '/OUT TIME_SERIALISATION:/ {print $3}')

  q=$(python -m advanced.dbg_indexer query -q "$QUERY" -i "$idx" -k "$k" -o "$res")
  td=$(echo "$q" | awk '/OUT TIME_DESERIALISATION:/ {print $3}')
  tq=$(echo "$q" | awk '/OUT TIME_QUERY:/ {print $3}')

  size=$(wc -c < "$idx" | tr -d ' ')
  echo -e "${k}\tadvanced\t${tb}\t${ts}\t${td}\t${tq}\t${size}" >> "$OUT"

done

echo
echo " Bench terminé. Résultats dans : $OUT"
