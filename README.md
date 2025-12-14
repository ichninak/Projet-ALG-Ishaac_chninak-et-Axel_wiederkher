# cDBG Indexer — Graphe de de Bruijn coloré

Ce projet implémente un outil d’indexation et de requêtage de génomes basé sur les **graphes de de Bruijn colorés (cDBG)**.
Il permet de comparer une ou plusieurs séquences requêtes à une collection de génomes par une approche *alignment-free*, en calculant des **ratios de k-mers partagés**.

Deux versions sont proposées :

* une version **naïve**, simple mais potentiellement coûteuse en mémoire ;
* une version **avancée**, exploitant la compaction du graphe en **unitigs** afin de réduire la redondance des informations de couleur.

---

## Organisation du dépôt

```
.
├── README.md
├── report.pdf
└── src
    ├── naive
    │   ├── build.py
    │   ├── query.py
    │   └── dbg_indexer.py
    └── advanced
        ├── build.py
        ├── query.py
        └── dbg_indexer.py
```

---

## Prérequis

* Python **3.8 ou plus**
* Bibliothèques standard uniquement (aucune dépendance externe)

---

## Données d’entrée

### Génomes

Un fichier texte contenant **un chemin vers un fichier FASTA par ligne**.

Exemple :

```
data/genome1.fa
data/genome2.fa
data/genome3.fa
```

### Requêtes

Un fichier FASTA contenant une ou plusieurs séquences à requêter.

---

## Construction de l’index

### Version naïve

```
python -m src.naive.dbg_indexer build \
  -i genomes_list.txt \
  -k 31 \
  -o naive_index.pkl
```

### Version avancée

```
python -m src.advanced.dbg_indexer build \
  -i genomes_list.txt \
  -k 31 \
  -o advanced_index.pkl
```

Lors de la construction, les temps suivants sont affichés dans la console :

* `OUT TIME_BUILD` : temps de construction de l’index
* `OUT TIME_SERIALISATION` : temps de sauvegarde sur disque

---

## Requêtage

### Version naïve

```
python -m src.naive.dbg_indexer query \
  -q queries.fa \
  -i naive_index.pkl \
  -k 31 \
  -o results_naive.txt
```

### Version avancée

```
python -m src.advanced.dbg_indexer query \
  -q queries.fa \
  -i advanced_index.pkl \
  -k 31 \
  -o results_advanced.txt
```

Lors du requêtage, les temps suivants sont affichés :

* `OUT TIME_DESERIALISATION` : temps de lecture de l’index
* `OUT TIME_QUERY` : temps de calcul des similarités

---

## Format de sortie

Le fichier de sortie contient **une ligne par séquence requête** :

```
<header>    sim_G1    sim_G2    ...    sim_GN
```

où `sim_Gi` est le **ratio de k-mers partagés** entre la requête et le génome `Gi`, arrondi à **4 décimales**.

Exemple :

```
query_1    0.1967    0.1964    0.1861
query_2    1.0000    0.0000    0.0000
```

---

## Remarques

* La valeur de `k` doit être **identique** lors de la construction et du requêtage.
* Les programmes ne sont pas interactifs : une fois lancés avec leurs arguments, ils s’exécutent et se terminent automatiquement.
* La version avancée exploite la factorisation de l’information de couleur au niveau des unitigs afin de réduire la taille de l’index.

---

## Auteurs

* **Ishaac Chninak**
* **Axel Wiederkher**

