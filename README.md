
# Superkmers

## Principe

Ce projet consiste à regrouper des sous-séquences successives d'ADN ensemble en superkmers, pour permettre d'en rechercher des fragments plus rapidement. 

Pour cela, on parcourt une séquence ADN en entrée (fichier .fasta). On crée toutes les sous-séquences de taille k (qui se chevauchent), appelées k-mers. Puis on crée toutes les sous-séquences des k-mers de taille m, appelées m-mers. On appelle "minimiseur d'un k-mer" son m-mer le plus petit, dans l'ordre alphabétique (A>C>G>T). 

Enfin on regroupe les k-mers successifs qui ont le même minimiseur ensemble, formant ainsi les superkmers. 


## Utilisation

#### Compilation
```bash
gcc -Wall main.c -o main -lm
```

#### Exécution

Le programme prend en argument : 
- file.fasta : un chemin vers un fichier en format fasta, contenant la séquence ADN à analyser
- directory : un chemin vers le répertoire où sera créé le fichier de sortie avec la liste des superkmers
- k : la taille des k-mers, comprise entre 1 et 31
- m : la taille des m-mers, comprise entre 1 et k

```bash
./main file.fasta directory k m
```

Par exemple : 

```
./main data/ecoli.fasta out/ 31 13
```

## Contact

Adèle DESMAZIERES
adesmaz.pro@gmail.com