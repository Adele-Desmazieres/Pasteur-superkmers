#ifndef KMER_H
#define KMER_H

struct kmer {
    int k; // taille du kmers en nbr de nucléotides
    long seq_val; // entier représentant la séquence ADN du kmer en base 4
    int *mmers; // tableau des entiers représentant la séquence ADN de chaque mmer
    // ordonné dans le même ordre que les mmers de la séquence
    int m; // taille des mmers
    int nbr_mmers; // le nombre de mmers
    int minimiseur; // entier représentant la séquence ADN du miniseur, présent aussi dans mmers
};
typedef struct kmer kmer;


#endif