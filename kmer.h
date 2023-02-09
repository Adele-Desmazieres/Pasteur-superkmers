#ifndef KMER_H
#define KMER_H

struct kmer {
    int sequence;
    int* mmers;
    int minimiseur;
};
typedef struct kmer kmer;


#endif