#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <math.h>
#include <string.h>

#include "kmer.h"

// input : un charactère nucléotide A, T, C ou G
// output : son code correspondant, A:0 C:1 T:2 G:3
int nucleotide_to_int(char n) {
	switch (n) {
	case 'A':
		return 0;
	case 'C' :
		return 1;
	case 'T' :
		return 2;
	case 'G' :
		return 3;
	default:
		perror("Nucléotide inconnu.\n");
		return -1;
	}
}

// input : un string représentant une séquence ADN
// output : l'entier codant cette séquence, en base 4
int seq_to_int(char* seq) {
	int result = 0;
	int len = strlen(seq);
	for (int i = 0; i < len; i++) {
		result += nucleotide_to_int(seq[i]) * pow(4, len-1-i);
	}
	return result;
}

// input : un string de la séquence ADN d'un kmer
// output : le kmer avec sa séquence convertie en int
// et ses autres attribus initialisés à une valeur nulle
kmer *seq_to_kmer(char* seq) {
	kmer *ret = malloc(sizeof(kmer));
	if (ret == NULL) { perror("malloc"); }
	ret->sequence = seq_to_int(seq);
	ret->mmers = NULL;
	ret->minimiseur = -1;
	return NULL;
}


// input : les entiers k et m
// output : -1 si k ou m ne respecte pas les conditions d'input, et 0 sinon
// il faut k <= 31 et m <= k
int check_args_k_m(int k, int m) {
	if (k > 31) {
		printf("k doit être strictement inférieur à 31.\n");
		return -1;
	}
	if (m > k) {
		printf("m doit être inférieur à k.\n");
		return -1;
	}
	return 0;
}

// fonction principale
int main(int argc, char *argv[]) {
	
	// récupération des arguments
	if (argc != 5) {
		printf("4 arguments requis.\n");
		return -1;
	}
	
	char* file_input_name = argv[1];
	char* dir_output_name = argv[2];
	int k = atoi(argv[3]);
	int m = atoi(argv[4]);
	
	if (check_args_k_m(k, m) != 0) { return -1; }
	
	FILE *finput = fopen(file_input_name, "r");
	if (finput == NULL) { printf("Fichier d'entrée introuvable.\n"); return -1; }
	
	DIR* dir_out = opendir(dir_output_name);
	if (! dir_out) { printf("Dossier de sortie introuvable.\n"); return -1; }
	
	
	// traitement des données
	printf("%d\n", seq_to_int("CTT"));
	
	
	
	
	
    closedir(dir_out);
	
	return 0;
}
