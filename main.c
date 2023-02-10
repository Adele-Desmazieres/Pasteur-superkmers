#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <math.h>
#include <string.h>

#include "kmer.h"

#define BASE 4

// input : un charactère de nucléotide A, T, C ou G
// output : son code correspondant, A:0 C:1 T:2 G:3
int nucleotide_to_val(char n) {
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

// input : un string représentant une séquence ADN et la longueur len
// output : l'entier codant cette séquence jusqu'à len, en base 4
int seq_to_val(char* seq, int len) {
	int result = 0;
	for (int i = 0; i < len; i++) {
		result += nucleotide_to_val(seq[i]) * pow(BASE, len-1-i);
	}
	return result;
}

// comparaison de seq1 et seq2 du début au charactère len
// output : -1 si seq1 précède seq2 selon l'ordre alphabetique
//          0 si elles sont égales
//          1 sinon
int cmp_seq(char *seq1, char *seq2, int len) {
	return strncmp(seq1, seq2, len);
}

// fonction appelée une fois, à l'initialisation
// input : un string de la séquence ADN d'un kmer 
//         et k la taille du kmer
//         et m la taille de ses mmers
// output : le kmer avec sa séquence convertie en int
// et ses autres attribus bien initialisés
kmer *seq_to_kmer(char *seq, int k, int m) {
	kmer *ret = malloc(sizeof(kmer));
	if (ret == NULL) { perror("malloc"); }
	ret->seq_val = seq_to_val(seq, k);
	
	// tableau de (k+1-m) strings, qui sont les mmers du kmer
	int nbr_mmers = k+1-m;
	int *mmers = malloc(sizeof(int) * nbr_mmers);
	char *seq_minimiseur = seq; // met a jour le minimieur
	
	for (int i = 0; i < k+1-m; i++) {
		char *seq_current = seq+i;
		int val = seq_to_val(seq_current, m);
		mmers[i] = val;
		
		if (cmp_seq(seq_current, seq_minimiseur, m) < 0) {
			seq_minimiseur = seq_current;
		}
	}
	ret->mmers = mmers;
	ret->minimiseur = seq_to_val(seq_minimiseur, m);
	ret->nbr_mmers = nbr_mmers;
	ret->k = k;
	ret->m = m;
	
	return ret;
}

// output : l'entier correspondant à la séquence de taille x, suivant celle de xmer
//          avec l'ajout du nucléotide new_nucl
int next_xmer_val(int xmer, int x, char new_nucl) {
	// suppression du coef de poids fort = suppression du nucléotide le plus à gauche
	int tmp = pow(BASE, x-1);
	while (xmer >= tmp) {
		xmer -= tmp; 
	}
	// décalage de bit
	xmer *= BASE; // TODO ? a remplacer avec un vrai bitshift <<
	// ajout nouveau nucléotide
	xmer += nucleotide_to_val(new_nucl);
	return xmer;
}


kmer *next_kmer(kmer *previous, char new_nucl) {
	kmer *ret = malloc(sizeof(kmer));
	if (ret == NULL) { perror("malloc"); }
	
	// TODO
	
	return NULL;
}

// output : le string de la séquence ADN de taille len, correspondant à l'entier en argument
char *val_to_seq(int val, int len) {
	char *ret = malloc(sizeof(char) * len);
	for (int i = 0; i < len; i++) {
		int nucl_val = 0;
		// convertisseur base 10 en base 4
		while (val >= pow(BASE, len-1-i)) {
			val -= pow(BASE, len-1-i);
			nucl_val += 1;
		}
		char c = '-';
		if (nucl_val == 0) { c = 'A'; }
		else if (nucl_val == 1) { c = 'C'; }
		else if (nucl_val == 2) { c = 'T'; }
		else if (nucl_val == 3) { c = 'G'; }
		
		ret[i] = c;
	}
	return ret;
}

// affiche le contenu d'un kmer, ses mmers et son minimiseur. Ne renvoie rien
void print_kmer(kmer *km) {
	printf("%s\n", val_to_seq(km->seq_val, km->k));
	printf("   mmers : ");
	for (int i = 0; i < km->nbr_mmers; i++) {
		printf("%s, ", val_to_seq(km->mmers[i], km->m));
	}
	printf("\n   minimiseur : %s\n", val_to_seq(km->minimiseur, km->m));
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
	/*
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
	*/
	
	// traitement des données
	int xmer1 = seq_to_val("CTT", 3); // 01 10 10 = 26
	printf("%d\n", xmer1);
	
	int xmer2 = next_xmer_val(xmer1, 3, 'A'); // 10 10 00 = 40
	printf("%d\n", xmer2);
	
	kmer *km1 = seq_to_kmer("TTAGGACAAA", 5, 3);
	print_kmer(km1);
	
	
	
    //closedir(dir_out);
	
	return 0;
}
