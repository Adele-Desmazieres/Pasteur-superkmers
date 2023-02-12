#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <math.h>
#include <string.h>
#include <limits.h>

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
		printf("Erreur : nucléotide inconnu.\n");
		exit(1);
	}
}

// input : un string représentant une séquence ADN et la longueur len
// output : l'entier codant cette séquence jusqu'à len, en base 4
long seq_to_val(char* seq, int len) {
	long result = 0;
	for (int i = 0; i < len; i++) {
		long tmp = ((long)nucleotide_to_val(seq[i])) * ((long)pow(BASE, len-1-i));
		result += tmp;
	}
	return result;
}

// output : le string de la séquence ADN de taille len, correspondant à l'entier en argument
char *val_to_seq(long val, int len) {
	char *ret = malloc(sizeof(char) * len + 1);
	for (int i = 0; i < len; i++) {
		int nucl_val = 0;
		// convertisseur base 10 en base 4
		while (val >= ((long)pow(BASE, len-1-i))) {
			val -= ((long)pow(BASE, len-1-i));
			nucl_val += 1;
		}
		char c = '-';
		if (nucl_val == 0) { c = 'A'; }
		else if (nucl_val == 1) { c = 'C'; }
		else if (nucl_val == 2) { c = 'T'; }
		else if (nucl_val == 3) { c = 'G'; }
		
		ret[i] = c;
	}
	ret[len] = '\0';
	return ret;
}

// affiche le contenu d'un kmer, ses mmers et son minimiseur. Ne renvoie rien
void print_kmer(kmer *km) {
	char *tmp = val_to_seq(km->seq_val, km->k);
	printf("%ld %s (", km->seq_val, tmp);
	free(tmp);
	
	tmp = val_to_seq(km->minimiseur, km->m);
	printf("%s) [", tmp);
	free(tmp);
	
	printf("%d", km->m);
	
	//for (int i = 0; i < km->nbr_mmers; i++) {
		//printf("%s ", val_to_seq(km->mmers[i], km->m));
	//}
	printf("]\n");
}

// comparaison de seq1 et seq2 du début au charactère len
// output : -1 si seq1 précède seq2 selon l'ordre alphabetique
//          0 si elles sont égales
//          1 sinon
int cmp_seq(char *seq1, char *seq2, int len) {
	return strncmp(seq1, seq2, len);
}

int cmp_val(int val1, int val2, int len) {
	char *tmp1 = val_to_seq(val1, len);
	char *tmp2 = val_to_seq(val2, len);
	int ret = cmp_seq(tmp1, tmp2, len);
	free(tmp1);
	free(tmp2);
	return ret;
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
	
	int nbr_mmers = k+1-m;
	long *mmers = malloc(sizeof(long) * nbr_mmers); // tableau des mmers du kmer en nombres
	char *seq_minimiseur = seq; // garde à jour le minimieur
	
	for (int i = 0; i < k+1-m; i++) {
		char *seq_current = seq+i;
		long val = seq_to_val(seq_current, m);
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
long next_xmer_val(long xmer, int x, char new_nucl) {
	// suppression du coef de poids fort = suppression du nucléotide le plus à gauche
	long tmp = (long)pow(BASE, x-1);
	while (xmer >= tmp) {
		xmer -= tmp; 
	}
	
	// décalage de bit
	xmer *= BASE;
	
	// ajout nouveau nucléotide
	xmer += nucleotide_to_val(new_nucl);
	
	return xmer;
}

// output : le kmer qui suit celui en argument, après l'ajout du nucléotide new_nucl
kmer *next_kmer(kmer *previous, char new_nucl) {
	kmer *ret = malloc(sizeof(kmer));
	if (ret == NULL) { perror("malloc"); }
	
	int k = previous->k;
	int m = previous->m;
	
	ret->k = k;
	ret->m = m;
	ret->nbr_mmers = previous->nbr_mmers;
	ret->seq_val = next_xmer_val(previous->seq_val, k, new_nucl); // nouveau kmer
	
	// nouveau mmer
	long new_mmer = next_xmer_val(previous->mmers[ret->nbr_mmers-1], m, new_nucl);
	
	// nouveau tableau des mmers à partir du précédent, en ajoutant le nouveau mmer et supprimant le 1er
	long *new_mmers = malloc(ret->nbr_mmers * sizeof(long));
	if (new_mmers == NULL) { perror("malloc"); }
	memcpy(new_mmers, previous->mmers + 1, ret->nbr_mmers * sizeof(long));
	new_mmers[ret->nbr_mmers-1] = new_mmer;
	ret->mmers = new_mmers;
	
	// nouveau minimiseur
	// check si le nouveau mmer et plus petit que le dernier minimiseur
	if (cmp_val(new_mmer, previous->minimiseur, m) < 0) {
		ret->minimiseur = new_mmer;
	}	
	// check si l'ancien minimiseur a été supprimé (vérifier s'il etait en premiere place)
	// dans ce cas, trouver le nouveau minimiseur en parcourant toute la liste
	else if (previous->minimiseur == previous->mmers[0]) {
		
		ret->minimiseur = ret->mmers[0];
		for (int i = 1; i < ret->nbr_mmers; i++) { // parcours de la liste des mmers
			
			if (cmp_val(ret->mmers[i], ret->minimiseur, m) < 0) {
				ret->minimiseur = ret->mmers[i];
			}
		}
	}
	// sinon on garde l'ancien minimiseur
	else { 
		ret->minimiseur = previous->minimiseur;
	}
	
	return ret;
}

// libère la mémoire d'un kmer
void free_kmer(kmer *kmer) {
	free(kmer);
}


int write_file_c(FILE *f, char *buff, int buffsize, int remplissage, char c) {
	buff[remplissage] = c;
	remplissage += 1;
	if (remplissage >= buffsize) {
		int written = fwrite(buff, sizeof(char), remplissage, f);
		if (written != remplissage) { perror("fwrite"); }
		remplissage = 0;
	}
	return remplissage;
}

int write_file_s(FILE *f, char *buff, int buffsize, int remplissage, char *s) {
	if (strlen(s) + remplissage < buffsize) {
		strncpy(buff+remplissage, s, strlen(s));
		remplissage += strlen(s);

	} else {
		int written = fwrite(buff, sizeof(char), remplissage, f);
		if (written != remplissage) { perror("fwrite"); }
		written = fwrite(s, sizeof(char), strlen(s)-1, f);
		if (written != strlen(s)-1) { perror("fwrite"); }
		remplissage = 0;
	}
	return remplissage;
}

void read_file(FILE *filein, FILE *fileout, int k, int m) {
	// ignore la première ligne
	int buffsize = 1024;
	char buff[buffsize];
	buff[buffsize] = '\0';
	if (fgets(buff, buffsize, filein) == NULL) {
		printf("Erreur : première ligne vide.\n");
		exit(1);
	}
	int remplissage = 0;
	
	// lit les k premiers caractères et crée le 1er kmer
	char first_seq[k];
	first_seq[k] = '\0';
	int nb_read = fread(first_seq, sizeof(char), k, filein);
	if (nb_read < k) {
		printf("Erreur : séquence trop courte, de taille inférieure à k.\n");
		exit(1);
	}
	kmer *current = seq_to_kmer(first_seq, k, m);
	kmer *next;
	
	char nucl;
	remplissage = write_file_s(fileout, buff, buffsize, remplissage, first_seq);
	
	// lit tous les autres caractères de la séquence et crée les kmers
	while ((nucl = fgetc(filein)) != EOF) {
		// ignore les nucléotides inconnus et les sauts de ligne
		if (nucl == 'N' || nucl == '\n') {
			continue;
		}
		
		next = next_kmer(current, nucl);

		// si le nouveau minimiseur est identique au précédent
		// alors on écrit le caractère à la suite du superkmer dans le fichier
		if (next->minimiseur == current->minimiseur) {
			remplissage = write_file_c(fileout, buff, buffsize, remplissage, nucl);
		
		// sinon on atteint la fin d'un superkmer, et on passe au suivant
		// en écrivant la séquence du kmer actuel
		} else {
			remplissage = write_file_c(fileout, buff, buffsize, remplissage, '\n');
			char *tmp = val_to_seq(next->seq_val, k);
			remplissage = write_file_s(fileout, buff, buffsize, remplissage, tmp);
			free(tmp);
		}
		
		free_kmer(current);
		current = next;
	}
	free_kmer(next);
	int written = fwrite(buff, sizeof(char), remplissage, fileout);
	if (written != remplissage) { perror("fwrite"); }
	remplissage = 0;
}


// input : les entiers k et m
// output : -1 si k ou m ne respecte pas les conditions d'input, et 0 sinon
// il faut k <= 31 et m <= k
int check_args_k_m(int k, int m) {
	if (k > 31) {
		printf("Erreur : k doit être strictement inférieur à 31.\n");
		return -1;
	}
	if (m > k) {
		printf("Erreur : m doit être inférieur à k.\n");
		return -1;
	}
	return 0;
}


// fonction principale
int main(int argc, char *argv[]) {
	
	printf("Running...\n");
	
	// récupération des arguments
	if (argc != 5) {
		printf("4 arguments requis.\n");
		exit(1);
	}
	
	char* file_input_name = argv[1];
	char* dir_output_name = argv[2];
	int k = atoi(argv[3]);
	int m = atoi(argv[4]);
	
	if (check_args_k_m(k, m) != 0) { exit(1); }
	
	// fichier d'entrée
	FILE *filein = fopen(file_input_name, "r");
	if (! filein) { 
		printf("Erreur : fichier d'entrée introuvable.\n"); 
		exit(1); 
	}
	
	// fichier de sortie
	int dlen = strlen(dir_output_name);
	char *nameout;
	if (dir_output_name[dlen-1] == '/') {
		nameout = malloc(sizeof(char) * (dlen + 8));
		if (! nameout) { perror("malloc"); }
		sprintf(nameout, "%sout.txt", dir_output_name);		
	} else {
		nameout = malloc(sizeof(char) * (dlen + 1 + 8));
		if (! nameout) { perror("malloc"); }
		sprintf(nameout, "%s/out.txt", dir_output_name);		
	}
		
	FILE *fileout = fopen(nameout, "w+");
	if (! fileout) {
		fclose(filein);
		printf("Erreur : impossible de créer %s.\n", nameout); 
		exit(1); 		
	}
	
	// traitement des données	
	read_file(filein, fileout, k, m);

	// fermeture des fichiers
	fclose(fileout);
	fclose(filein);
	
	printf("Done.\n");
	
	return 0;
}
