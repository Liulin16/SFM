#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <omp.h>

#include "LCS.h"

int main() {
	lcs_pubkey_t *pubkey;
	lcs_prvkey_t *prvkey;
	lcs_pub_para *PP;
	lcs_mkey *MK;

	mpz_t ztmp;
	mpz_init(ztmp);

	printf("\nGenerating 1024 bit key\n");
	// Generate 1024 bit keys
	lcs_PP_gen(1024, &PP, &MK, &lcs_get_rand_devurandom);
	//Gnenerate public and private key pairs
	lcs_keygen(PP, &pubkey, &prvkey);

	printf("PP N bit length is %ld\n", mpz_sizeinbase(PP->N , 2));

	/***********************
	Test for some functions
	***********************/

	//Test for load_dataset_S
	/*
	Dataset Data_A, Data_B;

	if(!load_dataset_S(&Data_A, &Data_B, PP)) {
		printf("Oops! I cannot read this file.\n" );
		return 0;
	}
	printf("%d\t%d\n", Data_A.M, NN);
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < NN; j++) {
			mpz_add(ztmp, Data_A.itemset[i].item[j]->m, Data_B.itemset[i].item[j]->m);
			if(mpz_cmp_ui(ztmp, 1) == 0) printf("%d\n", j+1);
		}
	}*/


	//Test for SF1I
	
	/*Dataset L1 = SF1I(&Data_A, &Data_B, pubkey, PP, MK);
	
	printf("%d\n", L1.M);
	for(int i = 0; i < L1.M; i++) {
		for(int j = 0; j < NN; j++) {
			if(mpz_cmp_ui(L1.itemset[i].item[j]->m, 1) == 0) printf("%d\n", j+1);
		}
	}*/

	
	//Test for SSC
	/*int num = 15; //the num for counting
	Itemset t;
	for(int i = 0; i < NN; i++) {
		t.item[i] = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
		pt_init(t.item[i]);
		mpz_set_ui(t.item[i]->m, 0);
	}
	mpz_set_ui(t.item[num-1]->m, 1);

	lcs_ciphertext_t** l1 = (lcs_ciphertext_t**)malloc(NN*sizeof(lcs_ciphertext_t*));
	Itemset G_A, G_B;
	for(int i = 0; i < NN; i++)  {
		l1[i] = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
		ct_init(l1[i]);
		lcs_enc(l1[i], PP, pubkey, t.item[i], &lcs_get_rand_devurandom);

		G_A.item[i] = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
		pt_init(G_A.item[i]);
		G_B.item[i] = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
		pt_init(G_B.item[i]);
		B2S(G_A.item[i], G_B.item[i], l1[i], pubkey, PP, MK);
	}
	G_A.lambda = 1;
	G_B.lambda = 1;

	lcs_ciphertext_t* P = SSC(&Data_A, &Data_B, &G_A, &G_B, pubkey, PP, MK);
	lcs_plaintext_t* res = lcs_dec(NULL, PP, pubkey, prvkey, P);
	gmp_printf("```````%Zd\n", res->m);*/

	printf("\n***ININ SOME ARGUEMENTS\n");
	clock_t start = clock(), diff;
	float wtime =  (float)omp_get_wtime(), diff_wtime;

	lcs_ciphertext_t*** Data = (lcs_ciphertext_t***)malloc(MM * sizeof(lcs_ciphertext_t**));
	for(int i = 0; i < MM; i++) {
		Data[i] = (lcs_ciphertext_t**)malloc(NN * sizeof(lcs_ciphertext_t*));
	}

	diff = clock() - start;
	int msec = diff * 1000 / CLOCKS_PER_SEC;
	printf("\tSUCCESS in %d miliseconds\n", msec);

	printf("\n***LOADING DATASETS\n");
	start = clock();

	if(!load_dataset_B(Data, pubkey, PP)) {
		printf("Oops! I cannot read this file.\n" );
		return 0;
	}
	
	diff = clock() - start;
	msec = diff * 1000 / CLOCKS_PER_SEC;
	diff_wtime = (float)omp_get_wtime() - wtime;
	printf("\tSUCCESS in %d miliseconds 【%d milisenconds】\n", msec, (int)(diff_wtime*1000));
	
	//Test for load_dataset_B
	/*FILE *fp;
	if((fp = fopen(FILENAME, "rt"))==NULL) {
		return 0;
	}
	int tmp;
	int data[3][NN] = {0};

	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < NUM; j++) {
			if(fscanf(fp,"%d",&tmp))
				data[i][tmp-1] = 1;
		}
	}

	fclose(fp);

	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < NN; j++) {
			lcs_plaintext_t* res = lcs_dec(NULL, PP, pubkey, prvkey, Data[i][j]);
			if(mpz_cmp_ui(res->m, 1) == 0 && data[i][j] != 1) printf("1, (%d, %d)\n", i+1, j+1);
			else if(mpz_cmp_ui(res->m, 0) == 0 && data[i][j] != 0) printf("0, (%d, %d)\n", i+1, j+1);
		}
	}*/

	int len;

	printf("\n***RUNNING SFM\n");
	start = clock();

	Dataset* D = SFM(Data, &len, pubkey, PP, MK);
	
	diff = clock() - start;
	msec = diff * 1000 / CLOCKS_PER_SEC;
	diff_wtime = (float)omp_get_wtime() - wtime;
	printf("\tSUCCESS in %d miliseconds 【%d milisenconds】\n", msec, (int)(diff_wtime*1000));

	printf("\n***PRINT L1 ~ L%d\n", len-1);
	for(int i = 0; i < len-1; i++) {
		printf("\t\t\t频繁项目集L%d如下(共%d个):\n", i+1, D[i].M);
		for(int j = 0; j < D[i].M; j++) {
			printf("{");
			for(int k = 0; k < NN; k++) {
				if(mpz_cmp_ui(D[i].itemset[j].item[k]->m, 1) == 0) {
					printf(" %d ", k+1);
				}
			}
			printf("}\n");
		}
	}

	return 0;

}