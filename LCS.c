#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#include "LCS.h"

void pk_init (lcs_pubkey_t* pk) {
	mpz_init(pk->h);
}

void pv_init (lcs_prvkey_t* pv) {
	mpz_init(pv->a);
}

void pp_init (lcs_pub_para* pp) {
	mpz_init(pp->N);
	mpz_init(pp->g);
	mpz_init(pp->Nsquare);
}

void mk_init (lcs_mkey* mk) {
	mpz_init(mk->lambda);
}

void pt_init (lcs_plaintext_t* pt) {
	mpz_init(pt->m);
}

void ct_init (lcs_ciphertext_t* ct) {
	mpz_init(ct->A);
	mpz_init(ct->B);
}

void init_rand(gmp_randstate_t rand, lcs_get_rand_t get_rand, int bytes) {
	void* buf;
	mpz_t s;

	buf = malloc(bytes);
	get_rand(buf, bytes);

	gmp_randinit_default(rand);
	mpz_init(s);
	mpz_import(s, bytes, 1, 1, 0, 0, buf);
	gmp_randseed(rand, s);
	mpz_clear(s);

	free(buf);
}

void lcs_PP_gen(int modulusbits, lcs_pub_para** PP, lcs_mkey** mk, lcs_get_rand_t get_rand) {
	mpz_t tem, tem1;
	mpz_t N_1;
	mpz_t p, q, p_0, q_0;
	gmp_randstate_t rand;

	/* allocate the new key structures */
	*PP  = (lcs_pub_para*)malloc(sizeof(lcs_pub_para));
	*mk  = (lcs_mkey*)malloc(sizeof(lcs_mkey));

	/* initialize our integers */
	pp_init(*PP);
	mk_init(*mk);
	mpz_init(q_0);
	mpz_init(p_0);
	mpz_init(p);
	mpz_init(q);
	mpz_init(N_1);
	mpz_init(tem);
	mpz_init(tem1);
	(*PP)->bits = modulusbits;

	/* pick random (modulusbits/2)-bit primes p and q */
	init_rand(rand, get_rand, modulusbits / 8 + 1);
	do {
		do
			mpz_urandomb(p, rand, modulusbits / 2);
		while(!mpz_probab_prime_p(p, 10));

		do
			mpz_urandomb(q, rand, modulusbits / 2);
		while(!mpz_probab_prime_p(q, 10));

		// compute the public modulus n = p q 
		mpz_mul((*PP)->N, p, q);
	} while(!mpz_tstbit((*PP)->N, modulusbits - 1));

	mpz_mul((*PP)->Nsquare, (*PP)->N, (*PP)->N);
	mpz_sub_ui(p_0, p, 1);
	mpz_sub_ui(q_0, q, 1);
	mpz_lcm((*mk)->lambda, p_0, q_0); //lambda = lcm(p-1, q-1)

	/*generate g of public parameter*/
	do
		mpz_urandomb(tem, rand, modulusbits);
	while(!mpz_probab_prime_p(tem, 20));
	mpz_mul_ui(N_1, (*PP)->N, 2);
	mpz_powm(tem1, tem, N_1, (*PP)->Nsquare);
	mpz_sub((*PP)->g, (*PP)->Nsquare, tem1); //g=N^2-tem^{2pq} (mod N^2)
	
	mpz_clear(p);
	mpz_clear(q);
	mpz_clear(tem);
	mpz_clear(tem1);
	mpz_clear(p_0);
	mpz_clear(q_0);
	mpz_clear(N_1);
  	gmp_randclear(rand);
}

void lcs_keygen(lcs_pub_para *PP, lcs_pubkey_t** pub, lcs_prvkey_t** prv) {
	mpz_t N_0;
	mpz_init(N_0);

	gmp_randstate_t rand; 
	*pub = (lcs_pubkey_t*) malloc(sizeof(lcs_pubkey_t));
	*prv = (lcs_prvkey_t*) malloc(sizeof(lcs_prvkey_t));
	pk_init(*pub);
	pv_init(*prv);
	gmp_randinit_default(rand); 
	gmp_randseed_ui(rand, time(NULL)); 
	mpz_urandomm((*prv)->a, rand, PP->N);
	mpz_powm((*pub)->h, PP->g, (*prv)->a, PP->Nsquare);

	mpz_clear(N_0);
	gmp_randclear(rand);
}

lcs_ciphertext_t* lcs_enc(lcs_ciphertext_t* res, lcs_pub_para* PP, lcs_pubkey_t* pub, lcs_plaintext_t* pt, lcs_get_rand_t get_rand) {
	mpz_t r;
	gmp_randstate_t rand;
	mpz_t x;

	/* pick random blinding factor */
	mpz_init(r);
	mpz_init(x);
 	init_rand(rand, get_rand, PP->bits / 8 + 1);
	mpz_urandomb(r, rand, PP->bits-3);

	/* compute ciphertext*/ 
	if(!res) {
		res = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
		mpz_init(res->A);
		mpz_init(res->B);
	}
	mpz_powm(res->A, PP->g, r, PP->Nsquare);
	mpz_mul(x, pt->m, PP->N);
	mpz_add_ui(x, x, 1);
	mpz_powm(res->B, pub->h, r, PP->Nsquare);
	mpz_mul(res->B, res->B, x);
	mpz_mod(res->B, res->B, PP->Nsquare);
	
	mpz_clear(x);
	mpz_clear(r);
  	gmp_randclear(rand);

	return res;
}

lcs_plaintext_t* lcs_dec(lcs_plaintext_t* res, lcs_pub_para* PP, lcs_pubkey_t* pub, lcs_prvkey_t* prv, lcs_ciphertext_t* ct) {
	if(!res) {
		res = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
		mpz_init(res->m);
	}

	mpz_t temp;
	mpz_init(temp);
	mpz_invert(temp, ct->A, PP->Nsquare);
	mpz_powm(temp, temp, prv->a, PP->Nsquare);
	mpz_mul(temp, ct->B, temp);
	mpz_sub_ui(temp, temp, 1);
	mpz_mod(temp, temp, PP->Nsquare);
	mpz_div(temp, temp, PP->N);
	mpz_mod(res->m , temp, PP->Nsquare);

	mpz_clear(temp);

	return res;
}

lcs_plaintext_t* lcs_mdec(lcs_plaintext_t* res, lcs_pub_para* PP, lcs_pubkey_t *pk, lcs_mkey* mkey, lcs_ciphertext_t* ct) {
	if(!res) {
		res = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
		mpz_init(res->m);
	}
	
	mpz_t pi, temp;
	mpz_init(pi);
	mpz_init(temp);	

	mpz_powm(temp, ct->B, mkey->lambda, PP->Nsquare);
	mpz_sub_ui(temp, temp, 1);
	mpz_cdiv_q(temp, temp, PP->N);
	mpz_invert(pi, mkey->lambda, PP->N);
	mpz_mul(res->m, temp, pi);
	mpz_mod(res->m, res->m, PP->N);

	mpz_clear(pi);
	mpz_clear(temp);

	return res;
}


void lcs_mul(lcs_pub_para* PP, lcs_ciphertext_t* res, lcs_ciphertext_t* ct0, lcs_ciphertext_t* ct1) {
	if(!res) {
		res = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
		ct_init(res);
	}

	mpz_mul(res->A, ct0->A, ct1->A);
	mpz_mod(res->A, res->A, PP->Nsquare);
	mpz_mul(res->B, ct0->B, ct1->B);
	mpz_mod(res->B, res->B, PP->Nsquare);
}

void lcs_exp(lcs_pub_para* PP, lcs_ciphertext_t* res, lcs_ciphertext_t* ct, lcs_plaintext_t* pt) {
	mpz_powm(res->A, ct->A, pt->m, PP->Nsquare);
	mpz_powm(res->B, ct->B, pt->m, PP->Nsquare);
}

lcs_plaintext_t* lcs_plaintext_from_ui(unsigned long int x) {
	lcs_plaintext_t* pt;

	pt = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	mpz_init_set_ui(pt->m, x);
	
	return pt;
}

lcs_plaintext_t* lcs_plaintext_from_bytes(void* m, int len) {
	lcs_plaintext_t* pt;

	pt = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	mpz_init(pt->m);
	mpz_import(pt->m, len, 1, 1, 0, 0, m);

	return pt;
}

void* lcs_plaintext_to_bytes(int len, lcs_plaintext_t* pt) {
	void* buf0;
	void* buf1;
	size_t written;

	buf0 = mpz_export(0, &written, 1, 1, 0, 0, pt->m);

 	if( written == len )
 		return buf0;

	buf1 = malloc(len);
	memset(buf1, 0, len);

	if( written == 0 )
		/* no need to copy anything, pt->m = 0 and buf0 was not allocated */
		return buf1;
	else if( written < len )
		/* pad with leading zeros */
		memcpy(buf1 + (len - written), buf0, written);
	else
		/* truncate leading garbage */
		memcpy(buf1, buf0 + (written - len), len);

	free(buf0);

	return buf1;
}

lcs_plaintext_t* lcs_plaintext_from_str(char* str) {
	return lcs_plaintext_from_bytes(str, strlen(str));
}

char* lcs_plaintext_to_str(lcs_plaintext_t* pt) {
	char* buf;
	size_t len;

	buf = (char*)mpz_export(0, &len, 1, 1, 0, 0, pt->m);
	buf = (char*)realloc(buf, len + 1);
	buf[len] = 0;

	return buf;
}

/*
lcs_ciphertext_t* lcs_ciphertext_from_bytes(void* c, int len) {
	lcs_ciphertext_t* ct;

	ct = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	mpz_init(ct->c);
	mpz_import(ct->c, len, 1, 1, 0, 0, c);

	return ct;
}

void*  lcs_ciphertext_to_bytes(int len, lcs_ciphertext_t* ct) {
	void* buf;
	int cur_len;

	cur_len = mpz_sizeinbase(ct->c, 2);
	cur_len = PAILLIER_BITS_TO_BYTES(cur_len);
	buf = malloc(len);
	memset(buf, 0, len);
	mpz_export(buf + (len - cur_len), 0, 1, 1, 0, 0, ct->c);

	return buf;
}

char* lcs_pubkey_to_hex(lcs_pubkey_t* pub) {
	return mpz_get_str(0, 16, pub->n);
}

char* lcs_prvkey_to_hex(lcs_prvkey_t* prv) {
	return mpz_get_str(0, 16, prv->lambda);
}

lcs_pubkey_t* lcs_pubkey_from_hex(char* str) {
	lcs_pubkey_t* pub;

	pub = (lcs_pubkey_t*)malloc(sizeof(lcs_pubkey_t));
	mpz_init_set_str(pub->n, str, 16);
	pub->bits = mpz_sizeinbase(pub->n, 2);
	mpz_init(pub->n_squared);
	mpz_init(pub->n_plusone);
	complete_pubkey(pub);

	return pub;
}


lcs_prvkey_t* lcs_prvkey_from_hex( char* str, lcs_pubkey_t* pub )
{
	lcs_prvkey_t* prv;

	prv = (lcs_prvkey_t*) malloc(sizeof(lcs_prvkey_t));
	mpz_init_set_str(prv->lambda, str, 16);
	mpz_init(prv->x);
	complete_prvkey(prv, pub);

	return prv;
}
*/

void lcs_freepubkey(lcs_pubkey_t* pub) {
	mpz_clear(pub->h);
	free(pub);
}

void lcs_freeprvkey(lcs_prvkey_t* prv) {
	mpz_clear(prv->a);
	free(prv);
}

void lcs_freeplaintext(lcs_plaintext_t* pt) {
	mpz_clear(pt->m);
	free(pt);
}

void lcs_freeciphertext(lcs_ciphertext_t* ct) {
	mpz_clear(ct->A);
	mpz_clear(ct->B);
	free(ct);
}

void lcs_freepubpara(lcs_pub_para* PP) {
	mpz_clear(PP->N);
	mpz_clear(PP->Nsquare);
	mpz_clear(PP->g);
	free(PP);
}

void lcs_freemkey(lcs_mkey* mk) {
	mpz_clear(mk->lambda);
	free(mk);
}


void lcs_get_rand_file(void* buf, int len, char* file) {
	FILE* fp;
	void* p;

	fp = fopen(file, "r");

	p = buf;
	while(len) {
		size_t s;
		s = fread(p, 1, len, fp);
		p += s;
		len -= s;
	}

	fclose(fp);
}

void lcs_get_rand_devrandom(void* buf, int len) {
	lcs_get_rand_file(buf, len, "/dev/random");
}

void lcs_get_rand_devurandom(void* buf, int len) {
	lcs_get_rand_file(buf, len, "/dev/urandom");
}
/*
lcs_ciphertext_t* 
lcs_create_enc_zero()
{
	lcs_ciphertext_t* ct;

	// make a NON-RERANDOMIZED encryption of zero for the purposes of homomorphic computation

	// note that this is just the number 1

	//ct = (lcs_ciphertext_t*) malloc(sizeof(lcs_ciphertext_t));
	//mpz_init_set_ui(ct->c, 1);

	//return ct;
//}*/

/***********************
SECRET SHARING BASIC OP
***********************/
void ss_add(lcs_pub_para* PP,  lcs_plaintext_t* res1, lcs_plaintext_t* res2, lcs_plaintext_t* x1, lcs_plaintext_t* x2, lcs_plaintext_t* y1, lcs_plaintext_t* y2) {
	mpz_t z1, z2;
	mpz_init(z1);
	mpz_init(z2);

	mpz_add(z1, x1->m, y1->m);
	mpz_mod(z1, z1, PP->N);
	mpz_add(z2, x2->m, y2->m);
	mpz_mod(z2, z2, PP->N);

	if(!res1) {
		res1 = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
		res2 = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
		pt_init(res1);
		pt_init(res2);
	}

	mpz_set(res1->m, z1);
	mpz_set(res2->m, z2);

	mpz_clear(z1);
	mpz_clear(z2);
}

void ss_mul(lcs_pub_para* PP, lcs_plaintext_t** z, lcs_plaintext_t* x1, lcs_plaintext_t* x2, lcs_plaintext_t* y1, lcs_plaintext_t* y2) {
	mpz_t a, a1, a2;
	mpz_t b, b1, b2;
	mpz_t c, c1, c2;
	mpz_init(a);
	mpz_init(a1);
	mpz_init(a2);
	mpz_init(b);
	mpz_init(b1);
	mpz_init(b2);
	mpz_init(c);
	mpz_init(c1);
	mpz_init(c2);

	gmp_randstate_t rand; 
	gmp_randinit_default(rand); 
	gmp_randseed_ui(rand, time(NULL));

	mpz_urandomb(a1, rand, mpz_sizeinbase(x1->m, 2)-1);
	mpz_urandomb(a2, rand, mpz_sizeinbase(x2->m, 2)-1);
	mpz_urandomb(b1, rand, mpz_sizeinbase(y1->m, 2)-1);
	mpz_urandomb(b2, rand, mpz_sizeinbase(y2->m, 2)-1);

	mpz_add(a, a1, a2);
	mpz_mod(a, a, PP->N);
	
	mpz_add(b, b1, b2);
	mpz_mod(b, b, PP->N);

	mpz_mul(c, a, b);
	mpz_mod(c, c, PP->N);
	mpz_urandomb(c1, rand, mpz_sizeinbase(c, 2)-1);
	mpz_sub(c2, c, c1);

	mpz_t e, e1, e2;
	mpz_t f, f1, f2;
	mpz_init(e);
	mpz_init(e1);
	mpz_init(e2);
	mpz_init(f);
	mpz_init(f1);
	mpz_init(f2);

	mpz_sub(e1, x1->m, a1);
	mpz_sub(e2, x2->m, a2);
	mpz_sub(f1, y1->m, b1);
	mpz_sub(f2, y2->m, b2);
	mpz_add(e, e1, e2);
	mpz_mod(e, e, PP->N);
	mpz_add(f, f1, f2);
	mpz_mod(f, f, PP->N);

	mpz_t z1, z2;
	mpz_init(z1);
	mpz_init(z2);
	
	mpz_mul(z1, f, a1);
	mpz_addmul(z1, e, b1);
	mpz_add(z1, z1, c1);
	mpz_mod(z1, z1, PP->N);

	mpz_mul(z2, f, e);
	mpz_addmul(z2, f, a2);
	mpz_addmul(z2, e, b2);
	mpz_add(z2, z2, c2);
	mpz_mod(z2, z2, PP->N);

	if(!z) {
		z = (lcs_plaintext_t**)malloc(2*sizeof(lcs_plaintext_t*));
		pt_init(z[0]);
		pt_init(z[1]);
	}
	mpz_set(z[0]->m, z1);
	mpz_set(z[1]->m, z2);

	mpz_clear(a);
	mpz_clear(a1);
	mpz_clear(a2);
	mpz_clear(b);
	mpz_clear(b1);
	mpz_clear(b2);
	mpz_clear(c);
	mpz_clear(c1);
	mpz_clear(c2);
	mpz_clear(e);
	mpz_clear(e1);
	mpz_clear(e2);
	mpz_clear(f);
	mpz_clear(f1);
	mpz_clear(f2);
	mpz_clear(z1);
	mpz_clear(z2);
	gmp_randclear(rand);
}

/***********
Section 4.2
***********/
void S2B(lcs_ciphertext_t* res, lcs_plaintext_t* x1, lcs_plaintext_t* x2, lcs_pubkey_t* pubkey, lcs_pub_para* PP) {
	lcs_ciphertext_t* c1 = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	lcs_ciphertext_t* c2 = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(c1);
	ct_init(c2);

	lcs_enc(c1, PP, pubkey, x1, &lcs_get_rand_devurandom);
	lcs_enc(c2, PP, pubkey, x2, &lcs_get_rand_devurandom);

	if(!res) {
		res = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
		ct_init(res);
	}

	lcs_mul(PP, res, c1, c2);

	lcs_freeciphertext(c1);
	lcs_freeciphertext(c2);
}

void B2S(lcs_plaintext_t* res1, lcs_plaintext_t* res2, lcs_ciphertext_t* c, lcs_pubkey_t* pubkey, lcs_pub_para* PP, lcs_mkey *MK) {
	lcs_plaintext_t* x1 = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	lcs_plaintext_t* x2 = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	pt_init(x1);
	pt_init(x2);

	gmp_randstate_t rand; 
	gmp_randinit_default(rand); 
	gmp_randseed_ui(rand, time(NULL));

	mpz_urandomb(x1->m, rand, mpz_sizeinbase(PP->N, 2)-1);

	lcs_ciphertext_t* c1 = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	lcs_ciphertext_t* c2 = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(c1);
	ct_init(c2);

	lcs_plaintext_t* tmp = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	pt_init(tmp);

	mpz_neg(tmp->m, x1->m);
	lcs_enc(c1, PP, pubkey, tmp, &lcs_get_rand_devurandom);
	lcs_mul(PP, c2, c, c1);

	lcs_mdec(x2, PP, pubkey, MK, c2);

	if(!res1) {
		res1 = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
		res2 = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
		pt_init(res1);
		pt_init(res2);
	}

	mpz_set(res1->m, x1->m);
	mpz_set(res2->m, x2->m);

	lcs_freeplaintext(x1);
	lcs_freeplaintext(x2);
	lcs_freeplaintext(tmp);
	gmp_randclear(rand);
	lcs_freeciphertext(c1);
	lcs_freeciphertext(c2);
}

void SIP(lcs_plaintext_t* res1,  lcs_plaintext_t* res2, lcs_plaintext_t** x1, lcs_plaintext_t** x2, lcs_plaintext_t** y1, lcs_plaintext_t** y2, lcs_pub_para* PP, int n) {
	if(!res1) {
		
		res1 = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
		res2 = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
		pt_init(res1);
		pt_init(res2);
	}
	mpz_set_ui(res1->m, 0);
	mpz_set_ui(res2->m, 0);
	//printf("Here\n");

	lcs_plaintext_t** tmp = (lcs_plaintext_t**)malloc(2*sizeof(lcs_plaintext_t*));
	tmp[0] = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	tmp[1] = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	pt_init(tmp[0]);
	pt_init(tmp[1]);

	int i;
	for(i = 0; i < n; i++) {
		ss_mul(PP, tmp, x1[i], x2[i], y1[i], y2[i]);
		ss_add(PP, res1, res2, res1, res2, tmp[0], tmp[1]);
	}

	lcs_freeplaintext(tmp[0]);
	lcs_freeplaintext(tmp[1]);
	free(tmp);
}

void TransDec(lcs_ciphertext_t* res, 	lcs_ciphertext_t* X, lcs_pubkey_t* pk1,
			  lcs_pubkey_t* pk2, lcs_pub_para* PP, 	lcs_mkey* MK)
{
	gmp_randstate_t rand; 
	gmp_randinit_default(rand); 
	gmp_randseed_ui(rand, time(NULL));
	mpz_t tmp;
	mpz_init(tmp);
	lcs_ciphertext_t* ctmp = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(ctmp);
	lcs_plaintext_t* ptmp = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	pt_init(ptmp);

	//Step 1
	mpz_t r1;
	mpz_init(r1);
	mpz_urandomb(r1, rand, mpz_sizeinbase(PP->N, 2)-2);
	mpz_set(ptmp->m, r1);
	lcs_enc(ctmp, PP, pk1, ptmp, &lcs_get_rand_devurandom); //[ra]_pkm
	lcs_mul(PP, ctmp, X, ctmp);

	//Step 2
	lcs_plaintext_t* z = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	pt_init(z); //
	lcs_mdec(z, PP, pk1, MK, ctmp);
	lcs_enc(ctmp, PP, pk2, z, &lcs_get_rand_devurandom);


	//Step 3
	lcs_ciphertext_t* CTMP = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(CTMP);
	lcs_ciphertext_t* CT = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(CT);
	lcs_enc(CTMP, PP, pk2, ptmp, &lcs_get_rand_devurandom);
	mpz_sub_ui(ptmp->m, PP->N, 1);
	lcs_exp(PP, CT, CTMP, ptmp);
	lcs_mul(PP, res, CT, ctmp);

	mpz_clear(tmp);
	mpz_clear(r1);
	lcs_freeplaintext(ptmp);
	lcs_freeplaintext(z);
	lcs_freeciphertext(ctmp);
	lcs_freeciphertext(CTMP);
	lcs_freeciphertext(CT);
}

int SCD(lcs_ciphertext_t* t, lcs_ciphertext_t* x, lcs_ciphertext_t* y,	lcs_pubkey_t* pk1, lcs_pubkey_t* pk2, lcs_pub_para* PP, lcs_mkey *MK) {
	//Init some arguments//
	mpz_t tmp;
	mpz_init(tmp);
	lcs_ciphertext_t* ctmp = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(ctmp);
	lcs_plaintext_t* ptmp = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	pt_init(ptmp);
	gmp_randstate_t rand; 
	gmp_randinit_default(rand); 
	gmp_randseed_ui(rand, time(NULL));

	//Step1-(1)//
	lcs_ciphertext_t* A = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(A);
	mpz_set_ui(ptmp->m, 2);
	lcs_exp(PP,	A, x, ptmp); //[x]_pkm^2
	mpz_set_ui(ptmp->m, 1);
	lcs_enc(ctmp, PP, pk1, ptmp, &lcs_get_rand_devurandom); //[1]_pkm
	lcs_mul(PP,	A, A, ctmp); //[x]_pkm^2 + [1]_pkm

	lcs_ciphertext_t* B = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(B);
	mpz_set_ui(ptmp->m, 2);
	lcs_exp(PP,	B, y, ptmp); //[y]_pks^2

	//Step1-(2)//
	mpz_t a;
	mpz_init(a);
	mpz_urandomb(a, rand, 1);

	lcs_ciphertext_t* C = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(C);
	lcs_ciphertext_t* D = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(D);

	if(mpz_cmp_ui(a, 0)  == 0) {
		mpz_set(C->A, A->A); mpz_set(C->B, A->B);
		mpz_sub_ui(ptmp->m, PP->N, 1);
		lcs_exp(PP, D, B, ptmp);
	} else {
		mpz_sub_ui(ptmp->m, PP->N, 1);
		lcs_exp(PP, C, A, ptmp);
		mpz_set(D->A, B->A); mpz_set(D->B, B->B);
	}

	//Step1-(3)//
	mpz_t r1;
	mpz_init(r1);
	mpz_urandomb(r1, rand, mpz_sizeinbase(PP->N, 2));
	mpz_set(ptmp->m, r1);
	lcs_enc(ctmp, PP, pk1, ptmp, &lcs_get_rand_devurandom); //[ra]_pkm
	lcs_mul(PP, C, C, ctmp); //C*[ra]_pkm

	mpz_t r2;
	mpz_init(r2);
	mpz_urandomb(r2, rand, mpz_sizeinbase(PP->N, 2));
	mpz_set(ptmp->m, r2);
	lcs_enc(ctmp, PP, pk2, ptmp, &lcs_get_rand_devurandom); //[rb]_pks
	lcs_mul(PP, D, D, ctmp); //D*[rb]_pks

	//Step2-(1)//
	lcs_plaintext_t* c_0 = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	pt_init(c_0); //c'
	lcs_mdec(c_0, PP, pk1, MK, C);

	lcs_plaintext_t* d_0 = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	pt_init(d_0); //d'
	lcs_mdec(d_0, PP, pk2, MK, D);

	//Step2-(2)//
	mpz_add(ptmp->m, c_0->m, d_0->m);

	lcs_ciphertext_t* E = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(E);
	lcs_enc(E, PP, pk2, ptmp, &lcs_get_rand_devurandom); //[c'+d']_pks

	//Step3-(1)//
	mpz_add(ptmp->m, r1, r2);
	mpz_neg(ptmp->m, ptmp->m); //-(ra+rb)
	lcs_enc(ctmp, PP, pk2, ptmp, &lcs_get_rand_devurandom); //[-(ra+rb)]_pks

	lcs_ciphertext_t* F = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(F);
	lcs_mul(PP, F, E, ctmp); //F = E*[-(ra+rb)]_pks

	//Step3-(2)//
	mpz_urandomb(ptmp->m, rand, mpz_sizeinbase(PP->N, 2)/2 - 1); //r1
	lcs_exp(PP, F, F, ptmp); //F^{r1}

	mpz_urandomb(ptmp->m, rand, mpz_sizeinbase(ptmp->m, 2)/2); //r2;
	lcs_enc(ctmp, PP, pk2, ptmp, &lcs_get_rand_devurandom); //[r2]_pks

	lcs_mul(PP, F, F, ctmp);

	//Step4-(1)//
	lcs_plaintext_t* z = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	pt_init(z);
	lcs_mdec(z, PP, pk2, MK, F);

	mpz_div_ui(tmp, PP->N, 2);
	if(mpz_cmp(z->m, tmp) < 0) mpz_set_ui(ptmp->m, 1);
	else mpz_set_ui(ptmp->m, 0);
	lcs_ciphertext_t* sigma = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(sigma);
	lcs_enc(sigma, PP, pk2, ptmp, &lcs_get_rand_devurandom);

	//Step5//
	if(!t) {
		t = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
		ct_init(t);
	}

	if(mpz_cmp_ui(a, 0) == 0) {
		mpz_set(t->A, sigma->A);
		mpz_set(t->B, sigma->B);
	} else {
		mpz_set_ui(ptmp->m, 1);
		lcs_enc(t, PP, pk2, ptmp, &lcs_get_rand_devurandom); //[1]_pks
		mpz_sub_ui(ptmp->m, PP->N, 1);
		lcs_exp(PP, ctmp, sigma, ptmp); //[sigma]_pks^{N-1}
		lcs_mul(PP, t, t, ctmp);
	}

	//Final//
	int res = 0;

	lcs_mdec(ptmp, PP, pk2, MK, t);
	if(mpz_cmp_ui(ptmp->m, 1) == 0) {
		res = 1;
	}

	mpz_clear(tmp);
	lcs_freeciphertext(ctmp);
	lcs_freeplaintext(ptmp);
	gmp_randclear(rand);

	lcs_freeciphertext(A);
	lcs_freeciphertext(B);
	mpz_clear(a);
	lcs_freeciphertext(C);
	lcs_freeciphertext(D);
	mpz_clear(r1);
	mpz_clear(r2);
	lcs_freeplaintext(c_0);
	lcs_freeplaintext(d_0);
	lcs_freeciphertext(E);
	lcs_freeciphertext(F);
	lcs_freeplaintext(z);
	lcs_freeciphertext(sigma);

	return res;
}

int SC(lcs_ciphertext_t* t, lcs_ciphertext_t* x, lcs_ciphertext_t* y, lcs_pubkey_t* pk, lcs_pub_para* PP, lcs_mkey *MK) {
	//Init some arguments//
	mpz_t tmp;
	mpz_init(tmp);
	lcs_ciphertext_t* ctmp = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(ctmp);
	lcs_plaintext_t* ptmp = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	pt_init(ptmp);
	gmp_randstate_t rand; 
	gmp_randinit_default(rand); 
	gmp_randseed_ui(rand, time(NULL));

	//Step1-(1)//
	lcs_ciphertext_t* A = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(A);
	mpz_set_ui(ptmp->m, 2);
	lcs_exp(PP,	A, x, ptmp); //[x]_pk^2
	mpz_set_ui(ptmp->m, 1);
	lcs_enc(ctmp, PP, pk, ptmp, &lcs_get_rand_devurandom); //[1]_pk
	lcs_mul(PP,	A, A, ctmp); //[x]_pk^2 + [1]_pk

	lcs_ciphertext_t* B = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(B);
	mpz_set_ui(ptmp->m, 2);
	lcs_exp(PP,	B, y, ptmp); //[y]_pk^2

	//Step1-(2)//
	mpz_t a;
	mpz_init(a);
	mpz_urandomb(a, rand, 1);

	lcs_ciphertext_t* C = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(C);

	if(mpz_cmp_ui(a, 1)  == 0) {
		mpz_sub_ui(ptmp->m, PP->N, 1);
		lcs_exp(PP, C, A, ptmp);
		lcs_mul(PP, C, C, B);
	} else {
		mpz_sub_ui(ptmp->m, PP->N, 1);
		lcs_exp(PP, C, B, ptmp);
		lcs_mul(PP, C, A, C);
	}

	//Step1-(3)//
	lcs_ciphertext_t* D = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(D);

	mpz_urandomb(ptmp->m, rand, mpz_sizeinbase(PP->N, 2)/2 - 1);
	lcs_exp(PP, D, C, ptmp);

	mpz_urandomb(ptmp->m, rand, mpz_sizeinbase(ptmp->m, 2)/2);
	lcs_enc(ctmp, PP, pk, ptmp, &lcs_get_rand_devurandom);

	lcs_mul(PP, D, D, ctmp);

	//Step2-(1)//
	lcs_mdec(ptmp, PP, pk, MK, D);
	mpz_div_ui(tmp, PP->N, 2);
	if(mpz_cmp(ptmp->m, tmp) < 0) mpz_set_ui(ptmp->m, 1);
	else mpz_set_ui(ptmp->m, 0);

	lcs_ciphertext_t* sigma = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(sigma);
	lcs_enc(sigma, PP, pk, ptmp, &lcs_get_rand_devurandom);

	//Step3//
	if(!t) {
		t = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
		ct_init(t);
	}

	if(mpz_cmp_ui(a, 0) == 0) {
		mpz_set(t->A, sigma->A);
		mpz_set(t->B, sigma->B);
	} else {
		mpz_set_ui(ptmp->m, 1);
		lcs_enc(t, PP, pk, ptmp, &lcs_get_rand_devurandom); //[1]_pk
		mpz_sub_ui(ptmp->m, PP->N, 1);
		lcs_exp(PP, ctmp, sigma, ptmp); //[sigma]_pk^{N-1}
		lcs_mul(PP, t, t, ctmp);
	}

	//Final//
	int res = 0;

	lcs_mdec(ptmp, PP, pk, MK, t);
	if(mpz_cmp_ui(ptmp->m, 1) == 0) {
		res = 1;
	}

	mpz_clear(tmp);
	lcs_freeciphertext(ctmp);
	lcs_freeplaintext(ptmp);
	gmp_randclear(rand);

	lcs_freeciphertext(A);
	lcs_freeciphertext(B);
	mpz_clear(a);
	lcs_freeciphertext(C);
	lcs_freeciphertext(D);
	lcs_freeciphertext(sigma);

	return res;
}

/***********
Section 4.3
***********/
int load_dataset_S(Dataset* Data_A, Dataset* Data_B, lcs_pub_para* PP) {
	clock_t start = clock(), diff;

	FILE *fp;
	if((fp = fopen(FILENAME, "rt"))==NULL) {
		return 0;
	}
	int tmp;
	int data[MM][NN] = {0};

	for(int i = 0; i < MM; i++) {
		for(int j = 0; j < NUM; j++) {
			if(fscanf(fp,"%d",&tmp))
				data[i][tmp-1] = 1;
		}
	}

	fclose(fp);

	Data_A->M = MM;
	Data_A->itemset = (Itemset*)malloc(MM * sizeof(Itemset));

	Data_B->M = MM;
	Data_B->itemset = (Itemset*)malloc(MM * sizeof(Itemset));
	
	gmp_randstate_t rand; 
	gmp_randinit_default(rand); 
	gmp_randseed_ui(rand, time(NULL));

	for(int i = 0; i < MM; i++) {
		diff = clock() - start;
		float msec = (float)diff / CLOCKS_PER_SEC;
		printf("\tloading %d in %d with %.3fs (left %.3fs)..\r", i+1, MM, msec, msec/(i+1) * (MM-1-i));

		for(int j = 0; j < NN; j++) {
			Data_A->itemset[i].item[j] = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
			Data_B->itemset[i].item[j] = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));

			pt_init(Data_A->itemset[i].item[j]);
			mpz_urandomb(Data_A->itemset[i].item[j]->m, rand, mpz_sizeinbase(PP->N, 2)-1);

			pt_init(Data_B->itemset[i].item[j]);
			mpz_neg(Data_B->itemset[i].item[j]->m, Data_A->itemset[i].item[j]->m);
			mpz_add_ui(Data_B->itemset[i].item[j]->m, Data_B->itemset[i].item[j]->m, data[i][j]);
		}

		fflush(stdout);
	}

	gmp_randclear(rand);

	printf("                                                                                \r");
	fflush(stdout);

	return 1;
}

int load_dataset_B(lcs_ciphertext_t*** Data, lcs_pubkey_t* pk, lcs_pub_para* PP) {
	clock_t start = clock(), diff;

	FILE *fp;
	if((fp = fopen(FILENAME, "rt"))==NULL) {
		return 0;
	}
	int tmp;
	int data[MM][NN] = {0};

	for(int i = 0; i < MM; i++) {
		for(int j = 0; j < NUM; j++) {
			if(fscanf(fp,"%d",&tmp))
				data[i][tmp-1] = 1;
		}
	}

	fclose(fp);

	double wtime =  omp_get_wtime();

	for(int i = 0; i < MM; i++) {
		diff = clock() - start;
		float msec = (float)diff / CLOCKS_PER_SEC;
		float diff_wtime = (float)omp_get_wtime() - wtime;

		printf("\tloading %d in %d with %.3fs【%.3fs (left %.3fs)】..\r", i+1, MM, msec, diff_wtime, diff_wtime/(i+1) * (MM-1-i));

		#pragma omp parallel for shared(Data) schedule(dynamic)
		for(int j = 0; j < NN; j++) {
			lcs_plaintext_t* ptmp = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
			pt_init(ptmp);

			mpz_set_ui(ptmp->m, data[i][j]);
			Data[i][j] = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
			ct_init(Data[i][j]);
			lcs_enc(Data[i][j], PP, pk, ptmp, &lcs_get_rand_devurandom);

			lcs_freeplaintext(ptmp);
		}

		fflush(stdout);
	}
	
	

	printf("                                                                                \r");
	fflush(stdout);

	return 1;
}

Dataset SF1I(Dataset* Data_A, Dataset* Data_B, lcs_pubkey_t* pk, lcs_pub_para* PP, lcs_mkey *MK) {
	clock_t start = clock(), diff;

	int count = 0;

	lcs_plaintext_t* ptmpA = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	pt_init(ptmpA);
	lcs_plaintext_t* ptmpB = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	pt_init(ptmpB);
	
	lcs_ciphertext_t* min_spt = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(min_spt);
	mpz_set_ui(ptmpA->m, ceil(MIN_SUPPORT * MM));
	lcs_enc(min_spt, PP, pk, ptmpA, &lcs_get_rand_devurandom);
	lcs_ciphertext_t* ctmp = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(ctmp);

	Dataset res;
	res.M = 0;

	double wtime =  omp_get_wtime();
	for(int j = 0; j < NN; j++) {
		diff = clock() - start;
		float msec = (float)diff / CLOCKS_PER_SEC;
		float diff_wtime = (float)omp_get_wtime() - wtime;
		printf("RUNNING SF1I: %d in %d with %.3fs【%.3fs (left %.3fs)】..\r", j+1, NN, msec, diff_wtime, diff_wtime/(j+1) * (NN-1-j));

		mpz_set_ui(ptmpA->m, 0);
		mpz_set_ui(ptmpB->m, 0);
		for(int i = 0; i < MM; i++) {
			ss_add(PP, ptmpA, ptmpB, ptmpA, ptmpB, Data_A->itemset[i].item[j], Data_B->itemset[i].item[j]);
		}
		S2B(ctmp, ptmpA, ptmpB, pk, PP);

		if(SCD(ctmp, ctmp, min_spt, pk, pk, PP, MK) == 1) {
			count++;

			if(count == 1)
				res.itemset = (Itemset*)malloc(sizeof(Itemset));
			else
				res.itemset = (Itemset*)realloc(res.itemset, count * sizeof(Itemset));

			res.itemset[count-1].lambda = 1;

			#pragma omp parallel for shared(res) schedule(dynamic)
			for(int k = 0; k < NN; k++) {
				res.itemset[count-1].item[k] = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
				pt_init(res.itemset[count-1].item[k]);
				mpz_set_ui(res.itemset[count-1].item[k]->m, 0);
			}
			mpz_set_ui(res.itemset[count-1].item[j]->m, 1);
		}

		fflush(stdout);
	}
	res.M = count;

	lcs_freeplaintext(ptmpA);
	lcs_freeplaintext(ptmpB);
	lcs_freeciphertext(min_spt);
	lcs_freeciphertext(ctmp);

	printf("                                                                                \r");
	fflush(stdout);

	return res;
}

lcs_ciphertext_t* SSC(Dataset* Data_A, Dataset* Data_B, Itemset* G_A, Itemset* G_B, lcs_pubkey_t* pk, lcs_pub_para* PP, lcs_mkey *MK) {
	lcs_plaintext_t* ptmp = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	pt_init(ptmp);
	mpz_set_ui(ptmp->m, 0);

	lcs_ciphertext_t* ctmp_c = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(ctmp_c);
	lcs_enc(ctmp_c, PP, pk, ptmp, &lcs_get_rand_devurandom);

	mpz_set_ui(ptmp->m, G_A->lambda);

	lcs_ciphertext_t* ctmp_lambda = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(ctmp_lambda);
	lcs_enc(ctmp_lambda, PP, pk, ptmp, &lcs_get_rand_devurandom);

	#pragma omp parallel for shared(ctmp_c) schedule(dynamic)
	for(int i = 0; i < Data_A->M; i++) {
		lcs_plaintext_t* omp_p1 = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
		pt_init(omp_p1);
		lcs_plaintext_t* omp_p2 = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
		pt_init(omp_p2);
		lcs_ciphertext_t* omp_c = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
		ct_init(omp_c);

		SIP(omp_p1, omp_p2, Data_A->itemset[i].item, Data_B->itemset[i].item, G_A->item, G_B->item, PP, NN);
		S2B(omp_c, omp_p1, omp_p2, pk, PP);
		SC(omp_c, omp_c, ctmp_lambda, pk, PP, MK);

		#pragma omp critical
		{
			lcs_mul(PP, ctmp_c, ctmp_c, omp_c);
		}

		lcs_freeplaintext(omp_p1);
		lcs_freeplaintext(omp_p2);
		lcs_freeciphertext(omp_c);
	}

	lcs_freeplaintext(ptmp);
	lcs_freeciphertext(ctmp_lambda);

	return ctmp_c;
}

Dataset* SFM(lcs_ciphertext_t*** Data, int* n, lcs_pubkey_t* pk, lcs_pub_para* PP, lcs_mkey* MK) {
	printf("\t\t*Step0-INIT\n");

	*n = 0;

	Dataset* Data_A = (Dataset*)malloc(sizeof(Dataset));
	Data_A->M = MM;

	Dataset* Data_B = (Dataset*)malloc(sizeof(Dataset));
	Data_B->M = MM;

	//1: //
	printf("\t\t*Step1-RUNNING B2S\n");

	clock_t start = clock(), diff;

	double wtime =  omp_get_wtime();	
	for(int i = 0; i < MM; i++) {
		diff = clock() - start;
		float msec = (float)diff / CLOCKS_PER_SEC;
		float diff_wtime = (float)omp_get_wtime() - wtime;
		printf("RUNNING B2S: %d in %d with %.3fs【%.3fs (left %.3fs)】..\r", i+1, MM, msec, diff_wtime, diff_wtime/(i+1) * (MM-1-i));

		if(i == 0) {
			Data_A->itemset = (Itemset*)malloc(sizeof(Itemset));
			Data_B->itemset = (Itemset*)malloc(sizeof(Itemset));
		} else {
			Data_A->itemset = (Itemset*)realloc(Data_A->itemset, (i+1) * sizeof(Itemset));
			Data_B->itemset = (Itemset*)realloc(Data_B->itemset, (i+1) * sizeof(Itemset));
		}

		#pragma omp parallel for shared(Data_A, Data_B) schedule(dynamic)
		for(int j = 0; j < NN; j++) {
			Data_A->itemset[i].item[j] = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
			pt_init(Data_A->itemset[i].item[j]);
			Data_B->itemset[i].item[j] = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
			pt_init(Data_B->itemset[i].item[j]);

			B2S(Data_A->itemset[i].item[j], Data_B->itemset[i].item[j],	Data[i][j], pk, PP, MK);
		}

		fflush(stdout);
	}

	printf("                                                                                \r");
	fflush(stdout);

	//2: //
	printf("\t\t*Step2-RUNNING SF1I\n");

	Dataset L = SF1I(Data_A, Data_B, pk, PP, MK);

	if(L.M == 0)
		return NULL;
	else
		(*n)++;

	Dataset* res = (Dataset*)malloc(sizeof(Dataset));
	res[0] = L;
	printf("\t\t\tget L1\n");

	printf("\t\t*Step3~11\n");

	lcs_plaintext_t* ptmp = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
	pt_init(ptmp);
	lcs_ciphertext_t* min_spt = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(min_spt);
	mpz_set_ui(ptmp->m, ceil(MIN_SUPPORT * MM));
	lcs_enc(min_spt, PP, pk, ptmp, &lcs_get_rand_devurandom);

	lcs_ciphertext_t* ctmp = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
	ct_init(ctmp);

	//3: //
	start = clock();

	while(res[(*n)-1].M != 0) {

		//4: //
		Dataset C = apriori_gen(res[(*n)-1]);
		(*n)++;
		res = (Dataset*)realloc(res, (*n) * sizeof(Dataset));
		res[(*n)-1].M = 0;

		Dataset C_A, C_B;
		C_A.M = C_B.M = C.M;
		C_A.itemset = (Itemset*)malloc(C_A.M * sizeof(Itemset));
		C_B.itemset = (Itemset*)malloc(C_B.M * sizeof(Itemset));

		//5: //
		double wtime =  omp_get_wtime();
		for(int i = 0; i < C.M; i++) {
			diff = clock() - start;
			float msec = (float)diff / CLOCKS_PER_SEC;
			float diff_wtime = (float)omp_get_wtime() - wtime;
			printf("RUNNING: %d in %d with %.3fs【%.3fs (left %.3fs)】..\r", i+1, C.M, msec, diff_wtime, diff_wtime/(i+1) * (C.M-1-i));

			//6: //
			C_A.itemset[i].lambda = C_B.itemset[i].lambda = C.itemset[i].lambda;

			#pragma omp parallel for shared(C_A, C_B) schedule(dynamic)
			for(int j = 0; j < NN; j++) {
				lcs_ciphertext_t* omp_ctmp = (lcs_ciphertext_t*)malloc(sizeof(lcs_ciphertext_t));
				ct_init(omp_ctmp);

				lcs_enc(omp_ctmp, PP, pk, C.itemset[i].item[j], &lcs_get_rand_devurandom);

				C_A.itemset[i].item[j] = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
				pt_init(C_A.itemset[i].item[j]);
				C_B.itemset[i].item[j] = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
				pt_init(C_B.itemset[i].item[j]);

				B2S(C_A.itemset[i].item[j], C_B.itemset[i].item[j], omp_ctmp, pk, PP, MK);

				lcs_freeciphertext(omp_ctmp);
			}

			//7: //
			lcs_ciphertext_t* c = SSC(Data_A, Data_B, &C_A.itemset[i], &C_B.itemset[i], pk, PP, MK);

			//8: //
			if(SCD(ctmp, c, min_spt, pk, pk, PP, MK) == 1) {
				//9: //
				res[(*n)-1].M++;
				if(res[(*n)-1].M == 1)
					res[(*n)-1].itemset = (Itemset*)malloc(sizeof(Itemset));
				else
					res[(*n)-1].itemset = (Itemset*)realloc(res[(*n)-1].itemset, res[(*n)-1].M * sizeof(Itemset));

				res[(*n)-1].itemset[res[(*n)-1].M - 1] = C.itemset[i];
			}

			lcs_freeciphertext(c);
			fflush(stdout);
		}

		printf("                                                                                \r");
		fflush(stdout);

		printf("\t\t\tget L%d\n", (*n));
	}

	for(int i = 0; i < Data_A->M; i++) {
		for(int j = 0; j < NN; j++) {
			lcs_freeplaintext(Data_A->itemset[i].item[j]);
			lcs_freeplaintext(Data_B->itemset[i].item[j]);
		}
	}
	free(Data_A->itemset);
	free(Data_B->itemset);
	free(Data_A);
	free(Data_B);
	lcs_freeplaintext(ptmp);
	lcs_freeciphertext(min_spt);
	lcs_freeciphertext(ctmp);

	return res;
}

Dataset apriori_gen(Dataset L) {
	Dataset res;
	int count = 0;

	#pragma omp parallel for shared(L,res, count) schedule(dynamic) ordered
	for(int i = 0; i < L.M-1; i++) {
		for(int j = i+1; j < L.M; j++) {
			int lambda = 0;

			if(L.itemset[0].lambda > 1) {
				 for(int k = 0; k < NN; k++) {
				 	if(mpz_cmp_ui(L.itemset[i].item[k]->m, 1) == 0) {
				 		if(mpz_cmp_ui(L.itemset[j].item[k]->m, 1) == 0) {
				 			lambda++;
				 		} else break;
				 	} else if(mpz_cmp_ui(L.itemset[j].item[k]->m, 1) == 0) break;
				}
			}
			
			#pragma omp critical
			{
				if(lambda == L.itemset[0].lambda-1) {
					Itemset itemset;
					itemset.lambda = L.itemset[0].lambda + 1;
					for(int k = 0; k < NN; k++) {
						itemset.item[k] = (lcs_plaintext_t*)malloc(sizeof(lcs_plaintext_t));
						pt_init(itemset.item[k]);
						mpz_add(itemset.item[k]->m, L.itemset[i].item[k]->m, L.itemset[j].item[k]->m);
						if(mpz_cmp_ui(itemset.item[k]->m, 0) != 0) mpz_set_ui(itemset.item[k]->m, 1);
					}

					count++;
					if(count == 1)
						res.itemset = (Itemset*)malloc(sizeof(Itemset));
					else
						res.itemset = (Itemset*)realloc(res.itemset, count * sizeof(Itemset));

					res.itemset[count-1] = itemset;
				}
			}
		}
	}

	res.M = count;

	return res;
}