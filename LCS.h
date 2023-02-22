/******
 TYPES
*******/ 
typedef struct {
	mpz_t h;
} lcs_pubkey_t;

typedef struct {
	mpz_t a;
} lcs_prvkey_t;

typedef struct {
	mpz_t N; // N=pq, p = 2p'+1, q = 2q'+1
	mpz_t g;
	mpz_t Nsquare;
	int bits;  /* e.g., 1024 */
} lcs_pub_para;

typedef struct {
	mpz_t lambda;
} lcs_mkey;

typedef struct {
	mpz_t m;
} lcs_plaintext_t;

typedef struct {
	mpz_t A;
	mpz_t B;
} lcs_ciphertext_t;

typedef void (*lcs_get_rand_t) (void* buf, int len);

/*****
 INIT
*****/
void pk_init (lcs_pubkey_t* pk);
void pv_init (lcs_prvkey_t* pv);
void pp_init (lcs_pub_para* pp);
void mk_init (lcs_mkey* mk);
void pt_init (lcs_plaintext_t* pt);
void ct_init (lcs_ciphertext_t* ct);

/*****************
 BASIC OPERATIONS
*****************/
void lcs_PP_gen (int modulusbits, lcs_pub_para** PP, lcs_mkey** mk, lcs_get_rand_t get_rand);

void lcs_keygen (lcs_pub_para *PP, lcs_pubkey_t** pub, lcs_prvkey_t** prv);

lcs_ciphertext_t* lcs_enc (lcs_ciphertext_t* res,
													lcs_pub_para* PP,
													lcs_pubkey_t* pub,
													lcs_plaintext_t* pt,
													lcs_get_rand_t get_rand);

lcs_plaintext_t* lcs_dec (lcs_plaintext_t* res,
													lcs_pub_para* PP,
													lcs_pubkey_t* pub,
													lcs_prvkey_t* prv,
													lcs_ciphertext_t* ct);

lcs_plaintext_t* lcs_mdec( lcs_plaintext_t* res,
													lcs_pub_para* PP,
													lcs_pubkey_t* pk,
													lcs_mkey* mkey,
													lcs_ciphertext_t* ct);


/*****************************
 USE OF ADDITIVE HOMOMORPHISM
*****************************/
void lcs_mul (lcs_pub_para* PP,
									lcs_ciphertext_t* res,
									lcs_ciphertext_t* ct0,
									lcs_ciphertext_t* ct1);

void lcs_exp (lcs_pub_para* PP,
									lcs_ciphertext_t* res,
									lcs_ciphertext_t* ct,
									lcs_plaintext_t* pt);

/****************************
 PLAINTEXT IMPORT AND EXPORT
****************************/
lcs_plaintext_t* lcs_plaintext_from_ui (unsigned long int x);
lcs_plaintext_t* lcs_plaintext_from_bytes (void* m, int len);
lcs_plaintext_t* lcs_plaintext_from_str (char* str);

char* lcs_plaintext_to_str (lcs_plaintext_t* pt);
void* lcs_plaintext_to_bytes (int len, lcs_plaintext_t* pt);

/*****************************
 CIPHERTEXT IMPORT AND EXPORT
*****************************/
lcs_ciphertext_t* lcs_ciphertext_from_bytes (void* c, int len);
void* lcs_ciphertext_to_bytes (int len, lcs_ciphertext_t* ct);

/**********************
 KEY IMPORT AND EXPORT
**********************/
char* lcs_pubkey_to_hex (lcs_pubkey_t* pub);
char* lcs_prvkey_to_hex (lcs_prvkey_t* prv);
lcs_pubkey_t* lcs_pubkey_from_hex (char* str);
lcs_prvkey_t* lcs_prvkey_from_hex (char* str, lcs_pubkey_t* pub);

/********
 CLEANUP
********/

void lcs_freepubpara (lcs_pub_para* PP);
void lcs_freemkey (lcs_mkey* mk);
void lcs_freepubkey (lcs_pubkey_t* pub);
void lcs_freeprvkey (lcs_prvkey_t* prv);
void lcs_freeplaintext (lcs_plaintext_t* pt);
void lcs_freeciphertext (lcs_ciphertext_t* ct);

/***********
 MISC STUFF
***********/
void lcs_get_rand_devrandom (void* buf, int len);
void lcs_get_rand_devurandom (void* buf, int len);

lcs_ciphertext_t* lcs_create_enc_zero();

#define PAILLIER_BITS_TO_BYTES(n) ((n) % 8 ? (n) / 8 + 1 : (n) / 8)


/***********************
SECRET SHARING BASIC OP
***********************/
void ss_add(lcs_pub_para* PP,		
									lcs_plaintext_t* res1, lcs_plaintext_t* res2,
									lcs_plaintext_t* x1, lcs_plaintext_t* x2,
									lcs_plaintext_t* y1, lcs_plaintext_t* y2);
void ss_mul(lcs_pub_para* PP,		
									lcs_plaintext_t** z,
									lcs_plaintext_t* x1, lcs_plaintext_t* x2,
									lcs_plaintext_t* y1, lcs_plaintext_t* y2);

/***
4.2
***/
//where c = [x], c1=[x1], c2=[x2]
void S2B(lcs_ciphertext_t* c,
									lcs_plaintext_t* x1, lcs_plaintext_t* x2,
									lcs_pubkey_t* pubkey, lcs_pub_para* PP);

void B2S(lcs_plaintext_t* x, lcs_plaintext_t* y,
									lcs_ciphertext_t* c, lcs_pubkey_t* pubkey,
									lcs_pub_para* PP, lcs_mkey* MK);

void SIP(lcs_plaintext_t* res1,		lcs_plaintext_t* res2, 
									lcs_plaintext_t** x1, lcs_plaintext_t** x2,
									lcs_plaintext_t** y1, lcs_plaintext_t** y2,
									lcs_pub_para* PP, int n);
int SCD(lcs_ciphertext_t* t,
									lcs_ciphertext_t* x, lcs_ciphertext_t* y,
									lcs_pubkey_t* pk1, lcs_pubkey_t* pk2,
									lcs_pub_para* PP, lcs_mkey* MK);
int SC(lcs_ciphertext_t* t,
									lcs_ciphertext_t* x, lcs_ciphertext_t* y,
									lcs_pubkey_t* pubkey,
									lcs_pub_para* PP, lcs_mkey* MK);

void TransDec(lcs_ciphertext_t* res,
									lcs_ciphertext_t* x,
									lcs_pubkey_t* pk1, lcs_pubkey_t* pk2,
									lcs_pub_para* PP, lcs_mkey* MK);

/***********
4.3 Apriori
***********/
#define NN 75				//数据集中的最大的数字
#define MM 3196				//读入的行数
#define NUM 37				//读入的列数
#define MIN_SUPPORT 0.97		//最小支持度
#define FILENAME "chess.dat"	//读入的文件名

typedef struct {
	int lambda;
	lcs_plaintext_t* item[NN];
} Itemset;

typedef struct {
	int M;
	Itemset* itemset;
} Dataset;

int load_dataset_S(Dataset* Data_A,
				   Dataset* Data_B,
										lcs_pub_para* PP);

int load_dataset_B(lcs_ciphertext_t*** Data,
										lcs_pubkey_t* pk,
										lcs_pub_para* PP);

int load_dataset_S(Dataset* Data_A,
				   Dataset* Data_B,
										lcs_pub_para* PP);

Dataset SF1I(Dataset* Data_A,
			 Dataset* Data_B,
										lcs_pubkey_t* pk,
			   							lcs_pub_para* PP, lcs_mkey *MK);

lcs_ciphertext_t* SSC(Dataset* A_Data,
					  Dataset* B_Data,
					  Itemset* G_A,
					  Itemset* G_B,
					  					lcs_pubkey_t* pk,
										lcs_pub_para* PP, lcs_mkey *MK);

Dataset* SFM(lcs_ciphertext_t*** Data,
			 int* n,
										lcs_pubkey_t* pk,
										lcs_pub_para* PP, lcs_mkey* MK);

Dataset apriori_gen(Dataset L);