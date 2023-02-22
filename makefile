SFM:SFM.c BCP.c BCP.h
	gcc -O3 -march=native -mtune=native -Wall -fopenmp SFM.c BCP.c -lgmp -o SFM
