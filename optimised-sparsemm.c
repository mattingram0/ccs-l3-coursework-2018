#include "utils.h"
#include <string.h>
#include <stdlib.h>
#include <glib.h>

#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

void basic_sparsemm(const COO, const COO, COO *);
void basic_sparsemm_sum(const COO, const COO, const COO,
		const COO, const COO, const COO,
		COO *);

void add_value_to_sparse(gpointer coords, gpointer value, gpointer sp){
	static int counter = 0;
	char* delim = ",";
	char* stri = strtok((char*)coords, delim); //Split up our string key "i,j" into i and j
	char* strj = strtok(NULL, delim);		   
	int i = atoi(stri); 					   //Convert out i, j strings into integers - sorry for using atoi NOTE
	int j = atoi(strj);
	
	((COO)sp)->coords[counter].i = i;
	((COO)sp)->coords[counter].j = j;
	((COO)sp)->data[counter] = *((double *)value);
	
	counter++;
}

/* Converts our hashmap to COO formulated array */
void convert_hashmap_to_sparse(GHashTable* hash, int m, int n, int NZ, COO* C){
	COO sp;
    alloc_sparse(m, n, NZ, &sp);

	//Loop through every key, value pair in our matrix, calling 'add_value_to_sparse' on each one, passing sp to the function
	g_hash_table_foreach(hash, add_value_to_sparse, sp);

	*C = sp;
}


/* Computes C = A*B.
 * C should be allocated by this routine.
 */

void optimised_sparsemm(const COO A, const COO B, COO *C)
{

	LIKWID_MARKER_START("Optimised Sparsemm - Whole Function");
	int a, b, nzA, nzB, nzC, m, n, estimate; 
	double product, *partial = NULL;
	char coords[40];	//NOTE - this may be vulnerable to buffer of if mat large
	GHashTable* sparseC = g_hash_table_new(g_str_hash, g_str_equal);

	gchar *key;
	gdouble *memPtr, *memPool;

	nzA = A->NZ;
	nzB = B->NZ;
	m = A->m;
	n = B->n;
	*C = NULL;
	nzC = product = 0;
	//density = (nzA / (A->m * A->n)) + (nzB / (B->m * B->n));
	//estimate = density * A->n * B->m; 
	estimate = (nzA * B->m / A->n) + (nzB * A->n / B->n);
	memPool = (double *) malloc(sizeof(double) * estimate); 
	memPtr = memPool;

	LIKWID_MARKER_START("Optimised Sparsemm - Loops");

	for(a = 0; a < nzA; a++){ //NOTE - time complexity is O(nzA . nzB) - much better than O(n^3)
		for(b = 0; b < nzB; b++){
			if(A->coords[a].j == B->coords[b].i){
				product = A->data[a] * B->data[b];
				sprintf(coords, "%d,%d", A->coords[a].i, B->coords[b].j); //Turn our coordinates into a string separated with a , - NOTE BO poss
//
				LIKWID_MARKER_START("Optimised Sparsemm - Hash Table Access");
				partial = (double *) g_hash_table_lookup(sparseC, coords);
			
				if(partial == NULL){ //NOTE - Optimised as only create the memory for new entries - space requirement is O(nzC) - linear
					if(nzC == estimate){ //If we have exceeded the number of non zero elements we estimated, then allocate more memory and realign pointers
						estimate += ((nzA - a)/nzA + 0.01) * estimate;
						memPool = (double *) realloc(memPool, sizeof(double) * estimate);
						memPtr = memPool + nzC; 
						//realloc more memory
					}
					key = g_strdup(coords);									  
					*memPtr = product;		
					
					g_hash_table_insert(sparseC, key, memPtr);
					nzC++;
					memPtr++;
				}else{
					*partial += product;
				}
			LIKWID_MARKER_STOP("Optimised Sparsemm - Hash Table Access");
			}
		}
	}
	LIKWID_MARKER_STOP("Optimised Sparsemm - Loops");
	convert_hashmap_to_sparse(sparseC, m, n, nzC, C);
	LIKWID_MARKER_STOP("Optimised Sparsemm - Whole Function");
	free(memPool);
	g_hash_table_destroy(sparseC);
	return;
}


/* Computes O = (A + B + C) (D + E + F).
 * O should be allocated by this routine.
 */
void optimised_sparsemm_sum(const COO A, const COO B, const COO C,
		const COO D, const COO E, const COO F,
		COO *O)
{
	return basic_sparsemm_sum(A, B, C, D, E, F, O);
}
