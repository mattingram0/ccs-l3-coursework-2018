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

	LIKWID_MARKER_START("OptimisedSMM");
	int a, b, nzA, nzB, nzC, m, n; 
	double product, *partial = NULL;
	char coords[40];	//NOTE - this may be vulnerable to buffer of if mat large
	GHashTable* sparseC = g_hash_table_new_full(g_str_hash, g_str_equal, free, free);

	gchar *key;
	gdouble *value;

	nzA = A->NZ;
	nzB = B->NZ;
	m = A->m;
	n = B->n;
	*C = NULL;
	nzC = product = 0;

	LIKWID_MARKER_START("Optimised dgemm");

	for(a = 0; a < nzA; a++){ //NOTE - time complexity is O(nzA . nzB) - much better than O(n^3)
		for(b = 0; b < nzB; b++){
			if(A->coords[a].j == B->coords[b].i){
				product = A->data[a] * B->data[b];
				sprintf(coords, "%d,%d", A->coords[a].i, B->coords[b].j); //Turn our coordinates into a string separated with a , - NOTE BO poss
				partial = (double *) g_hash_table_lookup(sparseC, key);
			
				if(partial == NULL){ //NOTE - Optimised as only create the memory for new entries - space requirement is O(nzC) - linear
					key = g_strdup(coords);									  
					value = (double *) malloc(sizeof(double));
					*value = product;
					
					g_hash_table_insert(sparseC, key, value);
					nzC++;
				}else{
					*partial += product;
				}
			}
		}
	}

	LIKWID_MARKER_STOP("Optimised dgemm");
	convert_hashmap_to_sparse(sparseC, m, n, nzC, C);
	g_hash_table_destroy(sparseC);
	LIKWID_MARKER_STOP("OptimisedSMM");
	return;
}


/* Computes O = (A + B + C) (D + E + F).
 * O should be allocated by this routine.
 */
void sparsemm_sum_2(const COO A, const COO B, COO *S)
{

	int a, b, nzA, nzB, nzS, m, n; 
	double *partial = NULL;
	char coords[40];	
	GHashTable* sparseS = g_hash_table_new_full(g_str_hash, g_str_equal, free, free);//Use the full version, specifying destructors

	gchar *key;
	gdouble *value;

	nzA = A->NZ;
	nzB = B->NZ;
	m = A->m;
	n = B->n;
	*S = NULL;
	nzS = 0;
	
	for(a = 0; a < nzA; a++){
		sprintf(coords, "%d,%d", A->coords[a].i, A->coords[a].j);
		value = (double *) malloc(sizeof(double));
		*value = A->data[a];
		key = g_strdup(coords);									  
		g_hash_table_insert(sparseS, key, value);	
		nzS++;
	}

	for(b = 0; b < nzB; b++){
		sprintf(coords, "%d,%d", B->coords[b].i, B->coords[b].j);
		partial = (double *) g_hash_table_lookup(sparseS, key);

		if(partial == NULL){
			value = (double *) malloc(sizeof(double));
			*value = B->data[b];
			key = g_strdup(coords);
			g_hash_table_insert(sparseS, key, value);	
			nzS++;
		} else {
			*pointer += B->data[b];
		}
	}

	convert_hashmap_to_sparse(sparseS, m, n, nzS, S);
	g_hash_table_destroy(sparseS);//Ensure that our potentially large hash table is freed from memory - NOTE
	return;
}

void sparsemm_sum_3(const COO A, const COO B, const COO C, COO *S)
{

	int a, b, c, nzA, nzB, nzC, nzS, m, n; 
	double *partial = NULL;
	char coords[40];	
	GHashTable* sparseS = g_hash_table_new_full(g_str_hash, g_str_equal, free, free);//Use the full version, specifying destructors

	gchar *key;
	gdouble *value;

	nzA = A->NZ;
	nzB = B->NZ;
	nzC = C->NZ;
	m = A->m;
	n = B->n;
	*S = NULL;
	nzS = 0;
	
	for(a = 0; a < nzA; a++){
		sprintf(coords, "%d,%d", A->coords[a].i, A->coords[a].j);
		value = (double *) malloc(sizeof(double));
		*value = A->data[a];
		key = g_strdup(coords);									  
		g_hash_table_insert(sparseS, key, value);	
		nzS++;
	}

	for(b = 0; b < nzB; b++){
		sprintf(coords, "%d,%d", B->coords[b].i, B->coords[b].j);
		partial = (double *) g_hash_table_lookup(sparseS, key);

		if(partial == NULL){
			value = (double *) malloc(sizeof(double));
			*value = B->data[b];
			key = g_strdup(coords);
			g_hash_table_insert(sparseS, key, value);	
			nzS++;
		} else {
			*pointer += B->data[b];
		}
	}
	for(c = 0; c < nzC; c++){
		sprintf(coords, "%d,%d", C->coords[c].i, C->coords[c].j);
		partial = (double *) g_hash_table_lookup(sparseS, key);

		if(partial == NULL){
			value = (double *) malloc(sizeof(double));
			*value = C->data[c];
			key = g_strdup(coords);
			g_hash_table_insert(sparseS, key, value);	
			nzS++;
		} else {
			*pointer += C->data[c];
		}
	}

	convert_hashmap_to_sparse(sparseS, m, n, nzS, S);
	g_hash_table_destroy(sparseS);//Ensure that our potentially large hash table is freed from memory - NOTE
	return;
}

void optimised_sparsemm_sum(const COO A, const COO B, const COO C,
		const COO D, const COO E, const COO F,
		COO *O)
{
	
	return basic_sparsemm_sum(A, B, C, D, E, F, O);
}
