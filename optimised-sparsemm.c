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

void free_memory(gpointer ptr, gpointer userData){
	free(ptr);
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
int col_sort(const void *a, const void *b){
	if(((int *)a)[1] < ((int *)b)[1]){
		return -1;
	}else if(((int *)a)[1] == ((int *)b)[1]){
		if(((int *)a)[0] < ((int *)b)[0]) 
			return -1;
		return 1;
	}else{
		return 1;
	}
}

int row_sort(const void *a, const void *b){
	if(((int *)a)[0] < ((int *)b)[0]){
		return -1;
	}else if(((int *)a)[0] == ((int *)b)[0]){
		if(((int *)a)[1] < ((int *)b)[1]) 
			return -1;
		return 1;
	}else{
		return 1;
	}
}

void sort_coo(const COO A, COO* S, int compare(const void *, const void *)){
	int nz, size, c;
	void *unsorted, *head;
	COO sp;
	nz = A->NZ;
	size = 2 * sizeof(int) + sizeof(double);
	c = 0;

	unsorted = malloc(nz * size);
	head = unsorted;

	//Create array of arrays to be sorted
	for(c; c < nz; c++){
		*(int *)head = A->coords[c].i;
		head = (int *)head + 1;
		*(int *)head = A->coords[c].j;
		head = (int *)head + 1;
		*(double *)head = A->data[c];
		head = (double *)head + 1;
	}

	//Sort using our col sort comparator
	qsort(unsorted, nz, size, compare);

	//Allocate a new sparse COO object
	alloc_sparse(A->m, A->n, nz, &sp);

	c = 0;
	head = unsorted;

	//Populate our sorted COO object
	for(c; c < nz; c++){
		sp->coords[c].i = *(int *)head;	
		head = (int *)head + 1;
		sp->coords[c].j = *(int *)head;	
		head = (int *)head + 1;
		sp->data[c] = *(double *)head;
		head = (double *)head + 1;
	}

	*S = sp;
	return;
}

void optimised_sparsemm(COO A, const COO B, COO *C)
{
	LIKWID_MARKER_START("Optimised Sparsemm - Whole Function");
	int a, b, nzA, nzB, nzC = 0, m, n, estimate, newEstimate, currentCol = 0, beginRow = 0, nextRow = 0;
	double product, errorProp, *partial = NULL;
	char coords[15];	//NOTE - this may be vulnerable to buffer of if mat large
	GHashTable* sparseC = g_hash_table_new(g_str_hash, g_str_equal);
	COO S;
	gchar *key;
	gdouble *memPtr, *newPtr;
	GPtrArray *memPtrs;
	
	//Array to keep track of memory pointers - could not simply use realloc as this would move the data pointed to by the pointers
	//in the hash table
	memPtrs = g_ptr_array_new();

	nzA = A->NZ;
	nzB = B->NZ;
	m = A->m;
	n = B->n;
	*C = NULL;
	estimate = (int)(((double)nzA * (double)B->m /(double)A->n) + ((double)nzB * (double)A->n / (double)B->n));
	memPtr = (double *) malloc(sizeof(double) * estimate); 
	g_ptr_array_add(memPtrs, (gpointer)memPtr);	

	//Sort by column, freeing our old matrix, then setting A to our sorted matrix to ensure it gets freed later in sparsemm.c
	sort_coo(A, &S, col_sort);
	free_sparse(&A);

	for(a = 0; a < nzA; a++){
		if(S->coords[a].j != currentCol){
			currentCol = S->coords[a].j;
			beginRow = nextRow;
		}

		b = beginRow;

		while(b < nzB && B->coords[b].i == currentCol){
			product = S->data[a] * B->data[b];
			sprintf(coords, "%d,%d", S->coords[a].i, B->coords[b].j);
			partial = (double *) g_hash_table_lookup(sparseC, coords);

			if(partial == NULL){ 
				if(nzC == estimate){ 
					errorProp = ((double)nzA * (double)nzB - (((double)a * (double)nzB) + (double)b))/(((double)a * (double)nzB) + (double)b);
					newEstimate = (int)(errorProp * estimate);
					memPtr = (double *) malloc(sizeof(double) * newEstimate);
					estimate += newEstimate;
					g_ptr_array_add(memPtrs, (gpointer)memPtr);		
				}

				key = g_strdup(coords);									  
				*memPtr = product;		

				g_hash_table_insert(sparseC, key, memPtr);
				nzC++;
				memPtr++;
			}else{
				*partial += product;
			}

			b++;	
		}

		nextRow = b;
	}

	//LIKWID_MARKER_STOP("Optimised Sparsemm - Loops");
	convert_hashmap_to_sparse(sparseC, m, n, nzC, C);
	LIKWID_MARKER_STOP("Optimised Sparsemm - Whole Function");
	g_ptr_array_foreach(memPtrs, free_memory, NULL);
	free_sparse(&S);
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
