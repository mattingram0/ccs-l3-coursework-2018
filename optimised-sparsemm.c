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

#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

void basic_sparsemm(const COO, const COO, COO *);
void basic_sparsemm_sum(const COO, const COO, const COO,
		const COO, const COO, const COO,
		COO *);

void add_value_to_sparse(gpointer coords, gpointer value, gpointer sp){
	char* delim = ",";
	char* stri = strtok((char*)coords, delim); //Split up our string key "i,j" into i and j
	char* strj = strtok(NULL, delim);

	int i = atoi(stri); 					   //Convert out i, j strings into integers - sorry for using atoi NOTE
	int j = atoi(strj);

	((COO)sp)->coords[((COO)sp)->NZ].i = i;
	((COO)sp)->coords[((COO)sp)->NZ].j = j;
	((COO)sp)->data[((COO)sp)->NZ] = *((double *)value);
	//printf("[%d, %d, %f]\n", i, j, *((double *)value));
	((COO)sp)->NZ++;
}

void print_ptr(gpointer ptr, gpointer userData){
	printf("Address: %p, Value: %p, Points to: %f\n", &ptr, ptr, *((double *)ptr));
}

void free_memory(gpointer ptr, gpointer userData){
	free(ptr);
}

/* Converts our hashmap to COO formulated array */
void convert_hashmap_to_sparse(GHashTable* hash, int m, int n, int NZ, COO* C){
	COO sp;
	alloc_sparse(m, n, NZ, &sp);

	//Set NZ = 0 to use as a counter for our add_value_to_sparse function
	sp->NZ = 0;
	//Loop through every key, value pair in our matrix, calling 'add_value_to_sparse' on each one, passing sp to the function
	g_hash_table_foreach(hash, add_value_to_sparse, sp);
	*C = sp;
	//printf("\nJust after convert_hashmap has *C = sp\n");
	//printf("Data 0: %f\n", (*C)->data[0]);
	//printf("Data 1: %f\n", (*C)->data[1]);
	//printf("Data 2: %f\n", (*C)->data[2]);
}


/* Computes C = A*B.
 * C should be allocated by this routine.
 */
//VEC?
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

	//VEC - Loop/Unravell
	//WASN'T VECTORISED
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

	//VEC - LoopUnravell
	//WAS VECTORIZED
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

void optimised_sparsemm(COO A, COO B, COO *C)
{
	LIKWID_MARKER_START("Optimised Sparsemm Sum - Multiplication");
	//LIKWID_MARKER_START("Optimised Sparsemm - Whole Function");
	//printf("\n In OPTS\n");
	//printf("A Data 1: %f\n", A->data[0]);
	//printf("A Data 2: %f\n", A->data[1]);
	//printf("A Data 3: %f\n", A->data[2]);
	//printf("A Data 4: %f\n", A->data[3]);
	//printf("A Data 5: %f\n", A->data[4]);
	//printf("A Data 6: %f\n", A->data[5]);
	//printf("A Data 7: %f\n", A->data[6]);
	//printf("B Data 1: %f\n", B->data[0]);
	//printf("B Data 2: %f\n", B->data[1]);
	//printf("B Data 3: %f\n", B->data[2]);
	//printf("B Data 4: %f\n", B->data[3]);
	//printf("B Data 5: %f\n", B->data[4]);
	//printf("B Data 6: %f\n", B->data[5]);
	//printf("B Data 7: %f\n", B->data[6]);
	int a, b, nzA, nzB, nzC = 0, m, n, estimate, newEstimate, currentCol = 0, beginRow = 0;
	double product, errorProp, *partial = NULL;
	char coords[15];	//NOTE - this may be vulnerable to buffer of if mat large
	GHashTable* sparseC = g_hash_table_new(g_str_hash, g_str_equal);

	//Create COO objects for the sorted A and B matrices
	COO AS, BS;
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
	sort_coo(A, &AS, col_sort);
	//Sort by row
	sort_coo(B, &BS, row_sort);

	//NOT VECTORISED - multiple nested loops
	//SLP doesn't divide the vector size. Unknown alignment for acess
	for(a = 0; a < nzA; a++){
		//Move our b back down to the beginning of the row which matches our column
		b = beginRow;

		//Skip values in our A matrix until our a column matches our b row
		//NOT VECTORISED - control flow in loop
		while(a < nzA && AS->coords[a].j < BS->coords[b].i){
			a++;
		}

		//Skip values in our B matrix until our b row matches our a column
		//NOT VECTORISED - control flow in loop
		while(b < nzB && BS->coords[b].i < AS->coords[a].j){
			b++;
		}

		beginRow = b;
		currentCol = AS->coords[a].j;

		//NOT VECTORISED - control flow in loop
		while(b < nzB && BS->coords[b].i == currentCol){
			product = AS->data[a] * BS->data[b];
			//printf("A Data: %f, B Data: %f\n", AS->data[a], BS->data[a]);
			sprintf(coords, "%d,%d", AS->coords[a].i, BS->coords[b].j);
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
	}

	//LIKWID_MARKER_STOP("Optimised Sparsemm - Loops");
	convert_hashmap_to_sparse(sparseC, m, n, nzC, C);
	//LIKWID_MARKER_STOP("Optimised Sparsemm - Whole Function");
	LIKWID_MARKER_STOP("Optimised Sparsemm Sum - Multiplication");
	g_ptr_array_foreach(memPtrs, free_memory, NULL);
	free_sparse(&AS);
	free_sparse(&BS);
	g_hash_table_destroy(sparseC);
	return;
}


/* Computes O = (A + B + C) (D + E + F).
 * O should be allocated by this routine.
 */
void add_3(const COO A, const COO B, const COO C, COO *S){

	LIKWID_MARKER_START("Optimised Sparsemm Sum - Addition");
	int a, b, c, nzA, nzB, nzC, nzS, estimate, newEstimate;
	double errorProp, value = 0.0, *partial = NULL;
	char coords[15];	//NOTE - this may be vulnerable to buffer of if mat large
	GHashTable* sparseS = g_hash_table_new(g_str_hash, g_str_equal);
	COO big, mid, sml, temp;
	COO order[3];
	gchar *key;
	gdouble *memPtr, *newPtr;
	GPtrArray *memPtrs;

	nzA = A->NZ;
	nzB = B->NZ;
	nzC = C->NZ;

	memPtrs = g_ptr_array_new();
	order[0] = A;
	order[1] = B;
	order[2] = C;
	//printf("\nnzA: %d, nzB: %d, nzC: %d\n", order[0]->NZ, order[1]->NZ, order[2]->NZ);

	if(order[0]->NZ > order[1]->NZ){
		order[0] = B;
		order[1] = A;
	}

	if(order[0]->NZ > order[2]->NZ){
		temp = order[0];
		order[0] = order[2];
		order[2] = temp;
	}

	if(order[1]->NZ > order[2]->NZ){
		temp = order[1];
		order[1] = order[2];
		order[2] = temp;
	}

	big = order[2];
	mid = order[1];
	sml = order[0];
	//printf("Big: %d, Medium: %d, Small: %d\n", big->NZ, mid->NZ, sml->NZ);
	*S = NULL;
	estimate = nzS = big->NZ;
	memPtr = (double *) malloc(sizeof(double) * estimate); 
	g_ptr_array_add(memPtrs, (gpointer)memPtr);	

	for(a = 0; a < big->NZ; a++){
		value =	big->data[a];
		sprintf(coords, "%d,%d", big->coords[a].i, big->coords[a].j);
		key = g_strdup(coords);									  
		*memPtr = value;
		g_hash_table_insert(sparseS, key, memPtr);
		memPtr++;
	}

	//NOT VEC - control flow in loop
	for(b = 0; b < mid->NZ; b++){
		value = mid->data[b];
		sprintf(coords, "%d,%d", mid->coords[b].i, mid->coords[b].j);
		partial = (double *) g_hash_table_lookup(sparseS, coords);

		if (partial == NULL){
			if(nzS == estimate){ 
				newEstimate = mid->NZ;
				memPtr = (double *) malloc(sizeof(double) * newEstimate);
				estimate += newEstimate;
				g_ptr_array_add(memPtrs, (gpointer)memPtr);		
			}

			key = g_strdup(coords);									  
			*memPtr = value;		
			g_hash_table_insert(sparseS, key, memPtr);
			
			nzS++;
			memPtr++;
		} else {
			*partial += value;
		}
	}

	//NOT VEC - control flow in loop
	for(c = 0; c < sml->NZ; c++){
		value = sml->data[c];
		sprintf(coords, "%d,%d", sml->coords[c].i, sml->coords[c].j);
		partial = (double *) g_hash_table_lookup(sparseS, coords);

		if (partial == NULL){
			if(nzS == estimate){ 
				newEstimate = sml->NZ;
				memPtr = (double *) malloc(sizeof(double) * newEstimate);
				estimate += newEstimate;
				g_ptr_array_add(memPtrs, (gpointer)memPtr);		
			}

			key = g_strdup(coords);									  
			*memPtr = value;		
			g_hash_table_insert(sparseS, key, memPtr);

			nzS++;
			memPtr++;
		} else {
			*partial += value;
		}
	}

	convert_hashmap_to_sparse(sparseS, A->m, A->n, nzS, S);
	//printf("\n Post Convert Data \n");
	//printf("S Data 0: %f\n", (*S)->data[0]);
	//printf("S Data 1: %f\n", (*S)->data[1]);

	g_ptr_array_foreach(memPtrs, free_memory, NULL);
	g_ptr_array_free(memPtrs, TRUE);
	g_hash_table_destroy(sparseS);
	LIKWID_MARKER_STOP("Optimised Sparsemm Sum - Addition");
	return;
}

void optimised_sparsemm_sum(const COO A, const COO B, const COO C,
		const COO D, const COO E, const COO F,
		COO *O)
{
	COO sum1, sum2;
	//printf("\nA + B + C");
	add_3(A, B, C, &sum1);
//	//printf("Sum1 Data 1: %f\n", sum1->data[0]);
//	printf("Sum1 Data 2: %f\n", sum1->data[1]);
//	printf("Sum1 Data 3: %f\n", sum1->data[2]);
//	printf("Sum1 Data 4: %f\n", sum1->data[3]);
//	printf("Sum1 Data 5: %f\n", sum1->data[4]);
//	printf("Sum1 Data 6: %f\n", sum1->data[5]);
//	printf("Sum1 Data 7: %f\n", sum1->data[6]);
//
	//printf("\nD + E + F");
	add_3(D, E, F, &sum2);
	//printf("Sum2 Data 1: %f\n", sum2->data[0]);
	//printf("Sum2 Data 2: %f\n", sum2->data[1]);
	//printf("Sum2 Data 3: %f\n", sum2->data[2]);
	//printf("Sum2 Data 4: %f\n", sum2->data[3]);
	//printf("Sum2 Data 5: %f\n", sum2->data[4]);
	//printf("Sum2 Data 6: %f\n", sum2->data[5]);
	//printf("Sum2 Data 7: %f\n", sum2->data[6]);
	//printf("\nD + E + F");
	optimised_sparsemm(sum1, sum2, O);
	return;
}
