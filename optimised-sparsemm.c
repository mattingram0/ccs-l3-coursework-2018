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

/* Computes C = A*B.
 * C should be allocated by this routine.
 */
void optimised_sparsemm(const COO A, const COO B, COO *C)
{

	//LIKWID_MARKER_START("OptimisedSMM");
	// A is a pointer to a struct
	//	printf("Variable Name | Memory Address | Value\n");
	//	printf("     A 	      | %p | %p \n", &A, A);
	//	printf("              | %p       | %d \n", A, *((int *)A)); //We must first cast our struct pointer to an int pointer and then deference - hacky 

	//	basic_sparsemm(A, B, C);

	// C is a pointer to a pointer to a struct
	//	printf("Variable Name | Memory Address | Value\n");
	//	printf("     C 	      | %p | %p \n", &C, C);
	//	printf("              | %p | %p \n", C, *C);
	//	printf("              | %p       | %d \n", *C, *((int *)*C)); //We deference our pointer which gives us another pointer (memory address), which we then must cast into an integer pointer before we can then derefernce it correctly and print it out using the %int

	int a, b, nzA, nzB, nzC, m, n; 
	double product, *partial, *c = NULL;
	char num1[40];
	char num2[20];
	GHashTable* sparseC = g_hash_table_new(g_str_hash, g_double_equal);

	nzA = A->NZ;
	nzB = B->NZ;
	m = A->m;
	n = B->n;
	*C = NULL;
	nzC = product = 0;

	alloc_dense(m, n, &c);
	zero_dense(m, n, c);

	for(a = 0; a < nzA; a++){
		for(b = 0; b < nzA; b++){
			if(A->coords[a].j == B->coords[b].i){
				c[B->coords[b].j * m + A->coords[a].i] = c[B->coords[b].j * m + A->coords[a].i] + A->data[a] * B->data[b];
				product = A->data[a] * B->data[b];
				sprintf(num1, "%d", "_", "%d", B->coords[b].j, A->coords[a].i); //Convert x coord to string
				//sprintf(num2, "%d", A->coords[a].i); //Convert y coord to string
				//strcat(num1, num2); //Concatenate together
				printf("%s \n", num1);

				partial = (double *) g_hash_table_lookup(sparseC, num1);
				if(partial == NULL){
					g_hash_table_insert(sparseC, num1, &product);
					nzC++;
				}else{
					*partial += product;
				}
			}
		}
	}
	printf("%d", (int) g_hash_table_size(sparseC));
	convert_dense_to_sparse(c, m, n, C);
	free_dense(&c);
	//LIKWID_MARKER_STOP("OptimisedSMM");
	return;
}

void add_value_to_sparse(gpointer coords, gpointer value, gpointer sp){
	//int index = (int)(strchr(coords, '_') - ((char *) coords)) //Get a pointer to where our _ is in our coord string, then subtract from it the pointer to the start of our coords string, thus giving us the index of _ in our string).
	char * j = (char)* coords;
	char * i = strtok(j, ',');
	sp->coords[NZ].i = i;
	sp->coords[NZ].j = j;
	sp->data[NZ] = val;
	
}

/* Converts our hashmap to COO formulated array */
void convert_hashmap_to_sparse(GHashTable* hash, int m, int n, int NZ, COO C){
	alloc_sparse(m, n, NZ, &sp);
	g_hash_table_foreach(hash, add_value_to_sparse, sp);
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
