#include "utils.h"
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

	LIKWID_MARKER_START("OptimisedSMM");
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

	int a, b, nzA, nzB, m, n;
	double *c = NULL;
	GHashTable* sparseC = g_hash_table_new(g_str_hash, g_str_equal);

	nzA = A->NZ;
	nzB = B->NZ;
	m = A->m;
	n = B->n;
	*C = NULL;

	alloc_dense(m, n, &c);
	zero_dense(m, n, c);

	for(a = 0; a < nzA; a++){
		for(b = 0; b < nzA; b++){
			if(A->coords[a].j == B->coords[b].i){
				c[B->coords[b].j * m + A->coords[a].i] = c[B->coords[b].j * m + A->coords[a].i] + A->data[a] * B->data[b];
			}
		}
	}
	convert_dense_to_sparse(c, m, n, C);
	free_dense(&c);
	LIKWID_MARKER_STOP("OptimisedSMM");
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
