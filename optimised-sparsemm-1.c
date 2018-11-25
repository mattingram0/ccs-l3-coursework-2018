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

	int a, b, nzA, nzB, m, n; 
	double *c = NULL;

	nzA = A->NZ;
	nzB = B->NZ;
	m = A->m;
	n = B->n;
	*C = NULL;
	nzC = 0;

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
	//LIKWID_MARKER_STOP("OptimisedSMM");
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
