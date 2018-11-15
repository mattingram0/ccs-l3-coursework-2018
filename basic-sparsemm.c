#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utils.h"


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

/* Compute C = C + A*B in dense, column major, format. */
static void dgemm(int m, int n, int k, const double *a, const double *b, double *c)
{

	LIKWID_MARKER_THREADINIT;
	LIKWID_MARKER_START("dgemm");
    int i, j, p;
    int lda = m;
    int ldb = k;
    int ldc = m;
    for (j = 0; j < n; j++) {
        for (p = 0; p < k; p++) {
            for (i = 0; i < m; i++) {
                c[j*ldc + i] = c[j*ldc + i] + a[p*lda + i] * b[j*ldb + p];
            }
        }
    }
    
	LIKWID_MARKER_STOP("dgemm");
	LIKWID_MARKER_CLOSE;  
}

/* Computes C = A*B by converting A and B to dense column major
 * format, performing the matrix-matrix multiplication, and then
 * converting the result back to sparse.
 * C will be allocated by this routine.
 */
//COO A is a pointer to a COO structure:
//A is a symbol referring to location 1234 holding: 50
//Location 50 holds: COO Structure
//COO *C is a pointer to a pointer to a COO structure
//C is a symbol referring to location 1235 holding: 800
//Location 800 holds: 55
//Location 55 holds: COO Structure
//*C = NULL means set the value at the memory address held at C to NULL
//i.e set the value at memory location 800 to NULL
//i.e C is a now a symbol referring to location 1234, holding 800, and location 800 now does not hold anything, so does not point to a COO struct anymore
void basic_sparsemm(const COO A, const COO B, COO *C)
{
	
    double *a = NULL;
    double *b = NULL;
    double *c = NULL;
    int m, n, k;
    //Converts the sparse COO format to dense format - an m*n two dimensional array
    convert_sparse_to_dense(A, &a);
    convert_sparse_to_dense(B, &b);

    *C = NULL;
    m = A->m;
    k = A->n;
    n = B->n;
    if (k != B->m) {
        fprintf(stderr, "Invalid matrix sizes, got %d x %d and %d x %d\n",
                A->m, A->n, B->m, B->n);
        free(a);
        free(b);
        exit(1);
    }
    alloc_dense(m, n, &c);
    zero_dense(m, n, c);


    dgemm(m, n, k, a, b, c); //Call the function that performs the actual matrix multiplication

    free_dense(&a);
    free_dense(&b);
    convert_dense_to_sparse(c, m, n, C);
    free_dense(&c);
}

/* Computes O = (A + B + C) (D + E + F) by converting to dense column
 * major, performing the matrix matrix multiplication, and converting
 * back to sparse.  This routine allocates O.*/
void basic_sparsemm_sum(const COO A, const COO B, const COO C,
                        const COO D, const COO E, const COO F,
                        COO *O)
{

    double *a = NULL;
    double *b = NULL;
    double *c = NULL;
    double *d = NULL;
    double *e = NULL;
    double *f = NULL;
    double *o = NULL;
    int i, j, m, n, k;

    m = A->m;
    k = A->n;
    n = D->n;
    if (A->m != B->m || A->n != B->n) {
        fprintf(stderr, "A (%d x %d) and B (%d x %d) are not the same shape\n",
                A->m, A->n, B->m, B->n);
        exit(1);
    }
    if (A->m != C->m || A->n != C->n) {
        fprintf(stderr, "A (%d x %d) and C (%d x %d) are not the same shape\n",
                A->m, A->n, C->m, C->n);
        exit(1);
    }
    if (D->m != E->m || D->n != E->n) {
        fprintf(stderr, "D (%d x %d) and E (%d x %d) are not the same shape\n",
                D->m, D->n, E->m, E->n);
        exit(1);
    }
    if (D->m != F->m || D->n != F->n) {
        fprintf(stderr, "D (%d x %d) and F (%d x %d) are not the same shape\n",
                D->m, D->n, F->m, F->n);
        exit(1);
    }

    if (A->n != D->m) {
        fprintf(stderr, "Invalid matrix sizes, got %d x %d and %d x %d\n",
                A->m, A->n, D->m, D->n);
        exit(1);
    }
        
    convert_sparse_to_dense(A, &a);
    convert_sparse_to_dense(B, &b);
    convert_sparse_to_dense(C, &c);
    convert_sparse_to_dense(D, &d);
    convert_sparse_to_dense(E, &e);
    convert_sparse_to_dense(F, &f);

    /* Compute sums */
    for (j = 0; j < k; j++) {
        for (i = 0; i < m; i++) {
            a[j*m + i] += b[j*m + i] + c[j*m + i];
        }
    }
    for (j = 0; j < n; j++) {
        for (i = 0; i < k; i++) {
            d[j*k + i] += e[j*k + i] + f[j*k + i];
        }
    }
    free_dense(&b);
    free_dense(&c);
    free_dense(&e);
    free_dense(&f);
    alloc_dense(m, n, &c);
    zero_dense(m, n, c);
    dgemm(m, n, k, a, d, c);
    free_dense(&a);
    free_dense(&d);
    convert_dense_to_sparse(c, m, n, O);
    free_dense(&c);

}
