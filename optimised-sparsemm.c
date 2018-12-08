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

void sort_coo(const COO A, COO* S, int compare(const void *, const void *), int** rowList){
	int nz, size, c, k, diff, currRow, smallestRowLength = INT_MAX;
	int	*IA = (int *)malloc((A->m + 1) * sizeof(int));
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

	c = currRow = 0;
	head = unsorted;
	IA[0] = 0;
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
		//printf("[%d, %d, %f]\n", sp->coords[c].i, sp->coords[c].j, sp->data[c]);

		//VECTORIZE??? NOTE - SPLIT INTO SEPARATE FUNCTION, may actually improve performance using SIMD
		diff = sp->coords[c].i - currRow;
		for(k = 0; k < diff; k++){
			IA[currRow + k + 1] = c;	
		}

		if(diff > 0){
			if(IA[currRow + diff] - IA[currRow] < smallestRowLength){
				smallestRowLength = IA[currRow + diff] - IA[currRow];
			}
		}
		currRow += diff;
	}

	IA[A->m] = nz;
	*rowList = IA;
	int z = 0;

	//printf("\n");

	for(z; z < A->m + 1; z++){
		//printf("%d ,", IA[z]);
	}
	//printf("\n");
	printf("Smallest Row Length: %d\n", smallestRowLength);
	free(unsorted);
	*S = sp;
	return;
}

void optimised_sparsemm(COO A, COO B, COO *C)
{
//	LIKWID_MARKER_START("Optimised Sparsemm Sum - Multiplication");
	LIKWID_MARKER_START("Optimised Sparsemm - Whole Function");
	//printf("\n In OPTS\n");
	int a, b, nzA, nzB, nzC = 0, m, n, estimate, newEstimate, currentCol = 0, beginRow = 0, *rowList; //TEST THIS - setting rowList = 0 is poor practice
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
	//printf("\nA + B + C: Column Sorted\n");
	sort_coo(A, &AS, col_sort, &rowList);
	//Sort by row, passing our row list (rowList) integer array pointer so that the sort function also returns our row list
	//printf("\nD + E: Row Sorted\n");
	sort_coo(B, &BS, row_sort, &rowList);

	//NOT VECTORISED - multiple nested loops
	//SLP doesn't divide the vector size. Unknown alignment for acess
	LIKWID_MARKER_START("Multiplication");
	
	//Recently added local variables
	int endRow = 0, l = 0;
	double aVal, *bProducts;
	////printf("%d\n", nzB + 1);

	for(a = 0; a < nzA; a++){
		//printf("a %d\n", a);
		//Move our b back down to the beginning of the row which matches our column

		b = beginRow;
		//printf("Before Corrections: a: %d, b: %d\n",a, b);
		//TEST
		while(BS->coords[b].i != AS->coords[a].j && b < nzB && a < nzA){
			while(BS->coords[b].i < AS->coords[a].j && b < nzB && a < nzA)
				b++;
			while(AS->coords[a].j < BS->coords[b].i && b < nzB && a < nzA)
				a++;
		}

		if(a == nzA || b == nzB)
			break;

		//Skip values in our B matrix until our b row matches our a column
		//NOT VECTORISED - control flow in loop
		//while(b < nzB && BS->coords[b].i < AS->coords[a].j){
//			b++;
//		}
//		//Skip values in our A matrix until our a column matches our b row
//		//NOT VECTORISED - control flow in loop
//		while(a < nzA && AS->coords[a].j < BS->coords[b].i){
//			a++;
//		}
//

		beginRow = b;
		aVal = AS->data[a];
		//printf("After Corrections: a: %d, b: %d\n",a, b);
		//printf("BS->coords[b].i + 1: %d\n", BS->coords[b].i + 1);
		endRow = rowList[BS->coords[b].i + 1];
		bProducts = malloc((endRow - beginRow) * sizeof(double));
		
		//Pad/Vec
		l = 0;

		for(b; b < endRow; b++){
			bProducts[l] = aVal * B->data[b];
			//printf("%f * %f = %f\n", aVal, B->data[b], bProducts[l]);
			l++;
		}

		l = 0;
		b = beginRow;

		//Store in Hash Table - ugly and inefficient
		for(b; b < endRow; b++){
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
				*memPtr = bProducts[l];		
				g_hash_table_insert(sparseC, key, memPtr);
				nzC++;
				memPtr++;
			}else{
				*partial += bProducts[l];
			}
			l++;
		}
		//printf("\n");
		free(bProducts);
		//printf("\n");	
		//If we have run out of memory, allocate more, proportionally to how far we are through our loops
	}
	//printf("Finished Loops");
	
	LIKWID_MARKER_STOP("Multiplication");
	convert_hashmap_to_sparse(sparseC, m, n, nzC, C);
	//LIKWID_MARKER_STOP("Optimised Sparsemm - Whole Function");
//	LIKWID_MARKER_STOP("Optimised Sparsemm Sum - Multiplication");
	g_ptr_array_foreach(memPtrs, free_memory, NULL);
	free_sparse(&AS);
	free_sparse(&BS);
	g_hash_table_destroy(sparseC);
	LIKWID_MARKER_STOP("Optimised Sparsemm - Whole Function");
	return;
}


/* Computes O = (A + B + C) (D + E + F).
 * O should be allocated by this routine.
 */
void add_3(const COO A, const COO B, const COO C, COO *S){

	LIKWID_MARKER_START("Addition");
	int a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, ar, ac, br, bc, cr, cc, dr, dc, er, ec, fr, fc, nzA, nzB, nzC, nzD, nzE, nzF, nzS = 0;
	double errorProp, estimate;
    struct coord *coords, *tempCoords;
    double *data, *tempData;
	COO D, E, F, sp;

	if(A->NZ > B->NZ && A->NZ > C->NZ){
		estimate = A->NZ * 1.1;
	}
	else if (B->NZ > A->NZ && B->NZ > C->NZ){
		estimate = B->NZ * 1.1;
	}
	else{
		estimate = C->NZ * 1.1;
	}	

	//printf("\nInitial Estimation: %d\n", (int)estimate);

	nzA = A->NZ;
	nzB = B->NZ;
	nzC = C->NZ;

    coords = calloc((int) estimate, sizeof(struct coord));
    data = calloc((int) estimate, sizeof(double));
    //printf("\nAddition:");
	
	//Insert values until all of one matrix has been copied/added
	//printf("Begin first loop\n");
	while(a < nzA && b < nzB && c < nzC){
		ar = A->coords[a].i;
		ac = A->coords[a].j;
		br = B->coords[b].i;
		bc = B->coords[b].j;
		cr = C->coords[c].i;
		cc = C->coords[c].j;

		if((ar < br && ar < cr) || (ar == br && ar < cr && ac < bc) || (ar == cr && ar < br && ac < cc) || (ar == br && ar == cr && ac < bc && ac < cc)){
			//printf("a: %d, b: %d, c: %d, a\n", a, b, c);
			coords[nzS].i = ar;
			coords[nzS].j = ac;
			data[nzS] = A->data[a];
			a++;
		}
		else if((br < ar && br < cr) || (br == ar && br < cr && bc < ac) || (br == cr && br < ar && bc < cc) || (ar == br && ar == cr && bc < ac && bc < cc)){
			//printf("a: %d, b: %d, c: %d, b\n", a, b, c);
			coords[nzS].i = br;
			coords[nzS].j = bc;
			data[nzS] = B->data[b];
			b++;
		}
		else if((cr < ar && cr < br) || (cr == ar && cr < br && cc < ac) || (cr == br && cr < ar && cc < bc) || (ar == br && ar == cr && cc < ac && cc < bc)){
			//printf("a: %d, b: %d, c: %d, c\n", a, b, c);
			coords[nzS].i = cr;
			coords[nzS].j = cc;
			data[nzS] = C->data[c];
			c++;
		}
		else if((ar == br && ac == bc ) && (ar < cr || ac < cc)){
			//printf("a: %d, b: %d, c: %d, a + b\n", a, b, c);
			coords[nzS].i = ar;
			coords[nzS].j = ac;
			data[nzS] = A->data[a] + B->data[b];
			a++;
			b++;
		}
		else if((ar == cr && ac == cc) && (ar < br || ac < bc)){
			//printf("a: %d, b: %d, c: %d, a + c\n", a, b, c);
			coords[nzS].i = ar;
			coords[nzS].j = ac;
			data[nzS] = A->data[a] + C->data[c];
			a++;
			c++;
		}
		else if((br == cr && bc == cc) && (br < ar || bc < ac)){
			//printf("a: %d, b: %d, c: %d, b + c\n", a, b, c);
			coords[nzS].i = br;
			coords[nzS].j = bc;
			data[nzS] = B->data[b] + C->data[c];
			b++;
			c++;
		}
		else{
			//printf("a: %d, b: %d, c: %d, a + b + c\n", a, b, c);
			coords[nzS].i = ar;
			coords[nzS].j = ac;
			data[nzS] = A->data[a] + B->data[b] + C->data[c];
			a++;
			b++;
			c++;
		}
		//printf("[%d, %d, %f]\n", coords[nzS].i, coords[nzS].j, data[nzS]);

		nzS++;

		if(nzS == (int)estimate){
			estimate /= (((double)((double)a + (double)b + (double)c))/((double)((double)nzA + (double)nzB + (double)nzC)) - 0.01); 
			//printf("\n1. Updated Estimation: %d\n", (int)estimate);
    		tempCoords = realloc(coords, (int)estimate * sizeof(struct coord));
    		tempData = realloc(data, (int)estimate * sizeof(double));
    		coords = tempCoords;
    		data = tempData;
		}	
	}
	
	if(a == A->NZ){
		D = B;
		d = b;
		E = C;
		e = c;
	}
	else if(b == B->NZ){
		D = A;
		d = a;
		E = C;
		e = c;
	}
	else{
		D = A;
		d = a;
		E = B;
		e = b;
	}

	nzD = D->NZ;
	nzE = E->NZ;

	//Insert values until all of the matrix of values have been inserted/added
	//printf("\nBegin second loop\n");
	while(d < nzD && e < nzE){
		dr = D->coords[d].i;
		dc = D->coords[d].j;
		er = E->coords[e].i;
		ec = E->coords[e].j;
		
		if((dr < er) || (dr == er && dc < ec)){
			coords[nzS].i = dr;
			coords[nzS].j = dc;
			data[nzS] = D->data[d];
			d++;
		}
		else if((er < dr) || (er == dr && ec < dc)){
			coords[nzS].i = er;
			coords[nzS].j = ec;
			data[nzS] = E->data[e];
			e++;
		}
		else{
			coords[nzS].i = dr;
			coords[nzS].j = dc;
			data[nzS] = D->data[d] + E->data[e];
			d++;
			e++;
		}
		
		//printf("[%d, %d, %f]\n", coords[nzS].i, coords[nzS].j, data[nzS]);
		nzS++;

		if(nzS == (int)estimate){
			estimate /= (((double)((double)d + (double)e))/((double)((double)nzD + (double)nzE)) - 0.01); 
			//printf("\n2. Updated Estimation: %d\n", (int)estimate);
    		tempCoords = realloc(coords, (int)estimate * sizeof(struct coord));
    		tempData = realloc(data, (int)estimate * sizeof(double));
    		coords = tempCoords;
    		data = tempData;
    	}
	}

	if(d == D->NZ){
		F = E;
		f = e;
	}
	else{
		F = D;
		f = d;
	}

	nzF = F->NZ;

	//Insert the remaining values from the last matrix
	//printf("\nBegin third loop\n");
	//NOTE - Vectorise definitely because we can calculate the number of non zeros left and check we have enough memory a priori
	while(f < nzF){
		coords[nzS].i = F->coords[f].i;
		coords[nzS].j = F->coords[f].j;
		data[nzS] = F->data[f];
		f++;
		
		//printf("[%d, %d, %f]\n", coords[nzS].i, coords[nzS].j, data[nzS]);

		nzS++;
		if(nzS == (int)estimate){
			estimate /= (((double)f)/((double)nzF) - 0.01); 
			//printf("\n3. Updated Estimation: %d\n", (int)estimate);
    		tempCoords = realloc(coords, (int)estimate * sizeof(struct coord));
    		tempData = realloc(data, (int)estimate * sizeof(double));
    		coords = tempCoords;
    		data = tempData;
		}
	}

	alloc_sparse(A->m, A->n, nzS, &sp);
	sp->m = A->m;
	sp->n = A->n;
	sp->NZ = nzS;
	sp->coords = coords;
	sp->data = data;
	
	*S = sp;
	//printf("\n");
	LIKWID_MARKER_STOP("Addition");
	return;
	
}

void optimised_sparsemm_sum(const COO A, const COO B, const COO C,
		const COO D, const COO E, const COO F,
		COO *O)
{
	LIKWID_MARKER_START("Optimised Sparsemm Sum - Whole Function");
	COO sum1, sum2;
	//printf("\nA + B + C");
//	LIKWID_MARKER_START("Optimised Sparsemm Sum - First Addition");
	add_3(A, B, C, &sum1);
//	LIKWID_MARKER_STOP("Optimised Sparsemm Sum - First Addition");
//	//printf("Sum1 Data 1: %f\n", sum1->data[0]);
//	printf("Sum1 Data 2: %f\n", sum1->data[1]);
//	printf("Sum1 Data 3: %f\n", sum1->data[2]);
//	printf("Sum1 Data 4: %f\n", sum1->data[3]);
//	printf("Sum1 Data 5: %f\n", sum1->data[4]);
//	printf("Sum1 Data 6: %f\n", sum1->data[5]);
//	printf("Sum1 Data 7: %f\n", sum1->data[6]);
//
	//printf("\nD + E + F");
//	LIKWID_MARKER_START("Optimised Sparsemm Sum - Second Addition");
	add_3(D, E, F, &sum2);
//	LIKWID_MARKER_STOP("Optimised Sparsemm Sum - Second Addition");
	//printf("Sum2 Data 1: %f\n", sum2->data[0]);
	//printf("Sum2 Data 2: %f\n", sum2->data[1]);
	//printf("Sum2 Data 3: %f\n", sum2->data[2]);
	//printf("Sum2 Data 4: %f\n", sum2->data[3]);
	//printf("Sum2 Data 5: %f\n", sum2->data[4]);
	//printf("Sum2 Data 6: %f\n", sum2->data[5]);
	//printf("Sum2 Data 7: %f\n", sum2->data[6]);
	//printf("\nD + E + F");
//	LIKWID_MARKER_START("Optimised Sparsemm Sum - Multiplication");
	optimised_sparsemm(sum1, sum2, O);
//	LIKWID_MARKER_STOP("Optimised Sparsemm Sum - Multiplication");
	LIKWID_MARKER_STOP("Optimised Sparsemm Sum - Whole Function");
	return;
}
