#ifndef _UTILS_H
#define _UTILS_H

#include <stdio.h>

//Coordinate Structure
struct coord {
    int i, j;
};

//Create a COO struct, where COO = "Coordinate List" = List of (row, column, value) tuples.
//m = no. row, n = no. columns, NZ = no. non zero. Coords pointer and data pointer
//struct coord *coords is a pointer to an array of coordinate structs
struct _p_COO {
    int m, n, NZ;
    struct coord *coords;
    double *data;
};

//Allows the creation of COO pointers - i.e variables which hold addresses pointing to COO structs,
//defined simply using COO *variable_name
//COO test defines a pointer called test, that will point to the _p_COO structs;
//Recall that -> is equivalent to (*). i.e it is used to access members of a struct pointer, and . is used to access the members of a struct
//Struct test1: test1 is a symbol a memory location on the stack, and all the structs members are located here sequentially
//Struct *test2: test2 is a symbol to a memory location on the stack, but this location holds a memory location, which will usually be a heap memory location
//and the struct itself will be located here. i.e test2 is a pointer to a struct
typedef struct _p_COO *COO;

//COO* in the function declaration means alloc_sparse takes a pointer to a pointer to a pointer to a COO struct. If you declare something of a COO type you get a pointer to a struct, if you then use COO* as an argument, you are telling C that you want a pointer to a pointer to a struct
void alloc_sparse(int, int, int, COO*);
void free_sparse(COO*);
void alloc_dense(int, int, double **);
void free_dense(double **);
void zero_dense(int, int, double *);

void convert_sparse_to_dense(const COO, double **);
void convert_dense_to_sparse(const double *, int, int, COO *);

void read_sparse(const char *, COO *);
void write_sparse(FILE *, COO);
void print_sparse(COO);
void random_matrix(int, int, double, COO *);

#endif
