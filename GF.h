#ifndef _GF_H_
#define _GF_H_

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "assert.h"
#include "time.h"

#define GFZERO 0
#define m 13
#define GFSIZE (1 << m)
#define primpoly 020033

#define VNAME(x) #x
typedef int GFNUM;

extern int   pow2vec[GFSIZE];
extern int   vec2pow[GFSIZE];
extern GFNUM v[m];
extern GFNUM omega[GFSIZE];
extern GFNUM s_beta[m][GFSIZE];
extern GFNUM p_inv[GFSIZE];
extern GFNUM para[m][GFSIZE];

extern long long mul_cnt;
extern long long add_cnt;
extern long long inv_cnt;
//extern int

/*
m       poly
m = 2   7
m = 3   013
m = 4   023
m = 5   045
m = 6   0103
m = 7   0211
m = 8   0435
m = 9   01021
m = 10  02011
m = 11  04005
m = 12  010123  
m = 13  020033
m = 14  042103
m = 15  0100003
m = 16  empty
*/

typedef int STATUS;
#define SUCCESS 1
#define FAIL    0

int init_GF();

GFNUM GF_ADD(GFNUM a, GFNUM b);

GFNUM GF_MUL(GFNUM a, GFNUM b);

GFNUM GF_POW(GFNUM a, int e);

GFNUM GF_INV(GFNUM a);

GFNUM substitution(GFNUM *poly, int deg, GFNUM beta);

GFNUM substition_bar_X(GFNUM *poly, int deg, GFNUM beta);

void gen_poly(GFNUM *roots, int num, GFNUM *poly);

void print_vec(GFNUM *v, int n);

void print_mat(GFNUM **mat, int r, int n);

void alloc_mat(GFNUM ***mat, int row, int col);

void init_mat(GFNUM **mat, int row, int col);

void free_mat(GFNUM **mat, int row);

void linear_transform(GFNUM *vec, GFNUM **mat, int row, int col, GFNUM *res);

void linear_transform_binary(GFNUM *vec, GFNUM **mat, int row, GFNUM *res);

void matrix_multiplication(GFNUM **mat1, int row1, int col1, GFNUM **mat2, int row2, int col2, GFNUM **mat3);

#endif