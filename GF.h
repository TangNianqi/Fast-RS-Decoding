/*
 * @Description: construct a finite field GF(2^m) = GF(2)[x] / primpoly. For example, m = 13, primpoly = 020033 = x^13 + x^4 + x^3 + x + 1.
                 Thus the finite field is GF(2)[x] / (x^13 + x^4 + x^3 + x + 1).
 * @Author: Nianqi Tang
 * @Date: 2025-01-02 00:37:06
 * @LastEditors: Nianqi Tang
 * @LastEditTime: 2025-01-03 22:48:27
 */
#ifndef _GF_H_
#define _GF_H_

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "assert.h"
#include "time.h"

//***********
#define m 8
#define primpoly 0435
//***********

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

//GFSIZE = 2 ^ m, GFZERO represents the addition identity
#define GFSIZE (1 << m)
#define GFZERO 0
typedef int GFNUM;

extern int   pow2vec[GFSIZE];
extern int   vec2pow[GFSIZE];

//the following parameters are used for Lin-Chung-Han FFT(LCH-FFT)
//where {v} is a basis of linear space GF(2^m) with respect to GF(2)
//omega[i] = i_0v_0 + i_1v_1 + ... + i_{m-1}v_{m-1}
//s_beta[i][j] = s_i(omega[j]), where s_i is the subspace polynomial
//p_inv[i] is the inverse of p[i], which is defined in LCH-FFT
//para[i][j]   = s_i(omega[j])/p[j], which is the parameter in LCH_FFT

extern GFNUM v[m];
extern GFNUM omega[GFSIZE];
extern GFNUM s_beta[m][GFSIZE];
extern GFNUM p_inv[GFSIZE];
extern GFNUM para[m][GFSIZE];

//the following parameters are used for counting the computation complexity
extern long long mul_cnt;
extern long long add_cnt;
extern long long inv_cnt;

typedef int STATUS;
#define SUCCESS 1
#define FAIL    0

//initalize the finite field
int init_GF();

//field addition
GFNUM GF_ADD(GFNUM a, GFNUM b);

//field multiplication
GFNUM GF_MUL(GFNUM a, GFNUM b);

//field power, return r = a ^ e.
GFNUM GF_POW(GFNUM a, int e);

//field inverse, if a = 0, return 0.
GFNUM GF_INV(GFNUM a);

//print the vector 'v', whose dimension is equal to 'n'.
void print_vec(GFNUM *v, int n);

#endif