#ifndef _POLY_H_
#define _POLY_H_

#include "GF.h"

#define DIVISIBLE    1
#define NONDIVISIBLE 0
#define DIVIDEZERO   2
#define NONCOPRIME   3
#define COPRIME      4

void poly_addition(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *abx, int *deg_abx);

int poly_multiplication(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *abx, int *deg_abx);

int poly_division(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *qx, int *deg_qx, GFNUM *rx, int *deg_rx);

int poly_modulo(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *rx, int *deg_rx);

void poly_normalization(GFNUM *ax, int deg_ax);

int Euclidean_algorithm(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *sx, int *deg_sx, GFNUM *tx, int *deg_tx);

void test_mul_div();

void test_Euclidean_algorithm();

#endif