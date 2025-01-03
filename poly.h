/*
 * @Description: This file provides the polynomial operations
 * @Author: Nianqi Tang
 * @Date: 2025-01-02 22:53:49
 * @LastEditors: Nianqi Tang
 * @LastEditTime: 2025-01-03 22:09:15
 */
#ifndef _POLY_H_
#define _POLY_H_

#include "GF.h"

#define DIVISIBLE    1
#define NONDIVISIBLE 0
#define DIVIDEZERO   2
#define NONCOPRIME   3
#define COPRIME      4

//the polynomial GFNUM *ax of degree deg_ax is equal to \sum_{i = 0}^{deg_ax}ax[i] * x^i.

//return abx = ax + bx
void poly_addition(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *abx, int *deg_abx);

//return abx = ax * bx
int poly_multiplication(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *abx, int *deg_abx);

//return ax = bx * qx + rx
int poly_division(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *qx, int *deg_qx, GFNUM *rx, int *deg_rx);

//return ax = bx * qx + rx, where qx is omitted.
int poly_modulo(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *rx, int *deg_rx);

//make the leading coefficient to be 1.
void poly_normalization(GFNUM *ax, int deg_ax);

//compute ax * sx + bx * tx = 1.
int Euclidean_algorithm(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *sx, int *deg_sx, GFNUM *tx, int *deg_tx);

//substitute 'beta' into 'poly'
GFNUM substitution(GFNUM *poly, int deg, GFNUM beta);

//substitute 'beta' into 'poly', where 'poly' is represented with respect to \bar{X}.
//here we assume that 'deg < 2 ^ m'.
GFNUM substition_bar_X(GFNUM *poly, int deg, GFNUM beta);

void test_mul_div();

void test_Euclidean_algorithm();

#endif