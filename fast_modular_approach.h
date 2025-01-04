/*
 * @Description: the implementation of fast modular approach
 * @Author: Nianqi Tang
 * @Date: 2025-01-04 19:10:16
 * @LastEditors: Nianqi Tang
 * @LastEditTime: 2025-01-04 22:41:27
 */
#ifndef _KEY_EQUATION_H_
#define _KEY_EQUATION_H_

#include "FFT.h"

struct mat_FMA{
    GFNUM *t11;
    GFNUM *t12;
    GFNUM *t21;
    GFNUM *t22;
    GFNUM *f11;
    GFNUM *f12;
    GFNUM *f21;
    GFNUM *f22;
    int deg_11;
    int deg_12;
    int deg_21;
    int deg_22;
};

void alloc_mat_FMA(struct mat_FMA *mat);

void init_mat_FMA(struct mat_FMA *mat);

void free_mat_FMA(struct mat_FMA *mat);

int fast_modular_approach(GFNUM *d, GFNUM *g, int mu, GFNUM beta, struct mat_FMA *mat, int *r1, int *r2);

void test_fast_modular_approach();



#endif