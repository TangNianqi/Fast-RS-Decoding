/*
 * @Description: This file implements LCH-FFT. All polynomials are represented in the basis bar{X}. More detailed background are referred to "https://doi.org/10.48550/arXiv.2207.11079".
 * @Author: Nianqi Tang
 * @Date: 2025-01-02 22:34:38
 * @LastEditors: Nianqi Tang
 * @LastEditTime: 2025-01-02 22:51:18
 */
#ifndef _FFT_H_
#define _FFT_H_

#include "GF.h"
#include "poly.h"

void FFT_bar_X(GFNUM *f, int mu, GFNUM beta, GFNUM *F);

void IFFT_bar_X(GFNUM *F, int mu, GFNUM beta, GFNUM *f);

void partial_IFFT(GFNUM *F, GFNUM *f, int epsilon, int mu, GFNUM beta);

void special_IFFT(GFNUM *F, GFNUM *f, int mu_epsilon, int mu, GFNUM beta);

void extended_IFFT_bar_X(GFNUM *F, int mu, GFNUM beta, GFNUM *f);

void extended_FFT_bar_X(GFNUM *f, int mu, GFNUM beta, GFNUM *F);

void test_FFT();

void test_partial_special_FFT();

void test_extended_IFFT();

void test_extended_FFT();

#endif