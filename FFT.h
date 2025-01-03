/*
 * @Description: This file implements LCH-FFT. All polynomials are represented in the basis bar{X}. More detailed background are referred to "https://doi.org/10.48550/arXiv.2207.11079".
 * @Author: Nianqi Tang
 * @Date: 2025-01-02 22:34:38
 * @LastEditors: Nianqi Tang
 * @LastEditTime: 2025-01-03 22:15:16
 */
#ifndef _FFT_H_
#define _FFT_H_

#include "GF.h"
#include "poly.h"

//FFT, compute F[i] = f(omega_i + beta) for i = 0, 1, ..., (2 ^ mu) - 1.
void FFT_bar_X(GFNUM *f, int mu, GFNUM beta, GFNUM *F);

//IFFT, compute f such that F[i] = f(omega_i + beta) for i = 0, 1, ..., (2 ^ mu) - 1.
void IFFT_bar_X(GFNUM *F, int mu, GFNUM beta, GFNUM *f);

//Given F[0], ..., F[epsilon - 1], compute f such that deg(f) < epsilon, and F[i] = f(omega_i + beta) for i = 0, 1, ..., (2 ^ mu) - 1.
void partial_IFFT(GFNUM *F, GFNUM *f, int epsilon, int mu, GFNUM beta);

//Given F[epsilon], ..., F[2^mu - 1], compute f such that deg(f) < 2 ^ mu - epsilon, and F[i] = f(omega_i + beta) for i = 0, 1, ..., (2 ^ mu) - 1.
void special_IFFT(GFNUM *F, GFNUM *f, int mu_epsilon, int mu, GFNUM beta);

//extended IFFT, compute f such that F[i] = f(omega_i + beta) for i = 0, 1, ..., (2 ^ mu).
void extended_IFFT_bar_X(GFNUM *F, int mu, GFNUM beta, GFNUM *f);

//extended FFT, compute F[i] = f(omega_i + beta) for i = 0, 1, ..., (2 ^ mu).
void extended_FFT_bar_X(GFNUM *f, int mu, GFNUM beta, GFNUM *F);

void test_FFT();

void test_partial_special_FFT();

void test_extended_IFFT();

void test_extended_FFT();

#endif