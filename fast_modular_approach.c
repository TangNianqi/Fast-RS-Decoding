#include "fast_modular_approach.h"

void alloc_mat_FMA(struct mat_FMA *mat){
    mat -> t11 = (GFNUM *) malloc (sizeof(GFNUM) * GFSIZE);
    mat -> t12 = (GFNUM *) malloc (sizeof(GFNUM) * GFSIZE);
    mat -> t21 = (GFNUM *) malloc (sizeof(GFNUM) * GFSIZE);
    mat -> t22 = (GFNUM *) malloc (sizeof(GFNUM) * GFSIZE);
    mat -> f11 = (GFNUM *) malloc (sizeof(GFNUM) * GFSIZE);
    mat -> f12 = (GFNUM *) malloc (sizeof(GFNUM) * GFSIZE);
    mat -> f21 = (GFNUM *) malloc (sizeof(GFNUM) * GFSIZE);
    mat -> f22 = (GFNUM *) malloc (sizeof(GFNUM) * GFSIZE);
}

void init_mat_FMA(struct mat_FMA *mat){
    memset(mat -> t11, 0, sizeof(GFNUM) * GFSIZE);
    memset(mat -> t12, 0, sizeof(GFNUM) * GFSIZE);
    memset(mat -> t21, 0, sizeof(GFNUM) * GFSIZE);
    memset(mat -> t22, 0, sizeof(GFNUM) * GFSIZE);
    memset(mat -> f11, 0, sizeof(GFNUM) * GFSIZE);
    memset(mat -> f12, 0, sizeof(GFNUM) * GFSIZE);
    memset(mat -> f21, 0, sizeof(GFNUM) * GFSIZE);
    memset(mat -> f22, 0, sizeof(GFNUM) * GFSIZE);
    mat -> deg_11 = 0;
    mat -> deg_12 = 0;
    mat -> deg_21 = 0;
    mat -> deg_22 = 0;
}

void free_mat_FMA(struct mat_FMA *mat){
    free(mat -> t11);
    free(mat -> t12);
    free(mat -> t21);
    free(mat -> t22);
    free(mat -> f11);
    free(mat -> f12);
    free(mat -> f21);
    free(mat -> f22);
}

void print_mat_FMA(struct mat_FMA *mat, int mu, GFNUM beta){
    printf("t11:");
    print_vec(mat -> t11, mat -> deg_11 + 1);
    printf("t12:");
    print_vec(mat -> t12, mat -> deg_12 + 1);
    printf("t21:");
    print_vec(mat -> t21, mat -> deg_21 + 1);
    printf("t22:");
    print_vec(mat -> t22, mat -> deg_22 + 1);

    printf("f11:");
    int i = 0;
    int size = 1 << mu;
    GFNUM pos = 0;
    for(i = 0 ; i <= size ; ++i){
        pos = GF_ADD(omega[i], beta);
        printf("%d ", mat -> f11[pos]);
    }
    printf("\n");
    printf("f12:");
    for(i = 0 ; i <= size ; ++i){
        pos = GF_ADD(omega[i], beta);
        printf("%d ", mat -> f12[pos]);
    }
    printf("\n");
    printf("f21:");
    for(i = 0 ; i <= size ; ++i){
        pos = GF_ADD(omega[i], beta);
        printf("%d ", mat -> f21[pos]);
    }
    printf("\n");
    printf("f22:");
    for(i = 0 ; i <= size ; ++i){
        pos = GF_ADD(omega[i], beta);
        printf("%d ", mat -> f22[pos]);
    }
    printf("\n");
}

void evaluate_mat_time(struct mat_FMA *mat, GFNUM beta, int mu){
    extended_FFT_bar_X(mat -> t11, mu, beta, mat -> f11 + beta);
    extended_FFT_bar_X(mat -> t12, mu, beta, mat -> f12 + beta);
    extended_FFT_bar_X(mat -> t21, mu, beta, mat -> f21 + beta);
    extended_FFT_bar_X(mat -> t22, mu, beta, mat -> f22 + beta);
    // printf("beta = %d, mu = %d\n", beta, mu);
    // print_vec(mat -> f21 + beta, 2);
}

void extended_IFFT_bar_X_FMA(GFNUM *F, int mu, GFNUM beta, GFNUM *f){
    int i = 0, j = 0;
    int size = 1 << mu;
    GFNUM *hat_f = (GFNUM *) malloc (sizeof(GFNUM) * (size));
    IFFT_bar_X(F + beta, mu, beta, hat_f);
    GFNUM pos = GF_ADD(beta, omega[size]);
    GFNUM temp = substition_bar_X(hat_f, size - 1, pos);
    f[size] = GF_ADD(F[pos], temp);
    memcpy(f, hat_f, sizeof(GFNUM) * size);
    f[0] = GF_ADD(f[0], GF_MUL(para[mu][beta], f[size]));
    free(hat_f);
}

int fast_modular_approach(GFNUM *d, GFNUM *g, int mu, GFNUM beta, struct mat_FMA *mat, int *r1, int *r2){
    int i = 0, j = 0;
    // printf("beta = %d, mu = %d\n", beta, mu);
    // getchar();
    init_mat_FMA(mat);
    int B = 0;
    if(mu == 0){
        GFNUM ext = GF_ADD(beta, omega[1]);

        B = (g[0] == GFZERO) || ((d[0] != GFZERO) && (*r1 < *r2)) ? 1 : 0;

        mat -> t11[0]    = g[0];
        mat -> f11[beta] = g[0];
        mat -> f11[ext]  = g[0];
        mat -> deg_11    = 0;

        mat -> t12[0]    = d[0];
        mat -> f12[beta] = d[0];
        mat -> f12[ext]  = d[0];
        mat -> deg_12    = 0;

        int temp_r1 = *r1, temp_r2 = *r2;
        if(B == 1){
            mat -> t21[0]    = beta;
            mat -> t21[1]    = omega[1];
            mat -> f21[beta] = 0;
            mat -> f21[ext]  = omega[1];
            mat -> deg_21    = 1;

            mat -> t22[0]    = 0;
            mat -> f22[beta] = 0;
            mat -> f22[ext]  = 0;
            mat -> deg_22    = 0;

            *r1 = temp_r2;
            *r2 = temp_r1 + 2;
        }
        else{
            mat -> t21[0]    = 0;
            mat -> f21[beta] = 0;
            mat -> f21[ext]  = 0;
            mat -> deg_21    = 0;

            mat -> t22[0]    = beta;
            mat -> t22[1]    = omega[1];
            mat -> f22[beta] = 0;
            mat -> f22[ext]  = omega[1];
            mat -> deg_22 = 1;

            *r1 = temp_r1;
            *r2 = temp_r2 + 2;
        }
        return 0;
    }
    else{
        struct mat_FMA mat1;
        alloc_mat_FMA(&mat1);
        init_mat_FMA(&mat1);
        int size = (1 << (mu - 1));
        GFNUM new_beta = GF_ADD(beta, omega[size]);
        fast_modular_approach(d, g, mu - 1, beta, &mat1, r1, r2);
        evaluate_mat_time(&mat1, new_beta, mu - 1);
        GFNUM temp = GF_ADD(beta, omega[1 << mu]);
        mat1.f11[temp] = substition_bar_X(mat1.t11, size, temp);
        mat1.f12[temp] = substition_bar_X(mat1.t12, size, temp);
        mat1.f21[temp] = substition_bar_X(mat1.t21, size, temp);
        mat1.f22[temp] = substition_bar_X(mat1.t22, size, temp);
        // print_mat_FMA(&mat1, mu, beta);
        // getchar();

        GFNUM pos  = 0;
        GFNUM *d1 = (GFNUM *) malloc (sizeof(GFNUM) * size);
        GFNUM *g1 = (GFNUM *) malloc (sizeof(GFNUM) * size);
        
        for(i = 0 ; i < size ; ++i){
            pos = GF_ADD(omega[i], new_beta);
            d1[i] = GF_ADD(GF_MUL(mat1.f11[pos], d[i + size]), GF_MUL(mat1.f12[pos], g[i + size]));
            g1[i] = GF_ADD(GF_MUL(mat1.f21[pos], d[i + size]), GF_MUL(mat1.f22[pos], g[i + size]));
        }
        // printf("d1, g1:\n");
        // print_vec(d1, size);
        // print_vec(g1, size);
        struct mat_FMA mat2;
        alloc_mat_FMA(&mat2);
        init_mat_FMA(&mat2);
        fast_modular_approach(d1, g1, mu - 1, new_beta, &mat2, r1, r2);
        evaluate_mat_time(&mat2, beta, mu - 1);
        temp = GF_ADD(beta, omega[1 << mu]);
        mat2.f11[temp] = substition_bar_X(mat2.t11, size, temp);
        mat2.f12[temp] = substition_bar_X(mat2.t12, size, temp);
        mat2.f21[temp] = substition_bar_X(mat2.t21, size, temp);
        mat2.f22[temp] = substition_bar_X(mat2.t22, size, temp);
        // print_mat_FMA(&mat2, mu, beta);
        // getchar();

        free(d1);
        free(g1);

        for(i = 0 ; i <= (1 << mu) ; ++i){
            pos = GF_ADD(omega[i], beta);
            mat -> f11[pos] = GF_ADD(GF_MUL(mat2.f11[pos], mat1.f11[pos]), GF_MUL(mat2.f12[pos], mat1.f21[pos]));
            mat -> f12[pos] = GF_ADD(GF_MUL(mat2.f11[pos], mat1.f12[pos]), GF_MUL(mat2.f12[pos], mat1.f22[pos]));
            mat -> f21[pos] = GF_ADD(GF_MUL(mat2.f21[pos], mat1.f11[pos]), GF_MUL(mat2.f22[pos], mat1.f21[pos]));
            mat -> f22[pos] = GF_ADD(GF_MUL(mat2.f21[pos], mat1.f12[pos]), GF_MUL(mat2.f22[pos], mat1.f22[pos]));
        }
        extended_IFFT_bar_X_FMA(mat -> f11, mu, beta, mat -> t11);
        extended_IFFT_bar_X_FMA(mat -> f12, mu, beta, mat -> t12);
        extended_IFFT_bar_X_FMA(mat -> f21, mu, beta, mat -> t21);
        extended_IFFT_bar_X_FMA(mat -> f22, mu, beta, mat -> t22);
        mat -> deg_11 = mat -> deg_12 = mat -> deg_21 = mat -> deg_22 = 0;
        for(i = 0 ; i <= (1 << mu) ; ++i){
            if(mat -> t11[i] != 0){
                mat -> deg_11 = i;
            }
            if(mat -> t12[i] != 0){
                mat -> deg_12 = i;
            }
            if(mat -> t21[i] != 0){
                mat -> deg_21 = i;
            }
            if(mat -> t22[i] != 0){
                mat -> deg_22 = i;
            }
        }
        free_mat_FMA(&mat1);
        free_mat_FMA(&mat2);
        // print_mat_FMA(mat, mu, beta);
        // getchar();
        return 0;
    }
}

void test_fast_modular_approach(){
    int i = 0, j = 0;
    int loops   = 0, maxloops = 50000;

    GFNUM d[GFSIZE] = {0};
    GFNUM g[GFSIZE] = {0};
    GFNUM wx[GFSIZE] = {0};
    GFNUM nx[GFSIZE] = {0};

    int mu = 0, size = 0;
    GFNUM beta = 0, pos = 0;
    int r1 = 0, r2 = 1;
    struct mat_FMA mat;
    alloc_mat_FMA(&mat);
    for(loops = 0 ; loops < maxloops ; ++loops){
        mu = 7;
        beta = 0;

        size = 1 << mu;
        for(i = 0 ; i < size ; ++i){
            d[i] = rand() % GFSIZE;
            g[i] = rand() % GFSIZE;
        }
        r1 = 0;
        r2 = 1;
        // print_vec(d, size);
        // print_vec(g, size);
        fast_modular_approach(d, g, mu, beta, &mat, &r1, &r2);
        for(i = 0 ; i < size ; ++i){
            pos = GF_ADD(beta, omega[i]);
            if(GF_ADD(GF_MUL(d[i], mat.f11[pos]), GF_MUL(g[i], mat.f12[pos])) != GFZERO){
                printf("FMA is wrong.\n");
            }

            if(GF_ADD(GF_MUL(d[i], mat.f21[pos]), GF_MUL(g[i], mat.f22[pos])) != GFZERO){
                printf("FMA is wrong.\n");
            }
        }
        printf("%.2f%%\r", (double) loops / maxloops * 100);
    }
    free_mat_FMA(&mat);
}