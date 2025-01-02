#include "FFT.h"

void FFT_bar_X(GFNUM *f, int mu, GFNUM beta, GFNUM *F){
    int i = 0 , j = 0, l = 0;
    
    int size = (1 << mu);
    memcpy(F, f, sizeof(GFNUM) * (size));
    int step = 0, step1 = 0;

    for(i = mu - 1 ; i >= 0 ; --i){
        step  = (1 << i);
        step1 = (1 << (i + 1));
        for(j = 0 ; j < size ; j += step1){
            for(l = j ; l < j + step ; ++l){
                F[l]        = GF_ADD(F[l], GF_MUL(para[i][beta ^ j], F[l + step]));
                F[l + step] = GF_ADD(F[l], F[l + step]);
            }
        }
    }
}

void IFFT_bar_X(GFNUM *F, int mu, GFNUM beta, GFNUM *f){
    int i = 0 , j = 0, l = 0;
    
    int size = (1 << mu);
    memcpy(f, F, sizeof(GFNUM) * (size));
    int step = 0, step1 = 0;

    for(i = 0 ; i < mu ; ++i){
        step  = (1 << i);
        step1 = (1 << (i + 1));
        for(j = 0 ; j < size ; j += step1){
            for(l = j ; l < j + step ; ++l){
                f[l + step] = GF_ADD(f[l], f[l + step]);
                f[l]        = GF_ADD(f[l], GF_MUL(para[i][beta ^ j], f[l + step]));
            }
        }
    }
}

void extended_IFFT_bar_X(GFNUM *F, int mu, GFNUM beta, GFNUM *f){
    int i = 0, j = 0;
    int size = 1 << mu;
    GFNUM *hat_f = (GFNUM *) malloc (sizeof(GFNUM) * (size));
    IFFT_bar_X(F, mu, beta, hat_f);
    GFNUM temp = substition_bar_X(hat_f, size - 1, GF_ADD(beta, omega[size]));
    f[size] = GF_ADD(F[size], temp);
    memcpy(f, hat_f, sizeof(GFNUM) * size);
    f[0] = GF_ADD(f[0], GF_MUL(para[mu][beta], f[size]));
    free(hat_f);
}

void extended_FFT_bar_X(GFNUM *f, int mu, GFNUM beta, GFNUM *F){
    int i = 0, j = 0;
    int size = 1 << mu;
    GFNUM *hat_f = (GFNUM *) malloc (sizeof(GFNUM) * (size));
    memcpy(hat_f, f, sizeof(GFNUM) * size);
    hat_f[0] = GF_ADD(hat_f[0], GF_MUL(f[size], para[mu][beta]));
    // print_vec(hat_f, size);
    FFT_bar_X(hat_f, mu, beta, F);
    // print_vec(F, size);
    free(hat_f);
}

void test_FFT(){
    int i = 0 , j = 0;

    long long loops = 0, maxloops = 50000;
    GFNUM f[GFSIZE] = {0};
    GFNUM F[GFSIZE] = {0};
    GFNUM beta      = 0;

    GFNUM t_f[GFSIZE] = {0};
    GFNUM t_F[GFSIZE] = {0};

    for(loops = 0 ; loops < maxloops ; ++loops){
        for(i = 0 ; i < GFSIZE ; ++i){
            f[i] = rand() % GFSIZE;
        }
        
        beta = rand() % GFSIZE;
        //fft
        FFT_bar_X(f, m, beta, F);

        for(i = 0 ; i < GFSIZE ; ++i){
            t_F[i] = substition_bar_X(f, GFSIZE - 1, GF_ADD(omega[i], beta));
        }

        IFFT_bar_X(F, m, beta, t_f);

        for(i = 0 ; i < GFSIZE ; ++i){
            if(f[i] != t_f[i]){
                printf("FFT/IFFT is wrong!\n");
                getchar();
            }
            if(F[i] != t_F[i]){
                print_vec(f, GFSIZE);
                print_vec(F, GFSIZE);
                print_vec(t_F, GFSIZE);
                printf("FFT/substitution_bar is wrong!\n");
                getchar();
            }
        }
        printf("%.2f%%\r", (double)(loops) / maxloops * 100);
    }
}

void test_extended_IFFT(){
    int i = 0 , j = 0;

    long long loops = 0, maxloops = 50000;
    GFNUM f[GFSIZE] = {0};
    GFNUM F[GFSIZE] = {0};
    GFNUM beta      = 0;
    int   mu        = 0;
    int   size      = 0;
    GFNUM t_f[GFSIZE] = {0};
    GFNUM t_F[GFSIZE] = {0};

    for(loops = 0 ; loops < maxloops ; ++loops){
        mu = rand() % m;
        size = 1 << mu;
        memset(f, 0, sizeof(GFNUM) * GFSIZE);
        for(i = 0 ; i <= size ; ++i){
            f[i] = rand() % GFSIZE;
        }
        
        beta = rand() % GFSIZE;
        //fft
        FFT_bar_X(f, m, beta, F);

        extended_IFFT_bar_X(F, mu, beta, t_f);

        for(i = 0 ; i <= size ; ++i){
            if(f[i] != t_f[i]){
                printf("extended IFFT is wrong!\n");
                getchar();
            }
        }
        printf("%.2f%%\r", (double)(loops) / maxloops * 100);
    }
}

void test_extended_FFT(){
    int i = 0 , j = 0;

    long long loops = 0, maxloops = 50000;
    GFNUM f[GFSIZE] = {0};
    GFNUM F[GFSIZE] = {0};
    GFNUM beta      = 0;
    int   mu        = 0;
    int   size      = 0;
    GFNUM t_f[GFSIZE] = {0};
    GFNUM t_F[GFSIZE] = {0};

    for(loops = 0 ; loops < maxloops ; ++loops){
        mu = rand() % m;
        size = 1 << mu;
        memset(f, 0, sizeof(GFNUM) * GFSIZE);
        for(i = 0 ; i <= size ; ++i){
            f[i] = rand() % GFSIZE;
        }
        
        beta = rand() % GFSIZE;
        //fft
        FFT_bar_X(f, m, beta, F);

        extended_FFT_bar_X(f, mu, beta, t_F);

        for(i = 0 ; i < size ; ++i){
            if(F[i] != t_F[i]){
                printf("extended FFT is wrong!\n");
                getchar();
            }
        }
        printf("%.2f%%\r", (double)(loops) / maxloops * 100);
    }
}

void partial_IFFT(GFNUM *F, GFNUM *f, int epsilon, int mu, GFNUM beta){
    if(mu == 0){
        f[0] = F[0];
        return;
    }
    if(epsilon == 0){
        memset(F, 0, sizeof(GFNUM) * (1 << mu));
        memset(f, 0, sizeof(GFNUM) * (1 << mu));
        return ;
    }
    int i      = 0;
    int size   = 1 << (mu - 1);
    GFNUM temp = GF_ADD(omega[size], beta);
    if(epsilon <= size){
        partial_IFFT(F, f, epsilon, mu - 1, beta);
        FFT_bar_X(f, mu - 1, temp, F + size);
        for(i = 0 ; i < size ; ++i){
            f[i + size] = 0;
        }
    }
    else{
        GFNUM *w  = (GFNUM *) malloc (sizeof(GFNUM) * size);
        GFNUM *w1 = (GFNUM *) malloc (sizeof(GFNUM) * size);
        GFNUM *F1 = (GFNUM *) malloc (sizeof(GFNUM) * size);

        IFFT_bar_X(F, mu - 1, beta, w);
        FFT_bar_X(w, mu - 1, temp, w1);

        for(i = 0 ; i < epsilon - size ; ++i){
            F1[i] = GF_ADD(w1[i], F[i + size]);
        }

        partial_IFFT(F1, f + size, epsilon - size, mu - 1, temp);

        for(i = 0 ; i < size ; ++i){
            f[i] = GF_ADD(w[i], GF_MUL(para[mu - 1][beta], f[i + size]));
            F[i + size] = GF_ADD(F1[i], w1[i]);
        }
        free(w);
        free(w1);
        free(F1);
    }
}

void special_IFFT(GFNUM *F, GFNUM *f, int mu_epsilon, int mu, GFNUM beta){
    // print_vec(F, (1 << mu));
    // print_vec(f, (1 << mu));
    // printf("%d %d\n", mu_epsilon, mu);
    // getchar();
    if(mu_epsilon == 0){
        memset(f, 0, sizeof(GFNUM) * (1 << mu));
        memset(F, 0, sizeof(GFNUM) * (1 << mu));
        return;
    }
    if(mu == 0){
        f[0] = F[0];
        return;
    }
    int i = 0, j = 0;
    int size = (1 << (mu - 1));
    int epsilon = (1 << mu) - mu_epsilon;
    GFNUM temp  = GF_ADD(omega[size], beta);
    GFNUM temp1 = 0;
    if(mu_epsilon <= size){
        special_IFFT(F + size, f, mu_epsilon, mu - 1, temp);
        FFT_bar_X(f, mu - 1, beta, F);
        memset(f + size, 0, sizeof(GFNUM) * size);
    }
    else{
        GFNUM *w  = (GFNUM *) malloc (sizeof(GFNUM) * size);
        GFNUM *w1 = (GFNUM *) malloc (sizeof(GFNUM) * size);
        GFNUM *F0 = (GFNUM *) malloc (sizeof(GFNUM) * size);

        IFFT_bar_X(F + size, mu - 1, temp, w);
        FFT_bar_X(w, mu - 1, beta, w1);
        for(i = epsilon ; i < size ; ++i){
            F0[i] = GF_ADD(w1[i], F[i]);
        }
        special_IFFT(F0, f + size, size - epsilon, mu - 1, beta);

        temp1 = GF_ADD(para[mu - 1][beta] , 1);
        for(i = 0 ; i < size ; ++i){
            f[i] = GF_ADD(w[i], GF_MUL(temp1, f[i + size]));
            F[i] = GF_ADD(F0[i], w1[i]);
        }

        free(w);
        free(w1);
        free(F0);
    }
}

void test_partial_special_FFT(){
    int i = 0 , j = 0 ;
    long long loops = 0, maxloops = 50000;

    GFNUM f[GFSIZE]   = {0};
    GFNUM F[GFSIZE]   = {0};
    GFNUM t_f[GFSIZE] = {0};
    GFNUM t_F[GFSIZE] = {0};

    int epsilon       = 0;
    int mu            = 0;
    GFNUM beta        = 0;

    for(loops = 0 ; loops < maxloops ; ++loops){
        mu = rand() % m + 1;
        epsilon = rand() % (1 << mu);//  + 1;
        beta = rand() % GFSIZE;
        for(i = 0 ; i < epsilon ; ++i){
            f[i] = rand() % GFSIZE;
        }

        for(i = epsilon ; i < (1 << mu) ; ++i){
            f[i] = 0;
        }

        FFT_bar_X(f, mu, beta, F);

        for(i = 0 ; i < epsilon ; ++i){
            t_F[i] = F[i];
        }
        //IFFT_bar_X(t_F, mu, beta, t_f);
        partial_IFFT(t_F, t_f, epsilon, mu, beta);

        for(i = 0 ; i < (1 << mu) ; ++i){
            if(t_f[i] != f[i] || t_F[i] != F[i]){
                printf("partial IFFT is wrong!\n");
            }
        }
        printf("%.2f%%\r", (double)loops / maxloops * 100);
    }


    for(loops = 0 ; loops < maxloops ; ++loops){
        mu = rand() % m + 1;
        epsilon = rand() % (1 << mu);
        beta = rand() % GFSIZE;
        for(i = 0 ; i < (1 << mu) - epsilon ; ++i){
            f[i] = rand() % GFSIZE;
        }

        for(i = (1 << mu) - epsilon ; i < (1 << mu) ; ++i){
            f[i] = 0;
        }
        // printf("FFT:\n");
        //printf("epsilon = %d\n", epsilon);
        FFT_bar_X(f, mu, beta, F);
        //getchar();
        for(i = epsilon ; i < (1 << mu) ; ++i){
            t_F[i] = F[i];
        }

        // printf("special FFT:\n");
        // print_vec(t_F, (1 << mu));
        // print_vec(f, (1 << mu));
        // getchar();
        //IFFT_bar_X(t_F, mu, beta, t_f);
        special_IFFT(t_F, t_f, (1 << mu) - epsilon, mu, beta);


        for(i = 0 ; i < (1 << mu) ; ++i){
            if(t_f[i] != f[i] || t_F[i] != F[i]){
                printf("special IFFT is wrong!\n");
            }
        }
        printf("%.2f%%\r", (double)loops / maxloops * 100);
    }
    
}