#include "GF.h"

int   pow2vec[GFSIZE] = {0};
int   vec2pow[GFSIZE] = {0};
GFNUM v[m] = {0};
GFNUM omega[GFSIZE] = {0};
GFNUM s_beta[m][GFSIZE] = {0};
GFNUM p_inv[GFSIZE] = {0};
GFNUM para[m][GFSIZE] = {0};

long long mul_cnt = 0;
long long add_cnt = 0;
long long inv_cnt = 0;

int init_GF(){
    int i = 0, j = 0;

    //first initialize the field
    int state    = 1;
    pow2vec[0]   = state;
    int mask1    = (1 << m);
    int mask2    = mask1 - 1;
    int feedback = primpoly & mask2;

    for(i = 1 ; i < GFSIZE ; ++i){
        state <<= 1;
        if((state & mask1) != 0){
            state ^= feedback;
        }
        state &= mask2;
        pow2vec[i]     = state;
        vec2pow[state] = i;
    }
    vec2pow[1] = 0;
    vec2pow[0] = -1;

    // for(int i = 0 ; i < GFSIZE ; ++i){
    //     printf("%d %d\n", i, pow2vec[i]);
    // }

    // for(int i = 0 ; i < GFSIZE ; ++i){
    //     printf("%d %d\n", i, vec2pow[i]);
    // }
    
    // getchar();

    //init the LCH-FFT parameters
    for(i = 0 ; i < m ; ++i) v[i] = (1 << i);

    for(i = 0 ; i < GFSIZE ; ++i){
        omega[i] = i;
    }

    for(i = 0 ; i < GFSIZE ; ++i){
        s_beta[0][i] = omega[i];
    }

    GFNUM beta = 0, temp = 0;

    for(i = 1 ; i < m ; ++i){
        for(j = 0 ; j < GFSIZE ; ++j){
            beta = j;
            temp = GF_ADD(beta, v[i - 1]);
            s_beta[i][beta] = GF_MUL(s_beta[i - 1][beta], s_beta[i - 1][temp]);
            //printf("%d, %d, %d \n", i, beta, s_beta[i][beta]);
        }
        //getchar();
    }

    for(i = 0 ; i < m ; ++i){
        for(j = 0 ; j < GFSIZE ; ++j){
            beta = j;
            para[i][beta] = GF_MUL(s_beta[i][beta], GF_INV(s_beta[i][v[i]]));
        }
    }

    for(i = 0 ; i < GFSIZE ; ++i){
        p_inv[i] = 1;
        for(j = 0 ; j < m ; ++j){
            if((i & (1 << j)) != 0){
                p_inv[i] = GF_MUL(p_inv[i], s_beta[j][v[j]]);
            }
        }
        p_inv[i] = GF_INV(p_inv[i]);
    }
}

GFNUM GF_ADD(GFNUM a, GFNUM b){
    assert(a >= GFZERO && a < GFSIZE);
    assert(b >= GFZERO && b < GFSIZE);
    GFNUM r = a ^ b;
    ++add_cnt;
    return r;
}

GFNUM GF_MUL(GFNUM a, GFNUM b){
    assert(a >= GFZERO && a < GFSIZE);
    assert(b >= GFZERO && b < GFSIZE);
    ++mul_cnt;
    if(a == GFZERO || b == GFZERO) return GFZERO;
    int e_a = vec2pow[a];
    int e_b = vec2pow[b];
    int e_r = (e_a + e_b) % (GFSIZE - 1);
    return pow2vec[e_r];
}

GFNUM GF_POW(GFNUM a, int e){
    assert(a >= GFZERO && a < GFSIZE);
    assert(e >= 0);
    if(e == 0) return 1;
    if(a == GFZERO) return GFZERO;
    long long e_a = vec2pow[a];
    long long e_e = e;
    long long e_r = e_a * e_e;
    e_r %= (GFSIZE - 1);
    int e_r_int = e_r;
    return pow2vec[e_r_int];
}

GFNUM GF_INV(GFNUM a){
    assert(0 <= a && a < GFSIZE);
    ++inv_cnt;
    if(a == GFZERO) return GFZERO;
    else{
        int pow = vec2pow[a];
        GFNUM res = pow2vec[GFSIZE - 1 - pow];
        return res;
    }
}

void print_vec(GFNUM *v, int n){
    int i = 0;
    for(i = 0 ; i < n ; ++i){
        printf("%d ", v[i]);
    }
    printf("\n");
}
