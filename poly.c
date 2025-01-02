/*
 * @Description: Polynomial operations. For a function with "_bar_X" in its name, the polynomial is represented with respect to bar{X}.
 * @Author: Nianqi Tang
 * @Date: 2025-01-02 22:53:49
 * @LastEditors: Nianqi Tang
 * @LastEditTime: 2025-01-02 22:57:46
 */
#include "poly.h"

void poly_addition(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *abx, int *deg_abx){
    int i = 0, j = 0;
    int max_deg = (deg_ax > deg_bx) ? deg_ax : deg_bx;

    (*deg_abx) = max_deg;

    for(i = 0 ; i <= max_deg ; ++i){
        abx[i] = GFZERO;

        if(i <= deg_ax){
            abx[i] = GF_ADD(abx[i], ax[i]);
        }

        if(i <= deg_bx){
            abx[i] = GF_ADD(abx[i], bx[i]);
        }
    }

    for(i = 0 ; i <= max_deg ; ++i){
        if(abx[i] != 0){
            *deg_abx = i;
        }
    }
}

int poly_multiplication(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *abx, int *deg_abx){
    int i = 0, j = 0;
    int t_ax = 0, t_bx = 0;
    //determine the degree of ax, bx
    for(i = 0 ; i <= deg_ax ; ++i){
        assert(ax[i] < GFSIZE && ax[i] >= 0);
        if(ax[i] != 0){
            t_ax = i;
        }
    }

    for(i = 0 ; i <= deg_bx ; ++i){
        assert(bx[i] < GFSIZE && bx[i] >= 0);
        if(bx[i] != 0){
            t_bx = i;
        }
    }

    if((t_ax == 0 && ax[0] == 0) || (t_bx == 0 && bx[0] == 0)){
        abx[0] = 0;
        (*deg_abx) = 0;
        return 0;
    }

    deg_ax = t_ax, deg_bx = t_bx;
    (*deg_abx) = deg_ax + deg_bx;
    memset(abx, 0, sizeof(GFNUM) * ((*deg_abx) + 1));
    GFNUM temp = 0;

    for(i = 0 ; i <= deg_ax ; ++i){
        for(j = 0 ; j <= deg_bx ; ++j){
            temp = GF_MUL(ax[i], bx[j]);
            abx[i + j] = GF_ADD(abx[i + j], temp);
        }
    }
    return 0;
}

int poly_division(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *qx, int *deg_qx, GFNUM *rx, int *deg_rx){
    int i = 0, j = 0;

    int t_a = 0, t_b = 0;
    //determine the degree of ax, bx
    for(i = 0 ; i <= deg_ax ; ++i){
        assert(ax[i] < GFSIZE && ax[i] >= 0);
        if(ax[i] != 0){
            t_a = i;
        }
    }

    for(i = 0 ; i <= deg_bx ; ++i){
        assert(bx[i] < GFSIZE && bx[i] >= 0);
        if(bx[i] != 0){
            t_b = i;
        }
    }

    if(t_b == 0 && bx[0] == 0){
        printf("The divisor is 0.\n");
        return DIVIDEZERO;
    }

    GFNUM *tx = (GFNUM *) malloc (sizeof(GFNUM) * (deg_ax + 1));
    memcpy(tx, ax, sizeof(GFNUM) * (deg_ax + 1));
    deg_ax = t_a;
    deg_bx = t_b;
    GFNUM coeff = GFZERO;

    memset(rx, 0, sizeof(GFNUM) * (deg_bx));
    if(deg_ax < deg_bx){
        (*deg_qx) = 0;
        qx[0] = 0;
        (*deg_rx) = deg_ax;
        memcpy(rx, ax, sizeof(GFNUM) * (deg_ax + 1));
    }else{

        (*deg_qx) = deg_ax - deg_bx;
        for(i = deg_ax ; i >= deg_bx ; --i){
            qx[i - deg_bx] = GF_MUL(tx[i], GF_INV(bx[deg_bx]));
            coeff = qx[i - deg_bx];
            for(j = i ; j >= (i - deg_bx) ; --j){
                tx[j] = GF_ADD(tx[j], GF_MUL(coeff, bx[j - i + deg_bx]));
            }
        }

        //
        //printf("deg_bx = %d\n", deg_bx);
        memcpy(rx, tx, sizeof(GFNUM) * (deg_bx));
        if(deg_bx == 0) rx[0] = 0;
        *deg_rx = 0;
        for(i = 0 ; i < deg_bx ; ++i){
            if(rx[i] != 0){
                (*deg_rx) = i;
            }
            //printf("%d, deg_rx = %d\n", i, *deg_rx);
        }
        //printf("%d deg_rx = %d\n", rx[0], *deg_rx);
    }
    free(tx);

    if((*deg_rx) == 0 && rx[0] == 0){
        return DIVISIBLE;
    }
    else{
        return NONDIVISIBLE;
    }
}

int poly_modulo(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *rx, int *deg_rx){
    GFNUM *qx  = (GFNUM *) malloc (sizeof(GFNUM) * (deg_ax + 1));
    int deg_qx = 0;
    int res = poly_division(ax, deg_ax, bx, deg_bx, qx, &deg_qx, rx, deg_rx);
    free(qx);
    return res;
}

void poly_normalization(GFNUM *ax, int deg_ax){
    int i = 0;
    if(deg_ax == 0 && ax[0] == 0) return;
    for(i = deg_ax ; i >= 0 ; --i){
        if(ax[i] != 0){
            deg_ax = i;
            break;
        }
    }
    GFNUM coeff = GF_INV(ax[deg_ax]);
    for(i = 0 ; i <= deg_ax ; ++i){
        ax[i] = GF_MUL(ax[i], coeff);
    }
}

//compute a(x)s(x) + b(x)t(x) = 1. if (a(x), b(x)) = 1
int Euclidean_algorithm(GFNUM *ax, int deg_ax, GFNUM *bx, int deg_bx, GFNUM *sx, int *deg_sx, GFNUM *tx, int *deg_tx){
    int i = 0, j = 0;

    if((deg_ax == 0 && ax[0] == 0) || (deg_bx == 0 && bx[0] == 0)){
        return NONCOPRIME;
    }

    //determine the actual degree
    int t_deg_ax = 0, t_deg_bx = 0;
    for(i = 0 ; i <= deg_ax ; ++i){
        if(ax[i] != 0){
            t_deg_ax = i;
        }
    }
    for(i = 0 ; i <= deg_bx ; ++i){
        if(bx[i] != 0){
            t_deg_bx = i;
        }
    }
    // printf("Euclidean algorithm:\n");
    int max_deg = (t_deg_ax >= t_deg_bx) ? t_deg_ax : t_deg_bx;
    int min_deg = (t_deg_ax <  t_deg_bx) ? t_deg_ax : t_deg_bx;
    GFNUM *t_ax = (GFNUM *) malloc (sizeof(GFNUM) * (max_deg + 1));
    GFNUM *t_bx = (GFNUM *) malloc (sizeof(GFNUM) * (min_deg + 1));

    if(t_deg_ax >= t_deg_bx){
        memcpy(t_ax, ax, sizeof(GFNUM) * (max_deg + 1));
        memcpy(t_bx, bx, sizeof(GFNUM) * (min_deg + 1));
    }
    else{
        memcpy(t_ax, bx, sizeof(GFNUM) * (max_deg + 1));
        memcpy(t_bx, ax, sizeof(GFNUM) * (min_deg + 1));
    }
    
    //initialization
    GFNUM *r_m1 = (GFNUM *) malloc (sizeof(GFNUM) * (max_deg + 1));
    GFNUM *r_0 = (GFNUM *) malloc (sizeof(GFNUM) * (max_deg + 1));
    GFNUM *r_1 = (GFNUM *) malloc (sizeof(GFNUM) * (max_deg + 1));
    GFNUM *s_m1 = (GFNUM *) malloc (sizeof(GFNUM) * (max_deg + 1));
    GFNUM *s_0 = (GFNUM *) malloc (sizeof(GFNUM) * (max_deg + 1));
    GFNUM *s_1 = (GFNUM *) malloc (sizeof(GFNUM) * (max_deg + 1));
    GFNUM *t_m1 = (GFNUM *) malloc (sizeof(GFNUM) * (max_deg + 1));
    GFNUM *t_0 = (GFNUM *) malloc (sizeof(GFNUM) * (max_deg + 1));
    GFNUM *t_1 = (GFNUM *) malloc (sizeof(GFNUM) * (max_deg + 1));
    GFNUM *q = (GFNUM *) malloc (sizeof(GFNUM) * (max_deg + 1));
    GFNUM *temp = (GFNUM *) malloc (sizeof(GFNUM) * (max_deg + 1));

    // printf("max_deg = %d, min_deg = %d\n", max_deg, min_deg);
    memcpy(r_m1, t_ax, sizeof(GFNUM) * (max_deg + 1));
    memcpy(r_0,  t_bx, sizeof(GFNUM) * (min_deg + 1)); 
    memset(s_m1, 0, sizeof(GFNUM) * (max_deg + 1));
    s_m1[0] = 1;
    memset(s_0, 0, sizeof(GFNUM) * (max_deg + 1));
    memset(t_m1, 0, sizeof(GFNUM) * (max_deg + 1));
    memset(t_0, 0, sizeof(GFNUM) * (max_deg + 1));
    t_0[0] = 1;
    int deg_r_m1 = max_deg, deg_r_0 = min_deg, deg_r_1 = 0;
    int deg_s_m1 = 0, deg_s_0 = 0, deg_s_1 = 0;
    int deg_t_m1 = 0, deg_t_0 = 0, deg_t_1 = 0;
    int deg_q    = 0;
    int deg_temp = 0;

    int flag = 0;
    while(1){
        if(deg_r_0 == 0){
            break;
        }

        // print_vec(r_m1, deg_r_m1 + 1);
        // print_vec(r_0,  deg_r_0 + 1);
        poly_division(r_m1, deg_r_m1, r_0, deg_r_0, q, &deg_q, r_1, &deg_r_1);
        // print_vec(q, deg_q + 1);
        // printf("deg_r_1 = %d\n", deg_r_1);
        // print_vec(r_1, deg_r_1 + 1);
        // getchar();

        // printf("deg_q = %d, deg_s_0 = %d\n", deg_q, deg_s_0);
        assert(deg_q + deg_s_0 <= max_deg);
        poly_multiplication(q, deg_q, s_0, deg_s_0, temp, &deg_temp);
        poly_addition(s_m1, deg_s_m1, temp, deg_temp, s_1, &deg_s_1);

        // printf("deg_q = %d, deg_t_0 = %d\n", deg_q, deg_t_0);
        assert(deg_q + deg_t_0 <= max_deg);
        poly_multiplication(q, deg_q, t_0, deg_t_0, temp, &deg_temp);
        poly_addition(t_m1, deg_t_m1, temp, deg_temp, t_1, &deg_t_1);

        memcpy(r_m1, r_0, sizeof(GFNUM) * (max_deg + 1));
        memcpy(r_0 , r_1, sizeof(GFNUM) * (max_deg + 1));
        deg_r_m1 = deg_r_0;
        deg_r_0  = deg_r_1;

        memcpy(s_m1, s_0, sizeof(GFNUM) * (max_deg + 1));
        memcpy(s_0 , s_1, sizeof(GFNUM) * (max_deg + 1));
        deg_s_m1 = deg_s_0;
        deg_s_0  = deg_s_1;

        memcpy(t_m1, t_0, sizeof(GFNUM) * (max_deg + 1));
        memcpy(t_0 , t_1, sizeof(GFNUM) * (max_deg + 1));
        deg_t_m1 = deg_t_0;
        deg_t_0  = deg_t_1;

    }
    // printf("end.\n");
    if(deg_r_0 == 0 && r_0[0] != 0){
        flag = COPRIME;
    }
    else{
        flag = NONCOPRIME;
    }

    GFNUM coeff = GF_INV(r_0[0]);
    // printf("r_0[0] = %d, coeff = %d\n", r_0[0], coeff);
    // print_vec(s_0, deg_s_0 + 1);
    // print_vec(t_0, deg_t_0 + 1);
    for(i = 0 ; i <= deg_s_0 ; ++i){
        s_0[i] = GF_MUL(s_0[i], coeff);
    }

    for(i = 0 ; i <= deg_t_0 ; ++i){
        t_0[i] = GF_MUL(t_0[i], coeff);
    }
    if(t_deg_ax >= t_deg_bx){
        memcpy(sx, s_0, sizeof(GFNUM) * (deg_s_0 + 1));
        memcpy(tx, t_0, sizeof(GFNUM) * (deg_t_0 + 1));
        *deg_sx = deg_s_0;
        *deg_tx = deg_t_0;
    }
    else{
        memcpy(sx, t_0, sizeof(GFNUM) * (deg_t_0 + 1));
        memcpy(tx, s_0, sizeof(GFNUM) * (deg_s_0 + 1));

        *deg_sx = deg_t_0;
        *deg_tx = deg_s_0;
    }

    free(t_ax);
    free(t_bx);
    free(r_m1);
    free(r_0);
    free(r_1);
    free(s_m1);
    free(s_0);
    free(s_1);
    free(t_m1);
    free(t_0);
    free(t_1);
    free(q);
    free(temp);
    return flag;
}

void test_mul_div(){
    int i = 0, j = 0;

    GFNUM ax[2000] = {0};
    GFNUM bx[2000] = {0};
    GFNUM abx[4000] = {0};
    GFNUM qx[4000] = {0};
    GFNUM rx[4000] = {0};
    GFNUM tx[4000] = {0};

    int deg_ax = 0, deg_bx = 0, deg_abx = 0;
    int deg_qx = 0, deg_rx = 0;

    int loops   = 0, maxloops = 50000;

    GFNUM beta = 0, temp = 0;
    int flag = 0;
    int deg_tx = 0;

    for(loops = 0 ; loops < maxloops ; ++loops){
        deg_ax = rand() % 500;
        deg_bx = rand() % 500;

        for(i = 0 ; i <= deg_ax ; ++i){
            ax[i] = rand() % GFSIZE;
        }
        for(i = 0 ; i <= deg_bx ; ++i){
            bx[i] = rand() % GFSIZE;
        }

        //print_vec(ax, deg_ax + 1);
        //print_vec(bx, deg_bx + 1);

        poly_multiplication(ax, deg_ax, bx, deg_bx, abx, &deg_abx);
        //print_vec(abx, deg_abx + 1);
        flag = poly_division(abx, deg_abx, bx, deg_bx, qx, &deg_qx, rx, &deg_rx);
        if(flag == DIVIDEZERO){
            continue;
        }
        // print_vec(qx, deg_qx + 1);
        // print_vec(rx, deg_rx + 1);
        // getchar();
        if(flag != DIVISIBLE){
            print_vec(ax, deg_ax + 1);
            print_vec(bx, deg_bx + 1);
            print_vec(abx, deg_abx + 1);
            print_vec(qx, deg_qx + 1);
            print_vec(rx, deg_rx + 1);
            printf("poly mul or div is wrong!\n");
            getchar();
        }

        // print_vec(ax, deg_ax + 1);
        // print_vec(bx, deg_bx + 1);
        // print_vec(qx, deg_qx + 1);

        flag = poly_division(ax, deg_ax, bx, deg_bx, qx, &deg_qx, rx, &deg_rx);
        if(flag == DIVIDEZERO){
            continue;
        }

        //print_vec(qx, deg_qx + 1);

        poly_multiplication(bx, deg_bx, qx, deg_qx, abx, &deg_abx);
        // printf("mul\n");
        
        // print_vec(abx, deg_abx + 1);
        // print_vec(rx, deg_rx + 1);
        // getchar();

        poly_addition(abx, deg_abx, rx, deg_rx, tx, &deg_tx);
        // printf("add\n");
        // getchar();
        // if(deg_tx != deg_ax){
        //     printf("poly mul or div is wrong!\n");
        //     print_vec(ax, deg_ax + 1);
        //     print_vec(bx, deg_bx + 1);
        //     print_vec(qx, deg_qx + 1);
        //     print_vec(rx, deg_rx + 1);
        //     print_vec(abx, deg_abx + 1);
        //     print_vec(tx, deg_tx + 1);
        //     getchar();
        // }
        for(i = 0 ; i <= deg_tx ; ++i){
            if(tx[i] != ax[i]){
                print_vec(ax, deg_ax + 1);
                print_vec(bx, deg_bx + 1);
                print_vec(abx, deg_abx + 1);
                print_vec(qx, deg_qx + 1);
                print_vec(rx, deg_rx + 1);
                printf("poly mul or div is wrong1!\n");
                getchar();
            }
        }

        printf("%.2f%%\r", (double) loops / maxloops * 100);
    }
}

void test_Euclidean_algorithm(){
    int i = 0, j = 0;

    GFNUM ax[2000] = {0};
    GFNUM bx[2000] = {0};
    GFNUM sx[4000] = {0};
    GFNUM tx[4000] = {0};
    GFNUM rx[4000] = {0};
    GFNUM asx[4000] = {0};
    GFNUM btx[4000] = {0};

    int deg_ax = 0, deg_bx = 0, deg_abx = 0;
    int deg_qx = 0, deg_rx = 0;
    int deg_sx = 0, deg_tx = 0;
    int deg_asx = 0, deg_btx = 0;

    int loops   = 0, maxloops = 50000;

    GFNUM beta = 0, temp = 0;
    int flag = 0;

    for(loops = 0 ; loops < maxloops ; ++loops){
        deg_ax = rand() % 500;
        deg_bx = rand() % 500;

        for(i = 0 ; i <= deg_ax ; ++i){
            ax[i] = rand() % GFSIZE;
        }
        for(i = 0 ; i <= deg_bx ; ++i){
            bx[i] = rand() % GFSIZE;
        }
        // deg_ax = 1;
        // deg_bx = 2;
        // ax[0] = 1, ax[1] = 2;
        // bx[0] = 1, bx[1] = 1; bx[2] = 1;

        //print_vec(ax, deg_ax + 1);
        //print_vec(bx, deg_bx + 1);
        flag = Euclidean_algorithm(ax, deg_ax, bx, deg_bx, sx, &deg_sx, tx, &deg_tx);
        // print_vec(sx, deg_sx + 1);
        // print_vec(tx, deg_tx + 1);
        // printf("deg_sx = %d, deg_tx = %d\n", deg_sx, deg_tx);
        // getchar();
        //a(x)s(x)
        poly_multiplication(ax, deg_ax, sx, deg_sx, asx, &deg_asx);
        poly_multiplication(bx, deg_bx, tx, deg_tx, btx, &deg_btx);
        poly_addition(asx, deg_asx, btx, deg_btx, rx, &deg_rx);
        // print_vec(asx, deg_asx + 1);
        // print_vec(btx, deg_btx + 1);
        // print_vec(rx, deg_rx + 1);
        // getchar();
        if(flag == NONCOPRIME){
            continue;
        }
        else if(deg_rx != 0 || rx[0] != 1){
            printf("Euclidean algorithm is wrong.\n");
        }
        
        printf("%.2f%%\r", (double) loops / maxloops * 100);
    }
}