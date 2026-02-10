
#include <stdlib.h>
#include <limits.h>
#include <omp.h>
#include <math.h>
#include <cstdio>


#define min(a,b) (((a)<(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))

int **c;
int **ck;
int **w;

int N;
int DIM;

#include "mem.h"
#define VEC

int main(int argc, char *argv[]){



    int num_proc=1;
    int i,j,k,ll,p,q;
    int c0, c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;

    int t1, t2, t3, t4, t5, t6;
    int lb, ub, lbp, ubp, lb2, ub2;
    int lbv, ubv;

    srand(time(NULL));


    if(argc > 1)
        num_proc = atoi(argv[1]);

    int kind=1;

    N = 8;
    DIM = 12;
    if(argc > 2)
        N = atoi(argv[2]);
    DIM = N+10;
    int n = N;
    int check = 1;
    if(argc > 3)
        kind = atoi(argv[3]);


    omp_set_num_threads(num_proc);


    c = mem();
    ck = mem();
    w = mem();

    for(i=0; i<DIM; i++)
        for(j=0; j<DIM; j++){
            ck[i][j] = i+j;
            c[i][j] = ck[i][j];
            w[i][j] = i-j;
        }


    double start = omp_get_wtime();
    if(kind==1 || check){
        printf("serial check\n");
#pragma scop
        for(i=n-1; i>=1; i--){
            if(!(i%100)) printf("%i\n", i);
            for(j=i+1; j<=n; j++)
                for(k=i+1; k<j; k++) {
                    if (i == 3 && j == 12)
                        printf("\033[0;32m[%d %d] [%d %d] : [%d] [%d] [w%d] : %d -> [%d] \n", i,k, k, j, ck[i][k], ck[k][j], w[i][j], ck[i][k]+ck[k][j]+w[i][j], ck[i][j]);
                    ck[i][j] = MIN(ck[i][j], w[i][j] + ck[i][k] + ck[k][j]);
                }
        }
#pragma endscop
    }

    if(kind==2)
    {
        int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
        int lb, ub, lbp, ubp, lb2, ub2;
        int lbv, ubv;
        /* Start of CLooG code */
        if (n >= 3) {
            for (t2=-1;t2<=floord(n-16,16);t2++) {
                lbp=t2+1;
                ubp=min(floord(n,16),floord(16*t2+n+13,16));
#pragma omp parallel for private(lbv,ubv,t5,t6,t7,t8,t9,t10,t4) shared(t2)
                for (t4=lbp;t4<=ubp;t4++) {
                    for (t5=max(max(-n+2,16*t2-16*t4),-16*t4-13);t5<=16*t2-16*t4+15;t5++) {
                        for (t7=max(16*t4,-t5+2);t7<=min(n,16*t4+15);t7++) {
                            for (t9=-t5+1;t9<=t7-1;t9++) {
                                c[-t5][t7] = MIN(c[-t5][t7], w[-t5][t7]+c[-t5][t9]+c[t9][t7]);;
                            }
                        }
                    }
                }
            }
        }
        /* End of CLooG code */
    }


    if(kind == 3)  //traco
    {

        {
            for( c1 = 0; c1 < n + floord(-3 * n - 3, 8); c1 += 1)
#pragma omp parallel for shared(c1) private(c3,c5,c9,c11) schedule(dynamic, 1)
                    for( c3 = max(0, c1 - (n + 6) / 8 + 1); c3 <= min(n / 2 - 1, c1 - (c1 + 6) / 5 + 1); c3 += 1)
                        for( c5 = 0; c5 <= c3 / 128; c5 += 1)
                            for( c7 = max(max(-n + 2 * c3 + 1, -n + 8 * c1 - 8 * c3 + 1), -n + c3 + 128 * c5 + 2); c7 <= min(-1, -n + 8 * c1 - 8 * c3 + 8); c7 += 1) {
                                if (n + 8 * c3 + c7 >= 8 * c1 + 2) {
                                    for( c11 = 256 * c5 - c7 + 1; c11 <= min(2 * c3 - c7, 256 * c5 - c7 + 256); c11 += 1)
                                        c[(-c7)][(2*c3-c7+1)] = MIN(c[(-c7)][(2*c3-c7+1)], w[(-c7)][(2*c3-c7+1)]+c[(-c7)][c11]+c[c11][(2*c3-c7+1)]);
                                    if (128 * c5 + 128 >= c3 && n + c7 >= 2 * c3 + 2) {
                                        if (c3 >= 128 * c5 + 1)
                                            for( c11 = -c7 + 1; c11 <= 256 * c5 - c7; c11 += 1)
                                                c[(-c7)][(2*c3-c7+2)] = MIN(c[(-c7)][(2*c3-c7+2)], w[(-c7)][(2*c3-c7+2)]+c[(-c7)][c11]+c[c11][(2*c3-c7+2)]);
                                        for( c11 = 256 * c5 - c7 + 1; c11 <= min(2 * c3 - c7 + 1, 256 * c5 - c7 + 256); c11 += 1)
                                            c[(-c7)][(2*c3-c7+2)] = MIN(c[(-c7)][(2*c3-c7+2)], w[(-c7)][(2*c3-c7+2)]+c[(-c7)][c11]+c[c11][(2*c3-c7+2)]);
                                    }
                                } else {
                                    for( c9 = max(n - 8 * c1 + 10 * c3, n - 8 * c1 + 8 * c3 + 256 * c5 + 1); c9 <= min(n, n - 8 * c1 + 10 * c3 + 1); c9 += 1)
                                        for( c11 = n - 8 * c1 + 8 * c3 + 256 * c5; c11 <= min(n - 8 * c1 + 8 * c3 + 256 * c5 + 255, c9 - 1); c11 += 1)
                                            c[(n-8*c1+8*c3-1)][c9] = MIN(c[(n-8*c1+8*c3-1)][c9], w[(n-8*c1+8*c3-1)][c9]+c[(n-8*c1+8*c3-1)][c11]+c[c11][c9]);
                                }
                            }
            if ((n - 2) % 8 == 0)
                for( c5 = 0; c5 <= floord(n - 10, 256); c5 += 1)
                    for( c11 = 256 * c5 + 2; c11 <= min(n - 1, 256 * c5 + 257); c11 += 1)
                        c[1][n] = MIN(c[1][n], w[1][n]+c[1][c11]+c[c11][n]);
        }


    }

    if(kind == 4)
    {
        for( c0 = 0; c0 <= floord(n - 2, 8); c0 += 1)
#pragma omp parallel for shared(c0) private(c1,c3,c4,c6,c10) schedule(dynamic, 1)
                for( c1 = (c0 + 1) / 2; c1 <= min(c0, (n - 2) / 16); c1 += 1)
                    for( c3 = max(2, 16 * c0 - 16 * c1 + 1); c3 <= min(n - 1, 16 * c0 - 16 * c1 + 16); c3 += 1)
                        for( c4 = max(0, -c1 + (n + 1) / 16 - 1); c4 <= min((n - 1) / 16, -c1 + (n + c3 - 2) / 16); c4 += 1)
                            for( c6 = max(max(-n + 16 * c1 + 1, -n + c3), -16 * c4 - 14); c6 <= min(min(-1, -n + 16 * c1 + 16), c3 - 16 * c4 - 1); c6 += 1)
                                for( c10 = max(16 * c4, -c6 + 1); c10 <= min(16 * c4 + 15, c3 - c6 - 1); c10 += 1)
                                    c[(-c6)][(c3-c6)] = MIN(c[(-c6)][(c3-c6)], w[(-c6)][(c3-c6)]+c[(-c6)][c10]+c[c10][(c3-c6)]);
    }


    // dapt
    if(kind == 5)
    {
        printf("dapt\n");
        int h0, i0, i1, i2, w0;
        for (int w0 = -1; w0 < floord(n, 16); w0 += 1) {
#pragma omp parallel for schedule(dynamic, 1) shared(w0) private(h0, i0, i1, i2)
            for (int h0 = max(w0 - (n + 16) / 16 + 1, -((n + 13) / 16)); h0 < 0; h0 += 1) {
                for (int i0 = max(max(-n + 2, -16 * w0 + 16 * h0 - 13), 16 * h0); i0 <= 16 * h0 + 15; i0 += 1) {
                    for (int i1 = max(16 * w0 - 16 * h0, -i0 + 2); i1 <= min(n, 16 * w0 - 16 * h0 + 15); i1 += 1) {
                        for (int i2 = -i0 + 1; i2 < i1; i2 += 1) {
                            ck[-i0][i1] = MIN(ck[-i0][i1], (w[-i0][i1] + ck[-i0][i2]) + ck[i2][i1]);
                        }
                    }
                }
            }
        }
    }

if(kind == 7)
{
    /* TC Optimizing Compiler 0.5.1 */
    /* ./tc ../examples/npdp/knt.scop.c --merge-tiling --free-scheduling --omp-for-codegen --iterative-tc --debug --align -b 16 */
    
    #pragma scop
    if (n >= 17) {
      for (register int k = 0; k <= (n - 1) / 16; k += 1) {
        #pragma omp parallel for
        for (register int ii0 = max(k - (n + 15) / 16, -((n + 12) / 16)); ii0 < 0; ii0 += 1) {
          if (n + 3 * ii0 >= 14) {
            for (register int i0 = max(max(-n + 3, 16 * ii0), -16 * k + 16 * ii0 + 3); i0 <= 16 * ii0 + 15; i0 += 1) {
              for (register int i1 = max(16 * k - 16 * ii0 - 16, -i0 + 2); i1 < min(n, 16 * k - 16 * ii0); i1 += 1) {
                for (register int i2 = -i0 + 1; i2 < i1; i2 += 1) {
                  ck[-i0][i1] = MIN(ck[-i0][i1], (w[-i0][i1] + ck[-i0][i2]) + ck[i2][i1]);
                }
              }
            }
          }
        }
      }
    }
    #pragma endscop


}


    
    if(kind == 6)  // ?
    {
        int bb = 32;
        for (int c0 = 0; c0 <= floord(N, bb); c0 += 1){
        #pragma omp parallel for
            for (int c1 = c0; c1 <= min(N / bb, (N + 2 * c0 - 2) / bb); c1 += 1) {
                int C[bb][bb];
                for (int i = 0; i < bb; i++)
                    for (int j = 0; j < bb; j++)
                        C[i][j] = INT_MAX;
                for (int c3 = -c0 + c1 + 1; c3 < c1; c3 += 1)
                    for (int c4 = max(1, -bb * c0 + bb * c1); c4 <= -bb * c0 + bb * c1 + bb-1; c4 += 1)
                        for (int c5 = bb * c1; c5 <= min(N, bb * c1 + bb-1); c5 += 1)
                            for (int c6 = bb * c3; c6 <= bb * c3 + bb-1; c6 += 1)
                           {
                               #ifdef VEC
                                C[c4 % bb][c5 % bb] = MIN(C[c4 % bb][c5 % bb], w[c4][c5] + c[c4][c6] + c[c6][c5]);  // c[c4][c5]
                               #else
                                c[c4][c5] = MIN(c[c4][c5], w[c4][c5] + c[c4][c6] + c[c6][c5]);
                               #endif
                            }
                #ifdef VEC
                if(c0 > 1) {
                    c[(-c0 + c1) * bb + bb-1][c1 * bb] = C[bb-1][0]; // nie ma problematycznych
                }
                #endif
                for (int c4 = max(2, bb * c0 - (bb-2)); c4 <= min(min(min(N - 1, bb * c0 + bb-1), N + bb * c0 - bb * c1), bb * c1 + bb - 2); c4 += 1) {
                    if (c0 >= 1) {
                        for (int c5 = max(max(bb * c1, -bb * c0 + bb * c1 + c4), c4 + 1); c5 <= min(min(N, bb * c1 + bb-1), -bb * c0 + bb * c1 + c4 + bb-1); c5 += 1) {
                            for (int c6 = -c4 + c5 + 1; c6 <= -bb * c0 + bb * c1 + bb-1; c6 += 1) {
                                  c[-c4 + c5][c5] = MIN(c[-c4 + c5][c5], w[-c4 + c5][c5] + c[-c4 + c5][c6] + c[c6][c5]);
                            }
                            #ifdef VEC
                            c[-c4 + c5][c5] = MIN(C[(-c4 + c5) % bb][c5 % bb], c[-c4 + c5][c5]);
                            #endif
                            for (int c6 = bb * c1; c6 < c5; c6 += 1) {
                                 c[-c4 + c5][c5] = MIN(c[-c4 + c5][c5], w[-c4 + c5][c5] + c[-c4 + c5][c6] + c[c6][c5]);
                            }
                        }
                    } else {
                        for (int c5 = max(bb * c1 + c4, c4 + 1); c5 <= min(N, bb * c1 + bb-1); c5 += 1)
                            for (int c6 = -c4 + c5 + 1; c6 < c5; c6 += 1)
                                c[-c4 + c5][c5] = MIN(c[-c4 + c5][c5], w[-c4 + c5][c5] + c[-c4 + c5][c6] + c[c6][c5]);
                    }
                }
            }



}


    }

    printf("\033[0m\n");


    double stop = omp_get_wtime();
    printf("%.4f\n",stop - start);


    if(check){
        for(i=0; i<DIM; i++)
            for(j=0; j<DIM; j++){
                if(c[i][j] != ck[i][j]){
                    if(i==3 && j==12)
                    printf("Bad! %d %d %i %i\n", i,j, c[i][j], ck[i][j]);
                    //exit(1);

                }
            }
    }


    return 0;
}
