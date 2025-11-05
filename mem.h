

#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int **mem()
{
    int i;
    int **S;
    S = (int **) malloc(DIM * sizeof(int*));

    for (i=0; i<DIM; i++)
        S[i] = (int*)malloc(DIM * sizeof(int));

    return S;
}

double **memd()
{
    int i;
    double **S = (double **) malloc(DIM * sizeof(double*));

    for (i=0; i<DIM; i++)
        S[i] = (double*)malloc(DIM * sizeof(double));




    return S;
}


void rand_seq(unsigned char*a, int N){
    int i, tmp;
    srand(time(NULL));
    for(i=0; i<N; i++)
    {
        tmp = i%4;

        switch(tmp){
            case 0 : a[i] = 'A'; break;
            case 1 : a[i] = 'G'; break;
            case 2 : a[i] = 'C'; break;
            case 3 : a[i] = 'U'; break;
        }

    }

}

void rna_array_init(double **S, int def, int def2){

    int i,j;

    for(i=0; i<=N+5; i++)
        for(j=0; j<=N+5; j++)
            if(i==j || i==0)
                S[i][j] = def;
            else
                S[i][j] = def2;



}

void rna_array_print(double **S){

    int i,j;

    for(i=0; i<N; i++){
        for(j=0; j<N; j++)
            if(i>j)
                printf("       ");
            else
                printf(" %5.3f ", S[i][j]);
        printf("\n");
    }
    printf("\n");
}
