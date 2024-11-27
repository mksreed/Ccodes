#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 1000

int main() {
    srand(0);

    double (*A)[SIZE] = malloc(SIZE * SIZE * sizeof(double));
    double (*B)[SIZE] = malloc(SIZE * SIZE * sizeof(double));
    double *C = malloc(SIZE * 128 * sizeof(double));
    double *D = malloc(SIZE * 128 * sizeof(double));
    double (*E)[SIZE / 4] = malloc((SIZE / 2) * (SIZE / 4) * sizeof(double));
    double (*F)[SIZE / 2] = malloc((SIZE / 2) * (SIZE / 2) * sizeof(double));
    double (*G)[SIZE / 2] = malloc((SIZE / 2) * (SIZE / 2) * sizeof(double));

    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            A[i][j] = (double)rand() / RAND_MAX;
            B[i][j] = (double)rand() / RAND_MAX;
        }
    }

    for (int i = 0; i < SIZE * 128; i++) {
        C[i] = (double)rand() / RAND_MAX;
        D[i] = (double)rand() / RAND_MAX;
    }

    for (int i = 0; i < SIZE / 2; i++) {
        for (int j = 0; j < SIZE / 4; j++) {
            E[i][j] = (double)rand() / RAND_MAX;
        }
    }

    for (int i = 0; i < SIZE / 2; i++) {
        for (int j = 0; j < SIZE / 2; j++) {
            F[i][j] = (double)rand() / RAND_MAX;
        }
    }

    // F = dot(F, F.T)
    for (int i = 0; i < SIZE / 2; i++) {
        for (int j = 0; j < SIZE / 2; j++) {
            double sum = 0.0;
            for (int k = 0; k < SIZE / 2; k++) {
                sum += F[i][k] * F[j][k];
            }
            F[i][j] = sum;
        }
    }

    for (int i = 0; i < SIZE / 2; i++) {
        for (int j = 0; j < SIZE / 2; j++) {
            G[i][j] = (double)rand() / RAND_MAX;
        }
    }

    A[1][1], B[20][20];

    // Matrix multiplication
    int N = 1;
    clock_t t = clock();
    for (int i = 0; i < N; i++) {
        double (*result)[SIZE] = malloc(SIZE * SIZE * sizeof(double));
        for (int j = 0; j < SIZE; j++) {
            for (int k = 0; k < SIZE; k++) {
                result[j][k] = 0;
                for (int l = 0; l < SIZE; l++) {
                    result[j][k] += A[j][l] * B[l][k];
                }
            }
        }
        free(result);
    }
    double delta = (double)(clock() - t) / CLOCKS_PER_SEC;
    printf("Dotted two %dx%d matrices in %0.2f s.\n", SIZE, SIZE, delta / N);

    t = clock();
    for (int i = 0; i < N; i++) {
        double (*result)[SIZE] = malloc(SIZE * SIZE * sizeof(double));
        for (int j = 0; j < SIZE; j++) {
            for (int k = 0; k < SIZE; k++) {
                result[j][k] = 0;
                for (int l = 0; l < SIZE; l++) {
                    result[j][k] += A[j][l] * B[l][k];
                }
            }
        }
        free(result);
    }
    delta = (double)(clock() - t) / CLOCKS_PER_SEC;
    printf("Dis matmul two %dx%d matrices in %0.2f s.\n", SIZE, SIZE, delta / N);

    free(A);
    free(B);
    free(C);
    free(D);
    free(E);
    free(F);
    free(G);

    return 0;
}

