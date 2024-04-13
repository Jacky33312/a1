#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define BLOCK_SIZE 4
#define LOOP_UNROLLING_FACTOR 4

void gaussian_elimination_revise(int n, double A[][n]) {
    int i, j;
    double temp;
    double epsilon = 1e-10; // Define epsilon for near-singularity checks

    for (i = 0; i < n - 1; i++) {
        int k = i;

        // Find the row with maximum absolute value in the column i
        for (j = i + 1; j < n; j++) {
            if (fabs(A[j][i]) > fabs(A[k][i])) {
                k = j;
            }
        }

        // Check if A is singular or nearly so
        if (fabs(A[k][i]) < epsilon) {
            printf("Matrix is singular or nearly singular\n");
            return;
        }

        // Swap rows i and k
        if (k != i) {
            for (j = 0; j < n; j++) {
                temp = A[i][j];
                A[i][j] = A[k][j];
                A[k][j] = temp;
            }
        }

        // Store multiplier in place of A(j,i)
        for (j = i + 1; j < n; j++) {
            A[j][i] /= A[i][i];
        }

        // Subtract multiple of row A(i,:) to zero out A(j,i)
        for (j = i + 1; j < n; j++) {
            int k_index = n - ((n-k) % 4);
            for (k = i + 1; k < k_index; k += 4) {
                A[j][k] -= A[j][i] * A[i][k];
                A[j][k + 1] -= A[j][i] * A[i][k + 1];
                A[j][k + 2] -= A[j][i] * A[i][k + 2];
                A[j][k + 3] -= A[j][i] * A[i][k + 3];
            }

            for (; k < n; k++) {
                A[j][k] -= A[j][i] * A[i][k];
            }

            A[j][i] = 0;
        }
    }
}
void gaussian_elimination(int n, double A[][n]) {
    int i, j;
    double temp;
    double epsilon = 1e-10; // Define epsilon for near-singularity checks

    for (i = 0; i < n - 1; i++) {
        int k = i;

        // Find the row with maximum absolute value in the column i
        for (j = i + 1; j < n; j++) {
            if (fabs(A[j][i]) > fabs(A[k][i])) {
                k = j;
            }
        }

        // Check if A is singular or nearly so
        if (fabs(A[k][i]) < epsilon) {
            printf("Matrix is singular or nearly singular\n");
            return;
        }

        // Swap rows i and k
        if (k != i) {
            for (j = 0; j < n; j++) {
                temp = A[i][j];
                A[i][j] = A[k][j];
                A[k][j] = temp;
            }
        }

        // Store multiplier in place of A(j,i)
        for (j = i + 1; j < n; j++) {
            A[j][i] /= A[i][i];
        }

        // Subtract multiple of row A(i,:) to zero out A(j,i)
        for (j = i + 1; j < n; j++) {
            for (k = i + 1; k < n; k++) {
                A[j][k] -= A[j][i] * A[i][k];
            }

            A[j][i] = 0;
        }
    }
}


int main() {
    int n;

    printf("Enter the size of the matrix A: ");
    int state = scanf("%d", &n);
    if (state != 1) {
        return 0;
    }

    double (*matrix)[n] = malloc(n * sizeof((*matrix)));
    double (*matrix2)[n] = malloc(n * sizeof((*matrix2)));

    if (matrix == NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }

    srand(time(NULL));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int number = rand() % 100;
            matrix[i][j] = number;
            matrix2[i][j] = number;
        }
    }


    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    gaussian_elimination_revise(n, matrix);
    clock_gettime(CLOCK_MONOTONIC, &end);

    /*for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%0.2f ", matrix[i][j]);
        }
        printf("\n");
    }*/


    double elapsedTime = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("Revise function execution time: %f seconds\n", elapsedTime);

    struct timespec start2, end2;
    clock_gettime(CLOCK_MONOTONIC, &start2);
    gaussian_elimination(n, matrix2);
    clock_gettime(CLOCK_MONOTONIC, &end2);

    /*for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%0.2f ", matrix2[i][j]);
        }
        printf("\n");
    }*/

    double elapsedTime2 = (end2.tv_sec - start2.tv_sec) + (end2.tv_nsec - start2.tv_nsec) / 1e9;

    printf("Function execution time: %f seconds\n", elapsedTime2);

    free(matrix);
    free(matrix2);

    return 0;
}