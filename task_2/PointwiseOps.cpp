#include "PointwiseOps.h"

void Copy(const float (&x)[XDIM][YDIM][ZDIM], float (&y)[XDIM][YDIM][ZDIM]) {
#pragma omp parallel for
    for (int i = 1; i < XDIM - 1; i++)
        for (int j = 1; j < YDIM - 1; j++)
            for (int k = 1; k < ZDIM - 1; k++) y[i][j][k] = x[i][j][k];
}

void Saxpy(const float (&x)[XDIM][YDIM][ZDIM],
           const float (&y)[XDIM][YDIM][ZDIM],
           float (&z)[XDIM][YDIM][ZDIM],
           const float scale) {
    // Should we use OpenMP parallel for here?
    for (int i = 1; i < XDIM - 1; i++)
        for (int j = 1; j < YDIM - 1; j++)
            for (int k = 1; k < ZDIM - 1; k++) z[i][j][k] = x[i][j][k] * scale + y[i][j][k];
}

void MergedSaxpy(float (&x)[XDIM][YDIM][ZDIM],
                 const float (&y)[XDIM][YDIM][ZDIM],
                 const float (&z)[XDIM][YDIM][ZDIM],
                 float (&p)[XDIM][YDIM][ZDIM],
                 const float scale_1,
                 const float scale_2) {
    for (int i = 1; i < XDIM - 1; i++)
        for (int j = 1; j < YDIM - 1; j++)
            for (int k = 1; k < ZDIM - 1; k++) {
                p[i][j][k] = x[i][j][k] * scale_1 + y[i][j][k];
                x[i][j][k] = z[i][j][k] + x[i][j][k] * scale_2;
            }
}
