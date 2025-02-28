#include "MergedOp.h"

float MergedComputeLaplacianInnerProduct(float (&u)[XDIM][YDIM][ZDIM], float (&Lu)[XDIM][YDIM][ZDIM]) {
    double result = 0.;
#pragma omp parallel for reduction(+ : result)
    for (int i = 1; i < XDIM - 1; i++)
        for (int j = 1; j < YDIM - 1; j++)
            for (int k = 1; k < ZDIM - 1; k++) {
                Lu[i][j][k] = -6 * u[i][j][k] + u[i + 1][j][k] + u[i - 1][j][k] + u[i][j + 1][k] + u[i][j - 1][k] +
                              u[i][j][k + 1] + u[i][j][k - 1];
                result += static_cast<double>(u[i][j][k]) * static_cast<double>(Lu[i][j][k]);
            }

    return static_cast<float>(result);
}
