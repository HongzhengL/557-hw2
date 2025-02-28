#include <iostream>

#include "Laplacian.h"
#include "MergedOp.h"
#include "Parameters.h"
#include "PointwiseOps.h"
#include "Reductions.h"
#include "Timer.h"
#include "Utilities.h"

extern Timer laplacian_timer_1;
extern Timer laplacian_timer_2;

extern Timer copy_timer_1;
extern Timer copy_timer_2;

extern Timer inner_product_timer_1;
extern Timer inner_product_timer_2;
extern Timer inner_product_timer_3;

extern Timer norm_timer_1;
extern Timer norm_timer_2;

extern Timer saxpy_timer_1;
extern Timer saxpy_timer_2;
extern Timer saxpy_timer_3;
extern Timer saxpy_timer_4;

void ConjugateGradients(float (&x)[XDIM][YDIM][ZDIM],
                        const float (&f)[XDIM][YDIM][ZDIM],
                        float (&p)[XDIM][YDIM][ZDIM],
                        float (&r)[XDIM][YDIM][ZDIM],
                        float (&z)[XDIM][YDIM][ZDIM],
                        const bool writeIterations) {
    // Algorithm : Line 2
    laplacian_timer_1.Restart();
    ComputeLaplacian(x, z);
    laplacian_timer_1.Pause();
    saxpy_timer_1.Restart();
    Saxpy(z, f, r, -1);
    saxpy_timer_1.Pause();
    norm_timer_1.Restart();
    float nu = Norm(r);
    norm_timer_1.Pause();

    // Algorithm : Line 3
    if (nu < nuMax) return;

    // Algorithm : Line 4
    copy_timer_1.Restart();
    Copy(r, p);
    copy_timer_1.Pause();
    inner_product_timer_1.Restart();
    float rho = InnerProduct(p, r);
    inner_product_timer_1.Pause();

    // Beginning of loop from Line 5
    for (int k = 0;; k++) {
        std::cout << "Residual norm (nu) after " << k << " iterations = " << nu << std::endl;

        // Algorithm : Line 6
        laplacian_timer_2.Restart();
        // ComputeLaplacian(p, z);
        float sigma = MergedComputeLaplacianInnerProduct(p, z);
        laplacian_timer_2.Pause();

        // Algorithm : Line 7
        float alpha = rho / sigma;

        // Algorithm : Line 8
        saxpy_timer_2.Restart();
        Saxpy(z, r, r, -alpha);
        saxpy_timer_2.Pause();
        norm_timer_2.Restart();
        nu = Norm(r);
        norm_timer_2.Pause();

        // Algorithm : Lines 9-12
        if (nu < nuMax || k == kMax) {
            saxpy_timer_3.Restart();
            Saxpy(p, x, x, alpha);
            saxpy_timer_3.Pause();
            std::cout << "Conjugate Gradients terminated after " << k << " iterations; residual norm (nu) = " << nu
                      << std::endl;
            if (writeIterations) WriteAsImage("x", x, k, 0, 127);
            return;
        }

        // Algorithm : Line 13
        copy_timer_2.Restart();
        Copy(r, z);
        copy_timer_2.Pause();
        inner_product_timer_3.Restart();
        float rho_new = InnerProduct(z, r);
        inner_product_timer_3.Pause();

        // Algorithm : Line 14
        float beta = rho_new / rho;

        // Algorithm : Line 15
        rho = rho_new;

        // Algorithm : Line 16
        saxpy_timer_4.Restart();
        MergedSaxpy(p, x, r, x, alpha, beta);
        saxpy_timer_4.Pause();

        // if (writeIterations) WriteAsImage("x", x, k, 0, 127);
    }
}
