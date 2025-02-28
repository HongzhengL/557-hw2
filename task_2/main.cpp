#include <string>
#include <vector>

#include "ConjugateGradients.h"
#include "Timer.h"
#include "Utilities.h"

#define ITERATIONS 256

Timer laplacian_timer_1;
Timer laplacian_inner_product_timer;

Timer copy_timer_1;
Timer copy_timer_2;

Timer inner_product_timer_1;
Timer inner_product_timer_3;

Timer norm_timer_1;
Timer norm_timer_2;

Timer saxpy_timer_1;
Timer saxpy_timer_2;
Timer saxpy_timer_3;
Timer saxpy_timer_4;

int main(int argc, char *argv[]) {
    using array_t = float(&)[XDIM][YDIM][ZDIM];

    float *xRaw = new float[XDIM * YDIM * ZDIM];
    float *fRaw = new float[XDIM * YDIM * ZDIM];
    float *pRaw = new float[XDIM * YDIM * ZDIM];
    float *rRaw = new float[XDIM * YDIM * ZDIM];
    float *zRaw = new float[XDIM * YDIM * ZDIM];

    array_t x = reinterpret_cast<array_t>(*xRaw);
    array_t f = reinterpret_cast<array_t>(*fRaw);
    array_t p = reinterpret_cast<array_t>(*pRaw);
    array_t r = reinterpret_cast<array_t>(*rRaw);
    array_t z = reinterpret_cast<array_t>(*zRaw);

    std::vector<Timer *> timers;
    timers.push_back(&laplacian_timer_1);
    timers.push_back(&laplacian_inner_product_timer);
    timers.push_back(&copy_timer_1);
    timers.push_back(&copy_timer_2);
    timers.push_back(&inner_product_timer_1);
    timers.push_back(&inner_product_timer_3);
    timers.push_back(&norm_timer_1);
    timers.push_back(&norm_timer_2);
    timers.push_back(&saxpy_timer_1);
    timers.push_back(&saxpy_timer_2);
    timers.push_back(&saxpy_timer_3);
    timers.push_back(&saxpy_timer_4);

    // clang-format off
    std::vector<std::string> timer_names = {
        "ComputeLaplacian(x, z)                   ",
        "MergedComputeLaplacianInnerProduct(p, z) ",
        "Copy(r, p)                               ",
        "Copy(r, z)                               ",
        "InnerProduct(p, r)                       ",
        "InnerProduct(z, r)                       ",
        "Norm(r)                                  ",
        "Norm(r)                                  ",
        "Saxpy(z, f, r, -1)                       ",
        "Saxpy(z, r, r, -alpha)                   ",
        "Saxpy(p, x, x, alpha)                    ",
        "MergedSaxpy(x, p, r, alpha, beta)        "
    };
    // clang-format on

    // Initialization
    {
        Timer timer;
        timer.Start();
        InitializeProblem(x, f);
        timer.Stop("Initialization : ");
    }

    for (auto &timer : timers) {
        timer->Reset();
    }

    // Call Conjugate Gradients algorithm
    {
        Timer timer;
        timer.Start();
        ConjugateGradients(x, f, p, r, z);
        timer.Stop("Conjugate Gradients : ");
        std::cout << std::endl;
    }

    std::cout << "----------------------------------" << std::endl;
    std::cout << "Per Kernel Cumulative:" << std::endl;
    for (size_t i = 0; i < timers.size(); i++) {
        timers[i]->Print(timer_names[i]);
    }

    // Per Iteration Average
    std::cout << "----------------------------------" << std::endl;
    std::cout << "Per Iteration Average:" << std::endl;
    for (size_t i = 0; i < timers.size(); i++) {
        std::cout << timer_names[i] << " : " << timers[i]->mElapsedTime.count() / (float)ITERATIONS << "ms"
                  << std::endl;
    }

    delete[] xRaw;
    delete[] fRaw;
    delete[] pRaw;
    delete[] rRaw;
    delete[] zRaw;

    return 0;
}
