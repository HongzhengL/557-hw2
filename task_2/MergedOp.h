#pragma once

#include "Parameters.h"

void MergedSaxpy(float (&x)[XDIM][YDIM][ZDIM],
                 const float (&y)[XDIM][YDIM][ZDIM],
                 const float (&z)[XDIM][YDIM][ZDIM],
                 float (&p)[XDIM][YDIM][ZDIM],
                 const float scale_1,
                 const float scale_2);

float MergedComputeLaplacianInnerProduct(float (&u)[XDIM][YDIM][ZDIM], float (&Lu)[XDIM][YDIM][ZDIM]);
