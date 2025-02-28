#pragma once

#include "Parameters.h"

void MergedSaxpy(float (&x)[XDIM][YDIM][ZDIM],
                 float (&p)[XDIM][YDIM][ZDIM],
                 const float (&z)[XDIM][YDIM][ZDIM],
                 const float alpha,
                 const float beta);

float MergedComputeLaplacianInnerProduct(float (&u)[XDIM][YDIM][ZDIM], float (&Lu)[XDIM][YDIM][ZDIM]);
