#pragma once

#include "Parameters.h"

float MergedComputeLaplacianInnerProduct(float (&u)[XDIM][YDIM][ZDIM], float (&Lu)[XDIM][YDIM][ZDIM]);
