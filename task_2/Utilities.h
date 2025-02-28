#pragma once

#include <string>

#include "Parameters.h"

void Clear(float (&x)[XDIM][YDIM][ZDIM]);
void InitializeProblem(float (&x)[XDIM][YDIM][ZDIM], float (&b)[XDIM][YDIM][ZDIM]);
void WriteAsImage(const std::string& filenamePrefix, const float (&x)[XDIM][YDIM][ZDIM], const int count,
                  const int axis, const int slice);
