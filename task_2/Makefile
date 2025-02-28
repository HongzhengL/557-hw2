CXX = g++
CXXFLAGS = -O3 -g -std=c++20 -Wall -fopenmp

SRCS = main.cpp ConjugateGradients.cpp Laplacian.cpp PointwiseOps.cpp Reductions.cpp Utilities.cpp MergedOp.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = solver
TEMP_FILES = *.pgm

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) -fopenmp

main.o: main.cpp ConjugateGradients.h Laplacian.h Parameters.h PointwiseOps.h Reductions.h Timer.h Utilities.h MergedOp.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

ConjugateGradients.o: ConjugateGradients.cpp ConjugateGradients.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

Laplacian.o: Laplacian.cpp Laplacian.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

PointwiseOps.o: PointwiseOps.cpp PointwiseOps.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

Reductions.o: Reductions.cpp Reductions.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

Utilities.o: Utilities.cpp Utilities.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

MergedOp.o: MergedOp.cpp MergedOp.h
	$(CXX) $(CXXFLAGS) -fno-strict-aliasing -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET) $(TEMP_FILES)

.PHONY: all clean
