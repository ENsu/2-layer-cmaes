CC = g++

# Flags for -c compilation and linking
FLAG = -c -g -I/usr/include -I/usr/local/include
LFLAGS = -lgsl -lgslcblas -lm 
# -L/usr/lib64 
CXXFLAGS = -O2 -Wall -march=native

OBJS  = main.o \
	benchmark.o cmaes.o group.o node.o global.o random.o \
        F01_shifted_sphere.o \
        F02_shifted_schwefel.o \
        F03_shifted_rotated_high_cond_elliptic.o \
        F04_shifted_schwefel_noise.o \
        F05_schwefel_global_opt_bound.o \
        F06_shifted_rosenbrock.o \
        F07_shifted_rotated_griewank.o \
        F08_shifted_rotated_ackley_global_opt_bound.o \
        F09_shifted_rastrigin.o \
        F10_shifted_rotated_rastrigin.o \
        F11_shifted_rotated_weierstrass.o \
        F12_schwefel.o \
        F13_shifted_expanded_griewank_rosenbrock.o \
        F14_shifted_rotated_expanded_scaffer.o \
        F15_hybrid_composition_1.o \
        F16_rotated_hybrid_composition_1.o \
        F17_rotated_hybrid_composition_1_noise.o \
        F18_rotated_hybrid_composition_2.o \
        F19_rotated_hybrid_composition_2_narrow_basin_global_opt.o \
        F20_rotated_hybrid_composition_2_global_opt_bound.o \
        F21_rotated_hybrid_composition_3.o \
        F22_rotated_hybrid_composition_3_high_cond_num_matrix.o \
        F23_noncontinuous_rotated_hybrid_composition_3.o \
        F24_rotated_hybrid_composition_4.o \
        F25_rotated_hybrid_composition_4_bound.o\

cmaes: $(OBJS)
	g++ -o cmaes $(OBJS) $(LFLAGS) 

clean:
	rm *.o cmaes

benchmark.o: benchmark.cpp benchmark.h random.hpp test_func.h HCJob.h 
	$(CC) $(FLAG) benchmark.cpp

cmaes.o: cmaes.cpp cmaes.h group.h global.h random.hpp
	$(CC) $(FLAG) cmaes.cpp

group.o: group.cpp group.h node.h global.h random.hpp
	$(CC) $(FLAG) group.cpp

node.o: node.cpp node.h global.h
	$(CC) $(FLAG) node.cpp
