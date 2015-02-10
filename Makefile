CC = g++

# Flags for -c compilation and linking
FLAG = -c -g -I/usr/include -I/usr/local/include
LFLAGS = -lgsl -lgslcblas -lm 
# -L/usr/lib64 
CXXFLAGS = -O2 -Wall -march=native

OBJS  = main.o \
	random.o benchmark.o \
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

ecga: $(OBJS)
	g++ -o rECGA_FHH $(OBJS) $(LFLAGS) 

clean:
	rm *.o rECGA_FHH
