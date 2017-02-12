GCC=gcc
MPICC=mpicc


all: serial

serial:
	$(GCC) seq_main.c seq_camedoids.c file_io.c -lm -O3 -w -o seq_camedoids

omp: 
	$(GCC) omp_main.c omp_camedoids.c file_io.c -fopenmp -lm -O3 -w -o omp_camedoids

mpi:
	$(MPICC) seq_main.c seq_camedoids.c  file_io.c -lm -O3 -w -o mpi_camedoids

hybrid:
	$(MPICC) hybrid_main.c hybrid_camedoids.c mpi_io.c file_io.c -fopenmp -lm -O3 -w -o hybrid_camedoids

clean:
	rm -f seq_camedoids omp_camedoids mpi_camedoids hybrid_camedoids
