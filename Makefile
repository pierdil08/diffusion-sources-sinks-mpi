all:
	mpicc -o diffusion main.c -lm

submit:
	sbatch myjob.sh

