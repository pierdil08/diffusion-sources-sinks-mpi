#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Define the source function g(x, y) used in the PDE
// It models two Gaussian sources, one positive (centered at x=1) and one negative (centered at x=-1)
double g(double x, double y, double lambda) {
    return (10 * lambda / sqrt(M_PI)) *
           (exp(-lambda*lambda*((x-1)*(x-1) + y*y)) -
            exp(-lambda*lambda*((x+1)*(x+1) + y*y)));
}

// Perform SOR updates on red or black points in the grid
void update_points(double *u_data, int start_j, int end_j, int Nx, double dx2, double omega,
                   double *g_data, int is_red, double *max_resid) {
    for (int j = start_j; j <= end_j; ++j) {
        for (int k = 1; k < Nx - 1; ++k) {
            if ((j + k) % 2 == is_red) {
                double *u = &u_data[j * Nx];
                double *u_top = &u_data[(j + 1) * Nx];
                double *u_bot = &u_data[(j - 1) * Nx];
                double old = u[k];
                double neighbor_sum = u_top[k] + u_bot[k] + u[k + 1] + u[k - 1];
                double new_val = old + omega * (0.25 * neighbor_sum - old - 0.25 * dx2 * g_data[j * Nx + k]);
                double R = new_val - old;
                u[k] = new_val;
                *max_resid = fmax(*max_resid, fabs(R));
            }
        }
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv); // Initialize MPI

    // Record timing info
    double precision = MPI_Wtick();
    double starttime = MPI_Wtime();

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Parse input arguments
    int N = atoi(argv[1]);              // Grid height
    double omega = atof(argv[2]);       // Relaxation parameter
    double tol = atof(argv[3]);         // Convergence tolerance
    int max_iters = atoi(argv[4]);      // Maximum number of iterations

    int Ny = N;
    int Nx = 2 * N - 1;                 // Grid width
    double dx = 4.0 / (Nx - 1);
    double dy = 2.0 / (Ny - 1);
    double dx2 = dx * dx;

    // Decompose domain into horizontal blocks
    int rows_per_proc = Ny / size;
    int start_j = rank * rows_per_proc + 1;
    int count_j = rows_per_proc;
    int local_Ny = count_j + 2; // Add halo rows

    // Allocate memory for grid and source term (contiguous layout)
    double *u_data = calloc(local_Ny * Nx, sizeof(double));
    double *g_data = calloc(local_Ny * Nx, sizeof(double));

    // Compute source term values for local subdomain
    for (int j = 0; j < local_Ny; ++j) {
        double y = -1.0 + (start_j - 1 + j) * dy;
        for (int k = 0; k < Nx; ++k) {
            double x = -2.0 + k * dx;
            g_data[j * Nx + k] = g(x, y, 100);
        }
    }

    double global_max, local_max;
    int iter = 0;

    // Main SOR iteration loop
    while (iter < max_iters) {
        global_max = 0.0;

        // Loop over red and black passes
        for (int color = 0; color <= 1; ++color) {
            // Exchange boundary rows with neighbors
            if (rank > 0)
                MPI_Sendrecv(&u_data[1 * Nx], Nx, MPI_DOUBLE, rank-1, 0,
                             &u_data[0 * Nx], Nx, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (rank < size - 1)
                MPI_Sendrecv(&u_data[(local_Ny-2) * Nx], Nx, MPI_DOUBLE, rank+1, 0,
                             &u_data[(local_Ny-1) * Nx], Nx, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Apply Neumann boundary conditions at left and right edges
            for (int j = 1; j < local_Ny - 1; ++j) {
                u_data[j * Nx + 0] = u_data[j * Nx + 1];
                u_data[j * Nx + (Nx - 1)] = u_data[j * Nx + (Nx - 2)];
            }

            // Apply Dirichlet boundary conditions at top and bottom
            if (rank == 0) {
                for (int k = 0; k < Nx; ++k)
                    u_data[1 * Nx + k] = 0.0;
            }
            if (rank == size - 1) {
                for (int k = 0; k < Nx; ++k)
                    u_data[(local_Ny-2) * Nx + k] = 0.0;
            }

            // Update red or black points
            local_max = 0.0;
            update_points(u_data, 1, local_Ny - 2, Nx, dx2, omega, g_data, color, &local_max);

            // Get global maximum residual
            MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        }

        // Check convergence
        if (global_max < tol) break;
        iter++;
    }

    // Get pointer to flat data block (excluding ghost rows)
    double *u_flat = &u_data[1 * Nx];
    double *final_u = NULL;
    if (rank == 0)
        final_u = malloc(Ny * Nx * sizeof(double));

    // Gather data from all processes to root
    MPI_Gather(u_flat, count_j * Nx, MPI_DOUBLE,
               final_u, count_j * Nx, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    // Root writes the output file in column-major order
    if (rank == 0) {
        FILE *fout = fopen("Sources.out", "wb");
        for (int k = 0; k < Nx; ++k) {
            for (int j = 0; j < Ny; ++j) {
                fwrite(&final_u[j * Nx + k], sizeof(double), 1, fout);
            }
        }
        fclose(fout);
        printf("Converged in %d iterations.\n", iter);
    }

    // Free allocated memory
    free(u_data);
    free(g_data);
    if (rank == 0) free(final_u);

    // Print timing summary
    double elapsedtime = MPI_Wtime() - starttime;
    if (rank == 0) {
        printf("Execution time = %le seconds\n", elapsedtime);
        printf("Precision of timing is %le seconds\n", precision);
    }

    MPI_Finalize();
    return 0;
}

