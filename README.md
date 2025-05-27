# Diffusion with Sources and Sinks (SOR, MPI, C)

This project solves a 2D steady-state diffusion equation with internal sources and sinks using the Successive Over-Relaxation (SOR) method. The simulation models heat distribution on a rectangular plate with insulated boundaries and two internal Dirac-type perturbations. Implemented in C and executed on Northwestern’s high-performance computing system using MPI parallelism.

## 🔬 Equation

The elliptic PDE solved is:

\[
\nabla^2 u = \frac{10\lambda}{\sqrt{\pi}} e^{-\lambda^2((x-1)^2+y^2)} - \frac{10\lambda}{\sqrt{\pi}} e^{-\lambda^2((x+1)^2+y^2)}
\]

Subject to boundary conditions:

\[
\frac{\partial u}{\partial x}(\pm2, y) = 0, \quad u(x, \pm1) = 0
\]

where:
- `λ` controls the width of the source/sink peaks
- The left and right boundaries are insulated
- The top and bottom are held at fixed temperature

## ⚙️ Tools Used
- **Language**: C
- **Libraries**: MPI
- **Numerical Method**: SOR
- **HPC System**: Quest @ Northwestern
- **Precision**: Double (error tolerance = 1e-9)

## 🚀 Features
- Parallel SOR solver using domain decomposition
- Adjustable grid size and convergence tolerance
- Performance timing for scalability analysis
- Output data in `Sources.out` with shape (2N-1) × N

## 🛠️ How to Run

```bash
make
mpirun -np 4 ./diffusion 128 1.8 1e-9 10000
