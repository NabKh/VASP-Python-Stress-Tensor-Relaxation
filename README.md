# VASP-Python-Stress-Tensor-Relaxation

# VASP Stress Tensor Relaxation Tool README

This repository hosts a custom tool designed to enhance the functionality of VASP (Vienna Ab initio Simulation Package) by enabling structure relaxation towards a predefined stress tensor, a feature not originally provided by standard VASP. The tool is capable of automatically modifying the POSCAR file and invoking VASP in a loop until the stress reaches the desired convergence criteria. It utilizes isotropic linear elastic theory for POSCAR adjustments. While suitable for many materials, performance may vary with highly anisotropic substances.

### Benchmarks

The tool has been successfully benchmarked on materials such as Tungsten (W), Molybdenum (Mo), Chromium (Cr), and alpha-Iron (α-Fe). However, caution is advised when dealing with highly anisotropic materials as convergence issues may arise. Implementing a full stiffness tensor (Cijkl) is a recommended solution for such cases.

## Usage Instructions

1. Insert the desired stress tensor values into the script.
2. Enter the estimated elastic modulus into the script.
3. In your VASP `INCAR` file, ensure that `ISIF=2`.
4. Duplicate your initial `POSCAR` file and rename the copy to `poscar.0`.
5. Edit the script to configure the `mpiexec` commands as per your VASP installation and job management system.
6. Navigate (`cd`) to your job directory.
7. Submit your job using the command `qsub fixpressure.py` (assuming the utilization of a PBS job system).

## Tips for Accurate Force Calculation

- Opt for a higher plane-wave cutoff energy (`ECUT`) and a dense k-point mesh for better convergence of forces and elastic constants, which tend to converge slowly with respect to k-points in VASP.
- Set `LREAL=.FALSE.` to avoid unrealistic forces that sometimes result from real space calculations.
- For large systems where reciprocal space calculations are computationally demanding, consider using `LREAL=Auto` coupled with `PREC=Accurate` and `ADDGRID=.TRUE.` to balance efficiency and accuracy.

## Recommendations for Computational Efficiency

- Employ `LWAVE=.TRUE.` along with the default settings for `ISTART` and `ICHARG` to enable VASP to reuse the `WAVECAR` file from the previous iteration, reducing the need to restart calculations from scratch.

## General Algorithm

1. Execute VASP with `ISIF=2` to determine the current stress tensor.
2. Apply generalized Hooke's law using the provided elastic modulus to calculate the necessary modifications to `POSCAR` for achieving the target pressure.
3. Repeat steps 1 and 2 until the discrepancy between the current and target pressure falls within the predefined convergence threshold.

---
