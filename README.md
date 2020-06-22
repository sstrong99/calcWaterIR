Code to calculate the IR spectrum of water according to [Yang&Skinner PCCP 2010](https://pubs.rsc.org/en/content/articlelanding/2010/CP/B918314K). 

# Compiling

To compile, you will need the [xdrfile library](ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.4.tar.gz).

To compile for GPUs, you will also need CUDA and the [MAGMA library](https://icl.cs.utk.edu/magma/)

Then do:
cd src
make cpu OR make gpu

# Syntax

To run the program, simply pass one argument: the name of the input script.
The input script is a list of key-value pairs. The keys are described below.
The keys are required unless a default is listed, or specified otherwise.
The "#" symbol is used for comments

*calc - specifies the type of calculation (default=0)
.*0 - calculate the QM transition dipole autocorrelation function "exactly" (within the mixed quantum-classical approximation) by numerical integration of the Schroedinger equation (NISE)
.*1 - calculate the distribution of eigenfrequencies
.*2 - calculate the spectrum using the time-averaging approximation (TAA), see e.g. [Auer&Skinner JCP 2008](https://aip.scitation.org/doi/10.1063/1.2925258)
*trajFile - name of the file with the molecular dynamics trajectory
*outPostfix - postfix for output files defaults to no postfix
*intMethod - method for integrating. Only used forcalc=0. (default=0)
.*0 - exact diagonalization
.*n - for a positive integer n, uses nth order Adams-Bashforth integration
*nSample - number of samples to average over (default=-1)
.*0 - Start a new sample every timestep. That is, take the max number of samples available according to nTCF and the length of the trajectory
.*-1 - Start a new sample every nTCF timesteps. That is, take the max number of samples without reusing data
.*n - for a positive integer n, specifies the number of samples
*timestep - The timestep at which to read data from the trajectory file. If set to 0, will use every timestep in the trajectory file. If a nonzero value is specified, it must be larger than the timestep of the data in the file, and be a multiple of it. (default=0) 
*avgF - the average frequency to subtract from each oscillator (in cm-1) to remove the fast osscillations in the integration. If specified as 0, will calculate the average frequency from a single snapshot and use that value. Only used for calc=0. (default=0)
*nTCF - the number of timesteps to sample the transition dipole autocorrelation function. Only used for calc=0. (default=6*T1 in timesteps)
*T1 - The phenomonological line broadening parameter in ps. Not used for calc=1. (default=0.26)
*Tavg - The averaging time in ps for the TAA approximation. Only used for calc=2 (default=0.076)

#Examples

The three example input scripts, named in.*, demonstrate the functionality of the program. The example trajectory provided "traj.xtc" is an 20 ps NVE simluation of 300 SPC/E water molecules at 300 K, sampled every 5 fs. It was generated using GROMACS 4.5.5.

Running the in.spec example with CPU code on an i7 processor with 8GB of RAM takes about 16 minutes. This example uses a single sample. With multiple samples, the code will use OpenMP to parallelize over all available cores. This can be tuned using OMP_NUM_THREADS.

Running the same example with GPU code on an nvidia K80 GPU takes about 15 seconds. This will also use OpenMP to parallelize each sample over the avaialable CPU cores.

# Notes

The code automatically discriminates between SPC/E, TIP4P, and TIP4P/2005 models based on the geometry of the water molecule (see traj.cpp). Other models will probably crash the program, and will definitely not give correct results. Based on the water model, the program automatically chooses which map to use. New maps can easily be added by adding a class in the map.h file and following the format in the other classes there. Getting the code to use your map requires editing the code in the constructor of calcW.cpp.

As is, the code will use the [Auer&Skinner JCP 2008](https://aip.scitation.org/doi/10.1063/1.2925258) for SPC/E water and the [Gruenbaum et al JCTC 2013](https://pubs.acs.org/doi/10.1021/ct400292q) map for TIP4P and TIP4P/2005 water.
