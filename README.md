# 2024_PNAS_Schlissel_etal
diffusion simulation code from 2024_PNAS_Schlissel_etal

This repositotory incldues MATLAB code to simulate diffusion of a particle on a 2D lattice

The code was written in MATLAB, but the output was captured and visualized using R.

there is a ZIP'd directory of the MATLAB output data in each simulation directory, which 
can be decompressed and used directly in the plotting script

Generallly... 
- use the slur_call_steadyState.sh scripts to initiate SLURM jobs for each parameter set
  - update the various filepaths as neccesssary in the slur_call_steadyState.sh script
- the SLURM caller script will call many instances of the steady_state_toy.m script (in MATLAB)
- the steady_state_toy.m script depends on the steady_diffusionModel_toy.m file to initalize functions
- The data will be dumped to an output folder called ./steady_state_test_data/ and the filenames will include some metadata
- the R script handle_matlab_data_toy_topo_pub.R will read in the file series from the MATLAB output and plot them
