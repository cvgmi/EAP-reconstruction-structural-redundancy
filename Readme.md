# EAP Reconstruction from Highly Undersampled (k,q)-space in diffusion MRI

Jiaqi Sun, Alireza Entezari, Baba C Vemuri 
Exploiting structural redundancy in q-space for improved EAP reconstruction from highly undersampled (k,q)-space in DMRI
Medical Image Analysis Volume 54 pgs. 122-137

This code is a MATLAB implementation of the above paper on EAP reconstruction from highly undersampled data by
exploiting q-space redundancy.

## Dependencies

* [Surfacelet Toolbox by Yue Lu](https://www.mathworks.com/matlabcentral/fileexchange/14485-surfacelet-toolbox)
* [FMINLBFGS: Fast Limited Memory Optimizer by Dirk-Jan Kroon](https://www.mathworks.com/matlabcentral/fileexchange/23245-fminlbfgs-fast-limited-memory-optimizer)

## Quickstart

Please put the `SurfBox` directory inside of `./Surfacelet` and the `fminlbfgs_version2c` directory in the main
directory. To compile the Surfacelet toolbox (highly recomended), modify and run `mexcompile.m` in the `SurfBox` directory. 

An example script is provided in `example.m` and can be run immediately. 
The main directory contains a synthetic sample signal `Sq_test_12.mat` and perfect 
reconstruction `Pr_test_12.mat` along with a set of sampling masks for varying sampling rates in the `masks` directory.
After running `example.m`, the reconstruction will be stored in the MATLAB variable `Pr_recon1`. 
The main function is contained in `mrics5_PLS_TV_init_515_objhalf` and a brief explanation of each of the parameters are
provided in that function. 
