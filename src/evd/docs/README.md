# Maximum Likelihood Estimator using Eigen Value Decomposition

This module estimates the maximum likelihood stack from a given stack of SLCs and a neighborhood map. The stack is expected to be provided as a GDAL VRT file and the neighborhood weights as a single GDAL dataset. 

The output is a complete folder with a single binary dataset for each SLC in the stack. In addition, two other outputs are included in the folder:

1) tcorr.bin - The posteriori coherence
2) compslc.bin - The compressed SLC corresponding to the entire stack


For full details on the processing approach, see references.

## Algorithm

1. For each pixel, a Covariance matrix in time is estimated using the local neighorhood.

2. The Covariance matrix is then used with the MLE Estimator to estimate a wrapped phase for each time-epoch in 3 steps
    - Coherence maxtrix is estimated as element-by-element absolute value of covariance matrix
    - The inverse of the coherence matrix is estimated
    - The Hadamard product of the inverse and Covariance matrix is used for further analysis

3. We have not yet implemented iterative improvements to the phase-linking algorithm. For now, we use the phase of the Eigen vector corresponding to the smallest Eigen value as the MLE estimate. Simple tests with Python show that the iterative improvements are fairly minor in most cases.

4. These MLE phase estimates are written out as individual layers. The other related measures like the temporal coherence and compressed SLC are also estimated on a pixel-by-pixel basis and written out.


### Matrix operations

1. We use armadillo cubes / matrices as containers for storing data.

2. For numerical operations, armadillo does not include an efficient way of generating a single eigen vector.

3. Hence, we bound Lapack's [zheevr](http://www.netlib.org/lapack/explore-html/d9/dd2/zheevr_8f.html) method to efficiently estimate single Eigen vectors of interst. 

4. For Hermitian matrix inverse, we use [zportf](http://www.netlib.org/lapack/explore-3.1.1-html/zpotrf.f.html) followed by [zpotri](http://www.netlib.org/lapack/explore-3.1.1-html/zpotri.f.html).

5. All the symbol bindings can be found under [include/fringe/EigenLapack.hpp](../../../include/fringe/EigenLapack.hpp) in the code base. 

6. EVD will need to be linked against lapack and blas/f77blas. On OS X, we rely on atlas to provide these packages. On linux, we have tested this with conda. 


### OpenMP issues

1. Note that most distributions of armadillo ship with OpenBLAS included and it is possible to use openMP optimizations. When armadillo linalg operations are called from within OpenMP loops, they raise warning / error messages. 

2. To turn these off, set environment variable OPENBLAS_NUM_THREADS=1. Note that setting this environment variable permanently may slow down your other armadillo dependent packages. Hence, suggested way to deal with this to include the environment variable as part of your command (see below for usage).

3. We are still working on a way to have python set this environment variable in a reliable way. Looks like the current way of setting the environment variable before importing the cython bound library works but needs verification on different computing platforms.


### References

1. Ansari, H., De Zan, F., Adam, N., Goel, K., & Bamler, R. (2016, July). Sequential estimator for distributed scatterer interferometry. In Geoscience and Remote Sensing Symposium (IGARSS), 2016 IEEE International (pp. 6859-6862). IEEE.

2. Ansari, H., De Zan, F., & Bamler, R. (2017). Sequential Estimator- A Proposal for High-Precision and Efficient Earth Deformation Monitoring with InSAR [FRINGE 2017 Presentation](http://fringe2017.esa.int/files/presentation324.pdf) 

3. A. Monti Guarnieri and S. Tebaldini, “On the Exploitation of Target Statistics for SAR Interferometry Applications,” IEEE Trans. Geosci. Remote Sens., vol. 46, no. 11, pp. 3436–3443, Nov. 2008.

4. A. M. Guarnieri and S. Tebaldini, “Hybrid Cramér-Rao Bounds for Crustal Displacement Field Estimators in SAR Interferometry,” IEEE Signal Process. Lett., vol. 14, no. 12, pp. 1012–1015, Dec. 2007.

5. A. Ferretti, A. Fumagalli, F. Novali, C. Prati, F. Rocca, and A. Rucci, “A New Algorithm for Processing Interferometric Data-Stacks: SqueeSAR,” IEEE Trans. Geosci. Remote Sens., vol. 49, no. 9, pp. 3460–3470, Sep. 2011.

6. Y. Wang and X. X. Zhu, “Robust Estimators for Multipass SAR Interferometry,” IEEE Trans. Geosci. Remote Sens., vol. 54, no. 2, pp. 968–980, Feb. 2016.


## evd.py 

Python-based script that uses Cython-bound version of nmap.
Python-based executable will always have more error checking and better handling of optional inputs. 

```
usage: evd.py [-h] -i INPUTDS -w WTSDS -o OUTPUTFOLDER [-l LINESPERBLOCK]
              [-r MEMORYSIZE] [-x HALFWINDOWX] [-y HALFWINDOWY]
              [-m MINNEIGHBORS]

Perform MLE-based phase-linking on a stack of coregistered SLCs

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTDS, --input INPUTDS
                        Input GDAL SLC stack VRT (default: None)
  -w WTSDS, --wts WTSDS
                        Input weights dataset (default: None)
  -o OUTPUTFOLDER, --output OUTPUTFOLDER
                        Output folder with phase-linked SLCs (default: None)
  -l LINESPERBLOCK, --linesperblock LINESPERBLOCK
                        Quantum for block of lines (default: 64)
  -r MEMORYSIZE, --ram MEMORYSIZE
                        Memory in Mb to use (default: 2048)
  -x HALFWINDOWX, --xhalf HALFWINDOWX
                        Half window size (range) (default: 5)
  -y HALFWINDOWY, --yhalf HALFWINDOWY
                        Half window size (azimuth) (default: 5)
  -m MINNEIGHBORS, --minneigh MINNEIGHBORS
                        Minimum number of neighbors for computation (default:
                        5)
```



## evd

C++ based executable. Not much error checking of inputs.

```
  evd {OPTIONS}

    Eigen value decomposition of a stack of SLCs

  OPTIONS:

      -h, --help                        Display this help menu
      -i[inputDS*]                      Input Stack VRT
      -w[wtsDS*]                        Input neighbor mask
      -o[output*]                       Output folder for new stack
      -l[linesperblock]                 Lines per block for processing
      -r[memorysize]                    Memory in Mb
      -x[xsize]                         Half window in x
      -y[ysize]                         Half window in y
      -n[minNeighbors]                  Minimum number of neighbors
```


Note that I often have to call the executable as shown below on conda-based environments:

```
> OPENBLAS_NUM_THREADS=1 LD_LIBRARY_PATH=$LD_RUN_PATH evd
```

### evdtest.py

The code base includes a script called "evdtest.py" to allow developers to compare data and the results of various matrix computatios in C++ and python. This is based on the MLE python scripts in the tests folder.


### [Return to main overview](../../../README.md)

