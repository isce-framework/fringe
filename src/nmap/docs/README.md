# Local neighborhood mapping module

This module estimates the local neighborhood of every pixel in a stack of SLCs.
The stack is expected to be provided as a GDAL VRT file. If an amplitude calibration constant is available in the metadata, it is applied before computing similarity between pixels.


## Algorithm

The algorithm estimates the  Goodness-of-Fit between amplitude distribution of a given pixel and every other pixel within a defined distance. Currently, Kolmogorov-Smirnov and Anderson-Darling 2-sample tests have been implemented. 

The code also accepts an optional input mask. Any pixel which does not have valid data in any of the coregistered SLC is also masked. The output consists of a bit-wise neighbor mask (stored within UInt32 layers) and a total count of self-similar neighbors for each pixel. Every valid pixel is labeled as its own neighbor. If this flag is set to zero, the pixel is equivalent to being masked out.


### Kolmogorov Smirnov 2-sample test

We compared numerous implementations of KS2 test available in the public domain. A brief description can be found [here](../../../tests/KS2/README.md)


### Anderson Darling 2-sample test

We compared numerous implementations of AD2 test available in the public domain. A brief description can be found [here](../../../tests/AD2/README.md)


### References

1. Ferretti, A., Fumagalli, A., Novali, F., Prati, C., Rocca, F. and Rucci, A., 2011. A new algorithm for processing interferometric data-stacks: SqueeSAR. IEEE Transactions on Geoscience and Remote Sensing, 49(9), pp.3460-3470.

2. Parizzi, A. and Brcic, R., 2011. Adaptive InSAR stack multilooking exploiting amplitude statistics: A comparison between different techniques and practical results. IEEE Geoscience and Remote Sensing Letters, 8(3), pp.441-445.



## nmap.py 

Python-based script that uses Cython-bound version of nmap.
Python-based executable will always have more error checking and better handling of optional inputs. 

```
usage: nmap.py [-h] -i INPUTDS -o OUTPUTDS -c COUNTDS [-m MASKDS]
               [-l LINESPERBLOCK] [-r MEMORYSIZE] [-x HALFWINDOWX]
               [-y HALFWINDOWY] [-p PVALUE] [-s METHOD]

Create neighborhood mask and count map using KS statistics for stack of coregistered SLCs

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTDS, --input INPUTDS
                        Input GDAL SLC stack VRT (default: None)
  -o OUTPUTDS, --output OUTPUTDS
                        Output neighborhood weights mask (default: None)
  -c COUNTDS, --count COUNTDS
                        Output count dataset (default: None)
  -m MASKDS, --mask MASKDS
                        Optional mask layer to speed up computation (default:
                        )
  -l LINESPERBLOCK, --linesperblock LINESPERBLOCK
                        Quantum for block of lines (default: 64)
  -r MEMORYSIZE, --ram MEMORYSIZE
                        Memory in Mb to use (default: 256)
  -x HALFWINDOWX, --xhalf HALFWINDOWX
                        Half window size (range) (default: 5)
  -y HALFWINDOWY, --yhalf HALFWINDOWY
                        Half window size (azimuth) (default: 5)
  -p PVALUE, --prob PVALUE
                        Minimum p-value for labeling neighbors. (default:
                        0.95)
  -s METHOD, --stat METHOD
                        Statistical test to use - KS2 or AD2 (default: KS2)
  --nogpu               Do not use GPU
```



## nmap

C++ based executable. Not much error checking of inputs.

```
  nmap {OPTIONS}

    Create neighborhood mask and count map using KS statistics

  OPTIONS:

      -h, --help                        Display this help menu
      -i[inputDS*]                      Input Stack VRT
      -o[outputDS*]                     Output neighbor mask
      -c[countDS*]                      Neighbor count dataset
      -m[maskDS*]                       Mask dataset
      -l[linesperblock]                 Lines per block for processing
      -r[memorysize]                    Memory in Mb
      -x[xsize]                         Half window in x
      -y[ysize]                         Half window in y
      -p[prob]                          P-value threshold
      -t[testMethod]                    Statistics test - KS2 or AD2
      --nogpu                           Do not use GPU
```


## Display output with mdx package
Note that neighbormask contains the neighborhood information stored in multiple bands as bit information.
It can therefore not easily be visualized. You can visualize the countdataset by using the WIDTH and precision.
This information is contained in the corresponding xml or hdr files. 

```
mdx -s WIDTH -i2 countdataset
```

### [Return to main overview](../../../README.md)

