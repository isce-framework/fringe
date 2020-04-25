# FRInGE Workflows
------

The process for generating InSAR time-series products with FRInGE can be broken down into two phases:

1. Coregistered SLC stack to PS+DS Wrapped Phase time-series
2. Wrapped time-series to deformation time-series


## Phase 1: Coregistered SLC stack to PS+DS Wrapped Phase time-series
-----------

FRInGE requires a coregistered stack of SLC images to work with. One may use different software to create the coregistered stack of SLCs. 

### Prepare the inputs for FRInGE 

Assuming ISCE2 stack processor has been used to create the coregistered stack of Sentinel-1 SLCs, the `tops2vrt.py` command can be used to create a VRT file which points to a subset of the dataset for further processing with FRInGE:

```
tops2vrt.py -i ../slc_stack/merged/ -s coreg_stack -g geometry -c slcs -B 37.85 38.034 -121.9 -122.15
```

### Finding the local neighborhood of each pixel

This step estimates the local neighborhood of every pixel in the stack of SLCs. The size of the local neighborhood in x and y directions can be determined. Basically the PDF of backscatter amplitude of the center pixel in the current window of x by y pixels is compared against the PDF of all other pixels in the current window and based on different statistical tests, the similiarity is decided: 

```
nmap.py -i coreg_stack/slcs_base.vrt -o KS2/nmap -c KS2/count -x 11 -y 5
```

### Estimate wrapped phase time-series for DS pixels

When the local neighborhood map is ready, one can estimate the wrapped phase time-series for DS pixels using the full covariance matrix of the entire stack using the `evd.py` command. Alternatively one may choose to break the stack to smaller mini-stacks and estimate a wrapped phase series for each mini-stack which can then be connected to each other using compressed SLCs: 

```
sequential.py -i ../slc_stack/merged/SLC -s 15 -o Sequential -w KS2/nmap -b 3200 4907 33084 38459 -x 11 -y 5
```

The wrapped phase series of mini-stacks can be then connected using the `adjustMiniStacks.py` command:

```
adjustMiniStacks.py -s slcs/ -m Sequential/miniStacks/ -d Sequential/Datum_connection/ -M 15 -o adjusted_wrapped_DS
```

### PS selection

To compute the amplitude dispersion of different pixels:

```
ampdispersion.py -i coreg_stack/slcs_base.vrt -o ampDispersion/ampdispersion -m ampDispersion/mean
```

To choose the PS pixels one can simply threshold the amplitude dispersion using a custom script. Here is an example using ISCE's imageMath function to threshold the PS pixels with a threshold of 0.4: 

```
imageMath.py -e="a<0.4" --a=ampDispersion/ampdispersion  -o ampDispersion/ps_pixels -t byte
```

### Integrate PS pixels to the wrapped phase series

For PS pixels we extract the wrapped phase series of the SLCs through time and integrate to the maps of wrapped phase time-series from the previous step. As a result of this step the wrapped phase time-series of PS and DS pixels are obtained in which the phase series of DS pixels were estimated from full covariance analysis over local neighborhoods as explained above and for PS pixels the wrapped phase series was extracted from SLCs.

```
integratePS.py -s coreg_stack/slcs_base.vrt -d adjusted_wrapped_DS/ -t Sequential/Datum_connection/EVD/tcorr.bin -p ampDispersion/ps_pixels -o PS_DS --unwrap_method snaphu
```

## Phase 2: Wrapped time-series to deformation time-series
-----------

FRInGE currently is focused on estimating wrapped phase time-series for PS and DS pixels as explained above. When the wrapped phase time-series is estimated, then user may use their favorite tools/algorithms to unwrap the wrapped phase time-series apply different corrections (atmospheric delay correction, tropospheric delay, DEM error, etc) and convert to deformation time-series. For example one may use Snaphu for unwrapping and then use regular time-series tools such as PyAPS and MintPy for atmospheric correction and displacement time-series analysis

### Unwrap the wrapped phase series

To prepare unwrapping commands for each epoch of time-series run `unwrapStack.py`. This will write the unwrapping command to a shell script _run_unwrap.sh_.

```
unwrapStack.py -s slcs -m Sequential/miniStacks/ -d Sequential/Datum_connection/ -M 15 -u 'unwrap_fringe.py -m snaphu'
```

Alternatively one may use `integratePS.py --unw_method snaphu` in previous step, which will write a shell script _run_unwrap_ps_ds.sh_ (recommended).

Each row of the run file is then the unwrapping command for each epoch which is indpendent from other epochs and therefore unwrapping different epochs can be run in parallel if computing resourses such as clusters are available.  

TODO: 
ingest FRInGE unwrapped phase series to MintPy ...
