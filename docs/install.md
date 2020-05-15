# Installation

## 1. Pre-requisites

1. cmake (>= 3.14)
2. C++-11 compatible compiler (gcc >= 5.0)
3. GDAL (>= 3.x)
4. [Armadillo](http://arma.sourceforge.net/) 
5. Lapack + Blas or similar implementation like openblas
6. Python (>= 3.7)
7. Cython


### 1.1 MacPorts on OSX

```
#clang compiler
> sudo port install clang-8.0
> sudo port select clang mp-clang-8.0

#Python and cython
> sudo port install python37 py37-cython

#Armadillo
> sudo port install armadillo

#gdal
> sudo port install gdal py37-gdal

#cmake
> sudo port install cmake
```

You will also need to install isce2 for using the workflows.

### 1.2 Anaconda on Linux / OSX

Armadillo installation should also install BLAS + Lapack.

Make sure the C++ compiler is new enough, otherwise, install the latest version using [conda](https://docs.conda.io/projects/conda-build/en/latest/resources/compiler-tools.html), i.e. `conda install gxx_linux-64` or `conda install clangxx_osx-64`. The latter will set an env variable `${CXX}` to the installed compiler.

Then the requirements file is as follows:

```
cmake
python
cython
gdal
libgdal
armadillo
isce2
```

## 2. Installation Process

### 2.1 Separate folders

The src, build and install folders need to be separate.

```
> mkdir install build src
```

### 2.2 Clone the repo

```
> cd src
> git clone https://github.com/isce-framework/fringe.git
```

### 2.3 Build and install

In the build folder

```
> cd build

#On OSX
> CXX=clang++ cmake -DCMAKE_INSTALL_PREFIX=../install ../src/fringe

#On Linux
> CXX=g++ cmake -DCMAKE_INSTALL_PREFIX=../install ../src/fringe

# with anaconda (on Linux) 
> CXX=g++ cmake -DCMAKE_INSTALL_PREFIX=../install ../src/fringe -DCMAKE_PREFIX_PATH=$CONDA_PREFIX

#If "conda install gxx_linux-64 / clangxx_osx-64"
> CXX=${CXX} cmake -DCMAKE_INSTALL_PREFIX=../install ../src/fringe

> make all
> make install
```

### 2.4 Set environment variables

Set the following variables after installation

```
export PATH=$PATH:$ISCE_HOME/bin:$ISCE_HOME/applications:path-to-install-folder/bin

export PYTHONPATH=$PYTHONPATH:path-to-install-folder/python

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path-to-install-folder/lib
```
If using within an Anaconda environment:
export GDAL_DIR=base-anaconda-directory/envs/fringe/

