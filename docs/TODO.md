# FRInGE To do list
-----

FRInGE is currently maintained by a very small team and was developed as a research prototype. There is a long laundry list of things to do and we attempt to capture it here. 

## 0. Testing framework
----

   - Convert exists tests to unit tests and add ctest support
   - Add CI capabilities


## 1. General software
-----

   - Code is still a research prototype
   - User needs to set PATH / PYTHONPATH like olden days
   - Package it better by building C++ code into a single lib and single cython extension
   - Single entry point instead of polluting the bin with multiple executables


## 2. nmap
-----

    - Generalize self-similarity comparison using class inheritance and virtual functions; to allow for easy inclusion of more methods. Currently AD2 and KS2 implementations available.
    - Add continuguous neighborhood option
    - CUDA code restructuring


## 3. evd
-----

    - CUDA support
    - Allow for using a small neighborhood for EVD than provided by NMAP. Currently, the neighborhood needs to be the same.
    - Ill-conditioned correltion matrices are now flagged during MLE and no attempt is made to regularize them. Consider implementing regularization.
    - Consider replacing amplitude of MLE solution with despeckled amplitude for each SLC. Implementation for this already available in despeck

