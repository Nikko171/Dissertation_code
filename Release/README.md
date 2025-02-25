#README For python_ccf_code 

If you discover any bugs or failure modes in this code, please email Kate (or Mouyuan) to let her (him) know!! As per the license in the code, you are using this code at your own risk (i.e. we are not responsible for any problems/bugs, though we do our best to find and fix them).

We provide two versions: a Pure Python Version; a Cython plus Python Version (5~10 times faster)

The usage of the Cython version is exactly the same as the pure Python version. The only additional thing you have to do is compile the Cython code once (which is super easy, see the end of this ReadMe)!

File descriptions:

In the Pure Python folder:
PYCCF.py is the actual python CCF software.

sample_runcode.py is a script that contains a demonstration of how to use PYCCF.py. It will output three data files and one plot to the current directory.

sample_lc1.dat and sample_lc2.dat are required to run the sample_runcode.py script.


In the CPython folder (from our benchmarks: this version is 5~10 times faster than the pure python version):
PYCCF.py is the actual python CCF software.

xcorspc.pyx is a Cython XCOR function; you have to install Cython and a C compiler (e.g., GCC)

setup.py is a configuration file; use this file to compile the Cython XCOR function

sample_runcode.py is a script that contains a demonstration of how to use CCF_interp.py. It will output three data files and one plot to the current directory.

sample_lc1.dat and sample_lc2.dat are required to run the sample_runcode.py script.

To compile the Cython code, please follow the following steps:
1. Change your terminal directory to the ICCF folder

2. running " python setup.py build_ext --inplace " and check if the compiling is finished without errors

3. copy the " xcorspc.so " and " PYCCF.py " to your work directory and enjoy the code (for illustration, see the sample_runcode.py file).

