If you discover any bugs or failure modes in this code, please email Victor (victoroknyansky@gmail.com) to let him know!! As per the license in the code, you are using this code at your own risk (i.e. we are not responsible for any problems/bugs, though we do our best to find and fix them).

If you use this code, please cite Gaskell & Spark 1986 (https://articles.adsabs.harvard.edu/pdf/1986ApJ...305..175G), Oknyanskii 1993 (https://articles.adsabs.harvard.edu/pdf/1993AstL...19..416O),
Peterson et al. 1998 (https://arxiv.org/abs/astro-ph/9802103),  http://ascl.net/code/v/1868 and https://ascl.net/code/v/3467
See details also here https://link.springer.com/content/pdf/10.1134/S1063773714090011.pdf

__Description of the code:__
PyMCCF emulates a PASCAL program written by V.L.Oknyansky for use with reverberation mapping (RM). The code is an updated version of the
Gaskell & Spark (1986) method ICCF. The code cross correlates two light curves that are unevenly sampled using linear interpolation, and
measures the peak and centroid of the cross-correlation function. The most significant improvement of the code is using a special parameter
MAX which reduces the number of used interpolated pointes just to those which are not further from the nearest real one than the MAX.
This way gives opportunities to significantly reduce noise from interpolation errors. Also in ICCF were introduced constant values at
the time space before the first and after the last points of the sets. The number of interpolated pairs in the ICCF were forever the same for
each time delay, but surely the signal/noise ratio drops down as the time delay grows. In our code the number correlated pairs of data are
not the same for each time delay. In the PyCCD version of the code published by Sun et al. also are used different numbers of pair for
each time delay due to not using any extrapolation of the data sets, but in the version are used interpolated points in the big gaps at
data sets which are very common for astronomical data. So we were based on the PyCCF code with some corrections which give us an
opportunity to introduce the new parameter MAX. Additionally we prefer to interpolate just the best one from the 2 data sets used for RM. The
method of getting errors for RM in our code PyMCCF were not changed and is exactly the same as in PyCCF. So it is possible, in addition,
to run Monte-Carlo iterations using flux randomization and random subset selection (RSS) to produce cross-correlation centroid
distributions to estimate the uncertainties in the cross correlation results.

__Python requirements__:

- Numpy (>1.4)
- Scipy (>0.1)
- Matplotlib (>1.0)

__File descriptions__:

PyMCCF.py is the actual python CCF software based on PYCCF project from [here](https://bitbucket.org/cgrier/python_ccf_code/src/master/) 
Peterson et al. 1998. We use this code to make our modernisation (Oknyanskii 1993) to the original method (Gaskell & Spark 1986) working as Python version.
The estimations of the errors were not changed and calculated the same way as in PYCCF.

`sample_lc1.dat` and `sample_lc2.dat` are required to run the sample_runcode.py script.

`PyMCCF.py` is the actual python CCF software.

`mccf_runner.py` is a script that contains a demonstration of how to use `PyMCCF.py`. It will output three data files and one plot to `results` directory based on input data from `input` folder.

`sample_lc1.dat` and `sample_lc2.dat` from `input` folder are required to run the `mccf_runner.py` script.

__Differences between `PyMCCF` and `PYCCF`__

- the most significant difference is the additional parameter MAX which can be set in `mccf_runner.py`. MAX is some compromise parameter which shows how far from the nearest real data point can the interpolated point be included in calculations 

- in `xcor(t1, y1, t2, y2, tlagmin, tlagmax, tunit, imode=0):` 

function block

``` 
while tau < tau_max:
    t2new = t1 + tau
    selin = np.where((t2new>=np.min(t2))&(t2new<=np.max(t2)), True, False)
```

was replaced with 

```
    while tau < tau_max:
        t2new = t1 + tau
        selin = np.where(
            [np.min(np.abs(t2new_mem - t2)) <= MAX and np.min(t2)-MAX <= t2new_mem <= np.max(t2)+MAX for t2new_mem in t2new],
            True, False)
```

- block 

```
while tau < tau_max:
    t1new = t2 - tau
    selin = np.where((t1new>=np.min(t1))&(t1new<=np.max(t1)), True, False)
```

was replaced with 

```
while tau < tau_max:
    t1new = t2 - tau
    selin = np.where(
        [np.min(np.abs(t1new_mem - t1)) <= MAX and np.min(t1)-MAX <= t1new_mem <= np.max(t1)+MAX for t1new_mem in t1new],
        True, False)
```
- in addition `xcor(t1, y1, t2, y2, tlagmin, tlagmax, tunit, imode=0)` gets new parameter `MAX=5` which is set in `mccf_runner.py`

- output folder `results` with inner folders `0`,`1`, `2` per `imode` which can be set in `mccf_runner.py`. Each such folder contains output files.
- output file `sample_ccf_yap.dat` contains one additional column. It is the number of selected pairs for getting correlation coefficient.
