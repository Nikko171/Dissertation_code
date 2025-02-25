'''
Version 1, by Mouyuan (Eric) Sun
email: ericsun88@live.com

-----------------------purpose--------------------------------------
The purpose of this code is to use Cython to speedup the cross
correlation calculation.

Avoid using numpy modules

Requirements: Cython (http://docs.cython.org/en/latest/index.html)
              a modern C compiler
--------------------------------------------------------------------
'''
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)

cdef double linearinterp(double x0, int n, double[:] x, double[:] y):
    '''
    general linear interpolation routine
    '''
    cdef double yinterp
    cdef int i
    if x0<=x[0]:
        yinterp = y[0] # lower boundary
    elif x0>=x[n-1]:
        yinterp = y[n-1] # upper boundary
    else:
        i = 1
        while i<n:
            if x0==x[i]:
                yinterp = y[i]
                i = n # exit the loop
            elif x0 < x[i] :
                yinterp = y[i-1] + (x0-x[i-1])*(y[i]-y[i-1])/(x[i]-x[i-1])
                i = n # exit the loop
            else:
                i = i + 1
    return yinterp
    


cpdef xcor1(int n1, double[:] t1, double[:] y1, int n2, double[:] t2, double[:] y2, double tlagmin,\
double tlagmax, double tunit, double[:] ccf, double[:] taulist, double[:] npts):
    '''
    Calculate cross-correlation function for unevenly 
    sampling data. Interpolate light curve 1 
    
    Inputs:
        n1, n2 -- number of data points for light curves 1, 2
        t1 -- time for light curve 1, assume increase;
        y1 -- flux for light curve 1;
        t2 -- time for light curve 2, assume increase;
        y2 -- flux for light curve 2;
        tlagmin -- minimum time lag;
        tlagmax -- maximum time lag;
        tunit -- tau step;        
    Outputs:
        ccf -- correlation coefficient;
        tlag -- time lag (t2 - t1); positive values mean second
                  light curve lags the first light curve, as per convention.
                 (edit by kate, march 2016)
        npts -- number of data points used;
    '''
    # declare some variables, like C
    cdef double safe, tau, fn, t2new, t1new
    cdef double rd1_sq, rd2_sq, y1sum, y1sqsum
    cdef double y2sum, y2sqsum, y1y2sum
    cdef double rd1, rd2, r
    cdef double tau_max
    cdef int knot, i, j, pointer
        
    safe = tunit*0.1
    tau_max = tlagmax+safe
    # now interpolate 1
    tau = tlagmin + 0.0
    pointer = 0
    while tau < tau_max:
        y1sum = 0.0
        y2sum = 0.0
        y1y2sum = 0.0
        y2sqsum = 0.0
        y1sqsum = 0.0
        knot = 0
        for i in xrange(n2):
            t1new = t2[i] - tau
            if t1new>=t1[0] and t1new<=t1[n1-1]:
                y1new = linearinterp(t1new, n1, t1, y1)
                y1sum = y1sum + y1new
                y1sqsum = y1sqsum + y1new*y1new
                y2sum = y2sum + y2[i]
                y2sqsum = y2sqsum + y2[i]*y2[i]
                y1y2sum = y1y2sum + y2[i]*y1new
                knot = knot + 1
        if knot>1:
            fn = float(knot)
            rd1_sq = fn*y2sqsum - y2sum*y2sum
            rd2_sq = fn*y1sqsum - y1sum*y1sum
            if rd1_sq>0.0:
                rd1 = rd1_sq**0.5
            else:
                rd1 = 0.0
            if rd2_sq>0.0:
                rd2 = rd2_sq**0.5
            else:
                rd2 = 0.0
            
            if rd1*rd2==0.0:
                r = 0.0
            else:
                r = (fn*y1y2sum - y2sum*y1sum)/(rd1*rd2)
            ccf[pointer] = r + 0.0
            taulist[pointer] = tau + 0.0
            npts[pointer] = knot + 0
            pointer = pointer + 1 # location of the ccf array
        tau = tau + tunit


cpdef xcor2(int n1, double[:] t1, double[:] y1, int n2, double[:] t2, double[:] y2, double tlagmin,\
double tlagmax, double tunit, double[:] ccf, double[:] taulist, double[:] npts):
    '''
    Calculate cross-correlation function for unevenly 
    sampling data. Interpolate light curve 2
    
    Inputs:
        n1, n2 -- number of data points for light curves 1, 2
        t1 -- time for light curve 1, assume increase;
        y1 -- flux for light curve 1;
        t2 -- time for light curve 2, assume increase;
        y2 -- flux for light curve 2;
        tlagmin -- minimum time lag;
        tlagmax -- maximum time lag;
        tunit -- tau step;        
    Outputs:
        ccf -- correlation coefficient;
        tlag -- time lag (t2 - t1); positive values mean second
                  light curve lags the first light curve, as per convention.
                 (edit by kate, march 2016)
        npts -- number of data points used;
    '''
    # declare some variables, like C
    cdef double safe, tau, fn, t2new, t1new
    cdef double rd1_sq, rd2_sq, y1sum, y1sqsum
    cdef double y2sum, y2sqsum, y1y2sum
    cdef double rd1, rd2, r
    cdef double tau_max
    cdef int knot, i, j, pointer
        
    safe = tunit*0.1
    tau_max = tlagmax+safe
    # now interpolate 1
    tau = tlagmin + 0.0
    pointer = 0
    while tau < tau_max:
        y1sum = 0.0
        y2sum = 0.0
        y1y2sum = 0.0
        y2sqsum = 0.0
        y1sqsum = 0.0
        knot = 0
        for i in xrange(n1):
            t2new = t1[i] + tau
            if t2new>=t2[0] and t2new<=t2[n2-1]:
                y2new = linearinterp(t2new, n2, t2, y2)
                y1sum = y1sum + y1[i]
                y1sqsum = y1sqsum + y1[i]*y1[i]
                y2sum = y2sum + y2new
                y2sqsum = y2sqsum + y2new*y2new
                y1y2sum = y1y2sum + y1[i]*y2new
                knot = knot + 1
        if knot>1:
            fn = float(knot)
            rd1_sq = fn*y2sqsum - y2sum*y2sum
            rd2_sq = fn*y1sqsum - y1sum*y1sum
            if rd1_sq>0.0:
                rd1 = rd1_sq**0.5
            else:
                rd1 = 0.0
            if rd2_sq>0.0:
                rd2 = rd2_sq**0.5
            else:
                rd2 = 0.0
            
            if rd1*rd2==0.0:
                r = 0.0
            else:
                r = (fn*y1y2sum - y2sum*y1sum)/(rd1*rd2)
            ccf[pointer] = r + 0.0
            taulist[pointer] = tau + 0.0
            npts[pointer] = knot + 0
            pointer = pointer + 1 # location of the ccf array
        tau = tau + tunit
