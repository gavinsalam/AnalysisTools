#!/usr/bin/env python3
""" module poly_fit.py/

given N pairs (xi,yi) of points, compute the coefficients of the
polynomial of degree N-1 going through the N points
"""

import numpy as np
import sys

        
#------------------------------------------------------------------------
def get_polynomial_coefficients(x, y, dy=None):
    """ 
given N points (xi,yi), compute the coefficients of the
polynomial of degree N-1 going through the N points

A vector (array) of the coefficients (starting with the constant term)
is returned [i.e. the ith coefficient is the coefficient of x^i

Note that the xi, yi can be multi-dimensional themselves. If x and y
are, sayt N-dimensional, then the first index is used to label the
series of points one needs to interpolate over and the interpolation is
done for each of the points in the remaining N-1-dimensional space.

In practice, if you have a vector of x values
  x = [x0, x1, ... x_{N-1}]
and for each xi an array of y values y0
  y = [y0, ... y_{N-1}]
where each of the y can be arrays themselves, calling
  get_polynomial_coefficients(x, y)
will give you the coefficients a0, ... a_{N-1} so that
  y_i = a0 + a1 xi + ... + a_{N-1} xi^{N-1}
If y is an arry, each of the ai's will be arrays as well.

If an uncertainty dy on each is specified, the function returns the a's
as well as their uncertainty.  

"""
    
    if isinstance(x,np.ndarray) and (x.shape[0] != y.shape[0] or (dy is not None and x.shape[0] != dy.shape[0])):
        print ('x, y and dy arrays must have the same shape', file=sys.stderr)
        return None

    n = x.shape[0]

    # case of a single point (-> constant polynomial)
    if n==1:
        return y.copy(), dy.copy()

    # case of 2 points (-> degree 1 polynomial)
    if n==2:
        d = x[1]-x[0]
        res = np.asarray([(y[0]*x[1]-y[1]*x[0])/d , (y[1]-y[0])/d])
        if dy is None: return res
        return res,\
            np.asarray([np.sqrt(np.square(dy[0]*x[1])+np.square(dy[1]*x[0]))/d , np.sqrt(np.square(dy[1])+np.square(dy[0]))/d]),

    # case of 3 points (-> degree 2 polynomial)
    t = np.copy(y)
    dt = np.copy(dy)
    for i in range(n):
        for j in range(n):
            if i==j: continue
            t[i]  =  t[i]/(x[j]-x[i])
            if dy is not None:
                dt[i] = dt[i]/(x[j]-x[i])

    if n==3:
        res = np.asarray([x[1]*x[2]*t[0] + x[0]*x[2]*t[1] + x[0]*x[1]*t[2],
                          -(x[1]+x[2])*t[0] - (x[0]+x[2])*t[1] - (x[0]+x[1])*t[2],
                          t[0] + t[1] + t[2]])
        if dy is None: return res
        return res,\
            np.asarray([np.sqrt(np.square(x[1]*x[2]*dt[0]) + np.square(x[0]*x[2]*dt[1]) + np.square(x[0]*x[1]*dt[2])),
                        np.sqrt(np.square((x[1]+x[2])*dt[0])+np.square((x[0]+x[2])*dt[1])+np.square((x[0]+x[1])*dt[2])),
                        np.sqrt(np.square(dt[0])+np.square(dt[1])+np.square(dt[2]))])
    
    # case of 4 points (-> degree 3 polynomial)
    if n==4:
        res = np.asarray([x[1]*x[2]*x[3]*t[0] + x[0]*x[2]*x[3]*t[1] + x[0]*x[1]*x[3]*t[2] + x[0]*x[1]*x[2]*t[3],
                          -(x[1]*x[2]+x[1]*x[3]+x[2]*x[3])*t[0] - (x[0]*x[2]+x[0]*x[3]+x[2]*x[3])*t[1] - (x[0]*x[1]+x[0]*x[3]+x[1]*x[3])*t[2] - (x[0]*x[1]+x[0]*x[2]+x[1]*x[2])*t[3],
                          (x[1]+x[2]+x[3])*t[0] + (x[0]+x[2]+x[3])*t[1] + (x[0]+x[1]+x[3])*t[2]+ (x[0]+x[1]+x[2])*t[3],
                          -(t[0] + t[1] + t[2]+ t[3])])
        if dy is None: return res
        return res,\
            np.asarray([np.sqrt(np.square(x[1]*x[2]*x[3]*dt[0]) + np.square(x[0]*x[2]*x[3]*dt[1]) + np.square(x[0]*x[1]*x[3]*dt[2]) + np.square(x[0]*x[1]*x[2]*dt[3])),
                        np.sqrt(np.square((x[1]*x[2]+x[1]*x[3]+x[2]*x[3])*dt[0]) + np.square((x[0]*x[2]+x[0]*x[3]+x[2]*x[3])*dt[1]) + np.square((x[0]*x[1]+x[0]*x[3]+x[1]*x[3])*dt[2]) + np.square((x[0]*x[1]+x[0]*x[2]+x[1]*x[2])*dt[3])),
                        np.sqrt(np.square((x[1]+x[2]+x[3])*dt[0]) + np.square((x[0]+x[2]+x[3])*dt[1]) + np.square((x[0]+x[1]+x[3])*dt[2])+ np.square((x[0]+x[1]+x[2])*dt[3])),
                        np.sqrt(np.square(dt[0]) + np.square(dt[1]) + np.square(dt[2])+ np.square(dt[3]))])

    # larger multiplicities are current;y not supported
    print ('only implemented up to 4 points so far', file=sys.stderr)
    return None

#------------------------------------------------------------------------
def get_poly_coeffs_scalar(x,y):
    """
    Simpler version of get_polynomial_coefficients for situation when
    the y values are a simple array of scalar objects rather than
    numpy arrays (as happens, for example, when using hfile.ValueAndError
    objects)
    """
    n = len(x)
    if (n != len(y)): raise ValueError("x and y lengths are not equal")

    if n == 1:
        return y

    if n == 2:
        d = x[1]-x[0]
        return [(y[0]*x[1]-y[1]*x[0])/d , (y[1]-y[0])/d]

    t = y * 1
    for i in range(n):
        for j in range(n):
            if i==j: continue
            t[i]  =  t[i]/(x[j]-x[i])

    if n == 3:
        res = [x[1]*x[2]*t[0] + x[0]*x[2]*t[1] + x[0]*x[1]*t[2],
                -(x[1]+x[2])*t[0] - (x[0]+x[2])*t[1] - (x[0]+x[1])*t[2],
                  t[0] + t[1] + t[2]]
        return res

#------------------------------------------------------------------------    
def apply_polynomial(coefs, x):
    """ 
return P(x) for a polynomial of given coefficients (starting w the constant term)
"""
    x0=1.0
    res=0.0
    for coef in coefs:
        res = res+coef*x0
        x0 = x0*x

    return res

#------------------------------------------------------------------------
def get_extrapolation_to_zero(x, y, dy=None):
    if dy is None:
        a = get_polynomial_coefficients(x,y)
        return a[0]
    a,da = get_polynomial_coefficients(x,y,dy)
    return a[0], da[0]

#------------------------------------------------------------------------
if __name__ == "__main__":
    # there are several formats supported by this script
    #
    # - the default one has lines of the form
    #      x  y dy
    #   extrapolating to x=0
    #   The output is then
    #      y(x->0)  dy(x->0)
    #   if only x and y are specified, no uncertaities are produced
    #   
    # - alternatively, all x values can be given on the 1st line
    #   starting the line by "x:"
    #      x: x1 x2 ... xn
    #   The following lines can be of one of the following 4 formats:
    #      format 1: y1 y2 ... yn
    #      format 2: y1 dy1 y2 dy2 ... yn dyn
    #      format 3: t y1 y2 ... yn
    #      format 4: t y1 dy1 y2 dy2 ... yn dyn
    #   In the last 2 cases, "t" is an auxiliary variable (useful if
    #   we want e.g. to extrapolate a full distribution).
    #   The variants w dy's include uncertainties.
    #   The output is as follows:
    #      format 1: y(x->0)
    #      format 2: y(x->0) dy(x->0)
    #      format 3: t y(x->0)
    #      format 4: t y(x->0) dy(x->0)
    #   Each input line is extrpaolated independently

    has_x_at_top = False

    xs=[]
    ts=[]
    ys=[]
    dys=[]

    print_all=False
    
    for l in sys.stdin:
        l=l.rstrip()
        if len(l)>0:
            # if a line says "all" we print all the coefficients (and their potentiak uncertainty)
            # [does not work for "has_x_at_top]
            if l=='all':
                print_all=True
                continue
                
            toks=l.split()
            if toks[0]=='x:':
                xs = np.asarray([float(t) for t in toks[1:]])
                has_x_at_top = True
                n = xs.shape[0]
            else:
                if has_x_at_top:
                    # format depends on number of tokens
                    if len(toks)==n:  # format 1
                        ys.append(np.asarray([float(t) for t in toks]))
                    elif len(toks)==2*n:  # format 2
                        ys .append(np.asarray([float(t) for t in toks[::2]]))
                        dys.append(np.asarray([float(t) for t in toks[1::2]]))
                    elif len(toks)==n+1:  # format 3
                        ts.append(float(toks[0]))
                        ys.append(np.asarray([float(t) for t in toks[1:]]))
                    elif len(toks)==2*n+1:  # format 4
                        ts.append(float(toks[0]))
                        ys .append(np.asarray([float(t) for t in toks[1::2]]))
                        dys.append(np.asarray([float(t) for t in toks[2::2]]))
                else:
                    # read x y dy
                    xs.append(float(toks[0]))
                    ys.append(float(toks[1]))
                    if len(toks)>2:
                        dys.append(float(toks[2]))
                        
    # now that we've read the inputs, process them
    if has_x_at_top:
        # ys and dys are in a format where each "row" (the 1st index)
        # is a value of "t". We want the opposite.
        ys  = np.transpose(np.asarray(ys))
        if len(dys)>0:
            dys = np.transpose(np.asarray(dys))
            a, da = get_extrapolation_to_zero(xs, ys, dys)
            if len(ts)>0:
                for t,y,dy in zip(ts, a, da):
                    print (t, y, dy)
            else:
                for y,dy in zip(a, da):
                    print (y, dy)
        else:
            a = get_extrapolation_to_zero(xs, ys)
            if len(ts)>0:
                for t,y in zip(ts, a):
                    print (t, y)
            else:
                for y in a:
                    print (y)
    else:
        if len(dys)>0:
            a, da = get_polynomial_coefficients(np.asarray(xs), np.asarray(ys), np.asarray(dys))
            if print_all:
                for i in range(len(a)): print (a[i],da[i],end=' ')
                print ('')
            else:
                print (a[0], da[0])
        else:
            a = get_polynomial_coefficients(np.asarray(xs), np.asarray(ys))
            if print_all:
                for i in range(len(a)): print (a[i],end=' ')
                print ('')
            else:
                print (a[0])
            
    
    
