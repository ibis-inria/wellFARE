"""

pycurves
---------

This module implements the Curve class. Curve instances store two
arrays `curve.x` and `curve.y` which represent a function y=f(x). Many
operations on the curves which generally require a few lines of code
can be done in one line using the methods of Curve instances:

Features
~~~~~~~~~

Here is what you get with pycurves:


- The usual operations +,-,/,*,**, are supported between
  curves, and between curves and constants.
- Curve instances behave like functions (``mycurve(x)``) using
  interpolation (linear by default) between its points.
- Find the maximum, minimum, argmaximum, and other caracteristics of
  curves.
- Analysis: differentiate, integrate, or inverse a curve.
- Advanced functions: find the time-shift between two curves
  (``curve.find_shift_gradient``). Find the intersection points
  between two curves (``curve.find_intersections``). Fit a custom
  model to a curve (``curve.fit``).
- Add a point to a curve (``curve.add_point``), merge all the points
  from different curves (``curve.merge_with(other)``).
- Remove a point from the curve (``curve.pop()``), or many points
  according to some condition (``curve.filter()``)
- Operations on several curves: compute the mean and standard
  deviation, or any percentile (e.g. median) of a family of curves.
- Save and load curves or collections of curves in data files.
- Plot the curves with Matplotlib.
- Many middle-level methods to make your own functions :D

"""
from __future__ import print_function

import pandas

import pickle
from copy import deepcopy

import decorator

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit, fmin
from scipy.integrate import odeint


try:
    import matplotlib.pyplot as plt
except ImportError:
    MATPLOTLIB_DETECTED = False
else:
    MATPLOTLIB_DETECTED = True




@decorator.decorator
def outplace(f, curve, *a, **k):
    new_curve = deepcopy(curve)
    f(new_curve, *a, **k)
    return new_curve




class Curve:

    """

    Parameters
    -----------

    x
      List or *np.array* containing values of x

    y
     List or *np.array* containing values of y

    xy
      List or *np.array* of the points, i.e. of the form
      ``[[x1,y1],[x2,y2],...]``. This parameter can be provided
      instead of parameters ``x`` and ``y`` to define the curve.

    left_value
      Default value to be returned when the curve is evaluated at a
      point ``x`` such that ``x<curve.x.min()``. Default is ``np.nan``.

    right_value
      Default value to be returned when the curve is evaluated at a
      point ``x`` such that ``x>curve.x.max()``. Default is ``np.nan``.

    interpolation
      Possible values are 'linear', 'quadratic', 'cubic'.


    Examples
    ----------
    ::

        from pylab import *
        from pycurves import Curve

        xx = arange(0,2,0.1)
        yy = sin(xx)
        curve_1 = Curve(xx, yy)

        x_y = [(0,1), (1.5, 3), (5, 2)]
        curve_2 = Curve(xy = xy)

        curve_3 = curve_1 - curve_2

        curve_3.plot()
    """



    def __init__(self, x=None, y=None, xy=None,
                 left_value= np.nan, right_value= np.nan,
                 interpolation='linear'):

        if xy is not None:
            x,y= zip(*xy)
        self.x= np.array( x )
        self.y= np.array(y)
        self.left_value = left_value
        self.right_value = right_value
        self.interpolation = interpolation
        self._update_interpolator()



    def _update_interpolator(self):
        """
        Special function to update the interpolator after some
        inplace changes have been made.
        Note that by default the left and right values are np.nan
        """
        itr = interp1d(self.x, self.y, kind=self.interpolation)

        def f(x):
            if hasattr(x,'__iter__'):
                return np.array(list(map(f,x)))
            else:
                return (self.left_value if (x<self.x[0]) else (
                        self.right_value if (x>self.x[-1]) else (
                        self.y[0] if x==self.x[0] else (
                        self.y[-1] if x==self.x[-1] else
                        float(itr(x) )))))

        self.interpolator = f


    def __call__(self,tt):
        """ This method enables the curve to be called like a function
        of the time. """
        return self.interpolator(tt)



    # =================================================================
    # PROPERTIES OF THE CURVE

    def __len__(self):
        return len(self.x)

    @property
    def xy(self):
        """ Array of the points of the curves.

        Returns a np.array([[x1,y1],[x2,y2]...]). makes it practical to
        iter over the points with

        >>> for x,y in curve.xy: # for each point
            do_something(x,y)
        """


        return np.array([self.x, self.y]).T



    def max(self):
        """ Returns the maximum value of the curve.

        >>> curve.max() == curve.y.max() # equivalent
        """
        return self.y.max()



    def argmax(self):
        """ Returns the value of `x` for which `y` is maximal.

        This funciton always returns one value only. If there are
        several values of `x` which maximize `y`, only the first one is
        returned. To obtain all the values that maximize the curve,
        use ``arg``:

        >>> curve.arg( curve.max() )
        """
        return self.x[np.argmax(self.y)]



    def min(self):
        """ Returns the minimum value of the curve.

        >>> curve.min() == curve.y.min() # equivalent
        """
        return self.y.min()



    def argmin(self):
        """ Returns the value of `x` for which `y` is minimal.

        This funciton always returns one value only. If there are
        several values of `x` which minimize `y`, only the first one is
        returned. To obtain all the values that minimize the curve,
        use ``arg``:

        >>> curve.arg( curve.min() )
        """
        return self.x[np.argmin(self.y)]



    def span_x(self):
        """ Returns x.max() - x.min() """

        return self.x[-1] - self.x[0]

    def xlim(self):
        """ Returns the (x.min, x.max) of the curve """
        return (self.x[0], self.x[-1])



    def span_y(self):
        """ Returns y.max() - y.min() """

        return self.y.max() - self.y.min()



    def mean(self):
        """ Mean value of the curve.

        This value is NOT curve.y.mean(), it is the temporal mean of
        the curve of the form I/(x.max - x.min) where I is the temporal
        integral of the curve, computed with the trapezoid rule.
        """
        return self.integral().y[-1] / self.span_x()



    def inv(self, check_bijective = True):
        """ Returns the inverse curve.

        The inverse of the curve x->y = f(x) is the curve f(x)-> x
        The original curve must be bijective (i.e. always strictly
        increasing, or always strictly decreasing). By default this is
        checked before creation of the inverse curve. You can disable
        this checking with ``check_bijective=False``.

        Examples
        ---------

        >>> tt = np.arange(0,10,.5)
        >>> inv_sin = Curve(tt, np.sin(tt)).inv()
        ValueError: Curve not bijective, cannot invert.
        >>> inv_x2 = Curve(tt, tt**2).inv()
        >>> inv_x2(16)
        4

        See Also
        ---------
        curve.arg(y)

        Notes
        --------
        ``curve.arg(y)`` retrieves all the x such that f(x)=y. Does
        not require the curve to be bijective.

        """

        x, y = self.x, self.y

        if check_bij:
            arr = y[1:] > y[:-1]
            if (True in arr) and (False in arr):
                raise ValueError("Curve not bijective, cannot invert.")

        if self.y[1]>self.y[0]:
            return Curve(y, x)
        else:
            return Curve(y[::-1], x[::-1])



    def arg(self, y):
        """ Returns a np.array containing all the values x such that
        ``curve(x)=y``. Note that `x` and `y` can be any real, not
        necessarily elements of curve.x and curve.y .

        Examples
        ---------

        >>> tt = np.arange(0,10,.01)
        >>> curve_sin = Curve(tt, np.sin(tt))
        >>> curve_sin.arg(0)
        [[ 0.          3.14159267  6.28318532  9.42477796]]
        """

        xx, yy = self.x, self.y
        inds = np.nonzero( ((yy[:-1]-y)*(yy[1:]-y)) <= 0 )
        return np.array([ xx[i] +
                          (xx[i+1]-xx[i])*(y-yy[i])/(yy[i+1]-yy[i])
                          for i in inds])



    # ================================================================
    # TRANSFORMATIONS


    @outplace
    def fy(self, f, *f_args, **f_kwargs):
        self.y = f(self.y, *f_args, **f_kwargs)
        self._update_interpolator()

    @outplace
    def fx(self,f, *f_args, **f_kwargs):
        self.x = f(self.x, *f_args, **f_kwargs)
        self._update_interpolator()



    def mul_y(self, factor):
        """ Multiplies all values of curve.y by some factor.
        Creates a new curve. """
        return self.fy(lambda y: factor*y)



    def add_y(self, constant):
        """ Adds a constant to all values of curve.y.
        Creates a new curve. """
        return self.fy(lambda y:y+constant)



    def mul_x(self, factor):
        """ Multiplies all values of curve.x by some factor.
        Creates a new curve """
        return self.fx(lambda x: factor*x)



    def add_x(self, constant):
        """ Adds a constant to all values of curve.x ('x-shift').
        Creates a new curve. """
        return self.fx(lambda x: x + constant)



    def apply(self, f, *f_args, **f_kwargs):
        return f( self, *f_args, **f_kwargs)



    def delete_ind(self, ind):
        """ Returns a curve where points at indices specified by `ind`
        are removed. `ind` can be a single index, or a list or array
        of indices. """
        return Curve(np.delete(self.x, ind), np.delete(self.y, ind))



    def normalized(self, t=None, region=None):
        """ Divides the curve by some number. Creates a new curve.

        >>> c2 = c1.normalize() # Divides by the temporal mean.
        >>> c2 = c1.normalize(a,b) # Divides by the mean on [a,b].
        >>> c2 = c1.normalize(t) # Divides by c1(t).
        """

        if t is not None:
            D = self(t)
        elif region is not None:
            D = self.crop(region).mean()
        else:
            D =  self.mean()

        return self/D



    def diff_win(self,deriv=1, win=2, order=1, sym=True):
        """
        Returns the discrete derivatives, or increases, of
        the measure
        """
        res = []
        l,r = (win/2,win/2) if sym else (0,win)
        lx = len(self.x)
        for i in range(0,lx):
            x = self.x[max(0,i-l) : min(lx, i+r)] -  self.x[i]
            y = self.y[max(0,i-l) : min(lx, i+r)]
            coef = np.polyfit(x,y,order)[-(deriv+1)]
            if deriv > 1:
                coef *= 1.0/np.array(range(1,deriv)).prod()
            res.append(coef)
        return Curve( self.x, np.array(res))



    def diff_sg(self, nx, hw= 1, deriv=1, order=3):
        """ Data differentiation using the SG-filter.

        Smooth (and opt. differentiate) data with a Savitzky-Golay
        filter. The SG-filter removes high frequency noise from data,
        while preserving the original shape and features of the signal
        better than other other moving average techniques.

        The original code with documentation and example here:
        http://www.scipy.org/Cookbook/SavitzkyGolay

        Parameters
        ----------

        nx
          Number of The first point will be curve.x[0], the last
          point will be curve.x[-1] and the points inbetween will be
          equally spaced.

        hw
            Number of points in the smoothing half-window. The more,
            the smoother the final curve will be

        order
            The order of the polynomial used in the filtering.
            The less, the smoother the curve will be. 1,2,3 or 4 are
            commonly used. ``order`` must be less than `2*hw`.

        deriv
            The order of the derivative to compute (default = 1).
            ``deriv=0`` will result in a smoothed version of the
            original curve.

        Returns
        -------
        curve :
            A curve instance with the algorithm

        Examples
        ---------

        >>> # g(x) = f(x) + f'(x), where f is a Curve
        >>> g = f + f.diff_sg(nt= 100, deriv=1)

        """


        xx = np.linspace(self.x[0], self.x[-1], nx)
        yy = self(xx)
        rate=1.0/(xx[1]-xx[0])
        order_range = np.arange(0,order+1)
        hw_range = np.arange(-hw,hw+1).reshape((2*hw+1,1))
        b = np.mat(hw_range**order_range)
        rate_corr = rate ** deriv * np.math.factorial(deriv)
        m = np.linalg.pinv(b).A[deriv] * rate
        firstvals = yy[0] - np.abs( yy[1:hw+1][::-1] - yy[0] )
        lastvals = yy[-1] + np.abs(yy[-hw-1:-1][::-1] - yy[-1])
        yy = np.concatenate((firstvals, yy, lastvals))
        yy_sg= np.convolve( m[::-1], yy, mode='valid')

        return Curve(xx, yy_sg)



    def smooth_sg(self, nx, hw= 2, order=2):
        """ Smoothing using the Savitsky-Golay filter """
        return self.diff_sg(nx, hw=hw, deriv=0, order=order)


    @outplace
    def integral(self):

        dx = np.diff(self.x)
        terms = 0.5 * (self.y[1:]+self.y[:-1]) *dx
        self.y = np.hstack([[0], np.cumsum(terms)])
        self.left_value = 0
        self.right_value = self.y[-1]




    # ================================================================
    # OPERATORS



    def _operate(self, operator, other, x=None):
        """ General operations on two curves.

        This is a function for applying an operator (addition,
        multiplication...)
        """

        if isinstance(other, Curve):
            self2, other2 = Curve._same_x([self, other], x=x)
            return Curve(self2.x, operator(self2.y, other2.y))
        else:
            return self.fy( lambda y: operator(y, other))



    def __add__(self, term):
        """ Adds, 2 measures, or one measure and a constant. """
        return self._operate(lambda y1, y2: y1+y2, term)



    def __sub__(self, term):
        """ Substracts, 2 measures, or one measure and a constant. """
        return self._operate(lambda y1, y2: y1-y2, term)



    def __truediv__(self, term):
        """ Divides, 2 measures, or one measure and a constant. """
        return self._operate(lambda y1, y2: y1/y2, term)



    def __mul__(self, term):
        """ Multiplies, 2 measures, or one measure and a constant. """
        return self._operate(lambda y1, y2: y1*y2, term)



    def __pow__(self, term):
        """ Elevates to the power 'term': y=f(x)^term """
        return self._operate(lambda y1, y2: y1**y2, term)



    def __radd__(self, term):
        return self+term



    def __rsub__(self, term):
        return (self - term)*(-1)



    def __rmul__(self, term):
        return self*term



    def __rdiv__(self, term):
        return (self**(-1.0)) * term



    # =================================================================
    # UTILITIES



    def fit(self, f, **kwargs):
        """ Fits a function to the Curve with Scipy's curve_fit.
            curve.fit(lambda x, c: c[0]*x + c[1])"""
        return curve_fit(f, self.x, self.y, **kwargs)



    @staticmethod
    def odeint(model, y0, x, args=(), **odeint_kw):
        """ Solves an ODE with Scipy's odeint, returns Curve(s).


        Parameters
        -----------

        model
          A function of the form ``Y,x -> (dY/dx)(x)`` or optionally
          with parameters ``Y,x,*args -> (dY/dx)(x)``.

        y0
          Initial conditions

        x
          Array of the values of x in which the model's solution
          must be evaluated.

        args
          Parameters to provide to the model (if the model accepts
          parameters). Must always be a tuple, even if there is
          only one parameter

        Examples
        ---------

        >>> import numpy as np
        >>> # Y'(t) = aY(t) + sin(t)
        >>> model = lambda Y,t : a*Y+sin(t)
        >>> tt =  np.arange(0,10,.01)
        >>> # case: Y has one variable only
        >>> curve_Y = Curve.odeint(model, 0, tt)
        >>> # case: Y has two variables
        >>> curve_Y1, curve_Y2 = Curve.odeint(model, [0,0], tt)
        """

        result = odeint(model, y0, x, args=args, **odeint_kw)

        n_x, n_vars = result.shape

        if n_vars == 1:
            return Curve(x, result.flatten())
        else:
            return [Curve(x, line) for line in result.T]




    # =================================================================
    # ADDING / REMOVING POINTS



    @outplace
    def add_point(self, x, y):
        xx =  np.hstack([self.x, np.array([x])])
        yy = np.hstack([self.y, np.array([y])])
        inds = np.argsort(xx)
        self.x = xx[inds]
        self.y = yy[inds]
        self._update_interpolator()


    @outplace
    def filter(self, cond):
        xx, yy = self.x, self.y
        keep = [i for i in range(len(self.x))
                   if cond(xx[i], yy[i])]
        self.x = xx[keep]
        self.y = yy[keep]
        self._update_interpolator()


    @outplace
    def pop(self,i=-1):
        L = len(self.x)
        if i<0: i = L-i
        indices = [ind for ind in range(L) if ind!=i]
        self.x = self.x[indices]
        self.y = self.y[indices]
        self._update_interpolator()


    @outplace
    def crop(self, xmin=None,xmax=None):
        """ Removes all measurements taken before xmin or after xmax.
        """

        xx = self.x
        if (xmin is None) or (xmin <= xx[0]):
            xmin= xx[0]

        if (xmax is None) or (xmax >= xx[-1]):
            xmax= xx[-1]

        x2 = xx[(xx >= xmin) & (xx <= xmax)]
        y2 = self(x2)

        new_x = np.hstack([[xmin],x2,[xmax]])
        new_y = np.hstack([[self(xmin)],y2,[self(xmax)]])
        self.x =  new_x
        self.y =  new_y

        self._update_interpolator()


    def resample(self, new_x, keep_xspan=True):
        """ Returns a new curve resampled at the times x """

        new_x = np.array(new_x)
        if keep_xspan:
            xmin, xmax = self.xlim()
            new_x = new_x[(new_x >= xmin) & (new_x <= xmax)]
        return Curve(new_x, self(new_x))

    #=================================================================
    # OPERATIONS ON SEVERAL CURVES



    @staticmethod
    def _same_x(curves, x=None):
        """
        Returns the interpolations of different curves on the times
        of all the curves combined.
        """

        xmin = max([c.x[0] for c in curves])
        xmax = min([c.x[-1] for c in curves])
        if x is None:
            x = sum([list(c.x) for c in curves], [])
            x = np.sort(np.array(list(set(x))))
        x = x[(x >= xmin) & (x <= xmax)]

        return [ c.resample(x) for c in curves]


    def corrcoef(self, other, time_aware=False):
        if time_aware:
            raise NotImplementedError()
            self_, other_ = Curve._same_x([self, other])
        else:
            other = other.resample(self.x)
            return np.corrcoef(self(other.x), other.y)

    @staticmethod
    def merge(curves):
        """
        merge([c1, c2, c3...]) returns one Curve containing all
        the points of the given measures.
        """
        xx = np.hstack([c.x for c in curves])
        yy = np.hstack([c.y for c in curves])
        xx_yy = np.vstack([xx, yy])
        xx_yy = xx_yy[:,xx_yy[0,:].argsort()]
        return Curve(xx_yy[0,:], xx_yy[1,:])



    def merge_with(self,curves):
        """
        c1.merge_with(c2) or c1.merge_with(c2,c3,c4)
        """
        if isinstance(curves, Curve):
            curves = [curves]
        return Curve.merge([self] + curves)



    @staticmethod
    def mean_std(curves, x=None):
        """ Mean and standard deviation of several curves.

        Parameters
        ----------

        curves
          A list of curve instances.

        x: np.array
          Times at which the final curve should be computed. If ``x``
          is None, the x of the first Curve in the list is used.

        Returns
        --------

        curve_mean, curve_std_dev
          Two Curves instances representing respectively the mean and
          the standard deviations of the family of curves, at the
          specified points ``x``.

        Examples
        ---------

        >>> # lets plot the mean and std deviations of 3 curves
        >>> means, stds = mean_std([cur1,cur2,cur3])
        >>> errorbars(cur1.x, means.y, stds.y) # matplotlib function
        """

        if x is None:
            curves = Curve._same_x(curves)
            x = curves[0].x
        else:
            curves = [c(x) for c in curves]

        yy = np.array([c.y for c in curves])
        means = yy.mean(axis=0)
        stds = yy.std(axis=0)
        return Curve(x,means), Curve(x,stds)



    @staticmethod
    def percentile(curves, p, x=None):
        """
        Returns a Curve containing the p-percentile of the curves,
        at different values x. If x is None, the x of the
        first Curve is used.
        """

        if x is None:
            x = curves[0].x
        yy_list = [np.array([c(px) for c in curves]) for px in x]
        percentiles = [np.percentile(yy,p) for yy in yy_list]
        return Curve(x,percentiles)

    @staticmethod
    def minimum(curves, x=None):
        """
        Returns a Curve containing the minimum of the curves,
        at different values x. If x is None, the x of the
        first Curve is used.
        """

        if x is None:
            x = curves[0].x
        yy = [min([c(px) for c in curves]) for px in x]
        return Curve(x,yy)

    @staticmethod
    def maximum(curves, x=None):
        """
        Returns a Curve containing the minimum of the curves,
        at different values x. If x is None, the x of the
        first Curve is used.
        """

        if x is None:
            x = curves[0].x
        yy = [max([c(px) for c in curves]) for px in x]
        return Curve(x,yy)



    def find_shift_fft(self, other, precision=None, cutoff_x=None):
        """ Shift-finding using cross-validation and the FFT


        Finds the x-shift to apply to the other curves so that they
        will be synchronized to the current curve. Uses the Fast
        Fourier Transform, and will be generally faster but less
        precise than Curve.find_shift_gradient.

        **WARNING** : this algorithm requires the signal to be
        wave-like, i.e. the beginning and the end should look alike,
        otherwise it will give meaningless results.

        Parameters
        -----------

        other
          A list of other Curves instances whose shifts with respect
          to the current curve are to be estimated.

        precision
          The wanted precision in the x-shift. Note that a curve with
          a x-span of Lx will be evaluated in Lx/precision points,
          so don't set it too low :). Default for precision will be
          Lx/Nx where Nx is the number of points in the curve.

        cutoff_x
          The elements of the curves spectrum corresponding to a
          period larger than cutoff_x will not be taken into account
          in the shift estimation. Therefore cutoff_x can be used to
          denoise or simplify the curves before synchronization, with
          no computational cost at all


        See Also
        ---------

        Curve.find_shift_gradient

        Examples
        ---------

        shifts = ref_curve.find_shift_fft(other_curves)
        synchronized_curves = [curve.add_x(shift)
                    for curve,shift in zip(other_curves, shifts)]
        """

        curves = [self]+other
        xstart = max(c.x[0] for c in curves)
        xend = min(c.x[-1] for c in curves)
        Lx = xend-xstart
        dx = 1.0*Lx/len(self.x) if precision is None else precision

        xx = np.arange(xstart, xend, dx)
        signals = [c(xx) for c in curves]


        rffts = [np.fft.rfft(s) for s in signals]

        Nx = len(signals[0])

        if cutoff_x is not None:
            # 'filter' the curves by cutting off high frequencies.

            cut = int(1.0*Lx/cutoff_x)
            if cut < Nx:
                for tf in rffts:
                    tf[ cut:] = 0

        ref_fft = rffts[0]
        argmaxs = [ np.argmax(np.fft.irfft(ref_fft * np.conj(f)))
                    for f in rffts[1:]]

        shifts = np.array([ ((a-Nx) if (a > Nx/2) else a)
                    for a in argmaxs ])

        return dx*shifts



    def find_shift_gradient(self, others, tt, shifts0=[0], **fmin_kw):
        """ Shift-finding using a gradient algorithm.

        Parameters
        ===========

        others
          A list of curves to be synchronized with respect to the
          current curve.


        tt
          time points used to compute the cross-correlation of the
          curves

        shifts0
          First shifts tried. The gradient algorithm is started from
          one of these shifts.

        fmin_kw
          keyword arguments for scipy.optimize.fmin

        Returns
        ========

        A list of shifts to apply to the curves in the list ``others``
        so that they will be maximally correlated (=synchronized)
        with the current curve.


        """


        self_y = self(tt)
        self_y = (self_y - self_y.mean())


        def make_obj(curve):

            def obj(s):
                curve_y = curve(tt-s)
                curve_y = curve_y - self_y.mean()
                return -(self_y*curve_y).sum() / curve_y.std()

            return obj

        def optimize(obj):

            shift0 = shifts0[np.argmin(list(map(obj, shifts0)))]
            return fmin(obj, x0=shift0, disp=0, **fmin_kw)[0]

        return [optimize( make_obj( c )) for c in others]



    def find_intersections(self, other):
        """ Finds all the `x` at which the curve intersects with the
        `other` curve. """

        return (self-other).arg(0)



    #=================================================================
    # I/O OPERATIONS



    def __getstate__(self):
        return {'x':self.x, 'y': self.y}



    def __setstate__(self, state):
        self.__init__(state['x'], state['y'])



    def save(self, filename=None):
        """ Saves the curve to a file.

        Examples
        ---------

        >>> curve.save("saved_curve.dat")
        >>> # Later...
        >>> curve = Curve.load("saved_curve.dat")
        """

        with open(filename, 'w+') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)



    def savetxt(self, filename, fmt='%.18e', delimiter=' ',
                newline='\n', header='', footer='', comments='# '):
        """ Writes the curve to a txt/csv file.

        Ideal for interfacing with other, less cool languages such as
        Matlab or PHP.

        See numpy.savetxt for the meaning of the different arguments.

        """

        np.savetxt(filename, X=self.xy, fmt=fmt, delimiter=delimiter,
                     newline=newline, header=header, footer=footer,
                     comments=comments)



    @staticmethod
    def save_object(obj, filename):
        """ Saves an object containing several curves.

        This function enables to save many curves at once into one
        file. The ``obj`` is any object (dictionnary, list) containing
        one or several curves.

        Example
        --------

        >>> curves_set = {'variable_1': curve_1,
                          'variable_2': curve_2,
                          'variable_3': curve_3}
        >>> curve_list = [curve_1, curve_2, curve_3]
        >>> Curve.save_object(curves_set, "curves_set.dat")
        >>> Curve.save_object(curves_list, "curves_list.dat")
        >>> # Later...
        >>> curves_set = Curves.load("curves_set.dat")
        >>> curves_list = Curves.load("curves_list.dat")
        """

        with open(filename, 'w+') as output:
            pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)



    @staticmethod
    def load(filename):
        """ Loads a curve or object, saved with ``curve.save()`` or
        with ``Curve.save_object()``.

        See these functions for details.
        """

        with open(filename, 'rb') as f:
            r = pickle.load(f)
        return r

    @staticmethod
    def loadtxt(filename, comments='#', delimiter=None,
                converters=None, skiprows=0, usecols=None,
                unpack=False, ndmin=0):
        """ Reads a curve from a txt/csv file.

        Ideal for interfacing with other, less cool languages such as
        Matlab or PHP.

        See numpy.loadtxt for the meaning of the different arguments.

        """

        xy = np.loadtxt(filename, comments=comments,
               delimiter=delimiter, converters=converters,
               skiprows=skiprows, usecols=usecols, unpack=unpack,
               ndmin=ndmin)

        return Curve(xy = xy)

    # ===   PRINTING AND PLOTTING ====================================



    def __str__(self):
        return str(self.xy)



    def plot(self,ax = None, hold = False, **plot_kw):
        """ Plots the Curve using Matplotlib

        ``curve.plot()`` is equivalent to (but shorter)

        >>> import matplotlib.pyplot as plt
        >>> plt.plot(curve.x, curve.y)

        Parameters
        -----------

        ax
          The Matplotlib `ax` on which to draw the curve

        hold
          Should be True in non-interactive sessions in order to be
          able to actually see the curve. Unnecessary in the IPython
          Notebook.

        **plot_kw
          Any of the (very many) options of Matplotlib's plot method:
          'color', 'marker', 'label' (for the legend),
          'lw' (linewidth), 'ls' (linestyle), etc. See Matplotlib's
          doc for more details.


        Examples
        ---------

        >>> fig, ax = plt.subplots(1,2)
        >>> curve.plot(ax[1], color='blue', marker='o', label='curve')

        """

        if not MATPLOTLIB_DETECTED:
            raise ImportError("You need to install Matplotlib in"
                               "order to plot curves.")

        if ax is None:
            fig,ax = plt.subplots(1)
        line = ax.plot(self.x, self.y, **plot_kw)

        if hold:
            plt.show()

        return line

    def hist(self,ax = None, hold = False, **plot_kw):
        """ Plots the histogram of the Curve using Matplotlib

        ``curve.hist()`` is equivalent to (but shorter)

        >>> import matplotlib.pyplot as plt
        >>> plt.hist(curve.y)

        Parameters
        -----------

        ax
          The Matplotlib `ax` on which to draw the curve

        hold
          Should be True in non-interactive sessions in order to be
          able to actually see the curve. Unnecessary in the IPython
          Notebook.

        **plot_kw
          Any of the (very many) options of Matplotlib's plot method:
          'color', 'marker', 'label' (for the legend),
          'lw' (linewidth), 'ls' (linestyle), etc. See Matplotlib's
          doc for more details.


        Examples
        ---------

        >>> fig, ax = plt.subplots(1,2)
        >>> curve.hist(ax[1], color='blue', bins=50, label='curve')

        """

        if not MATPLOTLIB_DETECTED:
            raise ImportError("You need to install Matplotlib in"
                               "order to plot curves.")

        if ax is None:
            fig,ax = plt.subplots(1)
        line = ax.hist(self.y, **plot_kw)

        if hold:
            plt.show()

        return line


# ====================================================================


if __name__ == '__main__':


    import numpy as np
    xx = np.linspace(0,100, 200)
    yy = np.sin(xx) + xx**2
    c =  Curve(xx, yy)
    print(c(0))
    print(c(50))

    c2 = c.crop(50)
    print(c2.x[0], c2.y[0])
    print("in 50", c2(50))
    print(c2(51))
    """
    from pylab import *

    x = arange(0,2,0.05)
    curve_1 = Curve(x, y=x**2 + 3*x)
    print curve_1.crop(0,1).x
    curve_2 = curve_1.diff_win()
    curve_3 = curve_1.fy( sin )
    curve_4 = (curve_1 + curve_3) / curve_2

    curve_4 = Curve.load('../../../test.dat')
    fig,ax = subplots(1)
    curve_3.plot(ax, label = 'curve 3')
    curve_4.plot(ax, label = 'curve 4')
    ax.legend()
    curves = [ Curve(x = x+1.0*i/100,
                     y = sin(x)+np.random.normal(0,5e-2,len(x)) )
               for i in range(8) ]
    fig, ax  = subplots(1)
    for c in curves:
        c.plot(ax, lw=0.5, alpha=0.65)
    mean, std = Curve.mean_std(curves)
    mean.plot(ax, lw=2, c='k')
    for c in (mean+std, mean-std):
        c.plot(ax, lw=2, c='k', ls='--')

    plt.show()
    """
