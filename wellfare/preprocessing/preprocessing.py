"""
This modules gathers function for the processing of data prior
to the inference of growth rate, gene activity, protein concentration...
For the moment, the only feature implemented  the removal of ouliers caused by beads.
"""


import numpy as np
from ..curves import Curve
from scipy.interpolate import UnivariateSpline


def remove_bumps(curve, side, percentile=50, niter=1, goal=0):
    """ Removes onesided outliers

    See also
    ----------

    filter_outliers

    Parameters
    -----------

    curve
       the Curve instance that will be filtered

    side
      +1 (default) for removing outliers above the curve, -1 to remove
      outliers below the curve


    percentile
       If side==1, of all the points with negative second derivative,
       only the percentile with higher values will be conserved, at
       each iteration.

    niter
      Maximal number of iterations (or actual number of iterations if
      goal==0)

    goal
      The iterations will stop if the lowest second derivative in the
      curve is larger than -goal.


    Returns
    --------

    filtered_curve
      A Curve instance without the outliers.


    """

    if side == -1:
        return (-1) * remove_bumps((-1) * curve, side=1,
                                   percentile=percentile,
                                   goal=goal, niter=niter)

    xx, yy = curve.x, curve.y

    for i in range(niter):

        dydx = np.diff(yy) / np.diff(xx)
        ddyddx = np.diff(dydx) / np.diff(xx[:-1])

        thr = np.percentile(-ddyddx[ddyddx < 0], percentile)
        inds = np.nonzero(ddyddx > -thr)

        if max(-ddyddx[ddyddx < 0]) < goal:
            break
        inds = [i + 1 for i in inds]
        xx, yy = xx[inds], yy[inds]

    return Curve(xx, yy)


def filter_outliers(curve, percentile_above=100, percentile_below=100,
                    niter_above=1, niter_below=1,
                    goal_above=0, goal_below=0,
                    smoothing_win=4, nstd=1,
                    above_first=True):
    """ Removes the ouliers in a curve.

    This is done in 4 steps: first, bumps are removed

    """

    curve_f = curve
    if not above_first:
        curve_f = remove_bumps(curve_f, side=-1,
                               percentile=percentile_below,
                               niter=niter_below, goal=goal_below)

    curve_f = remove_bumps(curve_f, side=1,
                           percentile=percentile_above,
                           niter=niter_above, goal=goal_above)

    if above_first:
        curve_f = remove_bumps(curve_f, side=-1,
                               percentile=percentile_below,
                               niter=niter_below, goal=goal_below)

    curve_s = curve_f.diff_win(0, win=smoothing_win)

    curve_std = (curve_f - curve_s).fy(abs).diff_win(
        0, win=smoothing_win)

    def fl(x, y): return abs(y - curve_s(x)) <= nstd * curve_std(x)

    return curve.filter(fl)


def calibration_curve(abscurve, fluocurve, smoothing, extrapoldistance, validinterval):
    """ Return raw and smoothed (polyfited) calibration curve for a background well. (Fluo = f(Abs))

    """

    absinterp = abscurve(fluocurve.x)
    autofluo = fluocurve(fluocurve.x)

    # sort and remove duplicates
    absinterp, sortedbyabs = np.unique(absinterp, return_index=True)
    autofluo = autofluo[sortedbyabs]

    # remove nan
    nonanargs = np.logical_not(np.logical_or(
        np.isnan(absinterp), np.isnan(autofluo)))
    absinterp = absinterp[nonanargs]
    autofluo = autofluo[nonanargs]

    rawcalibrationcurve = Curve(absinterp, autofluo)

    # keep valid interval data
    if (validinterval == None) or (not isinstance(validinterval, list)):
        validinterval = [None, None]
    if validinterval[0] == None:
        validinterval[0] = np.min(absinterp)
    if validinterval[1] == None:
        validinterval[1] = np.max(absinterp)

    validindexes = np.where(np.logical_and(
        absinterp >= validinterval[0], absinterp <= validinterval[1]))
    absinterp = absinterp[validindexes]
    autofluo = autofluo[validindexes]

    # smoothing
    if smoothing == None:
        smoothing = 10
    if smoothing > absinterp.shape[0]:
        print("smothing parameter too big")
        smoothing = absinterp.shape[0]

    calibrationcurve = Curve(absinterp, autofluo)
    smoothedcalcurve = calibrationcurve.smooth_sg(1000, smoothing, 3)

    # extrapolate
    spline = UnivariateSpline(smoothedcalcurve.x, smoothedcalcurve.y, k=1)
    #extrapoldistance = (absinterp.max()-absinterp.min())/10
    extendedabs = np.linspace(absinterp.min() - extrapoldistance,
                              absinterp.max() + extrapoldistance, 1000)

    return [rawcalibrationcurve, Curve(extendedabs, spline(extendedabs))]


def filter_outliersnew(curve, outliersnew_smoothing=20, outliersnew_nstd=2, outliersnew_iterations=1, outliersnew_uppenalty=1, outliersnew_downpenalty=1):
    """ Return outlier free curve using smoothed curve.
    """
    if outliersnew_uppenalty == None:
        outliersnew_uppenalty = 1
    if outliersnew_downpenalty == None:
        outliersnew_downpenalty = 1

    for i in range(outliersnew_iterations):
        # smoothing
        curvelength = curve.x.shape[0]
        if outliersnew_smoothing > curvelength:
            outliersnew_smoothing = curvelength - 1

        smoothedcurve = curve.smooth_sg(curvelength, outliersnew_smoothing, 3)

        std_cutoff = outliersnew_nstd * (curve - smoothedcurve).y.std()
        #percent_to_remove = np.percentile(np.abs((curve - smoothedcurve).y),100 - outliersnew_uppenalty)
        #std_cutoff = min(outliersnew_nstd*curve_std,percent_to_remove)
        # filter using nstd*std

        def fl(x, y):
            dist = y - smoothedcurve(x)
            if dist > 0:
                # Normalize cutoff with opposite penalty value
                return dist <= std_cutoff * outliersnew_downpenalty
            else:
                # Normalize cutoff with opposite penalty value
                return -dist <= std_cutoff * outliersnew_uppenalty
        #fl = lambda x,y : abs(y-smoothedcurve(x))<= std_cutoff

        curve = curve.filter(fl)
        if curvelength == curve.x.shape[0]:
            break

    return curve
