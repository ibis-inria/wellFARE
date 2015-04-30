"""
This modules gathers function for the processing of data prior
to the inference of growth rate, gene activity, protein concentration...
For the moment, the only feature implemented  the removal of ouliers caused by beads.
"""


import numpy as np
from ..curves import Curve

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
        return (-1)*remove_bumps((-1)*curve, side=1,
                                 percentile=percentile,
                                 goal=goal, niter=niter)

    xx, yy = curve.x, curve.y

    for i in range(niter):

        dydx = np.diff(yy)/np.diff(xx)
        ddyddx = np.diff(dydx)/np.diff(xx[:-1])
        thr = np.percentile(-ddyddx[ddyddx<0], percentile)

        if max(-ddyddx[ddyddx<0])< goal:
            break

        inds = np.nonzero(ddyddx > -thr)
        inds = [i+1 for i in inds]
        xx, yy = xx[inds],yy[inds]

    return Curve(xx,yy)



def filter_outliers(curve, percentile_above=100, percentile_below=100,
                             niter_above=1, niter_below=1,
                             goal_above=0, goal_below = 0,
                             smoothing_win = 4, nstd=1,
                             above_first = True):
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
                                            0,win=smoothing_win)
    fl = lambda x,y : abs(y-curve_s(x))<= nstd*curve_std(x)
    return curve.filter(fl)
