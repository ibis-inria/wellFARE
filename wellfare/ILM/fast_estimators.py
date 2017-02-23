"""

This module implements fast estimators for the time-profiles of
growth rate, promoter activity, and protein concentrations.

These estimators rely on a simple model in which gene expression
is modeled as a one-step process. This enables to compute the
observation matrix directly using an ad-hoc formula.

As a consequence these algorithms are faster and require less
parameters than their counterparts in module ``estimators``

Simple approximations are made to compute the observation matrix,
these are valid as long as the vector of estimation times (ttu) of
the different estimated input (growth rate, promoter actitivity,
protein concentration) has a fine time resolution. 

See also:
----------

estimators : collection of functions for the inference 

"""


from ..curves import Curve
from .methods import DEFAULT_ALPHAS, infer_control



def ilp_growth_rate(curve_volume, ttu, alphas=None, eps_L=.0001):
    """
    Returns
    --------
    
    mu, v_smoothed, model
      As described below.
      
    mu
      Vector of inferred mu.
    
    v_smoothed
      The predicted value of the observed volume at the same time
      points as the data. v_smoothed will appear smoothed compared to
      the measured volume.
    
    mod
      instance of sklearn.linear_model.RidgeCV, used for the Ridge
      regularization / cross-validation. Useful to get the value of
      the parameter alpha used etc.
    """
    
    if isinstance(curve_volume, list):

      results = [ilp_growth_rate(v, ttu,
                       alphas=alphas, eps_L=eps_L)
                 for v in curve_volume]
      return zip(*results)

    if alphas is None: alphas = DEFAULT_ALPHAS
    


    ttv = curve_volume.x
    dttu = 1.0*(ttu[1]-ttu[0])
    
    H_ic = np.ones((len(ttv),1))
    
    # dT is a Ny x Nu matrix with
    # dT[i,j] = ttv[i] - ttu[j]
    dT = np.array([ttv]).T - ttu
    
    H_u = ( np.maximum(0, np.minimum(dttu, dT))
            * curve_volume(ttu+ dttu/2))
        
    H = np.hstack([H_ic, H_u])
    
    growth_rate, v_smooth, ic, alpha, ascores = \
        infer_control(H, y= curve_volume.y, Nic= 1,
                      alphas= alphas, eps_L = eps_L)

    return ( Curve(ttu, growth_rate),
             Curve(ttv, v_smooth),
             ic, alpha, ascores )




def ilp_synthesis_rate(curve_fluo, curve_volume, ttu, degr,
                       alphas=None, eps_L=.0001):
    """

    
    dF/dt = s(t)V(t) - degr*F


    Parameters
    -----------

    curve_fluo
      A curve instance representing the (noisy) measured
      fluorescence

    curve_volume
      A curve instance representing the (noisy) measured
      volume

    ttu
      Times at which the control is


    Returns
    --------
    
    synth_rate, fluo_smoothed, ic, alpha, ascores
      As described below.
      
    synth_rate
      Vector. Inferred control.
      
    fluo_smoothed
      The predicted value of the observed data at the same time
      points as the data. y_smoothed will appear smoothed compared
      to y.
    
    mod
      instance of sklearn.linear_model.RidgeCV, used for the Ridge
      regularization / cross-validation. Useful to get the value
      of the parameter alpha used etc.
    """
    
    if isinstance(curve_fluo, list):

      results = [ilp_synthesis_rate(f, v, ttu, degr,
                       alphas=alphas, eps_L=eps_L)
                 for f, v in zip(curve_fluo, curve_volume)]
      return zip(*results)


    
    if alphas is None: alphas = DEFAULT_ALPHAS
    
    tt_fluo= curve_fluo.x
    H_ic = np.exp(-degr*tt_fluo).reshape((len(tt_fluo),1))
    model = lambda Y,t: 1 - degr*Y
    dtau = ttu[1]-ttu[0]
    m = odeint(model,0,[0,dtau]).flatten()[1]
    TT = (ttu-np.array([tt_fluo]).T)
    H_u = (m*np.exp(degr*TT)*(TT<0)) * curve_volume(ttu + dtau/2)
    
    H = np.hstack([H_ic, H_u])
    

    activity, fluo_smooth, ic, alpha, ascores = \
        infer_control(H, y= curve_fluo.y, Nic= 1, alphas= alphas,
                      eps_L = eps_L)

    return ( Curve(ttu, activity),
             Curve(tt_fluo, fluo_smooth),
             ic, alpha, ascores )



def ilp_concentration(curve_fluo, curve_volume, ttu, dR, dP,
                      alphas=None, eps_L=0.0001):
    """ Retrieves the concentration of a protein P, given
    the fluorescence of reporter R.

    
    Parameters
    -----------

    curve_fluo
      A curve instance representing the measured fluorescence
      (proportional to the quantities of reporter)

    curve_volume
      Volume of the population.

    dR
      Degradation rate of the reporter

    dP
      Degradation rate of the proteins.

    alphas
      Smoothing parameters to be tested.

    eps_L
      Negligible factor for the derivation matrix.

    """

    if isinstance(curve_fluo, list):

      results = [ilp_concentration(f, v, ttu, dR, dP,
                       alphas=alphas, eps_L=eps_L)
                 for f, v in zip(curve_fluo, curve_volume)]
      return zip(*results)




    tt = curve_fluo.x
    deltatau = ttu[1]-ttu[0]
    dT = np.array([tt]).T-ttu
    dTlz = dT >= 0 # ti-tj > 0
    dTlzsdtau = dTlz*(dT < deltatau) #  0 < ti-tj < delta_tau
    A = np.exp(dR*np.minimum(deltatau, dT)) - 1
    B = dTlz*np.exp(dT*(-dR))*(dP-dR)/dR
    
    Hu = (dTlzsdtau + A*B)*curve_volume(ttu+deltatau/2)
    
    Hic = np.array([np.exp(-dR*tt)]).reshape((len(tt),1))
    
    H = np.hstack([Hic, Hu])
    
    p_est, f_est, ic, a, ascores = infer_control(
                 H, curve_fluo.y, 1, alphas=alphas, eps_L=eps_L)
    
    return (Curve(ttu, p_est),
            Curve(tt, f_est),
            ic, a, ascores )