"""

This module implements estimators for the promoter activity and
protein concentrations, in which the expression of the reporter and
the gene of interest are modeled as multistep processes.

"""



import numpy as np
import scipy

from ..curves import Curve
from .methods import DEFAULT_ALPHAS, infer_control
from scipy.integrate import odeint

# ==== Kinetic models =============================================

# for the number of reporter protein and protein of interest.




# Observation "matrices"
obs_R = np.array([0, 0, 1])
obs_P = np.array([0, 1])



# === Helper functions ==========================================


def make_H(model, obs, ttu, tty):
    """ Fast method to make a "production" matrix H, representing
    the dependency of y_obs on the input u of a given model:
    
    H u = y_obs, where 
    dY/dt = model( Y, t, u ) 
    yobs = Y.dot(obs)
    """
    
    nvars = len(obs)
    
    def canon(i):
        a = np.zeros(nvars)
        a[i] = 1.0
        return a
    
    
    unull = lambda t:0
    tt = tty-ttu[0]
    Hic = np.vstack([odeint(model, canon(i), tt, args=(unull,)).dot(obs)
                     for i in range(nvars)] ).T
    
    dT = np.array([tty]).T - ttu
    
    dTf = dT.flatten()
    
    tt_sorted  = np.sort(list( set( dTf[dTf>=0] ) ))
    
    u0 = Curve( [ttu[0], ttu[1]] , [1.0, 1.0] )
    u0.left_value = u0.right_value = 0.0
    yy = odeint(model, nvars*[0.0], tt_sorted, args=(u0,)).dot(obs)
    
    D = dict(zip(tt_sorted,yy))
    values = [ D.get(t,0) for t in dTf]
    Hu = np.array(values).reshape(dT.shape)
    
    H = np.hstack([Hic, Hu])
    
    return H

# ===  Promact  =================================================


def infer_promact(curve_fluo, curve_volume, ttu, drna, kr, dR,
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
    
    ic
      Estimation of the initial concentration of rna/reporter/mature reporter

    alpha
      The value of alpha chosen by the algortihm.
    """

    if alphas is None:
      alphas = DEFAULT_ALPHAS

    def model_VR(Y, t, promact):
        """
        dRrna/dt = promact(t) - dRrna Rrna(t)
        dRi/dt = Rrna(t) - (dR+kR) Ri(t)
        dR/dt = kR Ri(t) - dR R(t)
        """
        kM = kru = 1.0
        Rrna, Ri, R = Y
        return [ kM*promact(t) - drna * Rrna,
                 kru*Rrna - (dR+kr) * Ri,
                 kr * Ri - dR * R ]
    
    


    tt_fluo = curve_fluo.x 
    ttu = tt_fluo[:-1]
    dttu = ttu[1]-ttu[0]
    H_F = make_H( model_VR, obs_R, ttu, tt_fluo )
    H_F[:,len(obs_R):] *= curve_volume(ttu+.5*dttu)

    synth_rate, fluo_smooth , ic, alpha, ascores = infer_control(
               H_F, curve_fluo.y, Nic=3, alphas=alphas, eps_L=1e-6)

    return (Curve(ttu, synth_rate), Curve(tt_fluo, fluo_smooth),
            ic, alpha, ascores)




def infer_prot_conc_multistep(curve_fluo, curve_volume, ttu, drna, kr, dR,
                    dP, alphas=None, eps_L=0.0001):
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
    
    if alphas is None:
        alphas = DEFAULT_ALPHAS

    tt_fluo = curve_fluo.x

    def model_VR(Y, t, promact):
        """
        dRrna/dt = promact(t) - dRrna Rrna(t)
        dRi/dt = Rrna(t) - (dR+kR) Ri(t)
        dR/dt = kR Ri(t) - dR R(t)
        """
        kM = kru = 1.0
        Rrna, Ri, R = Y
        return [ kM*promact(t) - drna * Rrna,
                 kru*Rrna - (dR+kr) * Ri,
                 kr * Ri - dR * R ]


    def model_VP(Y, t, promact):
        """ dPrna/dt = promact(t) - dPrna Prna(t)
            dP/dt = Prna - dP P(t) """
    
        Prna, P = Y
        return [ promact(t) - drna * Prna,
                 Prna - dP * P ]
    
    
    nvars_P, nvars_F = len(obs_P), len(obs_R)
    
    # ttu_promact represents the times of the promact control.
    # They must be in the same number as the ttu of the protein of
    # interest, but their maximum must be lower than the maximum of
    # the ttu of the protein of interest. Only then can the matrix
    # H_P be inverted.
    t_end = curve_volume.x.max()
    d = 2
    ttu_promact = np.linspace(0, t_end, len(ttu)-nvars_P+d)[:-d]

    H_F = make_H( model_VR, obs_R, ttu_promact, tt_fluo )
    H_P = make_H( model_VP, obs_P, ttu_promact, ttu )
    
    

    H_Pi = scipy.linalg.inv(H_P)
    
    H_Pi2V = ( H_Pi * curve_volume( ttu ) )[nvars_P:,:]
     
    Ny, Nx = H_Pi2V.shape
    Z_1 = 1.0* np.zeros((nvars_F, Nx))
    Z_2 = 1.0* np.zeros((Ny, nvars_F))
    I_F = 1.0* np.identity(nvars_F)

    H1 = np.vstack( [np.hstack( [ I_F , Z_1 ] ),
                     np.hstack( [ Z_2 , H_Pi2V ] )])
    
    HpF = H_F.dot( H1 )
    
    ######

    u, y , ic, alpha, scores = infer_control(HpF, curve_fluo.y, Nic=3,
                                             alphas=alphas, eps_L=1e-6)

    return Curve(ttu, u), Curve(tt_fluo, y), ic, alpha, scores