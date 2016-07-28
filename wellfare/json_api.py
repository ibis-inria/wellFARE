"""
This module contains helper functions to communicate with wellfare
using JSON dictionnaries. The main function is json_process, which
will call one of the subsequent functions depending on the TASK
to be performed.
"""

import numpy as np
from .curves import Curve

from .ILM import (infer_growth_rate,
                  infer_synthesis_rate,
                  infer_prot_conc_onestep,
                  infer_prot_conc_multistep)
                  
from .preprocessing import filter_outliers 



DEFAULTS = {
    'n_control_points':100,
    'dRNA': 1.0,
}
def get_var_with_default(data, var):
    if var in data:
        return data.get(var)
    elif var in DEFAULTS:
        return DEFAULTS[var]
    else:
        raise ValueError("Variable %s was not provided and no default value"%var
                        +" is known for this variable (check spelling ?)")

def check_noNaN(array, name, fun, additional_message=''):
    if np.isnan(np.sum(array)):
        raise AssertionError("Error: Array '%s' in function %s has NaNs ! %s"%(
                              name, fun, additional_message))

# THE MAIN FUNCTION, CALLED BY THE PYTHON/JS PROCESS:

def json_process(command, input_data):
    """ Calls the right function depending on the ``command``.

    This function is a 'hub': it will decide which function to
    apply to the data, depending on the command.
    For inputs and ouputs, see the doc of the different functions
    below.
    """

    return {'growth': wellfare_growth,
            'activity': wellfare_activity,
            'concentration': wellfare_concentration,
            'outliers': wellfare_outliers,
            'synchronize': wellfare_synchronize,
            'subtract': wellfare_subtract

           }[command](input_data)


# THE SPECIFIC FUNCTIONS, ONE FOR EACH TASK:


# === INFERENCE ============================================

def wellfare_growth(data):
    """ Computes the growth rate from volume data.

    Command : 'growth'

    Input :
      { 'times_volume': [...] ,
        'values_volume': [...],
        'n_control_points': 100 // optional, 100 is default
       }

    Output :
      { 'times_growth_rate': [...],
        'values_growth_rate': [...]
       }
    """
    print("Starting ---")
    curve_v = Curve(data['times_volume'],
                    data['values_volume'])

    check_noNaN(curve_v.y, "curve_v.y", "wellfare_growth")
    
    n_control_points = get_var_with_default(data, 'n_control_points')
    ttu = np.linspace(curve_v.x.min(), curve_v.x.max(), n_control_points+3)[:-3]
    
    print("Starting computations")
    alphas = 10.0**np.linspace(-5,8,1000)
    growth, volume, _, _, _ = infer_growth_rate(curve_v, ttu,
                                                alphas=alphas, eps_L=1e-6)
    print("finished computations")
    check_noNaN(growth.y, "growth.y", "wellfare_growth")
    
    return {'times_growth_rate': list(growth.x.astype(float)),
            'values_growth_rate': list(growth.y.astype(float))}



def wellfare_activity(data):
    """ Computes protein synthesis rate, or promoter activity,
    using a simple one-step model for the GFP synthesis. 

    Command : 'activity'
    
    Input:
      { 'times_volume': [...] ,
        'values_volume': [...],
        'times_fluo': [...],
        'values_fluo': [...],
        'dR': float, // degradation constant of the reporter
        'kR': float // (optional) folding constant of the reporter.
        'dRNA': float // (optional) degradation constant of the RNA.
        'n_control_points':100 // 100 is the default
       }

    Output:
      { 'times_activity': [...],
        'values_activity': [...]
       }
    """

    curve_v = Curve(data['times_volume'],
                    data['values_volume'])
    curve_f = Curve(data['times_fluo'],
                    data['values_fluo'])

    dR = data['dR']
    send_state = data['send_state']

    n_control_points = get_var_with_default(data, 'n_control_points')
    ttu = np.linspace(curve_v.x.min(), curve_v.x.max(), n_control_points+3)[:-3]

    if 'kR' in data:
        # use a two-step model of reporter expression
        # if no dRNA provided it is supposed to be very short-lived so that
        # the transcription step won't impact the dynamics of gene expression
        dRNA = data.get('dRNA', 1.0)
        synth_rate, _, _, _, _ = infer_promact(
            curve_v=curve_v,
            curve_f=curve_f,
            ttu = ttu,
            drna = dRNA,
            kr = kR,
            dR=dR)

    else:
        # use a one-step model of reporter expression
        synth_rate, _, _, _, _ = infer_synthesis_rate(
            curve_v=curve_v,
            curve_f=curve_f,
            ttu = ttu,
            degr=dR,
            send_state = send_state)

    return {'times_activity': list(synth_rate.x.astype(float)),
            'values_activity': list(synth_rate.y.astype(float))}



def wellfare_concentration(data):
    """ Computes the concentration of a protein from
    a fluorescence curve and an absorbance curve.

    Command: 'concentration'

    Input:
      { 'times_volume': [...] ,
        'values_volume': [...],
        'times_fluo: [...],
        'values_fluo: [...],
        'dR': float,
        'dP': float,
        'n_control_points': 100 // optional, 100 is default
       }

    Output:
      { 'times_concentration': [...],
        'values_concentration': [...]
       }
    """


    curve_v = Curve( data['times_volume'], 
                     data['values_volume'])
    curve_f = Curve( data['times_fluo'],
                     data['values_fluo'])
    dR = data['dR']
    dP = data['dP']
    
    n_control_points = get_var_with_default(data, 'n_control_points')
    ttu = np.linspace(curve_v.x.min(), curve_v.x.max(), n_control_points+3)[:-3]
    
    if 'kR' in data:
        # use a two-step model of reporter expression
        # if no dRNA provided it is supposed to be very short-lived so that
        # the transcription step won't impact the dynamics of gene expression
        dRNA = get_var_with_default(data, 'dRNA')
        kR = data['kR']
        concentration, _, _, _, _ = infer_prot_conc_multistep(
            curve_v=curve_v,
            curve_f=curve_f,
            ttu = ttu,
            drna= dRNA,
            kr=kR,
            dR=dR,
            dP=dP)
    else:
        concentration, _, _, _, _ = infer_prot_conc_onestep(
            curve_v=curve_v,
            curve_f=curve_f,
            ttu = ttu,
            dR=dR, dP = dP)

    return {'times_concentration': list(concentration.x.astype(float)),
            'values_concentration': list(concentration.y.astype(float))}



# ===  PREPROCESSING ==============================================

def wellfare_outliers(data):
    """ Removes outliers from a curve.
    

    Command: 'outliers'

    Input:
      {
        'times_curve': [...],
        'values_curve': [...],
        'percentile_above'= int/float,
        'percentile_below'=int/float,
        'niter_above'= int,
        'niter_below'= int,
        'goal_above'= float,
        'goal_below'= float,
        'smoothing_win'= int,
        'nstd'= int,
        'above_first': True or False (abso or fluo)
      }

    Output:
      { 'times_cleaned_curve': [...],
        'values_cleaned_curve': [...]}

    SUMMARY OF THE PARAMETERS
    --------------------------

    percentile_above -> 1-100, proportion of up-lying points to keep 
    percentile_below -> 1-100, proportion of up-lying points to keep 
    niter_above -> Number of times the up-filter is repeated
    niter_below -> Number of times the down-filter is repeated
    goal_above -> the algo will stop when second derivatives are below
    goal_below -> the algo will stop when -(ddt) are below
    smoothing_win -> window size (the more, the smoother)
    nstd=1 -> Keep points that are less than "nstd*std" from the smooth.
    above_first -> filter upliers first ? (yes for abso, no for fluo)

    INDICATIVE VALUES OF THE PARAMETERS (may vary)
    -----------------------------------------------

    OD:
    'percentile_above':50,
    'percentile_below':50,
    'niter_above':4,
    'niter_below':3,
    'goal_above':0.001,
    'goal_below':0.001,
    'smoothing_win':4,
    'nstd':1,
    'above_first':True

    Fluo:
    'percentile_above':90,
    'percentile_below':90,
    'niter_above':2,
    'niter_below':2,
    'goal_above':1,
    'goal_below':5,
    'smoothing_win':4,
    'nstd':1.5,
    'above_first':False

    """
    
    curve = Curve(data.pop('times_curve'),
                  data.pop('values_curve'))

    cleaned_curve = filter_outliers(curve, **data)


    return { 'times_cleaned_curve': list(cleaned_curve.x.astype(float)),
             'values_cleaned_curve': list(cleaned_curve.y.astype(float)) }



def wellfare_synchronize(data):
    """ Returns the lag between two curves.

    Command: 'synchronize'

    If [x1],[y1] and [x2],[y2] are two curves,

    returns the time shift such that

    [x1],[y1] ~~ [x2 + d],[y2]

    Will only find shifts smaller than the provided 'shift_max'

    Input:
     { 'times_curve_1': [...],
       'values_curve_1': [...],
       'times_curve_2': [...],
       'values_curve_2': [...],
       'max_shift': float}

     Output:
       { 'time_shift': float }

    """

    curve_1 = Curve(data['times_curve_1'],
                    data['values_curve_1'])
    curve_2 = Curve(data['times_curve_2'],
                    data['values_curve_2'])
    max_shift = data['max_shift']

    shifts0 = np.arange(-max_shift, max_shift, 10)
    t_min = max(curve_1.x[0],curve_2.x[0]) + max_shift + 1
    t_max = min(curve_1.x[-1],curve_2.x[-1]) - max_shift - 1
    tt = np.linspace(t_min, t_max, 50)

    time_shift = curve_1.find_shift_gradient([curve_2], tt,
                                           shifts0 = shifts0)[0]
    
    return {'time_shift': time_shift}



def wellfare_subtract(data):
    """ Returns the difference between two curves.

    Will return the curve corresponding to the difference
    ``(curve1 - curve2)``. The times of the returned curve are the
    the times of the curve ``curve1``.
    The values of the returned curve are computed using linear
    interpolation when necessary.
    
    Command: subtract

    Input:
     { 'times_curve_1': [...],
       'values_curve_1': [...],
       'times_curve_2': [...],
       'values_curve_2': [...] }

    Output:
     { 'times_subtraction'  : [...],
       'values_subtraction' : [...] }
    """

    curve_1 = Curve(data['times_curve_1'],
                    data['values_curve_1'])
    curve_2 = Curve(data['times_curve_2'],
                    data['values_curve_2'])
    
    subtraction = (curve_1 - curve_2)
    new_x = np.array([x for x in curve_1.x if x in subtraction.x])
    
    return { 'times_subtraction': list(new_x),
             'values_subtraction': list(subtraction(new_x)) }