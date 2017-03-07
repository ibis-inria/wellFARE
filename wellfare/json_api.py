"""
This module contains helper functions to communicate with wellfare
using JSON dictionnaries. The main function is json_process, which
will call one of the subsequent functions depending on the TASK
to be performed.
"""

import numpy as np
from .curves import Curve

from .ILM import (infer_growth_rate,
                  infer_synthesis_rate_onestep,
                  infer_synthesis_rate_multistep,
                  infer_prot_conc_onestep,
                  infer_prot_conc_multistep)

from .preprocessing import filter_outliers, filter_outliersnew, calibration_curve

from .parsing import parse_tecan, merge_wells_dicts


DEFAULTS = {
    'n_control_points':100,
    'dRNA': 1.0,
    'eps_L' : 0.000001,
    'alphalow' : -10,
    'alphahigh' : 10,
    'nalphastep' : 1000,
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
            'outliersnew': wellfare_outliersnew,
            'synchronize': wellfare_synchronize,
            'subtract': wellfare_subtract,
            'calibrationcurve': wellfare_calibrationcurve,
            'parsetecan': wellfare_parsetecan

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
        'positive' : boolean
        'alphalow' : -10,
        'alphahigh' : 10,
        'nalphastep' : 1000,
        'eps_L' : 0.000001
       }

    Output :
      { 'times_growth_rate': [...],
        'values_growth_rate': [...]
       }
    """

    curve_v = Curve(data['times_volume'],
                    data['values_volume'])

    positive = False
    if 'positive' in data:
        positive = data['positive']

    check_noNaN(curve_v.y, "curve_v.y", "wellfare_growth")

    n_control_points = get_var_with_default(data, 'n_control_points')
    ttu = np.linspace(curve_v.x.min(), curve_v.x.max(), n_control_points+3)[:-3]


    eps_L = get_var_with_default(data, 'eps_L')
    alphalow = get_var_with_default(data, 'alphalow')
    alphahigh = get_var_with_default(data, 'alphahigh')
    nalphastep = get_var_with_default(data, 'nalphastep')
    alphas = np.logspace(alphalow,alphahigh,nalphastep)

    growth, volume, _, _, _ = infer_growth_rate(curve_v, ttu,
                                                alphas=alphas, eps_L=eps_L, positive=positive)

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
        'n_control_points':100, // 100 is the default
        'alphalow' : -10,
        'alphahigh' : 10,
        'nalphastep' : 1000,
        'eps_L' : 0.000001
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

    n_control_points = get_var_with_default(data, 'n_control_points')
    ttu = np.linspace(curve_v.x.min(), curve_v.x.max(), n_control_points+3)[:-3]

    eps_L = get_var_with_default(data, 'eps_L')
    alphalow = get_var_with_default(data, 'alphalow')
    alphahigh = get_var_with_default(data, 'alphahigh')
    nalphastep = get_var_with_default(data, 'nalphastep')
    alphas = np.logspace(alphalow,alphahigh,nalphastep)

    if 'kR' in data:
        # use a two-step model of reporter expression
        # if no dRNA provided it is supposed to be very short-lived so that
        # the transcription step won't impact the dynamics of gene expression
        dRNA = get_var_with_default(data, 'dRNA')
        kR = data['kR']
        synth_rate, _, _, _, _ = infer_synthesis_rate_multistep(
            curve_v=curve_v,
            curve_f=curve_f,
            ttu = ttu,
            drna = dRNA,
            kr = kR,
            dR=dR,
            alphas=alphas,
            eps_L=eps_L)

    else:
        # use a one-step model of reporter expression
        synth_rate, _, _, _, _ = infer_synthesis_rate_onestep(
            curve_v=curve_v,
            curve_f=curve_f,
            ttu = ttu,
            degr=dR,
            alphas=alphas,
            eps_L=eps_L)

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
        'n_control_points': 100, // optional, 100 is default
        'alphalow' : -10,
        'alphahigh' : 10,
        'nalphastep' : 1000,
        'eps_L' : 0.000001
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

    eps_L = get_var_with_default(data, 'eps_L')
    alphalow = get_var_with_default(data, 'alphalow')
    alphahigh = get_var_with_default(data, 'alphahigh')
    nalphastep = get_var_with_default(data, 'nalphastep')
    alphas = np.logspace(alphalow,alphahigh,nalphastep)

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
            dP=dP,
            alphas=alphas,
            eps_L=eps_L)
    else:
        concentration, _, _, _, _ = infer_prot_conc_onestep(
            curve_v=curve_v,
            curve_f=curve_f,
            ttu = ttu,
            dR=dR, dP = dP,
            alphas=alphas,
            eps_L=eps_L)

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

def wellfare_outliersnew(data):
    """ Removes outliers from a curve using smooting spline.

    Command: 'outliersnew'

    Input:
      {
        'smoothing_win': int,
        'nstd': int/float,
        'iterations'= int,
      }

    Output:
      { 'times_cleaned_curve': [...],
        'values_cleaned_curve': [...]}

    SUMMARY OF THE PARAMETERS
    --------------------------
    smoothing_win -> window size (the more, the smoother)
    nstd -> Keep points that are less than "nstd*std" from the smooth.
    iterations --> number of time to do the smoothing and cut off using nstd
    """

    curve = Curve(data.pop('times_curve'),
                  data.pop('values_curve'))

    cleaned_curve = filter_outliersnew(curve, **data)

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


def wellfare_calibrationcurve(data):
    """ Returns the calibration curve (i.e. Fluo = f(Abs)) using polynomial fit.

    The returned curve gives to the autofluorescence of the well as a function of its absorbance.

    Command: calibrationcurve

    Input:
     { 'abs_time': [...],
       'abs_value': [...],
       'fluo_time': [...],
       'fluo_value': [...] }

    Output:
     { 'calibrationcurve_abs'  : [...],
       'calibrationcurve_fluo' : [...] }
    """

    abscurve = Curve(data['abs_time'],
                    data['abs_value'])
    fluocurve = Curve(data['fluo_time'],
                    data['fluo_value'])

    calibrationcurve, calibrationcurve_smoothextrapolation = calibration_curve(abscurve,fluocurve,data['smoothing'],data['extrapoldistance'],data['validinterval'])

    return {'calcurve_time': list(calibrationcurve.x.astype(float)),
            'calcurve_value': list(calibrationcurve.y.astype(float)),
            'calcurvesmoothed_time': list(calibrationcurve_smoothextrapolation.x.astype(float)),
            'calcurvesmoothed_value': list(calibrationcurve_smoothextrapolation.y.astype(float))}



def wellfare_parsetecan(data):
    """ Returns a dict containing parsed data


    Command: parsetecan

    Input:
     { 'inputfilename': str }

    Output:
     { 'parsedfile'  : dict }
    """

    filename = data['inputfilename']

    parsed_sheets, infodict = parse_tecan(filename, info=True)
    wells = merge_wells_dicts(parsed_sheets)

    jsonfile = {}
    if 'Start Time' in infodict:
        jsonfile['initialTime'] = infodict['Start Time']

    measureindex = 0
    jsonfile['measureTypes'] = {}
    jsonfile['programs'] = []
    while measureindex in infodict:
        programinfo = {}
        for info in infodict[measureindex]:
            if info[0] == 'Mode':
                mode = info[1]
            elif info[0] == 'Name':
                name = info[1]
            programinfo[info[0]] = info[1]

        measuretype = -1
        if 'Abs' in mode:
            measuretype = 0
        elif 'Fluo' in mode:
            measuretype = 1

        jsonfile['measureTypes'][name] = measuretype
        jsonfile['programs'].append(programinfo)

        measureindex += 1

    jsonfile['actions'] = infodict["actions"]

    jsonfile['wells'] = {}
    for wellname in wells:
        well = {}
        well['measures'] = {}
        for measurename in jsonfile['measureTypes']:
            if measurename in wells[wellname]:
                measure = {}
                measure["time"] = list(wells[wellname][measurename].x.astype(float))
                measure["originalSignal"] = list(wells[wellname][measurename].y.astype(float))
                well['measures'][measurename] = measure

        jsonfile['wells'][wellname] = well


    return  { 'parsedfile'  : jsonfile }
