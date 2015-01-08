"""
This module contains base functions to parse different layouts
of the excel files produced by the Fusion (add ref).


"""


import re
import datetime
import xlrd
import numpy as np
import csv

from ..curves import Curve

### ====   TIME  TOOLS =======

def to_minutes(seconds):
    return seconds * 1.0 / 60

def to_seconds(hours,minutes,secs):
    return 3600 * hours + 60 * minutes + secs

def str_to_seconds(s):
    """ s is a string like '4:5:2' or '04:05:02', etc... """
    return to_seconds(*map(int,s.rsplit(':')))

def seconds_delta(secs1,secs2):
    """ delta between to times in seconds.  The modulo ensures
        positive deltas in case of over-midnight experiments """
    SECONDS_PER_DAY = 86400
    return (secs1 - secs2) % SECONDS_PER_DAY

### ===================================================


def parse_fusion(dataFile, descriptionFile = None, verbose = False, **kwargs):
    """ Creates a new plate dataset from the raw data of the
        fusion plate reader.

        You can provide the vocable as keyword-argument:
        vocable = { 'lum':'Lum', 'fluo':'Fluo', 'Abs':'Abs' }
    """
    
    CSV_VOCABLE = kwargs.pop('vocable',
                  dict(zip(['fluo','lum','Abs'],['Fluo','Lum','Abs'])))
                  
    res = csv.reader(open(dataFile, "rb"))
    

    wells = wells_dict = { "%s%d"%(c,i) : dict() for c in "ABCDEFGH"
                                         for i in range(1,13) }
    
    # Transform all npArrays into lists to speed up appending
    
    
    # Read the file
    
    reference_time = None
    res = csv.reader(open(dataFile, "rb"))
    row = res.next()
    is_24h_passed = False # trick to handle 24h < t <48h experiments 
    
    try:
        cpt=0
        while(1):
            
            sRow = map(lambda s : s.strip(' '),row)
            
            if verbose :
                print sRow
                
            if len(sRow) and ( sRow[0] == 'Assay:'):
                
                parsed_variable = CSV_VOCABLE[sRow[1]]
                for well, dico in wells.items():
                    if parsed_variable not in dico.keys():
                        dico[parsed_variable] = {'times':[], 'values':[]}  
                
            elif len(sRow) and (sRow[0] == "Well Number"):
                
                r = res.next()
                
                while len(r):
                    
                    val_list = map(lambda s : s.strip(' '),r)
                    well_num = int(val_list[0])
                    #print well_num
                    well_name = "ABCDEFGH"[(well_num-1)/12] + str( ((well_num-1)%12)+1)
                    flag = val_list[-1]
                    
                    # It is possible that one line of the csv
                    # contains more than one measurement, at which
                    # case it is of the form
                    # num, t1,t2,t3, val1,val2,val3, descr, flag
                    
                    n_measures = ( len(val_list) - 3 ) / 2
                    
                    for i in range(n_measures) :
                        
                        hour = val_list[1+i]
                        value = val_list[1+i+n_measures]
                        value = float(value)        
                    
                        current_time = str_to_seconds(hour)
                        
                        if (reference_time == None ):
                            reference_time = current_time
                           
                        t = to_minutes(seconds_delta(current_time,reference_time))
                        
                        if (t > 60*24):
                            is_24h_passed = True    
                        elif is_24h_passed:
                            
                            t = t + 60*24
                            
                        wells[well_name][parsed_variable]['times'].append(t)
                        wells[well_name][parsed_variable]['values'].append(float(value))
                    
                    r = res.next()  
                    if verbose:
                        print(r)
                    
            row = res.next()
            
    except StopIteration:
        pass
    
    # Transform back all the data into proper Curves
    for well in wells.values():
        for var, dico in well.items():
            well[var] = Curve( dico['times'], dico['values'])

    return wells
