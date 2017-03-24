"""
This module contains base functions to parse different layouts
of the excel files produced by the TECAN infinite pro.


"""


import re
import datetime
import xlrd
import numpy as np

from ..curves import Curve

def parse_tecan(filename, sheet_index=None, info=False):
    """ Parses a .xlsx file from a cinetic experiment

    File specifications:


    """
    sheets = workbook2numpy(filename, sheet_index=sheet_index)

    if info :
        info_dict = get_info(sheets)

    if isinstance(sheets, list):
        starts = map(find_start_in_sheet, sheets)
        t0 = min([ s for s in starts if (s is not None)])
        if info:
            return [[parse_sheet(sheet, t0=t0)[1] for sheet in sheets], info_dict]
        else:
            return [parse_sheet(sheet, t0=t0)[1] for sheet in sheets]
    else:
        if info:
            return [parse_sheet(sheets), info_dict]
        else:
            return parse_sheet(sheets)


def get_info(sheets):

    info_dict = {}

    if isinstance(sheets, list):
        i = 0
        while len(sheets[i][0]) == 0:
            i += 1
        sheet = sheets[i]

    i = 0
    print("SHEET", sheet.shape)
    modeindex = 0
    nameindex = 0
    while i < sheet.shape[0]:
        if sheet[i][0].startswith('List of actions'):
            print("ACTIONS",i)
            info_dict["actions"] = []
            i+=1
            while not ('Label' in sheet[i][0]):
                if len(sheet[i][0]) != 0 :
                    linelist = [var for var in sheet[i] if var]
                    info_dict["actions"].append([ linelist[0] , ' '.join(linelist[1:]) ])
                i += 1

        if sheet[i][0].startswith('Mode'):
            print("MODE", i)
            linelist = [var for var in sheet[i] if var]
            info_dict[modeindex] = [[ linelist[0] , ' '.join(linelist[1:]) ]]
            i += 1
            while not (sheet[i][0].startswith('Mode') or len(sheet[i][0]) == 0 or sheet[i][0].startswith('Start Time') ):
                linelist = [var for var in sheet[i] if var]
                info_dict[modeindex].append([ linelist[0] , ' '.join(linelist[1:]) ])
                i += 1

            if len(sheet[i][0]) != 0:
                i -= 1
            modeindex += 1

        if sheet[i][0].startswith('Start Time'):
            linelist = [var for var in sheet[i] if var]
            info_dict["Start Time"] = ' '.join(linelist[1:])

        if sheet[i][0].startswith('Cycle Nr'):
            info_dict[nameindex].append(['Name',sheet[i-1][0]])
            nameindex += 1

        i += 1

    return info_dict


def workbook2numpy(filename, sheet_index=None):
    """ loads the xlsx file as a (Numpy) array, or list of
        numpy arrays if there are several sheets.
        If `sheetindex` is None, """

    book = xlrd.open_workbook(filename)
    sheets = np.array(book.sheets())
    if sheet_index is None:
        sheet_index = range(len(sheets))
    if np.isscalar(sheet_index):
        return sheet2numpy(sheets[0])
    else:
        res = []
        for sh in sheets[sheet_index]:
            try:
                res.append(sheet2numpy(sh))
            except:
                pass

    return res[0] if len(res)==1 else res



def find_start_in_sheet(sheet):
    for line in sheet:
        if len(line) == 0:
            pass
        elif line[0] == "Start Time:":
            return date2seconds(line[1])
    return None



def sheet2numpy(sheet):
    """ Conversts a xlread Excel sheet to a numpy array """
    X,Y = sheet.ncols,  sheet.nrows
    arr = [[sheet.cell(y,x).value for x in range(X)]
              for y in range(Y)]
    return np.array(arr)



def parse_sheet(sheet, t0 = None):

    wells_dict = { "%s%d"%(c,i) : dict() for c in "ABCDEFGH"
                                         for i in range(1,13) }
    start_time = 0
    for i,line in enumerate(sheet):
        if len(line) == 0:
            pass
        elif line[0] == "Start Time:":
            start_time = date2seconds(line[1])
            if t0 is None:
                t0 = start_time
            start_time = start_time-t0
            parse_labels(sheet,i, wells_dict, start_time)

    return t0, wells_dict


def parse_labels(sheet,i, wells_dict, start_time):
    """
    Parses the different labels encountered (and fills the given
    plate until an "End Time:" cell is found.
    """
    j = i
    while sheet[j][0] != "End Time:":
        if sheet[j][0] == "Cycle Nr.":
            parse_label(sheet,j, wells_dict,  start_time)
        j +=1



def parse_label(sheet,i, wells_dict,  start_time=0,
                timePerWell=True, over_replace = -1,
                per_column=False):
    """
    Parses an array of measurements, supposing that
    line i of arr in the first line of an array of the form:
    Cycle Nr,  1,  2
    Time [s],  0,  34.5
    Temp. [C], 23.5, 23.5
    A1, 0.3174999952, 0.3181999922
    t   00            30 <- time per well activated
    A2, 0.3980999887, 0.4104000032
    t   02            32
    """
    label = sheet[i-1,0]


    if per_column:
        sheet = sheet[i:,:].T
        i=0


    try:
        xmax = list(sheet[i]).index(u'') - 1
    except:
        xmax = len(list(sheet[i]))

    if sheet[i+1][1] == '':
        return #return if the first data element is empty (meaning all data should be empty)

    if not timePerWell:
        # read the times once and for all
        tt = sheet[i+1, 1:xmax].astype(float)/60000 + start_time
        j = i+3
    else:
        j=i+2

    while (j<sheet.shape[0]) and (sheet[j,0] != u''):

        try:
            xmax = list(sheet[j]).index(u'') - 1
        except:
            xmax = len(list(sheet[j]))

        try:
            well = sheet[j,0]

            if timePerWell:
                tt = sheet[j+1, 1:xmax].astype(float)/60000 + start_time

            yy = sheet[j,1:xmax]
            yy[yy == 'OVER'] = over_replace


            curve = Curve(tt.astype(float), yy.astype(float))

            if not (label in wells_dict[well].keys()):
                wells_dict[well][label] = curve
            else:
                wells_dict[well][label] = wells[well][label].merge_with(curve)

            j += 2 if timePerWell else 1
        except:
            j += 2 if timePerWell else 1
            continue
            pass


def merge_wells_dicts(wells_dicts):
    """
    Merges the dictionnaries
    """

    result = { "%s%d"%(c,i) : dict() for c in "ABCDEFGH"
                                         for i in range(1,13) }
    for wells_dict in wells_dicts:

        for well, curves_dict in wells_dict.items():

            for label,curve in curves_dict.items():

                if not (label in result[well].keys()):
                    result[well][label] = curve
                else:
                    result[well][label] =\
                        result[well][label].merge_with(curve)
    return result


def to_coords(s):
    """ Converts "A5", "C11", ... into (0,5), (2,11) ... """
    return  ( "ABCDEFGH".index(s[0]), int(s[1:])-1 )



def date2seconds(timeString):
    """
    Converts a Tecan date string ("08/10/2013 17:44:24")
    into seconds since 1970/1/1
    """
    template = "(\d+)/(\d+)/(\d+) (\d+):(\d+):(\d+)"
    matching = re.match(template,timeString)
    day,mth,yr,hr,mn,sec = map(int,matching.groups())
    t1 = datetime.datetime(yr,mth,day,hr,mn,sec)
    t0 = datetime.datetime(1970,1,1)
    return (t1-t0).total_seconds()
