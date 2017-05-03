"""
This module contains a function load_WI_JSON_export to import WellInverter exported data into wellFARE.
load_WI_JSON_export returns a wells dictonary.

The individual well data is stored in wells[wellname][measurename]

The available measure names are stored in wells["measures"]
The group list is stored in wells["groups"].keys()

The group mean, std and sem are stored in wells["groups"][groupname]["SD" or "Mean" or "SEM"][measurename]
The list of wells in groups is stored in wells["groups"][groupname]["WellList"]
"""

from ..curves import Curve
import json


def load_WI_JSON(filename):
    with open(filename) as data_file:
        data = json.load(data_file)

    wells = {}
    wells["groups"] = {}
    wells["measures"] = []

    for measurename in data:
        for datatype in data[measurename]:
            dataname = measurename + '-' + datatype
            if dataname not in wells["measures"]:
                wells["measures"].append(dataname)
            for wellorgroupname in data[measurename][datatype]:
                if (data[measurename][datatype][wellorgroupname].get("WellList") == None):
                    # A unique well
                    if wells.get(wellorgroupname) == None:
                        wells[wellorgroupname] = {}
                    wells[wellorgroupname][dataname] = Curve(data[measurename][datatype][wellorgroupname].get(
                        'time'), data[measurename][datatype][wellorgroupname].get('signal'))
                else:
                    # A group
                    if wells["groups"].get(wellorgroupname) == None:
                        wells["groups"][wellorgroupname] = {}
                    for stat in data[measurename][datatype][wellorgroupname]:
                        if wells["groups"][wellorgroupname].get(stat) == None:
                            wells["groups"][wellorgroupname][stat] = {}
                        if stat == "WellList":
                            wells["groups"][wellorgroupname][stat] = data[measurename][datatype][wellorgroupname][stat]
                        else:
                            wells["groups"][wellorgroupname][stat][dataname] = Curve(data[measurename][datatype][wellorgroupname][stat].get(
                                'time'), data[measurename][datatype][wellorgroupname][stat].get('signal'))

    return wells
