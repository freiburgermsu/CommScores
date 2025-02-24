from __future__ import absolute_import

# import non-local packages
from pandas import DataFrame, read_csv, concat, set_option
from cobra.io import write_sbml_model, read_sbml_model
from configparser import ConfigParser
from numpy import load, nan, ndarray
from zipfile import ZipFile
from typing import Iterable
from math import isclose
from pathlib import Path
import json, time, copy
import sys, os, re
import hashlib
import sigfig
import logging
import uuid

# define a local pointer to import package versions
path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
path2 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site_packages = "/home/afreiburger/.local/lib/python3.11/site-packages"
kbmodules_path = "/scratch/shared/code/KBBaseModules"
modelseed_path = "/scratch/shared/ModelSEEDpy_APF"
mscommunity_path = "/scratch/shared/code/MSCommunity"
utilsModule_path = "/scratch/shared/code/chenry_utility_module/lib"
msrecon_path = "/scratch/shared/code/KB-ModelSEEDReconstruction/lib"
commscores_path = "/scratch/shared/code/CommScores"
for p in [path, path2, site_packages, kbmodules_path, modelseed_path, mscommunity_path, utilsModule_path, msrecon_path, commscores_path]:
    sys.path.insert(0, p)
# print(sys.path) 

## cobrakbase imports
from cobrakbase.core.kbasefba import FBAModel

## modelseedpy imports
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core import MSMinimalMedia
from modelseedpy.helpers import get_template

## mscommunity imports
from mscommunity.commhelper import build_from_species_models

class BadModels(Exception):
    def __init__(self, message):
        print(message)

        
class CommScoresUtil():
    @staticmethod
    def _categorize_mets(metIDs):
        # load compound categories
        package_dir = os.path.abspath(os.path.dirname(__file__))
        categories_dir = os.path.join(package_dir, "data", "categories")
        sugars, aminoacids = (
            load(os.path.join(categories_dir, "sugars.npy")),
            load(os.path.join(categories_dir, "aminoAcids.npy")),
        )
        vitamins, minerals = (
            load(os.path.join(categories_dir, "vitamins.npy")),
            load(os.path.join(categories_dir, "minerals.npy")),
        )
        energy_compounds = load(os.path.join(categories_dir, "energy_compounds.npy"))
        energy_compounds_dic = dict(zip(energy_compounds[:, 0], energy_compounds[:, 1]))
        sugars_dic, aminoacids_dic = (
            dict(zip(sugars[:, 0], sugars[:, 1])),
            dict(zip(aminoacids[:, 0], aminoacids[:, 1])),
        )
        vitamins_dic, minerals_dic = (
            dict(zip(vitamins[:, 0], vitamins[:, 1])),
            dict(zip(minerals[:, 0], minerals[:, 1])),
        )
        
        # partition the compounds into their respective directories
        met_sugars, met_aminoAcids, met_vitamins, met_minerals, met_energy, met_other = ([], [], [], [], [], [])
        for metID in metIDs:
            if metID in sugars[:, 0]:    met_sugars.append(f"{sugars_dic[metID]} ({metID})")
            elif metID in aminoacids[:, 0]:  met_aminoAcids.append(f"{aminoacids_dic[metID]} ({metID})")
            elif metID in vitamins[:, 0]:    met_vitamins.append(f"{vitamins_dic[metID]} ({metID})")
            elif metID in minerals[:, 0]:    met_minerals.append(f"{minerals_dic[metID]} ({metID})")
            elif metID in energy_compounds[:, 0]:   met_energy.append(f"{energy_compounds_dic[metID]} ({metID})")
            else:    met_other.append(metID)
        return met_sugars, met_aminoAcids, met_vitamins, met_minerals, met_energy, met_other


    @staticmethod
    def remove_metadata(element):
        rm_costless = re.compile("(\s\(.+\))")
        try:
            element = float(rm_costless.sub("", str(element)).replace("%", ""))
        except:
            pass
        return element

    @staticmethod
    def convert_to_int(element):
        try:
            element = int(element)
        except:
            pass
        return element

    @staticmethod
    def _process_mets(metIDs):
        return [", ".join(lst) for lst in CommScoresUtil._categorize_mets(metIDs)]

    @staticmethod
    def _compatibilize(member_models: Iterable, printing=False):
        return member_models

    @staticmethod
    def _load_models(member_models: Iterable, com_model=None, compatibilize=True, commID=None, printing=False):
        # ic(member_models, com_model, compatibilize)
        if not com_model and member_models:
            return member_models, build_from_species_models(member_models, commID, "CommScores_community")  # (model, names=names, abundances=abundances)
        # elif com_model and not member_models:
        #     return com_model.members, com_model  # TODO the individual models of a community model can be parsed
        if compatibilize:
            return (
                CommScoresUtil._compatibilize(member_models, printing),
                CommScoresUtil._compatibilize([com_model], printing)[0],
            )
        return member_models, com_model


    @staticmethod
    def _get_media(
        media=None,
        com_model=None,
        model_s_=None,
        min_growth=0.1,
        environment=None,
        interacting=True,
        printing=True,
        minimization_method="minFlux",
        climit=None,
        o2limit=None
    ):
        assert com_model is not None or model_s_ is not None, "com_model or model_s_ must be parameterized."
        if com_model is True and isinstance(model_s_, (set, list, tuple)):
            com_model = build_from_species_models(model_s_)
        if media is not None:
            if model_s_ is not None and not isinstance(model_s_, (list, set, tuple)):
                return media["members"][model_s_.id]["media"]
            elif com_model is not None:
                return media["community_media"]
            return media
        com_media = None
        if com_model is not None:
            minGrowth = min_growth
            count = 0
            while com_media is None or com_media == {}:
                if count > 0:  minGrowth *= 1.1
                com_media, media_sol = MSMinimalMedia.determine_min_media(
                    com_model, minimization_method, minGrowth, environment, interacting, climit, o2limit, printing)
                count += 1
            if not isclose(minGrowth, min_growth, abs_tol=0.01):   print(f"{com_model.id} needed {minGrowth} v {min_growth}\n")
            if model_s_ is None:
                return com_media, media_sol
        if model_s_ is not None:
            min_media = None
            if not isinstance(model_s_, (list, set, tuple, ndarray)):
                minGrowth = min_growth
                count = 0
                while min_media is None or min_media == {}:
                    if count > 0:  minGrowth *= 1.1
                    min_media, media_sol = MSMinimalMedia.determine_min_media(
                        model_s_, minimization_method, minGrowth, environment, interacting, climit, o2limit, printing)
                    count += 1
                if not isclose(minGrowth, min_growth, abs_tol=0.01):   print(f"{model_s_.id} needed {minGrowth} v {min_growth}\n")
                return min_media, media_sol
            members_media = {}
            for model in model_s_:
                minGrowth = min_growth
                count = 0
                while min_media is None or min_media == {}:
                    if count > 0:  minGrowth *= 1.1
                    min_media, media_sol = MSMinimalMedia.determine_min_media(
                        model, minimization_method, minGrowth, environment, interacting, climit, o2limit, printing)
                    count += 1
                members_media[model.id] = {"media": min_media, "solution": media_sol}
                min_media = None
                if not isclose(minGrowth, min_growth, abs_tol=0.01):  print(f"{model.id} needed {minGrowth} v {min_growth}")
            if com_model is None:
                return members_media
            return {"community_media": com_media, "members": members_media}
        raise BadModels(f"The parameterized community model of type {type(com_model)} and member models {model_s_} are not properly captured.")

    @staticmethod
    def _sigfig_check(value, sigfigs, default):
        if str(value) in ["inf", "nan"]:
            value = ""
        if FBAHelper.isnumber(value):
            return sigfig.round(value, sigfigs)
        else:
            return default

    @staticmethod
    def _calculate_jaccard_score(set1, set2):
        if set1 == set2:
            print(f"The sets are identical, with a length of {len(set1)}.")
        if len(set1.union(set2)) == 0:
            return (None, None)
        return (
            set1.intersection(set2),
            len(set1.intersection(set2)) / len(set1.union(set2)),
        )


    @staticmethod
    def _check_model(model_util, media, model_str=None):
        default_media = model_util.model.medium
        # print("test")
        if media is not None:
            model_util.add_medium(media)
        obj_val = model_util.model.slim_optimize()
        model_str = model_str or model_util.model.id
        # print(model_str, obj_val)
        if isclose(obj_val, 0, abs_tol=1e-6) or not FBAHelper.isnumber(obj_val):
            print(f"The {model_str} model is not operational on the provided media")
            model_util.add_medium(default_media)
        return model_util.model


    @staticmethod
    def _determine_growths(modelUtils, environ, sigfig=5):
        obj_vals = []
        for util in modelUtils:
            util.add_medium(environ)
            val = util.model.slim_optimize()
            sigfiggedVal = CommScoresUtil._sigfig_check(val, sigfig, "")
            obj_vals.append(sigfiggedVal)
        return obj_vals

    @staticmethod
    def nanFilter(value, string=True):
        if isinstance(value, str) or value is None:
            if string:
                return value
            else:
                return nan
        if any([value < 0, value > 1e5]):
            return "" if string else nan
        return value

    
util = CommScoresUtil()