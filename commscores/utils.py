from modelseedpy.community.commhelper import build_from_species_models
from modelseedpy.core.msminimalmedia import MSMinimalMedia
from modelseedpy.core.fbahelper import FBAHelper
from numpy import load, nan, ndarray
from typing import Iterable
from math import isclose
import sigfig
import os
import re


class BadModels(Exception):
    def __init__(self, message):
        print(message)


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
sugars_dic, aminoacids_dic = (
    dict(zip(sugars[:, 0], sugars[:, 1])),
    dict(zip(aminoacids[:, 0], aminoacids[:, 1])),
)
vitamins_dic, minerals_dic = (
    dict(zip(vitamins[:, 0], vitamins[:, 1])),
    dict(zip(minerals[:, 0], minerals[:, 1])),
)
energy_compounds_dic = dict(zip(energy_compounds[:, 0], energy_compounds[:, 1]))


def _categorize_mets(metIDs):
    met_sugars, met_aminoAcids, met_vitamins, met_minerals, met_energy, met_other = (
        [],
        [],
        [],
        [],
        [],
        [],
    )
    for metID in metIDs:
        if metID in sugars[:, 0]:
            met_sugars.append(f"{sugars_dic[metID]} ({metID})")
        elif metID in aminoacids[:, 0]:
            met_aminoAcids.append(f"{aminoacids_dic[metID]} ({metID})")
        elif metID in vitamins[:, 0]:
            met_vitamins.append(f"{vitamins_dic[metID]} ({metID})")
        elif metID in minerals[:, 0]:
            met_minerals.append(f"{minerals_dic[metID]} ({metID})")
        elif metID in energy_compounds[:, 0]:
            met_energy.append(f"{energy_compounds_dic[metID]} ({metID})")
        else:
            met_other.append(metID)
    return met_sugars, met_aminoAcids, met_vitamins, met_minerals, met_energy, met_other


rm_costless = re.compile("(\s\(.+\))")


def remove_metadata(element):
    try:
        element = float(rm_costless.sub("", str(element)).replace("%", ""))
    except:
        pass
    return element


def convert_to_int(element):
    try:
        element = int(element)
    except:
        pass
    return element


def _process_mets(metIDs):
    return [", ".join(lst) for lst in _categorize_mets(metIDs)]


def _compatibilize(member_models: Iterable, printing=False):
    return member_models


def _load_models(
    member_models: Iterable, com_model=None, compatibilize=True, commID=None, printing=False
):
    # ic(member_models, com_model, compatibilize)
    if not com_model and member_models:
        return member_models, build_from_species_models(member_models, commID, "CommScores_community")  # (model, names=names, abundances=abundances)
    # elif com_model and not member_models:
    #     return com_model.members, com_model  # TODO the individual models of a community model can be parsed
    if compatibilize:
        return (
            _compatibilize(member_models, printing),
            _compatibilize([com_model], printing)[0],
        )
    return member_models, com_model


def _get_media(
    media=None,
    com_model=None,
    model_s_=None,
    min_growth=None,
    environment=None,
    interacting=True,
    printing=False,
    minimization_method="minFlux",
    skip_bad_media=False,
):
    if com_model is None and model_s_ is None:
        raise TypeError("< com_model > or < model_s_ > must be parameterized.")
    if media is not None:
        if model_s_ is not None and not isinstance(model_s_, (list, set, tuple)):
            return media["members"][model_s_.id]["media"]
        elif com_model is not None:
            return media["community_media"]
        return media
    com_media = None
    if com_model is not None:
        while com_media is None:
            com_media, media_sol = MSMinimalMedia.determine_min_media(
                com_model,
                minimization_method,
                min_growth,
                None,
                interacting,
                5,
                printing,
            )
            min_growth *= 1.1
        if model_s_ is None:
            return com_media, media_sol
    if model_s_ is not None:
        if not isinstance(model_s_, (list, set, tuple, ndarray)):
            while com_media is None:
                com_media, media_sol = MSMinimalMedia.determine_min_media(
                    com_model,
                    minimization_method,
                    min_growth,
                    None,
                    interacting,
                    5,
                    printing,
                )
                min_growth *= 1.1
            return com_media, media_sol
        members_media = {}
        for model in model_s_:
            # print(model.id)
            while com_media is None:
                com_media, media_sol = MSMinimalMedia.determine_min_media(
                        model,
                        minimization_method,
                        min_growth,
                        environment,
                        interacting,
                        printing,
                    )[0]
                minimal_growth *= 1.1
            members_media[model.id] = {"media": (com_media, media_sol)}
        # print(members_media)
        if com_model is None:
            return members_media
        return {"community_media": com_media, "members": members_media}
    raise BadModels(f"The parameterized community model of type {type(com_model)} and member models {model_s_} are not properly captured.")


def _sigfig_check(value, sigfigs, default):
    if str(value) in ["inf", "nan"]:
        value = ""
    if FBAHelper.isnumber(value):
        return sigfig.round(value, sigfigs)
    else:
        return default


def _calculate_jaccard_score(set1, set2):
    if set1 == set2:
        print(f"The sets are identical, with a length of {len(set1)}.")
    if len(set1.union(set2)) == 0:
        return (None, None)
    return (
        set1.intersection(set2),
        len(set1.intersection(set2)) / len(set1.union(set2)),
    )



def _check_model(model_util, media, model_str=None, skip_bad_media=True):
    default_media = model_util.model.medium
    # print("test")
    if media is not None:
        model_util.add_medium(media)
    obj_val = model_util.model.slim_optimize()
    model_str = model_str or model_util.model.id
    # print(model_str, obj_val)
    if isclose(obj_val, 0, abs_tol=1e-6) or not FBAHelper.isnumber(obj_val):
        print(f"The {model_str} model is not operational on the provided media")
        if not skip_bad_media:
            print(" and will be gapfilled.")
            return MSGapfill.gapfill(model_util.model, media)
        model_util.add_medium(default_media)
    return model_util.model


def _determine_growths(modelUtils, environ):
    obj_vals = []
    for util in modelUtils:
        util.add_medium(environ)
        obj_vals.append(util.model.slim_optimize())
    return obj_vals


def nanFilter(value, string=True):
    if isinstance(value, str) or value is None:
        if string:
            return value
        else:
            return nan
    if any([value < 0, value > 1e5]):
        return "" if string else nan
    return value
