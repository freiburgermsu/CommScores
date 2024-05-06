from modelseedpy.community.commhelper import build_from_species_models
from modelseedpy.core.msminimalmedia import MSMinimalMedia
from modelseedpy.core.fbahelper import FBAHelper
from numpy import load, nan, ndarray
from typing import Iterable
import sigfig
import os
import re


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
    member_models: Iterable, com_model=None, compatibilize=True, printing=False
):
    # ic(member_models, com_model, compatibilize)
    if not com_model and member_models:
        model = build_from_species_models(member_models, name="SMETANA_pair")
        return member_models, model  # (model, names=names, abundances=abundances)
    # models = PARSING_FUNCTION(community_model)  # TODO the individual models of a community model can be parsed
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
    # print(media, com_model, model_s_)
    if com_model is None and model_s_ is None:
        raise TypeError("< com_model > or < model_s_ > must be parameterized.")
    if media is not None:
        if model_s_ is not None and not isinstance(model_s_, (list, set, tuple)):
            return media["members"][model_s_.id]["media"]
        elif com_model is not None:
            return media["community_media"]
        return media
    # model_s_ is either a singular model or a list of models
    if com_model is not None:
        try:
            com_media, media_sol = MSMinimalMedia.determine_min_media(
                com_model,
                minimization_method,
                min_growth,
                None,
                interacting,
                5,
                printing,
            )
        except Exception as e:
            if skip_bad_media:
                com_media, media_sol = None, None
            else:
                print(e)
    if model_s_ is not None:
        if not isinstance(model_s_, (list, set, tuple, ndarray)):
            try:
                return MSMinimalMedia.determine_min_media(
                    model_s_,
                    minimization_method,
                    min_growth,
                    environment,
                    interacting,
                    printing,
                )
            except Exception as e:
                if not skip_bad_media:
                    print(e)
                return None
        members_media = {}
        for model in model_s_:
            try:
                members_media[model.id] = {
                    "media": MSMinimalMedia.determine_min_media(
                        model,
                        minimization_method,
                        min_growth,
                        environment,
                        interacting,
                        printing,
                    )[0]
                }
                continue
            except Exception as e:
                if skip_bad_media:
                    continue
                else:
                    print(e)
        # print(members_media)
        if com_model is None:
            return members_media
    else:
        return com_media, media_sol
    return {"community_media": com_media, "members": members_media}


def _sigfig_check(value, sigfigs, default):
    if str(value) in ["inf", "nan"]:
        value = ""
    if FBAHelper.isnumber(value):
        return sigfig.round(value, sigfigs)
    else:
        return default


def nanFilter(value, string=True):
    if isinstance(value, str) or value is None:
        if string:
            return value
        else:
            return nan
    if any([value < 0, value > 1e5]):
        return "" if string else nan
    return value
