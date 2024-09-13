import re
from itertools import combinations
from typing import Iterable
from math import isclose

from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core.msmodelutl import MSModelUtil

# allows to singular execution of this script, besides loading CommScores as an entire package
import sys
from pathlib import Path
# if __name__ == "__main__" and (__package__ is None or __package__ == ''):
parent_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(parent_dir))
from commscoresutil import CommScoresUtil
# else:   from ..commscoresutil import CommScoresUtil

rm_comp = FBAHelper.remove_compartment

def remove_comp(string):
    return re.sub(r"(\_\w\d+$)", "", string.replace("EX_", ""))


def compute_score(minMedia, model_utils, environ, index=0):
    scores = {}
    model1_util, model2_util = model_utils
    sol_growths = []
    for util in model_utils:
        ogMedia = util.model.medium
        util.add_medium(environ)
        sol_growths.append(util.run_fba(None, True).fluxes[util.biomass_objective])
        util.add_medium(ogMedia)
    min_growth = min(sol_growths)
    # if isclose(0, min_growth):
    #     return {}
    minMedia = minMedia or CommScoresUtil._get_media(
        model_s_=[model1_util.model, model2_util.model], min_growth=min_growth, environment=environ)
    model1_media = set(list(map(remove_comp, list(minMedia[model1_util.id]["media"].keys()))))
    model2_media = set(list(map(remove_comp, list(minMedia[model2_util.id]["media"].keys()))))
    model1_internal = {rm_comp(met.id) for rxn in model1_util.internal_list() for met in rxn.products}
    model2_internal = {rm_comp(met.id) for rxn in model2_util.internal_list() for met in rxn.products}
    if len(model2_media) > 0:
        scores[f"{model1_util.id} supporting {model2_util.id} in media{index}"] = (
            list(model1_internal), len(model2_media & model1_internal) / len(model2_media))
    if len(model1_media) > 0:
        scores[f"{model2_util.id} supporting {model1_util.id} in media{index}"] = (
            list(model2_internal), len(model1_media & model2_internal) / len(model1_media))
    return scores


def bss(
    member_models: Iterable = None,
    model_utils: Iterable = None,
    environment=None,
    minMedia=None,
):
    bss_scores = {}
    for combination in combinations(member_models if model_utils is None else model_utils, 2):
        model1_util = MSModelUtil(combination[0], True) if model_utils is None else combination[0]
        model2_util = MSModelUtil(combination[1], True) if model_utils is None else combination[1]
        comb_utils = [model1_util, model2_util]
        if isinstance(environment, (tuple, list, set)):
            for index, environ in enumerate(environment):
                bss_scores.update(compute_score(minMedia, comb_utils, environ, index))
        else:
            bss_scores.update(compute_score(minMedia, comb_utils, environment))
    return bss_scores
