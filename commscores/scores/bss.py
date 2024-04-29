import re
from itertools import combinations
from typing import Iterable

from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core.msmodelutl import MSModelUtil

from ..utils import _get_media

rm_comp = FBAHelper.remove_compartment


def bss(
    member_models: Iterable = None,
    model_utils: Iterable = None,
    environments=None,
    minMedia=None,
    skip_bad_media=False,
):
    def compute_score(minMedia, environment=None, index=0):
        minMedia = minMedia or _get_media(
            model_s_=[modelUtil.model for modelUtil in model_utils],
            environment=environment,
            skip_bad_media=skip_bad_media,
        )
        model1_media = set(
            [
                re.sub(r"(\_\w\d+$)", "", rxnID.replace("EX_", ""))
                for rxnID in minMedia[model1_util.id]["media"].keys()
            ]
        )
        model2_media = set(
            [
                re.sub(r"(\_\w\d+$)", "", rxnID.replace("EX_", ""))
                for rxnID in minMedia[model2_util.id]["media"].keys()
            ]
        )
        model1_internal = {
            rm_comp(met.id)
            for rxn in model1_util.internal_list()
            for met in rxn.products
        }
        model2_internal = {
            rm_comp(met.id)
            for rxn in model2_util.internal_list()
            for met in rxn.products
        }
        bss_scores[f"{model1_util.id} supporting {model2_util.id} in media{index}"] = (
            model1_internal,
            len(model2_media.intersection(model1_internal)) / len(model2_media),
        )
        bss_scores[f"{model2_util.id} supporting {model1_util.id} in media{index}"] = (
            model2_internal,
            len(model1_media.intersection(model2_internal)) / len(model1_media),
        )

    bss_scores = {}
    for combination in combinations(model_utils or member_models, 2):
        if model_utils is None:
            model1_util = MSModelUtil(combination[0], True)
            model2_util = MSModelUtil(combination[1], True)
            model_utils = [model1_util, model2_util]
        else:
            model1_util = combination[0]
            model2_util = combination[1]
        if environments:
            for index, environment in enumerate(environments):
                compute_score(minMedia, environment, index)
        else:
            compute_score(minMedia)
    return bss_scores
