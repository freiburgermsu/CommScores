from itertools import combinations
from typing import Iterable

from modelseedpy.community.mscommunity import MSCommunity
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.msmodelutl import MSModelUtil
from numpy import isclose


def _check_model(model_util, media, model_str, skip_bad_media=True):
    default_media = model_util.model.medium
    # print("test")
    if media is not None:
        model_util.add_medium(media)
    obj_val = model_util.model.slim_optimize()
    print(model_util.model.id, obj_val)
    if isclose(obj_val, 0, abs_tol=1e-6) or not FBAHelper.isnumber(obj_val):
        print(f"The {model_str} model is not operational")
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


def gyd(
    member_models: Iterable = None,
    model_utils: Iterable = None,
    environment=None,
    coculture_growth=False,
    community=None,
    check_models=True,
):
    gyds = {}
    for combination in combinations(model_utils or member_models, 2):
        if model_utils is None:
            model1_util = MSModelUtil(combination[0], True)
            model2_util = MSModelUtil(combination[1], True)
            print(
                f"{model1_util.model.id} ++ {model2_util.model.id}",
                model1_util.model.slim_optimize(),
                model2_util.model.slim_optimize(),
            )
            if check_models:
                model1_util = _check_model(model1_util, environment)
                model2_util = _check_model(model2_util, environment)
        else:
            model1_util = combination[0]
            model2_util = combination[1]
        if not coculture_growth:
            G_m1, G_m2 = _determine_growths([model1_util, model2_util], environment)
            G_m1, G_m2 = (
                G_m1 if FBAHelper.isnumber(str(G_m1)) else 0,
                (G_m2 if FBAHelper.isnumber(str(G_m2)) else 0),
            )
        else:
            community = community or MSCommunity(
                member_models=[model1_util.model, model2_util.model],
                ids=[mem.id for mem in member_models],
            )
            community.run_fba()
            member_growths = community.parse_member_growths()
            G_m1, G_m2 = (
                member_growths[model1_util.model.id],
                member_growths[model2_util.model.id],
            )
        if G_m2 <= 0 or G_m1 <= 0:
            gyds[f"{model1_util.model.id} ++ {model2_util.model.id}"] = (
                "",
                "",
                G_m1,
                G_m2,
            )
            continue
        gyds[f"{model1_util.model.id} ++ {model2_util.model.id}"] = (
            abs(G_m1 - G_m2) / G_m1,
            abs(G_m2 - G_m1) / G_m2,
            G_m1,
            G_m2,
        )
    return gyds
