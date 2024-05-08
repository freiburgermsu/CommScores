import sigfig
from modelseedpy.community.commhelper import build_from_species_models
from modelseedpy.community.mscommunity import MSCommunity
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core.msmodelutl import MSModelUtil
from multiprocess import current_process
from collections import Iterable

from ..logger import logger
from ..utils import _get_media, _sigfig_check
from .bss import bss
from .cip import cip
from .fs import fs
from .gyd import _check_model, _determine_growths, gyd
from .mip import mip
from .mro import mro
from .pc import pc


def _load(model, kbase_obj):
    model_str = model
    if not isinstance(model, Iterable):   model = kbase_obj.get_from_ws(model)
    else:    model = kbase_obj.get_from_ws(*model)
    return model, model_str


def calculate_scores(pairs, models_media=None, environments=None, annotated_genomes=True,
                     lazy_load=False, kbase_obj=None, cip_score=True, costless=True,
                     skip_bad_media=False, check_models=True, print_progress=False):
    from pandas import Series

    if isinstance(pairs, list):
        (pairs, models_media, environments, annotated_genomes, lazy_load, kbase_obj) = pairs
    series, mets = [], []
    if not isinstance(environments, (list, tuple)):   environments = [environments]
    if isinstance(environments, (list, tuple)) and hasattr(environments[0], "name"):
        environments = {m.name: FBAHelper.convert_kbase_media(m, 1000) for m in environments}
    elif not isinstance(environments, dict):
        environments = {f"media{i}": m for i, m in enumerate(environments)}
    pid = current_process().name
    model_utils = {}
    count = 0
    for model1, models in pairs.items():
        if model1.id == "":    model1.id = "model1"
        if lazy_load:    model1, model1_str = _load(model1, kbase_obj)
        else:    model1_str = model1.id
        if model1.id not in models_media:
            models_media[model1.id] = {
                "media": _get_media(model_s_=model1, skip_bad_media=skip_bad_media)
            }
            if models_media[model1.id] is None:   continue
        if model1.id not in model_utils:
            model_utils[model1.id] = MSModelUtil(model1, True)
        # print(pid, model1)
        for model_index, model2 in enumerate(models):
            if model2.id == "":     model2.id = "model2"
            if lazy_load:     model2, model2_str = _load(model2, kbase_obj)
            else:     model2_str = model2.id
            if model2.id not in models_media:
                models_media[model2.id] = {
                    "media": _get_media(model_s_=model2, skip_bad_media=skip_bad_media)
                }
                if models_media[model2.id] is None:    continue
            if model2.id not in model_utils:   model_utils[model2.id] = MSModelUtil(model2, True)
            grouping = [model1, model2]
            grouping_utils = [model_utils[model1.id], model_utils[model2.id]]
            modelIDs = [model.id for model in grouping]
            comm_model = build_from_species_models(grouping)
            community = MSCommunity(comm_model, ids=modelIDs)
            comm_sol = comm_model.optimize()
            print(f"{pid}~~{count}\t{modelIDs}\t{comm_sol.objective_value}")
            for environName, environ in environments.items():
                if print_progress:
                    print(f"\tEnvironment\t{environName}", end="\t")
                if check_models:  # check that the models grow in the environment
                    _check_model(model_utils[model1.id], environ, model1_str, skip_bad_media)
                    _check_model(model_utils[model2.id], environ, model2_str, skip_bad_media)
                # initiate the KBase output
                report_dic = {f"model{i+1}": modelID for i, modelID in enumerate(modelIDs)}
                # the model growths are determined and the environmental media is parameterized for each of the members
                g1, g2, comm = [
                    _sigfig_check(val, 5, "")
                    for val in _determine_growths(
                        [model_utils[model1.id], model_utils[model2.id], community.util,],  environ
                        )
                ]
                coculture_growths = {memID: abundance * comm
                                     for memID, abundance in community.predict_abundances().items()}
                report_dic.update(
                    {"media": environName, "monoculture growth model1": g1, "monoculture growth model2": g2}
                )
                report_dic.update(
                    {f"coculture growth model{modelIDs.index(memID)+1}": growth
                     for memID, growth in coculture_growths.items()}
                )
                report_dic.update({"community growth": comm})
                # define the MRO content
                mro_values = mro(grouping, models_media, raw_content=True, environment=environ)
                report_dic.update(
                    {f"MRO_model{modelIDs.index(models_string.split('--')[0])+1}": f"{100*len(intersection)/len(memMedia):.3f}% ({len(intersection)}/{len(memMedia)})"
                        for models_string, (intersection, memMedia) in mro_values.items()
                    }
                )
                mets.append({"MRO metabolites": list(mro_values.values())[0][0]})
                if print_progress:   print("MRO done", end="\t")
                # define the CIP content
                if cip_score:
                    cip_values = cip(modelutils=[model_utils[mem.id] for mem in grouping])
                    report_dic.update({"CIP": cip_values[1]})
                    mets[-1].update({"CIP metabolites": list(cip_values[0])})
                    if print_progress:
                        print("CIP done", end="\t")
                # define the MIP content
                mip_values = mip(grouping, comm_model, 0.1, None, None, environ, print_progress,
                                 True, costless, costless, skip_bad_media)
                # print(mip_values)
                if mip_values is not None:
                    report_dic.update(
                        {f"MIP_model{modelIDs.index(models_name)+1}": str(len(received))
                         for models_name, received in mip_values[0].items()}
                    )
                    mets[-1].update(
                        {"MIP model1 metabolites": list(mip_values[0].values())[0],
                         "MIP model2 metabolites": list(mip_values[0].values())[1]}
                    )
                    if costless:
                        for models_name, received in mip_values[1].items():
                            report_dic[f"MIP_model{modelIDs.index(models_name)+1} (costless)"] = (
                                report_dic[f"MIP_model{modelIDs.index(models_name)+1}"]
                                + f" ({len(received)})"
                            )
                            del report_dic[f"MIP_model{modelIDs.index(models_name)+1}"]
                        if print_progress:
                            print("costless_MIP  done", end="\t")
                else:
                    report_dic.update(
                        {"MIP_model1 (costless)": "", "MIP_model2 (costless)": ""}
                    )
                    mets[-1].update({"MIP model1 metabolites": [None], "MIP model2 metabolites": [None]})
                if print_progress:
                    print("MIP done", end="\t")
                # define the BSS content
                bss_values = bss(grouping, grouping_utils, environments, models_media, skip_bad_media)
                report_dic.update(
                    {f"BSS_model{modelIDs.index(name.split(' supporting ')[0])+1}": f"{_sigfig_check(100*val, 5, '')}%"
                     for name, (mets, val) in bss_values.items()}
                )
                mets[-1].update(
                    {
                        "BSS model1 metabolites": [
                            met_set for met_set, val in bss_values.values()
                        ][0],
                        "BSS model2 metabolites": [
                            met_set for met_set, val in bss_values.values()
                        ][1],
                    }
                )
                # mets[-1].update({"bss_mets": list(bss_values[0].values())})
                if print_progress:
                    print("BSS done", end="\t")
                # define the PC content
                pc_values = pc(
                    grouping,
                    grouping_utils,
                    comm_model,
                    None,
                    comm_sol,
                    environ,
                    True,
                    community,
                )
                report_dic.update(
                    {
                        "PC_model1": _sigfig_check(
                            list(pc_values[1].values())[0], 5, ""
                        ),
                        "PC_model2": _sigfig_check(
                            list(pc_values[1].values())[1], 5, ""
                        ),
                        "PC_comm": _sigfig_check(pc_values[0], 5, ""),
                        "BIT": pc_values[3],
                    }
                )
                if print_progress:
                    print("PC  done\tBIT done", end="\t")
                # print([mem.slim_optimize() for mem in grouping])
                # define the GYD content
                logger.debug(
                    list(
                        gyd(
                            grouping,
                            grouping_utils,
                            environ,
                            False,
                            community,
                            check_models,
                        ).values()
                    )
                )
                gyd1, gyd2, g1, g2 = list(
                    gyd(
                        grouping,
                        grouping_utils,
                        environ,
                        False,
                        community,
                        check_models,
                    ).values()
                )[0]
                report_dic.update(
                    {
                        "GYD1": _sigfig_check(gyd1, 5, ""),
                        "GYD2": _sigfig_check(gyd2, 5, ""),
                    }
                )
                if print_progress:
                    print("GYD done\t\t", end="\t" if annotated_genomes else "\n")
                # define the FS content

                if kbase_obj is not None and annotated_genomes:
                    fs_values = list(
                        fs(
                            grouping, kbase_obj, annotated_genomes=annotated_genomes
                        ).values()
                    )[0]
                    print(
                        len(fs_values[0]) if fs_values[0] is not None else "NaN",
                        fs_values[1],
                    )
                    report_dic.update({"FS": sigfig.round(fs_values[1], 5)})
                    if fs_values is not None:
                        mets[-1].update({"FS features": fs_values[0]})
                    if print_progress:
                        print("FS done\t\t")
                # return a pandas Series, which can be easily aggregated with other results into a DataFrame
                series.append(Series(report_dic))
            count += 1
    return series, mets
