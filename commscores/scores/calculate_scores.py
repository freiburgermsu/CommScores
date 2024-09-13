from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core.msmodelutl import MSModelUtil

from mscommunity.commhelper import build_from_species_models
from mscommunity.mscommsim import MSCommunity

from multiprocess import current_process
from collections.abc import Iterable
from pandas import concat, Series
import sigfig
import signal

import sys
from pathlib import Path
# if __name__ == "__main__" and (__package__ is None or __package__ == ''):
parent_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(parent_dir))
from commscoresutil import CommScoresUtil
from logger import logger
from scores.bss import bss
from scores.cip import cip
from scores.fs import fs
from scores.gyd import gyd
from scores.mip import mip
from scores.mro import mro
from scores.pc import pc
# else:
#     from ..logger import logger
#     from ..commscoresutil import CommScoresUtil 
#     from .bss import bss
#     from .cip import cip
#     from .fs import fs
#     from .gyd import gyd
#     from .mip import mip
#     from .mro import mro
#     from .pc import pc

def timeout_handler(signum, frame):
    raise TimeoutError("Operation took too long!")

def _load(model, kbase_obj):
    model_str = model
    if not isinstance(model, Iterable):   model = kbase_obj.get_from_ws(model)
    else:    model = kbase_obj.get_from_ws(*model)
    return model, model_str


def calculate_scores(pairs, member_media=None, environments=None, annotated_genomes=True, lazy_load=False,
                     kbase_obj=None, costless=True, check_models=True, print_progress=False):    
    # process the arguments
    if isinstance(pairs, list):  (pairs, member_media, environments, annotated_genomes, lazy_load, kbase_obj) = pairs
    series, mets = [], []
    if isinstance(environments, (list, tuple)) and hasattr(environments[0], "name"):
        environments = {m.name: FBAHelper.convert_kbase_media(m, 1000) for m in environments}
    elif isinstance(environments, (list, tuple)) and not isinstance(environments[0], dict):
        environments = {f"media{i}": m for i, m in enumerate(environments)}
    # print(environments)
        
    # compute the scores
    pid = current_process().name
    member_media = member_media or {}
    print(member_media)
    model_utils = {}
    count = 0
    for model1, models in pairs.items():
        # load and process model1
        print(model1, models)
        if model1.id == "":    model1.id = "model1"
        model1_str = model1.id
        if lazy_load:    model1, model1_str = _load(model1, kbase_obj)
        # member_media[model1.id] = member_media.get(model1.id, CommScoresUtil._get_media(model_s_=model1))
        if member_media[model1.id]["media"] is None:   continue
        if model1.id not in model_utils:    model_utils[model1.id] = MSModelUtil(model1, True)
        for model_index, model2 in enumerate(models):
            # load and process model2
            if model2.id == "":     model2.id = "model2"
            model2_str = model2.id
            if lazy_load:     model2, model2_str = _load(model2, kbase_obj)
            # member_media[model2.id] = member_media.get(model2.id, CommScoresUtil._get_media(model_s_=model2))
            if member_media[model2.id]["media"] is None:   continue
            if model2.id not in model_utils:   model_utils[model2.id] = MSModelUtil(model2, True)
            
            # define group the model1 and model2 pair
            grouping = [model1, model2]
            grouping_utils = [model_utils[model1.id], model_utils[model2.id]]
            modelIDs = [model.id for model in grouping]
            
            ## construct a community model
            comm_model = build_from_species_models(grouping)
            community = MSCommunity(comm_model, ids=modelIDs)
            comm_sol = comm_model.optimize()
            print(f"{pid}~~{count}\t{modelIDs}\t{type(community.util.model.solver)}\t{comm_sol.objective_value}")
            
            # test every given environment
            for environName, environ in environments.items():
                if print_progress:   print(f"\tEnvironment\t{environName}", end="\t")
                ## check that the models grow in the environment
                if check_models:
                    CommScoresUtil._check_model(model_utils[model1.id], environ, model1_str)
                    CommScoresUtil._check_model(model_utils[model2.id], environ, model2_str)
                ## initiate the KBase output
                report_dic = {f"model{i+1}": modelID for i, modelID in enumerate(modelIDs)}
                ### the mono- and co-cultural growths are determined in the environmental media
                groupUtils = [model_utils[model1.id], model_utils[model2.id], community.util]
                g1, g2, comm = CommScoresUtil._determine_growths(groupUtils, environ, 5)
                isolate_growths = {model1.id: g1, model2.id: g2}
                abundances = community.predict_abundances(environ, True, 10)
                if abundances is None:
                    print(f"The {community.id} failed to compute abundances")
                    continue
                print(abundances)
                coculture_growths = {memID: abundance * comm for memID, abundance in abundances.items()}
                report_dic.update({"media": environName, "monoculture growth model1": g1, "monoculture growth model2": g2})
                report_dic.update({f"coculture growth model{modelIDs.index(memID)+1}": growth for memID, growth in coculture_growths.items()})
                report_dic.update({"community growth": comm})
                
                ### add the MRO score
                mro_values = mro(grouping, member_media, raw_content=True, environment=environ)
                report_dic.update(
                    {f"MRO_model{modelIDs.index(models_string.split('--')[0])+1}": f"{100*len(intersection)/len(memMedia):.3f}% ({len(intersection)}/{len(memMedia)})"
                        for models_string, (intersection, memMedia) in mro_values.items()})
                mets.append({"MRO metabolites": list(mro_values.values())[0][0]})
                if print_progress:   print("MRO done", end="\t")
                
                ### add the CIP score
                cip_mets, cipVal = cip(modelutils=[model_utils[mem.id] for mem in grouping])
                report_dic.update({"CIP": cipVal})
                mets[-1].update({"CIP metabolites": list(cip_mets)})
                if print_progress:
                    print("CIP done", end="\t")
                
                ### add the MIP score
                mip_mets = mip(grouping, comm_model, 0.1, None, None, environ, print_progress, True, costless, cip_mets)
                if mip_mets is not None:
                    report_dic.update(
                        {f"MIP_model{modelIDs.index(models_name)+1}": str(len(received))
                         for models_name, received in mip_mets.items() if models_name != "costless"}
                    )
                    mets[-1].update({"MIP model1 metabolites": list(mip_mets.values())[0], "MIP model2 metabolites": list(mip_mets.values())[1]})
                    if costless:
                        for models_name, received in mip_mets["costless"].items():
                            report_dic[f"MIP_model{modelIDs.index(models_name)+1} (costless)"] = (
                                report_dic[f"MIP_model{modelIDs.index(models_name)+1}"] + f" ({len(received)})")
                            del report_dic[f"MIP_model{modelIDs.index(models_name)+1}"]
                        if print_progress:
                            print("costless_MIP  done", end="\t")
                else:
                    report_dic.update({"MIP_model1 (costless)": "", "MIP_model2 (costless)": ""})
                    mets[-1].update({"MIP model1 metabolites": [None], "MIP model2 metabolites": [None]})
                if print_progress:  print("MIP done", end="\t")
                
                ### add the BSS score
                bss_values = bss(grouping, grouping_utils, environments, member_media)
                report_dic.update(
                    {f"BSS_model{modelIDs.index(name.split(' supporting ')[0])+1}": f"{CommScoresUtil._sigfig_check(100*val, 5, '')}%"
                     for name, (mets, val) in bss_values.items()}
                )
                mets[-1].update(
                    {"BSS model1 metabolites": [met_set for met_set, val in bss_values.values()][0],
                     "BSS model2 metabolites": [met_set for met_set, val in bss_values.values()][1]}
                )
                if print_progress:
                    print("BSS done", end="\t")
                
                ### add the PC score
                pc_values = pc(grouping, grouping_utils, comm_model, isolate_growths, comm_sol, environ, True, community)
                if pc_values is not None:
                    report_dic.update(
                        {"PC_model1": CommScoresUtil._sigfig_check(list(pc_values[1].values())[0], 5, ""),
                         "PC_model2": CommScoresUtil._sigfig_check(list(pc_values[1].values())[1], 5, ""),
                         "PC_comm": CommScoresUtil._sigfig_check(pc_values[0], 5, ""),
                         "BIT": pc_values[3]}
                    )
                    if print_progress:  print("PC  done\tBIT done", end="\t")
                
                ### add the GYD score
                if g1 > 0 and g1 > 0:
                    report_dic.update({"GYD1": CommScoresUtil._sigfig_check(abs(g1-g2)/g1, 5, ""),
                                       "GYD2": CommScoresUtil._sigfig_check(abs(g2-g1)/g2, 5, "")})
                    if print_progress:
                        print("GYD done\t\t", end="\t" if annotated_genomes else "\n")
                
                ### add the FS score
                if isinstance(annotated_genomes, (list, dict, tuple)) or kbase_obj is not None and annotated_genomes is True:
                    print(kbase_obj, annotated_genomes)
                    fs_values = list(fs(grouping, kbase_obj, annotated_genomes=annotated_genomes).values())[0]
                    print(len(fs_values[0]) if fs_values[0] is not None else "NaN", fs_values[1])
                    report_dic.update({"FS": sigfig.round(fs_values[1], 5)})
                    if fs_values is not None:  mets[-1].update({"FS features": fs_values[0]})
                    if print_progress:  print("FS done\t\t")
                series.append(Series(report_dic))
                count += 1
                if count % 100 == 0:
                    df = concat(series, axis=1).T
                    df.to_csv(f"{pid}_CommScores.csv")
    # return a pandas Series, which can be easily aggregated with other results into a DataFrame
    return series, mets
