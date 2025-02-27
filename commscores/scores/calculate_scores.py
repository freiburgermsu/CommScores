from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core.msmodelutl import MSModelUtil

from mscommunity.commhelper import build_from_species_models
from mscommunity.mscommsim import MSCommunity

from multiprocess import current_process
from collections.abc import Iterable
from pandas import concat, Series
from os import path
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
from scores.pc import pc, bit
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


def calculate_scores(pairs, member_media=None, environments=None, annotated_genomes=True, lazy_load=False, kbase_obj=None,
                     climit=120, o2limit=120/3, kinCoef=500, costless=True, check_models=True, print_progress=False): 
    # process the arguments
    if isinstance(pairs, list):
        (pairs, member_media, environments, annotated_genomes, lazy_load, kbase_obj, climit, o2limit, kinCoef) = pairs
    series, mets = [], []
    if isinstance(environments, (list, tuple)) and hasattr(environments[0], "name"):
        environments = {m.name: FBAHelper.convert_kbase_media(m, 1000) for m in environments}
    elif isinstance(environments, (list, tuple)) and not isinstance(environments[0], dict):
        environments = {f"media{i}": m for i, m in enumerate(environments)}
    singleclimit = climit/2 if FBAHelper.isnumber(climit) else None
    singleo2limit = o2limit/2 if FBAHelper.isnumber(o2limit) else None
        
    # compute the scores
    pid = current_process().name
    member_media = member_media or {}
    model_utils = {}
    count = 0
    for model1, models in pairs.items():
        # load and process model1
        # print(model1, models)
        if model1.id == "":    model1.id = "model1"
        model1_str = model1.id
        if lazy_load:    model1, model1_str = _load(model1, kbase_obj)
        # member_media[model1.id] = member_media.get(model1.id, CommScoresUtil._get_media(model_s_=model1))
        if member_media[model1.id] is None:  print(f"skipping {model1.id}") ; continue
        if model1.id not in model_utils:   model_utils[model1.id] = MSModelUtil(model1, True, None, singleclimit, singleo2limit)
        for model_index, model2 in enumerate(models):
            # load and process model2
            if model2.id == "":     model2.id = "model2"
            model2_str = model2.id
            if lazy_load:     model2, model2_str = _load(model2, kbase_obj)
            # member_media[model2.id] = member_media.get(model2.id, CommScoresUtil._get_media(model_s_=model2))
            if member_media[model2.id] is None:  print(f"skipping {model1.id}") ; continue
            if model2.id not in model_utils:   model_utils[model2.id] = MSModelUtil(model2, True, None, singleclimit, singleo2limit)
            
            # define group the model1 and model2 pair
            grouping = [model1, model2]
            grouping_utils = [model_utils[model1.id], model_utils[model2.id]]
            modelIDs = [model.id for model in grouping]
            
            ## construct a community model
            comm_model = build_from_species_models(grouping, climit=climit, o2limit=o2limit)
            print("raw model", comm_model.slim_optimize())
            community = MSCommunity(comm_model, ids=modelIDs, kinetic_coeff=kinCoef, climit=climit, o2limit=o2limit)
            print(f"{pid}~~{count}\t{type(community.util.model.solver)}\t{comm_model.slim_optimize()}")
            
            # test every given environment
            groupMedia = {k:v for k,v in member_media.items() if k in modelIDs}
            modelSeries = []
            for environName, environ in environments.items():
                print(f"\t{modelIDs}\tEnvironment\t{environName}", end="\t")
                ## check that the models grow in the environment
                if check_models:
                    CommScoresUtil._check_model(model_utils[model1.id], environ, model1_str)
                    CommScoresUtil._check_model(model_utils[model2.id], environ, model2_str)
                ## initiate the KBase output
                report_dic = {"model1": model1.id, "model2": model2.id}
                ### the mono- and co-cultural growths are determined in the environmental media
                groupUtils = [model_utils[model1.id], model_utils[model2.id], community.util]
                g1, g2, comm = CommScoresUtil._determine_growths(groupUtils, environ, 5)
                isolate_growths = {model1.id: g1, model2.id: g2}
                abundances = community.predict_abundances(environ, True, 10, environName)
                if abundances is None:
                    print(f"FAILED: Abundances of {community.id} were not computed")
                    continue
                # print(abundances)
                coculture_growths = {memID: abundance * comm for memID, abundance in abundances.items()}
                report_dic.update({"media": environName,
                                   "monoculture growth model1": CommScoresUtil._sigfig_check(g1, 5, ""),
                                   "monoculture growth model2": CommScoresUtil._sigfig_check(g2, 5, "")})
                report_dic.update({f"coculture growth model{modelIDs.index(memID)+1}": CommScoresUtil._sigfig_check(growth, 5, "")
                                   for memID, growth in coculture_growths.items()})
                report_dic.update({"community growth": comm})
                
                ### add the MRO score
                mro_values = mro(grouping, groupMedia, raw_content=True, environment=environ)
                if len(list(mro_values.values())[0]) == 2:
                    for commName, (inter, m2Media) in mro_values.items():
                        mroName = f"MRO_model{modelIDs.index(commName.split('--')[0])+1}"
                        mroScore = f"{100*len(inter)/len(m2Media):.3f}% ({len(inter)}/{len(m2Media)})"
                        report_dic.update({mroName: mroScore})
                        mets.append({"MRO metabolites": inter})
                    if print_progress:   print("MRO done", end="\t")
                
                ### add the CIP score
                cip_mets, cipVal = cip(modelutils=[model_utils[mem.id] for mem in grouping])
                report_dic.update({"CIP": cipVal})
                mets[-1].update({"CIP metabolites": list(cip_mets)})
                if print_progress:   print("CIP done", end="\t")
                
                ### add the MIP score
                mip_mets = mip(grouping_utils, community.util, 0.1, None, None, environ, print_progress, True, costless, cip_mets, climit, o2limit)
                if mip_mets is not None:
                    report_dic.update({f"MIP_model{modelIDs.index(models_name)+1}": str(len(received))
                                       for models_name, received in mip_mets.items() if models_name != "costless"})
                    mets[-1].update({"MIP model1 metabolites": list(mip_mets.values())[0], "MIP model2 metabolites": list(mip_mets.values())[1]})
                    if costless:
                        for models_name, received in mip_mets["costless"].items():
                            report_dic[f"MIP_model{modelIDs.index(models_name)+1} (costless)"] = (
                                report_dic[f"MIP_model{modelIDs.index(models_name)+1}"] + f" ({len(received)})")
                            del report_dic[f"MIP_model{modelIDs.index(models_name)+1}"]
                        if print_progress:    print("costless_MIP  done", end="\t")
                else:
                    report_dic.update({"MIP_model1 (costless)": "", "MIP_model2 (costless)": ""})
                    mets[-1].update({"MIP model1 metabolites": [None], "MIP model2 metabolites": [None]})
                if print_progress:  print("MIP done", end="\t")
                
                ### add the BSS score
                bss_values = bss(grouping, grouping_utils, environments, groupMedia)
                if len(list(bss_values.values())[0]) == 2:
                    report_dic.update(
                        {f"BSS_model{modelIDs.index(name.split(' supporting ')[0])+1}": f"{CommScoresUtil._sigfig_check(100*val, 5, '')}%"
                         for name, (mets, val) in bss_values.items()}
                    )
                    mets[-1].update(
                        {"BSS model1 metabolites": [met_set for met_set, val in bss_values.values()][0],
                         "BSS model2 metabolites": [met_set for met_set, val in bss_values.values()][1]}
                    )
                    if print_progress:   print("BSS done", end="\t")
                
                ### add the PC score
                memberGrowth = sum(list(isolate_growths.values()))
                if memberGrowth > 0:
                    comm_growth_effect = {}
                    pc_score = comm / memberGrowth
                    print([mem.primary_biomass.id for mem in community.members])
                    for mem in community.members:
                        if isolate_growths[mem.id] > 0:  comm_growth_effect[mem.id] = coculture_growths[mem.id] / isolate_growths[mem.id]
                        else:     comm_growth_effect[mem.id] = ""
                    pc1, pc2 = [comm_growth_effect[m] for m in modelIDs]
                    report_dic.update({"PC_model1": CommScoresUtil._sigfig_check(pc1, 5, ""),
                                       "PC_model2": CommScoresUtil._sigfig_check(pc2, 5, ""),
                                       "PC_comm": CommScoresUtil._sigfig_check(pc_score, 5, "")})
                    if print_progress:  print("PC  done\t")
                    if not isinstance(pc1, str) and not isinstance(pc2, str):
                        bit_score = bit(comm_growth_effect, 1.1, 0.9, coculture_growths, isolate_growths)
                        report_dic.update({"BIT": bit_score})
                        if print_progress:  print("BIT done", end="\t")
                
                ### add the GYD score
                if g1 > 0:  report_dic.update({"GYD1": CommScoresUtil._sigfig_check((g1-g2)/g1, 5, "")})
                if g2 > 0:  report_dic.update({"GYD2": CommScoresUtil._sigfig_check((g2-g1)/g2, 5, "")})
                if (g1 > 0 or g2 > 0) and print_progress:   print("GYD done\t\t", end="\t" if annotated_genomes else "\n")
                
                ### add the FS score
                if isinstance(annotated_genomes, (list, dict, tuple)) or kbase_obj is not None and annotated_genomes is True:
                    print(kbase_obj, annotated_genomes)
                    fs_values = list(fs(grouping, kbase_obj, annotated_genomes=annotated_genomes).values())[0]
                    print(len(fs_values[0]) if fs_values[0] is not None else "NaN", fs_values[1])
                    report_dic.update({"FS": sigfig.round(fs_values[1], 5)})
                    if fs_values is not None:  mets[-1].update({"FS features": fs_values[0]})
                    if print_progress:  print("FS done\t\t")
                series.append(Series(report_dic))
                modelSeries.append(Series(report_dic))
            df = concat(modelSeries, axis=1).T
            df.to_csv(f"{model1.id}_{model2.id}_CommScores.csv")
            count += 1
    # return a pandas Series, which can be easily aggregated with other results into a DataFrame
    return series, mets
