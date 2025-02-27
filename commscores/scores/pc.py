from mscommunity.mscommsim import MSCommunity
from modelseedpy.core.msmodelutl import MSModelUtil
from numpy import array

# allows to singular execution of this script, besides loading CommScores as an entire package
import sys
from pathlib import Path
# if __name__ == "__main__" and (__package__ is None or __package__ == ''):
parent_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(parent_dir))
from commscoresutil import CommScoresUtil
# else:   from ..commscoresutil import CommScoresUtil


def bit(comm_growth_effect, th_pos, th_neg, comm_member_growths, isolate_growths):
    growth_diffs = array([CommScoresUtil.nanFilter(x, False) for x in list(comm_growth_effect.values())])
    if all(growth_diffs > th_pos):
        bit = "mutualism"
    elif all(growth_diffs < th_neg):
        bit = "competitive"
    elif ((th_pos > growth_diffs) & (growth_diffs > th_neg)).all():
        bit = "neutral"
    elif all(growth_diffs > th_neg) and any(growth_diffs > th_pos):
        bit = "commensalism"
    elif all(growth_diffs < th_pos) and any(growth_diffs < th_neg):
        bit = "amensalism"
    elif any(growth_diffs > th_pos) and any(growth_diffs < th_neg):
        bit = "parasitism"
    else:
        print(
            f"The relative growths {comm_growth_effect} from {comm_member_growths} coculture and"
            f" {isolate_growths} monoculture are not captured."
        )
        bit = ""
    return bit

def pc(
    member_models=None,
    modelutils=None,
    com_model=None,
    isolate_growths=None,
    comm_sol=None,
    environment=None,
    comm_effects=True,
    community=None,
    interaction_threshold=0.1,
    compatibilized=False,
    pfba=True,
    climit=None,
    o2limit=None
):
    assert member_models or modelutils or community, "Either member_models or modelutils or community must be defined."
    member_models = member_models or [mem.model for mem in modelutils] or community.members
    if com_model is None:
        member_models, com_model = CommScoresUtil._load_models(member_models, None, not compatibilized, None, climit, o2limit, False)
    community = community or MSCommunity(com_model, member_models, climit=climit, o2limit=o2limit)
    comm_sol = comm_sol or community.util.model.optimize()
    modelutils = modelutils or [MSModelUtil(mem, True, None, climit, o2limit) for mem in member_models]
    if isolate_growths is None:
        for mem in modelutils:
            mem.add_medium(environment)
        isolate_growths = {mem.id: mem.model.slim_optimize() for mem in modelutils}
    memberGrowth = sum(list(isolate_growths.values()))
    if memberGrowth == 0:  return 
    pc_score = comm_sol.objective_value / memberGrowth
    if not comm_effects:   return pc_score

    # compute the community growth effects of each member
    comm_member_growths, comm_growth_effect = {}, {}
    for mem in community.members:
        comm_member_growths[mem.id] = float(comm_sol.fluxes[mem.primary_biomass.id])
        comm_growth_effect[mem.id] = float(comm_member_growths[mem.id] / isolate_growths[mem.id])
    
    th_pos, th_neg = 1 + interaction_threshold, 1 - interaction_threshold    
    return (pc_score, comm_growth_effect, comm_member_growths, bit(comm_growth_effect, th_pos, th_neg, comm_member_growths, isolate_growths))
