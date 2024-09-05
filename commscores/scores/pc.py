from mscommunity.mscommunity import MSCommunity
from modelseedpy.core.msmodelutl import MSModelUtil
from numpy import array

# allows to singular execution of this script, besides loading CommScores as an entire package
import sys
from pathlib import Path
if __name__ == "__main__" and (__package__ is None or __package__ == ''):
    parent_dir = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(parent_dir))
    from commscoresutil import CommScoresUtil
else:   from ..commscoresutil import CommScoresUtil


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
    pfba=True
):
    assert member_models or modelutils or community, "Either member_models or modelutils or community must be defined."
    member_models = member_models or [mem.model for mem in modelutils] or community.members
    if com_model is None:
        member_models, com_model = CommScoresUtil._load_models(member_models, None, not compatibilized, printing=False)
    community = community or MSCommunity(com_model, member_models)
    if comm_sol is None:
        abundances = community.predict_abundances(environment, False)
        comm_sol = community.solution
    model_utils = modelutils or [MSModelUtil(mem, True) for mem in member_models]
    modelutils = []
    for mem in model_utils:
        mem.add_medium(environment)
        modelutils.append(mem)
    if isolate_growths is None:
        isolate_growths = {mem.id: mem.model.slim_optimize() for mem in modelutils}
    pc_score = comm_sol.objective_value / sum(list(isolate_growths.values()))
    if not comm_effects:
        return pc_score

    # compute the community growth effects of each member
    comm_member_growths, comm_growth_effect = {}, {}
    for mem in community.members:
        comm_member_growths[mem.id] = comm_sol.fluxes[mem.primary_biomass.id]
        comm_growth_effect[mem.id] = CommScoresUtil.nanFilter(comm_member_growths[mem.id] / isolate_growths[mem.id])
    
    growth_diffs = array([CommScoresUtil.nanFilter(x, False) for x in list(comm_growth_effect.values())])
    th_pos, th_neg = 1 + interaction_threshold, 1 - interaction_threshold
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
    return (pc_score, comm_growth_effect, comm_member_growths, bit)
