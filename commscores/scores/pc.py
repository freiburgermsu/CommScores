from modelseedpy.community.mscommunity import MSCommunity
from modelseedpy.core.msmodelutl import MSModelUtil
from numpy import array

from ..utils import _load_models, nanFilter


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
):
    assert member_models or modelutils or community, (
        "Members must be defined through either < member_models >"
        " or < modelutils > or < community >."
    )
    member_models = (
        member_models or [mem.model for mem in modelutils] or community.members
    )
    if com_model is None:
        member_models, com_model = _load_models(
            member_models, None, not compatibilized, printing=False
        )
    community = community or MSCommunity(com_model, member_models)
    if comm_sol is None:
        community.util.add_medium(environment)
        comm_sol = community.util.model.optimize()
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

    comm_member_growths = {
        mem.id: comm_sol.fluxes[mem.primary_biomass.id] for mem in community.members
    }
    comm_growth_effect = {
        memID: nanFilter(comm_environ / isolate_growths[memID])
        for memID, comm_environ in comm_member_growths.items()
    }
    growth_diffs = array(
        [nanFilter(x, False) for x in list(comm_growth_effect.values())]
    )
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
