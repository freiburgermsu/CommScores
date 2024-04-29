import re
from typing import Iterable

from deepdiff import DeepDiff  # (old, new)
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core.msmodelutl import MSModelUtil

from ..utils import _get_media, _load_models
from .cip import cip


def mip(
    member_models: Iterable,
    com_model=None,
    min_growth=0.1,
    interacting_media_dict=None,
    noninteracting_media_dict=None,
    environment=None,
    printing=False,
    compatibilized=False,
    costless=False,
    multi_output=False,
    skip_bad_media=False,
):
    """Determine the quantity of nutrients that can be potentially sourced through syntrophy"""
    member_models, community = _load_models(
        member_models, com_model, not compatibilized, printing=printing
    )
    # determine the interacting and non-interacting media for the specified community  .util.model
    noninteracting_medium, noninteracting_sol = _get_media(
        noninteracting_media_dict,
        community,
        None,
        min_growth,
        environment,
        False,
        skip_bad_media=skip_bad_media,
    )
    if noninteracting_medium is None:
        return None
    if "community_media" in noninteracting_medium:
        noninteracting_medium = noninteracting_medium["community_media"]
    interacting_medium, interacting_sol = _get_media(
        interacting_media_dict,
        community,
        None,
        min_growth,
        environment,
        True,
        skip_bad_media=skip_bad_media,
    )
    if interacting_medium is None:
        return None
    if "community_media" in interacting_medium:
        interacting_medium = interacting_medium["community_media"]
    interact_diff = DeepDiff(noninteracting_medium, interacting_medium)
    if "dictionary_item_removed" not in interact_diff:
        return None
    cross_fed_exIDs = [
        re.sub("(root\['|'\])", "", x) for x in interact_diff["dictionary_item_removed"]
    ]
    # Determine each direction of the MIP score interactions
    comm_util = MSModelUtil(community)
    cross_fed_metIDs = [
        ex.replace("EX_", "").replace("_e0", "") for ex in cross_fed_exIDs
    ]
    cross_fed_copy = cross_fed_metIDs[:]
    directionalMIP = {mem.id: [] for mem in member_models}
    for rxn in comm_util.transport_list():
        # print(rxn.reaction, "\t", [met.id for met in rxn.metabolites if "_e0" in met.id])
        metIDs = list(
            set([met.id.split("_")[0] for met in rxn.reactants]).intersection(
                set([met.id.split("_")[0] for met in rxn.products])
            )
        )
        if len(metIDs) == 1:
            metID = metIDs[0]
        else:
            if "cpd00067" in metIDs:
                metIDs.remove("cpd00067")
            metID = metIDs[0]
        if metID not in cross_fed_metIDs:
            continue
        rxn_index = FBAHelper.compartment_index(rxn.id.split("_")[-1])
        if rxn_index == 0:
            continue
        mets = [met for met in rxn.metabolites if met.id == f"{metID}_c{rxn_index}"]
        if mets == []:
            print(f"The {metID}_c{rxn_index} is missing in {rxn.reaction}.")
            continue
        rxn_model = member_models[rxn_index - 1]
        # comm_trans[metID] = comm_trans.get(f"{metID}_c{rxn_index}", {})
        if (
            rxn.metabolites[mets[0]] > 0
            and interacting_sol.fluxes[rxn.id] > 0
            or rxn.metabolites[mets[0]] < 0
            and interacting_sol.fluxes[rxn.id] < 0
        ):  # donor
            directionalMIP[rxn_model.id].append(metID)
            if metID in cross_fed_copy:
                cross_fed_copy.remove(metID)
                continue
        # if printing:  print(f"{mets[0]} in {rxn.id} ({rxn.reaction}) is not assigned a receiving member.")
    if cross_fed_copy != [] and printing:
        print(f"Missing directions for the {cross_fed_copy} cross-fed metabolites")
    outputs = [directionalMIP]
    # TODO categorize all of the cross-fed substrates to examine potential associations of specific compounds
    if costless:
        costless_mets, numExs = cip(member_models=member_models)
        # print(list(directionalMIP.values()), costless_mets)
        costlessDirectionalMIP = {
            member_name: set(receive_mets).intersection(costless_mets)
            for member_name, receive_mets in directionalMIP.items()
        }
        if not multi_output:
            return costlessDirectionalMIP
        outputs.append(costlessDirectionalMIP)
    return outputs
