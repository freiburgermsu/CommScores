from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.fbahelper import FBAHelper
from deepdiff import DeepDiff  # (old, new)
from typing import Iterable
from itertools import chain
import re


class NoMedia:
    def __init__(message):
        print(message)


# allows to singular execution of this script, besides loading CommScores as an entire package
import sys
from pathlib import Path
# if __name__ == "__main__" and (__package__ is None or __package__ == ''):
    # Uses directory of script as starting point
parent_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(parent_dir))
from commscoresutil import CommScoresUtil
# else:   from ..commscoresutil import CommScoresUtil


def mip(member_models: Iterable, com_model=None, min_growth=0.1, interacting_media_dict=None,
        noninteracting_media_dict=None, environment=None, printing=False, compatibilized=False, costless=False, cip_mets=None):
    """Determine the quantity of nutrients that can be potentially sourced through syntrophy"""
    member_models, community = CommScoresUtil._load_models(member_models, com_model, not compatibilized, "MIP_comm", printing)

    ## non-interacting media
    print("Non-interacting community, minimize transporters", end="\t")
    noninteracting_medium, noninteracting_sol = CommScoresUtil._get_media(noninteracting_media_dict, community, None, min_growth, environment, False)
    if noninteracting_medium is None:   raise NoMedia("There is no non-interacting media.")
    if "community_media" in noninteracting_medium:
        noninteracting_medium = noninteracting_medium["community_media"]
    
    ## interacting media
    print("Interacting community, minimize exchanges", end="\t")
    interacting_medium, interacting_sol = CommScoresUtil._get_media(interacting_media_dict, community, None, min_growth, environment, True)
    if interacting_medium is None:      raise NoMedia("There is no Interacting media.")
    if "community_media" in interacting_medium:
        interacting_medium = interacting_medium["community_media"]
    
    ## assess the media difference, the interacting relative to the non-interacting
    directionalMIP = {mem.id: [] for mem in member_models}
    interact_diff = DeepDiff(noninteracting_medium, interacting_medium)
    if "dictionary_item_removed" not in interact_diff:
        print(f"The {community.id} model has this media difference: {interact_diff}")
        directionalMIP.update({"costless": {memName: [] for memName in directionalMIP}})
        return directionalMIP
    cross_fed_exIDs = [re.sub("(root\['|'\])", "", x) for x in interact_diff["dictionary_item_removed"]]
    
    # Determine each direction of the MIP score interactions
    comm_util = MSModelUtil(community)
    cross_fed_metIDs = [ex.replace("EX_", "").replace("_e0", "") for ex in cross_fed_exIDs]
    cross_fed_copy = cross_fed_metIDs[:]
    for rxn in comm_util.transport_list():
        metIDs = list(set([met.id.split("_")[0] for met in rxn.reactants]) & set([met.id.split("_")[0] for met in rxn.products]))
        if len(metIDs) == 1:   metID = metIDs[0]
        else:
            ## not select protons
            if "cpd00067" in metIDs: metIDs.remove("cpd00067")
            metID = metIDs[0]
        if metID not in cross_fed_metIDs:   continue
        rxn_index = FBAHelper.compartment_index(rxn.id.split("_")[-1])
        if rxn_index == 0:   continue
        mets = [met for met in rxn.metabolites if met.id == f"{metID}_c{rxn_index}"]
        # if mets == []:
        #     print(f"The {metID}_c{rxn_index} is missing in {rxn.reaction}.")
        #     continue
        rxn_model = member_models[rxn_index - 1]
        ## Adding syntrophic exchanges from the donor/production perspective
        metStoich = rxn.metabolites[mets[0]]
        commFlux = interacting_sol.fluxes[rxn.id]
        if metStoich > 0 and commFlux > 0 or metStoich < 0 and commFlux < 0:
            directionalMIP[rxn_model.id].append(metID)
            if metID in cross_fed_copy:
                cross_fed_copy.remove(metID)
                continue
        # if printing:  print(f"{mets[0]} in {rxn.id} ({rxn.reaction}) is not assigned a receiving member.")
    if cross_fed_copy != []:   print(f"Missing directions for the {cross_fed_copy} cross-fed metabolites")
    outputs = directionalMIP
    # TODO categorize all of the cross-fed substrates to examine potential associations of specific compounds
    if costless or cip_mets:
        if cip_mets is None:
            modelutils = {MSModelUtil(model) for model in member_models}
            cip_mets = chain.from_iterable([modelutil.costless_excreta() for modelutil in modelutils])
        outputs.update({"costless":{memName: set(received_mets) & set(cip_mets) for memName, received_mets in directionalMIP.items()}})
    return outputs
