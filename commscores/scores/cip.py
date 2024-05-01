from itertools import chain

from modelseedpy.core.msmodelutl import MSModelUtil


def cip(modelutils=None, member_models=None):  # costless interaction potential
    if not modelutils:
        modelutils = {MSModelUtil(model) for model in member_models}
    costless_mets = set(
        chain.from_iterable([modelutil.costless_excreta() for modelutil in modelutils])
    )
    return costless_mets, len(costless_mets)
