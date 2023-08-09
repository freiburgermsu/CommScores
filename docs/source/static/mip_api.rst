Metabolic Interaction Potential (MIP)
----------------------------------------------------------------------------

----------------------
mip()
----------------------

The metabolic capacity for a member to support the other member is quantified:

.. code-block:: python

 CommScores.mip(member_models: Iterable, com_model=None, min_growth=0.1, interacting_media_dict=None,
                noninteracting_media_dict=None, environment=None, printing=False, compatibilized=False,
                costless=False, multi_output=False, skip_bad_media=False)

- *member_models* ``list|set``: The models from which the MSModelUtil objects and the community model can be constructed.
- *com_model* ``cobrakbase.core.kbasefba.fbamodel.FBAModel``: The community model that will be simulated, which can alternatively be constructed from the *member_models* or *model_utils* arguments when *comm_model* remains ``None``.
- *min_growth* ``float``: The model minimum growth in the specified media, which is used as the lower boundary when computing the minimal media for each model.
- *interacting_media_dict* & *noninteracting_media_dict* ``dict``: The media dictionaries of exchange IDs and fluxes when the members are allowed and disallowed to interact through cross-feeding, respectively.
- *environment* ``dict|cobrakbase.core.kbasebiochem.media.Media``: The media environment in which the member models will be simulated.
- *printing* ``bool``: specifies whether progress and errors during computation are printed.
- *compatibilized* ``bool``: specifies whether the member models have been aligned with ModelSEED conventions.
- *costless* ``bool``: specifies whether the costless MIP subscore is computed and reported, which is the number of cross-fed compounds that are also costlessly excreted.
- *multi_output* ``bool``: specifies whether the costless MIP is reported in addition to, or in lieu of, the standard MIP score.
- *skip_bad_media* ``bool``: specifies whether media in which the members do not grow are skipped without error or throw an error.

**Returns** ``list|dict``: A list is returned unless *multi_output* is ``False``, in which case a dictionary of the costlessly cross-fed compounds donated from each model is returned. The first element of the returned list is a dictionary of cross-fed compounds donated from each model. A second list element is added when *costless* is ``True``, which is the set of costlessly cross-fed compounds donated from each model.