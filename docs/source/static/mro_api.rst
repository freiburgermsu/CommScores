Metabolic Resource Overlap (MRO)
----------------------------------------------------------------------------

----------------------
mro()
----------------------

The proportion of nutritional requirements for each member that are also required by the other member:

.. code-block:: python

 CommScores.mro(member_models:Iterable=None, mem_media:dict=None, min_growth=0.1, media_dict=None,
                raw_content=False, environment=None, skip_bad_media=False, printing=False, compatibilized=False)

- *member_models* ``list|set``: The models from which the MSModelUtil objects and the community model can be constructed.
- *mem_media* ``dict``: The minimal media of each member, which follows the structure: `"members": {<member ID>: {"media": {< exchange ID> : < flux >}}}` or `<member ID>: {"media": {< exchange ID> : < flux >}}`.
- *min_growth* ``float``: The model minimum growth in the specified media, which is used as the lower boundary when computing the minimal media for each model.
- *media_dict* ``dict``: A full dictionary of minimal media that can be parsed to acquire the member minimal media where *mem_media* is not defined.
- *raw_content* ``bool``: specifies whether the metabolites under competition are returned with the MRO score as a tuple, otherwise the competed metabolites are returned in a separate key of the returned dictionary. 
- *environment* ``dict|cobrakbase.core.kbasebiochem.media.Media``: The media environment in which the member models will be simulated.
- *skip_bad_media* ``bool``: specifies whether media in which the members do not grow are skipped without error or throw an error.
- *printing* ``bool``: specifies whether progress and errors during computation are printed.
- *compatibilized* ``bool``: specifies whether the member models have been aligned with ModelSEED conventions.

**Returns** ``dict``: The output when *raw_content* is ``True`` is a tuple of the competed metabolites and the minimal media of the respective member for each pairwise direction. The output for *raw_content* of ``True`` is a tuple of the MRO score, the number of competed compounds, and the number of nutritionally required compounds for each pairwise direction, with the set of competed metabolites stored in a separate key.