Growth Yield Difference (GYD)
--------------------------------------

----------------------
gyd()
----------------------

The relative difference in the maximal growth yields between two members in a given environment is calculated:

.. code-block:: python

 CommScores.gyd(member_models:Iterable=None, model_utils:Iterable=None, environment=None, coculture_growth=False,
                community=None, anme_comm=False)

- *member_models* ``list|set``: The models from which pairwise combinations will be created and then evaluated.
- *model_utils* ``list|set``: The MSModeulUtils objects of the member models that can be provided as an alternative to *member_models* and offers quicker execution time.
- *environment* ``dict|cobrakbase.core.kbasebiochem.media.Media``: The media environment in which the member models will be simulated.
- *coculture_growth* ``bool``: specifies whether the evaluated growth yields are for their growths in the community or as isolates.
- *community* ``cobrakbase.core.kbasefba.fbamodel.FBAModel``: The community that is used when *coculture_growth* is ``True`` and the member growths within the community must be determined. A fresh Community model can be constructed from the *member_models* and *model_utils* when *community* remains ``None``.
- *anme_comm* ``bool``: specifies whether an environment is parameterized to the models, which may be undesirable for some communities that fail to growth in isolation, such as syntrophic ANME members.

**Returns** ``dict``: A tuple of the GYD score followed by each of the two member growths, with keys of the model IDs delimited by " ++ ".