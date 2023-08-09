Biosynthetic Support Score (BSS)
----------------------------------------------------------------------------

----------------------
bss()
----------------------

The metabolic capacity for a member to support the other member is quantified:

.. code-block:: python

 CommScores.bss(member_models:Iterable=None, model_utils:Iterable=None, environments=None, minMedia=None, skip_bad_media=False)

- *member_models* ``list|set``: The models from which the MSModelUtil objects and the community model can be constructed.
- *model_utils* ``list|set``: The MSModeulUtils objects of the member models that can be provided as an alternative to *member_models* for quicker execution time or to conserve additional constraints in the models.
- *environments* ``list<dict|cobrakbase.core.kbasebiochem.media.Media>``: The media environments in which the member models will be simulated.
- *minMedia* ``dict``: The minimal media of the members, with the structure of `< member ID >: {"media": {< exchange ID> : < flux >}}`.
- *skip_bad_media* ``bool``: specifies whether media in which the members do not grow are skipped without error or throw an error.

**Returns** ``dict``: Keys of each directional interaction -- Model1 supporting Model2 and visa versa -- and values of tuples that contain the compounds that are used to support the other member and the BSS score.