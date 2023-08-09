Potential Cooperation (PC) & Biological Interaction Type (BIT)
----------------------------------------------------------------------------

----------------------
pc()
----------------------

The relative growth of the members between their isolate growth and their growth in the pairwise community, which can optionally generate the member-level PC scores and the BIT score, depending on the specified arguments:

.. code-block:: python

 CommScores.pc(member_models=None, modelutils=None, com_model=None, isolate_growths=None, comm_sol=None,
                environment=None, comm_effects=True, community=None, interaction_threshold=0.1, compatibilized=False)

- *member_models* ``list|set``: The models from which the MSModelUtil objects and the community model can be constructed.
- *modelutils* ``list|set``: The MSModeulUtils objects of the member models that can be provided as an alternative to *member_models* for quicker execution time or to conserve additional constraints in the models.
- *com_model* ``cobrakbase.core.kbasefba.fbamodel.FBAModel``: The community model that will be simulated, which can alternatively be constructed from the *member_models* or *model_utils* arguments when *comm_model* remains ``None``.
- *isolate_growth* ``dict``: The maximal growth of the member isolates, keyed by their IDs, which can be directly provided through this argument to circumvent external simulation.
- *com_sol* ``cobra.core.solution.Solution``: The flux solution of the community at maximal growth, which is computed when *com_sol* is left as ``None``.
- *environment* ``dict|cobrakbase.core.kbasebiochem.media.Media``: The media environment in which the member models will be simulated.
- *comm_effects* ``bool``: specifies whether the biological interaction type will be approximated based on the PC score for each member, otherwise the community PC score is returned.
- *community* ``modelseedpy.community.mscommunity.MSCommunity``: The MSCommunity object that encapsulates the *com_model*.
- *interaction_threshold* ``float``: The deviation threshold from a PC score of 1 is that used to identify a positive or negative effect on a given member, which is used to compute the BIT score.
- *compatibilized* ``bool``: specifies whether the member models have been compatibilized to ModelSEED conventions, which is necessary to appropriately simulate the community growth and behaviors.

**Returns** ``float|tuple``: The function returns the community PC score when *comm_effects* is ``False``. A *comm_effects* of ``True`` returns instead a tuple of the community PC score, a dictionary with the PC score of each member, and the string BIT score prediction based on the member PC scores.