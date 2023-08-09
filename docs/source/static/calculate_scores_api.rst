Calculate Scores
----------------------------------------------------------------------------

----------------------
calculate_scores()
----------------------

A high-level function that simulates all of the scores for a range of models and environments, and returns a DataFrame and list of metabolites for concise interpretation:

.. code-block:: python

 CommScores.calculate_scores(pairs, models_media=None, environments=None, annotated_genomes=True, lazy_load=False,
                             kbase_obj=None, cip_score=True, costless=True, skip_bad_media=False, anme_comm=False,
                             print_progress=False)

- *pairs* ``list|dict``: A provided list denotes all of the argument inputs from parallelization -- ``pairs``, ``models_media``, ``environments``, ``annotated_genomes``, ``lazy_load``, ``kbase_obj`` -- and is accordingly unpacked. A provided dictionary specifies all of the models that are paired with each given model, as a concise means of simulating only specified pairs.
- *models_media* ``dict``: The minimal media of each member, which follows the structure: `<member ID>: {"media": {< exchange ID> : < flux >}}`.
- *environments* ``list<dict|cobrakbase.core.kbasebiochem.media.Media>``: The media environments in which the member models will be simulated.
- *annotated_genomes* ``dict``: The collection of annotated genomes that will be compared, as an alternative to acquiring the model genomes via *kbase_object*.



- *modelutils* ``list|set``: The MSModeulUtils objects of the member models that can be provided as an alternative to *member_models* for quicker execution time or to conserve additional constraints in the models.
- *com_model* ``cobrakbase.core.kbasefba.fbamodel.FBAModel``: The community model that will be simulated, which can alternatively be constructed from the *member_models* or *model_utils* arguments when *comm_model* remains ``None``.
- *isolate_growth* ``dict``: The maximal growth of the member isolates, keyed by their IDs, which can be directly provided through this argument to circumvent external simulation.
- *com_sol* ``cobra.core.solution.Solution``: The flux solution of the community at maximal growth, which is computed when *com_sol* is left as ``None``.
- *environment* ``dict|cobrakbase.core.kbasebiochem.media.Media``: The media environment in which the member models will be simulated.
- *comm_effects* ``bool``: specifies whether the biological interaction type will be approximated based on the PC score for each member, otherwise the community PC score is returned.
- *community* ``modelseedpy.community.mscommunity.MSCommunity``: The MSCommunity object that encapsulates the *com_model*.
- *interaction_threshold* ``float``: The deviation threshold from a PC score of 1 is that used to identify a positive or negative effect on a given member, which is used to compute the BIT score.
- *compatibilized* ``bool``: specifies whether the member models have been compatibilized to ModelSEED conventions, which is necessary to appropriately simulate the community growth and behaviors.

**Returns** ``float|tuple``: The function returns the community PC score when *comm_effects* is ``False``. A *comm_effects* of ``True`` returns instead a tuple of the community PC score, a dictionary with the PC score of each member, and the string BIT score prediction based on the .