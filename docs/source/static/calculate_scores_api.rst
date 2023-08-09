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
- *lazy_load* ``bool``: specifies whether only models that are necessary for the current comparison are loaded, to save RAM for analyses of many models.
- *kbase_object* ``cobrakbase.kbaseapi.KBaseAPI``: The KBase API object that allows the corresponding genomes for each model to acquired.
- *cip_score* ``bool``: specifies whether the CIP score will be computed and reported with the other scores.
- *costless* ``bool``: specifies whether the costless MIP subscore is computed and reported, which is the number of cross-fed compounds that are also costlessly excreted.
- *skip_bad_media* ``bool``: specifies whether media in which the members do not grow are skipped without error or throw an error.
- *anme_comm* ``bool``: specifies whether an environment is parameterized to the models, which may be undesirable for some communities that fail to growth in isolation, such as syntrophic ANME members.
- *print_progress* ``bool``: specifies whether progress and auxillary information is printed with each loop over all pair and environment combinations.

**Returns** ``tuple``: The first tuple element is a list of Pandas Series objects that represent the results of all scores for a single pair in a single environment, which can be seamlessly converted into a Pandas Dataframe via `pandas.concat(series_list, axis=1).T`. The second tuple element is the list of dictionaries that detail the metabolites that are involved with each respective score.