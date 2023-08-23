Report generation
----------------------------------------------------------------------------

----------------------
report_generation()
----------------------

The highest-level function in CommScores that executes **calculate_scores**, possibly with parallelization, and returns a Pandas Dataframe and all involved metabolites for succinct analysis:

.. code-block:: python

 CommScores.calculate_scores(all_models:iter=None, pairs:dict=None, mem_media:dict=None, pair_limit:int=None, exclude_pairs:list=None, 
                             kbase_obj=None, annotated_genomes:dict=True, see_media=True, environments:iter=None, pool_size:int=None, 
                             cip_score=True, costless=True, skip_bad_media=False, anme_comm=False, print_progress=False)

- *all_models* ``list``: All of the models that will be simulated. This can be a list of models, where all combinations of the models will be simulated. This argument can alternatively be a list of lists of models, where only pairs between the various lists are examined. The set of all models can also be deduced from *pairs* if that argument is provided instead of *all_models*.
- *pairs* ``dict``: A dictionary of the models that are paired with each given model, as a concise means of simulating only specified pairs.
- *mem_media* ``dict``: The minimal media of each member, which follows the structure: `<member ID>: {"media": {< exchange ID> : < flux >}}`.
- *pair_limit* ``int``: Defines a limit for the number of pairs that are examined for a given system, where a limit that is less than the total number of combinations randomly samples the possible pairwise combinations.
- *excluded_pairs* ``list|set``: The collection of model pairs that will not be simulated, even if they can be constructed from the provided model sets.
- *kbase_obj* ``cobrakbase.kbaseapi.KBaseAPI``: The KBase API object that allows the corresponding genomes for each model to acquired.
- *annotated_genomes* ``dict``: The collection of annotated genomes that will be compared, as an alternative to acquiring the model genomes via *kbase_object*.
- *see_media* ``bool``: specifies whether the computed minimal media for the member models will be printed, which can facilitate storing these media and providing them through the `mem_media` argument in future simulations to save computation time.
- *environments* ``list<dict|cobrakbase.core.kbasebiochem.media.Media>``: The media environments in which the member models will be simulated.
- *pool_size* ``int``: Defines the number of instances that will be employed to parallelize the computation, where parallelization only occurs when this argument is not ``None``.
- *cip_score* ``bool``: specifies whether the CIP score will be computed and reported with the other scores.
- *costless* ``bool``: specifies whether the costless MIP subscore is computed and reported, which is the number of cross-fed compounds that are also costlessly excreted.
- *skip_bad_media* ``bool``: specifies whether media in which the members do not grow are skipped without error or throw an error.
- *anme_comm* ``bool``: specifies whether an environment is parameterized to the models, which may be undesirable for some communities that fail to growth in isolation, such as syntrophic ANME members.
- *print_progress* ``bool``: specifies whether progress and auxillary information is printed with each loop over all pair and environment combinations.

**Returns** ``tuple``: The first tuple element is a Pandas DataFrame object where each row presents all of the score values for a single pair in a single environment. The second tuple element is the list of dictionaries that detail the metabolites that are involved with each respective score. These outputs can be directly fed into the **commscores_report** function to construct an HTML report of the results.