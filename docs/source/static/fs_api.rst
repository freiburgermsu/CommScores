Functional Similarity (FS)
--------------------------------------

----------------------
fs()
----------------------

The Jiccard Index between two organism's genomes is determined to assess their functional similarity:

.. code-block:: python

 CommScores.fs(models:Iterable=None, kbase_object=None, cobrakbase_repo_path:str=None,
               kbase_token_path:str=None, annotated_genomes:dict=None, printing=False)

- *models* ``list|set``: The models from which pairwise combinations will be created and then evaluated.
- *kbase_object* ``cobrakbase.kbaseapi.KBaseAPI``: The KBase API object that allows the corresponding genomes for each model to acquired.
- *cobrakbase_repo_path* ``str``: The path to the COBRA-KBase GitHub repository, from which a *kbase_object* object will be created where it is not provided and where the genomes are not explicitly given by *annotated_genomes*.
- *kbase_token_path* ``str``: The path to a kbase user token, which is necessary to access and acquire content from KBase.
- *annotated_genomes* ``dict``: The collection of annotated genomes that will be compared, as an alternative to acquiring the model genomes via *kbase_object*.
- *printing* ``bool``: specifies whether progress and auxillary information is printed during the anaysis.

**Returns** ``dict``: The FS score, as the Jiccard Index of the gene ontologies between each member's genome, with keys of the two member's genome IDs delimited by " ++ ".