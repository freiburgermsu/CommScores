Costless Interaction Potential (CIP)
--------------------------------------

----------------------
cip()
----------------------

The set of costless excretions between two evaluated members:

.. code-block:: python

 CommScores.cip(modelutils=None, member_models=None)

- *modelutils* ``list|set``: The MSModeulUtils objects of the member models that can be provided as a quicker alternative to *member_models*.
- *member_models* ``list|set``: The models from which MSModelUtil objects will combinations will be created and then evaluated.

**Returns** ``tuple``: A tuple of the costless excreta among the two members followed by the number of these compounds, which is the CIP score.