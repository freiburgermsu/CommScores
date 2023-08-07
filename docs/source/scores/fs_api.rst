Functional Similarity (FS)
--------------------------------------

----------------------
fs()
----------------------

A list of compounds are ascribed a value to construct a media:

.. code-block:: python

 bioplate.add_base(compounds, value)

- *compounds* ``list``: The compounds that will constitute a media.
- *value* ``float``: The value that will be assigned to each compound.


----------------------
Accessible content
----------------------

The ``BilevelPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *plates* ``dict``: The collection of plates (``values``) for all plate IDs (``keys``), will is updated with simulation results from each well media.
