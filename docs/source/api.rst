CommScores API
________________________________________________________________________

The detailed documentation of all user-operable classes and functions in the CommScores package.


Static Methods
---------

The ``core`` sub-library permits parsing and manipulating metabolic models::

 from modelseedpy.core import *   
 

.. toctree::
   :includehidden:
   
   static/fs_api
   static/gyd_api
   static/cip_api
   static/mro_api
   static/mip_api
   static/bss_api
   static/pc_api
   static/calculate_scores_api
   static/kbase_report_api
   
   
Class object
------------
   
The aforementioned static methods can also be used synergistically via the ``CommScores`` class object, which allows the parameterization of a model once and the the excution of each of score on the single set of parameters. The score results are subsequently stored in the Class object for succinct reference.
