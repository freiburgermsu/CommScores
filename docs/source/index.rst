Quantifying microbial interactions within microbial communities
________________________________________________________________________

|PyPI version| |Downloads| |License|

.. |Supported Python Versions| image:: https://img.shields.io/pypi/pyversions/commscores)
   :target: https://pypi.org/project/commscores/
   :alt: Python versions

.. |PyPI version| image:: https://img.shields.io/pypi/v/modelseedpy.svg?logo=PyPI&logoColor=brightgreen
   :target: https://pypi.org/project/commscores/
   :alt: PyPI version

.. |Actions Status| image:: https://github.com/freiburgermsu/modelseedpy/workflows/Test%20modelseedpy/badge.svg
   :target: https://github.com/freiburgermsu/commscores/actions
   :alt: Actions Status

.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

.. |Downloads| image:: https://pepy.tech/badge/commscores
   :target: https://pepy.tech/project/commscores
   :alt: Downloads

.. Metabolic modeling is an pivotal method for computational research in synthetic biology and precision medicine. The metabolic models, such as the constrint-based flux balance analysis (FBA) algorithm, are improved with comprehensive datasets that capture more metabolic chemistry in the model and improve the accuracy of simulation predictions. We therefore developed ModelSEEDpy as a comprehensive suite of packages that bootstrap metabolic modeling with the ModelSEED Database (`Seaver et al., 2021 <https://academic.oup.com/nar/article/49/D1/D575/5912569?login=true>`_ ). These packages parse and manipulate (e.g. gapfill missing reactions or calculated chemical properties of metabolites), constrain (with kinetic, thermodynamics, and nutrient uptake), and simulate cobrakbase models (both individual models and communities). This is achieved by standardizing COBRA models through the   ``cobrakbase`` module into a form that is amenable with the KBase/ModelSEED ecosystem. These functionalities are exemplified in `Python Notebooks <https://github.com/ModelSEED/ModelSEEDpy/examples>`_ . Please submit errors, inquiries, or suggestions as `GitHub issues <https://github.com/ModelSEED/ModelSEEDpy/issues>`_ where they can be addressed by our developers.


.. note::

   This project is under active development, and may be subject to losing back-compatibility.

----------------------
Installation
----------------------

**Windows users** must separately install the ``pyeda`` module: 1) download the appropriate wheel for your Python version from `this website <https://www.lfd.uci.edu/~gohlke/pythonlibs/#pyeda>`_ ; and 2) install the wheel through the following commands in a command prompt/powershell console::

 cd path/to/pyeda/wheel
 pip install pyeda_wheel_name.whl

The primary dependency of CommScores is the `Freiburgermsu` fork of the ModelSEEDpy library, which must be installed via ``pip`` in a terminal/command prompt::

 pip install git+https://github.com/freiburgermsu/ModelSEEDpy.git

Every other dependency should automatically install from the command::

 pip install commscores 

.. toctree::
   :hidden:
   
   scores   
   api
