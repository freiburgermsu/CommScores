# Quantifying microbial interactions within microbial communities

[![PyPI version](https://img.shields.io/pypi/v/modelseedpy.svg?logo=PyPI&logoColor=brightgreen)](https://pypi.org/project/commscores/)
[![Downloads](https://pepy.tech/badge/commscores)](https://pepy.tech/project/commscores)
[![Documentation Status](https://readthedocs.org/projects/commscores/badge/?version=latest)](https://commscores.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Microbial communities predicate most biological systems on Earth, yet
the interaction dynamics between community members remains opaque.
Quantitative metrics offer a means of isolating these complex,
multi-dimensional, interactions into single biological dimensions;
although, several ostensibly important biological dimensions evade
existing metrics and extant metrics have moreover not be consolidated
into a single operable package. We therefore developed CommScores as a
comprehensive package for quantifying microbial interaction dimensions
within a microbial community, and thereby elucidating the dynamics that
govern the given community. CommScores leverages
[ModelSEEDpy](https://github.com/ModelSEED/ModelSEEDpy) and
[COBRApy](https://github.com/opencobra/cobrapy) packages for metabolic
modeling, and the [COBRA-KBase](https://github.com/fliu/cobrakbase)
package for acquiring genomic information from KBase. CommScores should
accelerate fundamental discoveries in microbial ecology and the rational
design of microbial communities for diverse applications in medicine,
ecology, and industry.

::: note
::: title
Note
:::

This project is under active development, and may be subject to losing
back-compatibility.
:::

## Installation

CommScores and all of its dependencies should automatically install when
it is installed from PyPI:

    pip install commscores
