from cobra.io import read_sbml_model
from shutil import rmtree
from glob import glob
import commscores
import re, os


test_models = [read_sbml_model(model) for model in glob("*.xml")]


def test_mro():
    # validate the pure scores
    mro_output = commscores.mro(test_models)
    for key, val in mro_output.items():
        assert isinstance(val[0], list), "the MRO compounds under competition"
        assert isinstance(val[1], dict), "the MRO minimal media is faulty"
        if key == "model1---model2":   assert val == (48.38709677419355, 15, 31), "the MRO computation has changed"
        if key == "model2---model1":   assert val == (71.42857142857143, 15, 21), "the MRO computation has changed"
        if key == "model1---model3":   assert val == (58.064516129032256, 18, 31), "the MRO computation has changed"
        if key == "model3---model1":   assert val == (60.0, 18, 30), "the MRO computation has changed"
        if key == "model1---model4":   assert val == (51.61290322580645, 16, 31), "the MRO computation has changed"
        if key == "model4---model1":   assert val == (69.56521739130434, 16, 23), "the MRO computation has changed"
        if key == "model2---model3":   assert val == (71.42857142857143, 15, 21), "the MRO computation has changed"
        if key == "model3---model2":   assert val == (50.0, 15, 30), "the MRO computation has changed"
        if key == "model2---model4":   assert val == (85.71428571428571, 18, 21), "the MRO computation has changed"
        if key == "model4---model2":   assert val == (78.26086956521739, 18, 23), "the MRO computation has changed"
        if key == "model3---model4":   assert val == (53.333333333333336, 16, 30), "the MRO computation has changed"
        if key == "model4---model3":   assert val == (69.56521739130434, 16, 23), "the MRO computation has changed"

    # validate the metabolite predictions and compounds under competition
    mro_output = commscores.mro(test_models, raw_content=True)
    for key, val in mro_output.items():
        if "mets" not in key:   assert isinstance(val, list), "the MRO mets are faulty"
        else:    assert all([isinstance(ele, float) or isinstance(ele, int) for ele in val]),  "the processed MRO scores are faulty"


def test_mip():
    # validate the predicted cross-fed compounds
    mip_output = commscores.mip(test_models),
    for modelID, val in mip_output.items():
        assert isinstance(val[0], list), "the MIP syntrophic compounds are faulty"
        if modelID == "model1":
            assert val == ['cpd00128', 'cpd00136', 'cpd01080', 'cpd00081', 'cpd00129', 'cpd00048', 'cpd00039', 'cpd00066'], "the MIP computation has changed"
        if modelID == "model2":
            assert val == ['cpd00324', 'cpd00071', 'cpd00048', 'cpd00013'], "the MIP computation has changed"
        if modelID == "model3":
            assert val == ['cpd00100', 'cpd00281', 'cpd00129', 'cpd00363', 'cpd00071', 'cpd00239', 'cpd00065'], "the MIP computation has changed"
        if modelID == "model4":
            assert val == ['cpd00023', 'cpd00324', 'cpd00048', 'cpd00363', 'cpd00081', 'cpd00180'], "the MIP computation has changed"


