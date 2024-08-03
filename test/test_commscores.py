from cobra.io import read_sbml_model
from shutil import rmtree
from math import isclose
from glob import glob
from json import load
import commscores
import re, os


# load the test models and media
test_models, test_genomes = [], []
for model_path in glob("*.xml"):
    model = read_sbml_model(model_path)
    model.id = model_path.replace(".xml", '')
    test_models.append(model)

with open("carbonDGlucose.json", 'r') as jsonIn:
    media = load(jsonIn)



# works
# def test_mro():
#     # load the existing output
#     with open("mro_results.json", 'r') as jsonIn:
#         mro_results = load(jsonIn)

#     # validate the pure scores
#     mro_output = commscores.mro(test_models, environment=media)
#     # print(mro_output)
#     for key, val in mro_output.items():
#         if "met" in key:
#             assert set(val) == set(mro_results[key]), f"The {key} mets are different"
#             continue
#         # print(key, val, mro_results[key])
#         for index, ele in enumerate(val):
#             print(mro_results[key][index])
#             assert isclose(ele, mro_results[key][index]), f"The {key} results have changed: old {mro_results[key][index]} new {ele}"

#     # validate the metabolite predictions and compounds under competition
#     with open("mro_results2.json", 'r') as jsonIn:
#         mro_results2 = load(jsonIn)
#     mro_output = commscores.mro(test_models, environment=media, raw_content=True)
#     print(mro_output)
#     for key, val in mro_output.items():
#         assert set(val[0]) == set(mro_results2[key][0]), f"The {key} metabolites lists are different: old {mro_results2[key][0]} new {val[0]}"
#         assert val[1] == mro_results2[key][1], f"The {key} minimal media are different"


# def test_mip():
#     # affirm that the outputs are sensible
#     with open("mip_results.json", 'r') as jsonIn:
#         mip_results = load(jsonIn)
#     mip_output = commscores.mip(test_models, environment=media, costless=True)
#     for key, val in mip_output.items():
#         if key == "costless":
#             for key2, val2 in mip_output["costless"].items():
#                 # assertion bounds are necessary because there is inherent variability in the flux profile
#                 assert len(val2) - len(mip_results["costless"][key2]) <= 3, f"The costless {key2} is incorrect"
#                 assert len(val2) <= len(mip_output[key2]), f"The costless metabolites for {key2} is greater than the total metabolites"
#             continue
#         assert len(val) - len(mip_results[key]) <= 5, f"The {key} is incorrect"
    


# def test_grd():
#     # affirm that the outputs are consistent
#     with open("gyd_results.json", 'r') as jsonIn:
#         gyd_results = load(jsonIn)
#     gyd_output = commscores.gyd(test_models)
#     print(gyd_output)
#     for memberIDs, val in gyd_output.items():
#         assert set(val) == set(gyd_results[memberIDs]), f"The GYD {memberIDs} results have changed"



# def test_pc():
#     # affirm that the outputs are consistent
#     with open("pc_results.json", 'r') as jsonIn:
#         pc_results = load(jsonIn)
#     pc_output = commscores.pc(test_models)
#     print(pc_output)
#     assert isclose(pc_results[0], pc_output[0], abs_tol=1e-5), f"The PC score of {pc_output[0]} differs from the accepted value {pc_results[0]}"
#     assert pc_results[3] == pc_output[3], f"The BIT of {pc_output[0]} differs from the accepted value {pc_results[0]}"
#     for mem, CommGrowth in pc_output[1].items():
#         assert isclose(CommGrowth, pc_results[1][mem], abs_tol=1e-5), f"The growth of {mem} in the community {CommGrowth} differs from the accepted value {pc_results[1][mem]}"
#     for mem, IndGrowth in pc_output[2].items():
#         assert isclose(IndGrowth, pc_results[2][mem], abs_tol=1e-2), f"The growth of {mem} in the community {IndGrowth} differs from the accepted value {pc_results[2][mem]}"


# def test_cip():
#     # affirm that the outputs are consistent
#     with open("cip_results.json", 'r') as jsonIn:
#         cip_results = load(jsonIn)
#     cip_output = commscores.cip(member_models=test_models, environment=media)
#     assert set(cip_output[0]) == set(cip_results[0]), f"The CIP metabolites {cip_output[0]} is inconsistent with the accepted metabolites {cip_results[0]}"
#     assert cip_output[1] == cip_results[1], f"The CIP count {cip_output[1]} is inconsistent with the accepted {cip_results[1]}"
        



# Developing
# def test_bss():
#     # affirm that the outputs are consistent
#     with open("cip_results.json", 'r') as jsonIn:
#         cip_results = load(jsonIn)
#     mip_output = commscores.mip(test_models)
#     for modelID, val in mip_output.items():



# def test_fs():
#     # affirm that the outputs are consistent
#     with open("gyd_results.json", 'r') as jsonIn:
#         gyd_results = load(jsonIn)
#     fs_output = commscores.fs(test_models)
#     for modelID, val in mip_output.items():
#         assert isinstance(val[0], list), "the MIP syntrophic compounds are faulty"
