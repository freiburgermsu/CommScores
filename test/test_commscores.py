from cobra.io import read_sbml_model
from itertools import combinations
from shutil import rmtree
from math import isclose
from glob import glob
from json import load
import commscores
import re, os


relDir = os.path.dirname(__file__)

# load the test models and media
test_models, test_genomes = [], []
for model_path in glob(f"{relDir}/*.xml"):
    model = read_sbml_model(model_path)
    model.id = os.path.basename(model_path).replace(".xml", '')
    print(model.slim_optimize())
    test_models.append(model)
    print([mdl.id for mdl in test_models])

with open(f"{relDir}/carbonDGlucose.json", 'r') as jsonIn:
    media = load(jsonIn)



# works
def test_mro():
    # load the existing output
    with open(f"{relDir}/mro_results.json", 'r') as jsonIn:
        mro_results = load(jsonIn)

    # validate the pure scores
    mro_output = commscores.mro(test_models, environment=media)
    # print(mro_output)
    for key, val in mro_output.items():
        if "met" in key:
            assert set(val) == set(mro_results[key]), f"The {key} mets are different"
            continue
        # print(key, val, mro_results[key])
        for index, ele in enumerate(val):
            assert isclose(ele, mro_results[key][index]), f"The {key} results have changed: old {mro_results[key][index]} new {ele}"

    # validate the metabolite predictions and compounds under competition
    with open(f"{relDir}/mro_results2.json", 'r') as jsonIn:
        mro_results2 = load(jsonIn)
    mro_output = commscores.mro(test_models, environment=media, raw_content=True)
    print(mro_output)
    for key, val in mro_output.items():
        assert set(val[0]) == set(mro_results2[key][0]), f"The {key} metabolites lists are different: old {mro_results2[key][0]} new {val[0]}"
        assert val[1] == mro_results2[key][1], f"The {key} minimal media are different"



def test_mip():
    # affirm that the outputs are sensible
    with open(f"{relDir}/mip_results.json", 'r') as jsonIn:
        mip_results = load(jsonIn)
    mip_output = commscores.mip(test_models, environment=media, costless=True)
    for key, val in mip_output.items():
        if key == "costless":
            for key2, val2 in mip_output["costless"].items():
                # assertion bounds are necessary because there is inherent variability in the flux profile
                assert abs(len(val2) - len(mip_results["costless"][key2])) <= 3, f"The costless {key2} is incorrect"
                assert len(val2) <= len(mip_output[key2]), f"The costless metabolites for {key2} is greater than the total metabolites"
            continue
        assert abs(len(val) - len(mip_results[key])) <= 5 and val != [], f"The {key} is incorrect"
    


def test_grd():
    # affirm that the outputs are consistent
    with open(f"{relDir}/gyd_results.json", 'r') as jsonIn:
        gyd_results = load(jsonIn)
    gyd_output = commscores.gyd(test_models)
    print(gyd_output)
    for memberIDs, val in gyd_output.items():
        assert set(val) == set(gyd_results[memberIDs]), f"The GYD {memberIDs} results have changed"



def test_pc():
    # affirm that the outputs are consistent
    with open(f"{relDir}/pc_results.json", 'r') as jsonIn:
        pc_results = load(jsonIn)
    pc_output = commscores.pc(test_models)
    print(pc_output)
    assert isclose(pc_results[0], pc_output[0], abs_tol=1e-5), f"The PC score of {pc_output[0]} differs from the accepted value {pc_results[0]}"
    assert pc_results[3] == pc_output[3], f"The BIT of {pc_output[0]} differs from the accepted value {pc_results[0]}"
    for mem, CommGrowth in pc_output[1].items():
        assert isclose(CommGrowth, pc_results[1][mem], abs_tol=1e-5), f"The growth of {mem} in the community {CommGrowth} differs from the accepted value {pc_results[1][mem]}"
    for mem, IndGrowth in pc_output[2].items():
        assert isclose(IndGrowth, pc_results[2][mem], abs_tol=1e-2), f"The growth of {mem} in the community {IndGrowth} differs from the accepted value {pc_results[2][mem]}"



def test_cip():
    # affirm that the outputs are consistent
    with open(f"{relDir}/cip_results.json", 'r') as jsonIn:
        cip_results = load(jsonIn)
    cip_output = commscores.cip(member_models=test_models, environment=media)
    assert set(cip_output[0]) == set(cip_results[0]), f"The CIP metabolites {cip_output[0]} is inconsistent with the accepted metabolites {cip_results[0]}"
    assert cip_output[1] == cip_results[1], f"The CIP count {cip_output[1]} is inconsistent with the accepted {cip_results[1]}"
        


def test_bss():
    # affirm that the outputs are consistent
    with open(f"{relDir}/bss_results2.json", 'r') as jsonIn:
        bss_results = load(jsonIn)
    bss_output = commscores.bss(test_models, environment=media)
    for modelID, val in bss_output.items():
        assert set(bss_results[modelID][0]) == set(val[0]), f"The {len(val)} parasitized metabolites of {modelID} is inconsistent with the accepted value of {len(bss_results[modelID])}"
        assert bss_results[modelID][1] == val[1], f"The {val[1]} parasitized metabolites of {modelID} is inconsistent with the accepted value of {bss_results[modelID][1]}"

        
def test_fs():
    # affirm that the outputs are consistent
    from glob import glob
    from json import load
    from os import path
    
    with open(f"{relDir}/fs_results.json", 'r') as jsonIn:
        fs_results = load(jsonIn)
        
    genomes = {}
    for genetics in glob(f"{relDir}/*.fna.RAST.json"):
        baseName = path.basename(genetics).replace(".json", '.mdl')
        genomes[baseName] = load(open(genetics, 'r'))
    
    fs_output = commscores.fs(test_models, annotated_genomes=genomes)
    for combo, vals in fs_output.items():
        if combo not in fs_results:
            orgs = combo.split(" ++ ")
            combo = f"{orgs[1]} ++ {orgs[0]}"
        SSOs = set(fs_results[combo][0])
        fs = float(fs_results[combo][1])
        assert vals[0] == SSOs, f"The {combo} computed set is not identical to the established SSOs: {vals[0]} ; {SSOs}"
        assert isclose(vals[1], fs, abs_tol=1e-5), f"The {combo} computed FS {vals[1]} is not identical to the established FS {fs}"

        
        

# Developing
# def test_calculate_scores():
#     # affirm that the outputs are consistent
#     from glob import glob
#     from json import load
#     from os import path
    
#     with open(f"{relDir}/calculate_scores_results.json", 'r') as jsonIn:
#         fs_results = load(jsonIn)
        
#     # load the genomes
#     genomes = {}
#     for genetics in glob(f"{relDir}/*.fna.RAST.json"):
#         baseName = path.basename(genetics).replace(".json", '.mdl')
#         genomes[baseName] = load(open(genetics, 'r'))
    
#     # calculate the scores
#     scores_output = commscores.calculate_scores(pairs=list(combinations(all_models, 2)), member_media=test_models,
#                                             environments=media, annotated_genomes=genomes)
#     for combo, vals in fs_output.items():
#         SSOs = set(fs_results[combo][0])
#         fs = float(fs_results[combo][1])
#         assert vals[0] == SSOs, f"The {combo} computed set is not identical to the established SSOs: {vals[0]} ; {SSOs}"
#         assert isclose(vals[1], fs, abs_tol=1e-5), f"The {combo} computed FS {vals[1]} is not identical to the established FS {fs}"