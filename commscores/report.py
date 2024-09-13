import os
import re
from itertools import chain, combinations, permutations

from numpy import array, sort, unique, where
from numpy.random import shuffle
from pandas import concat
from logger import logger
from scores.calculate_scores import calculate_scores
from commscoresutil import CommScoresUtil

package_dir = os.path.abspath(os.path.dirname(__file__))


def create_pairs(all_models, pair_limit: int = None):
    if not isinstance(all_models[0], list):
        all_models = list(set(all_models))
        model_pairs = array(list(combinations(all_models, 2)))
    else:
        model_pairs = []
        for models1, models2 in combinations(all_models, 2):
            models1 = set(models1)
            models2 = set(models2)
            if len(models1) > len(models2):
                larger_list = models1
                smaller_list = models2
            else:
                larger_list = models2
                smaller_list = models1
            model_pairs.append([list(zip(combin, smaller_list))
                                for combin in permutations(larger_list, len(smaller_list))])
        # flatten the assembled pairs and filter duplicates
        model_pairs = array([x for x in set(tuple(x) for x in [
            i for y in list(chain.from_iterable(model_pairs)) for i in y])])
        all_models = list(chain.from_iterable(all_models))
    if pair_limit is not None:
        shuffle(model_pairs)
        new_pairs = []
        for index, pair in enumerate(model_pairs):
            if set(pair) not in exclude_pairs and index < pair_limit:
                new_pairs.append(pair)
            elif index >= pair_limit:  break
        model_pairs = array(new_pairs)
    if isinstance(model_pairs[0], str):
        model_pairs = unique(sort(model_pairs, axis=1))
    pairs = {first: model_pairs[where(model_pairs[:, 0] == first)][:, 1]
             for first in model_pairs[:, 0]}
    return model_pairs, pairs


def report_generation(
    all_models: iter = None,  # a list of distinct lists is provided for specifying exclusive groups
    pairs: dict = None,
    mem_media: dict = None,
    pair_limit: int = None,
    exclude_pairs: list = None,
    kbase_obj=None,
    annotated_genomes: dict = True,  # True triggers internal acquisition of the genomes, where None skips
    see_media=True,
    environments: iter = None,  # a collection of environment dicts or KBase media objects
    pool_size: int = None,
    cip_score=True,
    costless=True,
    check_models=True,
    ids = None,
    print_progress=False,
):
    # define the model pairs
    if pairs:
        model_pairs = unique([{model1, model2} for model1, models in pairs.items() for model2 in models])
    elif all_models is not None:
        model_pairs, pairs = create_pairs(all_models, pair_limit)
    else:
        raise ValueError("Either < all_models > or < pairs > must be defined to simulate interactions.")
    if not all_models:
        all_models = list(chain(*[list(values) for values in pairs.values()])) + list(pairs.keys())
        
    # assign IDs to the models
    new_models = []
    for index, model in enumerate(all_models):
        if model.id == "":
            if ids is not None:     model.id = ids[index]
            else:     model.id = f"model{index}"
        new_models.append(model)
    all_models = new_models[:]
    models_media = mem_media.copy()
    # define the minimal media of each model
    if mem_media is None:
        models_media = {}
        for model in all_models:
            print(model.id)
            models_media[model.id] = {}
            for name, environ in environments:
                print(name)
                models_media[model.id][name] = CommScoresUtil._get_media(model_s_=model, environment=environ)
    else:
        models_media = mem_media.copy()
        missing_models = set()
        missing_modelID = []
        for model in all_models:
            if model is not None and model.id not in models_media:
                missing_models.add(model)
                missing_modelID.append(model if not hasattr(model, "id") else model.id)
        if missing_models != set():
            logger.error(f"Media of the {missing_modelID} models are not defined, and will be calculated separately.")
            models_media.update(CommScoresUtil._get_media(model_s_=missing_models))
    if see_media:
        print(f"The minimal media of all members:")
        display(models_media)
    print(f"\nExamining the {len(list(model_pairs))} model pairs")
    
    # run the simulation, and paralellize or not
    lazy_load = len(model_pairs) > 1000  # all_models[0], (list,set,tuple))
    if lazy_load:
        assert kbase_obj is not None, "The < kbase_obj > argument must be provided to lazy load models."
    if pool_size is not None:
        from datetime import datetime
        from multiprocess import Pool

        pool_size = min([pool_size, len(pairs)*5])
        print(f"Loading {int(pool_size)} workers and computing the scores", datetime.now())
        pool = Pool(int(pool_size))  # .map(calculate_scores, [{k: v} for k,v in pairs.items()])
        args = [[dict([pair]), models_media, environments, annotated_genomes, lazy_load, kbase_obj, check_models, print_progress]
                for pair in list(pairs.items())]
        output = pool.map(calculate_scores, args)
        series = chain.from_iterable([ele[0] for ele in output])
        mets = chain.from_iterable([ele[1] for ele in output])
    else:
        series, mets = calculate_scores(pairs, models_media, environments, annotated_genomes,
                                        lazy_load, kbase_obj, costless, check_models, print_progress)
    return concat(series, axis=1).T, mets


def html_report(df, mets, export_html_path="commscores_report.html", msdb_path=None):
    import jinja2
    from pandas import to_numeric

    def names_updateCPD(metIDs, update_cpdNames):
        names = []
        for metID in metIDs:
            if metID not in cpdNames:
                if "msdb" not in locals().keys():
                    from modelseedpy.biochem import from_local

                    assert msdb_path, "An MSDB path is needed to load compound names"
                    msdb = from_local(msdb_path)
                name = msdb.compounds.get(metID, None)
                if name is None:   name = metID
                else:   name = name["name"]
                update_cpdNames[metID] = name
            else:  name = cpdNames[metID]
            names.append(name)
        return names, update_cpdNames

    # construct a heatmap
    df.index.name = "Community_index"
    heatmap_df = df.copy(deep=True)  # takes some time
    heatmap_df_index = zip(heatmap_df["model1"].to_numpy(), heatmap_df["model2"].to_numpy())
    heatmap_df.index = [" ++ ".join(index) for index in heatmap_df_index]
    heatmap_df.index.name = "model1 ++ model2"
    if "media" in heatmap_df.columns:
        media_list = heatmap_df["media"].tolist()
        new_index = [f"{models} in {media_list[i]}" for i, models in enumerate(heatmap_df.index)]
        heatmap_df.index = new_index
        heatmap_df.index.name = "model1 ++ model2 in Media"
    heatmap_df = heatmap_df.loc[~heatmap_df.index.duplicated(), :]
    heatmap_df = heatmap_df.drop(["model1", "model2"], axis=1)
    if "media" in heatmap_df:
        heatmap_df = heatmap_df.drop(["media"], axis=1)
    costless = re.compile(r"(?<=\s\()(\d)(?=\))")
    if "MIP_model1 (costless)" in heatmap_df.columns:
        mip_model1, mip_model2 = [], []
        for e in heatmap_df["MIP_model1 (costless)"]:
            if e == "":
                mip_model1.append("")
                continue
            mip_model1.append(costless.search(str(e)).group() if e not in [0, "0"] else "")
        for e in heatmap_df["MIP_model2 (costless)"]:
            if e == "":
                mip_model2.append("")
                continue
            mip_model2.append(costless.search(str(e)).group() if e not in [0, "0"] else "")
        for col, lis in {
            "c_MIP1": mip_model1,
            "c_MIP2": mip_model2,
            "MIP_model1": heatmap_df["MIP_model1 (costless)"].apply(CommScoresUtil.remove_metadata),
            "MIP_model2": heatmap_df["MIP_model2 (costless)"].apply(CommScoresUtil.remove_metadata),
        }.items():
            heatmap_df[col] = to_numeric(lis, errors="coerce")
    for col in ["MRO_model1", "MRO_model2", "BSS_model1", "BSS_model2", "PC_model1",
                "PC_model2", "FS", "GYD"]:
        if col not in heatmap_df:
            print(f"The {col} is not computed")
            continue
        heatmap_df[col] = to_numeric(heatmap_df[col].apply(CommScoresUtil.remove_metadata), errors="coerce")
     # TODO colorize the BIT entries as well
    del (heatmap_df["BIT"], heatmap_df["MIP_model1 (costless)"], heatmap_df["MIP_model2 (costless)"])
    heatmap_df = heatmap_df.astype(float)
    int_cols = ["CIP", "MIP_model1", "MIP_model2"]
    if "costless_MIP_model1" in heatmap_df.columns:  int_cols.extend(["c_MIP1", "c_MIP2"])
    for col in int_cols:
        heatmap_df[col] = heatmap_df[col].apply(CommScoresutil.convert_to_int)

    # construct a metabolites table
    from pandas import DataFrame

    # from pandas import set_option
    # set_option("display.max_colwidth", None, 'display.width', 1500)
    ## Process the score metabolites
    mro_mets, mro_mets_names, mip_model1_mets, mip_model1_mets_names = ([],[],[],[])
    mip_model2_mets, mip_model2_mets_names, cip_mets, cip_mets_names = ([],[],[],[])
    from json import dump, load

    cpdNames_path = os.path.join(package_dir, "data", "compoundNames.json")
    with open(cpdNames_path, "r") as jsonIn:  cpdNames = load(jsonIn)
    update_cpdNames = {}
    for met in mets:
        # MRO metabolites
        mro_metIDs = [metID for metID in map(str, met["MRO metabolites"]) if metID not in ["None", None]]
        mro_mets.append(", ".join(mro_metIDs))
        names, update_cpdNames = names_updateCPD(mro_metIDs, update_cpdNames)
        mro_mets_names.append(", ".join(names))
        # MIP metabolites
        mip_model1_metIDs = [metID for metID in map(str, met["MIP model1 metabolites"]) if metID not in ["None", None]]
        mip_model1_mets.append(", ".join(mip_model1_metIDs))
        names, update_cpdNames = names_updateCPD(mip_model1_metIDs, update_cpdNames)
        mip_model1_mets_names.append(", ".join(names))
        ## model2 MIP metabolites
        mip_model2_metIDs = [metID for metID in map(str, met["MIP model2 metabolites"]) if metID not in ["None", None]]
        mip_model2_mets.append(", ".join(mip_model2_metIDs))
        names, update_cpdNames = names_updateCPD(mip_model2_metIDs, update_cpdNames)
        mip_model2_mets_names.append(", ".join(names))
        # CIP metabolites
        cip_metIDs = [metID for metID in map(str, met["CIP metabolites"]) if metID not in ["None", None]]
        cip_mets.append(", ".join(cip_metIDs))
        names, update_cpdNames = names_updateCPD(cip_metIDs, update_cpdNames)
        cip_mets_names.append(", ".join(names))
    df_content = {
        "MRO metabolite names": mro_mets_names,
        "MRO metabolite IDs": mro_mets,
        "MIP model1 metabolite names": mip_model1_mets_names,
        "MIP model1 metabolite IDs": mip_model1_mets,
        "MIP model2 metabolite names": mip_model2_mets_names,
        "MIP model2 metabolite IDs": mip_model2_mets,
        "CIP metabolite names": cip_mets_names,
        "CIP metabolite IDs": cip_mets,
    }
    if update_cpdNames != {}:
        cpdNames.update(update_cpdNames)
        with open(cpdNames_path, "w") as jsonOut:
            dump(cpdNames, jsonOut, indent=3)
    # print(list(map(len, df_content.values())))
    mets_table = DataFrame(data=df_content)
    mets_table.index.name = "Community_index"

    # populate the HTML template with the assembled simulation data from the DataFrame -> HTML conversion
    content = {
        "table": df.to_html(table_id="main", classes="display"),
        "mets_table": mets_table.to_html(),
        "heatmap": heatmap_df.applymap(lambda x: round(x, 3))
        .style.background_gradient()
        .to_html(table_id="heat", classes="display"),
    }
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(os.path.join(package_dir, "community")),
        autoescape=jinja2.select_autoescape(["html", "xml"]),
    )
    html_report = env.get_template("commscores_template.html").render(content)
    with open(export_html_path, "w") as out:
        out.writelines(html_report)
    return html_report
