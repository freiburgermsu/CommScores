from itertools import chain, combinations, permutations

from numpy import array, sort, unique, where
from numpy.random import shuffle

from .logger import logger
from .scores.calculate_scores import calculate_scores
from .utils import _get_media


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
    skip_bad_media=False,
    check_models=True,
    print_progress=False,
):
    from pandas import concat

    if pairs:
        model_pairs = unique(
            [{model1, model2} for model1, models in pairs.items() for model2 in models]
        )
    elif all_models is not None:
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
                model_pairs.append(
                    [
                        list(zip(combin, smaller_list))
                        for combin in permutations(larger_list, len(smaller_list))
                    ]
                )
            # flatten the assembled pairs and filter duplicates
            model_pairs = array(
                [
                    x
                    for x in set(
                        tuple(x)
                        for x in [
                            i for y in list(chain.from_iterable(model_pairs)) for i in y
                        ]
                    )
                ]
            )
            all_models = list(chain.from_iterable(all_models))
        if pair_limit is not None:
            shuffle(model_pairs)
            new_pairs = []
            for index, pair in enumerate(model_pairs):
                if set(pair) not in exclude_pairs and index < pair_limit:
                    new_pairs.append(pair)
                elif index >= pair_limit:
                    break
            model_pairs = array(new_pairs)
        if isinstance(model_pairs[0], str):
            model_pairs = unique(sort(model_pairs, axis=1))
        pairs = {
            first: model_pairs[where(model_pairs[:, 0] == first)][:, 1]
            for first in model_pairs[:, 0]
        }
    else:
        raise ValueError(
            "Either < all_models > or < pairs > must be defined to simulate interactions."
        )
    if not all_models:
        all_models = list(chain(*[list(values) for values in pairs.values()])) + list(
            pairs.keys()
        )
    lazy_load = len(model_pairs) > 10000  # all_models[0], (list,set,tuple))
    if lazy_load and not kbase_obj:
        ValueError("The < kbase_obj > argument must be provided to lazy load models.")
    new_models = []
    for index, model in enumerate(all_models):
        if model.id == "":
            model.id = f"model_index{index}"
        new_models.append(model)
    all_models = new_models[:]
    if not mem_media:
        models_media = _get_media(model_s_=all_models, skip_bad_media=skip_bad_media)
    else:
        models_media = mem_media.copy()
        missing_models = set()
        missing_modelID = []
        for model in all_models:
            if model is not None and model.id not in models_media:
                missing_models.add(model)
                missing_modelID.append(model if not hasattr(model, "id") else model.id)
        if missing_models != set():
            logger.error(
                f"Media of the {missing_modelID} models are not defined, and will be calculated separately."
            )
            models_media.update(
                _get_media(model_s_=missing_models), skip_bad_media=skip_bad_media
            )
    if see_media:
        print(f"The minimal media of all members:\n{models_media}")
    print(f"\nExamining the {len(list(model_pairs))} model pairs")
    if pool_size is not None:
        from datetime import datetime

        from multiprocess import Pool

        print(
            f"Loading {int(pool_size)} workers and computing the scores",
            datetime.now(),
        )
        pool = Pool(
            int(pool_size)
        )  # .map(calculate_scores, [{k: v} for k,v in pairs.items()])
        args = [
            [
                dict([pair]),
                models_media,
                environments,
                annotated_genomes,
                lazy_load,
                kbase_obj,
            ]
            for pair in list(pairs.items())
        ]
        output = pool.map(calculate_scores, args)
        series = chain.from_iterable([ele[0] for ele in output])
        mets = chain.from_iterable([ele[1] for ele in output])
    else:
        series, mets = calculate_scores(
            pairs,
            models_media,
            environments,
            annotated_genomes,
            lazy_load,
            kbase_obj,
            cip_score,
            costless,
            skip_bad_media,
            check_models,
            print_progress,
        )
    return concat(series, axis=1).T, mets
