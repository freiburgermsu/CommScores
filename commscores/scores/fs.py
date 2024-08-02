from itertools import combinations
from typing import Iterable


# allows to singular execution of this script, besides loading CommScores as an entire package
import sys
from pathlib import Path
if __name__ == "__main__" and (__package__ is None or __package__ == ''):
    parent_dir = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(parent_dir))
    from utils import _calculate_jaccard_score
else:   from ..utils import _calculate_jaccard_score


def fs(
    models: Iterable = None,
    kbase_object=None,
    token_string: str = None,
    kbase_token_path: str = None,
    annotated_genomes: dict = None,
    printing=False,
):
    if not isinstance(annotated_genomes, dict):
        if not kbase_object:
            import cobrakbase  # ; os.environ["HOME"] = cobrakbase_repo_path ; import cobrakbase

            if token_string is not None:
                kbase_object = cobrakbase.KBaseAPI(token_string)
            else:
                with open(kbase_token_path) as token_file:
                    kbase_object = cobrakbase.KBaseAPI(token_file.readline())
        annotated_genomes = {
            model.id: kbase_object.get_from_ws(model.genome_ref)
            for model in models
            if hasattr(model, "genome_ref")
        }
    elif isinstance(annotated_genomes, list):
        annotated_genomes = dict(zip([model.id for model in models], annotated_genomes))
    elif models is not None:
        annotated_genomes = {
            k: v
            for k, v in annotated_genomes.items()
            if k in [model.id for model in models]
        }
    genome_combinations = list(combinations(annotated_genomes.keys(), 2))
    if printing:
        print(
            f"The Functionality Score (FS) will be calculated for {len(genome_combinations)} pairs."
        )
    if not isinstance(list(annotated_genomes.values())[0], dict):
        genome1_set, genome2_set = set(), set()
        distances = {}
        for genome1, genome2 in genome_combinations:
            for j in annotated_genomes[genome1].features:
                for key, val in j.ontology_terms.items():
                    if key == "SSO":
                        genome1_set.update(val)
            for j in annotated_genomes[genome2].features:
                for key, val in j.ontology_terms.items():
                    if key == "SSO":
                        genome2_set.update(val)
            distances[f"{genome1} ++ {genome2}"] = _calculate_jaccard_score(
                genome1_set, genome2_set
            )
    else:
        distances = {
            f"{genome1} ++ {genome2}": _calculate_jaccard_score(
                set(
                    list(content["SSO"].keys())[0]
                    for dic in annotated_genomes[genome1]["cdss"]
                    for x, content in dic.items()
                    if x == "ontology_terms" and len(content["SSO"].keys()) > 0
                ),
                set(
                    list(content["SSO"].keys())[0]
                    for dic in annotated_genomes[genome2]["cdss"]
                    for x, content in dic.items()
                    if x == "ontology_terms" and len(content["SSO"].keys()) > 0
                ),
            )
            for genome1, genome2 in combinations(annotated_genomes.keys(), 2)
        }
    return distances
