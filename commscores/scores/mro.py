from modelseedpy.core.exceptions import ParameterError
from itertools import combinations
from typing import Iterable

# allows to singular execution of this script, besides loading CommScores as an entire package
import sys
from pathlib import Path
# if __name__ == "__main__" and (__package__ is None or __package__ == ''):
# Uses directory of script as starting point
parent_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(parent_dir))
from commscoresutil import CommScoresUtil
# else:   from ..commscoresutil import CommScoresUtil

def mro(member_models: Iterable = None, mem_media: dict = None, min_growth=0.1, media_dict=None,
        raw_content=False, environment=None, printing=False, compatibilized=False):
    """Determine the overlap of nutritional requirements (minimal media) between member organisms."""
    # determine the member minimal media if they are not parameterized
    if not mem_media:
        if not member_models:
            raise ParameterError("The either member_models or minimal_media parameter must be defined.")
        member_models = member_models if compatibilized else CommScoresUtil._compatibilize(member_models, printing)
        mem_media = CommScoresUtil._get_media(media_dict, None, member_models, min_growth, environment, True, printing)
        # print(mem_media)
        # print(CommScoresUtil._get_media(model_s_=member_models))
    # MROs = array(list(map(len, pairs.values()))) / array(list(map(len, mem_media.values())))
    mro_values = {}
    for model1, model2 in combinations(member_models, 2):
        intersection = set(mem_media[model1.id]["media"].keys()) & set(mem_media[model2.id]["media"].keys())
        inter = [ex.replace("EX_", "").replace("_e0", "") for ex in intersection]
        m1_media = mem_media[model1.id]["media"]
        m2_media = mem_media[model2.id]["media"]
        if raw_content:
            mro_values.update({f"{model1.id}---{model2.id}": (inter, m1_media), f"{model2.id}---{model1.id}": (inter, m2_media)})
        else:
            mro_values.update({f"{model1.id}---{model2.id}": (100 * len(inter) / len(m1_media), len(inter), len(m1_media)),
                               f"{model2.id}---{model1.id}": (100 * len(inter) / len(m2_media), len(inter), len(m2_media)),
                               f"{model2.id}---{model1.id}--mets": inter})
    return mro_values
    # return mean(list(map(len, pairs.values()))) / mean(list(map(len, mem_media.values())))
