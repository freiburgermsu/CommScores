from modelseedpy.core.exceptions import ParameterError
from numpy import array


def antiSMASH(json_path=None, zip_path=None):
    # TODO Scores 2, 4, and 5 are being explored for relevance to community formation and reveal specific member interactions/targets
    # load the antiSMASH report from either the JSON or the raw ZIP, or both
    from json import load
    from os import listdir, mkdir, path
    from zipfile import ZipFile

    if json_path:
        cwd_files = listdir()
        if json_path not in cwd_files and zip_path:
            with ZipFile(zip_path, "r") as zip_file:
                zip_file.extract(json_path)
        with open(json_path, "r") as json_file:
            data = load(json_file)
    elif zip_path:
        mkdir("extracted_antiSMASH")
        with ZipFile(zip_path, "r") as zip_file:
            zip_file.extractall("extracted_antiSMASH")
        json_files = [x for x in listdir("extracted_antiSMASH") if x.endswith("json")]
        if len(json_files) > 1:
            print(
                f"The antiSMASH report describes {len(json_files)} JSON files, the first of which is selected "
                f"{json_files[0]} for analysis, otherwise explicitly identify the desired JSON file in the json_path parameter."
            )
        with open(path.join("extracted_antiSMASH", json_files[0]), "r") as json_file:
            data = load(json_file)
    else:
        raise ParameterError(
            "Either the json_path or zip_path from the antiSMASH analysis must be provided,"
            " for these scores to be determined."
        )
    # Parse data and scores from the antiSMASH report
    biosynthetic_areas = data["records"][0]["areas"]
    BGCs = set(
        array(
            [
                data["records"][0]["areas"][i]["products"]
                for i in range(biosynthetic_areas)
            ]
        ).flatten()
    )
    len_proteins = len(
        data["records"][0]["modules"]["antismash.modules.clusterblast"]["knowncluster"][
            "proteins"
        ]
    )
    protein_annotations = [
        data["records"][0]["modules"]["antismash.modules.clusterblast"]["knowncluster"][
            "proteins"
        ][i]["annotations"]
        for i in range(len_proteins)
    ]
    clusterBlast = [s for s in protein_annotations if "resistance" in s]
    num_clusterBlast = sum([item.count("resistance") for item in protein_annotations])

    return (
        biosynthetic_areas,
        BGCs,
        protein_annotations,
        clusterBlast,
        num_clusterBlast,
    )
