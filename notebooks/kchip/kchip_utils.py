""" Utilities functions for kchip data analysis """

from pathlib import Path

import cobra
import pandas as pd
import seaborn as sns
from cobra.io import read_sbml_model


def read_experimental_yield(file_path: Path) -> pd.DataFrame:
    met_profile = pd.read_csv(file_path, index_col=0).T
    met_profile = met_profile.rename(index={"SF1": "SF"})
    met_profile = met_profile.rename(
        columns={"ArabinoseD": "D-Arabinose", "ArabinoseL": "L-Arabinose"}
    )
    # labels extracted from heatmap image
    row_labels = [
        "PAl",
        "PAg2",
        "PAg1",
        "PAg3",
        "LA",
        "EL",
        "RP2",
        "RP1",
        "SF",
        "BI",
        "KA",
        "CF",
        "EC",
        "EA",
        "PK",
        "PR1",
        "PAr",
        "PR2",
        "PH",
        "PP",
    ]
    # ignore the diluted media conditions
    column_labels = [
        "Glycerol",
        "Glutamine",
        "Cellobiose",
        "Rhamnose",
        "Maltose",
        "Mannose",
        "GlcNAc",
        "Trehalose",
        "Glucose",
        "Mix",
        "Pyruvate",
        "Alanine",
        "Fructose",
        "Galactose",
        "Ribose",
        "Xylose",
        "Mannitol",
        "L-Arabinose",
        "Sorbitol",
        "Lactose",
        "Sucrose",
        "Raffinose",
        "Uridine",
        "Arabinogalactan",
        "Melezitose",
        "Water",
        "D-Arabinose",
        "Serine",
        "Isoleucine",
        "Arginine",
        "Acetate",
        "Citrate",
        "Fumarate",
        "Succinate",
        "Proline",
    ]
    yield_df: pd.DataFrame = met_profile.loc[row_labels, column_labels].apply(
        lambda row: row / row.max(), axis=1
    )  # type: ignore
    yield_df[yield_df < 0] = 0
    return yield_df


def load_models_and_media(
    models_folder: Path, media_folder: Path
) -> tuple[dict[str, cobra.Model], dict[str, list[str]]]:
    model_dict = {}
    media_dict = {}
    for model_file in models_folder.glob("*.sbml"):
        label = model_file.stem.removesuffix(".genome.mdl")
        model_dict[label] = read_sbml_model(model_file)
    for media_file in media_folder.glob("*.tsv"):
        label = media_file.stem
        media_dict[label] = list(pd.read_csv(media_file, sep="\t")["compounds"])
    return model_dict, media_dict


def simulate_predicted_growth(
    model_dict: dict[str, cobra.Model],
    media_dict: dict[str, list[str]],
    carbon_sources: dict[str, list[str]],
) -> pd.DataFrame:
    predicted_growth_data = []
    for model_name, model in model_dict.items():
        for media_name, media_compounds in media_dict.items():
            curr_model = model.copy()
            # update the medium
            curr_media = curr_model.medium
            compound_exchanges = [f"EX_{m}_e0" for m in media_compounds]
            for m in curr_media:
                if m in compound_exchanges:
                    if m in carbon_sources[media_name]:
                        curr_media[m] = 10
                    else:
                        curr_media[m] = 1000
                else:
                    curr_media[m] = 0
            # Constrain the O2 exchange flux
            curr_media[f"EX_cpd00007_e0"] = 10 / 5
            curr_model.medium = curr_media
            # simulate growth
            growth = curr_model.slim_optimize()
            data_item = {"model": model_name, "media": media_name, "growth": growth}
            predicted_growth_data.append(data_item)
    predicted_growth = pd.DataFrame(predicted_growth_data).pivot_table(
        index="model", columns="media", values="growth"
    )
    return predicted_growth


def plot_confusion_matrix(pred_data: pd.DataFrame, exp_data: pd.DataFrame) -> None:
    tp = ((pred_data > 0) & (exp_data > 0)).sum(axis=1)
    fp = ((pred_data > 0) & (exp_data == 0)).sum(axis=1)
    tn = ((pred_data == 0) & (exp_data == 0)).sum(axis=1)
    fn = ((pred_data == 0) & (exp_data > 0)).sum(axis=1)
    confusion_df = pd.DataFrame({"TP": tp, "FP": fp, "TN": tn, "FN": fn})
    confusion_df.index.name = "model"
    plot_data = confusion_df.reset_index().melt(
        id_vars="model",
        value_vars=["TP", "FP", "TN", "FN"],
        var_name="Type",
        value_name="Count",
    )
    g = sns.catplot(
        data=plot_data,
        x="model",
        y="Count",
        col="Type",
        kind="bar",
        aspect=2,
        col_wrap=2,
    )
    # change ylimit for each axis
    for ax in g.axes:
        ax.set_ylim(0, pred_data.shape[1])
