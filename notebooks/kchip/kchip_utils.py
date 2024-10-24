""" Utilities functions for kchip data analysis """

from pathlib import Path

import pandas as pd


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
    )
    yield_df[yield_df < 0] = 0
    return yield_df
