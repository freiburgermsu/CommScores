"""
Utilities and functions for the exploration of CommScores metrics using the AT leaf microbiome
"""

import pandas as pd
import numpy as np


def read_monoculture_data(df: pd.DataFrame, media: set[str]) -> pd.DataFrame:
    monoculture_experiment_data = []
    for _, row in df.iterrows():
        # Skip if co-culture
        if "/" in row.Treatment_name:
            continue
        # Skip if unsupported media condition
        if row.Medium not in media:
            continue
        microbe = row.Treatment_name
        colonization = row["Colonization (CFU g-1 plant weight)"]
        data_item = {
            "microbe": microbe,
            "colonization": colonization,
            "dilution": row.Dilution,
        }
        monoculture_experiment_data.append(data_item)
    monoculture_experiment_all = pd.DataFrame(monoculture_experiment_data)
    return monoculture_experiment_all


def make_pair(m1: str, m2: str, measured_m: str) -> str:
    if measured_m == m1:
        p1 = m2
        p2 = m1
    elif measured_m == m2:
        p1 = m1
        p2 = m2
    else:
        raise ValueError(f"Unknown microbe {measured_m}")
    pair = "->".join([p1, p2])
    return pair


def make_pairindex(df: pd.DataFrame, i1: str, i2: str, i: str) -> pd.Index:
    inds = []
    for p1, p2, p in zip(df[i1], df[i2], df[i]):
        ind = make_pair(p1, p2, p)
        inds.append(ind)
    return pd.Index(inds)


def read_paircoculture_data(df: pd.DataFrame, media: set[str]) -> pd.DataFrame:
    pairwise_experiment_data = []
    for _, row in df.iterrows():
        name = row.Treatment_name
        # Skip if not co-culture
        if "/" not in name:
            continue
        # Skip is co-culture with more than 2 strains
        if len(name.split("/")) > 2:
            continue
        # Skip if unsupported media condition
        if row.Medium not in media:
            continue
        microbe1, microbe2 = sorted(name.split("/"))
        strain = "L" + row.Strain.strip("Leaf")
        pair = make_pair(microbe1, microbe2, strain)
        colonization = row["Colonization (CFU g-1 plant weight)"]
        data_item = {
            "microbe1": microbe1,
            "microbe2": microbe2,
            "strain": strain,
            "pair": pair,
            "dilution": row.Dilution,
            "colonization": colonization,
        }
        pairwise_experiment_data.append(data_item)
    pairwise_experiment_all = pd.DataFrame(pairwise_experiment_data)
    return pairwise_experiment_all


def get_colonization_values(
    row: pd.Series, monoculture_exp: pd.DataFrame
) -> tuple[float, float]:
    m = row.strain
    y_pair = row.colonization
    y_m = monoculture_exp.colonization[monoculture_exp.microbe == m].values[0]
    return y_pair, y_m


def get_logfc(row: pd.Series, monoculture_exp: pd.DataFrame) -> float:
    y_pair, y_m = get_colonization_values(row, monoculture_exp)
    return np.log2(y_pair / y_m)


def get_col_perc(row: pd.Series, monoculture_exp: pd.DataFrame) -> float:
    y_pair, y_m = get_colonization_values(row, monoculture_exp)
    return (y_pair / y_m) * 100


def identify_interaction(
    row: pd.Series, monoculture_exp: pd.DataFrame, thres: float
) -> str:
    val = get_logfc(row, monoculture_exp)
    weak_lb, weak_ub = -1, 1
    thres_lb, thres_ub = -thres / 2, thres / 2
    # weak_lb, weak_ub = -0.3, 0.3  # log 10
    # thres_lb, thres_ub = -0.1 / 2, 0.1 / 2  # log 10
    threshold = 0
    if thres_lb <= val <= thres_ub:
        interaction = "neutral"
    elif weak_lb <= val <= threshold:
        interaction = "weak_negative"
    elif threshold <= val <= weak_ub:
        interaction = "weak_positive"
    elif val >= threshold:
        interaction = "positive"
    else:
        interaction = "negative"
    return interaction


def get_leafid(full_gcf_id: str, gcf_leafid_map: dict[str, str]) -> str:
    gcf_id = full_gcf_id.split(".")[0] + ".1"
    # if ".1" in full_gcf_id:
    #     gcf_id = full_gcf_id.split(".")[0] + ".1"
    # else:
    #     gcf_id = "GCF_" + full_gcf_id.split("_")[1] + ".1"
    leaf_id = gcf_leafid_map.get(gcf_id, "missing")
    return leaf_id


def get_mro(mro_string: str) -> float:
    mro = float(mro_string.split(" ")[0].strip("%"))
    return mro


def get_mip(mip_string: str) -> int:
    mip = int(mip_string.split(" ")[0])
    return mip


def get_mip_costless(mip_string: str) -> int:
    mip_costless = int(mip_string.split(" ")[-1].strip("()"))
    return mip_costless


def parse_commscores_data(
    raw_commscores: pd.DataFrame, gcf_leafid_map: dict[str, str], dtype="undirected"
) -> pd.DataFrame:
    commscores_data = []
    for _, row in raw_commscores.iterrows():
        model1 = get_leafid(row.model1, gcf_leafid_map)
        model2 = get_leafid(row.model2, gcf_leafid_map)
        gr1 = row["model1 growth"]
        gr2 = row["model2 growth"]
        gr_comm = row["community growth"]
        if np.isclose(0.0, gr1) or np.isclose(0.0, gr2) or np.isclose(0.0, gr_comm):
            continue
        mip1, mip2 = map(
            get_mip, [row["MIP_model1 (costless)"], row["MIP_model2 (costless)"]]
        )
        mro1, mro2 = map(get_mro, [row.MRO_model1, row.MRO_model2])
        mip_c1, mip_c2 = map(
            get_mip_costless,
            [row["MIP_model1 (costless)"], row["MIP_model2 (costless)"]],
        )
        if dtype == "undirected":
            data_item = {
                "model1": model1,
                "model2": model2,
                "gr_comm": gr_comm,
                "mro": (mro1 + mro2) / 2,
                "cip": row.CIP,
                "mip": (mip1 + mip2) // 2,
                "mip_c": (mip_c1 + mip_c2) // 2,
                "bss": (row.BSS_model1 + row.BSS_model2) / 2,
                "pc_comm": row.PC_comm,
                "bit": row.BIT,
                "gyd": 1 if row.GYD > 1 else row.GYD,
                "fs": row.FS,
            }
            commscores_data.append(data_item)
        elif dtype == "directed":
            data_item = {
                "model1": model1,
                "model2": model2,
                "gr1": gr1,
                "gr2": gr2,
                "gr_comm": gr_comm,
                "gyd": (gr1 - gr2) / gr1,
                "mro1": mro1,
                "mro2": mro2,
                "cip": row.CIP,
                "mip1": mip1,
                "mip2": mip2,
                "mip_c1": mip_c1,
                "mip_c2": mip_c2,
                "bss1": row.BSS_model1,
                "bss2": row.BSS_model2,
                "pc1": row.PC_model1,
                "pc2": row.PC_model2,
                "pc_comm": row.PC_comm,
                "fs": row.FS,
            }
            commscores_data.append(data_item)
            data_item = {
                "model1": model2,
                "model2": model1,
                "gr1": gr2,
                "gr2": gr1,
                "gr_comm": gr_comm,
                "gyd": (gr2 - gr1) / gr2,
                "mro1": mro2,
                "mro2": mro1,
                "cip": row.CIP,
                "mip1": mip2,
                "mip2": mip1,
                "mip_c1": mip_c2,
                "mip_c2": mip_c1,
                "bss1": row.BSS_model2,
                "bss2": row.BSS_model1,
                "pc1": row.PC_model2,
                "pc2": row.PC_model1,
                "pc_comm": row.PC_comm,
                "fs": row.FS,
            }
            commscores_data.append(data_item)
        else:
            raise ValueError(f"Unsupported dtype {dtype}")

    commscores = pd.DataFrame(commscores_data)
    commscores.index = make_pairindex(commscores, "model1", "model2", "model1")
    return commscores
