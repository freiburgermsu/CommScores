#!/usr/bin/env python3

import pathlib

import pandas as pd


def main(
    media_w_vit_folder: pathlib.Path,
    media_wo_vit_folder: pathlib.Path,
    vit_file: pathlib.Path,
):
    vitamins = list(pd.read_csv(vit_file, index_col=0)["cpd_ids"])
    for media_file in media_w_vit_folder.glob("*.tsv"):
        print(f"Processing {media_file.name}")
        media_w_vit = pd.read_csv(media_file, index_col=0, sep="\t")
        # remove vitamins from media based on index (cpd_ids)
        media_wo_vit = media_w_vit.drop(vitamins, errors="ignore")
        media_wo_vit_file = media_wo_vit_folder / media_file.name
        media_wo_vit.to_csv(media_wo_vit_file, sep="\t", index=True)


if __name__ == "__main__":
    # Media folders
    media_w_vit_folder = pathlib.Path("../../data/processed/at_leaf/media_w_vit")
    media_wo_vit_folder = pathlib.Path("../../data/processed/at_leaf/media_wo_vit")
    # Vitamin file
    vit_file = pathlib.Path("../../data/raw/at_leaf/media/vitamins.csv")
    main(media_w_vit_folder, media_wo_vit_folder, vit_file)