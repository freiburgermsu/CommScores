#!/usr/bin/env python3

from pathlib import Path
from typing import Iterable

import pandas as pd

CPDS = {
    "co2": {"compounds": "cpd00011", "name": "CO2"},
    "D-Xylose": {"compounds": "cpd00154", "name": "Xylose"},
    "D-Raffinose pentahydrate": {"compounds": "cpd00382", "name": "Melitose"},
}


def main(media_files: Iterable[Path], output_folder: Path) -> None:
    output_folder.mkdir(parents=True, exist_ok=True)
    for media_file in media_files:
        media_df = pd.read_excel(media_file, sheet_name=0, index_col=0)
        # 1. Remove CO2
        media_df.drop(CPDS["co2"]["compounds"], inplace=True)
        # TODO: Manually replace carbon sources
        media_file_name = media_file.stem.removeprefix("M9_") + ".tsv"
        output_file = output_folder / media_file_name
        media_df.to_csv(output_file, sep="\t", index=True)


if __name__ == "__main__":
    media_files = Path("media_old").glob("*.xls")
    output_folder = Path("media")
    main(media_files, output_folder)
