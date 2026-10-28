import numpy
from pathlib import Path
from unittest import TestCase

import pandas

from .utils import parse_sample_df

module_dir = Path(__file__).parent


class TestUtils(TestCase):
    def test_parse_sample_df(self):
        config_dir = module_dir / ".." / "configs"
        sample_path = config_dir / "sample_metadata.csv"

        sample_df = parse_sample_df(sample_path)

        simple_rows = sample_df.loc[sample_df.well_type != "Multiplexed"].index
        self.assertTrue(
            numpy.all(pandas.isnull(sample_df.loc[simple_rows, "mult_genotype"]))
        )

        multiplexed_rows = sample_df.loc[sample_df.well_type == "Multiplexed"].index
        self.assertTrue(
            numpy.all(
                sample_df.loc[multiplexed_rows, "mult_genotype"].apply(
                    lambda x: "_" in x
                )
            )
        )
