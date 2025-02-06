from domaps import SpatialDataSetComparison
import json
import os
import pkg_resources
from pytest import fixture
import pandas as pd
from unittest.mock import patch

filenames_HeLa = [
    os.path.join(os.path.dirname(__file__), "../domaps/referencedata/" + el)
    for el in pkg_resources.resource_listdir("domaps", "referencedata")
    if "HeLa" in el
]

filenames_minimal = [
    os.path.join(
        os.path.dirname(__file__), "data/DDA_Borner_short_AnalysedDatasets.json"
    ),
]


def SDC(filenames: list) -> SpatialDataSetComparison:
    comp = SpatialDataSetComparison()
    data = dict()
    for file in filenames:
        with open(file) as file_input:
            content = json.load(file_input)
        data.update(content)
    multi_choice = [k for k in data.keys()]
    comp.json_dict = {k: data[k] for k in multi_choice}
    comp.read_jsonFile()
    return comp


@fixture
def SDC_HeLa() -> SpatialDataSetComparison:
    return SDC(filenames_HeLa)


@fixture
def SDC_minimal() -> SpatialDataSetComparison:
    return SDC(filenames_minimal)


valid_df_01_filtered_combined = pd.DataFrame(
    SDC(filenames_minimal).df_01_filtered_combined.copy()
)


@patch("domaps.SpatialDataSetComparison._calc_biological_precision")
def test_calc_biological_precision_attributes(
    mock_calc_biological_precision, SDC_minimal
):
    """Test that the function updates the attribute and does not add new ones."""
    mock_calc_biological_precision.return_value = pd.DataFrame(columns=["test"])
    SDC = SDC_minimal
    attr_before = SDC.__dict__.keys()
    df_returned = SDC.calc_biological_precision()
    attr_after = SDC.__dict__.keys()

    # No new attributes written
    assert attr_before == attr_after

    # Attribute is updated
    pd.testing.assert_frame_equal(SDC.df_distance_comp, df_returned)


def test_calc_biological_precision_values():
    """Test that the private function calculates distances correclty."""
    df_returned = SpatialDataSetComparison._calc_biological_precision(
        df01=valid_df_01_filtered_combined,
        experiments=["DDA_Borner"],
        markerproteins={"Dummy": ["O00116", "O00299"]},
        clusters=["Dummy"],
    )

    # Column names are correct
    assert list(df_returned.columns) == [
        "Exp_Map",
        "Gene names",
        "Protein IDs",
        "merge type",
        "Compartment",
        "Cluster",
        "Experiment",
        "Map",
        "distance",
    ]

    # Values are correct
    assert df_returned.iloc[0, :].to_list() == [
        "DDA_Borner_Map1",
        "AGPS",
        "O00116",
        "primary id",
        "Peroxisome",
        "Dummy",
        "DDA_Borner",
        "Map1",
        0.31006103920000005,
    ]

    d_1 = df_returned.loc[
        df_returned["Protein IDs"] == "O00116", "distance"
    ].reset_index(drop=True)
    d_2 = df_returned.loc[
        df_returned["Protein IDs"] == "O00299", "distance"
    ].reset_index(drop=True)

    subtracted_distance = d_1 - d_2

    pd.testing.assert_series_equal(
        subtracted_distance,
        pd.Series([0.0, 0.0, 0.0], name="distance"),
        check_exact=False,
    )
