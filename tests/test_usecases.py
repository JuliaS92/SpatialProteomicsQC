from pandas import read_csv
from pandas.testing import assert_frame_equal
from domaps import SpatialDataSet, SpatialDataSetComparison
import os
import pkg_resources
import json


class TestAnalysis:
    filenames = {
        "HCC827": "HCC827-gef.csv",
        "HL_U2OS": "hyperLOPIT_U2OS_short.csv",
        "DOM_Kozik": "Kozik_Dom_Data_short.csv",
        "DOM_DDA": "DDA_Borner_short.txt",
    }

    def run_analysis(_, filename):
        # set up filenames
        basename = filename.split(".")[0]
        datafile = os.path.join(os.path.dirname(__file__), "data/" + filename)
        settings = os.path.join(
            os.path.dirname(__file__), "data/" + basename + "_settings.json"
        )
        overview = os.path.join(
            os.path.dirname(__file__), "data/" + basename + "_overview.csv"
        )
        dfoutput = os.path.join(
            os.path.dirname(__file__), "data/" + basename + "_dfoutput.csv"
        )

        # set up gui
        with open(settings, "r") as file:
            settings = json.load(file)
            settings["filename"] = datafile

        # run analysis
        data = SpatialDataSet.from_settings(settings)
        data.run_pipeline()

        assert_frame_equal(
            data.reframe_df_01ORlog_for_Perseus(data.df_01_stacked).reset_index(),
            read_csv(dfoutput),
            check_dtype=False,
        )

        assert_frame_equal(
            data.results_overview_table().reset_index(),
            read_csv(overview),
            check_dtype=False,
        )

    def test_HCC827(self):
        self.run_analysis(self.filenames["HCC827"])

    def test_HL_U2OS(self):
        self.run_analysis(self.filenames["HL_U2OS"])

    def test_DOM_Kozik(self):
        self.run_analysis(self.filenames["DOM_Kozik"])

    def test_DOM_DDA(self):
        self.run_analysis(self.filenames["DOM_DDA"])


class TestBenchmark:
    def run_benchmark(_, filenames, svm=False):
        # read files
        data = dict()
        for file in filenames:
            with open(file) as file_input:
                content = json.load(file_input)
            data.update(content)

        # Set up class
        multi_choice = [k for k in data.keys()]

        comp = SpatialDataSetComparison()

        # align data
        comp.json_dict = {k: data[k] for k in multi_choice}
        comp.read_jsonFile()

        # analyse data
        comp.calc_biological_precision()
        comp.get_complex_coverage()
        comp.perform_pca_comparison()
        comp.calculate_global_scatter(
            metric="manhattan distance to average profile", consolidation="average"
        )
        if svm:
            comp.svm_processing()

        # run plotting functions
        comp.quantity_pr_pg_barplot_comparison(multi_choice)
        comp.coverage_comparison(multi_choice)
        comp.venn_sections(multi_choice, omit_size=50)
        comp.get_complex_coverage(5)
        comp.plot_intramap_scatter(
            multi_choice=multi_choice,
            clusters_for_ranking=list(comp.coverage_lists[0].keys()),
        )
        comp.plot_reproducibility_distribution(x_cut=0.2, show_full=False, q=0.7)
        if svm:
            comp.plot_svm_detail(
                multi_choice=multi_choice, orientation="h", score="F1 score"
            )

    def test_referencedata_HeLa(self):
        filenames = [
            os.path.join(os.path.dirname(__file__), "../domaps/referencedata/" + el)
            for el in pkg_resources.resource_listdir("domaps", "referencedata")
            if "HeLa" in el
        ]
        self.run_benchmark(filenames, svm=False)
