from pandas import read_csv
from pandas.testing import assert_frame_equal
from domaps import SpatialDataSet
import domaps.gui as gui
import os

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
        datafile = os.path.join(os.path.dirname(__file__), "data/"+filename)
        settings = os.path.join(os.path.dirname(__file__), "data/"+basename+"_settings.json")
        overview = os.path.join(os.path.dirname(__file__), "data/"+basename+"_overview.csv")
        dfoutput = os.path.join(os.path.dirname(__file__), "data/"+basename+"_dfoutput.csv")
        
        # set up gui
        interface = gui.ConfigureSingleFile()
        with open(datafile, "br") as file:
            interface._content.file.value = file.read()
            interface._content.file.filename = file.name
        with open(settings, "br") as file:
            interface._btn_load.value = file.read()
        
        # run analysis
        data = SpatialDataSet.from_settings(interface.get_settings(), legacy=False)
        data.run_pipeline()
        
        assert_frame_equal(
            data.reframe_df_01ORlog_for_Perseus(data.df_01_stacked).reset_index(),
            read_csv(dfoutput),
            check_dtype=False
        )
        
        assert_frame_equal(
            data.results_overview_table().reset_index(),
            read_csv(overview),
            check_dtype=False
        )
    
    def test_HCC827(self):
        self.run_analysis(self.filenames["HCC827"])
    
    def test_HL_U2OS(self):
        self.run_analysis(self.filenames["HL_U2OS"])
    
    def test_DOM_Kozik(self):
        self.run_analysis(self.filenames["DOM_Kozik"])
    
    def test_DOM_DDA(self):
        self.run_analysis(self.filenames["DOM_DDA"])