import natsort
import numpy as np
import pandas as pd
import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import re
import traceback
from io import BytesIO
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise as pw
import json
import statistics 
import matplotlib.pyplot as plt
import matplotlib_venn as venn
from matplotlib_venn import venn2, venn3, venn3_circles
from PIL import Image
from upsetplot import from_memberships
from upsetplot import plot as upplot
import pkg_resources
from scipy.stats import zscore
import urllib.parse
import urllib.request
import bisect

def natsort_index_keys(x):
    order = natsort.natsorted(np.unique(x.values))
    return pd.Index([order.index(el) for el in x], name=x.name)
    
def natsort_list_keys(x):
    order = natsort.natsorted(np.unique(x))
    return [order.index(el) for el in x]

class SpatialDataSet:
    
    regex = {
        "imported_columns": "^[Rr]atio H/L (?!normalized|type|is.*|variability|count)[^ ]+|^Ratio H/L variability.... .+|^Ratio H/L count .+|id$|[Mm][Ss].*[cC]ount.+$|[Ll][Ff][Qq].*|.*[nN]ames.*|.*[Pp][rR]otein.[Ii][Dd]s.*|[Pp]otential.[cC]ontaminant|[Oo]nly.[iI]dentified.[bB]y.[sS]ite|[Rr]everse|[Ss]core|[Qq]-[Vv]alue|R.Condition|PG.Genes|PG.ProteinGroups|PG.Cscore|PG.Qvalue|PG.RunEvidenceCount|PG.Quantity|^Proteins$|^Sequence$"
    }
    
    acquisition_set_dict = {
        "LFQ6 - Spectronaut" : ["LFQ intensity", "MS/MS count"],
        "LFQ5 - Spectronaut" : ["LFQ intensity", "MS/MS count"],
        "LFQ5 - MQ" : ["[Ll][Ff][Qq].[Ii]ntensity", "[Mm][Ss]/[Mm][Ss].[cC]ount", "[Ii]ntensity"],
        "LFQ6 - MQ" : ["[Ll][Ff][Qq].[Ii]ntensity", "[Mm][Ss]/[Mm][Ss].[cC]ount", "[Ii]ntensity"],
        "SILAC - MQ" : [ "[Rr]atio.[Hh]/[Ll](?!.[Vv]aria|.[Cc]ount)","[Rr]atio.[Hh]/[Ll].[Vv]ariability.\[%\]", "[Rr]atio.[Hh]/[Ll].[cC]ount"],
        "Custom": ["(?!Protein IDs|Gene names)"]
    }
    
    Spectronaut_columnRenaming = {
        "R.Condition": "Map", "PG.Genes" : "Gene names", "PG.Qvalue": "Q-value", "PG.Cscore":"C-Score", 
        "PG.ProteinGroups" : "Original Protein IDs", "PG.RunEvidenceCount" : "MS/MS count", "PG.Quantity" : "LFQ intensity"
    }
    
    css_color = ["#b2df8a", "#6a3d9a", "#e31a1c", "#b15928", "#fdbf6f", "#ff7f00", "#cab2d6", "#fb9a99", "#1f78b4", "#ffff99", "#a6cee3", 
                      "#33a02c", "blue", "orange", "goldenrod", "lightcoral", "magenta", "brown", "lightpink", "red", "turquoise",
                      "khaki", "darkgoldenrod","darkturquoise", "darkviolet", "greenyellow", "darksalmon", "hotpink", "indianred", "indigo","darkolivegreen", 
                      "coral", "aqua", "beige", "bisque", "black", "blanchedalmond", "blueviolet", "burlywood", "cadetblue", "yellowgreen", "chartreuse",
                      "chocolate", "cornflowerblue", "cornsilk", "darkblue", "darkcyan", "darkgray", "darkgrey", "darkgreen", "darkkhaki", "darkmagenta", 
                      "darkorange", "darkorchid", "darkred", "darkseagreen", "darkslateblue", "snow", "springgreen", "darkslategrey", "mediumpurple", "oldlace", 
                      "olive", "lightseagreen", "deeppink", "deepskyblue", "dimgray", "dimgrey", "dodgerblue", "firebrick", "floralwhite", "forestgreen", 
                      "fuchsia", "gainsboro", "ghostwhite", "gold", "gray", "ivory", "lavenderblush", "lawngreen", "lemonchiffon", "lightblue", "lightcyan",
                      "fuchsia", "gainsboro", "ghostwhite", "gold", "gray", "ivory", "lavenderblush", "lawngreen", "lemonchiffon", "lightblue", "lightcyan",
                      "lightgoldenrodyellow", "lightgray", "lightgrey", "lightgreen", "lightsalmon", "lightskyblue", "lightslategray", "lightslategrey",
                      "lightsteelblue", "lightyellow", "lime", "limegreen", "linen", "maroon", "mediumaquamarine", "mediumblue", "mediumseagreen",
                      "mediumslateblue", "mediumspringgreen", "mediumturquoise", "mediumvioletred", "midnightblue", "mintcream", "mistyrose", "moccasin",
                      "olivedrab", "orangered", "orchid", "palegoldenrod", "palegreen", "paleturquoise", "palevioletred", "papayawhip", "peachpuff", "peru",
                      "pink", "plum", "powderblue", "rosybrown", "royalblue", "saddlebrown", "salmon", "sandybrown", "seagreen", "seashell", "sienna", "silver",
                      "skyblue", "slateblue", "steelblue", "teal", "thistle", "tomato", "violet", "wheat", "white", "whitesmoke", "slategray", "slategrey",
                      "aquamarine", "azure","crimson", "cyan", "darkslategray", "grey","mediumorchid","navajowhite", "navy"]
    

    analysed_datasets_dict = {}
    
    df_organellarMarkerSet = pd.read_csv(pkg_resources.resource_stream(__name__, 'annotations/organellemarkers/{}.csv'.format("Homo sapiens - Uniprot")),
                                       usecols=lambda x: bool(re.match("Compartment|Protein ID", x)))

    def __init__(self, filename, expname, acquisition, comment, name_pattern="e.g.:.* (?P<cond>.*)_(?P<rep>.*)_(?P<frac>.*)", reannotate_genes=False, **kwargs):
        
        self.filename = filename
        self.expname = expname
        self.acquisition = acquisition
        self.name_pattern = name_pattern
        self.comment = comment
        self.imported_columns = self.regex["imported_columns"]
        
        self.fractions, self.map_names = [], []
        self.df_01_stacked, self.df_log_stacked = pd.DataFrame(), pd.DataFrame()
        self.reannotate_genes = reannotate_genes
        
        if reannotate_genes:
            self.reannotation_source = kwargs["reannotation_source"]
        
        if acquisition == "SILAC - MQ":
            if "RatioHLcount" not in kwargs.keys():
                self.RatioHLcount = 2
            else:
                self.RatioHLcount = kwargs["RatioHLcount"]
                del kwargs["RatioHLcount"]
            if "RatioVariability" not in kwargs.keys():
                self.RatioVariability = 30
            else:
                self.RatioVariability = kwargs["RatioVariability"]
                del kwargs["RatioVariability"]
        
        elif acquisition == "Custom":
            self.custom_columns = kwargs["custom_columns"]
            self.custom_normalized = kwargs["custom_normalized"]
            self.imported_columns = "^"+"$|^".join(["$|^".join(el) if type(el) == list else el for el in self.custom_columns.values() if el not in [[], None, ""]])+"$"
            
        #elif acquisition == "LFQ5 - MQ" or acquisition == "LFQ6 - MQ" or acquisition == "LFQ6 - Spectronaut" or acquisition == "LFQ5 - Spectronaut":
        else:
            if "summed_MSMS_counts" not in kwargs.keys():
                self.summed_MSMS_counts = 2
            else:
                self.summed_MSMS_counts = kwargs["summed_MSMS_counts"]
                del kwargs["summed_MSMS_counts"]
            if "consecutiveLFQi" not in kwargs.keys():
                self.consecutiveLFQi = 4
            else:
                self.consecutiveLFQi = kwargs["consecutiveLFQi"]
                del kwargs["consecutiveLFQi"]
        
        #self.markerset_or_cluster = False if "markerset_or_cluster" not in kwargs.keys() else kwargs["markerset_or_cluster"]
        if "organism" not in kwargs.keys():
            marker_table = pd.read_csv(pkg_resources.resource_stream(__name__, 'annotations/complexes/{}.csv'.format("Homo sapiens - Uniprot")))
            self.markerproteins = {k: v.replace(" ", "").split(",") for k,v in zip(marker_table["Cluster"], marker_table["Members - Protein IDs"])}
            df_organellarMarkerSet = pd.read_csv(
                                                 pkg_resources.resource_stream(__name__,'annotations/organellemarkers/{}.csv'.format("Homo sapiens - Uniprot")),
                                                 usecols=lambda x: bool(re.match("Compartment|Protein ID", x)))
            self.df_organellarMarkerSet = df_organellarMarkerSet
        else:
            assert kwargs["organism"]+".csv" in pkg_resources.resource_listdir(__name__, "annotations/complexes")
            marker_table = pd.read_csv(pkg_resources.resource_stream(__name__, 'annotations/complexes/{}.csv'.format(kwargs["organism"])))
            self.markerproteins = {k: v.replace(" ", "").split(",") for k,v in zip(marker_table["Cluster"], marker_table["Members - Protein IDs"])}
            df_organellarMarkerSet = pd.read_csv(
                                                 pkg_resources.resource_stream(__name__,'annotations/organellemarkers/{}.csv'.format(kwargs["organism"])),
                                                 usecols=lambda x: bool(re.match("Compartment|Protein ID", x)))
            self.df_organellarMarkerSet = df_organellarMarkerSet
            self.organism = kwargs["organism"]
            del kwargs["organism"]
        
        self.analysed_datasets_dict = {}
        self.analysis_summary_dict = {}
    
    
    def run_pipeline(self, content=None, progressbar=None):
        """
        Run the whole processing pipeline without generating any visual output.
        
        Args:
            self: needs to be fully configured by init arguments
            content: None or file like object. Passed on to data_reading.
            progressbar: None or progress indicator, "STDOUT" print, panel.pane object and panel.widget value setting implemented.
        """
        if str(type(progressbar)).startswith("panel.pane"):
            def setprogress(x, progressbar=progressbar):
                progressbar.object = x
        elif str(type(progressbar)).startswith("panel.widget"):
            def setprogress(x, progressbar=progressbar):
                progressbar.value = x
        elif progressbar == "STDOUT":
            setprogress = lambda x: print(x)
        else:
            setprogress = lambda x: None
        
        setprogress("Data Reading ...")
        self.data_reading(content=content)
        setprogress("Data Processing ...")
        self.processingdf()
        self.quantity_profiles_proteinGroups()
        setprogress("PCA ...")
        self.perform_pca()
        setprogress("Calculating Manhattan Distance ...")
        self.calc_biological_precision()
        setprogress("Assembling analysis output ...")
        self.results_overview_table()
        setprogress("Analysis finished.")
    
    
    def data_reading(self, filename=None, content=None):
        """
        Data import. Can read the df_original from a file or buffer. 
        df_original contains all information of the raw file; tab separated file is imported,

        Args:
            self:
                filename: string
                imported_columns : dictionry; columns that correspond to this regular expression will be imported
            filename: default None, to use the class attribute. Otherwise overwrites the class attribute upon success.
            content: default None, to use the filename. Any valid input to pd.read_csv can be provided, e.g. a StringIO buffer.

        Returns:
            self.df_orginal: raw, unprocessed dataframe, single level column index
        """
        
        # use instance attribute if no filename is provided
        if filename is None:
            filename = self.filename
        # if no buffer is provided for the content read straight from the file
        if content is None:
            content = filename

        if filename.endswith("xls") or filename.endswith("txt"):
            self.df_original = pd.read_csv(content, sep="\t", comment="#", usecols=lambda x: bool(re.match(self.imported_columns, x)), low_memory = True)
        else: #assuming csv file
            self.df_original = pd.read_csv(content, sep=",", comment="#", usecols=lambda x: bool(re.match(self.imported_columns, x)), low_memory = True)
        assert self.df_original.shape[0]>10 and self.df_original.shape[1]>5
        
        self.filename = filename

        return self.df_original
    

    def processingdf(self, name_pattern=None, summed_MSMS_counts=None, consecutiveLFQi=None, RatioHLcount=None, RatioVariability=None, custom_columns=None, custom_normalized=None):
        """
        Analysis of the SILAC/LFQ-MQ/LFQ-Spectronaut data will be performed. The dataframe will be filtered, normalized, and converted into a dataframe, 
        characterized by a flat column index. These tasks is performed by following functions:
            indexingdf(df_original, acquisition_set_dict, acquisition, fraction_dict, name_pattern)
            spectronaut_LFQ_indexingdf(df_original, Spectronaut_columnRenaming, acquisition_set_dict, acquisition, fraction_dict, name_pattern)
            stringency_silac(df_index)
            normalization_01_silac(df_stringency_mapfracstacked):
            logarithmization_silac(df_stringency_mapfracstacked):
            stringency_lfq(df_index):
            normalization_01_lfq(df_stringency_mapfracstacked):
            logarithmization_lfq(df_stringency_mapfracstacked):

        Args:
            self.acquisition: string, "LFQ6 - Spectronaut", "LFQ5 - Spectronaut", "LFQ5 - MQ", "LFQ6 - MQ", "SILAC - MQ"
            additional arguments can be used to override the value set by the class init function

        Returns:
            self:
                map_names: list of Map names
                df_01_stacked: df; 0-1 normalized data with "normalized profile" as column name
                df_log_stacked: df; log transformed data
                analysis_summary_dict["0/1 normalized data - mean"] : 0/1 normalized data across all maps by calculating the mean
                                     ["changes in shape after filtering"]
                                     ["Analysis parameters"] : {"acquisition" : ..., 
                                                                "filename" : ...,
                                                            #SILAC#
                                                                "Ratio H/L count 1 (>=X)" : ...,
                                                                "Ratio H/L count 2 (>=Y, var<Z)" : ...,
                                                                "Ratio variability (<Z, count>=Y)" : ...
                                                            #LFQ# 
                                                                "consecutive data points" : ...,
                                                                "summed MS/MS counts" : ...
                                                               }
        """
        
        if name_pattern is None:
            name_pattern = self.name_pattern
        
        if self.acquisition == "SILAC - MQ":
            if RatioHLcount is None:
                RatioHLcount = self.RatioHLcount
            if RatioVariability is None:
                RatioVariability = self.RatioVariability
        elif self.acquisition == "Custom":
            if custom_columns is None:
                custom_columns = self.custom_columns
            if custom_normalized is None:
                custom_normalized = self.custom_normalized
        else:
            if summed_MSMS_counts is None:
                summed_MSMS_counts = self.summed_MSMS_counts
            if consecutiveLFQi is None:
                consecutiveLFQi = self.consecutiveLFQi
            
        shape_dict = {}
        df_organellarMarkerSet = self.df_organellarMarkerSet.copy()
        df_organellarMarkerSet = df_organellarMarkerSet.rename({"Protein ID": "Protein IDs"}, axis=1).set_index("Protein IDs")
        

        
        def indexingdf():
            """
            For data output from MaxQuant, all columns - except of "MS/MS count" and "LFQ intensity" (LFQ) | "Ratio H/L count", "Ratio H/L variability [%]" 
            (SILAC) - will be set as index. A multiindex will be generated, containing "Set" ("MS/MS count", "LFQ intensity"|  "Ratio H/L count", "Ratio H/L
            variability [%]"), "Fraction" (= defined via "name_pattern") and "Map" (= defined via "name_pattern") as level names, allowing the stacking and 
            unstacking of the dataframe. The dataframe will be filtered by removing matches to the reverse database, matches only identified by site, and 
            potential contaminants.
            
            Args:
                self:
                    df_original: dataframe, columns defined through self.imported_columns
                    acquisition_set_dict: dictionary, all columns will be set as index, except of those that are listed in acquisition_set_dict
                    acquisition: string, one of "LFQ6 - Spectronaut", "LFQ5 - Spectronaut", "LFQ5 - MQ", "LFQ6 - MQ", "SILAC - MQ"
                    fraction_dict: "Fraction" is part of the multiindex; fraction_dict allows the renaming of the fractions e.g. 3K -> 03K
                    name_pattern: regular expression, to identify Map-Fraction-(Replicate)

            Returns:
                self:
                    df_index: mutliindex dataframe, which contains 3 level labels: Map, Fraction, Type
                    shape_dict["Original size"] of df_original
                    shape_dict["Shape after categorical filtering"] of df_index
                    fractions: list of fractions e.g. ["01K", "03K", ...]
            """
            
            df_original = self.df_original.copy()
            df_original.drop("Protein IDs", axis=1, inplace=True)
            df_original.rename({"Majority protein IDs": "Original Protein IDs"}, axis=1, inplace=True)
            df_original.insert(0, "Protein IDs", [split_ids(el) for el in df_original["Original Protein IDs"]])
            if self.reannotate_genes:
                df_original["Gene names"] = reannotate_genes_uniprot(df_original["Protein IDs"], **self.reannotation_source)
            df_original = df_original.set_index([col for col in df_original.columns
                                                 if any([re.match(s, col) for s in self.acquisition_set_dict[self.acquisition]]) == False])
    
            # multindex will be generated, by extracting the information about the Map, Fraction and Type from each individual column name
            multiindex = pd.MultiIndex.from_arrays(
                arrays=[
                    [[re.findall(s, col)[0] for s in self.acquisition_set_dict[self.acquisition] if re.match(s,col)][0]
                    for col in df_original.columns],
                    [re.match(self.name_pattern, col).group("rep") for col in df_original.columns] if not "<cond>" in self.name_pattern
                     else ["_".join(re.match(self.name_pattern, col).group("cond", "rep")) for col in df_original.columns],
                    [re.match(self.name_pattern, col).group("frac") for col in df_original.columns],
                ],
                names=["Set", "Map", "Fraction"]
            )
            
            df_original.columns = multiindex
            df_original.sort_index(1, inplace=True)
            
            shape_dict["Original size"] = df_original.shape
            
            try:
                df_index = df_original.xs(
                    np.nan, 0, "Reverse")
            except:
                pass
            
            try:
                df_index = df_index.xs(
                np.nan, 0, "Potential contaminant")
            except:
                pass
            
            try:
                df_index = df_index.xs(
                np.nan, 0, "Only identified by site")
            except:
                pass
            
            df_index.replace(0, np.nan, inplace=True)
            shape_dict["Shape after categorical filtering"] = df_index.shape
            
            df_index.rename(columns={"MS/MS Count":"MS/MS count"}, inplace=True)

            fraction_wCyt = list(df_index.columns.get_level_values("Fraction").unique())
            
            ##############Cyt should get only be removed if it is not an NMC split
            if "Cyt" in fraction_wCyt and len(fraction_wCyt) >= 4:
                df_index.drop("Cyt", axis=1, level="Fraction", inplace=True)
            try:
                if self.acquisition == "LFQ5 - MQ":
                    df_index.drop("01K", axis=1, level="Fraction", inplace=True)
            except:
                pass
            
            self.fractions = natsort.natsorted(list(df_index.columns.get_level_values("Fraction").unique()))
            
            # merge with markerset
            df_index.columns = df_index.columns.values
            df_index = df_index.join(df_organellarMarkerSet, how="left", on="Protein IDs").set_index("Compartment", append=True)
            df_index.columns = pd.MultiIndex.from_tuples(df_index.columns, names=["Set", "Map", "Fraction"])
            self.df_index = df_index
            
            return df_index
        
        
        
        def custom_indexing_and_normalization():
            df_original = self.df_original.copy()
            df_original.rename({custom_columns["ids"]: "Protein IDs", custom_columns["genes"]: "Gene names"}, axis=1, inplace=True)
            if self.reannotate_genes:
                df_original["Gene names"] = reannotate_genes_uniprot(df_original["Protein IDs"], **self.reannotation_source)
            df_original = df_original.set_index([col for col in df_original.columns
                                                 if any([re.match(s, col) for s in self.acquisition_set_dict[self.acquisition]]) == False])
    
            # multindex will be generated, by extracting the information about the Map, Fraction and Type from each individual column name
            multiindex = pd.MultiIndex.from_arrays(
                arrays=[
                    ["normalized profile" for col in df_original.columns],
                    [re.match(self.name_pattern, col).group("rep") for col in df_original.columns] if not "<cond>" in self.name_pattern
                     else ["_".join(re.match(self.name_pattern, col).group("cond", "rep")) for col in df_original.columns],
                    [re.match(self.name_pattern, col).group("frac") for col in df_original.columns],
                ],
                names=["Set", "Map", "Fraction"]
            )
            
            df_original.columns = multiindex
            df_original.sort_index(1, inplace=True)
            
            shape_dict["Original size"] = df_original.shape
            
            # for custom upload assume full normalization for now. this should be extended to valid value filtering and 0-1 normalization later
            df_index = df_original.copy()
            self.fractions = natsort.natsorted(list(df_index.columns.get_level_values("Fraction").unique()))
            
            # merge with markerset
            df_index.columns = df_index.columns.values
            df_index = df_index.join(df_organellarMarkerSet, how="left", on="Protein IDs").set_index("Compartment", append=True)
            df_index.columns = pd.MultiIndex.from_tuples(df_index.columns, names=["Set", "Map", "Fraction"])
            self.df_index = df_index
            
            return df_index
        
        
        def spectronaut_LFQ_indexingdf():
            """
            For data generated from the Spectronaut software, columns will be renamed, such it fits in the scheme of MaxQuant output data.  Subsequently, all
            columns - except of "MS/MS count" and "LFQ intensity" will be set as index. A multiindex will be generated, containing "Set" ("MS/MS count" and 
            "LFQ intensity"), Fraction" and "Map" (= defined via "name_pattern"; both based on the column name R.condition - equivalent to the column name "Map" 
            in df_renamed["Map"]) as level labels.
            !!!
            !!!It is very important to define R.Fraction, R.condition already during the setup of Spectronaut!!!
            !!!
            
            Args:
                self:
                    df_original: dataframe, columns defined through self.imported_columns
                    Spectronaut_columnRenaming
                    acquisition_set_dict: dictionary, all columns will be set as index, except of those that are listed in acquisition_set_dict
                    acquisition: string, "LFQ6 - Spectronaut", "LFQ5 - Spectronaut"
                    fraction_dict: "Fraction" is part of the multiindex; fraction_dict allows the renaming of the fractions e.g. 3K -> 03K
                    name_pattern: regular expression, to identify Map-Fraction-(Replicate)
                
            Returns:
                self:
                    df_index: mutliindex dataframe, which contains 3 level labels: Map, Fraction, Type
                    shape_dict["Original size"] of df_index
                    fractions: list of fractions e.g. ["01K", "03K", ...]
            """
            
            df_original = self.df_original.copy()
            
            df_renamed = df_original.rename(columns=self.Spectronaut_columnRenaming)
            df_renamed.insert(0, "Protein IDs", [split_ids(el) for el in df_renamed["Original Protein IDs"]])
            if self.reannotate_genes:
                df_renamed["Gene names"] = reannotate_genes_uniprot(df_renamed["Protein IDs"], **self.reannotation_source)
            
            df_renamed["Fraction"] = [re.match(self.name_pattern, i).group("frac") for i in df_renamed["Map"]]
            df_renamed["Map"] = [re.match(self.name_pattern, i).group("rep") for i in df_renamed["Map"]] if not "<cond>" in self.name_pattern else ["_".join(
                        re.match(self.name_pattern, i).group("cond", "rep")) for i in df_renamed["Map"]]
            
            df_index = df_renamed.set_index([col for col in df_renamed.columns if any([re.match(s, col) for s in self.acquisition_set_dict[self.acquisition]])==False])
            
            df_index.columns.names = ["Set"]
            
            # In case fractionated data was used this needs to be catched and aggregated
            try:
                df_index = df_index.unstack(["Map", "Fraction"])
            except ValueError:
                df_index = df_index.groupby(by=df_index.index.names).agg(np.nansum, axis=0)
                df_index = df_index.unstack(["Map", "Fraction"])
            
            df_index.replace(0, np.nan, inplace=True)
            shape_dict["Original size"]=df_index.shape
            
            fraction_wCyt = list(df_index.columns.get_level_values("Fraction").unique())
            #Cyt is removed only if it is not an NMC split
            if "Cyt" in fraction_wCyt and len(fraction_wCyt) >= 4:
                df_index.drop("Cyt", axis=1, level="Fraction", inplace=True)
            try:
                if self.acquisition == "LFQ5 - Spectronaut":
                    df_index.drop("01K", axis=1, level="Fraction", inplace=True)
            except:
                pass
            
            self.fractions = natsort.natsorted(list(df_index.columns.get_level_values("Fraction").unique()))
            
            # merge with markerset
            df_index.columns = df_index.columns.values
            df_index = df_index.join(df_organellarMarkerSet, how="left", on="Protein IDs").set_index("Compartment", append=True)
            df_index.columns = pd.MultiIndex.from_tuples(df_index.columns, names=["Set", "Map", "Fraction"])
            
            self.df_index = df_index
            
            return df_index
        
        
        def stringency_silac(df_index):
            """
            The multiindex dataframe is subjected to stringency filtering. Only Proteins with complete profiles are considered (a set of f.e. 5 SILAC ratios 
            in case you have 5 fractions / any proteins with missing values were rejected). Proteins were retained with 3 or more quantifications in each 
            subfraction (=count). Furthermore, proteins with only 2 quantification events in one or more subfraction were retained, if their ratio variability for 
            ratios obtained with 2 quantification events was below 30% (=var). SILAC ratios were linearly normalized by division through the fraction median. 
            Subsequently normalization to SILAC loading was performed.Data is annotated based on specified marker set e.g. eLife.

            Args:
                df_index: multiindex dataframe, which contains 3 level labels: MAP, Fraction, Type
                RatioHLcount: int, 2
                RatioVariability: int, 30 
                df_organellarMarkerSet: df, columns: "Protein ID", "Compartment", no index 
                fractions: list of fractions e.g. ["01K", "03K", ...]

            Returns:
                df_stringency_mapfracstacked: dataframe, in which "MAP" and "Fraction" are stacked;
                                              columns "Ratio H/L count", "Ratio H/L variability [%]", and "Ratio H/L" stored as single level indices
                shape_dict["Shape after Ratio H/L count (>=3)/var (count>=2, var<30) filtering"] of df_countvarfiltered_stacked
                shape_dict["Shape after filtering for complete profiles"] of df_stringency_mapfracstacked
            """

            # Fraction and Map will be stacked
            df_stack = df_index.stack(["Fraction", "Map"])

            # filtering for sufficient number of quantifications (count in "Ratio H/L count"), taken variability (var in Ratio H/L variability [%]) into account
            # zip: allows direct comparison of count and var
            # only if the filtering parameters are fulfilled the data will be introduced into df_countvarfiltered_stacked
            #default setting: RatioHLcount = 2 ; RatioVariability = 30
            
            df_countvarfiltered_stacked = df_stack.loc[[count>RatioHLcount or (count==RatioHLcount and var<RatioVariability) 
                                            for var, count in zip(df_stack["Ratio H/L variability [%]"], df_stack["Ratio H/L count"])]]
            
            shape_dict["Shape after Ratio H/L count (>=3)/var (count==2, var<30) filtering"] = df_countvarfiltered_stacked.unstack(["Fraction", "Map"]).shape

            # "Ratio H/L":normalization to SILAC loading, each individual experiment (FractionXMap) will be divided by its median
            # np.median([...]): only entries, that are not NANs are considered
            df_normsilac_stacked = df_countvarfiltered_stacked["Ratio H/L"]\
                .unstack(["Fraction", "Map"])\
                .apply(lambda x: x/np.nanmedian(x), axis=0)\
                .stack(["Map", "Fraction"])

            df_stringency_mapfracstacked = df_countvarfiltered_stacked[["Ratio H/L count", "Ratio H/L variability [%]"]].join(
                pd.DataFrame(df_normsilac_stacked, columns=["Ratio H/L"]))

            # dataframe is grouped (Map, id), that allows the filtering for complete profiles
            df_stringency_mapfracstacked = df_stringency_mapfracstacked.groupby(["Map", "id"]).filter(lambda x: len(x)>=len(self.fractions))
            
            shape_dict["Shape after filtering for complete profiles"]=df_stringency_mapfracstacked.unstack(["Fraction", "Map"]).shape
            
            # Ratio H/L is converted into Ratio L/H
            df_stringency_mapfracstacked["Ratio H/L"] = df_stringency_mapfracstacked["Ratio H/L"].transform(lambda x: 1/x)
            
            #Annotation with marker genes
            df_organellarMarkerSet = self.df_organellarMarkerSet
            
            df_stringency_mapfracstacked.reset_index(inplace=True)
            df_stringency_mapfracstacked.set_index([c for c in df_stringency_mapfracstacked.columns
                                                    if c not in ["Ratio H/L count","Ratio H/L variability [%]","Ratio H/L"]], inplace=True)
            df_stringency_mapfracstacked.rename(index={np.nan:"undefined"}, level="Compartment", inplace=True)

            return df_stringency_mapfracstacked


        def normalization_01_silac(df_stringency_mapfracstacked):
            """
            The multiindex dataframe, that was subjected to stringency filtering, is 0-1 normalized ("Ratio H/L").

            Args:
                df_stringency_mapfracstacked: dataframe, in which "Map" and "Fraction" are stacked; 
                                              columns "Ratio H/L count", "Ratio H/L variability [%]", and "Ratio H/L" stored as single level indices
                self:
                    fractions: list of fractions e.g. ["01K", "03K", ...]
                    data_completeness: series, for each individual map, as well as combined maps: 1 - (percentage of NANs)

            Returns:
                df_01_stacked: dataframe, in which "MAP" and "Fraction" are stacked; data in the column "Ratio H/L" is 0-1 normalized and renamed to "normalized
                               profile"; the columns "Ratio H/L count", "Ratio H/L variability [%]", and "normalized profile" stored as single level indices; 
                               plotting is possible now
                self:
                    analysis_summary_dict["Data/Profile Completeness"] : df, with information about Data/Profile Completeness
                                        column: "Experiment", "Map", "Data completeness", "Profile completeness"
                                        no row index
            """

            df_01norm_unstacked = df_stringency_mapfracstacked["Ratio H/L"].unstack("Fraction")

            # 0:1 normalization of Ratio L/H
            df_01norm_unstacked = df_01norm_unstacked.div(df_01norm_unstacked.sum(axis=1), axis=0)

            df_01_stacked = df_stringency_mapfracstacked[["Ratio H/L count", "Ratio H/L variability [%]"]].join(pd.DataFrame
                (df_01norm_unstacked.stack("Fraction"),columns=["Ratio H/L"]))

            # "Ratio H/L" will be renamed to "normalized profile"
            df_01_stacked.columns = [col if col!="Ratio H/L" else "normalized profile" for col in df_01_stacked.columns]

            return df_01_stacked


        def logarithmization_silac(df_stringency_mapfracstacked):
            """
            The multiindex dataframe, that was subjected to stringency filtering, is logarithmized ("Ratio H/L").

            Args:
                df_stringency_mapfracstacked: dataframe, in which "MAP" and "Fraction" are stacked; the columns "Ratio H/L count", "Ratio H/L variability [%]", 
                                              and "Ratio H/L" stored as single level indices

            Returns:
                df_log_stacked: dataframe, in which "MAP" and "Fraction" are stacked; data in the column "log profile" originates from logarithmized "Ratio H/L" 
                                data; the columns "Ratio H/L count", "Ratio H/L variability [%]" and  "log profile" are stored as single level indices; 
                                PCA is possible now

            """
            
            # logarithmizing, basis of 2
            df_lognorm_ratio_stacked = df_stringency_mapfracstacked["Ratio H/L"].transform(np.log2)
            df_log_stacked = df_stringency_mapfracstacked[["Ratio H/L count", "Ratio H/L variability [%]"]].join(
                pd.DataFrame(df_lognorm_ratio_stacked, columns=["Ratio H/L"]))

            # "Ratio H/L" will be renamed to "log profile"
            df_log_stacked.columns = [col if col !="Ratio H/L" else "log profile" for col in df_log_stacked.columns]

            return df_log_stacked

        
        def stringency_lfq(df_index):
            """
            The multiindex dataframe is subjected to stringency filtering. Only Proteins which were identified with
            at least [4] consecutive data points regarding the "LFQ intensity", and if summed MS/MS counts >= n(fractions)*[2]
            (LFQ5: min 10 and LFQ6: min 12, respectively; coverage filtering) were included.
            Data is annotated based on specified marker set e.g. eLife.

            Args:
                df_index: multiindex dataframe, which contains 3 level labels: MAP, Fraction, Typ
                self:
                    df_organellarMarkerSet: df, columns: "Protein ID", "Compartment", no index
                    fractions: list of fractions e.g. ["01K", "03K", ...]
                    summed_MSMS_counts: int, 2
                    consecutiveLFQi: int, 4

            Returns:
                df_stringency_mapfracstacked: dataframe, in which "Map" and "Fraction" is stacked; "LFQ intensity" and "MS/MS count" define a 
                                              single-level column index
                
                self:
                    shape_dict["Shape after MS/MS value filtering"] of df_mscount_mapstacked
                    shape_dict["Shape after consecutive value filtering"] of df_stringency_mapfracstacked
                """

            df_index = df_index.stack("Map")

            # sorting the level 0, in order to have LFQ intensity -    MS/MS count instead of continuous alternation
            df_index.sort_index(axis=1, level=0, inplace=True)
            
            # "MS/MS count"-column: take the sum over the fractions; if the sum is larger than n[fraction]*2, it will be stored in the new dataframe
            minms = (len(self.fractions) * self.summed_MSMS_counts)
            
            if minms > 0:
                df_mscount_mapstacked = df_index.loc[df_index[("MS/MS count")].apply(np.sum, axis=1) >= minms]
    
                shape_dict["Shape after MS/MS value filtering"]=df_mscount_mapstacked.unstack("Map").shape
                df_stringency_mapfracstacked = df_mscount_mapstacked.copy()
            else:
                df_stringency_mapfracstacked = df_index.copy()

            # series no dataframe is generated; if there are at least i.e. 4 consecutive non-NANs, data will be retained
            df_stringency_mapfracstacked.sort_index(level="Fraction", axis=1, key=natsort_index_keys, inplace=True)
            df_stringency_mapfracstacked = df_stringency_mapfracstacked.loc[
                df_stringency_mapfracstacked[("LFQ intensity")]\
                .apply(lambda x: np.isfinite(x), axis=0)\
                .apply(lambda x: sum(x) >= self.consecutiveLFQi and any(x.rolling(window=self.consecutiveLFQi).sum() >= self.consecutiveLFQi), axis=1)]
            
            shape_dict["Shape after consecutive value filtering"]=df_stringency_mapfracstacked.unstack("Map").shape

            df_stringency_mapfracstacked = df_stringency_mapfracstacked.copy().stack("Fraction")
            
            #Annotation with marker genes
            df_organellarMarkerSet = self.df_organellarMarkerSet
            
            df_stringency_mapfracstacked.reset_index(inplace=True)
            df_stringency_mapfracstacked.set_index([c for c in df_stringency_mapfracstacked.columns
                                                    if c!="MS/MS count" and c!="LFQ intensity"], inplace=True)
            df_stringency_mapfracstacked.rename(index={np.nan : "undefined"}, level="Compartment", inplace=True)
            
            return df_stringency_mapfracstacked


        def normalization_01_lfq(df_stringency_mapfracstacked):
            """
            The multiindex dataframe, that was subjected to stringency filtering, is 0-1 normalized ("LFQ intensity").

            Args:
                df_stringency_mapfracstacked: dataframe, in which "Map" and "Fraction" is stacked, "LFQ intensity" and "MS/MS count" define a 
                                              single-level column index
                self:
                    fractions: list of fractions e.g. ["01K", "03K", ...]


            Returns:
                df_01_stacked: dataframe, in which "MAP" and "Fraction" are stacked; data in the column "LFQ intensity" is 0-1 normalized and renamed to 
                               "normalized profile"; the columns "normalized profile" and "MS/MS count" are stored as single level indices; plotting is possible now

            """
            df_01norm_mapstacked = df_stringency_mapfracstacked["LFQ intensity"].unstack("Fraction")
            
            # 0:1 normalization of Ratio L/H
            df_01norm_unstacked = df_01norm_mapstacked.div(df_01norm_mapstacked.sum(axis=1), axis=0)
            
            df_rest = df_stringency_mapfracstacked.drop("LFQ intensity", axis=1)
            df_01_stacked = df_rest.join(pd.DataFrame(df_01norm_unstacked.stack(
                   "Fraction"),columns=["LFQ intensity"]))
            
            # rename columns: "LFQ intensity" into "normalized profile"
            df_01_stacked.columns = [col if col!="LFQ intensity" else "normalized profile" for col in
                                     df_01_stacked.columns]
            
            #imputation 
            df_01_stacked = df_01_stacked.unstack("Fraction").replace(np.NaN, 0).stack("Fraction")
            
            df_01_stacked = df_01_stacked.sort_index()
            
            return df_01_stacked


        def logarithmization_lfq(df_stringency_mapfracstacked):
            """The multiindex dataframe, that was subjected to stringency filtering, is logarithmized ("LFQ intensity").

            Args:
                df_stringency_mapfracstacked: dataframe, in which "Map" and "Fraction" is stacked; "LFQ intensity" and "MS/MS count" define a 
                                              single-level column index

            Returns:
                df_log_stacked: dataframe, in which "MAP" and "Fraction" are stacked; data in the column "log profile" originates from logarithmized 
                                "LFQ intensity"; the columns "log profile" and "MS/MS count" are stored as single level indices; PCA is possible now
            """

            df_lognorm_ratio_stacked = df_stringency_mapfracstacked["LFQ intensity"].transform(np.log2)
            
            df_rest = df_stringency_mapfracstacked.drop("LFQ intensity", axis=1)
            df_log_stacked = df_rest.join(pd.DataFrame(df_lognorm_ratio_stacked, columns=["LFQ intensity"]))
            
            # "LFQ intensity" will be renamed to "log profile"
            df_log_stacked.columns = [col if col!="LFQ intensity" else "log profile" for col in df_log_stacked.columns]

            return df_log_stacked
        
        if self.acquisition == "SILAC - MQ":
            
            # Index data
            df_index = indexingdf()
            map_names = df_index.columns.get_level_values("Map").unique()
            self.map_names = map_names
            
            # Run stringency filtering and normalization
            df_stringency_mapfracstacked = stringency_silac(df_index)
            self.df_stringencyFiltered = df_stringency_mapfracstacked
            self.df_01_stacked = normalization_01_silac(df_stringency_mapfracstacked)
            self.df_log_stacked = logarithmization_silac(df_stringency_mapfracstacked)
            
            # format and reduce 0-1 normalized data for comparison with other experiments
            df_01_comparison = self.df_01_stacked.copy()
            df_01_comparison.drop(["Ratio H/L count", "Ratio H/L variability [%]"], inplace=True, axis=1)
            df_01_comparison = df_01_comparison.unstack(["Map", "Fraction"])
            df_01_comparison.columns = ["?".join(el) for el in df_01_comparison.columns.values]
            df_01_comparison = df_01_comparison.copy().reset_index().drop(["C-Score", "Q-value", "Score", "Majority protein IDs", "Protein names", "id"], axis=1, errors="ignore")
            
            # poopulate analysis summary dictionary with (meta)data
            self.analysis_summary_dict["0/1 normalized data"] = df_01_comparison.to_json()
            self.analysis_summary_dict["changes in shape after filtering"] = shape_dict.copy() 
            analysis_parameters = {"acquisition" : self.acquisition, 
                                   "filename" : self.filename,
                                   "comment" : self.comment,
                                   "Ratio H/L count" : self.RatioHLcount,
                                   "Ratio variability" : self.RatioVariability,
                                   "organism" : self.organism,
                                  }
            self.analysis_summary_dict["Analysis parameters"] = analysis_parameters.copy()
            
            # TODO this line needs to be removed.
            self.analysed_datasets_dict[self.expname] = self.analysis_summary_dict.copy()


        elif self.acquisition == "LFQ5 - MQ" or self.acquisition == "LFQ6 - MQ" or self.acquisition == "LFQ5 - Spectronaut" or self.acquisition == "LFQ6 - Spectronaut":
        
            #if not summed_MS_counts:
            #    summed_MS_counts = self.summed_MS_counts
            #if not consecutiveLFQi:
            #    consecutiveLFQi = self.consecutiveLFQi
        
            if self.acquisition == "LFQ5 - MQ" or self.acquisition == "LFQ6 - MQ":
                df_index = indexingdf()
            elif self.acquisition == "LFQ5 - Spectronaut" or self.acquisition == "LFQ6 - Spectronaut":
                df_index = spectronaut_LFQ_indexingdf()
            
            map_names = df_index.columns.get_level_values("Map").unique()
            self.map_names = map_names
            
            df_stringency_mapfracstacked = stringency_lfq(df_index)
            self.df_stringencyFiltered = df_stringency_mapfracstacked
            self.df_log_stacked = logarithmization_lfq(df_stringency_mapfracstacked)
            self.df_01_stacked = normalization_01_lfq(df_stringency_mapfracstacked)
            
            df_01_comparison = self.df_01_stacked.copy()
            df_01_comparison.drop("MS/MS count", inplace=True, axis=1, errors="ignore")
            df_01_comparison = df_01_comparison.unstack(["Map", "Fraction"])
            df_01_comparison.columns = ["?".join(el) for el in df_01_comparison.columns.values]
            df_01_comparison = df_01_comparison.copy().reset_index().drop(["C-Score", "Q-value", "Score", "Majority protein IDs", "Protein names", "id"], axis=1, errors="ignore")
            self.analysis_summary_dict["0/1 normalized data"] = df_01_comparison.to_json()#double_precision=4) #.reset_index()
            
            self.analysis_summary_dict["changes in shape after filtering"] = shape_dict.copy() 
            analysis_parameters = {"acquisition" : self.acquisition, 
                                   "filename" : self.filename,
                                   "comment" : self.comment,
                                   "consecutive data points" : self.consecutiveLFQi,
                                   "summed MS/MS counts" : self.summed_MSMS_counts,
                                   "organism" : self.organism,
                                  }
            self.analysis_summary_dict["Analysis parameters"] = analysis_parameters.copy() 
            self.analysed_datasets_dict[self.expname] = self.analysis_summary_dict.copy()
            #return self.df_01_stacked
        elif self.acquisition == "Custom":
            df_index = custom_indexing_and_normalization()
            map_names = df_index.columns.get_level_values("Map").unique()
            self.map_names = map_names
            df_01_stacked = df_index.stack(["Map", "Fraction"]).reset_index()
            df_01_stacked.set_index([c for c in df_01_stacked.columns if c not in ["normalized profile"]], inplace=True)
            df_01_stacked.rename(index={np.nan:"undefined"}, level="Compartment", inplace=True)
            self.df_01_stacked = df_01_stacked
            
            df_01_comparison = self.df_01_stacked.copy()
            df_01_comparison.drop("MS/MS count", inplace=True, axis=1, errors="ignore")
            df_01_comparison = df_01_comparison.unstack(["Map", "Fraction"])
            df_01_comparison.columns = ["?".join(el) for el in df_01_comparison.columns.values]
            df_01_comparison = df_01_comparison.copy().reset_index().drop(["C-Score", "Q-value", "Score", "Majority protein IDs", "Protein names", "id"], axis=1, errors="ignore")
            self.analysis_summary_dict["0/1 normalized data"] = df_01_comparison.to_json()#double_precision=4) #.reset_index()
            
            self.analysis_summary_dict["changes in shape after filtering"] = shape_dict.copy() 
            analysis_parameters = {"acquisition" : self.acquisition, 
                                   "filename" : self.filename,
                                   "comment" : self.comment,
                                   "organism" : self.organism,
                                  }
            self.analysis_summary_dict["Analysis parameters"] = analysis_parameters.copy() 
            self.analysed_datasets_dict[self.expname] = self.analysis_summary_dict.copy()
            
        else:
            return "I do not know this"
    
    
    def plot_log_data(self):     
        """
        
        Args:
            self.df_log_stacked
        
        
        Returns:
            log_histogram: Histogram of log transformed data
        
        """  
      
        log_histogram = px.histogram(self.df_log_stacked.reset_index().sort_values(["Map", "Fraction"], key=natsort_list_keys), 
                                     x="log profile",
                                     facet_col="Fraction",
                                     facet_row="Map",
                                     template="simple_white",
                                     labels={"log profile": "log tranformed data ({})".format("LFQ intenisty" if self.acquisition != "SILAC - MQ" else "Ratio H/L")}
                                    )
        
        log_histogram.for_each_xaxis(lambda axis: axis.update(title={"text":""}))
        log_histogram.for_each_yaxis(lambda axis: axis.update(title={"text":""}))
        log_histogram.add_annotation(x=0.5, y=0, yshift=-50, xref="paper",showarrow=False, yref="paper",
                    text="log2(LFQ intensity)")
        log_histogram.add_annotation(x=0, y=0.5, textangle=270, xref="paper",showarrow=False, yref="paper", xshift=-50,
                    text="count")
        log_histogram.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        
        return log_histogram
    
    
    def quantity_profiles_proteinGroups(self):
        """
        Number of profiles, protein groups per experiment, and the data completness of profiles (total quantity, intersection) is calculated.
        
        Args:
            self:
                acquisition: string, "LFQ6 - Spectronaut", "LFQ5 - Spectronaut", "LFQ5 - MQ", "LFQ6 - MQ", "SILAC - MQ"
                df_index: multiindex dataframe, which contains 3 level labels: MAP, Fraction, Typ
                df_01_stacked: df; 0-1 normalized data with "normalized profile" as column name
            
        Returns:
            self:
                df_quantity_pr_pg: df; no index, columns: "filtering", "type", "npg", "npr", "npr_dc"; containign following information:
                    npg_t: protein groups per experiment total quantity
                    npgf_t = groups with valid profiles per experiment total quanitity
                    npr_t: profiles with any valid values
                    nprf_t = total number of valid profiles
                    
                    npg_i: protein groups per experiment intersection
                    npgf_i = groups with valid profiles per experiment intersection
                    npr_i: profiles with any valid values in the intersection
                    nprf_i = total number of valid profiles in the intersection
                    
                    npr_t_dc: profiles, % values != nan
                    nprf_t_dc = profiles, total, filtered, % values != nan
                    npr_i_dc: profiles, intersection, % values != nan
                    nprf_i_dc = profiles, intersection, filtered, % values != nan
                    
                df_npg | df_npgf: index: maps e.g. "Map1", "Map2",..., columns: fractions e.g. "03K", "06K", ...
                    npg_f = protein groups, per fraction
                    or npgf_f = protein groups, filtered, per fraction
                df_npg_dc | df_npgf_dc: index: maps e.g. "Map1", "Map2",..., columns: fractions e.g. "03K", "06K", ...
                    npg_f_dc = protein groups, per fraction, % values != nan 
                    or npgf_f_dc = protein groups, filtered, per fraction, % values != nan                        
                    
        """
    
        if self.acquisition == "SILAC - MQ":
            df_index = self.df_index["Ratio H/L"]
            df_01_stacked = self.df_01_stacked["normalized profile"]
        elif self.acquisition.startswith("LFQ"):
            df_index = self.df_index["LFQ intensity"]
            df_01_stacked = self.df_01_stacked["normalized profile"].replace(0, np.nan)
        elif self.acquisition == "Custom":
            df_index = self.df_index["normalized profile"]
            df_01_stacked = self.df_01_stacked["normalized profile"].replace(0, np.nan)
        
        #unfiltered
        npg_t = df_index.shape[0]
        df_index_MapStacked = df_index.stack("Map")
        npr_t = df_index_MapStacked.shape[0]/len(self.map_names)
        npr_t_dc = 1-df_index_MapStacked.isna().sum().sum()/np.prod(df_index_MapStacked.shape)
        
        #filtered
        npgf_t = df_01_stacked.unstack(["Map", "Fraction"]).shape[0]
        df_01_MapStacked = df_01_stacked.unstack("Fraction")
        nprf_t = df_01_MapStacked.shape[0]/len(self.map_names)
        nprf_t_dc = 1-df_01_MapStacked.isna().sum().sum()/np.prod(df_01_MapStacked.shape)
        
        #unfiltered intersection
        try:
            df_index_intersection = df_index_MapStacked.groupby(level="Sequence").filter(lambda x : len(x)==len(self.map_names))
        except:
            df_index_intersection = df_index_MapStacked.groupby(level="Protein IDs").filter(lambda x : len(x)==len(self.map_names))
        npr_i = df_index_intersection.shape[0]/len(self.map_names)
        npr_i_dc = 1-df_index_intersection.isna().sum().sum()/np.prod(df_index_intersection.shape)
        npg_i = df_index_intersection.unstack("Map").shape[0]
        
        #filtered intersection
        try:
            df_01_intersection = df_01_MapStacked.groupby(level = "Sequence").filter(lambda x : len(x)==len(self.map_names))
        except:
            df_01_intersection = df_01_MapStacked.groupby(level = "Protein IDs").filter(lambda x : len(x)==len(self.map_names))
        nprf_i = df_01_intersection.shape[0]/len(self.map_names)
        nprf_i_dc = 1-df_01_intersection.isna().sum().sum()/np.prod(df_01_intersection.shape)
        npgf_i = df_01_intersection.unstack("Map").shape[0]
        
        # summarize in dataframe and save to attribute
        df_quantity_pr_pg = pd.DataFrame(
            {
            "filtering": pd.Series(["before filtering", "before filtering", "after filtering", "after filtering"], dtype=np.dtype("O")),
            "type": pd.Series(["total", "intersection", "total", "intersection"], dtype=np.dtype("O")),
            "number of protein groups": pd.Series([npg_t, npg_i, npgf_t, npgf_i], dtype=np.dtype("float")),
            "number of profiles": pd.Series([npr_t, npr_i, nprf_t, nprf_i], dtype=np.dtype("float")),
            "data completeness of profiles": pd.Series([npr_t_dc, npr_i_dc, nprf_t_dc, nprf_i_dc], dtype=np.dtype("float"))})
        
        self.df_quantity_pr_pg = df_quantity_pr_pg.reset_index()
        self.analysis_summary_dict["quantity: profiles/protein groups"] = self.df_quantity_pr_pg.to_json() 
        
        #additional depth assessment per fraction
        dict_npgf = {}
        dict_npg = {}
        list_npg_dc = []
        list_npgf_dc = []
        
        for df_intersection in [df_index_intersection, df_01_intersection]:
            for fraction in self.fractions:
                df_intersection_frac = df_intersection[fraction]
                npgF_f_dc = 1-df_intersection_frac.isna().sum()/len(df_intersection_frac)
                npgF_f = df_intersection_frac.unstack("Map").isnull().sum(axis=1).value_counts()
                if fraction not in dict_npg.keys():
                    dict_npg[fraction] = npgF_f
                    list_npg_dc.append(npgF_f_dc)
                else:
                    dict_npgf[fraction] = npgF_f
                    list_npgf_dc.append(npgF_f_dc)
        
        df_npg = pd.DataFrame(dict_npg)
        df_npg.index.name =  "Protein Groups present in:"
        df_npg.rename_axis("Fraction", axis=1, inplace=True)
        df_npg = df_npg.stack("Fraction").reset_index()
        df_npg = df_npg.rename({0: "Protein Groups"}, axis=1)
        df_npg.sort_values(["Fraction", "Protein Groups present in:"], inplace=True, key=natsort_list_keys)     
        
        df_npgf = pd.DataFrame(dict_npgf)
        df_npgf.index.name =  "Protein Groups present in:"
        df_npgf.rename_axis("Fraction", axis=1, inplace=True)
        df_npgf = df_npgf.stack("Fraction").reset_index()
        df_npgf = df_npgf.rename({0: "Protein Groups"}, axis=1)
        df_npgf.sort_values(["Fraction", "Protein Groups present in:"], inplace=True, key=natsort_list_keys)
        
        max_df_npg = df_npg["Protein Groups present in:"].max()
        min_df_npg = df_npg["Protein Groups present in:"].min()
        rename_numOFnans =  {}
        for x, y in zip(range(max_df_npg,min_df_npg-1, -1), range(max_df_npg+1)):
            if y == 1:
                rename_numOFnans[x] = "{} Map".format(y)
            elif y == 0:
                rename_numOFnans[x] = "PG not identified".format(y)
            else:
                rename_numOFnans[x] = "{} Maps".format(y)
        for keys in rename_numOFnans.keys():
            df_npg.loc[df_npg["Protein Groups present in:"] ==keys, "Protein Groups present in:"] = rename_numOFnans[keys]
            df_npgf.loc[df_npgf["Protein Groups present in:"] ==keys, "Protein Groups present in:"] = rename_numOFnans[keys]
        
        # summarize in dataframe and save to attributes
        self.df_npg_dc = pd.DataFrame(
            {
            "Fraction" : pd.Series(self.fractions),
            "Data completeness before filtering": pd.Series(list_npg_dc),
            "Data completeness after filtering": pd.Series(list_npgf_dc),
            })
        
        self.df_npg = df_npg
        self.df_npgf = df_npgf
    
    
    def plot_quantity_profiles_proteinGroups(self):
        """
        
        Args:
            self:
                df_quantity_pr_pg: df; no index, columns: "filtering", "type", "npg", "npr", "npr_dc"; further information: see above
                
        Returns:
            
        
        """
        df_quantity_pr_pg = self.df_quantity_pr_pg
        
        layout = go.Layout(barmode="overlay", 
                           xaxis_tickangle=90, 
                           autosize=False,
                           width=300,
                           height=500,
                           xaxis=go.layout.XAxis(linecolor="black",
                                                 linewidth=1,
                                                 #title="Map",
                                                 mirror=True),
                           yaxis=go.layout.YAxis(linecolor="black",
                                                 linewidth=1,
                                                 mirror=True),
                           template="simple_white")
        
        fig_npg = go.Figure()
        for t in df_quantity_pr_pg["type"].unique():
            plot_df = df_quantity_pr_pg[df_quantity_pr_pg["type"] == t]
            fig_npg.add_trace(go.Bar(
                x=plot_df["filtering"],
                y=plot_df["number of protein groups"],
                name=t))
        fig_npg.update_layout(layout, title="Number of Protein Groups", yaxis=go.layout.YAxis(title="Protein Groups"))
        
                
        fig_npr = go.Figure()
        for t in df_quantity_pr_pg["type"].unique():
            plot_df = df_quantity_pr_pg[df_quantity_pr_pg["type"] == t]
            fig_npr.add_trace(go.Bar(
                x=plot_df["filtering"],
                y=plot_df["number of profiles"],
                name=t))
        fig_npr.update_layout(layout, title="Number of Profiles")
        
        df_quantity_pr_pg = df_quantity_pr_pg.sort_values("filtering")
        fig_npr_dc = go.Figure()
        for t in df_quantity_pr_pg["filtering"].unique():
            plot_df = df_quantity_pr_pg[df_quantity_pr_pg["filtering"] == t]
            fig_npr_dc.add_trace(go.Bar(
                x=plot_df["type"],
                y=plot_df["data completeness of profiles"],
                name=t))
        fig_npr_dc.update_layout(layout, title="Coverage", yaxis=go.layout.YAxis(title="Data completness"))
        #fig_npr_dc.update_xaxes(tickangle=30)
        
        fig_npg_F = px.bar(self.df_npg,
                          x="Fraction", 
                          y="Protein Groups", 
                          color="Protein Groups present in:",
                          template="simple_white",
                          title = "Protein groups per fraction - before filtering",
                          width=500)
        
        fig_npgf_F = px.bar(self.df_npgf,
                  x="Fraction", 
                  y="Protein Groups", 
                  color="Protein Groups present in:",
                  template="simple_white",
                  title = "Protein groups per fraction - after filtering",
                  width=500)
        
        fig_npg_F_dc = go.Figure()
        for data_type in ["Data completeness after filtering", "Data completeness before filtering"]:
            fig_npg_F_dc.add_trace(go.Bar(
                    x=self.df_npg_dc["Fraction"],
                    y=self.df_npg_dc[data_type],
                    name=data_type))
        fig_npg_F_dc.update_layout(layout, barmode="overlay", title="Data completeness per fraction", yaxis=go.layout.YAxis(title=""), height=450, width=600)
        
        
        return fig_npg, fig_npr, fig_npr_dc, fig_npg_F, fig_npgf_F, fig_npg_F_dc
                                                  
                                                                                
    def perform_pca(self):
        """
        PCA will be performed, using logarithmized data.

        Args:
            self:
                markerproteins: dictionary, key: cluster name, value: gene names (e.g. {"Proteasome" : ["PSMA1", "PSMA2",...], ...}
            "V-type proton ATP
                df_log_stacked: dataframe, in which "MAP" and "Fraction" are stacked; data in the column "log profile" originates from logarithmized "LFQ intensity" 
                                and "Ratio H/L", respectively; additionally the columns "MS/MS count" and "Ratio H/L count|Ratio H/L variability [%]" are stored 
                                as single level indices
                df_01_stacked: dataframe, in which "MAP" and "Fraction" are stacked; data in the column "LFQ intensity" is 0-1 normalized and renamed to "normalized 
                               profile"; the columns "normalized profile"" and "MS/MS count" are stored as single level indices; plotting is possible now

        Returns:
            self:
                df_pca: df, PCA was performed, while keeping the information of the Maps 
                            columns: "PC1", "PC2", "PC3"
                            index: "Protein IDs", "Majority protein IDs", "Protein names", "Gene names", "Q-value", "Score", "id", "Map" "Compartment"
                df_pca_combined: df, PCA was performed across the Maps 
                            columns: "PC1", "PC2", "PC3"
                            index: "Protein IDs", "Majority protein IDs", "Protein names", "Gene names", "Q-value", "Score", "id", "Compartment"
                df_pca_all_marker_cluster_maps: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3", filtered for marker genes, that are consistent 
                                                throughout all maps / coverage filtering.
        """

        markerproteins = self.markerproteins

        if self.acquisition == "SILAC - MQ":
            df_01orlog_fracunstacked = self.df_log_stacked["log profile"].unstack("Fraction").dropna()
            df_01orlog_MapFracUnstacked = self.df_log_stacked["log profile"].unstack(["Fraction", "Map"]).dropna()  
            
        elif self.acquisition.startswith("LFQ") or self.acquisition == "Custom":
            df_01orlog_fracunstacked = self.df_01_stacked["normalized profile"].unstack("Fraction").dropna()
            df_01orlog_MapFracUnstacked = self.df_01_stacked["normalized profile"].unstack(["Fraction", "Map"]).dropna()
            
            
        pca = PCA(n_components=3)

        # df_pca: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3"
        df_pca = pd.DataFrame(pca.fit_transform(df_01orlog_fracunstacked.apply(zscore, axis=0).replace(np.nan, 0)))
        df_pca.columns = ["PC1", "PC2", "PC3"]
        df_pca.index = df_01orlog_fracunstacked.index
        self.df_pca = df_pca.sort_index(level=["Protein IDs", "Compartment"])
        
        # df_pca: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3"
        df_pca_combined = pd.DataFrame(pca.fit_transform(df_01orlog_MapFracUnstacked.apply(zscore, axis=0).replace(np.nan, 0)))
        df_pca_combined.columns = ["PC1", "PC2", "PC3"]
        df_pca_combined.index = df_01orlog_MapFracUnstacked.index
        self.df_pca_combined = df_pca_combined.sort_index(level=["Protein IDs", "Compartment"])
        
        map_names = self.map_names
        df_pca_all_marker_cluster_maps = pd.DataFrame()
        df_pca_filtered = df_pca.unstack("Map").dropna()
        for clusters in markerproteins:
            for marker in markerproteins[clusters]:
                try:
                    plot_try_pca = df_pca_filtered.xs(marker, level="Protein IDs", drop_level=False)
                except KeyError:
                    continue
                df_pca_all_marker_cluster_maps = df_pca_all_marker_cluster_maps.append(
                    plot_try_pca)
        if len(df_pca_all_marker_cluster_maps) == 0:
            df_pca_all_marker_cluster_maps = df_pca_filtered.stack("Map")
        else:
            df_pca_all_marker_cluster_maps = df_pca_all_marker_cluster_maps.stack("Map")
        self.df_pca_all_marker_cluster_maps = df_pca_all_marker_cluster_maps.sort_index(level=["Protein IDs", "Compartment"])

        
    def plot_global_pca(self, map_of_interest="Map1", cluster_of_interest="Proteasome", x_PCA="PC1", y_PCA="PC3", collapse_maps=False):
        """"
        PCA plot will be generated

        Args:
            self:
                df_organellarMarkerSet: df, columns: "Protein ID", "Compartment", no index
                df_pca: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3",
                    index: "Gene names", "Protein IDs", "C-Score", "Q-value", "Map", "Compartment", 

        Returns:
            pca_figure: global PCA plot
        """
        
        if collapse_maps == False:
            df_global_pca = self.df_pca.unstack("Map").swaplevel(0,1, axis=1)[map_of_interest].reset_index()
        else:
            df_global_pca = self.df_pca_combined.reset_index()
            
        for i in self.markerproteins[cluster_of_interest]:
            df_global_pca.loc[df_global_pca["Protein IDs"] == i, "Compartment"] = "Selection"

        compartments = self.df_organellarMarkerSet["Compartment"].unique()
        compartment_color = dict(zip(compartments, self.css_color))
        compartment_color["Selection"] = "black"
        compartment_color["undefined"] = "lightgrey"
        
        fig_global_pca = px.scatter(data_frame=df_global_pca,
                                    x=x_PCA,
                                    y=y_PCA,
                                    color="Compartment",
                                    color_discrete_map=compartment_color,
                                    title= "Protein subcellular localization by PCA for {}".format(map_of_interest) 
                                        if collapse_maps == False else "Protein subcellular localization by PCA of combined maps", 
                                    hover_data=["Protein IDs", "Gene names", "Compartment"],
                                    template="simple_white",
                                    opacity=0.9
                                    )
        return fig_global_pca                         
            
            
    def plot_cluster_pca(self, cluster_of_interest="Proteasome"):
        """
        PCA plot will be generated

        Args:
            self:
                markerproteins: dictionary, key: cluster name, value: gene names (e.g. {"Proteasome" : ["PSMA1", "PSMA2",...], ...}
                df_pca_all_marker_cluster_maps: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3", filtered for marker genes, that are 
                                                consistent throughout all maps / coverage filtering.

        Returns:
            pca_figure: PCA plot, for one protein cluster all maps are plotted
        """

        df_pca_all_marker_cluster_maps = self.df_pca_all_marker_cluster_maps
        map_names = self.map_names
        markerproteins = self.markerproteins
        
        try:
            for maps in map_names:
                df_setofproteins_PCA = pd.DataFrame()
                for marker in markerproteins[cluster_of_interest]:
                    try:
                        plot_try_pca = df_pca_all_marker_cluster_maps.xs((marker, maps), level=["Protein IDs", "Map"],
                                                                         drop_level=False)
                    except KeyError:
                        continue
                    df_setofproteins_PCA = df_setofproteins_PCA.append(plot_try_pca)
    
                df_setofproteins_PCA.reset_index(inplace=True)
                if maps == map_names[0]:
                    pca_figure = go.Figure(
                        data=[go.Scatter3d(x=df_setofproteins_PCA.PC1, 
                                           y=df_setofproteins_PCA.PC2, 
                                           z=df_setofproteins_PCA.PC3,
                                           hovertext=df_setofproteins_PCA["Gene names"], 
                                           mode="markers", 
                                           name=maps
                                           )])
                else:
                    pca_figure.add_trace(go.Scatter3d(x=df_setofproteins_PCA.PC1, 
                                                      y=df_setofproteins_PCA.PC2, 
                                                      z=df_setofproteins_PCA.PC3,
                                                      hovertext=df_setofproteins_PCA["Gene names"], 
                                                      mode="markers", 
                                                      name=maps
                                                     ))

            pca_figure.update_layout(autosize=False, width=500, height=500,
                                  title="PCA plot for <br>the protein cluster: {}".format(cluster_of_interest), 
                                  template="simple_white")
            return pca_figure
        
        except:
            return "This protein cluster was not quantified"
    
    
    def calc_biological_precision(self):
        """
        This function calculates the biological precision of all quantified protein clusters. It provides access to the data slice for all marker proteins, the distance profiles and the aggregated distances. It repeatedly applies the methods get_marker_proteins_unfiltered and calc_cluster_distances.
        
        TODO: integrate optional arguments for calc_cluster_distances: complex_profile, distance_measure.
        TODO: replace compatibiliy attributes with function return values and adjust attribute usage in downstream plotting functions.
        
        Args:
            self attributes:
                markerproteins: dict, contains marker protein assignments
                df_01_stacked: df, contains 0-1 nromalized data, required for execution of get_marker_proteins_unfiltered
        Returns:
            df_alldistances_individual_mapfracunstacked: df, distance profiles, fully unstacked
            df_alldistances_aggregated_mapunstacked: df, profile distances (manhattan distance by default), fully unstacked
            df_allclusters_01_unfiltered_mapfracunstacked: df, collected marker protein data
            self attributes:
                df_distance_noindex: compatibility version of df_alldistances_aggregated_mapunstacked
                df_allclusters_01_unfiltered_mapfracunstacked
                df_allclusters_clusterdist_fracunstacked_unfiltered: compatibility version of df_allclusters_01_unfiltered_mapfracunstacked (only used by quantificaiton_overview)
                df_allclusters_clusterdist_fracunstacked: compatibility version of df_alldistances_individual_mapfracunstacked
                genenames_sortedout_list = list of gene names with incomplete coverage
                analysis_summary_dict entries:
                    "Manhattan distances" = df_distance_noindex
                    "Distances to the median profile": df_allclusters_clusterdist_fracunstacked, sorted and melted
        """
        df_alldistances_individual_mapfracunstacked = pd.DataFrame()
        df_alldistances_aggregated_mapunstacked = pd.DataFrame()
        df_allclusters_01_unfiltered_mapfracunstacked = pd.DataFrame()
        
        for cluster in self.markerproteins.keys():
            # collect data irrespective of coverage
            df_cluster_unfiltered = self.get_marker_proteins_unfiltered(cluster)
            df_allclusters_01_unfiltered_mapfracunstacked = df_allclusters_01_unfiltered_mapfracunstacked.append(df_cluster_unfiltered)
            
            # filter for coverage and calculate distances
            df_cluster = df_cluster_unfiltered.dropna()
            if len(df_cluster) == 0:
                continue
            df_distances_aggregated, df_distances_individual = self.calc_cluster_distances(df_cluster)
            df_alldistances_individual_mapfracunstacked = df_alldistances_individual_mapfracunstacked.append(df_distances_individual)
            df_alldistances_aggregated_mapunstacked = df_alldistances_aggregated_mapunstacked.append(df_distances_aggregated)
        if len(df_alldistances_individual_mapfracunstacked) == 0:
            self.df_distance_noindex = pd.DataFrame(columns = ["Protein IDs", "Gene names", "Map", "Cluster", "distance"])
            self.df_allclusters_01_unfiltered_mapfracunstacked = pd.DataFrame(columns = ["Protein IDs", "Gene names", "Map", "Cluster", "distance"])
            self.df_allclusters_clusterdist_fracunstacked_unfiltered = pd.DataFrame(columns = ["Fraction"])
            self.df_allclusters_clusterdist_fracunstacked = pd.DataFrame(columns = ["Fraction"])
            self.genenames_sortedout_list = "No clusters found"
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        else:
            df_alldistances_aggregated_mapunstacked.columns.name = "Map"
            ## Get compatibility with plotting functions, by mimicking assignment of old functions:
            # old output of distance_calculation
            self.df_distance_noindex = df_alldistances_aggregated_mapunstacked.stack("Map").reset_index().rename({0: "distance"}, axis=1) 
            self.analysis_summary_dict["Manhattan distances"] = self.df_distance_noindex.to_json()
            # old output of multiple_iterations
            # self.df_allclusters_clusterdist_fracunstacked_unfiltered --> this won't exist anymore, replaced by:
            self.df_allclusters_01_unfiltered_mapfracunstacked = df_allclusters_01_unfiltered_mapfracunstacked
            # kept for testing of quantification table:
            self.df_allclusters_clusterdist_fracunstacked_unfiltered = df_allclusters_01_unfiltered_mapfracunstacked.stack("Map")
            # same as before, but now already abs
            self.df_allclusters_clusterdist_fracunstacked = df_alldistances_individual_mapfracunstacked.stack("Map")
            df_dist_to_median = self.df_allclusters_clusterdist_fracunstacked.stack("Fraction")
            df_dist_to_median.name = "distance"
            df_dist_to_median = df_dist_to_median.reindex(index=natsort.natsorted(df_dist_to_median.index))
            self.analysis_summary_dict["Distances to the median profile"] = df_dist_to_median.reset_index().to_json()
            self.genenames_sortedout_list = [el for el in df_allclusters_01_unfiltered_mapfracunstacked.index.get_level_values("Gene names")
                                            if el not in df_alldistances_individual_mapfracunstacked.index.get_level_values("Gene names")]
                                            
            return df_alldistances_individual_mapfracunstacked, df_alldistances_aggregated_mapunstacked, df_allclusters_01_unfiltered_mapfracunstacked
    
    
    def get_marker_proteins_unfiltered(self, cluster):
        """
        This funciton retrieves the 0-1 normalized data for any given protein cluster, unfiltered for coverage.
        
        Args:
            cluster: str, cluster name, should be one of self.markerproteins.keys()
            self attributes:
                df_01_stacked: df, contains the fully stacked 0-1 normalized data
                markerproteins: dict, contains marker protein assignments
        
        Returns:
            df_cluster_unfiltered: df, unfiltered data for the selected cluster, maps and fractions are unstacked.
            self attribtues:
                None
        """
        df_in = self.df_01_stacked["normalized profile"].unstack("Fraction")
        markers = self.markerproteins[cluster]
        
        # retrieve marker proteins
        df_cluster_unfiltered = pd.DataFrame()
        for marker in markers:
            try:
                df_p = df_in.xs(marker, level="Protein IDs", axis=0, drop_level=False)
            except:
                continue
            df_cluster_unfiltered = df_cluster_unfiltered.append(df_p)
        if len(df_cluster_unfiltered) == 0:
            return df_cluster_unfiltered
        
        # Unstack maps and add Cluster to index
        df_cluster_unfiltered = df_cluster_unfiltered.unstack("Map")
        df_cluster_unfiltered.set_index(pd.Index(np.repeat(cluster, len(df_cluster_unfiltered)), name="Cluster"), append=True, inplace=True)
        
        return df_cluster_unfiltered
    
    
    def calc_cluster_distances(self, df_cluster, complex_profile=np.median, distance_measure="manhattan"):
        """
        Calculates the absolute differences in each fraction and the profile distances relative to the center of a cluster.
        Per default this is the manhattan distance to the median profile.
        
        Args:
            df_cluster: df, 0-1 normalized profiles of cluster members, should already be filtered for full coverage and be in full wide format.
            complex_profile: fun, function provided to apply for calculating the reference profile, default: np.median.
            distance_measure: str, selected distance measure to calculate. Currently only 'manhattan' is supported, everything else raises a ValueError.
            self attributes:
                None
        Returns:
            df_distances_aggregated: df, proteins x maps, if stacked distance column is currently named 0 but contains manhattan distances.
            df_distances_individual: df, same shape as df_cluster, but now with absolute differences to the reference.
            self attribtues:
                None
        """
        df_distances_aggregated = pd.DataFrame()
        
        ref_profile = pd.DataFrame(df_cluster.apply(complex_profile, axis=0, result_type="expand")).T
        df_distances_individual = df_cluster.apply(lambda x: np.abs(x-ref_profile.iloc[0,:]), axis=1)
            
        # loop over maps
        maps = set(df_cluster.columns.get_level_values("Map"))
        for m in maps:
            if distance_measure == "manhattan":
                d_m = pw.manhattan_distances(df_cluster.xs(m, level="Map", axis=1), ref_profile.xs(m, level="Map", axis=1))
            else:
                raise ValueError(distance_measure)
            d_m = pd.DataFrame(d_m, columns=[m], index=df_cluster.index)
            df_distances_aggregated = pd.concat([df_distances_aggregated, d_m], axis=1)
        
        df_distances_aggregated.columns.set_names(names="Map", inplace=True)
        return df_distances_aggregated, df_distances_individual
    
    
    def profiles_plot(self, map_of_interest="Map1", cluster_of_interest="Proteasome"):
        """
        The function allows the plotting of filtered and normalized spatial proteomic data using plotly.express.
        The median profile is also calculated based on the overlapping proteins. Profiles of proteins that are not quantified in all maps are dashed.

        Args:
            map_of_interest: str, must be in self.map_names
            cluster_of_interest: str, must be in self.markerproteins.keys()
            self attribtues:
                df_allclusters_01_unfiltered_mapfracunstacked: df, contains 0-1 normalized profiles for all markerproteins detected in any map

        Returns:
            abundance_profiles_and_median_figure: plotly line plot, displaying the relative abundance profiles.
        """
        
        try:
            df_setofproteins = self.df_allclusters_01_unfiltered_mapfracunstacked.xs(cluster_of_interest, level="Cluster", axis=0)
            df_setofproteins_median = df_setofproteins.dropna().xs(map_of_interest, level="Map", axis=1).median(axis=0)
    
            # fractions get sorted
            df_setofproteins = df_setofproteins.xs(map_of_interest, level="Map", axis=1).stack("Fraction")
            df_setofproteins = df_setofproteins.reindex(index=natsort.natsorted(df_setofproteins.index))
            df_setofproteins.name = "normalized profile"
    
            # make it available for plotting
            df_setofproteins = df_setofproteins.reindex(index=natsort.natsorted(df_setofproteins.index))
            df_setofproteins = df_setofproteins.reset_index()
            abundance_profiles_figure = px.line(df_setofproteins, 
                                                x="Fraction", 
                                                y="normalized profile",
                                                color="Gene names",
                                                line_group="Sequence" if "Sequence" in df_setofproteins.columns else "Gene names",
                                                template="simple_white",
                                                hover_data=["Protein IDs", "Gene names"],
                                                title="Relative abundance profile for {} of <br>the protein cluster: {}".format(map_of_interest, cluster_of_interest)
                                               )
    
            df_setofproteins_median.name = "normalized profile"
    
            #fractions get sorted
            df_setofproteins_median = df_setofproteins_median.reindex(index=natsort.natsorted(df_setofproteins_median.index))
    
            # make it available for plotting
            df_setofproteins_median = df_setofproteins_median.reset_index()
            df_setofproteins_median.insert(0, "Gene names", np.repeat("Median profile", len(df_setofproteins_median)))
    
            abundance_profiles_and_median_figure = abundance_profiles_figure.add_scatter(x=df_setofproteins_median["Fraction"],
                                                                                         y=df_setofproteins_median["normalized profile"],
                                                                                         name="Median profile"
                                                                                        )
            # dash lines for proteins that have insufficient coverage across maps
            abundance_profiles_and_median_figure.for_each_trace(lambda x: x.update(line={"dash":"dash"}),
                                                                selector=lambda x: x.name in self.genenames_sortedout_list)
    
            return abundance_profiles_and_median_figure
        
        except:
            return "This protein cluster was not quantified"


    def quantification_overview(self, cluster_of_interest="Proteasome"):
        """
        
        Args:
            self.df_allclusters_clusterdist_fracunstacked_unfiltered 
                columns: 01K, 03K, 06K, 12K, 24K, 80K
                index: Gene names, Protein IDs, C-Score, Q-value, Map, Compartment, Cluster
            
        Returns:
            df
        """
        
        df_quantification_overview = self.df_allclusters_clusterdist_fracunstacked_unfiltered.xs(cluster_of_interest, level="Cluster", axis=0)\
                                                                                             [self.fractions[0]].unstack("Map")
        if "Sequence" in df_quantification_overview.index.names:
            df_quantification_overview = df_quantification_overview.droplevel([i for i in df_quantification_overview.index.names if not i in ["Sequence","Gene names"]])
        else:
            df_quantification_overview = df_quantification_overview.droplevel([i for i in df_quantification_overview.index.names if not i=="Gene names"])
        df_quantification_overview = df_quantification_overview.notnull().replace({True: "x", False: "-"})
        
        return df_quantification_overview

    
    def distance_boxplot(self, cluster_of_interest="Proteasome"):
        """
        A box plot for 1 desired cluster, and across all maps is generated displaying the distribution of the e.g.
        Manhattan distance.

        Args:
            self:
            df_distance_noindex: stored as attribute (self.df_distance_noindex),index is reset.
            It contains the column name "distance", in which the e.g. Manhattan distances for each individual protein
            of the specified clusters (see self.markerproteins) are stored
            map_names: individual map names are stored as an index

        Returns:
            distance_boxplot_figure: boxplot. Along the x-axis the maps, along the y-axis the distances are shown
        """

        map_names = self.map_names

        df_distance_noindex = self.df_distance_noindex

        # "Gene names", "Map", "Cluster" and transferred into the index
        df_distance_map_cluster_gene_in_index = df_distance_noindex.set_index(["Gene names", "Map", "Cluster"])
        if "Sequence" in df_distance_map_cluster_gene_in_index.columns:
            df_distance_map_cluster_gene_in_index.set_index("Sequence", append=True, inplace=True)

        df_cluster_xmaps_distance_with_index = pd.DataFrame()
        
        try:
        # for each individual map and a defined cluster data will be extracted from the dataframe
        # "df_distance_map_cluster_gene_in_index" and appended to the new dataframe df_cluster_xmaps_distance_with_index
            for maps in map_names:
                plot_try = df_distance_map_cluster_gene_in_index.xs((cluster_of_interest, maps),
                                                                    level=["Cluster", "Map"], drop_level=False)
                df_cluster_xmaps_distance_with_index = df_cluster_xmaps_distance_with_index.append(plot_try)
                
            df_cluster_xmaps_distance_with_index["Combined Maps"] = "Combined Maps"
            
            #number of proteins within one cluster
            self.proteins_quantified_across_all_maps = df_cluster_xmaps_distance_with_index.unstack("Map").shape[0]
        
            # index will be reset, required by px.box
            df_cluster_xmaps_distance = df_cluster_xmaps_distance_with_index.reset_index()
    
            distance_boxplot_figure = go.Figure()
            distance_boxplot_figure.add_trace(go.Box(
                x=df_cluster_xmaps_distance["Map"],
                y=df_cluster_xmaps_distance["distance"],
                boxpoints="all",
                whiskerwidth=0.2,
                marker_size=2,
                hovertext=df_cluster_xmaps_distance["Gene names"]
            ))
            
            distance_boxplot_figure.add_trace(go.Box(
                x=df_cluster_xmaps_distance["Combined Maps"],
                y=df_cluster_xmaps_distance["distance"],
                boxpoints="all",
                whiskerwidth=0.2,
                marker_size=2,
                hovertext=df_cluster_xmaps_distance["Gene names"]
            ))  
    
            distance_boxplot_figure.update_layout(
                title="Manhattan distance distribution for <br>the protein cluster: {}".format(cluster_of_interest),
                autosize=False,
                showlegend=False,
                width=500,
                height=500,
    
                # black box around the graph
                xaxis=go.layout.XAxis(linecolor="black",
                                      linewidth=1,
                                      title="Map",
                                      mirror=True),
    
                yaxis=go.layout.YAxis(linecolor="black",
                                      linewidth=1,
                                      title="distance",
                                      mirror=True),
                template="simple_white"
            )
    
            return distance_boxplot_figure
    
        except:
            self.cache_cluster_quantified = False
            
    
    def distance_to_median_boxplot(self, cluster_of_interest="Proteasome"):
        """
        A box plot for 1 desired cluster, across all maps and fractions is generated displaying the
        distribution of the distance to the median. For each fraction, one box plot will be displayed.

        Args:
            self:
            df_allclusters_clusterdist_fracunstacked, dataframe with single level column, stored as attribute
            (self.allclusters_clusterdist_fracunstacked), in which "Fraction" is unstacked. It contains only the
            normalized data of individual protein clusters substracted by the median of the respective protein cluster
            for each fraction.
            map_names: individual map names are stored as an index

        Returns:
            distance_to_median_boxplot_figure: Box plot. Along the x-axis, the maps are shown, along the y-axis
            the distances is plotted
        """
        
        df_boxplot_manymaps = pd.DataFrame()

        try:
        # for each individual map and a defined cluster data will be extracted from the dataframe
        # "df_allclusters_clusterdist_fracunstacked" and appended to the new dataframe df_boxplot_manymaps
            for maps in self.map_names:
                plot_try = self.df_allclusters_clusterdist_fracunstacked.xs((cluster_of_interest, maps), level=["Cluster", "Map"], drop_level=False)
                df_boxplot_manymaps = df_boxplot_manymaps.append(plot_try)
            
            self.df_boxplot_manymaps = df_boxplot_manymaps
        
            # index will be reset, required by px.violin
            df_boxplot_manymaps = abs(df_boxplot_manymaps.stack("Fraction"))
    
            df_boxplot_manymaps.name = "distance"
    
            df_boxplot_manymaps = df_boxplot_manymaps.reindex(index=natsort.natsorted(df_boxplot_manymaps.index))
    
            df_boxplot_manymaps = df_boxplot_manymaps.reset_index()
            
            # box plot will be generated, every fraction will be displayed in a single plot
            distance_to_median_boxplot_figure = px.box(df_boxplot_manymaps, 
                                                       x="Map", 
                                                       y="distance", 
                                                       facet_col="Fraction",
                                                       facet_col_wrap=2,
                                                       boxmode="overlay", height=900, width=700, points="all",
                                                       hover_name="Gene names",
                                                       template="simple_white",
                                                       title="Distribution of the distance to the median for <br>the protein cluster: {}".format(cluster_of_interest))
    
            return distance_to_median_boxplot_figure
            
        except:
            return "This protein cluster was not quantified" 
    
    
    def results_overview_table(self):
        """
        Dataframe will be created, that provides information about "range", "mean" and "standardeviation",
        given as the column names, based on the data given in df_distance_noindex

        Args:
            self:
                df_distance_noindex: stored as attribute (self.df_distance_noindex),index is reset. It contains the column name "distance", 
                                     in which the e.g. Manhattan distances for each individual protein of the specified clusters (see self.markerproteins) 
                                     are stored
                markerproteins: dictionary, key: cluster name, value: gene names (e.g. {"Proteasome" : ["PSMA1", "PSMA2",...], ...}
                
        """
        
        df_distance_noindex = self.df_distance_noindex
        df_distance_map_cluster_gene_in_index = df_distance_noindex.set_index(["Gene names", "Map", "Cluster"])
        map_names = self.map_names

        df_overview = pd.DataFrame()

        for clusters in self.markerproteins:
            #if a certain cluster is not available in the dataset at all
            try:
                for maps in map_names:
                    df_dist_map_cluster = df_distance_map_cluster_gene_in_index.xs((clusters, maps), level=["Cluster", "Map"], drop_level=False)
                    statistic_table = {"range": (df_dist_map_cluster["distance"].max(axis=0)) - (df_dist_map_cluster["distance"].min(axis=0)),
                                       "median": df_dist_map_cluster["distance"].median(axis=0),
                                       "standardeviation": df_dist_map_cluster["distance"].std(axis=0),
                                       "Cluster": clusters,
                                       "Map": maps
                                      }
                    statistic_series = pd.Series(data=statistic_table)
                    df_statistic_table_individual_cluster = pd.DataFrame(statistic_series).T
                    df_overview = df_overview.append(df_statistic_table_individual_cluster)
                
                df_dist_cluster = df_distance_map_cluster_gene_in_index.xs(clusters, level="Cluster")    
                statistic_table_combined = {
                    "range": (df_dist_cluster["distance"].max(axis=0)) - (df_dist_cluster["distance"].min(axis=0)),
                    "median": df_dist_cluster["distance"].median(axis=0),
                    "standardeviation": df_dist_cluster["distance"].std(axis=0),
                    "Cluster": clusters,
                    "Map": "combined maps"
                }
                statistic_series_combined = pd.Series(data=statistic_table_combined)
                df_statistic_table_individual_cluster = pd.DataFrame(statistic_series_combined).T
                df_overview = df_overview.append(df_statistic_table_individual_cluster)
                
            except:
                continue
        try:
            df_overview.set_index(["Cluster", "Map"], inplace=True)
            df_overview.sort_index(axis=0, level=0, inplace=True)
        except:
            df_overview = pd.DataFrame()

        self.analysis_summary_dict["Overview table"] = df_overview.reset_index().to_json()
        self.analysed_datasets_dict[self.expname] = self.analysis_summary_dict.copy() 
#self.analysis_summary_dict.clear()
        
        return df_overview


    def reframe_df_01ORlog_for_Perseus(self, df_01ORlog):
        """"
        To be available for Perseus df_01_stacked needs to be reframed.

        Args:
            df_01ORlog:
            df_distance_noindex: stored as attribute (self.df_distance_noindex),index is reset.
            It contains the column name "distance", in which the e.g. Manhattan distances for each individual protein
            of the specified clusters (see self.markerproteins) are stored
            map_names: individual map names are stored as an index

        Returns:
            df_01ORlog_svm: 
                LFQ:
                columns: "MS/MS count_Map1_01K", "normalized profile_Map1_01K"
                index: "Gene names", "Protein IDs", "C-Score", "Q-value", "Compartment"
                
                SILAC:
                columns: e.g. "Ratio H/L count_MAP2_80K", "Ratio H/L variability [%]_MAP1_03K", "normalized profile_MAP5_03K" 
                index: "Q-value", "Score", "Protein IDs", "Majority protein IDs", "Protein names", "Gene names", "id", "Compartment"
                
        """
        
        df_01ORlog_svm = df_01ORlog.copy()
        
        #df_01_filtered_combined = df_01_filtered_combined.stack(["Experiment", "Map"]).swaplevel(0,1, axis=0).dropna(axis=1)
        index_ExpMap = df_01ORlog_svm.index.get_level_values("Map")+"_"+df_01ORlog_svm.index.get_level_values("Fraction")
        index_ExpMap.name = "Map_Frac"
        df_01ORlog_svm.set_index(index_ExpMap, append=True, inplace=True)      
        
        df_01ORlog_svm.index = df_01ORlog_svm.index.droplevel(["Map", "Fraction"])
        df_01ORlog_svm = df_01ORlog_svm.unstack("Map_Frac")
        #df_01ORlog_svm = df_01ORlog_svm.dropna(axis=0, subset=df_01ORlog_svm.loc[[], ["normalized profile"]].columns)
        df_01ORlog_svm.columns = ["_".join(col) for col in df_01ORlog_svm.columns.values]
        df_01ORlog_svm.rename(index={"undefined" : np.nan}, level="Compartment", inplace=True)
        
        return df_01ORlog_svm


class SpatialDataSetComparison:
    
        
    analysed_datasets_dict = {}
    css_color = SpatialDataSet.css_color
    #cache_stored_SVM = True


    def __init__(self, ref_exp="Exp2", **kwargs): #clusters_for_ranking=["Proteasome", "Lysosome"]
        
        #self.clusters_for_ranking = clusters_for_ranking
        self.ref_exp = ref_exp
        self.json_dict = {}
        #self.fractions, self.map_names = [], []  #self.df_01_stacked, self.df_log_stacked = pd.DataFrame(), pd.DataFrame()
        #collapse_maps,collapse_cluster,  cluster_of_interest_comparison, multi_choice, multi_choice_venn, x_PCA_comp, y_PCA_comp
        
        
        #if "organism" not in kwargs.keys():
        #    self.markerproteins = self.markerproteins_set["Human - Swissprot"]
        #else:
        #    assert kwargs["organism"] in self.markerproteins_set.keys()
        #    self.markerproteins = self.markerproteins_set[kwargs["organism"]]
        #    del kwargs["organism"]
        
        self.exp_names, self.exp_map_names = [], []
        
        self.df_01_filtered_combined, self.df_distance_comp = pd.DataFrame(), pd.DataFrame()
        self.df_quantity_pr_pg_combined = pd.DataFrame()
        self.svm_results = dict()
        
        self.parameters = dict()
        

    def read_jsonFile(self): #, content=None
        """
        Read-out of the JSON-file and currently analysed dataset, stored in "analysed_datasets_dict". It wil create df_distances_combined ("Gene
        names", "Cluster" are stacked; "Map" and Experiment names (are not stored in an additional level name) are unstacked. Layout will be 
        adjusted for distance-plotting.
    
        Args: 
            self.json_dict: contains the dictionary stored in AnalysedDatasets.json
            
            {"Experiment name" : {
                "changes in shape after filtering" : {
                
                  ##SILAC##
                    "Original size" : tuple,
                    "Shape after categorical filtering" : tuple,
                    "Shape after Ratio H/L count (>= 3)/var (count>=2, var<30) filtering" : tuple,
                    "Shape after filtering for complete profiles" : tuple,
                    
                  ##LFQ/spectronaut##
                    "Original size" : tuple,
                    "Shape after MS/MS value filtering" : tuple,
                    "Shape after consecutive value filtering" : tuple,
                },
                   
                "quantity: profiles/protein groups" : df - number of protein groups | number of profiles | data completeness of profiles
                "Analysis parameters" : {
                    "acquisition" : str,
                    "filename" : str,
                    
                  ##SILAC##
                    "Ratio H/L count 1 (>= X)" : int,
                    "Ratio H/L count 2 (>=Y, var<Z)" : int,
                    "Ratio variability (<Z, count>=Y)" : int,
                    
                  ##LFQ/spectronaut##
                    "consecutive data points" : int,
                    "summed MS/MS counts" : int,
                },

                "0/1 normalized data - mean" : df - mean of all datapoints,
                "0/1 normalized data" : df - individual cluster,
                "Distances to the median profile" : df - individual cluster,
                "Manhattan distances" : df - individual cluster,
                "Overview table" : df - individual cluster,

               ##if user perform the Misclassification Analysis befor downloading the dictionary AnalysedDatasets.json##
                {"Misclassification Analysis": {
                    "True: ER" : {
                        "Recall": int,
                        "FDR": int,
                        "Precision": int,
                        "F1": int
                        }
                    "True: NPC" : {...}
                    ...
                    "Summary": {
                        "Total - Recall": int, 
                        "Membrane - Recall" : int, 
                        "Av per organelle - Recall": int,
                        "Median per organelle - Recall" : int,
                        "Av precision organelles" : int,
                        "Av F1 organelles" : int,
                        "Av F1 all clusters" :  int,
                        }
                    }   
                }
            }
            
        
        Returns: 
            self:
                df_01_filtered_combined: df, "Fraction" is unstacked; "Experiment", "Gene names", "Map", "Exp_Map" are stacked 
                df_distance_comp: df, no index, column names: "Gene names", "Cluster", "Protein IDs", "Compartment", "Experiment", "Map", "Exp_Map", "distance"
                            "distance": Manhattan distances for each individual protein of the specified clusters (see self.markerproteins) are stored
                df_quantity_pr_pg_combined: df, no index, column names: "filtering", "type", "number of protein groups", "number of profiles", 
                                            "data completeness of profiles", "Experiment"
                exp_map_names: list of unique Exp_Map - fusions e.g. LFQ_Map1
                exp_names: list of unique Experiment names - e.g. LFQ
        """
        
        
        json_dict = self.json_dict
    
        self.analysis_parameters_total = {}
        self.exp_names = list(json_dict.keys())
        
        df_01_combined = pd.DataFrame()
        for exp_name in json_dict.keys():
            for data_type in json_dict[exp_name].keys():
                if data_type == "0/1 normalized data":
                    df_01_toadd = pd.read_json(json_dict[exp_name][data_type])
                    df_01_toadd.insert(0,"Experiment",np.repeat(exp_name, len(df_01_toadd)))
                    if len(df_01_combined) == 0:
                        df_01_combined = df_01_toadd.copy()
                    else:
                        df_01_combined = pd.concat([df_01_combined,df_01_toadd], sort=False, axis=0)
                        
                elif data_type == "quantity: profiles/protein groups" and exp_name == list(json_dict.keys())[0]:
                    df_quantity_pr_pg_combined = pd.read_json(json_dict[exp_name][data_type])
                    df_quantity_pr_pg_combined["Experiment"] = exp_name
                    
                elif data_type == "quantity: profiles/protein groups" and exp_name != list(json_dict.keys())[0]:
                    df_quantity_pr_pg_toadd = pd.read_json(json_dict[exp_name][data_type])
                    df_quantity_pr_pg_toadd["Experiment"] = exp_name
                    df_quantity_pr_pg_combined = pd.concat([df_quantity_pr_pg_combined, df_quantity_pr_pg_toadd])
                
                elif data_type == "Analysis parameters":
                    self.analysis_parameters_total[exp_name] = json_dict[exp_name][data_type]
                
                elif data_type == "SVM results":
                    for res in self.json_dict[exp_name]["SVM results"].keys():
                        self.add_svm_result(exp_name,
                                             pd.read_json(self.json_dict[exp_name]["SVM results"][res]["misclassification"]),
                                             name=res,
                                             prediction=pd.read_json(self.json_dict[exp_name]["SVM results"][res]["prediction"]),
                                             comment=self.json_dict[exp_name]["SVM results"][res]["comment"],
                                             overwrite=True # supposed to always exceed upload from old format
                                             )
                
                elif data_type == "Misclassification Matrix":
                    try:
                        self.add_svm_result(exp_name, pd.read_json(self.json_dict[exp_name]["Misclassification Matrix"]),
                                            comment="read from old json file version", overwrite=False) 
                    except:
                        # should only happen if this has been loaded before and both SVM results and Misclassification Matrix are present
                        pass
                    
                #try:
                #    for paramters in json_dict[exp_name][data_type].keys():
                #        if paramters=="acquisition":
                #            acquisition_loaded.append(json_dict[exp_name][data_type][paramters])
                #        #elif parameters=="Non valid profiles":
                #except:
                #    continue
                #
        ### New code for alignment:
        df_01_combined = df_01_combined.set_index(["Experiment", "Gene names", "Compartment"])
        if "Original Protein IDs" in df_01_combined.columns and "Protein IDs" in df_01_combined.columns:
            df_01_combined.set_index(["Original Protein IDs"], append=True, inplace=True)
            df_01_combined.drop("Protein IDs", axis=1, inplace=True)
        elif "Original Protein IDs" in df_01_combined.columns:
            df_01_combined.set_index("Original Protein IDs", append=True, inplace=True)
        elif "Protein IDs" in df_01_combined.columns:
            df_01_combined.set_index(pd.Index(df_01_combined["Protein IDs"], name="Original Protein IDs"), append=True, inplace=True)
            df_01_combined.drop("Protein IDs", axis=1, inplace=True)
        else:
            raise KeyError("No Protein IDs were found while loading the json file")
            
        df_01_combined.columns = pd.MultiIndex.from_tuples([el.split("?")[1::] for el in df_01_combined.columns], names=["Map", "Fraction"])
        df_01_combined = df_01_combined.stack("Map").dropna().unstack("Map")
        df_01_filtered_combined, id_alignment = align_datasets(df_01_combined)
        self.id_alignment = id_alignment
        index_ExpMap = df_01_filtered_combined.index.get_level_values("Experiment")+"_"+df_01_filtered_combined.index.get_level_values("Map")
        index_ExpMap.name = "Exp_Map"
        df_01_filtered_combined.set_index(index_ExpMap, append=True, inplace=True)
        
        self.exp_map_names = list(index_ExpMap.unique())
        
        self.df_01_filtered_combined = df_01_filtered_combined
        
        self.df_quantity_pr_pg_combined = df_quantity_pr_pg_combined
        
        try:
            organism = json_dict[list(json_dict.keys())[0]]["Analysis parameters"]['organism']
        except:
            organism = "Homo sapiens - Uniprot"
        
        marker_table = pd.read_csv(pkg_resources.resource_stream(__name__, 'annotations/complexes/{}.csv'.format(organism)))
        self.markerproteins = {k: v.replace(" ", "").split(",") for k,v in zip(marker_table["Cluster"], marker_table["Members - Protein IDs"])}
        
        self.clusters_for_ranking = self.markerproteins.keys()        
           

    def add_svm_result(self, experiment, misclassification, name="default", prediction=pd.DataFrame(), comment="", overwrite=True):
        """
        Add SVM output (either only misclassification matrix or whole prediction) to the data.
        
        Args:
            self:
                svm_results: dict, here all SVM output is stored
            experiment: str
            misclassification: pd.DataFrame
            name: str, default=default
            prediction: pd.DataFrame, default=empty dataframe
            comment: str, default=empty string
            overwrite: bool, default=True
        Returns:
            self.svm_results is being set
        """
        if experiment not in self.exp_names:
            raise KeyError(f"Experiment {experiment} not found during SVM matrix storage")
        
        if experiment not in self.svm_results.keys():
            self.svm_results[experiment] = dict()
        
        if name in self.svm_results[experiment].keys() and not overwrite:
            raise KeyError(f"SVM result named {name} already exists for experiment {experiment}. Rename it or set overwrite=True.")
        
        self.svm_results[experiment][name] = {
            "comment": comment,
            "misclassification": misclassification,
            "prediction": prediction
        }
    
    
    def perform_pca_comparison(self):
        """
        PCA will be performed, using logarithmized data.

        Args:
            self:
                df_01_filtered_combined: df, which contains 0/1 normalized data for each map - for all experiments
                    columns: Fractions, e.g. "03K", "06K", "12K", "24K", "80K"
                    index: "Protein IDs", "Gene names", "Compartment", "Experiment", "Map", "Exp_Map"    
                df_01_mean_filtered_combined: df, which contains (global) 0/1 normalized data across all maps (mean) - for all experiments and for all protein IDs, 
                    that are consistent throughout all experiments
                    columns: Fractions, e.g. "03K", "06K", "12K", "24K", "80K"
                    index: "Gene names", "Protein IDs", "Compartment", "Experiment"
                    
        Returns:
            self:
                df_pca_for_plotting: PCA processed dataframe
                    index: "Experiment", "Gene names", "Map", "Exp_Map"
                    columns: "PC1", "PC2", "PC3"
                    contains only marker genes, that are consistent throughout all maps / experiments
                df_global_pca: PCA processed dataframe
                    index: "Gene names", "Protein IDs", "Compartment", "Experiment", 
                    columns: "PC1", "PC2", "PC3"
                    contains all protein IDs, that are consistent throughout all experiments
        """

        markerproteins = self.markerproteins.copy()
        
        #df_01_filtered_combined = self.df_01_filtered_combined
        #df_01_filtered_combined = self.df_01_filtered_combined 
        
        
        df_mean = pd.DataFrame()
        for exp in self.exp_names:
            df_exp = self.df_01_filtered_combined.stack("Fraction").unstack(["Experiment", "Map","Exp_Map"])[exp].mean(axis=1).to_frame(name=exp)
            df_mean = pd.concat([df_mean, df_exp], axis=1)
        df_mean = df_mean.rename_axis("Experiment", axis="columns").stack("Experiment").unstack("Fraction")
        
        pca = PCA(n_components=3)
        
        df_pca = pd.DataFrame(pca.fit_transform(df_mean.apply(zscore, axis=0).replace(np.nan, 0)))
        df_pca.columns = ["PC1", "PC2", "PC3"]
        df_pca.index = df_mean.index
        
        ###only one df, make annotation at that time
        df_cluster = pd.DataFrame([(k, i) for k, l in markerproteins.items() for i in l], columns=["Cluster", "Protein IDs"])
        df_global_pca = df_pca.reset_index().merge(df_cluster, how="left", on="Protein IDs")
        df_global_pca.Cluster.replace(np.NaN, "Undefined", inplace=True)
        
        self.df_pca = df_pca 
        self.df_global_pca = df_global_pca
            
            
    def plot_pca_comparison(self, cluster_of_interest_comparison="Proteasome", multi_choice=["Exp1", "Exp2"]):
        """
        A PCA plot for desired experiments (multi_choice) and 1 desired cluster is generated.
        Either the maps for every single experiment are displayed individually or in a combined manner
    
        Args:
            self:
                markerproteins: dictionary, key: cluster name, value: gene names (e.g. {"Proteasome" : ["PSMA1", "PSMA2",...], ...}
                multi_choice: list of experiment names
                cluster_of_interest_comparison: string, protein cluster (key in markerproteins, e.g. "Proteasome")
                df_pca: PCA processed dataframe
                    index: "Experiment", "Gene names", "Map", "Exp_Map"
                    columns: "PC1", "PC2", "PC3"
                    contains only marker genes, that are consistent throughout all maps / experiments
    
        Returns:
            pca_figure: PCA plot for a specified protein cluster.
        """
    
        df_pca = self.df_pca.copy()
        markerproteins = self.markerproteins
 
        try:
            df_setofproteins_PCA = pd.DataFrame()
            for map_or_exp in multi_choice:
                    for marker in markerproteins[cluster_of_interest_comparison]:
                        try:
                            plot_try_pca = df_pca.xs((marker, map_or_exp), level=["Protein IDs", "Experiment"], drop_level=False)
                        except KeyError:
                            continue
                        df_setofproteins_PCA = df_setofproteins_PCA.append(plot_try_pca)
            df_setofproteins_PCA.reset_index(inplace=True)
            
            df_setofproteins_PCA = df_setofproteins_PCA.assign(Experiment_lexicographic_sort=pd.Categorical(df_setofproteins_PCA["Experiment"], categories=multi_choice,
                                                                                                              ordered=True))
            df_setofproteins_PCA.sort_values("Experiment_lexicographic_sort", inplace=True)
                    
            pca_figure = px.scatter_3d(df_setofproteins_PCA, 
                                       x="PC1", 
                                       y="PC2", 
                                       z="PC3", 
                                       color="Experiment",
                                       template="simple_white",
                                       hover_data=["Gene names"]
                                      )
                    
            pca_figure.update_layout(autosize=False, 
                                     width=700, 
                                     height=500,
                                     title="PCA plot for <br>the protein cluster: {}".format(cluster_of_interest_comparison),
                                     template="simple_white"
                                    )
    
            return pca_figure
        except:
            return "This protein cluster was not identified in all experiments"
            
                
    def plot_global_pca_comparison(self, cluster_of_interest_comparison="Proteasome", x_PCA="PC1", y_PCA="PC3", 
        markerset_or_cluster=False, multi_choice=["Exp1", "Exp2"]):
        """"
        PCA plot will be generated
    
        Args:
            self:
                df_organellarMarkerSet: df, columns: "Gene names", "Compartment", no index
                multi_choice: list of experiment names
                css_color: list of colors
                df_global_pca: PCA processed dataframe
                    index: "Gene names", "Protein IDs", "Compartment", "Experiment", 
                    columns: "PC1", "PC2", "PC3"
                    contains all protein IDs, that are consistent throughout all experiments    
    
        Returns:
            pca_figure: global PCA plot, clusters based on the markerset based (df_organellarMarkerSet) are color coded. 
        """
        
        
        df_global_pca_exp = self.df_global_pca.loc[self.df_global_pca["Experiment"].isin(multi_choice)]
        df_global_pca_exp.reset_index(inplace=True)

        compartments = list(SpatialDataSet.df_organellarMarkerSet["Compartment"].unique())
        compartment_color = dict(zip(compartments, self.css_color))
        compartment_color["Selection"] = "black"
        compartment_color["undefined"] = "lightgrey"
        compartments.insert(0, "undefined")
        compartments.insert(len(compartments), "Selection")
            
        cluster = self.markerproteins.keys()
        cluster_color = dict(zip(cluster, self.css_color))
        cluster_color["Undefined"] = "lightgrey"
                
        
        if markerset_or_cluster == True:
            df_global_pca = df_global_pca_exp[df_global_pca_exp.Cluster!="Undefined"].sort_values(by="Cluster")
            df_global_pca = df_global_pca_exp[df_global_pca_exp.Cluster=="Undefined"].append(df_global_pca)
        else:
            for i in self.markerproteins[cluster_of_interest_comparison]:
                df_global_pca_exp.loc[df_global_pca_exp["Protein IDs"] == i, "Compartment"] = "Selection"
            df_global_pca = df_global_pca_exp.assign(Compartment_lexicographic_sort = pd.Categorical(df_global_pca_exp["Compartment"], 
                                                                                                     categories=[x for x in compartments], 
                                                                                                     ordered=True))
            df_global_pca.sort_values(["Compartment_lexicographic_sort", "Experiment"], inplace=True)
            
        fig_global_pca = px.scatter(data_frame=df_global_pca,
                                    x=x_PCA,
                                    y=y_PCA,
                                    color="Compartment" if markerset_or_cluster == False else "Cluster",
                                    color_discrete_map=compartment_color if markerset_or_cluster == False else cluster_color,
                                    title="Protein subcellular localization by PCA",
                                    hover_data=["Protein IDs", "Gene names", "Compartment"],
                                    facet_col="Experiment",
                                    facet_col_wrap=2,
                                    opacity=0.9, 
                                    template="simple_white"
                                    )
        
        fig_global_pca.update_layout(autosize=False, 
                                     width=1800 if markerset_or_cluster == False else 1600, 
                                     height=400*(int(len(multi_choice) / 2) + (len(multi_choice) % 2 > 0)),
                                     template="simple_white"
                                    )
        
        return fig_global_pca
    
    
    def get_marker_proteins(self, experiments, cluster):
        df_in = self.df_01_filtered_combined.copy()
        markers = self.markerproteins[cluster]
        
        # retrieve marker proteins
        df_cluster = pd.DataFrame()
        for marker in markers:
            try:
                df_p = df_in.xs(marker, level="Protein IDs", axis=0, drop_level=False)
            except:
                continue
            df_cluster = df_cluster.append(df_p)
        if len(df_cluster) == 0:
            return df_cluster
        
        # filter for all selected experiments
        df_cluster = df_cluster.droplevel("Exp_Map", axis=0)
        df_cluster = df_cluster.unstack(["Experiment", "Map"])
        if any([el not in df_cluster.columns.get_level_values("Experiment") for el in experiments]):
            return pd.DataFrame()
        drop_experiments = [el for el in df_cluster.columns.get_level_values("Experiment") if el not in experiments]
        if len(drop_experiments) > 0:
            df_cluster.drop([el for el in df_cluster.columns.get_level_values("Experiment") if el not in experiments],
                            level="Experiment", axis=1, inplace=True)
        df_cluster.dropna(inplace=True)
        if len(df_cluster) == 0:
            return df_cluster
        df_cluster.set_index(pd.Index(np.repeat(cluster, len(df_cluster)), name="Cluster"), append=True, inplace=True)
        
        return df_cluster
    
    
    def calc_cluster_distances(self, df_cluster, complex_profile=np.median, distance_measure="manhattan"):
        df_distances = pd.DataFrame()
        
        # loop over experiments
        experiments = set(df_cluster.columns.get_level_values("Experiment"))
        for exp in experiments:
            df_exp = df_cluster.xs(exp, level="Experiment", axis=1)
            ref_profile = pd.DataFrame(df_exp.apply(complex_profile, axis=0, result_type="expand")).T
            
            # loop over maps
            maps = set(df_exp.columns.get_level_values("Map"))
            for m in maps:
                if distance_measure == "manhattan":
                    d_m = pw.manhattan_distances(df_exp.xs(m, level="Map", axis=1), ref_profile.xs(m, level="Map", axis=1))
                else:
                    raise ValueError(distance_measure)
                d_m = pd.DataFrame(d_m, columns=[(exp, m)], index=df_exp.index)
                df_distances = pd.concat([df_distances, d_m], axis=1)
        
        df_distances.columns = pd.MultiIndex.from_tuples(df_distances.columns, names=["Experiment", "Map"])
        return df_distances
    
    
    def calc_biological_precision(self, experiments=None, clusters=None):
        """
        Method to calculate the distance table for assessing biological precision
        """
        
        df_distances = pd.DataFrame()
        if experiments is None:
            experiments = self.exp_names
        if clusters is None:
            clusters = self.markerproteins.keys()
        
        for cluster in clusters:
            df_cluster = self.get_marker_proteins(experiments, cluster)
            if len(df_cluster) == 0:
                continue
            dists_cluster = self.calc_cluster_distances(df_cluster)
            df_distances = df_distances.append(dists_cluster)
        df_distances = df_distances.stack(["Experiment", "Map"]).reset_index()\
            .sort_values(["Experiment","Gene names"]).rename({0: "distance"}, axis=1)
        df_distances.insert(0, "Exp_Map", ["_".join([e,m]) for e,m in zip(df_distances["Experiment"], df_distances["Map"])])
        
        self.df_distance_comp = df_distances
        
        return df_distances
    
    
    def get_complex_coverage(self, min_n=5):
        full_coverage = {}
        for complx in self.markerproteins.keys():
            df = self.get_marker_proteins(self.exp_names, complx)
            if len(df) >= min_n:
                full_coverage[complx] = len(df)
        partial_coverage = {}
        for exp in self.exp_names:
            for complx in self.markerproteins.keys():
                if complx in full_coverage.keys():
                    continue
                df = self.get_marker_proteins([exp], complx)
                #print(df)
                if complx in partial_coverage.keys():
                    partial_coverage[complx].append(len(df))
                else:
                    partial_coverage[complx] = [len(df)]
        no_coverage = {}
        for k in partial_coverage.keys():
            if all([el < min_n for el in partial_coverage[k]]):
                no_coverage[k] = partial_coverage[k]
        for k in no_coverage.keys():
            del partial_coverage[k]
        self.coverage_lists = [full_coverage, partial_coverage, no_coverage]
        return full_coverage, partial_coverage, no_coverage
    
    
    def distance_boxplot_comparison(self, cluster_of_interest_comparison="Proteasome", collapse_maps=False, multi_choice=["Exp1", "Exp2"]):
        """
        A box plot for desired experiments (multi_choice) and 1 desired cluster is generated displaying the distribution of the e.g.
        Manhattan distance. Either the maps for every single experiment are displayed individually or in a combined manner.

        Args:
            self:
                multi_choice: list of experiment names
                collapse_maps: boolean
                cluster_of_interest_comparison: string, protein cluster (key in markerproteins, e.g. "Proteasome")
                map_names: individual map names are stored as an index
                df_distance_comp: df_distance_comp: no index, column names: "Gene names", "Cluster", "Protein IDs", "Compartment", "Experiment", "Map", 
                                 "Exp_Map", "distance"
                                 "distance": Manhattan distances for each individual protein of the specified clusters (see self.markerproteins) are stored

        Returns:
            distance_boxplot_figure: boxplot. Along the x-axis the maps, along the y-axis the distances are shown
        """
        
        #an error massage, if no Experiments are selected, will be displayed already, that is why: return ""
        if len(multi_choice)>=1:
            pass
        else:
            return ("")
        
        df_distance_comp = self.df_distance_comp.copy()
        #set categroical column, allowing lexicographic sorting
        df_distance_comp["Experiment_lexicographic_sort"] = pd.Categorical(df_distance_comp["Experiment"], categories=multi_choice, ordered=True)
        df_distance_comp.sort_values(["Experiment_lexicographic_sort", "Map"], inplace=True)
        
        if collapse_maps == False:
            #get only values form experiment of interest
            df_distance_selectedExp = df_distance_comp.loc[df_distance_comp["Experiment"].isin(multi_choice)]
            #get only values form cluster of interest
            df_distance_selectedExp = df_distance_selectedExp.loc[df_distance_selectedExp["Cluster"]==cluster_of_interest_comparison]
            
            if df_distance_selectedExp.shape[0] == 0:
                self.cache_cluster_quantified = False
                
            else:
                individual_distance_boxplot_figure=go.Figure()
                for i, exp in enumerate(multi_choice):
                    df_plot=df_distance_selectedExp[df_distance_selectedExp["Experiment"]==exp]
                    individual_distance_boxplot_figure.add_trace(go.Box(
                        x=[df_plot["Experiment"], df_plot["Map"]],
                        y=df_plot["distance"],
                        #line=dict(color=pio.templates["simple_white"].layout["colorway"][i]),
                        boxpoints="all",
                        whiskerwidth=0.2,
                        marker_size=2,
                        name=exp, 
                        hovertext=df_plot["Gene names"]
                    ))
                    
                individual_distance_boxplot_figure.update_layout(boxmode="group", 
                                                                 xaxis_tickangle=90, 
                                                                 title="Manhattan distance distribution for <br>the protein cluster: {}".format(cluster_of_interest_comparison),
                                                                 autosize=False,
                                                                 width=350*len(multi_choice),
                                                                 height=500,
                                                                 xaxis=go.layout.XAxis(linecolor="black",
                                                                                       linewidth=1,
                                                                                       title="Experiment",
                                                                                       mirror=True),
                                                                 yaxis=go.layout.YAxis(linecolor="black",
                                                                                       linewidth=1,
                                                                                       title="Distance",
                                                                                       mirror=True),
                                                                 template="simple_white")
                
                return individual_distance_boxplot_figure
        
        else: 
            map_or_exp_names = multi_choice
            level_of_interest = "Experiment"
            boxplot_color = "Experiment"
            df_distance_selectedExp_global = df_distance_comp
    
            # "Gene names", "Map", "Cluster" and transferred into the index
            df_distance_selectedExp_global.set_index(["Gene names", level_of_interest, "Cluster"], inplace=True) 
    
            df_cluster_xmaps_distance_global = pd.DataFrame()
    
            # for each individual map and a defined cluster data will be extracted from the dataframe
            # "df_distance_selectedExp_global" and appended to the new dataframe df_cluster_xmaps_distance_global
            for map_or_exp in map_or_exp_names:
                plot_try = df_distance_selectedExp_global.xs((cluster_of_interest_comparison, map_or_exp), level=["Cluster", 
                                                                                                                       level_of_interest], drop_level=False)
                df_cluster_xmaps_distance_global = df_cluster_xmaps_distance_global.append(plot_try)
            
            df_cluster_xmaps_distance_global.sort_values("Experiment_lexicographic_sort", inplace=True)
            df_cluster_xmaps_distance_global.reset_index(inplace=True)
        
            distance_boxplot_figure = px.box(df_cluster_xmaps_distance_global, 
                                             x=level_of_interest, 
                                             y="distance", 
                                             points="all",
                                             hover_name="Gene names", 
                                             color=boxplot_color,
                                             template="simple_white",
                                             title="Global Manhattan distance distribution for the protein cluster: {}".format(cluster_of_interest_comparison)
                                             )
    
            distance_boxplot_figure.update_layout(autosize=False,
                                                  width=250*len(multi_choice),
                                                  height=500,
                                                  xaxis=go.layout.XAxis(linecolor="black",
                                                                        linewidth=1,
                                                                        title="Map",
                                                                        mirror=True),
                                                  yaxis=go.layout.YAxis(linecolor="black",
                                                                        linewidth=1,
                                                                        title="distance",
                                                                        mirror=True),
                                                  template="simple_white"
                                                 )
            
            return distance_boxplot_figure
    
    
    def plot_biological_precision(self, multi_choice=None, clusters_for_ranking=None, min_members=5, reference=""):
        if multi_choice is None:
            multi_choice = self.exp_names
        if clusters_for_ranking is None:
            clusters_for_ranking = self.clusters_for_ranking
        if len(multi_choice) == 0 or len(clusters_for_ranking) == 0:
            return("Please provide at least one experiment and one cluster for ranking")
        
        df = self.df_distance_comp.copy()
        df = df[df["Experiment"].isin(multi_choice)]
        df = df[df["Cluster"].isin(clusters_for_ranking)]
        
        df_m = df.groupby(["Cluster", "Experiment", "Map"]).filter(lambda x: len(x)>=min_members)
        df_c = df_m.groupby(["Cluster", "Experiment"]).median().reset_index()
        df_m = df_m.groupby(["Cluster", "Experiment", "Map"]).median().reset_index()
        
        df_m = df_m.assign(Experiment_lexicographic_sort = pd.Categorical(df_m["Experiment"], categories=multi_choice, ordered=True))
        df_m = df_m.sort_values("Experiment_lexicographic_sort").drop("Experiment_lexicographic_sort", axis=1)\
               .groupby("Experiment", as_index=False, group_keys=False, sort=False).apply(lambda x: x.sort_values("distance", ascending=False))
        df_c = df_c.assign(Experiment_lexicographic_sort = pd.Categorical(df_c["Experiment"], categories=multi_choice, ordered=True))
        df_c = df_c.sort_values("Experiment_lexicographic_sort").drop("Experiment_lexicographic_sort", axis=1)\
               .groupby("Experiment", as_index=False, group_keys=False, sort=False).apply(lambda x: x.sort_values("distance", ascending=False))
        
        bp_stacked_bar = px.bar(df_m, x="Experiment", y="distance", color="Cluster", hover_data=["Map"],
                                width=400+80*len(multi_choice), template="simple_white", height=100+30*len(clusters_for_ranking)).update_layout(legend_traceorder="reversed")
        
        bp_box_minus_min = px.box(df_m.set_index(["Experiment", "Cluster", "Map"]).unstack(["Experiment", "Map"])\
                                      .apply(lambda x: x-x.min(), axis=1).stack(["Experiment", "Map"]).reset_index()\
                                      .sort_values(["Experiment"], key=lambda x: [multi_choice.index(el) for el in x]),
                                  x="Experiment", y="distance", color="Experiment", hover_data=["Cluster", "Map"],
                                  width=200+100*len(multi_choice), template="simple_white", height=400, points="all")\
                                  .update_yaxes(title="distance - cluster offset (minimum)")
        bp_box_minus_ref = px.box(df_c.set_index(["Experiment", "Cluster"]).unstack(["Experiment"])\
                                      .apply(lambda x: x/x[("distance", reference)], axis=1).stack(["Experiment"]).reset_index()\
                                      .sort_values(["Experiment"], key=lambda x: [multi_choice.index(el) for el in x])\
                                      .loc[lambda x: x.Experiment != reference],
                                  x="Experiment", y="distance", color="Experiment", hover_data=["Cluster"],
                                  color_discrete_sequence=[px.colors.qualitative.D3[multi_choice.index(el)]
                                      for el in multi_choice if el != reference],
                                  width=200+100*len(multi_choice), template="simple_white", height=400, points="all")\
                                  .update_yaxes(title="distance relative to {}".format(reference))
        
        return bp_stacked_bar, bp_box_minus_min, bp_box_minus_ref
        
        
    
    def distance_ranking_barplot_comparison(self, collapse_cluster=False, multi_choice=["Exp1", "Exp2"], clusters_for_ranking=None, ranking_boxPlot="Box plot"):#, toggle_sumORmedian=False):
    #ref_exp="Exp1", 
        if clusters_for_ranking is None:
            clusters_for_ranking = self.clusters_for_ranking
            
            #an error massage, if no Experiments are selected, will be displayed already, that is why: return ""
        if len(multi_choice)>=1:
            pass
        else:
            return ("")
        
    #dict_cluster_normalizedMedian = {}
        #multi_choice = i_multi_choice.value
        #clusters_for_ranking =  i_clusters_for_ranking.value
        df_distance_comp = self.df_distance_comp.copy()
        df_distance_comp = df_distance_comp[df_distance_comp["Experiment"].isin(multi_choice)]
        df_distance_comp = df_distance_comp[df_distance_comp["Cluster"].isin(clusters_for_ranking)]
        
        df_quantified_cluster = df_distance_comp.reset_index()
        df_quantified_cluster = df_distance_comp.drop_duplicates(subset=["Cluster", "Experiment"]).set_index(["Cluster", 
                                                                                                                "Experiment"])["distance"].unstack("Cluster")
        self.df_quantified_cluster = df_quantified_cluster.notnull().replace({True: "x", False: "-"})
        
        
        dict_quantified_cluster = {}
        dict_cluster_normalizedMedian_ref = {}
        dict_median_distance_ranking = {}
        for cluster in clusters_for_ranking:
            try:
                df_cluster = df_distance_comp[df_distance_comp["Cluster"]==cluster]
                cluster_quantitity = df_cluster["Protein IDs"].unique().size
                if  cluster_quantitity>= 5:
                    dict_quantified_cluster[cluster] = cluster_quantitity
                    all_median_one_cluster_several_exp = {}
                    #ref = df_cluster["distance"].median()
                    for exp in multi_choice:
                        median = df_cluster[df_cluster["Experiment"]==exp]["distance"].median()
                        all_median_one_cluster_several_exp[exp] = float(median)
                        #new
                        #if exp == ref_exp:
                        #    ref = median
                    ref = np.median(list(all_median_one_cluster_several_exp.values()))
                    dict_median_distance_ranking[cluster] = all_median_one_cluster_several_exp
                    
                    median_ranking_ref = {exp: median/ref for exp, median in all_median_one_cluster_several_exp.items()}
                    dict_cluster_normalizedMedian_ref[cluster] = median_ranking_ref
                else:
                    continue
            except:
                continue
        
        self.cluster_above_treshold = dict_quantified_cluster.keys()
        self.df_quantified_cluster2 = pd.DataFrame.from_dict({"Number of PG per Cluster":dict_quantified_cluster}).T
        
        df_cluster_normalizedMedian_ref = pd.DataFrame(dict_cluster_normalizedMedian_ref)
        df_cluster_normalizedMedian_ref.index.name="Experiment"
        df_cluster_normalizedMedian_ref.rename_axis("Cluster", axis=1, inplace=True)
        
        #median makes a huge differnece, improves result of DIA, MQ, libary
        df_RelDistanceRanking = pd.concat([df_cluster_normalizedMedian_ref.median(axis=1), df_cluster_normalizedMedian_ref.sem(axis=1)], axis=1, 
                                        keys=["Distance Ranking (rel, median)", "SEM"]).reset_index().sort_values("Distance Ranking (rel, median)")
        
        ranking_sum = df_cluster_normalizedMedian_ref.sum(axis=1).round(2)
        ranking_sum.name = "Normalized Median - Sum"
        df_ranking_sum = ranking_sum.reset_index()
        
        #ranking_product = df_cluster_normalizedMedian.product(axis=1).round(2)
        #ranking_product.name = "Normalized Median - Product"
        #df_globalRanking = pd.concat([pd.DataFrame(ranking_sum), pd.DataFrame(ranking_product)], axis=1).reset_index()
        
        df_cluster_normalizedMedian_ref = df_cluster_normalizedMedian_ref.stack("Cluster")
        df_cluster_normalizedMedian_ref.name="Normalized Median"
        df_cluster_normalizedMedian_ref = df_cluster_normalizedMedian_ref.reset_index()
        self.df_cluster_normalizedMedian_ref = df_cluster_normalizedMedian_ref
        df_cluster_normalizedMedian_ref = df_cluster_normalizedMedian_ref.assign(Experiment_lexicographic_sort = pd.Categorical(df_cluster_normalizedMedian_ref["Experiment"], categories=multi_choice, ordered=True))
        df_cluster_normalizedMedian_ref.sort_values("Experiment_lexicographic_sort", inplace=True)
        
        
        if collapse_cluster == False:
            
            fig_ranking = px.bar(df_cluster_normalizedMedian_ref, 
                                x="Cluster", 
                                y="Normalized Median", 
                                color="Experiment", 
                                barmode="group", 
                                title="Ranking - normalization to reference experiments the median across all experiments for each cluster",
                                template="simple_white"
                                )
            
            fig_ranking.update_xaxes(categoryorder="total ascending")
                
            fig_ranking.update_layout(autosize=False,
                                    width=1200 if len(multi_choice)<=3 else 300*len(multi_choice),
                                    height=500,
                                    template="simple_white"
                                    )
            return fig_ranking
 
        else:
            if ranking_boxPlot == "Bar plot - median":
                fig_globalRanking = px.bar(df_RelDistanceRanking.sort_values("Distance Ranking (rel, median)"), 
                                            x="Experiment",
                                            y="Distance Ranking (rel, median)", 
                                            title="Median manhattan distance distribution for <br>all protein clusters (n>=5 per cluster)",# - median of all individual normalized medians - reference experiment is the median across all experiments for each cluster",
                                            error_x="SEM", error_y="SEM", 
                                            color="Experiment", 
                                            template="simple_white")
                    
        
                                                
            if ranking_boxPlot == "Box plot": 
                fig_globalRanking = px.box(df_cluster_normalizedMedian_ref,
                                        x="Experiment",
                                        y="Normalized Median", 
                                        title="Median manhattan distance distribution for <br>all protein clusters (n>=5 per cluster)",# "Ranking - median of all individual normalized medians - reference is the median across all experiments for each cluster",
                                        color="Experiment",
                                        points="all",
                                        template="simple_white",
                                        hover_name="Cluster")
                #return pn.Column(pn.Row(fig_globalRanking), pn.Row(fig_globalRanking2))
            else:
                fig_globalRanking = px.bar(df_ranking_sum.sort_values("Normalized Median - Sum"), 
                    x="Experiment",
                    template="simple_white",
                    y="Normalized Median - Sum", 
                    title="Ranking - median of all individual normalized medians - reference is the median across all experiments for each cluster",
                    color="Experiment")
            
            fig_globalRanking.update_layout(autosize=False,
                                            width=250*len(multi_choice),
                                            height=500,
                                            template="simple_white"
                                            )
            
            return fig_globalRanking
        
    
    def quantity_pr_pg_barplot_comparison(self, multi_choice=["Exp1", "Exp2"]):
        """
        Barplot, showing number of protein groups/profiles. 
        
        Args:
            self:
                df_quantity_pr_pg_combined: df, no index, column names: "filtering", "type", "number of protein groups", "number of profiles", 
                                            "data completeness of profiles", "Experiment"    
                multi_choice: list of experiment names
                
        Returns: 
            fig_quantity_pr_pg: barplot, number of protein groups/profiles before/after filtering of the intersection/total quantity
        """
        
        df_quantity_pr_pg_combined = self.df_quantity_pr_pg_combined.copy()
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined[df_quantity_pr_pg_combined["Experiment"].isin(multi_choice)]
        
        df_quantity_pr_pg_combined.insert(0,"Expxfiltering",[" ".join([e,f]) for e,f in zip(
            df_quantity_pr_pg_combined.Experiment, df_quantity_pr_pg_combined.filtering)])
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined.assign(
            Experiment_lexicographic_sort = pd.Categorical(df_quantity_pr_pg_combined["Experiment"], categories=multi_choice, ordered=True))

        df_quantity_pr_pg_combined.sort_values(["Experiment_lexicographic_sort", "type"], ascending=[True, False], inplace=True)
        
        layout = go.Layout(barmode="overlay", 
          #xaxis_tickangle=90, 
          autosize=False,
          width=100*len(multi_choice)+150,
          height=400,
          template="simple_white")
        filtered = list(np.tile(["id","profile"],len(multi_choice)))
        
        fig_quantity_pg = px.bar(df_quantity_pr_pg_combined, x="Expxfiltering", y="number of protein groups",
                                 color="Experiment", barmode="overlay", hover_data=["type"],
                                 opacity=0.8, color_discrete_sequence=px.colors.qualitative.D3)
        fig_quantity_pg.update_layout(layout, title="Number of Protein Groups",
                                      xaxis={"tickmode":"array", "tickvals":[el for el in range(len(multi_choice)*2)],
                                             "ticktext":filtered, "title": {"text": None}})
        
         
         
        fig_quantity_pr = px.bar(df_quantity_pr_pg_combined, x="filtering", y="number of profiles",
                                 color="type", barmode="overlay", labels={"Experiment":"", "filtering":""}, 
                                 facet_col="Experiment",template="simple_white", opacity=1)\
                                .for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig_quantity_pr.update_layout(layout, title="Number of Profiles" )
        
        return fig_quantity_pg, fig_quantity_pr
    
    
    def coverage_comparison(self, multi_choice=["Exp1", "Exp2"]):
        """
        Barplot, showing data completeness of profiles. 
        
        Args:
            self:
                df_quantity_pr_pg_combined: df, no index, column names: "filtering", "type", "number of protein groups", "number of profiles", 
                                            "data completeness of profiles", "Experiment"
                multi_choice: list of experiment names
                
        Returns: 
            fig_pr_dc: barplot, data completeness of profiles before/after filtering of intersection/total qunatity
        """
        
        df_quantity_pr_pg_combined = self.df_quantity_pr_pg_combined.copy()
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined[df_quantity_pr_pg_combined["Experiment"].isin(multi_choice)]

        df_quantity_pr_pg_combined.insert(0,"Expxfiltering",[" ".join([e,f]) for e,f in zip(
            df_quantity_pr_pg_combined.Experiment, df_quantity_pr_pg_combined.filtering)])
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined.assign(
            Experiment_lexicographic_sort = pd.Categorical(df_quantity_pr_pg_combined["Experiment"], categories=multi_choice, ordered=True))

        #df_quantity_pr_pg_combined.sort_values("Experiment_lexicographic_sort", inplace=True)
        df_quantity_pr_pg_combined.sort_values(["Experiment_lexicographic_sort", "filtering"], ascending=[True, False], inplace=True)
        filtered = list(np.tile(["id","profile"],len(multi_choice)))
        
        fig_pr_dc = px.bar(df_quantity_pr_pg_combined.loc[df_quantity_pr_pg_combined.type=="total"], x="Expxfiltering", y="data completeness of profiles",
                           color="Experiment", barmode="overlay", hover_data=["filtering"],
                           template="simple_white", opacity=0.8, color_discrete_sequence=px.colors.qualitative.D3)
        
        fig_pr_dc.update_layout(
            title="Profile completeness of all<br>identified protein groups", autosize=False,
            width=100*len(multi_choice)+150, height=400,
            template="simple_white", xaxis={"tickmode":"array", "tickvals":[el for el in range(len(multi_choice)*2)],
                                            "ticktext":filtered, "title": {"text": None}})
        
        return fig_pr_dc
    
    
    def venn_sections(self, multi_choice_venn=["Exp1"], omit_size=300):
        """
        UpsetPlot is created based on list of experiments. If 2/3 experiments are given, the Upsetlot displays all possible
        mutually exclusive overlapping combinations of these experiments. Additionally a Venn Diagram is created using matplotlib. 
        Latter figure has to be transformed from matplotlib object to jpg, to make it available for the webinterface via panel/holoviz.
        If more than 3 experiments are given, the UpsetPlot will be calculated only for those combinations of these experiments with at least 300 entries.
        
        Another way to think of this is the mutually exclusive sections of a venn diagram of the sets.  If the original list has N sets, 
        the returned list will have (2**N)-1 sets.
        
        Args:
            multi_choice_venn: list of experiment names 
            self:
                df_01_filtered_combined: pd.DataFrame, Uses Protein IDs index level
        
        Returns:
            im: Venn diagram, made availabe flor plotly/webinterface
            figure_UpSetPlot: Upsetplot figure
            
        combinations : list of tuple
            tag : str
                Binary string representing which sets are included / excluded in
                the combination.
            set : set
                The set formed by the overlapping input sets.
        """
        
        def create_upsetplot(sets, multi_choice):
            num_combinations = 2 ** len(sets)
            bit_flags = [2 ** n for n in range(len(sets))]
            flags_zip_sets = [z for z in zip(bit_flags, sets)]
            combo_sets = []
            overlapping_ids = []
            experiments = []
            #dictio = {}
            for bits in range(num_combinations - 1, 0, -1):
                include_sets = [s for flag, s in flags_zip_sets if bits & flag]
                exclude_sets = [s for flag, s in flags_zip_sets if not bits & flag]
                combo = set.intersection(*include_sets)
                combo = set.difference(combo, *exclude_sets)
                tag = "".join([str(int((bits & flag) > 0)) for flag in bit_flags])
                
                experiment_decoded = []
                for digit, exp in zip(list(tag), multi_choice):
                    if digit=="0":
                        continue
                    else:
                        experiment_decoded.append(exp)
                #dictio[len(combo)] = experiment_decoded
                if len(multi_choice)>3:
                    if len(combo)>omit_size:
                        overlapping_ids.append(len(combo))
                        experiments.append(experiment_decoded)
                else:
                    if len(combo)>0:
                        overlapping_ids.append(len(combo))
                        experiments.append(experiment_decoded)
                #combo_sets.append((tag, len(combo)))
                
            fig_UpSetPlot = plt.Figure()
            series_UpSetPlot = from_memberships(experiments, data=overlapping_ids)
            upplot(series_UpSetPlot, fig=fig_UpSetPlot, show_counts="%d")
            return fig_UpSetPlot
            
        
        if "Sequence" not in self.df_01_filtered_combined.index.names:
            sets_proteins_total = [set(self.df_01_filtered_combined.xs(i, axis=0, level="Experiment").index.get_level_values("Protein IDs"))
                                for i in multi_choice_venn]
            sets_proteins_intersection = [set(self.df_01_filtered_combined.xs(i, axis=0, level="Experiment").unstack(["Map", "Exp_Map"]).dropna()\
                                        .index.get_level_values("Protein IDs")) for i in multi_choice_venn]
        else:
            sets_proteins_total = [set(self.df_01_filtered_combined.xs(i, axis=0, level="Experiment").index.get_level_values("Sequence"))
                                for i in multi_choice_venn]
            sets_proteins_intersection = [set(self.df_01_filtered_combined.xs(i, axis=0, level="Experiment").unstack(["Map", "Exp_Map"]).dropna()\
                                        .index.get_level_values("Sequence")) for i in multi_choice_venn]
        figure_UpSetPlot_total = create_upsetplot(sets_proteins_total, multi_choice_venn)
        figure_UpSetPlot_int = create_upsetplot(sets_proteins_intersection, multi_choice_venn)
        
        #make matplot figure available for plotly
        def convert_venn_jpg(vd):
            vd = vd.figure
            out_img = BytesIO()
            plt.savefig(out_img, bbox_inches="tight",format="jpg", dpi=72)
            out_img.seek(0)  # rewind file
            im = Image.open(out_img)
            plt.clf()
            return im
        
        if len(multi_choice_venn) == 2:
            vd_t = venn2(sets_proteins_total, set_labels=([i for i in multi_choice_venn]),
                         set_colors=px.colors.qualitative.D3[0:2], alpha=0.8)
            vd_t = plt.title("in at least one map")
            im_t = convert_venn_jpg(vd_t)
            vd_i = venn2(sets_proteins_intersection, set_labels=([i for i in multi_choice_venn]),
                         set_colors=px.colors.qualitative.D3[0:2], alpha=0.8)
            vd_i = plt.title("in all maps")
            im_i = convert_venn_jpg(vd_i)
        elif len(multi_choice_venn) == 3:
            vd_t = venn3(sets_proteins_total, set_labels=([i for i in multi_choice_venn]),
                         set_colors=px.colors.qualitative.D3[0:3], alpha=0.8)
            vd_t = plt.title("in at least one map")
            im_t = convert_venn_jpg(vd_t)
            vd_i = venn3(sets_proteins_intersection, set_labels=([i for i in multi_choice_venn]),
                         set_colors=px.colors.qualitative.D3[0:3], alpha=0.8)
            vd_i = plt.title("in all maps")
            im_i = convert_venn_jpg(vd_i)
        
        else:
            im = "Venn diagram can be displayed for 3 Experiments or less"
            return im,im, figure_UpSetPlot_total, figure_UpSetPlot_int
        
        return im_t, im_i, figure_UpSetPlot_total, figure_UpSetPlot_int
    
    
    def calculate_global_scatter(self, metric, consolidation):
        """
        A distribution plot of the profile scatter in each experiment is generated, with variable distance metric and consolidation of replicates.
        
        Args:
            self:
                df_01_filtered_combined: df, indexed
            multi_choice: list of experiment names
            metric: distance metric, one of 'euclidean distance', 'manhattan distance', '1 - cosine correlation', '1 - pearson correlation'
            consolidation: method to consolidate replicate distances, one of 'median', 'average', 'sum'
        
        Returns:
            plot: plotly.figure_factory.displot, shows kernel density estiamtion in the main pane and a rug plot underneath. Traces are sorted by ascending median of the distribution.
        """
        
        # Option dictionaries
        cons_functions = {
            "median": np.median,
            "average": np.mean,
            "sum": np.sum
        }
        metrics = {
            "euclidean distance": "euclidean",
            "manhattan distance": "manhattan",
            "1 - cosine correlation": "cosine",
            "1 - pearson correlation": lambda x,y: 1-np.corrcoef(x,y)[0][1],
            "manhattan distance to average profile": [np.mean, pw.paired_manhattan_distances],
            "manhattan distance to median profile": [np.median, pw.paired_manhattan_distances]
        }
        
        # Option assertion
        assert consolidation in cons_functions.keys()
        assert metric in metrics.keys()
        
        self.parameters["reproducibility metric"] = metric
        self.parameters["reproducibility consolidation"] = consolidation
        
        df = self.df_01_filtered_combined.copy()
        df.index = df.index.droplevel(["Exp_Map", "Gene names", "Compartment"])
        if "Sequence" in df.index.names:
            df.index = df.index.droplevel(["Protein IDs"])
        
        # Calculate and consolidate distances
        distances = pd.DataFrame()
        for exp in self.exp_names:
            df_m = df.xs(exp, level="Experiment", axis=0).unstack("Map").dropna().stack("Map")
            maps = list(set(df_m.index.get_level_values("Map")))
            
            # this if clause switches between pairwise comparisons of profiles (else) and comparisons to an average/median profile
            if " to " in metric:
                df_m = df_m.unstack("Map")
                
                # calculate reference profiles
                df_profiles = df_m.stack("Fraction").apply(metrics[metric][0], axis=1).unstack("Fraction")
                
                # calculate the distance for every map
                distances_m = pd.DataFrame()
                for m in maps:
                    dist_m = pd.DataFrame(metrics[metric][1](df_m.xs(m, level="Map", axis=1), df_profiles), columns = [m])
                    distances_m = pd.concat([distances_m, dist_m], axis=1)
                
                distances_m.index = df_m.index
                
            else:
                distances_m = pd.DataFrame()
                
                # loop over pairs of maps
                for i,mapi in enumerate(maps):
                    for j,mapj in enumerate(maps):
                        # only look at each comparison once
                        if j <= i:
                            continue
                        dist = pw.paired_distances(df_m.xs(mapi, level="Map", axis=0).values,
                                                df_m.xs(mapj, level="Map", axis=0).values,
                                                metric = metrics[metric])
                        dist = pd.Series(dist, name="_".join([mapi,mapj]))
                        distances_m = pd.concat([distances_m, dist], axis=1)
                distances_m.index = df_m.xs(maps[0], level="Map", axis=0).index
            
            distances = pd.concat([distances,
                                   pd.Series(distances_m.apply(cons_functions[consolidation], axis=1), name=exp, index=distances_m.index)],
                                  axis=1)
        
        self.distances = distances
    
    
    def plot_reproducibility_distribution(self, multi_choice=[], q=0.5,
                                          show_full=True, show_rug=False, x_cut=None):
        
        if len(multi_choice) == 0:
            multi_choice = self.exp_names
        
        distances = self.distances[multi_choice].copy()
        if x_cut == None:
            x_cut = 1.05*distances.max().max()
        
        # get data and labels set up for figure factory
        if show_full:
            plotdata = distances.join(distances.dropna(), lsuffix=" full", rsuffix=" overlap")
        else:
            plotdata = distances.dropna()
        plotlabels = list(plotdata.columns)
        plotdata = [vals.dropna() for k,vals in plotdata.T.iterrows()]
        
        # calculate quantiles
        quantiles = dict()
        for i, (name, data) in enumerate(zip(plotlabels, plotdata)):
            quantiles[name] = np.quantile(data, q)
        
        # create plot
        plot = ff.create_distplot(plotdata, plotlabels, show_hist=False, show_rug=show_rug, show_curve=True)
        plot.update_layout(title="Distribution of {} {}, overlap = {}".format(self.parameters["reproducibility consolidation"],
                                                                              self.parameters["reproducibility metric"],
                                                                              len(plotdata[-1])),
                           width=1000, height=600, template="simple_white",
                           xaxis={"range": (0,x_cut), "rangemode": "nonnegative"},
                           )
        
        # get color dictionary and assign same color to overlap and full coverage traces, assign legend groups
        colors = []
        plot.for_each_trace(lambda x: colors.append(x.marker.color),
                            selector=lambda x: not x.name.endswith(" overlap"))
        if show_full:
            colors[distances.shape[1]::] = colors[:distances.shape[1]]
        colors = {k: v for k,v in zip(plotlabels, colors)}
        plot.for_each_trace(lambda x: x.update(legendgroup=x.name, marker_color=colors[x.name]))
        
        # get density estimates at quantile and add highlighting lines
        densities = []
        plot.for_each_trace(lambda x: densities.append(x.y[bisect.bisect(x.x, quantiles[x.name])]),
                            selector=dict(type="scatter"))
        for i, (name, data) in enumerate(zip(plotlabels, plotdata)):
            plot.add_scatter(x=[quantiles[name], quantiles[name]], y=[0,densities[i]],
                             text=["", "{:.4f}".format(quantiles[name])],
                             name=name, legendgroup=name, showlegend=False,
                             line_color=colors[name],
                             mode="lines+text", textposition='top right')
        
        # dash full distribution traces
        if show_full:
            plot.for_each_trace(lambda x: x.update(line_dash="dash"),
                                selector=lambda x: not x.name.endswith(" overlap"))
        
        return plot
    
    
    def plot_overview(self, multi_choice, clusters, quantile):
        dists = self.df_distance_comp.query('Cluster in @clusters').query('Experiment in @multi_choice')\
            .groupby(["Experiment", "Map", "Cluster"]).median()\
            .groupby("Experiment").sum().rename({"distance": "complex scatter"}, axis=1)
        
        rep = self.distances[multi_choice].dropna().apply(lambda x: np.quantile(x, quantile), axis=0)
        rep.name = "intermap scatter"
        
        depth = self.df_quantity_pr_pg_combined.query('type == "intersection"').query('filtering == "after filtering"').query('Experiment in @multi_choice')\
            [["Experiment", "number of protein groups"]].set_index("Experiment")\
            .rename({"number of protein groups": "profiled depth"}, axis=1)
        
        df_plot = dists.join(rep).join(depth).melt(ignore_index=False).reset_index()
        
        fig = make_subplots(
            1,3,
            subplot_titles=["profiled depth", "complex scatter", "intermap scatter"],
            horizontal_spacing=0.17
        )
        
        for i,y in enumerate(["profiled depth", "complex scatter", "intermap scatter"]):
            for j,e in enumerate(df_plot.Experiment.unique()):
                fig.add_trace(go.Bar(x=df_plot.query('variable == @y and Experiment == @e').Experiment,
                                     y=df_plot.query('variable == @y and Experiment == @e').value, name=e,
                                     marker=dict(color=px.colors.qualitative.D3[j], coloraxis="coloraxis"),
                                     legendgroup=e,
                                     hovertemplate="%{y:.5r}"
                                    ), 1, i+1)
        fig.update_xaxes(matches='x').update_layout(template="simple_white", showlegend=False)
        fig.for_each_yaxis(lambda x: x.update(title="full coverage protein groups" if x.anchor=="x"\
                                              else "∑ median intracomplex distances" if x.anchor=="x2"\
                                              else f"{quantile*100}% quantile of shared protein groups",
                                              tickformat=".0s" if x.anchor == "x" else ".2r",
                                              nticks=8, title_standoff=8 if x.anchor=="x" else 0
                                             ))
        fig.update_layout(width=320+(200*(j+1)), height=500, margin=dict(t=70,b=0,r=0,l=0),
                          title="Overview of benchmarking output")
        return fig
    
    def svm_processing(self):
        """
        The misclassification matrix, generated by Perseus, will be used for Recall/Precision calculation of each individual cluster and on a global level. 
        Data will be stored in a local dictionary that will be assigned to the global dictionary. 
        
        Args:
            self.df_SVM: dataframe, provided by Perseus, no index; 
                        Column names: e.g. "Predicted: ER", "Predicted: NPC"
                        Rows: e.g. "True: ER", "True: NPC"
        
        Returns:
            self.analysed_datasets_dict:
                local dictionary (SVM_dict) will be assigned to the global dictionary self.analysed_datasets_dict, that is available for downloading
                {"Experiment name" : {see read_jsonFile(self) [below]}
                                     {"Misclassification Analysis": 
                                         {
                                          "True: ER" : {
                                                       "Recall": int,
                                                       "FDR": int,
                                                       "Precision": int,
                                                       "F1": int
                                                       }
                                           "True: NPC" : {...}
                                            ...          
                                           "Summary": {...}
                                           }
                                        }
                }
        """
        
        global_SVM_dict_total = {}
        global_SVM_dict = {}
        for exp in self.svm_results.keys():
            try:
                df_SVM = self.svm_results[exp]["default"]["misclassification"]
                df_SVM["T: True group"] = df_SVM["T: True group"].str.replace(r'True: ', '')
            except KeyError:
                continue
            SVM_dict = {}
            all_correct = np.diag(df_SVM)
            members = df_SVM.sum(axis=1)
            total_members = 0
            membrame_members = 0
            membrane_correct = 0
            all_organelle_recall = []
            all_organelle_precision = []
            all_organelle_f1 = []
            F1_all_cluster = []
            no_of_membrane_clusters = 0
            total_correct = sum(all_correct)
            predicted_one_organelle = df_SVM.sum(axis=0)
            
            for i in range(len(df_SVM)):
                total_members = total_members + members[i]
                recall = all_correct[i]/members[i]
                fdr = (predicted_one_organelle[i]-all_correct[i])/predicted_one_organelle[i]
                precision = 1-fdr
                F1 = statistics.harmonic_mean([recall, precision])
                F1_all_cluster.append(F1)
                SVM_dict[df_SVM["T: True group"][i]] = {"Recall": recall, "FDR": fdr, "Precision": precision, "F1": F1} 
                if df_SVM["T: True group"][i]!="Nuclear pore complex" and df_SVM["T: True group"][i]!="Large Protein Complex" and df_SVM["T: True group"][i]!="Actin binding proteins" :
                    no_of_membrane_clusters = no_of_membrane_clusters+1
                    membrame_members = membrame_members + members[i]
                    membrane_correct = membrane_correct + all_correct[i]
                    all_organelle_f1.append(F1)
                    all_organelle_recall.append(recall)
                    all_organelle_precision.append(precision)
                
            total_recall = total_correct/total_members
            membrane_recall = membrane_correct/membrame_members
            av_per_organelle_recall = statistics.mean(all_organelle_recall)
            median_per_organelle_recall = statistics.median(all_organelle_recall)
            av_per_organelle_precision = statistics.mean(all_organelle_precision)
            avg_organelle_f1 = statistics.mean(all_organelle_f1)
            avg_F1_all_cluster = statistics.mean(F1_all_cluster)
            
            SVM_dict_total = {}
            SVM_dict_total["Avg. all clusters"] = {"Recall": total_recall, "F1": avg_F1_all_cluster} #total recall = marker prediction accuracy
            SVM_dict_total["Avg. all organelles"] = {"Recall": av_per_organelle_recall, "F1": avg_organelle_f1, "Precision": av_per_organelle_precision}
            SVM_dict_total["Membrane"] = {"Recall": membrane_recall}
            SVM_dict_total["Median. per organelle"] = {"Recall": median_per_organelle_recall}
            
            global_SVM_dict[exp] = SVM_dict
            global_SVM_dict_total[exp] = SVM_dict_total
            self.global_SVM_dict = global_SVM_dict
            self.global_SVM_dict_total = global_SVM_dict_total
            
        #f global_SVM_dict=={}:
        #   self.cache_stored_SVM = False
        #   return
        #lse:
        df_clusterPerformance_global = pd.DataFrame.from_dict({(i,j): global_SVM_dict[i][j] 
                                for i in global_SVM_dict.keys() 
                                for j in global_SVM_dict[i].keys()},
                            orient='index')
        df_clusterPerformance_global.index.names = ["Experiment", "Type"]
        self.df_clusterPerformance_global = df_clusterPerformance_global.T 
        
        df_AvgClusterPerformance_global = pd.DataFrame.from_dict({(i,j): global_SVM_dict_total[i][j] 
                                for i in global_SVM_dict_total.keys() 
                                for j in global_SVM_dict_total[i].keys()},
                            orient='index')
        df_AvgClusterPerformance_global.index.names = ["Experiment", "Type"]
        self.df_AvgClusterPerformance_global = df_AvgClusterPerformance_global.T
        #elf.cache_stored_SVM = True
        return 
            
            
    def svm_plotting(self, multi_choice):
        """
        The markerperformance (line/scatter plot) as well as marker prediction accuracy (bar plot) is visuaized.
        
        Args:
            self: df_AvgClusterPerformance_global 
                  df_clusterPerformance_global
            multi_choice: list of experiment names
        """
        df_clusterPerformance_global = self.df_clusterPerformance_global
        df_AvgClusterPerformance_global = self.df_AvgClusterPerformance_global
        
        df_AvgAllCluster = df_AvgClusterPerformance_global.xs("Avg. all clusters", level='Type', axis=1)
        fig_markerPredictionAccuracy = go.Figure()#data=[go.Bar(x=df_test.columns, y=df_test.loc["Recall"])])
        for exp in multi_choice:
            fig_markerPredictionAccuracy.add_trace(go.Bar(x=[exp], y=[df_AvgAllCluster[exp].loc["Recall"]], name=exp))
        fig_markerPredictionAccuracy.update_layout(template="simple_white", #showlegend=False, 
                    title="Marker prediction accuracy - Overall recall",
                    xaxis=go.layout.XAxis(linecolor="black",
                                        linewidth=1,
                                        mirror=True),
                    yaxis=go.layout.YAxis(linecolor="black",
                                        linewidth=1,
                                        title="Marker prediction accuracy [%]",
                                        mirror=True),
                    )
        
        fig_clusterPerformance = go.Figure()
        list_data_type = ["Avg. all clusters", "Avg. all organelles"]
        for i,exp in enumerate(multi_choice):
            df_clusterPerformance = df_clusterPerformance_global.xs(exp, level='Experiment', axis=1).sort_index(axis=1)
            df_AvgClusterPerformance = df_AvgClusterPerformance_global.xs(exp, level='Experiment', axis=1)
            fig_clusterPerformance.add_trace(go.Scatter(x=df_clusterPerformance.columns, y=df_clusterPerformance.loc["F1"], 
                                                    marker=dict(color=pio.templates["simple_white"].layout["colorway"][i]), name=exp))
            for  data_type in list_data_type:
                fig_clusterPerformance.add_trace(go.Scatter(x=[data_type], y=[df_AvgClusterPerformance[data_type].loc["F1"]],
                            mode="markers",
                            showlegend=False,
                            marker=dict(color=pio.templates["simple_white"].layout["colorway"][i])
                                ))
            fig_clusterPerformance.update_layout(template="simple_white", #showlegend=False, 
                            title="Cluster wise SVM analysis",
                            xaxis=go.layout.XAxis(linecolor="black",
                                                linewidth=1,
                                                mirror=True),
                            yaxis=go.layout.YAxis(linecolor="black",
                                                linewidth=1,
                                                title="F1 score", #- harmonic mean of recall and precision
                                                mirror=True),
                            )
        
        return fig_markerPredictionAccuracy, fig_clusterPerformance
        

        
    def __repr__(self):
        return str(self.__dict__)
        #return "This is a spatial dataset with {} lines.".format(len(self.df_original))

def svm_heatmap(df_SVM):
    """
    The misclassification matrix, generated by Perseus, will be displayed as a heatmap. 
    
    Args:
        self.df_SVM: dataframe, provided by Perseus, no index; 
                    Column names: e.g. "Predicted: ER", "Predicted: NPC"
                    Rows: e.g. "True: ER", "True: NPC"
    
    Returns:
        fig_SVMheatmap: heatmap of the misclassification matrix 
    """
    
    #df_SVM = self.df_SVM.copy()
    #if hasattr(df_SVM, "keys") == True: 
    try:
        df_SVM = pd.read_json(df_SVM["Misclassification Matrix"])
        df_SVM = df_SVM.set_index("T: True group")[::-1]
    
    except:
        df_SVM = df_SVM.set_index("T: True group")[::-1]
    
    y_axis_label = df_SVM.index
    x_axis_label = df_SVM.columns
    data_svm = df_SVM.values
    fig_SVMheatmap = go.Figure()
    
    fig_SVMheatmap.add_trace(go.Heatmap(
                             z=data_svm,
                             x = x_axis_label,
                             y = y_axis_label,
                             colorscale=[
                                 [0.0, "green"],
                                 [0.01, "white"],
                                 [1.0, "red"]
                             ],
                         ))
               
    return fig_SVMheatmap
    
    
def reframe_df_01_fromJson_for_Perseus(json_dict):
    """
    Make 0-1 normalized data from all experiments available for Perseus
    
    Args:
        json: dictionary, json file uploaded in manage dataset tab.
    
    Return: 
        df: 0-1 normlaized data (globally normalized), with Gene names, Protein IDs, Comaprtment as columns
            Pattern for Column data: Exp_Map_Fraction
    """
    for exp_name in json_dict.keys():
        for data_type in json_dict[exp_name].keys():                   
            if data_type == "0/1 normalized data" and exp_name == list(json_dict.keys())[0]:
                df_01_combined = pd.read_json(json_dict[exp_name][data_type])
                df_01_combined = df_01_combined.set_index(["Gene names", "Protein IDs", "Compartment"]).copy()
                df_01_combined.drop([col for col in df_01_combined.columns if not col.startswith("normalized profile")])
                df_01_combined.columns = pd.MultiIndex.from_tuples([el.split("?") for el in df_01_combined.columns], names=["Set", "Map", "Fraction"])
                df_01_combined.rename(columns = {"normalized profile":exp_name}, inplace=True)
    
            elif data_type == "0/1 normalized data" and exp_name != list(json_dict.keys())[0]:
                df_01_toadd = pd.read_json(json_dict[exp_name][data_type])
                df_01_toadd = df_01_toadd.set_index(["Gene names", "Protein IDs", "Compartment"]).copy()
                df_01_toadd.drop([col for col in df_01_toadd.columns if not col.startswith("normalized profile")])
                df_01_toadd.columns = pd.MultiIndex.from_tuples([el.split("?") for el in df_01_toadd.columns], names=["Set", "Map", "Fraction"])
                df_01_toadd.rename(columns = {"normalized profile":exp_name}, inplace=True)
                df_01_combined = pd.concat([df_01_combined, df_01_toadd], axis=1)
                                   
    df_01_combined.columns.names = ["Experiment", "Map", "Fraction"]
    df = df_01_combined.stack(["Experiment", "Map"]).dropna(axis=0)
    df = df.div(df.sum(axis=1), axis=0)
    index_ExpMap = df.index.get_level_values("Experiment")+"_"+df.index.get_level_values("Map")
    index_ExpMap.name = "Exp_Map"
    df.set_index(index_ExpMap, append=True, inplace=True)
    
    df.index = df.index.droplevel(["Map", "Experiment"])
    df = df.stack("Fraction").unstack(["Exp_Map", "Fraction"])
    df.columns = ["_".join(col) for col in df.columns.values]
    
    return df
   
def reannotate_genes_uniprot(proteingroups, source="uniprot", idmapping=[]):
    """
    Annotate protein groups with gene names from a specified source for uniprot annotations.
    
    Args:
        proteingroups: listlike, each protein group is split at ';' for multiple entries and isoforms are shortened at '-'
        source: str, one of 'uniprot', 'fasta', 'tsv'. If uniprot (default) annotations will be loaded live from uniprot.org
        idmapping: listlike if source=fasta, pandas DataFrame if source=tsv.
                   source=fasta: All pairs of protein ids and gene names (GN key) are retrieved from each element.
                                 Also works with a list of ";" separated fasta headers.
                   source=tsv: Must contain columns 'Entry' and 'Gene names  (primary )'.
    Returns:
        genes: list in same shape as proteingroups. Contains unique gene names, ';'-separated, by order of occurence in the protein group
    """
    
    protein_ids = []
    split_ids = [[i.split("-")[0] for i in el.split(";")] for el in proteingroups]
    protein_gene = {}
    
    if source == "uniprot":
        for ids in split_ids:
            for i in set(ids):
                if i not in protein_ids:
                    protein_ids.append(i)
        
        url = 'https://www.uniprot.org/uploadlists/'
        
        params = {
            'from': 'ACC+ID',
            'to': 'GENENAME',
            'format': 'tab',
            'query': "\t".join(protein_ids)
        }
        
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
            response = f.read()
        for l in response.decode('utf-8').split("\n")[1:-1]:
            l = l.split("\t")
            if l[0] in protein_gene.keys():
                protein_gene[l[0]] = protein_gene[l[0]]+"/"+l[1]
            else:
                protein_gene[l[0]] = l[1]
    
    elif source == "fasta_headers":
        for header in idmapping:
            p_g = re.findall("\|([^-\|]+)\|[^;]+GN=([^ ]+) ", header)
            for el in p_g:
                protein_gene[el[0]] = el[1]
    
    elif source == "tsv":
        for p,g in zip(idmapping["Entry"], idmapping["Gene names  (primary )"]):
            if type(g) == str:
                protein_gene[p] = g.replace("; ", "/")
    
    else:
        raise ValueError(f"Unkown source for gene reannotation: {source}")
    
    genes = list()
    for ids in split_ids:
        ids = [ids[el] for el in sorted([ids.index(i) for i in set(ids)])]
        genes.append(";".join([p if p not in protein_gene.keys() else protein_gene[p] for p in ids]))
    
    return genes


def split_ids(el, primary=";", isoform="-"):
    """
    This finds the primary canoncial protein ID in the protein group.
    If multiple canonical ids are present, it leaves the protein group unchanged.
    If no canonical ID is present it returns the first isoform ID.
    """
    p1 = el.split(primary)
    p2 = [p.split(isoform)[0] for p in p1]
    
    if len(set(p2)) > 1:
        return el
    elif p2[0] in p1:
        return p2[0]
    else:
        return p1[0]


def relabel_groups(exp, pgs, genes, sep=";"):
    
    assert len(exp) == len(pgs) == len(genes)
    
    src = np.array([
        exp,
        [el.split(sep) for el in pgs],
        [el.count(sep) for el in pgs],
        [el.split(sep) for el in genes],
        [el.count(sep) for el in genes],
        [i for i in range(len(genes))]
    ], dtype=object)
    
    src_single_gene = src[:,src[4] == 0]
    src_single_gene = src_single_gene[:,src_single_gene[2].argsort()[::-1]]
    src_multi_gene = src[:,src[4] != 0]
    src_multi_gene = src_multi_gene[:,np.lexsort((src_multi_gene[2], src_multi_gene[4]))[::-1]]
    
    out = np.empty([5,0])
    
    # Process multi gene entries
    while src_multi_gene.shape[1] > 0:
        group_max = src_multi_gene[3,0]
        pg_max = src_multi_gene[1,0]
        
        # find members from groups with multiple gene names
        members_multi = [0]
        for i,pg in enumerate(src_multi_gene[1]):
            if all(el in pg_max for el in pg) \
            and src_multi_gene[0,i] not in src_multi_gene[0, members_multi]:
                members_multi.append(i)
        group = src_multi_gene[:, members_multi]
        src_multi_gene = np.delete(src_multi_gene, members_multi, axis=1)
        
        # find members from groups with single gene names
        members_single = []
        for i,pg in enumerate(src_single_gene[1]):
            if all(el in pg_max for el in pg) \
            and src_single_gene[0,i] not in group[0] \
            and src_single_gene[0,i] not in src_single_gene[0, members_single]:
                members_single.append(i)
        group = np.append(group, src_single_gene[:, members_single], axis=1)
        src_single_gene = np.delete(src_single_gene, members_single, axis=1)
        
        # sort group by number of genes, then by number of protein groups and lastly 
        group[1] = [sep.join(el) for el in group[1]]
        group[3] = [sep.join(el) for el in group[3]]
        pg_counts = {k:v for k,v in zip(*np.unique(group[1], return_counts=True))}
        group = np.append(group, [[pg_counts[el] for el in group[1]]], axis=0)
        group = group[:,np.lexsort((group[6], group[2], group[4]))]
        #print(group_max, group)
        
        group_out = np.array([group[0],
                              np.repeat(group[1,-1], group.shape[1]),
                              np.repeat(group[3,-1], group.shape[1]),
                              np.repeat("multiple genes" if len(np.unique(group[3])) == 1 else "gene level conflict", group.shape[1]),
                              group[5]])
        out = np.append(out, group_out, axis=1)
    
    # Process leftover single gene entries
    single_out = np.array([src_single_gene[0],
                           [split_ids(sep.join(el)) for el in src_single_gene[1]],
                           [el[0] for el in src_single_gene[3]],
                           ["primary id" for el in src_single_gene[2]],
                           src_single_gene[5]])
    out = np.append(out, single_out, axis=1)
    out = out[0:4, out[-1].argsort()].T
    return pd.DataFrame(out, columns = [
        "Experiment", "Protein IDs", "Gene names",
        "merge type"
    ])


def align_datasets(df):
    """
    Align protein groups for profiling datasets from different experiments.
    
    Args:
    - df_01: pd.DataFrame
             index levels required:
              - Experiment
              - Gene names
              - Compartment -> This will be reassigned based on the aligned protein ID
              - Original Protein IDs -> This is the original ID collection for traceability
             column levels required:
              - Map
              - Fraction
    """
    
    assert all([el in df.index.names for el in
                ["Experiment","Gene names","Compartment","Original Protein IDs"]])
    assert all([el in df.columns.names for el in ["Map", "Fraction"]])
    
    df_01 = df.copy()
    n_exp = len(set(df_01.index.get_level_values("Experiment")))
    
    # 2. Retrieve data with simple full coverage alignment so it doesn't get torn during the more complicated matching
    df_01.set_index(pd.Index(
        [split_ids(el) for el in df_01.index.get_level_values("Original Protein IDs")],
        name="Protein IDs"), append=True, inplace=True)
    full_cov = df_01.groupby("Protein IDs").filter(lambda x: x.shape[0] >= n_exp)
    
    # 3. Relabel more complicated cases
    compl = df_01.drop(full_cov.index)
    new_labels = relabel_groups(compl.index.get_level_values("Experiment"),
                         compl.index.get_level_values("Original Protein IDs"),
                         [str(el) for el in compl.index.get_level_values("Gene names")])
    compl_relabeled = compl.copy()
    for i in new_labels.columns:
        if i in compl_relabeled.index.names:
            compl_relabeled.index = compl_relabeled.index.droplevel(i)
        compl_relabeled.set_index(pd.Index(new_labels[i], name=i), append=True, inplace=True)
    full_cov.set_index(pd.Index(np.repeat("primary id", len(full_cov)), name="merge type"), append=True, inplace=True)
    compl_relabeled.index = compl_relabeled.index.reorder_levels(full_cov.index.names)
    
    # 4. Get index mapping
    index_mapping = pd.concat([
        full_cov.index.to_frame(index=False),
        compl_relabeled.index.to_frame(index=False)], axis=0)
    index_mapping.set_index(["Experiment", "Gene names", "Protein IDs", "merge type"], inplace=True)
    index_mapping = index_mapping.unstack("Experiment")
    
    # 5. Rejoin datasets and relabel Compartments
    df_01_aligned = pd.concat([full_cov, compl_relabeled], axis=0)
    df_01_aligned.index = df_01_aligned.index.droplevel(["Original Protein IDs", "Compartment"])
    df_01_aligned = df_01_aligned.unstack("Experiment") # rejoin index
    c = []
    for el in df_01_aligned.index.get_level_values("Protein IDs"):
        cs = [el for el in index_mapping.xs(el, level="Protein IDs", axis=0).Compartment.iloc[0,:].values if type(el) == str]
        if len(set(cs)) == 1:
            c.append(cs[0])
        else:
            c.append("undefined")
    df_01_aligned.set_index(pd.Index(c, name="Compartment"), append=True, inplace=True)
    df_01_aligned = df_01_aligned.stack(["Experiment", "Map"])
    
    return df_01_aligned, index_mapping