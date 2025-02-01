import natsort
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import datashader as ds
import re
from io import BytesIO, StringIO
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise as pw
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn import svm
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.covariance import MinCovDet
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests
from scipy.stats import combine_pvalues
import statistics
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from PIL import Image
from upsetplot import from_memberships
from upsetplot import plot as upplot
import pkg_resources
from scipy.stats import zscore
import urllib.parse
import urllib.request
import bisect
import inspect
import warnings
from itertools import cycle, combinations

from domaps.constants import (
    DataFrameStrings,
    DefaultAcquisitionSettings,
    SettingStrings,
)

from domaps.__init__ import __version__ as version


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

    # fmt: off
    css_color = ["#b2df8a","#6a3d9a","#e31a1c","#b15928","#fdbf6f","#ff7f00","#cab2d6","#fb9a99","#1f78b4","#ffff99","#a6cee3","#33a02c","blue","orange","goldenrod","lightcoral","magenta","brown","lightpink","red","turquoise","khaki","darkgoldenrod","darkturquoise","darkviolet","greenyellow","darksalmon","hotpink","indianred","indigo","darkolivegreen","coral","aqua","beige","bisque","black","blanchedalmond","blueviolet","burlywood","cadetblue","yellowgreen","chartreuse","chocolate","cornflowerblue","cornsilk","darkblue","darkcyan","darkgray","darkgrey","darkgreen","darkkhaki","darkmagenta","darkorange","darkorchid","darkred","darkseagreen","darkslateblue","snow","springgreen","darkslategrey","mediumpurple","oldlace","olive","lightseagreen","deeppink","deepskyblue","dimgray","dimgrey","dodgerblue","firebrick","floralwhite","forestgreen","fuchsia","gainsboro","ghostwhite","gold","gray","ivory","lavenderblush","lawngreen","lemonchiffon","lightblue","lightcyan","fuchsia","gainsboro","ghostwhite","gold","gray","ivory","lavenderblush","lawngreen","lemonchiffon","lightblue","lightcyan","lightgoldenrodyellow","lightgray","lightgrey","lightgreen","lightsalmon","lightskyblue","lightslategray","lightslategrey","lightsteelblue","lightyellow","lime","limegreen","linen","maroon","mediumaquamarine","mediumblue","mediumseagreen","mediumslateblue","mediumspringgreen","mediumturquoise","mediumvioletred","midnightblue","mintcream","mistyrose","moccasin","olivedrab","orangered","orchid","palegoldenrod","palegreen","paleturquoise","palevioletred","papayawhip","peachpuff","peru","pink","plum","powderblue","rosybrown","royalblue","saddlebrown","salmon","sandybrown","seagreen","seashell","sienna","silver","skyblue","slateblue","steelblue","teal","thistle","tomato","violet","wheat","white","whitesmoke","slategray","slategrey","aquamarine","azure","crimson","cyan","darkslategray","grey","mediumorchid","navajowhite","navy"]
    # fmt: on

    def __init__(
        self,
        filename,
        expname,
        acquisition="LFQ",
        source="MaxQuant_proteins_pivot",  # future
        orientation="pivot",  # future
        comment="",
        name_pattern="e.g.:.* (?P<cond>.*)_(?P<rep>.*)_(?P<frac>.*)",
        reannotate_genes=False,
        reannotation_source={"source": False, "idmapping": None},
        RatioCount=2,
        RatioVariability=30,
        average_MSMS_counts=2,  # future
        consecutive=4,  # future
        organism="Homo sapiens - Uniprot",  # future
        organelles="Homo sapiens - Uniprot",  # future
        complexes="Homo sapiens - Uniprot",  # future
        fractions=[],
        **kwargs,
    ):
        self.filename = filename
        self.expname = expname
        self.acquisition = acquisition
        self.name_pattern = name_pattern
        self.comment = comment
        self.fractions = fractions

        self.source = source
        self.orientation = orientation

        self.map_names = []
        self.df_01_stacked, self.df_log_stacked = pd.DataFrame(), pd.DataFrame()

        if SettingStrings.REANNOTATE in kwargs.keys():
            reannotate_genes = kwargs[SettingStrings.REANNOTATE]

        self.reannotate_genes = reannotate_genes
        if reannotate_genes:
            self.reannotation_source = reannotation_source

        if acquisition == SettingStrings.SILAC_ACQUISITION:
            self.RatioCount = RatioCount
            self.RatioVariability = RatioVariability
            self.mainset = DataFrameStrings.RATIO

        elif acquisition != SettingStrings.SILAC_ACQUISITION:
            self.average_MSMS_counts = average_MSMS_counts
            self.consecutive = consecutive
            self.mainset = DefaultAcquisitionSettings[acquisition][SettingStrings.SETS][
                0
            ]

        self.transformed_mainset = self.mainset

        self.organism = organism
        if complexes + ".csv" in pkg_resources.resource_listdir(
            __name__, "annotations/complexes"
        ):
            marker_table = pd.read_csv(
                pkg_resources.resource_stream(
                    __name__, "annotations/complexes/{}.csv".format(complexes)
                )
            )
        else:
            marker_table = pd.read_csv(StringIO(complexes))
        self.markerproteins = {
            k: v.replace(" ", "").split(",")
            for k, v in zip(
                marker_table[DataFrameStrings.CLUSTER],
                marker_table["Members - Protein IDs"],
            )
        }
        if organelles + ".csv" in pkg_resources.resource_listdir(
            __name__, "annotations/organellemarkers"
        ):
            df_organellarMarkerSet = pd.read_csv(
                pkg_resources.resource_stream(
                    __name__,
                    "annotations/organellemarkers/{}.csv".format(organelles),
                ),
                usecols=lambda x: bool(re.match("Compartment|Protein ID", x)),
            ).drop_duplicates()
        else:
            df_organellarMarkerSet = pd.read_csv(
                StringIO(organelles),
                usecols=lambda x: bool(re.match("Compartment|Protein ID", x)),
            ).drop_duplicates()
        self.df_organellarMarkerSet = df_organellarMarkerSet

        self.kwargs = kwargs

        self.imported_columns = generate_usecols_regex(
            sets=kwargs[SettingStrings.SETS],
            original_protein_ids=kwargs[SettingStrings.ORIGINAL_PROTEIN_IDS],
            genes=kwargs[SettingStrings.GENES],
            samples=kwargs.get(SettingStrings.SAMPLES, None),
            column_filters=kwargs.get(SettingStrings.COLUMN_FILTERS, None),
            annotation_columns=kwargs.get(SettingStrings.ANNOTATION_COLUMNS, None),
        )

        self.analysed_datasets_dict = {}
        self.analysis_summary_dict = {}

    def from_settings(settings):
        settings[SettingStrings.REANNOTATION_SOURCE] = parse_reannotation_source(
            settings[SettingStrings.REANNOTATE],
            settings[SettingStrings.REANNOTATION_SOURCE],
        )
        ds = SpatialDataSet(**settings)
        return ds

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
        self.process_input()
        self.quantity_profiles_proteinGroups()
        setprogress("PCA ...")
        self.perform_pca()
        setprogress("Calculating Manhattan Distance ...")
        self.calc_biological_precision()
        setprogress("Calculating replicate correlations ... ")
        self.calculate_pairwise_correlations()
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

        if (
            filename.endswith("xls")
            or filename.endswith("txt")
            or filename.endswith("tsv")
        ):
            self.df_original = pd.read_csv(
                content,
                sep="\t",
                comment="#",
                usecols=lambda x: bool(re.match(self.imported_columns, x)),
                low_memory=True,
            )
        else:  # assuming csv file
            self.df_original = pd.read_csv(
                content,
                sep=",",
                comment="#",
                usecols=lambda x: bool(re.match(self.imported_columns, x)),
                low_memory=True,
            )
        assert self.df_original.shape[0] > 10 and self.df_original.shape[1] > 2

        self.filename = filename

        return self.df_original

    def process_input(self):
        """ """

        processing_steps = []

        orientation = self.orientation

        ## Format data based on input columns
        if self.orientation == "long":
            df_index = format_data_long(
                df=self.df_original,
                name_pattern=self.name_pattern,
                index_cols=self.kwargs[SettingStrings.ANNOTATION_COLUMNS]
                + [k for k in self.kwargs["column_filters"].keys()],
                **{
                    k: v
                    for k, v in self.kwargs.items()
                    if k in inspect.getfullargspec(format_data_long).args
                },
            )
            processing_steps.append("Formatting long data")
        elif self.orientation == "pivot":
            df_index = format_data_pivot(
                df=self.df_original,
                name_pattern=self.name_pattern,
                index_cols=self.kwargs[SettingStrings.ANNOTATION_COLUMNS]
                + [k for k in self.kwargs["column_filters"].keys()],
                **{
                    k: v
                    for k, v in self.kwargs.items()
                    if k in inspect.getfullargspec(format_data_pivot).args
                },
            )
            processing_steps.append("Formatting pivot data")
        else:
            raise ValueError("Unknown input format")

        self.map_names = df_index.columns.get_level_values(
            DataFrameStrings.MAP
        ).unique()
        if "<cond>" in self.name_pattern:
            self.conditions = list(set([el.split("_")[0] for el in self.map_names]))
        else:
            self.conditions = [""]

        ## Reannotate genes
        if self.reannotate_genes:
            df_index.reset_index(DataFrameStrings.GENE_NAMES, inplace=True, drop=True)
            df_index.set_index(
                pd.Index(
                    name=DataFrameStrings.GENE_NAMES,
                    data=reannotate_genes_uniprot(
                        df_index.index.get_level_values(DataFrameStrings.PROTEIN_IDS),
                        **self.reannotation_source,
                    ),
                ),
                inplace=True,
                append=True,
            )
            processing_steps.append("Reannotating genes")
        else:
            genes = [
                g if str(g) != "nan" else p
                for g, p in zip(
                    df_index.index.get_level_values(DataFrameStrings.GENE_NAMES),
                    df_index.index.get_level_values(DataFrameStrings.PROTEIN_IDS),
                )
            ]
            df_index.reset_index(DataFrameStrings.GENE_NAMES, inplace=True, drop=True)
            df_index.set_index(
                pd.Index(name=DataFrameStrings.GENE_NAMES, data=genes),
                inplace=True,
                append=True,
            )
            processing_steps.append("Filling missing gene names with Protein IDs")

        self.df_index = df_index

        ## Annotate marker proteins
        df_index.columns = df_index.columns.values
        df_index = df_index.join(
            self.df_organellarMarkerSet.set_index("Protein ID"),
            how="left",
            on=DataFrameStrings.PROTEIN_IDS,
        ).set_index(DataFrameStrings.COMPARTMENT, append=True)
        df_index.columns = pd.MultiIndex.from_tuples(
            df_index.columns,
            names=[
                DataFrameStrings.SET,
                DataFrameStrings.MAP,
                DataFrameStrings.FRACTION,
            ],
        )
        df_index.rename(
            index={np.nan: "undefined"},
            level=DataFrameStrings.COMPARTMENT,
            inplace=True,
        )
        processing_steps.append("Annotating marker proteins")

        ## Filter data quality as applicable
        df_filtered = df_index.copy()
        for k, v in self.kwargs["column_filters"].items():
            df_filtered = filter_singlecolumn_keep(
                df_filtered, column=k, operator=v[0], value=v[1]
            )
            processing_steps.append(f"Filtering {k}")
        if len(self.kwargs[SettingStrings.QUALITY_FILTER]) == 0:
            processing_steps.append("No further quality filters applied.")
        elif all(
            [
                el
                in [
                    SettingStrings.FILTER_SILAC_COUNTVAR,
                    SettingStrings.FILTER_MSMS_COUNT,
                    SettingStrings.FILTER_CONSECUTIVE,
                ]
                for el in self.kwargs[SettingStrings.QUALITY_FILTER]
            ]
        ):
            if (
                SettingStrings.FILTER_SILAC_COUNTVAR
                in self.kwargs[SettingStrings.QUALITY_FILTER]
            ):
                df_filtered = filter_SILAC_countvar(
                    df_filtered, self.RatioCount, self.RatioVariability
                )
                processing_steps.append("Filtering SILAC countvar")
            if (
                SettingStrings.FILTER_MSMS_COUNT
                in self.kwargs[SettingStrings.QUALITY_FILTER]
            ):
                df_filtered = filter_msms_count(df_filtered, self.average_MSMS_counts)
                processing_steps.append("Filtering MSMS count")
            if (
                SettingStrings.FILTER_CONSECUTIVE
                in self.kwargs[SettingStrings.QUALITY_FILTER]
            ):
                df_filtered = filter_consecutive(
                    df_filtered,
                    self.consecutive,
                    sets={DataFrameStrings.ABUNDANCE: self.mainset},
                )
                processing_steps.append("Filtering consecutive values")
        else:
            raise ValueError("Unknown quality filter")

        ## Run data tranformations
        if self.kwargs[SettingStrings.INPUT_LOGGED]:
            df_transformed = transform_data(
                df_filtered[[self.mainset]],
                invert=False,
                unlog=self.kwargs[SettingStrings.INPUT_LOGGED],
                log=False,
            )
            processing_steps.append("Unlogging input data")
        else:
            df_transformed = df_filtered[[self.mainset]].copy()
        if self.kwargs[SettingStrings.INPUT_SAMPLENORMALIZATION]:
            df_transformed = normalize_samples(
                df_transformed,
                method=self.kwargs[SettingStrings.INPUT_SAMPLENORMALIZATION],
            )
            df_transformed.rename(
                {self.transformed_mainset: "Normalized " + self.transformed_mainset},
                axis=1,
                inplace=True,
            )
            self.transformed_mainset = "Normalized " + self.transformed_mainset
            processing_steps.append(
                f"Normalizing samples by {self.kwargs[SettingStrings.INPUT_SAMPLENORMALIZATION]}"
            )
        if self.kwargs[SettingStrings.INPUT_INVERT]:
            df_transformed = transform_data(
                df_transformed,
                invert=self.kwargs[SettingStrings.INPUT_INVERT],
                unlog=False,
                log=False,
            )
            df_transformed.rename(
                {self.transformed_mainset: "Inverted " + self.transformed_mainset},
                axis=1,
                inplace=True,
            )
            self.transformed_mainset = "Inverted " + self.transformed_mainset
            processing_steps.append("Inverting data")

        # if len(self.kwargs["yields"]) != 0:
        #    df_transformed = domaps.weigh_yields(df_transformed)

        df_transformed = df_filtered.drop(self.mainset, axis=1).join(df_transformed)

        if self.acquisition == "SILAC":
            # complete profiles
            df_transformed = filter_consecutive(
                df_transformed,
                len(self.fractions),
                sets={DataFrameStrings.ABUNDANCE: self.transformed_mainset},
            )
            processing_steps.append("Filtering for complete ratio profiles")
        if self.acquisition == "LFQ":
            df_transformed = (
                df_transformed.stack(DataFrameStrings.MAP)
                .replace(np.NaN, 0)
                .unstack(DataFrameStrings.MAP)
            )

        ## Normalize profiles
        df_01_stacked = normalize_sum1(
            df_transformed, sets={DataFrameStrings.ABUNDANCE: self.transformed_mainset}
        ).stack(DataFrameStrings.FRACTION)
        processing_steps.append("Normalizing to % profiles")

        ## Store dataset versions
        self.df_index = df_index
        self.df_filtered = df_filtered
        df_log_stacked = df_01_stacked.drop(DataFrameStrings.NORMALIZED_PROFILE, axis=1)
        df_log_stacked = df_log_stacked.drop(self.transformed_mainset, axis=1).join(
            df_log_stacked[[self.transformed_mainset]]
            .replace(0, np.nan)
            .transform(np.log2)
            .rename({self.transformed_mainset: DataFrameStrings.LOG_PROFILE}, axis=1)
        )
        self.df_log_stacked = df_log_stacked
        self.df_01_stacked = df_01_stacked

        # format and reduce 0-1 normalized data for comparison with other experiments
        df_01_comparison = df_01_stacked[[DataFrameStrings.NORMALIZED_PROFILE]].copy()
        df_01_comparison = df_01_comparison.unstack(
            [DataFrameStrings.MAP, DataFrameStrings.FRACTION]
        )
        df_01_comparison.columns = [
            "?".join(el) for el in df_01_comparison.columns.values
        ]
        df_01_comparison = df_01_comparison.reset_index(
            level=[
                DataFrameStrings.GENE_NAMES,
                DataFrameStrings.ORIGINAL_PROTEIN_IDS,
                DataFrameStrings.PROTEIN_IDS,
                DataFrameStrings.COMPARTMENT,
            ]
        ).reset_index(drop=True)

        # poopulate analysis summary dictionary with (meta)data
        self.analysis_summary_dict["0/1 normalized data"] = df_01_comparison.to_json()
        self.analysis_summary_dict[SettingStrings.ORGANELLES] = (
            self.df_organellarMarkerSet.to_json()
        )
        self.analysis_summary_dict[SettingStrings.COMPLEXES] = self.markerproteins
        # self.analysis_summary_dict["changes in shape after filtering"] = shape_dict.copy()
        # add more settings here
        analysis_parameters = {
            SettingStrings.ACQUISITION: self.acquisition,
            SettingStrings.FILENAME: self.filename,
            SettingStrings.COMMENT: self.comment,
            SettingStrings.ORGANISM: self.organism,
            "processing steps": "\n".join(processing_steps),
            "legacy": False,  # TODO: remove this with 2.0 (only required for reading new files with versions < 1.0.4)
            "domaps version": version,
        }
        self.analysis_summary_dict["Analysis parameters"] = analysis_parameters.copy()

        return df_01_stacked

    def plot_log_data(self):
        """

        Args:
            self.df_log_stacked


        Returns:
            log_histogram: Histogram of log transformed data

        """

        log_histogram = px.histogram(
            self.df_log_stacked.reset_index().sort_values(
                [DataFrameStrings.MAP, DataFrameStrings.FRACTION], key=natsort_list_keys
            ),
            x=DataFrameStrings.LOG_PROFILE,
            facet_col=DataFrameStrings.FRACTION,
            facet_row=DataFrameStrings.MAP,
            template="simple_white",
            labels={
                DataFrameStrings.LOG_PROFILE: "log tranformed data ({})".format(
                    "LFQ intenisty" if self.acquisition != "SILAC - MQ" else "Ratio H/L"
                )
            },
        )

        log_histogram.for_each_xaxis(lambda axis: axis.update(title={"text": ""}))
        log_histogram.for_each_yaxis(lambda axis: axis.update(title={"text": ""}))
        log_histogram.add_annotation(
            x=0.5,
            y=0,
            yshift=-50,
            xref="paper",
            showarrow=False,
            yref="paper",
            text="log2(LFQ intensity)",
        )
        log_histogram.add_annotation(
            x=0,
            y=0.5,
            textangle=270,
            xref="paper",
            showarrow=False,
            yref="paper",
            xshift=-50,
            text="count",
        )
        log_histogram.for_each_annotation(
            lambda a: a.update(text=a.text.split("=")[-1])
        )

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

        df_01_stacked = self.df_01_stacked[DataFrameStrings.NORMALIZED_PROFILE].replace(
            0, np.nan
        )
        df_index = self.df_index[self.mainset]

        # unfiltered
        npg_t = df_index.shape[0]
        df_index_MapStacked = df_index.stack(DataFrameStrings.MAP)
        npr_t = df_index_MapStacked.shape[0] / len(self.map_names)
        npr_t_dc = 1 - df_index_MapStacked.isna().sum().sum() / np.prod(
            df_index_MapStacked.shape
        )

        # filtered
        npgf_t = df_01_stacked.unstack(
            [DataFrameStrings.MAP, DataFrameStrings.FRACTION]
        ).shape[0]
        df_01_MapStacked = df_01_stacked.unstack(DataFrameStrings.FRACTION)
        nprf_t = df_01_MapStacked.shape[0] / len(self.map_names)
        nprf_t_dc = 1 - df_01_MapStacked.isna().sum().sum() / np.prod(
            df_01_MapStacked.shape
        )

        # unfiltered intersection
        try:
            df_index_intersection = df_index_MapStacked.groupby(
                level="Sequence"
            ).filter(lambda x: len(x) == len(self.map_names))
        except:
            df_index_intersection = df_index_MapStacked.groupby(
                level=DataFrameStrings.PROTEIN_IDS
            ).filter(lambda x: len(x) == len(self.map_names))
        npr_i = df_index_intersection.shape[0] / len(self.map_names)
        npr_i_dc = 1 - df_index_intersection.isna().sum().sum() / np.prod(
            df_index_intersection.shape
        )
        npg_i = df_index_intersection.unstack(DataFrameStrings.MAP).shape[0]

        # filtered intersection
        try:
            df_01_intersection = df_01_MapStacked.groupby(level="Sequence").filter(
                lambda x: len(x) == len(self.map_names)
            )
        except:
            df_01_intersection = df_01_MapStacked.groupby(
                level=DataFrameStrings.PROTEIN_IDS
            ).filter(lambda x: len(x) == len(self.map_names))
        nprf_i = df_01_intersection.shape[0] / len(self.map_names)
        nprf_i_dc = 1 - df_01_intersection.isna().sum().sum() / np.prod(
            df_01_intersection.shape
        )
        npgf_i = df_01_intersection.unstack(DataFrameStrings.MAP).shape[0]

        # summarize in dataframe and save to attribute
        df_quantity_pr_pg = pd.DataFrame(
            {
                "filtering": pd.Series(
                    [
                        "before filtering",
                        "before filtering",
                        "after filtering",
                        "after filtering",
                    ],
                    dtype=np.dtype("O"),
                ),
                "type": pd.Series(
                    ["total", "intersection", "total", "intersection"],
                    dtype=np.dtype("O"),
                ),
                "number of protein groups": pd.Series(
                    [npg_t, npg_i, npgf_t, npgf_i], dtype=np.dtype("float")
                ),
                "number of profiles": pd.Series(
                    [npr_t, npr_i, nprf_t, nprf_i], dtype=np.dtype("float")
                ),
                "data completeness of profiles": pd.Series(
                    [npr_t_dc, npr_i_dc, nprf_t_dc, nprf_i_dc], dtype=np.dtype("float")
                ),
            }
        )

        self.df_quantity_pr_pg = df_quantity_pr_pg.reset_index()
        self.analysis_summary_dict["quantity: profiles/protein groups"] = (
            self.df_quantity_pr_pg.to_json()
        )

        # additional depth assessment per fraction
        dict_npgf = {}
        dict_npg = {}
        list_npg_dc = []
        list_npgf_dc = []

        for df_intersection in [df_index_intersection, df_01_intersection]:
            for fraction in self.fractions:
                df_intersection_frac = df_intersection[fraction]
                npgF_f_dc = 1 - df_intersection_frac.isna().sum() / len(
                    df_intersection_frac
                )
                npgF_f = (
                    df_intersection_frac.unstack(DataFrameStrings.MAP)
                    .isnull()
                    .sum(axis=1)
                    .value_counts()
                )
                if fraction not in dict_npg.keys():
                    dict_npg[fraction] = npgF_f
                    list_npg_dc.append(npgF_f_dc)
                else:
                    dict_npgf[fraction] = npgF_f
                    list_npgf_dc.append(npgF_f_dc)

        df_npg = pd.DataFrame(dict_npg)
        df_npg.index.name = "Protein Groups present in:"
        df_npg.rename_axis(DataFrameStrings.FRACTION, axis=1, inplace=True)
        df_npg = df_npg.stack(DataFrameStrings.FRACTION).reset_index()
        df_npg = df_npg.rename({0: "Protein Groups"}, axis=1)
        df_npg.sort_values(
            [DataFrameStrings.FRACTION, "Protein Groups present in:"],
            inplace=True,
            key=natsort_list_keys,
        )

        df_npgf = pd.DataFrame(dict_npgf)
        df_npgf.index.name = "Protein Groups present in:"
        df_npgf.rename_axis(DataFrameStrings.FRACTION, axis=1, inplace=True)
        df_npgf = df_npgf.stack(DataFrameStrings.FRACTION).reset_index()
        df_npgf = df_npgf.rename({0: "Protein Groups"}, axis=1)
        df_npgf.sort_values(
            [DataFrameStrings.FRACTION, "Protein Groups present in:"],
            inplace=True,
            key=natsort_list_keys,
        )

        max_df_npg = df_npg["Protein Groups present in:"].max()
        min_df_npg = df_npg["Protein Groups present in:"].min()
        rename_numOFnans = {}
        for x, y in zip(range(max_df_npg, min_df_npg - 1, -1), range(max_df_npg + 1)):
            if y == 1:
                rename_numOFnans[x] = "{} Map".format(y)
            elif y == 0:
                rename_numOFnans[x] = "PG not identified"
            else:
                rename_numOFnans[x] = "{} Maps".format(y)
        for keys in rename_numOFnans.keys():
            df_npg.loc[
                df_npg["Protein Groups present in:"] == keys,
                "Protein Groups present in:",
            ] = rename_numOFnans[keys]
            df_npgf.loc[
                df_npgf["Protein Groups present in:"] == keys,
                "Protein Groups present in:",
            ] = rename_numOFnans[keys]

        # summarize in dataframe and save to attributes
        self.df_npg_dc = pd.DataFrame(
            {
                DataFrameStrings.FRACTION: pd.Series(self.fractions),
                "Data completeness before filtering": pd.Series(list_npg_dc),
                "Data completeness after filtering": pd.Series(list_npgf_dc),
            }
        )

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

        layout = go.Layout(
            barmode="overlay",
            xaxis_tickangle=90,
            autosize=False,
            width=300,
            height=500,
            xaxis=go.layout.XAxis(
                linecolor="black",
                linewidth=1,
                # title="Map",
                mirror=True,
            ),
            yaxis=go.layout.YAxis(linecolor="black", linewidth=1, mirror=True),
            template="simple_white",
        )

        fig_npg = go.Figure()
        for t in df_quantity_pr_pg["type"].unique():
            plot_df = df_quantity_pr_pg[df_quantity_pr_pg["type"] == t]
            fig_npg.add_trace(
                go.Bar(
                    x=plot_df["filtering"],
                    y=plot_df["number of protein groups"],
                    name=t,
                )
            )
        fig_npg.update_layout(
            layout,
            title="Number of Protein Groups",
            yaxis=go.layout.YAxis(title="Protein Groups"),
        )
        fig_npg.update_layout(**calc_width_categorical(fig_npg))

        fig_npr = go.Figure()
        for t in df_quantity_pr_pg["type"].unique():
            plot_df = df_quantity_pr_pg[df_quantity_pr_pg["type"] == t]
            fig_npr.add_trace(
                go.Bar(x=plot_df["filtering"], y=plot_df["number of profiles"], name=t)
            )
        fig_npr.update_layout(layout, title="Number of Profiles")
        fig_npr.update_layout(**calc_width_categorical(fig_npr))

        df_quantity_pr_pg = df_quantity_pr_pg.sort_values("filtering")
        fig_npr_dc = go.Figure()
        for t in df_quantity_pr_pg["filtering"].unique():
            plot_df = df_quantity_pr_pg[df_quantity_pr_pg["filtering"] == t]
            fig_npr_dc.add_trace(
                go.Bar(
                    x=plot_df["type"],
                    y=plot_df["data completeness of profiles"],
                    name=t,
                )
            )
        fig_npr_dc.update_layout(
            layout, title="Coverage", yaxis=go.layout.YAxis(title="Data completness")
        )
        fig_npr_dc.update_layout(**calc_width_categorical(fig_npr_dc))

        fig_npg_F = px.bar(
            self.df_npg,
            x=DataFrameStrings.FRACTION,
            y="Protein Groups",
            color="Protein Groups present in:",
            template="simple_white",
            title="Protein groups per fraction - before filtering",
            width=500,
        )
        fig_npg_F.update_layout(**calc_width_categorical(fig_npg_F))

        fig_npgf_F = px.bar(
            self.df_npgf,
            x=DataFrameStrings.FRACTION,
            y="Protein Groups",
            color="Protein Groups present in:",
            template="simple_white",
            title="Protein groups per fraction - after filtering",
            width=500,
        )
        fig_npgf_F.update_layout(**calc_width_categorical(fig_npgf_F))

        fig_npg_F_dc = go.Figure()
        for data_type in [
            "Data completeness after filtering",
            "Data completeness before filtering",
        ]:
            fig_npg_F_dc.add_trace(
                go.Bar(
                    x=self.df_npg_dc[DataFrameStrings.FRACTION],
                    y=self.df_npg_dc[data_type],
                    name=data_type,
                )
            )
        fig_npg_F_dc.update_layout(
            layout,
            barmode="overlay",
            title="Data completeness per fraction",
            yaxis=go.layout.YAxis(title=""),
            height=450,
            width=600,
        )
        fig_npg_F_dc.update_layout(**calc_width_categorical(fig_npg_F_dc))

        return fig_npg, fig_npr, fig_npr_dc, fig_npg_F, fig_npgf_F, fig_npg_F_dc

    def perform_pca(self, n=3):
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
        df_01orlog_fracunstacked = (
            self.df_01_stacked[DataFrameStrings.NORMALIZED_PROFILE]
            .unstack(DataFrameStrings.FRACTION)
            .dropna()
        )
        df_01orlog_MapFracUnstacked = (
            self.df_01_stacked[DataFrameStrings.NORMALIZED_PROFILE]
            .unstack([DataFrameStrings.FRACTION, DataFrameStrings.MAP])
            .dropna()
        )

        pca = PCA(n_components=n)

        # df_pca: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3"
        df_pca = pd.DataFrame(
            pca.fit_transform(
                df_01orlog_fracunstacked.apply(zscore, axis=0).replace(np.nan, 0)
            )
        )
        df_pca.columns = [f"PC{el + 1}" for el in range(n)]
        df_pca.index = df_01orlog_fracunstacked.index
        self.df_pca = df_pca.sort_index(
            level=[DataFrameStrings.PROTEIN_IDS, DataFrameStrings.COMPARTMENT]
        )
        self.df_pca_loadings = pd.DataFrame(
            pca.components_,
            columns=df_01orlog_fracunstacked.columns,
            index=df_pca.columns,
        ).T.reset_index()
        self.df_pca_var = (
            pd.DataFrame(pca.explained_variance_ratio_, index=df_pca.columns)
            .reset_index()
            .rename({"index": "Component", 0: "variance explained"}, axis=1)
        )

        # df_pca: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3"
        df_pca_combined = pd.DataFrame(
            pca.fit_transform(
                df_01orlog_MapFracUnstacked.apply(zscore, axis=0).replace(np.nan, 0)
            )
        )
        df_pca_combined.columns = [f"PC{el + 1}" for el in range(n)]
        df_pca_combined.index = df_01orlog_MapFracUnstacked.index
        self.df_pca_combined = df_pca_combined.sort_index(
            level=[DataFrameStrings.PROTEIN_IDS, DataFrameStrings.COMPARTMENT]
        )
        self.df_pca_combined_loadings = pd.DataFrame(
            pca.components_,
            columns=df_01orlog_MapFracUnstacked.columns,
            index=df_pca_combined.columns,
        ).T.reset_index()
        self.df_pca_combined_var = (
            pd.DataFrame(pca.explained_variance_ratio_, index=df_pca_combined.columns)
            .reset_index()
            .rename({"index": "Component", 0: "variance explained"}, axis=1)
        )

        map_names = self.map_names
        df_pca_all_marker_cluster_maps = pd.DataFrame()
        df_pca_filtered = df_pca.unstack(DataFrameStrings.MAP)
        for clusters in markerproteins:
            for marker in markerproteins[clusters]:
                try:
                    plot_try_pca = df_pca_filtered.xs(
                        marker, level=DataFrameStrings.PROTEIN_IDS, drop_level=False
                    ).set_index(
                        pd.Index([clusters], name=DataFrameStrings.CLUSTER), append=True
                    )
                except KeyError:
                    continue
                df_pca_all_marker_cluster_maps = pd.concat(
                    [df_pca_all_marker_cluster_maps, plot_try_pca],
                    axis=0,
                )
        if len(df_pca_all_marker_cluster_maps) == 0:
            df_pca_all_marker_cluster_maps = df_pca_filtered.stack(DataFrameStrings.MAP)
        else:
            df_pca_all_marker_cluster_maps = df_pca_all_marker_cluster_maps.stack(
                DataFrameStrings.MAP
            )
        self.df_pca_all_marker_cluster_maps = df_pca_all_marker_cluster_maps.sort_index(
            level=[DataFrameStrings.PROTEIN_IDS, DataFrameStrings.COMPARTMENT]
        )

    def plot_global_pca(
        self,
        map_of_interest="Map1",
        cluster_of_interest="Proteasome",
        x_PCA="PC1",
        y_PCA="PC3",
        collapse_maps=False,
    ):
        """ "
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
            df_pca = (
                self.df_pca.unstack(DataFrameStrings.MAP)
                .swaplevel(0, 1, axis=1)[map_of_interest]
                .reset_index()
            )
        else:
            df_pca = self.df_pca_combined.reset_index()

        for i in self.markerproteins[cluster_of_interest]:
            df_pca.loc[
                df_pca[DataFrameStrings.PROTEIN_IDS] == i, DataFrameStrings.COMPARTMENT
            ] = "Selection"

        compartments = self.df_organellarMarkerSet[
            DataFrameStrings.COMPARTMENT
        ].unique()
        compartment_color = dict(zip(compartments, self.css_color))
        compartment_color["Selection"] = "black"
        compartment_color["undefined"] = "lightgrey"

        fig_global_pca = px.scatter(
            data_frame=df_pca,
            x=x_PCA,
            y=y_PCA,
            color=DataFrameStrings.COMPARTMENT,
            color_discrete_map=compartment_color,
            title="Protein subcellular localization by PCA for {}".format(
                map_of_interest
            )
            if collapse_maps == False
            else "Protein subcellular localization by PCA of combined maps",
            hover_data=[
                DataFrameStrings.PROTEIN_IDS,
                DataFrameStrings.GENE_NAMES,
                DataFrameStrings.COMPARTMENT,
            ],
            template="simple_white",
            opacity=0.9,
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
                        plot_try_pca = df_pca_all_marker_cluster_maps.xs(
                            (marker, maps),
                            level=[DataFrameStrings.PROTEIN_IDS, DataFrameStrings.MAP],
                            drop_level=False,
                        )
                    except KeyError:
                        continue
                    df_setofproteins_PCA = pd.concat(
                        [df_setofproteins_PCA, plot_try_pca]
                    )

                df_setofproteins_PCA.reset_index(inplace=True)
                if maps == map_names[0]:
                    pca_figure = go.Figure(
                        data=[
                            go.Scatter3d(
                                x=df_setofproteins_PCA.PC1,
                                y=df_setofproteins_PCA.PC2,
                                z=df_setofproteins_PCA.PC3,
                                hovertext=df_setofproteins_PCA[
                                    DataFrameStrings.GENE_NAMES
                                ],
                                mode="markers",
                                name=maps,
                            )
                        ]
                    )
                else:
                    pca_figure.add_trace(
                        go.Scatter3d(
                            x=df_setofproteins_PCA.PC1,
                            y=df_setofproteins_PCA.PC2,
                            z=df_setofproteins_PCA.PC3,
                            hovertext=df_setofproteins_PCA[DataFrameStrings.GENE_NAMES],
                            mode="markers",
                            name=maps,
                        )
                    )

            pca_figure.update_layout(
                autosize=False,
                width=500,
                height=500,
                title="PCA plot for <br>the protein cluster: {}".format(
                    cluster_of_interest
                ),
                template="simple_white",
            )
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
            df_cluster_unfiltered = self._get_marker_proteins_unfiltered(cluster)
            df_allclusters_01_unfiltered_mapfracunstacked = pd.concat(
                [df_allclusters_01_unfiltered_mapfracunstacked, df_cluster_unfiltered]
            )

            # filter for coverage and calculate distances
            df_cluster = df_cluster_unfiltered.dropna()
            if len(df_cluster) == 0:
                continue
            df_distances_aggregated, df_distances_individual = (
                self._calc_cluster_distances(df_cluster)
            )
            df_alldistances_individual_mapfracunstacked = pd.concat(
                [df_alldistances_individual_mapfracunstacked, df_distances_individual]
            )
            df_alldistances_aggregated_mapunstacked = pd.concat(
                [df_alldistances_aggregated_mapunstacked, df_distances_aggregated]
            )
        if len(df_alldistances_individual_mapfracunstacked) == 0:
            self.df_distance_noindex = pd.DataFrame(
                columns=[
                    DataFrameStrings.PROTEIN_IDS,
                    DataFrameStrings.GENE_NAMES,
                    DataFrameStrings.MAP,
                    DataFrameStrings.CLUSTER,
                    "distance",
                ]
            )
            self.df_allclusters_01_unfiltered_mapfracunstacked = pd.DataFrame(
                columns=[
                    DataFrameStrings.PROTEIN_IDS,
                    DataFrameStrings.GENE_NAMES,
                    DataFrameStrings.MAP,
                    DataFrameStrings.CLUSTER,
                    "distance",
                ]
            )
            self.df_allclusters_clusterdist_fracunstacked_unfiltered = pd.DataFrame(
                columns=[DataFrameStrings.FRACTION]
            )
            self.df_allclusters_clusterdist_fracunstacked = pd.DataFrame(
                columns=[DataFrameStrings.FRACTION]
            )
            self.genenames_sortedout_list = "No clusters found"
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        else:
            df_alldistances_aggregated_mapunstacked.columns.name = DataFrameStrings.MAP
            ## Get compatibility with plotting functions, by mimicking assignment of old functions:
            # old output of distance_calculation
            self.df_distance_noindex = (
                df_alldistances_aggregated_mapunstacked.stack(DataFrameStrings.MAP)
                .reset_index()
                .rename({0: "distance"}, axis=1)
            )
            self.analysis_summary_dict["Manhattan distances"] = (
                self.df_distance_noindex.to_json()
            )
            # old output of multiple_iterations
            # self.df_allclusters_clusterdist_fracunstacked_unfiltered --> this won't exist anymore, replaced by:
            self.df_allclusters_01_unfiltered_mapfracunstacked = (
                df_allclusters_01_unfiltered_mapfracunstacked
            )
            # kept for testing of quantification table:
            self.df_allclusters_clusterdist_fracunstacked_unfiltered = (
                df_allclusters_01_unfiltered_mapfracunstacked.stack(
                    DataFrameStrings.MAP
                )
            )
            # same as before, but now already abs
            self.df_allclusters_clusterdist_fracunstacked = (
                df_alldistances_individual_mapfracunstacked.stack(DataFrameStrings.MAP)
            )
            df_dist_to_median = self.df_allclusters_clusterdist_fracunstacked.stack(
                DataFrameStrings.FRACTION
            )
            df_dist_to_median.name = "distance"
            df_dist_to_median = df_dist_to_median.reindex(
                index=natsort.natsorted(df_dist_to_median.index)
            )
            self.analysis_summary_dict["Distances to the median profile"] = (
                df_dist_to_median.reset_index().to_json()
            )
            self.genenames_sortedout_list = [
                el
                for el in df_allclusters_01_unfiltered_mapfracunstacked.index.get_level_values(
                    DataFrameStrings.GENE_NAMES
                )
                if el
                not in df_alldistances_individual_mapfracunstacked.index.get_level_values(
                    DataFrameStrings.GENE_NAMES
                )
            ]

            return (
                df_alldistances_individual_mapfracunstacked,
                df_alldistances_aggregated_mapunstacked,
                df_allclusters_01_unfiltered_mapfracunstacked,
            )

    def _get_marker_proteins_unfiltered(self, cluster):
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
        df_in = self.df_01_stacked[DataFrameStrings.NORMALIZED_PROFILE].unstack(
            DataFrameStrings.FRACTION
        )
        markers = self.markerproteins[cluster]

        # retrieve marker proteins
        df_cluster_unfiltered = pd.DataFrame()
        for marker in markers:
            try:
                df_p = df_in.xs(
                    marker, level=DataFrameStrings.PROTEIN_IDS, axis=0, drop_level=False
                )
            except:
                continue
            df_cluster_unfiltered = pd.concat([df_cluster_unfiltered, df_p])
        if len(df_cluster_unfiltered) == 0:
            return df_cluster_unfiltered

        # Unstack maps and add Cluster to index
        df_cluster_unfiltered = df_cluster_unfiltered.unstack(DataFrameStrings.MAP)
        df_cluster_unfiltered.set_index(
            pd.Index(
                np.repeat(cluster, len(df_cluster_unfiltered)),
                name=DataFrameStrings.CLUSTER,
            ),
            append=True,
            inplace=True,
        )

        return df_cluster_unfiltered

    def _calc_cluster_distances(
        self, df_cluster, complex_profile=np.median, distance_measure="manhattan"
    ):
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

        ref_profile = pd.DataFrame(
            df_cluster.apply(complex_profile, axis=0, result_type="expand")
        ).T
        df_distances_individual = df_cluster.apply(
            lambda x: np.abs(x - ref_profile.iloc[0, :]), axis=1
        )

        # loop over maps
        maps = set(df_cluster.columns.get_level_values(DataFrameStrings.MAP))
        for m in maps:
            if distance_measure == "manhattan":
                d_m = pw.manhattan_distances(
                    df_cluster.xs(m, level=DataFrameStrings.MAP, axis=1),
                    ref_profile.xs(m, level=DataFrameStrings.MAP, axis=1),
                )
            else:
                raise ValueError(distance_measure)
            d_m = pd.DataFrame(d_m, columns=[m], index=df_cluster.index)
            df_distances_aggregated = pd.concat([df_distances_aggregated, d_m], axis=1)

        df_distances_aggregated.columns.set_names(
            names=DataFrameStrings.MAP, inplace=True
        )
        return df_distances_aggregated, df_distances_individual

    def run_outliertest(
        self,
        cond_1="",
        cond_2="",
        pairs=[],
        proportion=0,
        iterations=1,
        rmode="median",
        stop_at_95_05=True,
        prefilter=0.9,
        canvas=None,
    ):
        mrhash = f"{cond_1}{cond_2}{''.join([''.join(p) for p in pairs])}{proportion}{iterations}{stop_at_95_05}{rmode}{prefilter}"
        try:
            self.mr_results.keys()
        except:
            self.mr_results = dict()
        if mrhash in self.mr_results.keys():
            mr = self.mr_results[mrhash]
            if canvas is not None:
                canvas.objects = [
                    "This analysis has already been run and is shown below. No temporary data collected during processing can be shown."
                ]
            return mr

        cols = []
        for col in self.correlations.columns:
            colel = col.split("_")
            if colel[0] == cond_1:
                reps = [el[0] for el in pairs]
                if colel[1] in reps and colel[3] in reps:
                    cols.append(col)
            elif colel[0] == cond_2:
                reps = [el[1] for el in pairs]
                if colel[1] in reps and colel[3] in reps:
                    cols.append(col)

        min_correl = self.correlations[cols].min(axis=1)

        df01 = (
            self.df_01_stacked[DataFrameStrings.NORMALIZED_PROFILE]
            .unstack([DataFrameStrings.MAP, DataFrameStrings.FRACTION])[
                [
                    el
                    for el in self.map_names
                    if el.startswith(cond_1 + "_") or el.startswith(cond_2 + "_")
                ]
            ]
            .dropna()
            .stack([DataFrameStrings.FRACTION])
        )
        df01 = (
            pd.DataFrame(df01)
            .join(pd.Series(min_correl[min_correl > 0.8], name="filter"), how="left")
            .dropna()
            .query("filter >= @prefilter")
            .drop("filter", axis=1)
            .unstack([DataFrameStrings.FRACTION])
        )

        pvals = pd.DataFrame(
            index=pd.MultiIndex.from_tuples([], names=df01.index.names)
        )
        deltas = pd.DataFrame(
            index=pd.MultiIndex.from_tuples([], names=df01.index.names)
        )
        plot_it = px.scatter(
            template="simple_white", title="Approximation over iterations"
        ).update_layout(
            xaxis_title="Iteration", yaxis_title="95th percentile of delta(median)"
        )
        if canvas is not None:
            canvas.objects = ["Status", f"Number of proteins is {len(df01)}", plot_it]

        for p in pairs:
            m1 = cond_1 + "_" + p[0]
            m2 = cond_2 + "_" + p[1]
            delta = df01[[m1, m2]].stack(DataFrameStrings.FRACTION)
            delta = delta[m2] - delta[m1]

            if canvas is not None:
                canvas.append(
                    px.histogram(
                        delta.reset_index(),
                        x=0,
                        facet_col=DataFrameStrings.FRACTION,
                        title=f"Delta profile distribution {m2}-{m1}",
                        template="simple_white",
                        nbins=50,
                    )
                    .update_layout(
                        width=800, height=200, margin=dict(l=0, r=0, b=0, t=50)
                    )
                    .update_xaxes(title="delta")
                )

            delta = delta.unstack(DataFrameStrings.FRACTION)
            delta.columns = [f"delta {m2}-{m1} {el}" for el in delta.columns]
            if canvas is not None:
                canvas.objects[0].object = f"**Calculating {m2}-{m1}**"

            pv, med = outlier_test(
                delta.values,
                iterations=iterations,
                support_fraction=proportion,
                stop_at_95_05=stop_at_95_05,
            )

            pvals = pvals.join(
                pd.DataFrame(
                    np.array(pv).T, index=delta.index, columns=[f"p-value {m2}-{m1}"]
                ),
                how="outer",
            )
            deltas = deltas.join(delta, how="outer")

            plot_it.add_scatter(
                x=[i for i in range(2, len(med) + 1)],
                y=[
                    np.quantile(abs(med[i] - med[i - 1]) / med[i - 1], 0.95)
                    for i in range(1, len(med))
                ],
                name=f"{m2}-{m1}",
            )

        fisherp = pvals.apply(lambda x: combine_pvalues(x, method="fisher")[1], axis=1)
        mscore = -np.log10(multipletests(fisherp, method="fdr_bh")[1])
        pvals.insert(0, "M", mscore)

        if canvas is not None:
            canvas.objects[0].object = "**Calculating R-scores**"
        correlations = rscores_multi(
            df01, pairs=[[cond_1 + "_" + p[0], cond_2 + "_" + p[1]] for p in pairs]
        )
        correlations.columns = [f"correlation {col}" for col in correlations.columns]

        rscore = correlations.apply(
            lambda x: x.min()
            if rmode == "smallest correlation"
            else x.max()
            if rmode == "largest correlation"
            else x.median(),
            axis=1,
        )

        correlations.insert(0, "R", rscore)

        mr = pvals.join(correlations).join(deltas)

        if canvas is not None:
            canvas.objects[0].object = "**Done**"

        self.mr_results[mrhash] = mr.copy()

        return mr

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
            df_setofproteins = self.df_allclusters_01_unfiltered_mapfracunstacked.xs(
                cluster_of_interest, level=DataFrameStrings.CLUSTER, axis=0
            )
            df_setofproteins_median = (
                df_setofproteins.dropna()
                .xs(map_of_interest, level=DataFrameStrings.MAP, axis=1)
                .median(axis=0)
            )

            # fractions get sorted
            df_setofproteins = df_setofproteins.xs(
                map_of_interest, level=DataFrameStrings.MAP, axis=1
            ).stack(DataFrameStrings.FRACTION)
            df_setofproteins = df_setofproteins.reindex(
                index=natsort.natsorted(df_setofproteins.index)
            )
            df_setofproteins.name = DataFrameStrings.NORMALIZED_PROFILE

            # make it available for plotting
            df_setofproteins = df_setofproteins.reindex(
                index=natsort.natsorted(df_setofproteins.index)
            )
            df_setofproteins = df_setofproteins.reset_index()
            abundance_profiles_figure = px.line(
                df_setofproteins,
                x=DataFrameStrings.FRACTION,
                y=DataFrameStrings.NORMALIZED_PROFILE,
                color=DataFrameStrings.GENE_NAMES,
                line_group="Sequence"
                if "Sequence" in df_setofproteins.columns
                else DataFrameStrings.GENE_NAMES,
                template="simple_white",
                hover_data=[DataFrameStrings.PROTEIN_IDS, DataFrameStrings.GENE_NAMES],
                title="Relative abundance profile for {} of <br>the protein cluster: {}".format(
                    map_of_interest, cluster_of_interest
                ),
            )

            df_setofproteins_median.name = DataFrameStrings.NORMALIZED_PROFILE

            # fractions get sorted
            df_setofproteins_median = df_setofproteins_median.reindex(
                index=natsort.natsorted(df_setofproteins_median.index)
            )

            # make it available for plotting
            df_setofproteins_median = df_setofproteins_median.reset_index()
            df_setofproteins_median.insert(
                0,
                DataFrameStrings.GENE_NAMES,
                np.repeat("Median profile", len(df_setofproteins_median)),
            )

            abundance_profiles_and_median_figure = (
                abundance_profiles_figure.add_scatter(
                    x=df_setofproteins_median[DataFrameStrings.FRACTION],
                    y=df_setofproteins_median[DataFrameStrings.NORMALIZED_PROFILE],
                    name="Median profile",
                )
            )
            # dash lines for proteins that have insufficient coverage across maps
            abundance_profiles_and_median_figure.for_each_trace(
                lambda x: x.update(line={"dash": "dash"}),
                selector=lambda x: x.name in self.genenames_sortedout_list,
            )

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

        df_quantification_overview = (
            self.df_allclusters_clusterdist_fracunstacked_unfiltered.xs(
                cluster_of_interest, level=DataFrameStrings.CLUSTER, axis=0
            )[self.fractions[0]].unstack(DataFrameStrings.MAP)
        )
        if "Sequence" in df_quantification_overview.index.names:
            df_quantification_overview = df_quantification_overview.droplevel(
                [
                    i
                    for i in df_quantification_overview.index.names
                    if i not in ["Sequence", DataFrameStrings.GENE_NAMES]
                ]
            )
        else:
            df_quantification_overview = df_quantification_overview.droplevel(
                [
                    i
                    for i in df_quantification_overview.index.names
                    if not i == DataFrameStrings.GENE_NAMES
                ]
            )
        df_quantification_overview = df_quantification_overview.notnull().replace(
            {True: "x", False: "-"}
        )

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
        df_distance_map_cluster_gene_in_index = df_distance_noindex.set_index(
            [
                DataFrameStrings.GENE_NAMES,
                DataFrameStrings.MAP,
                DataFrameStrings.CLUSTER,
            ]
        )
        if "Sequence" in df_distance_map_cluster_gene_in_index.columns:
            df_distance_map_cluster_gene_in_index.set_index(
                "Sequence", append=True, inplace=True
            )

        df_cluster_xmaps_distance_with_index = pd.DataFrame()

        try:
            # for each individual map and a defined cluster data will be extracted from the dataframe
            # "df_distance_map_cluster_gene_in_index" and appended to the new dataframe df_cluster_xmaps_distance_with_index
            for maps in map_names:
                plot_try = df_distance_map_cluster_gene_in_index.xs(
                    (cluster_of_interest, maps),
                    level=[DataFrameStrings.CLUSTER, DataFrameStrings.MAP],
                    drop_level=False,
                )
                df_cluster_xmaps_distance_with_index = pd.concat(
                    [df_cluster_xmaps_distance_with_index, plot_try]
                )

            df_cluster_xmaps_distance_with_index["Combined Maps"] = "Combined Maps"

            # number of proteins within one cluster
            self.proteins_quantified_across_all_maps = (
                df_cluster_xmaps_distance_with_index.unstack(
                    DataFrameStrings.MAP
                ).shape[0]
            )

            # index will be reset, required by px.box
            df_cluster_xmaps_distance = (
                df_cluster_xmaps_distance_with_index.reset_index()
            )

            distance_boxplot_figure = go.Figure()
            distance_boxplot_figure.add_trace(
                go.Box(
                    x=df_cluster_xmaps_distance[DataFrameStrings.MAP],
                    y=df_cluster_xmaps_distance["distance"],
                    boxpoints="all",
                    whiskerwidth=0.2,
                    marker_size=2,
                    hovertext=df_cluster_xmaps_distance[DataFrameStrings.GENE_NAMES],
                )
            )

            distance_boxplot_figure.add_trace(
                go.Box(
                    x=df_cluster_xmaps_distance["Combined Maps"],
                    y=df_cluster_xmaps_distance["distance"],
                    boxpoints="all",
                    whiskerwidth=0.2,
                    marker_size=2,
                    hovertext=df_cluster_xmaps_distance[DataFrameStrings.GENE_NAMES],
                )
            )

            distance_boxplot_figure.update_layout(
                title="Manhattan distance distribution for <br>the protein cluster: {}".format(
                    cluster_of_interest
                ),
                autosize=False,
                showlegend=False,
                width=500,
                height=500,
                # black box around the graph
                xaxis=go.layout.XAxis(
                    linecolor="black",
                    linewidth=1,
                    title=DataFrameStrings.MAP,
                    mirror=True,
                ),
                yaxis=go.layout.YAxis(
                    linecolor="black", linewidth=1, title="distance", mirror=True
                ),
                template="simple_white",
            )

            return distance_boxplot_figure

        except:
            self.cache_cluster_quantified = False

    def calculate_pairwise_correlations(self, metric="cosine correlation"):
        """
        A distribution plot of the profile scatter in each experiment is generated, with variable distance metric and consolidation of replicates.

        Args:
            self:
                df_01_stacked: df, indexed
            metric: 'cosine correlation', 'pearson correlation'

        Returns:
            self attribute correlations with one column per within condition pair of correlation values
        """

        # Option dictionaries

        metrics = {
            "cosine correlation": "cosine",
            "pearson correlation": lambda x, y: np.corrcoef(x, y)[0][1],
        }

        # Option assertion
        assert metric in metrics.keys()

        df = (
            self.df_01_stacked.unstack(DataFrameStrings.MAP)
            .xs(DataFrameStrings.NORMALIZED_PROFILE, axis=1, level=DataFrameStrings.SET)
            .copy()
        )
        df.index = df.index.droplevel(
            [DataFrameStrings.GENE_NAMES, DataFrameStrings.COMPARTMENT]
        )

        # Calculate and consolidate distances
        distances = pd.DataFrame()
        for cond in self.conditions:
            cond += "_" if cond != "" else ""
            df_m = (
                df[[el for el in self.map_names if el.startswith(cond)]]
                .dropna()
                .stack([DataFrameStrings.MAP])
                .unstack(DataFrameStrings.FRACTION)
            )
            maps = list(set(df_m.index.get_level_values(DataFrameStrings.MAP)))

            distances_m = pd.DataFrame()

            # loop over pairs of maps
            for i, mapi in enumerate(maps):
                for j, mapj in enumerate(maps):
                    # only look at each comparison once
                    if j <= i:
                        continue
                    dist = pw.paired_distances(
                        df_m.xs(mapi, level=DataFrameStrings.MAP, axis=0).values,
                        df_m.xs(mapj, level=DataFrameStrings.MAP, axis=0).values,
                        metric=metrics[metric],
                    )
                    dist = pd.Series(dist, name="_".join([mapi, mapj]))
                    if metric == "cosine correlation":
                        dist = 1 - dist
                    distances_m = pd.concat([distances_m, dist], axis=1)
            distances_m.index = df_m.xs(
                maps[0], level=DataFrameStrings.MAP, axis=0
            ).index

            distances = pd.concat([distances, distances_m], axis=1)

        self.correlations = distances

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
                plot_try = self.df_allclusters_clusterdist_fracunstacked.xs(
                    (cluster_of_interest, maps),
                    level=[DataFrameStrings.CLUSTER, DataFrameStrings.MAP],
                    drop_level=False,
                )
                df_boxplot_manymaps = pd.concat([df_boxplot_manymaps, plot_try])

            self.df_boxplot_manymaps = df_boxplot_manymaps

            # index will be reset, required by px.violin
            df_boxplot_manymaps = abs(
                df_boxplot_manymaps.stack(DataFrameStrings.FRACTION)
            )

            df_boxplot_manymaps.name = "distance"

            df_boxplot_manymaps = df_boxplot_manymaps.reindex(
                index=natsort.natsorted(df_boxplot_manymaps.index)
            )

            df_boxplot_manymaps = df_boxplot_manymaps.reset_index()

            # box plot will be generated, every fraction will be displayed in a single plot
            distance_to_median_boxplot_figure = px.box(
                df_boxplot_manymaps,
                x=DataFrameStrings.MAP,
                y="distance",
                facet_col=DataFrameStrings.FRACTION,
                facet_col_wrap=2,
                boxmode="overlay",
                height=900,
                width=700,
                points="all",
                hover_name=DataFrameStrings.GENE_NAMES,
                template="simple_white",
                title="Distribution of the distance to the median for <br>the protein cluster: {}".format(
                    cluster_of_interest
                ),
            )

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
        df_distance_map_cluster_gene_in_index = df_distance_noindex.set_index(
            [
                DataFrameStrings.GENE_NAMES,
                DataFrameStrings.MAP,
                DataFrameStrings.CLUSTER,
            ]
        )
        map_names = self.map_names

        df_overview = pd.DataFrame()

        for clusters in self.markerproteins:
            # if a certain cluster is not available in the dataset at all
            try:
                for maps in map_names:
                    df_dist_map_cluster = df_distance_map_cluster_gene_in_index.xs(
                        (clusters, maps),
                        level=[DataFrameStrings.CLUSTER, DataFrameStrings.MAP],
                        drop_level=False,
                    )
                    statistic_table = {
                        "range": (df_dist_map_cluster["distance"].max(axis=0))
                        - (df_dist_map_cluster["distance"].min(axis=0)),
                        "median": df_dist_map_cluster["distance"].median(axis=0),
                        "standardeviation": df_dist_map_cluster["distance"].std(axis=0),
                        DataFrameStrings.CLUSTER: clusters,
                        DataFrameStrings.MAP: maps,
                    }
                    statistic_series = pd.Series(data=statistic_table)
                    df_statistic_table_individual_cluster = pd.DataFrame(
                        statistic_series
                    ).T
                    df_overview = pd.concat(
                        [df_overview, df_statistic_table_individual_cluster]
                    )

                df_dist_cluster = df_distance_map_cluster_gene_in_index.xs(
                    clusters, level=DataFrameStrings.CLUSTER
                )
                statistic_table_combined = {
                    "range": (df_dist_cluster["distance"].max(axis=0))
                    - (df_dist_cluster["distance"].min(axis=0)),
                    "median": df_dist_cluster["distance"].median(axis=0),
                    "standardeviation": df_dist_cluster["distance"].std(axis=0),
                    DataFrameStrings.CLUSTER: clusters,
                    DataFrameStrings.MAP: "combined maps",
                }
                statistic_series_combined = pd.Series(data=statistic_table_combined)
                df_statistic_table_individual_cluster = pd.DataFrame(
                    statistic_series_combined
                ).T
                df_overview = pd.concat(
                    [df_overview, df_statistic_table_individual_cluster]
                )

            except:
                continue
        try:
            df_overview.set_index(
                [DataFrameStrings.CLUSTER, DataFrameStrings.MAP], inplace=True
            )
            df_overview.sort_index(axis=0, level=0, inplace=True)
        except:
            df_overview = pd.DataFrame()

        self.analysis_summary_dict["Overview table"] = (
            df_overview.reset_index().to_json()
        )
        self.analysed_datasets_dict[self.expname] = self.analysis_summary_dict.copy()
        # self.analysis_summary_dict.clear()

        return df_overview

    def reframe_df_01ORlog_for_Perseus(self, df_01ORlog):
        """ "
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

        # df_01_filtered_combined = df_01_filtered_combined.stack(["Experiment", "Map"]).swaplevel(0,1, axis=0).dropna(axis=1)
        index_ExpMap = (
            df_01ORlog_svm.index.get_level_values(DataFrameStrings.MAP)
            + "_"
            + df_01ORlog_svm.index.get_level_values(DataFrameStrings.FRACTION)
        )
        index_ExpMap.name = "Map_Frac"
        df_01ORlog_svm.set_index(index_ExpMap, append=True, inplace=True)

        df_01ORlog_svm.index = df_01ORlog_svm.index.droplevel(
            [DataFrameStrings.MAP, DataFrameStrings.FRACTION]
        )
        df_01ORlog_svm = df_01ORlog_svm.unstack("Map_Frac")
        # df_01ORlog_svm = df_01ORlog_svm.dropna(axis=0, subset=df_01ORlog_svm.loc[[], [DataFrameStrings.NORMALIZED_PROFILE]].columns)
        df_01ORlog_svm.columns = [
            "_".join(col) for col in df_01ORlog_svm.columns.values
        ]
        df_01ORlog_svm.rename(
            index={"undefined": np.nan},
            level=DataFrameStrings.COMPARTMENT,
            inplace=True,
        )

        return df_01ORlog_svm


class SpatialDataSetComparison:
    analysed_datasets_dict = {}
    css_color = SpatialDataSet.css_color
    # cache_stored_SVM = True

    def __init__(
        self, ref_exp=None, **kwargs
    ):  # clusters_for_ranking=["Proteasome", "Lysosome"]
        # self.clusters_for_ranking = clusters_for_ranking
        self.ref_exp = ref_exp
        self.json_dict = {}
        # self.fractions, self.map_names = [], []  #self.df_01_stacked, self.df_log_stacked = pd.DataFrame(), pd.DataFrame()
        # collapse_maps,collapse_cluster,  cluster_of_interest_comparison, multi_choice, multi_choice_venn, x_PCA_comp, y_PCA_comp

        # if "organism" not in kwargs.keys():
        #    self.markerproteins = self.markerproteins_set["Human - Swissprot"]
        # else:
        #    assert kwargs["organism"] in self.markerproteins_set.keys()
        #    self.markerproteins = self.markerproteins_set[kwargs["organism"]]
        #    del kwargs["organism"]

        self.exp_names, self.exp_map_names = [], []

        self.df_01_filtered_combined, self.df_distance_comp = (
            pd.DataFrame(),
            pd.DataFrame(),
        )
        self.df_quantity_pr_pg_combined = pd.DataFrame()
        self.df_svm_performance = pd.DataFrame()
        self.svm_results = dict()
        self.svm_runs = dict()
        self.parameters = dict()

    def read_jsonFile(self):  # , content=None
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
        if self.ref_exp == None:
            self.ref_exp = self.exp_names[0]

        df_01_combined = pd.DataFrame()
        for exp_name in json_dict.keys():
            for data_type in json_dict[exp_name].keys():
                if data_type == "0/1 normalized data":
                    df_01_toadd = pd.read_json(json_dict[exp_name][data_type])
                    df_01_toadd.insert(
                        0,
                        DataFrameStrings.EXPERIMENT,
                        np.repeat(exp_name, len(df_01_toadd)),
                    )
                    if len(df_01_combined) == 0:
                        df_01_combined = df_01_toadd.copy()
                    else:
                        df_01_combined = pd.concat(
                            [df_01_combined, df_01_toadd], sort=False, axis=0
                        )

                elif (
                    data_type == "quantity: profiles/protein groups"
                    and exp_name == list(json_dict.keys())[0]
                ):
                    df_quantity_pr_pg_combined = pd.read_json(
                        json_dict[exp_name][data_type]
                    )
                    df_quantity_pr_pg_combined[DataFrameStrings.EXPERIMENT] = exp_name

                elif (
                    data_type == "quantity: profiles/protein groups"
                    and exp_name != list(json_dict.keys())[0]
                ):
                    df_quantity_pr_pg_toadd = pd.read_json(
                        json_dict[exp_name][data_type]
                    )
                    df_quantity_pr_pg_toadd[DataFrameStrings.EXPERIMENT] = exp_name
                    df_quantity_pr_pg_combined = pd.concat(
                        [df_quantity_pr_pg_combined, df_quantity_pr_pg_toadd]
                    )

                elif data_type == "Analysis parameters":
                    self.analysis_parameters_total[exp_name] = json_dict[exp_name][
                        data_type
                    ]

                elif data_type == "SVM results":
                    for res in self.json_dict[exp_name]["SVM results"].keys():
                        misclassification = pd.read_json(
                            self.json_dict[exp_name]["SVM results"][res][
                                "misclassification"
                            ]
                        )
                        if "T: True group" in misclassification.columns:
                            misclassification.set_index("T: True group", inplace=True)
                            misclassification.columns = [
                                el.split(": ")[1] for el in misclassification.columns
                            ]
                            misclassification.index = [
                                el.split(": ")[1] for el in misclassification.index
                            ]
                        self.add_svm_result(
                            exp_name,
                            misclassification,
                            name=res,
                            source="direct",
                            prediction=pd.read_json(
                                self.json_dict[exp_name]["SVM results"][res][
                                    "prediction"
                                ]
                            ),
                            comment=self.json_dict[exp_name]["SVM results"][res][
                                "comment"
                            ],
                            overwrite=True,  # supposed to always exceed upload from old format
                        )

                elif data_type == "Misclassification Matrix":
                    try:
                        self.add_svm_result(
                            exp_name,
                            pd.read_json(
                                self.json_dict[exp_name]["Misclassification Matrix"]
                            ),
                            comment="read from old json file version",
                            overwrite=False,
                        )
                    except:
                        # should only happen if this has been loaded before and both SVM results and Misclassification Matrix are present
                        pass

                # try:
                #    for paramters in json_dict[exp_name][data_type].keys():
                #        if paramters=="acquisition":
                #            acquisition_loaded.append(json_dict[exp_name][data_type][paramters])
                #        #elif parameters=="Non valid profiles":
                # except:
                #    continue
                #
        ### New code for alignment:
        df_01_combined = df_01_combined.set_index(
            [
                DataFrameStrings.EXPERIMENT,
                DataFrameStrings.GENE_NAMES,
                DataFrameStrings.COMPARTMENT,
            ]
        )
        if (
            DataFrameStrings.ORIGINAL_PROTEIN_IDS in df_01_combined.columns
            and DataFrameStrings.PROTEIN_IDS in df_01_combined.columns
        ):
            df_01_combined.set_index(
                [DataFrameStrings.ORIGINAL_PROTEIN_IDS], append=True, inplace=True
            )
            df_01_combined.drop(DataFrameStrings.PROTEIN_IDS, axis=1, inplace=True)
        elif DataFrameStrings.ORIGINAL_PROTEIN_IDS in df_01_combined.columns:
            df_01_combined.set_index(
                DataFrameStrings.ORIGINAL_PROTEIN_IDS, append=True, inplace=True
            )
        elif DataFrameStrings.PROTEIN_IDS in df_01_combined.columns:
            df_01_combined.set_index(
                pd.Index(
                    df_01_combined[DataFrameStrings.PROTEIN_IDS],
                    name=DataFrameStrings.ORIGINAL_PROTEIN_IDS,
                ),
                append=True,
                inplace=True,
            )
            df_01_combined.drop(DataFrameStrings.PROTEIN_IDS, axis=1, inplace=True)
        else:
            raise KeyError("No Protein IDs were found while loading the json file")

        df_01_combined.columns = pd.MultiIndex.from_tuples(
            [el.split("?")[1::] for el in df_01_combined.columns],
            names=[DataFrameStrings.MAP, DataFrameStrings.FRACTION],
        )
        # df_01_combined = df_01_combined.stack("Map").dropna().unstack("Map")
        if len(df_01_combined) == 0:
            raise ValueError(
                "Merged dataframe is empty. This is most likely due to unequal naming of fractions."
            )
        if len(
            set(df_01_combined.index.get_level_values(DataFrameStrings.EXPERIMENT))
        ) != len(self.exp_names):
            raise ValueError(
                "At least one of the datsets was removed completely, likely because it is lacking one of the fractions in the other experiments."
            )
        df_01_filtered_combined, id_alignment = align_datasets(df_01_combined)
        self.id_alignment = id_alignment
        index_ExpMap = (
            df_01_filtered_combined.index.get_level_values(DataFrameStrings.EXPERIMENT)
            + "_"
            + df_01_filtered_combined.index.get_level_values(DataFrameStrings.MAP)
        )
        index_ExpMap.name = DataFrameStrings.EXP_MAP
        df_01_filtered_combined.set_index(index_ExpMap, append=True, inplace=True)

        self.exp_map_names = list(index_ExpMap.unique())

        self.df_01_filtered_combined = df_01_filtered_combined

        self.df_quantity_pr_pg_combined = df_quantity_pr_pg_combined

        # find out whether mixed fractions were loaded
        columns = (
            df_01_filtered_combined.unstack([DataFrameStrings.EXPERIMENT])
            .dropna(axis=1, how="all")
            .columns
        )
        fractions = []
        for f, _ in columns:
            fractions.append(f)
        fractions = list(set(fractions))
        mixed = False
        for f in fractions:
            for e in self.exp_names:
                if (f, e) not in columns:
                    mixed = True
                    continue
        self.mixed = mixed
        self.fractions = fractions

        # This is to support files written between 1.0.0 and 1.0.3 (legacy) and files written from 1.0.4 onwards (version)
        if (
            "legacy" in json_dict[self.ref_exp]["Analysis parameters"]
            or "domaps version" in json_dict[self.ref_exp]["Analysis parameters"]
        ):
            self.markerproteins = json_dict[self.ref_exp][SettingStrings.COMPLEXES]
        # This is to support files written before 1.0
        else:
            try:
                organism = json_dict[list(json_dict.keys())[0]]["Analysis parameters"][
                    SettingStrings.ORGANISM
                ]
            except:
                organism = "Homo sapiens - Uniprot"

            marker_table = pd.read_csv(
                pkg_resources.resource_stream(
                    __name__, "annotations/complexes/{}.csv".format(organism)
                )
            )
            self.markerproteins = {
                k: v.replace(" ", "").split(",")
                for k, v in zip(
                    marker_table[DataFrameStrings.CLUSTER],
                    marker_table["Members - Protein IDs"],
                )
            }

        self.clusters_for_ranking = self.markerproteins.keys()

        self.color_maps = dict()
        compartments = sorted(
            set(
                df_01_filtered_combined.index.get_level_values(
                    DataFrameStrings.COMPARTMENT
                )
            )
        )
        compartments.pop(compartments.index("undefined"))
        self.color_maps["Compartments"] = dict(zip(compartments, self.css_color))
        self.color_maps["Compartments"]["undefined"] = "lightgrey"

        self.color_maps["Clusters"] = dict(
            zip(self.clusters_for_ranking, self.css_color)
        )
        self.color_maps["Clusters"]["Undefined"] = "lightgrey"

    def add_svm_result(
        self,
        experiment,
        misclassification,
        source="Perseus",
        name="default",
        prediction=pd.DataFrame(),
        comment="",
        overwrite=True,
    ):
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
            raise KeyError(
                f"Experiment {experiment} not found during SVM matrix storage"
            )

        if experiment not in self.svm_results.keys():
            self.svm_results[experiment] = dict()

        if name in self.svm_results[experiment].keys() and not overwrite:
            raise KeyError(
                f"SVM result named {name} already exists for experiment {experiment}. Rename it or set overwrite=True."
            )

        if source == "Perseus":
            misclassification.set_index("T: True group", inplace=True)
            misclassification.columns = [
                el.split(": ")[1] for el in misclassification.columns
            ]
            misclassification.index = [
                el.split(": ")[1] for el in misclassification.index
            ]
        elif source == "MetaMass":
            misclassification.index = misclassification.columns
            misclassification = misclassification.T
        elif source == "direct":
            misclassification.index = misclassification.columns
        else:
            raise ValueError(
                f"Source {source} for SVM results is not defined. Choose 'direct' and have true classes in rows and predicted classes in columns."
            )

        self.svm_results[experiment][name] = {
            "comment": comment,
            "misclassification": misclassification,
            "prediction": prediction,
        }

        try:
            self.df_svm_performance = self.df_svm_performance.query(
                "Experiment != @experiment or Set != @name"
            )
        except:
            pass
        self.df_svm_performance = pd.concat(
            [
                self.df_svm_performance,
                process_mc_errorbars(misclassification).set_index(
                    pd.MultiIndex.from_arrays(
                        [
                            np.repeat(experiment, 3 * (len(misclassification) + 3)),
                            np.repeat(name, 3 * (len(misclassification) + 3)),
                        ],
                        names=[DataFrameStrings.EXPERIMENT, DataFrameStrings.SET],
                    ),
                    append=True,
                ),
            ]
        )

    def perform_pca_comparison(self, n=3):
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
                df_pca: PCA processed dataframe
                    index: "Gene names", "Protein IDs", "Compartment", "Experiment",
                    columns: "PC1", "PC2", "PC3"
                    contains all protein IDs, that are consistent throughout all experiments
        """

        markerproteins = self.markerproteins.copy()

        df_mean = pd.DataFrame()
        for exp in self.exp_names:
            df_exp = (
                self.df_01_filtered_combined.stack(DataFrameStrings.FRACTION)
                .unstack(
                    [
                        DataFrameStrings.EXPERIMENT,
                        DataFrameStrings.MAP,
                        DataFrameStrings.EXP_MAP,
                    ]
                )[exp]
                .mean(axis=1)
                .to_frame(name=exp)
            )
            df_mean = pd.concat([df_mean, df_exp], axis=1)
        df_mean = (
            df_mean.rename_axis(DataFrameStrings.EXPERIMENT, axis="columns")
            .stack(DataFrameStrings.EXPERIMENT)
            .unstack(DataFrameStrings.FRACTION)
        )
        df_zscore = df_mean.apply(zscore, axis=0).replace(np.nan, 0)

        pca = PCA(n_components=n)

        df_pca = pd.DataFrame(
            pca.fit_transform(df_zscore),
            columns=[f"PC{el + 1}" for el in range(n)],
            index=df_zscore.index,
        )
        self.df_pca_loadings = pd.DataFrame(
            pca.components_, columns=df_zscore.columns, index=df_pca.columns
        ).T.reset_index()
        self.df_pca_var = (
            pd.DataFrame(pca.explained_variance_ratio_, index=df_pca.columns)
            .reset_index()
            .rename({"index": "Component", 0: "variance explained"}, axis=1)
        )

        ###only one df, make annotation at that time
        df_cluster = pd.DataFrame(
            [(k, i) for k, l in markerproteins.items() for i in l],
            columns=[DataFrameStrings.CLUSTER, DataFrameStrings.PROTEIN_IDS],
        )
        df_pca = df_pca.reset_index().merge(
            df_cluster, how="left", on=DataFrameStrings.PROTEIN_IDS
        )
        df_pca.Cluster.replace(np.NaN, "Undefined", inplace=True)

        df_cluster_pca = (
            df_zscore.reset_index()
            .merge(df_cluster, how="right", on=DataFrameStrings.PROTEIN_IDS)
            .dropna()
        )
        df_cluster_pca.set_index(
            df_zscore.index.names + [DataFrameStrings.CLUSTER], inplace=True
        )
        df_cluster_pca = pd.DataFrame(
            pca.transform(df_cluster_pca),
            columns=[f"PC{el + 1}" for el in range(n)],
            index=df_cluster_pca.index,
        )

        self.df_pca = df_pca
        self.df_cluster_pca = df_cluster_pca

    def plot_pca_comparison(
        self, cluster_of_interest_comparison="Proteasome", multi_choice=["Exp1", "Exp2"]
    ):
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
            df_setofproteins_PCA = df_pca.loc[
                df_pca.Cluster == cluster_of_interest_comparison, :
            ]

            df_setofproteins_PCA = df_setofproteins_PCA.assign(
                Experiment_lexicographic_sort=pd.Categorical(
                    df_setofproteins_PCA[DataFrameStrings.EXPERIMENT],
                    categories=multi_choice,
                    ordered=True,
                )
            )
            df_setofproteins_PCA.sort_values(
                "Experiment_lexicographic_sort", inplace=True
            )

            pca_figure = px.scatter_3d(
                df_setofproteins_PCA,
                x="PC1",
                y="PC2",
                z="PC3",
                color=DataFrameStrings.EXPERIMENT,
                template="simple_white",
                hover_data=[DataFrameStrings.GENE_NAMES],
            )

            pca_figure.update_layout(
                autosize=False,
                width=700,
                height=500,
                title="PCA plot for <br>the protein cluster: {}".format(
                    cluster_of_interest_comparison
                ),
                template="simple_white",
            )

            return pca_figure
        except:
            return "This protein cluster was not identified in all experiments"

    def plot_global_pca_comparison(
        self,
        cluster_of_interest_comparison="Proteasome",
        x_PCA="PC1",
        y_PCA="PC3",
        markerset_or_cluster=False,
        multi_choice=["Exp1", "Exp2"],
    ):
        """ "
        PCA plot will be generated

        Args:
            self:
                df_organellarMarkerSet: df, columns: "Gene names", "Compartment", no index
                multi_choice: list of experiment names
                css_color: list of colors
                df_pca: PCA processed dataframe
                    index: "Gene names", "Protein IDs", "Compartment", "Experiment",
                    columns: "PC1", "PC2", "PC3"
                    contains all protein IDs, that are consistent throughout all experiments

        Returns:
            pca_figure: global PCA plot, clusters based on the markerset based (df_organellarMarkerSet) are color coded.
        """

        df_pca_exp = self.df_pca.loc[
            self.df_pca[DataFrameStrings.EXPERIMENT].isin(multi_choice)
        ]
        df_pca_exp.reset_index(inplace=True)

        compartment_color = self.color_maps["Compartments"]
        compartments = list(compartment_color.keys())

        cluster = self.markerproteins.keys()
        cluster_color = dict(zip(cluster, self.css_color))
        cluster_color["Undefined"] = "lightgrey"

        if markerset_or_cluster == True:
            df_pca = df_pca_exp[df_pca_exp.Cluster != "Undefined"].sort_values(
                by=DataFrameStrings.CLUSTER
            )
            df_pca = pd.concat([df_pca_exp[df_pca_exp.Cluster == "Undefined"], df_pca])
        else:
            for i in self.markerproteins[cluster_of_interest_comparison]:
                df_pca_exp.loc[
                    df_pca_exp[DataFrameStrings.PROTEIN_IDS] == i,
                    DataFrameStrings.COMPARTMENT,
                ] = "Selection"
            df_pca = df_pca_exp.assign(
                Compartment_lexicographic_sort=pd.Categorical(
                    df_pca_exp[DataFrameStrings.COMPARTMENT],
                    categories=[x for x in compartments],
                    ordered=True,
                )
            )
            df_pca.sort_values(
                ["Compartment_lexicographic_sort", DataFrameStrings.EXPERIMENT],
                inplace=True,
            )

        fig_global_pca = px.scatter(
            data_frame=df_pca,
            x=x_PCA,
            y=y_PCA,
            color=DataFrameStrings.COMPARTMENT
            if markerset_or_cluster == False
            else DataFrameStrings.CLUSTER,
            color_discrete_map=compartment_color
            if markerset_or_cluster == False
            else cluster_color,
            title="Protein subcellular localization by PCA",
            hover_data=[
                DataFrameStrings.PROTEIN_IDS,
                DataFrameStrings.GENE_NAMES,
                DataFrameStrings.COMPARTMENT,
            ],
            facet_col=DataFrameStrings.EXPERIMENT,
            facet_col_wrap=2,
            opacity=0.9,
            template="simple_white",
        )

        fig_global_pca.update_layout(
            autosize=False,
            width=1800 if markerset_or_cluster == False else 1600,
            height=400 * (int(len(multi_choice) / 2) + (len(multi_choice) % 2 > 0)),
            template="simple_white",
        )

        return fig_global_pca

    def _get_marker_proteins(self, experiments, cluster):
        df_in = self.df_01_filtered_combined.copy()
        markers = self.markerproteins[cluster]

        # retrieve marker proteins
        df_cluster = pd.DataFrame()
        for marker in markers:
            try:
                df_p = df_in.xs(
                    marker, level=DataFrameStrings.PROTEIN_IDS, axis=0, drop_level=False
                )
            except:
                continue
            df_cluster = pd.concat([df_cluster, df_p], axis=0)
        if len(df_cluster) == 0:
            return df_cluster

        # filter for all selected experiments
        df_cluster = df_cluster.droplevel(DataFrameStrings.EXP_MAP, axis=0)
        df_cluster = df_cluster.unstack(
            [DataFrameStrings.EXPERIMENT, DataFrameStrings.MAP]
        )
        if any(
            [
                el
                not in df_cluster.columns.get_level_values(DataFrameStrings.EXPERIMENT)
                for el in experiments
            ]
        ):
            return pd.DataFrame()
        drop_experiments = [
            el
            for el in df_cluster.columns.get_level_values(DataFrameStrings.EXPERIMENT)
            if el not in experiments
        ]
        if len(drop_experiments) > 0:
            df_cluster.drop(
                [
                    el
                    for el in df_cluster.columns.get_level_values(
                        DataFrameStrings.EXPERIMENT
                    )
                    if el not in experiments
                ],
                level=DataFrameStrings.EXPERIMENT,
                axis=1,
                inplace=True,
            )
        df_cluster = df_cluster.dropna(axis=1, how="all").dropna()
        if len(df_cluster) == 0:
            return df_cluster
        df_cluster.set_index(
            pd.Index(
                np.repeat(cluster, len(df_cluster)), name=DataFrameStrings.CLUSTER
            ),
            append=True,
            inplace=True,
        )

        return df_cluster

    def _calc_cluster_distances(
        self, df_cluster, complex_profile=np.median, distance_measure="manhattan"
    ):
        df_distances = pd.DataFrame()

        # loop over experiments
        experiments = set(
            df_cluster.columns.get_level_values(DataFrameStrings.EXPERIMENT)
        )
        for exp in experiments:
            df_exp = df_cluster.xs(exp, level=DataFrameStrings.EXPERIMENT, axis=1)
            ref_profile = pd.DataFrame(
                df_exp.apply(complex_profile, axis=0, result_type="expand")
            ).T

            # loop over maps
            maps = set(df_exp.columns.get_level_values(DataFrameStrings.MAP))
            for m in maps:
                if distance_measure == "manhattan":
                    d_m = pw.manhattan_distances(
                        df_exp.xs(m, level=DataFrameStrings.MAP, axis=1),
                        ref_profile.xs(m, level=DataFrameStrings.MAP, axis=1),
                    )
                else:
                    raise ValueError(distance_measure)
                d_m = pd.DataFrame(d_m, columns=[(exp, m)], index=df_exp.index)
                df_distances = pd.concat([df_distances, d_m], axis=1)

        df_distances.columns = pd.MultiIndex.from_tuples(
            df_distances.columns,
            names=[DataFrameStrings.EXPERIMENT, DataFrameStrings.MAP],
        )
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
            df_cluster = self._get_marker_proteins(experiments, cluster)
            if len(df_cluster) == 0:
                continue
            dists_cluster = self._calc_cluster_distances(df_cluster)
            df_distances = pd.concat([df_distances, dists_cluster])
        if len(df_distances) == 0:
            raise ValueError(
                "Could not calculate biological precision, because no complexes could be extracted"
            )
        df_distances = (
            df_distances.stack([DataFrameStrings.EXPERIMENT, DataFrameStrings.MAP])
            .reset_index()
            .sort_values([DataFrameStrings.EXPERIMENT, DataFrameStrings.GENE_NAMES])
            .rename({0: "distance"}, axis=1)
        )
        df_distances.insert(
            0,
            DataFrameStrings.EXP_MAP,
            [
                "_".join([e, m])
                for e, m in zip(
                    df_distances[DataFrameStrings.EXPERIMENT],
                    df_distances[DataFrameStrings.MAP],
                )
            ],
        )

        self.df_distance_comp = df_distances

        return df_distances

    def get_complex_coverage(self, min_n=5):
        full_coverage = {}
        for complx in self.markerproteins.keys():
            df = self._get_marker_proteins(self.exp_names, complx)
            if len(df) >= min_n:
                full_coverage[complx] = len(df)
        partial_coverage = {}
        for exp in self.exp_names:
            for complx in self.markerproteins.keys():
                if complx in full_coverage.keys():
                    continue
                df = self._get_marker_proteins([exp], complx)
                # print(df)
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

    def plot_intramap_scatter_cluster(
        self,
        cluster_of_interest_comparison="Proteasome",
        collapse_maps=False,
        multi_choice=["Exp1", "Exp2"],
    ):
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

        # an error massage, if no Experiments are selected, will be displayed already, that is why: return ""
        if len(multi_choice) >= 1:
            pass
        else:
            return ""

        df_distance_comp = self.df_distance_comp.copy()
        # set categroical column, allowing lexicographic sorting
        df_distance_comp["Experiment_lexicographic_sort"] = pd.Categorical(
            df_distance_comp[DataFrameStrings.EXPERIMENT],
            categories=multi_choice,
            ordered=True,
        )
        df_distance_comp.sort_values(
            ["Experiment_lexicographic_sort", DataFrameStrings.MAP], inplace=True
        )

        if collapse_maps == False:
            # get only values form experiment of interest
            df_distance_selectedExp = df_distance_comp.loc[
                df_distance_comp[DataFrameStrings.EXPERIMENT].isin(multi_choice)
            ]
            # get only values form cluster of interest
            df_distance_selectedExp = df_distance_selectedExp.loc[
                df_distance_selectedExp[DataFrameStrings.CLUSTER]
                == cluster_of_interest_comparison
            ]

            if df_distance_selectedExp.shape[0] == 0:
                self.cache_cluster_quantified = False

            else:
                individual_distance_boxplot_figure = go.Figure()
                for i, exp in enumerate(multi_choice):
                    df_plot = df_distance_selectedExp[
                        df_distance_selectedExp[DataFrameStrings.EXPERIMENT] == exp
                    ]
                    individual_distance_boxplot_figure.add_trace(
                        go.Box(
                            x=[
                                df_plot[DataFrameStrings.EXPERIMENT],
                                df_plot[DataFrameStrings.MAP],
                            ],
                            y=df_plot["distance"],
                            # line=dict(color=pio.templates["simple_white"].layout["colorway"][i]),
                            boxpoints="all",
                            whiskerwidth=0.2,
                            marker_size=2,
                            name=exp,
                            hovertext=df_plot[DataFrameStrings.GENE_NAMES],
                        )
                    )

                individual_distance_boxplot_figure.update_layout(
                    boxmode="group",
                    xaxis_tickangle=90,
                    title="Manhattan distance distribution for <br>the protein cluster: {}".format(
                        cluster_of_interest_comparison
                    ),
                    autosize=False,
                    width=350 * len(multi_choice),
                    height=500,
                    xaxis=go.layout.XAxis(
                        linecolor="black",
                        linewidth=1,
                        title=DataFrameStrings.EXPERIMENT,
                        mirror=True,
                    ),
                    yaxis=go.layout.YAxis(
                        linecolor="black", linewidth=1, title="Distance", mirror=True
                    ),
                    template="simple_white",
                )
                individual_distance_boxplot_figure.update_layout(
                    **calc_width_categorical(individual_distance_boxplot_figure)
                )

                return individual_distance_boxplot_figure

        else:
            map_or_exp_names = multi_choice
            level_of_interest = DataFrameStrings.EXPERIMENT
            boxplot_color = DataFrameStrings.EXPERIMENT
            df_distance_selectedExp_global = df_distance_comp

            # "Gene names", "Map", "Cluster" and transferred into the index
            df_distance_selectedExp_global.set_index(
                [
                    DataFrameStrings.GENE_NAMES,
                    level_of_interest,
                    DataFrameStrings.CLUSTER,
                ],
                inplace=True,
            )

            df_cluster_xmaps_distance_global = pd.DataFrame()

            # for each individual map and a defined cluster data will be extracted from the dataframe
            # "df_distance_selectedExp_global" and appended to the new dataframe df_cluster_xmaps_distance_global
            for map_or_exp in map_or_exp_names:
                plot_try = df_distance_selectedExp_global.xs(
                    (cluster_of_interest_comparison, map_or_exp),
                    level=[DataFrameStrings.CLUSTER, level_of_interest],
                    drop_level=False,
                )
                df_cluster_xmaps_distance_global = pd.concat(
                    [df_cluster_xmaps_distance_global, plot_try]
                )

            df_cluster_xmaps_distance_global.sort_values(
                "Experiment_lexicographic_sort", inplace=True
            )
            df_cluster_xmaps_distance_global.reset_index(inplace=True)

            distance_boxplot_figure = px.box(
                df_cluster_xmaps_distance_global,
                x=level_of_interest,
                y="distance",
                points="all",
                hover_name=DataFrameStrings.GENE_NAMES,
                color=boxplot_color,
                template="simple_white",
                title="Global Manhattan distance distribution for the protein cluster: {}".format(
                    cluster_of_interest_comparison
                ),
            )

            distance_boxplot_figure.update_layout(
                autosize=False,
                width=250 * len(multi_choice),
                height=500,
                xaxis=go.layout.XAxis(
                    linecolor="black",
                    linewidth=1,
                    title=DataFrameStrings.MAP,
                    mirror=True,
                ),
                yaxis=go.layout.YAxis(
                    linecolor="black", linewidth=1, title="distance", mirror=True
                ),
                template="simple_white",
            )
            distance_boxplot_figure.update_layout(
                **calc_width_categorical(distance_boxplot_figure)
            )

            return distance_boxplot_figure

    def plot_intramap_scatter(
        self,
        normalization=np.median,
        aggregate_proteins=True,
        aggregate_maps=False,
        plot_type="strip",
        min_size=5,
        multi_choice=None,
        clusters_for_ranking=None,
        highlight=None,
    ):
        if multi_choice is None:
            multi_choice = self.exp_names
        if clusters_for_ranking is None:
            clusters_for_ranking = self.clusters_for_ranking
        if len(multi_choice) == 0 or len(clusters_for_ranking) == 0:
            return "Please provide at least one experiment and one cluster for ranking"

        df = self.df_distance_comp.copy()
        df.drop(
            [DataFrameStrings.EXP_MAP, "merge type", DataFrameStrings.COMPARTMENT],
            axis=1,
            inplace=True,
        )
        df = df[df[DataFrameStrings.EXPERIMENT].isin(multi_choice)]
        df = df[df[DataFrameStrings.CLUSTER].isin(clusters_for_ranking)]
        df = df.groupby(
            [
                DataFrameStrings.CLUSTER,
                DataFrameStrings.EXPERIMENT,
                DataFrameStrings.MAP,
            ]
        ).filter(lambda x: len(x) >= min_size)
        df.set_index([col for col in df.columns if col != "distance"], inplace=True)

        if aggregate_maps is not False:
            if aggregate_maps is True:
                aggregate_maps = np.nanmean
            df = pd.DataFrame(
                df.unstack(DataFrameStrings.MAP).apply(aggregate_maps, axis=1),
                columns=["distance"],
            )

        if normalization == np.median:
            df = df.groupby([DataFrameStrings.CLUSTER]).apply(
                lambda x: x / normalization(x)
            )
            df.reset_index(
                level=df.index.names.index(DataFrameStrings.CLUSTER),
                inplace=True,
                drop=True,
            )
        elif normalization in multi_choice:
            df = (
                df.unstack(DataFrameStrings.EXPERIMENT)
                .apply(lambda x: x / x[("distance", normalization)], axis=1)
                .stack(DataFrameStrings.EXPERIMENT)
            )
        else:
            pass

        if aggregate_proteins is not False:
            if aggregate_proteins is True:
                aggregate_proteins = np.nanmedian
            df = pd.DataFrame(
                df.unstack(
                    [DataFrameStrings.PROTEIN_IDS, DataFrameStrings.GENE_NAMES]
                ).apply(aggregate_proteins, axis=1),
                columns=["distance"],
            )

        plotargs = dict(
            x=DataFrameStrings.EXPERIMENT,
            y="distance",
            hover_data=df.index.names,
            template="simple_white",
        )

        df = df.sort_values("distance", ascending=False).sort_index(
            level=DataFrameStrings.EXPERIMENT,
            key=lambda x: pd.Index(
                [multi_choice.index(el) for el in x], name=DataFrameStrings.EXPERIMENT
            ),
        )

        medians = df.groupby(DataFrameStrings.EXPERIMENT).median()

        df_plot = df.reset_index()

        if plot_type == "strip":
            plot = px.strip(
                df_plot,
                color=DataFrameStrings.EXPERIMENT,
                stripmode="overlay",
                **plotargs,
            )
            plot.update_traces(width=2.3)
            plot.update_xaxes(range=(-0.6, len(multi_choice) - 0.4))
            for index, m in medians.iterrows():
                i = multi_choice.index(index)
                plot.add_shape(
                    x0=i - 0.45,
                    x1=i + 0.45,
                    y0=m[0],
                    y1=m[0],
                    line_color=px.colors.DEFAULT_PLOTLY_COLORS[i],
                    line_width=3,
                    opacity=0.8,
                )
        elif plot_type == "box":
            plot = px.box(df_plot, color=DataFrameStrings.EXPERIMENT, **plotargs)
        elif plot_type == "violin":
            plot = px.violin(df_plot, color=DataFrameStrings.EXPERIMENT, **plotargs)
        elif plot_type == "stacked":
            plot = px.bar(
                df_plot,
                color=DataFrameStrings.CLUSTER,
                barmode="stack",
                **plotargs,
            )
            plot.update_layout(
                legend_traceorder="reversed",
                height=30 * len(df_plot[DataFrameStrings.CLUSTER].unique()),
            ).update_traces(marker_line_color="black", marker_line_width=1)
        elif plot_type == "histogram":
            plot = px.histogram(
                df_plot,
                x="distance",
                color=DataFrameStrings.EXPERIMENT,
                barmode="overlay",
                **{k: v for k, v in plotargs.items() if k != "x"},
            )
        else:
            raise ValueError("Unknown plot type")

        if (
            highlight is not None
            and highlight in df.index.get_level_values(DataFrameStrings.CLUSTER)
            and plot_type not in ["histogram", "stacked"]
        ):
            df_highlight = df.xs(
                highlight, level=DataFrameStrings.CLUSTER, axis=0, drop_level=False
            )
            for index, v in df_highlight.iterrows():
                i = multi_choice.index(
                    index[df.index.names.index(DataFrameStrings.EXPERIMENT)]
                )
                plot.add_annotation(
                    x=i - 0.45,
                    y=v[0],
                    ax=-10,
                    ay=0,
                    showarrow=True,
                    arrowside="end",
                    arrowhead=1,
                    arrowwidth=2,
                    arrowsize=1,
                    hovertext="<br>".join(index),
                )

        plot.update_layout(**calc_width_categorical(plot), xaxis_tickangle=45)

        return medians, plot

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
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined[
            df_quantity_pr_pg_combined[DataFrameStrings.EXPERIMENT].isin(multi_choice)
        ]

        df_quantity_pr_pg_combined.insert(
            0,
            "Expxfiltering",
            [
                " ".join([e, f])
                for e, f in zip(
                    df_quantity_pr_pg_combined.Experiment,
                    df_quantity_pr_pg_combined.filtering,
                )
            ],
        )
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined.assign(
            Experiment_lexicographic_sort=pd.Categorical(
                df_quantity_pr_pg_combined[DataFrameStrings.EXPERIMENT],
                categories=multi_choice,
                ordered=True,
            )
        )

        df_quantity_pr_pg_combined.sort_values(
            ["Experiment_lexicographic_sort", "type"],
            ascending=[True, False],
            inplace=True,
        )

        layout = go.Layout(
            barmode="overlay",
            # xaxis_tickangle=90,
            autosize=False,
            width=100 * len(multi_choice) + 150,
            height=400,
            template="simple_white",
        )
        filtered = list(np.tile(["id", "profile"], len(multi_choice)))

        fig_quantity_pg = px.bar(
            df_quantity_pr_pg_combined,
            x="Expxfiltering",
            y="number of protein groups",
            color=DataFrameStrings.EXPERIMENT,
            barmode="overlay",
            hover_data=["type"],
            opacity=0.8,
            color_discrete_sequence=px.colors.qualitative.D3,
        )
        fig_quantity_pg.update_layout(
            layout,
            title="Number of Protein Groups",
            xaxis={
                "tickmode": "array",
                "tickvals": [el for el in range(len(multi_choice) * 2)],
                "ticktext": filtered,
                "title": {"text": None},
            },
        )
        fig_quantity_pg.update_layout(**calc_width_categorical(fig_quantity_pg))

        fig_quantity_pr = px.bar(
            df_quantity_pr_pg_combined,
            x="filtering",
            y="number of profiles",
            color="type",
            barmode="overlay",
            labels={DataFrameStrings.EXPERIMENT: "", "filtering": ""},
            facet_col=DataFrameStrings.EXPERIMENT,
            template="simple_white",
            opacity=1,
        ).for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig_quantity_pr.update_layout(layout, title="Number of Profiles")
        fig_quantity_pr.update_layout(**calc_width_categorical(fig_quantity_pr))

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
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined[
            df_quantity_pr_pg_combined[DataFrameStrings.EXPERIMENT].isin(multi_choice)
        ]

        df_quantity_pr_pg_combined.insert(
            0,
            "Expxfiltering",
            [
                " ".join([e, f])
                for e, f in zip(
                    df_quantity_pr_pg_combined.Experiment,
                    df_quantity_pr_pg_combined.filtering,
                )
            ],
        )
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined.assign(
            Experiment_lexicographic_sort=pd.Categorical(
                df_quantity_pr_pg_combined[DataFrameStrings.EXPERIMENT],
                categories=multi_choice,
                ordered=True,
            )
        )

        # df_quantity_pr_pg_combined.sort_values("Experiment_lexicographic_sort", inplace=True)
        df_quantity_pr_pg_combined.sort_values(
            ["Experiment_lexicographic_sort", "filtering"],
            ascending=[True, False],
            inplace=True,
        )
        filtered = list(np.tile(["id", "profile"], len(multi_choice)))

        fig_pr_dc = px.bar(
            df_quantity_pr_pg_combined.loc[df_quantity_pr_pg_combined.type == "total"],
            x="Expxfiltering",
            y="data completeness of profiles",
            color=DataFrameStrings.EXPERIMENT,
            barmode="overlay",
            hover_data=["filtering"],
            template="simple_white",
            opacity=0.8,
            color_discrete_sequence=px.colors.qualitative.D3,
        )

        fig_pr_dc.update_layout(
            title="Profile completeness of all<br>identified protein groups",
            autosize=False,
            width=100 * len(multi_choice) + 150,
            height=400,
            template="simple_white",
            xaxis={
                "tickmode": "array",
                "tickvals": [el for el in range(len(multi_choice) * 2)],
                "ticktext": filtered,
                "title": {"text": None},
            },
        )
        fig_pr_dc.update_layout(**calc_width_categorical(fig_pr_dc))

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
            bit_flags = [2**n for n in range(len(sets))]
            flags_zip_sets = [z for z in zip(bit_flags, sets)]
            combo_sets = []
            overlapping_ids = []
            experiments = []
            # dictio = {}
            for bits in range(num_combinations - 1, 0, -1):
                include_sets = [s for flag, s in flags_zip_sets if bits & flag]
                exclude_sets = [s for flag, s in flags_zip_sets if not bits & flag]
                combo = set.intersection(*include_sets)
                combo = set.difference(combo, *exclude_sets)
                tag = "".join([str(int((bits & flag) > 0)) for flag in bit_flags])

                experiment_decoded = []
                for digit, exp in zip(list(tag), multi_choice):
                    if digit == "0":
                        continue
                    else:
                        experiment_decoded.append(exp)
                # dictio[len(combo)] = experiment_decoded
                if len(multi_choice) > 3:
                    if len(combo) > omit_size:
                        overlapping_ids.append(len(combo))
                        experiments.append(experiment_decoded)
                else:
                    if len(combo) > 0:
                        overlapping_ids.append(len(combo))
                        experiments.append(experiment_decoded)
                # combo_sets.append((tag, len(combo)))

            fig_UpSetPlot = plt.Figure()
            series_UpSetPlot = from_memberships(experiments, data=overlapping_ids)
            upplot(series_UpSetPlot, fig=fig_UpSetPlot, show_counts="%d")
            return fig_UpSetPlot

        if "Sequence" not in self.df_01_filtered_combined.index.names:
            sets_proteins_total = [
                set(
                    self.df_01_filtered_combined.xs(
                        i, axis=0, level=DataFrameStrings.EXPERIMENT
                    ).index.get_level_values(DataFrameStrings.PROTEIN_IDS)
                )
                for i in multi_choice_venn
            ]
            sets_proteins_intersection = [
                set(
                    self.df_01_filtered_combined.xs(
                        i, axis=0, level=DataFrameStrings.EXPERIMENT
                    )
                    .unstack([DataFrameStrings.MAP, DataFrameStrings.EXP_MAP])
                    .dropna()
                    .index.get_level_values(DataFrameStrings.PROTEIN_IDS)
                )
                for i in multi_choice_venn
            ]
        else:
            sets_proteins_total = [
                set(
                    self.df_01_filtered_combined.xs(
                        i, axis=0, level=DataFrameStrings.EXPERIMENT
                    ).index.get_level_values("Sequence")
                )
                for i in multi_choice_venn
            ]
            sets_proteins_intersection = [
                set(
                    self.df_01_filtered_combined.xs(
                        i, axis=0, level=DataFrameStrings.EXPERIMENT
                    )
                    .unstack([DataFrameStrings.MAP, DataFrameStrings.EXP_MAP])
                    .dropna()
                    .index.get_level_values("Sequence")
                )
                for i in multi_choice_venn
            ]
        figure_UpSetPlot_total = create_upsetplot(
            sets_proteins_total, multi_choice_venn
        )
        figure_UpSetPlot_int = create_upsetplot(
            sets_proteins_intersection, multi_choice_venn
        )

        # make matplot figure available for plotly
        def convert_venn_jpg(vd):
            vd = vd.figure
            out_img = BytesIO()
            plt.savefig(out_img, bbox_inches="tight", format="jpg", dpi=72)
            out_img.seek(0)  # rewind file
            im = Image.open(out_img)
            plt.clf()
            return im

        if len(multi_choice_venn) == 2:
            vd_t = venn2(
                sets_proteins_total,
                set_labels=([i for i in multi_choice_venn]),
                set_colors=px.colors.qualitative.D3[0:2],
                alpha=0.8,
            )
            vd_t = plt.title("in at least one map")
            im_t = convert_venn_jpg(vd_t)
            vd_i = venn2(
                sets_proteins_intersection,
                set_labels=([i for i in multi_choice_venn]),
                set_colors=px.colors.qualitative.D3[0:2],
                alpha=0.8,
            )
            vd_i = plt.title("in all maps")
            im_i = convert_venn_jpg(vd_i)
        elif len(multi_choice_venn) == 3:
            vd_t = venn3(
                sets_proteins_total,
                set_labels=([i for i in multi_choice_venn]),
                set_colors=px.colors.qualitative.D3[0:3],
                alpha=0.8,
            )
            vd_t = plt.title("in at least one map")
            im_t = convert_venn_jpg(vd_t)
            vd_i = venn3(
                sets_proteins_intersection,
                set_labels=([i for i in multi_choice_venn]),
                set_colors=px.colors.qualitative.D3[0:3],
                alpha=0.8,
            )
            vd_i = plt.title("in all maps")
            im_i = convert_venn_jpg(vd_i)

        else:
            im = "Venn diagram can be displayed for 3 Experiments or less"
            return im, im, figure_UpSetPlot_total, figure_UpSetPlot_int

        return im_t, im_i, figure_UpSetPlot_total, figure_UpSetPlot_int

    def calculate_global_scatter(
        self, metric="manhattan distance to average profile", consolidation="average"
    ):
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
        cons_functions = {"median": np.median, "average": np.mean, "sum": np.sum}
        metrics = {
            "euclidean distance": "euclidean",
            "manhattan distance": "manhattan",
            "1 - cosine correlation": "cosine",
            "1 - pearson correlation": lambda x, y: 1 - np.corrcoef(x, y)[0][1],
            "manhattan distance to average profile": [
                np.mean,
                pw.paired_manhattan_distances,
            ],
            "manhattan distance to median profile": [
                np.median,
                pw.paired_manhattan_distances,
            ],
        }

        # Option assertion
        assert consolidation in cons_functions.keys()
        assert metric in metrics.keys()

        self.parameters["reproducibility metric"] = metric
        self.parameters["reproducibility consolidation"] = consolidation

        df = self.df_01_filtered_combined.copy()
        df.index = df.index.droplevel(
            [
                DataFrameStrings.EXP_MAP,
                DataFrameStrings.GENE_NAMES,
                DataFrameStrings.COMPARTMENT,
            ]
        )
        if "Sequence" in df.index.names:
            df.index = df.index.droplevel([DataFrameStrings.PROTEIN_IDS])

        # Calculate and consolidate distances
        distances = pd.DataFrame()
        for exp in self.exp_names:
            df_m = (
                df.xs(exp, level=DataFrameStrings.EXPERIMENT, axis=0)
                .unstack(DataFrameStrings.MAP)
                .dropna(axis=1, how="all")
                .dropna()
                .stack(DataFrameStrings.MAP)
            )
            maps = list(set(df_m.index.get_level_values(DataFrameStrings.MAP)))

            # this if clause switches between pairwise comparisons of profiles (else) and comparisons to an average/median profile
            if " to " in metric:
                df_m = df_m.unstack(DataFrameStrings.MAP)

                # calculate reference profiles
                df_profiles = (
                    df_m.stack(DataFrameStrings.FRACTION)
                    .apply(metrics[metric][0], axis=1)
                    .unstack(DataFrameStrings.FRACTION)
                )

                # calculate the distance for every map
                distances_m = pd.DataFrame()
                for m in maps:
                    dist_m = pd.DataFrame(
                        metrics[metric][1](
                            df_m.xs(m, level=DataFrameStrings.MAP, axis=1), df_profiles
                        ),
                        columns=[m],
                    )
                    distances_m = pd.concat([distances_m, dist_m], axis=1)

                distances_m.index = df_m.index

            else:
                distances_m = pd.DataFrame()

                # loop over pairs of maps
                for i, mapi in enumerate(maps):
                    for j, mapj in enumerate(maps):
                        # only look at each comparison once
                        if j <= i:
                            continue
                        dist = pw.paired_distances(
                            df_m.xs(mapi, level=DataFrameStrings.MAP, axis=0).values,
                            df_m.xs(mapj, level=DataFrameStrings.MAP, axis=0).values,
                            metric=metrics[metric],
                        )
                        dist = pd.Series(dist, name="_".join([mapi, mapj]))
                        distances_m = pd.concat([distances_m, dist], axis=1)
                distances_m.index = df_m.xs(
                    maps[0], level=DataFrameStrings.MAP, axis=0
                ).index

            distances = pd.concat(
                [
                    distances,
                    pd.Series(
                        distances_m.apply(cons_functions[consolidation], axis=1),
                        name=exp,
                        index=distances_m.index,
                    ),
                ],
                axis=1,
            )

        self.distances = distances

    def plot_reproducibility_distribution(
        self, multi_choice=[], q=0.5, show_full=True, show_rug=False, x_cut=None
    ):
        if len(multi_choice) == 0:
            multi_choice = self.exp_names

        distances = self.distances[multi_choice].copy()
        if x_cut == None:
            x_cut = 1.05 * distances.max().max()

        # get data and labels set up for figure factory
        if show_full:
            plotdata = distances.join(
                distances.dropna(), lsuffix=" full", rsuffix=" overlap"
            )
        else:
            plotdata = distances.dropna()
        plotlabels = list(plotdata.columns)
        plotdata = [vals.dropna() for k, vals in plotdata.T.iterrows()]

        # calculate quantiles
        quantiles = dict()
        for i, (name, data) in enumerate(zip(plotlabels, plotdata)):
            quantiles[name] = np.quantile(data, q)

        # create plot
        plot = ff.create_distplot(
            plotdata, plotlabels, show_hist=False, show_rug=show_rug, show_curve=True
        )
        plot.update_layout(
            title="Distribution of {} {}, overlap = {}".format(
                self.parameters["reproducibility consolidation"],
                self.parameters["reproducibility metric"],
                len(plotdata[-1]),
            ),
            width=1000,
            height=600,
            template="simple_white",
            xaxis={"range": (0, x_cut), "rangemode": "nonnegative"},
        )

        # get color dictionary and assign same color to overlap and full coverage traces, assign legend groups
        colors = []
        plot.for_each_trace(
            lambda x: colors.append(x.marker.color),
            selector=lambda x: not x.name.endswith(" overlap"),
        )
        if show_full:
            colors[distances.shape[1] : :] = colors[: distances.shape[1]]
        colors = {k: v for k, v in zip(plotlabels, colors)}
        plot.for_each_trace(
            lambda x: x.update(legendgroup=x.name, marker_color=colors[x.name])
        )

        # get density estimates at quantile and add highlighting lines
        densities = []
        plot.for_each_trace(
            lambda x: densities.append(x.y[bisect.bisect(x.x, quantiles[x.name])]),
            selector=dict(type="scatter"),
        )
        for i, (name, data) in enumerate(zip(plotlabels, plotdata)):
            plot.add_scatter(
                x=[quantiles[name], quantiles[name]],
                y=[0, densities[i]],
                text=["", "{:.4f}".format(quantiles[name])],
                name=name,
                legendgroup=name,
                showlegend=False,
                line_color=colors[name],
                mode="lines+text",
                textposition="top right",
            )

        # dash full distribution traces
        if show_full:
            plot.for_each_trace(
                lambda x: x.update(line_dash="dash"),
                selector=lambda x: not x.name.endswith(" overlap"),
            )

        return plot

    def train_svm(
        self,
        experiments,
        random_state,
        test_split,
        classes,
        output,
        canvas,
        rounds,
        C0,
        C1,
        g0,
        g1,
        overwrite=False,
    ):
        hash_data = SVMComp._construct_hash_data(
            experiments=experiments,
            random_state=random_state,
            test_split=test_split,
            classes=classes,
        )
        if hash_data in self.svm_runs.keys() and overwrite == False:
            raise RuntimeError("Confirm retraining")

        svmcomp = SVMComp(hash_data)
        svmcomp.set_df(self.df_01_filtered_combined)
        svmcomp.train_svm(
            C0=C0, C1=C1, g0=g0, g1=g1, canvas=canvas, output=output, rounds=rounds
        )
        svmcomp.df = pd.DataFrame()
        self.svm_runs[hash_data] = svmcomp

        return svmcomp.C, svmcomp.gamma

    def predict_svm(
        self,
        experiments,
        random_state,
        test_split,
        classes,
        C,
        gamma,
        min_p,
        min_diff,
        svmset,
    ):
        hash_data = SVMComp._construct_hash_data(
            experiments=experiments,
            random_state=random_state,
            test_split=test_split,
            classes=classes,
        )
        if hash_data in self.svm_runs.keys():
            svmcomp = self.svm_runs[hash_data]
        else:
            svmcomp = SVMComp(hash_data)
        svmcomp.set_df(self.df_01_filtered_combined)
        svmcomp.set_parameters(C=C, gamma=gamma)
        svmcomp._predict_svm()
        svmcomp._score_predictions(min_p=min_p, min_diff=min_diff)
        for exp in svmcomp.experiments:
            self.add_svm_result(
                misclassification=svmcomp._get_confusion(exp),
                source="direct",
                experiment=exp,
                comment=svmcomp.hash_data + "___" + svmcomp.hash_parameters,
                name=svmset,
                prediction=svmcomp.predictions[exp]
                .join(svmcomp.probabilities[exp])
                .dropna(),
            )
        svmcomp.df = pd.DataFrame()
        self.svm_runs[hash_data] = svmcomp
        return svmcomp.predictions

    def plot_overview(self, multi_choice, clusters, quantile):
        dists, _ = self.plot_intramap_scatter(
            normalization=np.median,
            min_size=1,
            multi_choice=multi_choice,
            clusters_for_ranking=clusters,
        )
        dists = dists.rename({"distance": "complex scatter"}, axis=1)

        rep = (
            self.distances[multi_choice]
            .dropna()
            .apply(lambda x: np.quantile(x, quantile), axis=0)
        )
        rep.name = "intermap scatter"

        depth = (
            self.df_quantity_pr_pg_combined.query('type == "intersection"')
            .query('filtering == "after filtering"')
            .query("Experiment in @multi_choice")[
                [DataFrameStrings.EXPERIMENT, "number of protein groups"]
            ]
            .set_index(DataFrameStrings.EXPERIMENT)
            .rename({"number of protein groups": "profiled depth"}, axis=1)
        )

        df_plot = dists.join(rep).join(depth).melt(ignore_index=False).reset_index()

        fig = make_subplots(
            1,
            3,
            subplot_titles=["profiled depth", "complex scatter", "intermap scatter"],
            horizontal_spacing=0.17,
        )

        for i, y in enumerate(
            ["profiled depth", "complex scatter", "intermap scatter"]
        ):
            for j, e in enumerate(df_plot.Experiment.unique()):
                fig.add_trace(
                    go.Bar(
                        x=df_plot.query(
                            "variable == @y and Experiment == @e"
                        ).Experiment,
                        y=df_plot.query("variable == @y and Experiment == @e").value,
                        name=e,
                        marker=dict(
                            color=px.colors.qualitative.D3[j], coloraxis="coloraxis"
                        ),
                        legendgroup=e,
                        hovertemplate="%{y:.5r}",
                    ),
                    1,
                    i + 1,
                )
        fig.update_xaxes(matches="x").update_layout(
            template="simple_white", showlegend=False
        )
        fig.for_each_yaxis(
            lambda x: x.update(
                title="full coverage protein groups"
                if x.anchor == "x"
                else "median normalized intracomplex scatter"
                if x.anchor == "x2"
                else f"{quantile * 100}% quantile of shared protein groups",
                tickformat=".0s" if x.anchor == "x" else ".2r",
                nticks=8,
                title_standoff=8 if x.anchor == "x" else 0,
            )
        )
        fig.update_layout(
            width=320 + (200 * (j + 1)),
            height=500,
            margin=dict(t=70, b=0, r=0, l=0),
            title="Overview of benchmarking output",
        )
        return fig

    def plot_svm_summary(
        self, svmset="default", multi_choice=[], score="F1 score", orientation="v"
    ):
        if multi_choice == []:
            multi_choice = self.exp_names

        try:
            df = self.df_svm_performance[score].xs(
                svmset, level=DataFrameStrings.SET, axis=0
            )
        except KeyError:
            raise KeyError(
                f"{score} is not among the criteria calculated for classification performance."
            )

        df = df.unstack("value").reset_index()
        df = df.loc[df.Experiment.isin(multi_choice)]
        df = df.loc[
            df.Compartment.isin(
                [
                    "Overall performance",
                    "Average all classes",
                    "Average membraneous organelles",
                ]
            )
        ]
        df.rename({DataFrameStrings.COMPARTMENT: "Aggregation"}, axis=1, inplace=True)
        df = df.sort_values(
            DataFrameStrings.EXPERIMENT,
            key=lambda x: [multi_choice.index(el) for el in x],
        )

        if score != "Class size":
            plot = px.bar(
                df,
                x="Aggregation" if orientation == "v" else "mean",
                y="mean" if orientation == "v" else "Aggregation",
                error_x=None if orientation == "v" else "std",
                error_y="std" if orientation == "v" else None,
                orientation=orientation,
                color=DataFrameStrings.EXPERIMENT,
                template="simple_white",
                barmode="group",
            ).update_yaxes(title=score)
        else:
            plot = px.bar(
                df,
                x="Aggregation" if orientation == "v" else "direct",
                y="direct" if orientation == "v" else "Aggregation",
                orientation=orientation,
                color=DataFrameStrings.EXPERIMENT,
                template="simple_white",
                barmode="group",
            ).update_yaxes(title=score)

        return plot

    def plot_svm_detail(
        self, svmset="default", multi_choice=[], score="F1 score", orientation="v"
    ):
        if multi_choice == []:
            multi_choice = self.exp_names

        try:
            df = (
                self.df_svm_performance[score]
                .xs(svmset, level=DataFrameStrings.SET, axis=0)
                .drop(
                    [
                        "Overall performance",
                        "Average all classes",
                        "Average membraneous organelles",
                    ],
                    level=DataFrameStrings.COMPARTMENT,
                    axis=0,
                )
            )
        except KeyError:
            raise KeyError(
                f"{score} is not among the criteria calculated for classification performance."
            )

        df = df.unstack("value").reset_index()
        df = df.loc[df.Experiment.isin(multi_choice)]
        df = df.sort_values(
            DataFrameStrings.EXPERIMENT,
            key=lambda x: [multi_choice.index(el) for el in x],
        )
        if score != "Class size":
            plot = px.bar(
                df,
                x=DataFrameStrings.COMPARTMENT if orientation == "v" else "mean",
                y="mean" if orientation == "v" else DataFrameStrings.COMPARTMENT,
                orientation=orientation,
                error_x=None if orientation == "v" else "std",
                error_y="std" if orientation == "v" else None,
                color=DataFrameStrings.EXPERIMENT,
                template="simple_white",
                barmode="group",
            ).update_yaxes(title=score)
        else:
            plot = px.bar(
                df,
                x=DataFrameStrings.COMPARTMENT if orientation == "v" else "direct",
                y="direct" if orientation == "v" else DataFrameStrings.COMPARTMENT,
                orientation=orientation,
                color=DataFrameStrings.EXPERIMENT,
                template="simple_white",
                barmode="group",
            ).update_yaxes(title=score)

        return plot

    def __repr__(self):
        return str(self.__dict__)


class SVMComp:
    """
    This is a class meant to operate on the data from a SpatialDataSetComparison object and holds all functionality for configuring and training SVMs for orgenellar assignments.
    """

    def __init__(self, hash_data, **kwargs):
        """
        hash_data can be constructed using _construct_hash_data.
        """
        # extract main parameters from the data hash
        self.hash_data = hash_data
        self.experiments = hash_data.split("__")[0].split(";")
        self.random_state = int(hash_data.split("__")[1])
        self.test_split = float(hash_data.split("__")[2])
        self.classes = hash_data.split("__")[3].split(";")

        # If accuracies have oreviously determined and are passed, set them.
        if "accuracies" in kwargs.keys():
            self.accuracies = kwargs["accuracies"]
        else:
            self.accuracies = pd.DataFrame(
                columns=pd.MultiIndex.from_tuples(
                    [
                        (el, s)
                        for el in self.experiments
                        for s in ["C", "gamma", "Accuracy"]
                    ],
                    names=[DataFrameStrings.EXPERIMENT, "Score"],
                )
            )

        # Currently not used, but additional interfaceparameters like the configuration for training could be saved here.
        if "interfaceparameters" in kwargs.keys():
            self.interfaceparameters = kwargs["interfaceparameters"]
        else:
            self.interfaceparameters = dict()

        # Set default parameters C and gamma
        self.set_parameters(C=5, gamma=25)

        # Initialize empty dataframes for probabilities, predictions and input data
        self.probabilities, self.predictions, self.df = (
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
        )

    def _construct_hash_data(
        experiments: list,
        random_state: int = 42,
        test_split: float = 0.2,
        classes: list = [],
    ):
        """
        This constructs the data hash for a SVMComp class, which can then be used to initialize it.
        """
        return "__".join(
            [
                ";".join(experiments),
                str(random_state),
                str(test_split),
                ";".join(classes),
            ]
        )

    def _construct_hash_parameters(C: float = 5, gamma: float = 25):
        """
        This constructs the parameter hash for a SVMComp class, which can be used to uniquely identify the SVM
        """
        return "__".join([str(C), str(gamma)])

    def set_df(self, df_01_filtered_combined):
        """
        Reduce the input dataframe from the SpatialDatasetComparison class to the required experiments.
        """
        df = df_01_filtered_combined.reset_index(
            [
                el
                for el in df_01_filtered_combined.index.names
                if el
                not in [
                    DataFrameStrings.EXPERIMENT,
                    DataFrameStrings.MAP,
                    DataFrameStrings.PROTEIN_IDS,
                    DataFrameStrings.COMPARTMENT,
                ]
            ],
            drop=True,
        )
        df = df[
            df.index.get_level_values(DataFrameStrings.EXPERIMENT).isin(
                self.experiments
            )
        ]
        self.df = df

    def set_parameters(self, C: float = 5, gamma: float = 25):
        """
        Set parameters for the SVM and store the constructed hash
        """
        self.C = C
        self.gamma = gamma
        self.hash_parameters = SVMComp._construct_hash_parameters(C=C, gamma=gamma)

    def __repr__(self):
        text = """This is an SVM object for comparing several spatial datasets.

 Experiments:  {}
 Test split:   {}
 Random state: {}
 Classes:      {}
 
 C:     {}
 gamma: {}
 
 DataFrame dimensions are:
 Input data:          {}
 Grid accuracies:     {}
 Class probabilities: {}
 Class assignments:   {}""".format(
            ", ".join(self.experiments),
            self.test_split,
            self.random_state,
            ", ".join(self.classes),
            self.C,
            self.gamma,
            " x ".join([str(el) for el in self.df.shape]),
            " x ".join([str(el) for el in self.accuracies.shape]),
            " x ".join([str(el) for el in self.probabilities.shape]),
            " x ".join([str(el) for el in self.predictions.shape]),
        )
        return text

    def _get_markers(self, which="shared"):
        """
        Retrieve marker subsets across experiments, for classes specified in the classes attribute.
        Arguments:
            which: one of shared, unshared
        """

        if which == "shared":
            shared = (
                self.df.unstack([DataFrameStrings.EXPERIMENT, DataFrameStrings.MAP])
                .dropna()
                .index.to_frame(index=False)[
                    [DataFrameStrings.PROTEIN_IDS, DataFrameStrings.COMPARTMENT]
                ]
            )
            return shared[shared[DataFrameStrings.COMPARTMENT].isin(self.classes)]
        else:
            unshared = (
                self.df.unstack([DataFrameStrings.EXPERIMENT, DataFrameStrings.MAP])
                .drop(
                    self.df_01_filtered_combined.unstack(
                        [DataFrameStrings.EXPERIMENT, DataFrameStrings.MAP]
                    )
                    .dropna()
                    .index,
                    axis=0,
                )
                .index.to_frame(index=False)[
                    [DataFrameStrings.PROTEIN_IDS, DataFrameStrings.COMPARTMENT]
                ]
            )
            return unshared[unshared[DataFrameStrings.COMPARTMENT].isin(self.classes)]

    def _train_test_split(
        self, test_unshared=False, test_percent=None, random_state=None
    ):
        """

        Arguments:
            test_unshared: Whether to include unshared markers in the test sets
            test_percent: Percentage of shared markers to set aside as a testset
            random_state: Random state passed to train_test_split from sklearn

        Returns:
            train: pd.DataFrame, columns = ["Protein IDs", "Comaprtment"]
            test: pd.DataFrame, columns = ["Protein IDs", "Comaprtment"]
        """

        shared = self._get_markers(which="shared")
        if test_percent is None:
            test_percent = self.test_split
        if random_state is None:
            random_state = self.random_state

        if test_percent != 0:
            train, test = train_test_split(
                shared,
                test_size=test_percent,
                random_state=random_state,
                stratify=shared[DataFrameStrings.COMPARTMENT],
            )
        else:
            train = shared
            test = shared

        if test_unshared == True:
            test = pd.concat(
                [test, self._get_markers(which="unshared")], axis=0, ignore_index=True
            )

        return train, test

    def _get_grid(C0=1, C1=60, Cn=4, g0=0.01, g1=60, gn=4):
        """
        Return dictionary suitable as input for GridSearchCV, using np.geomspace for spacing.
        """
        return {
            "C": np.geomspace(C0, C1, Cn).round(5),
            "gamma": np.geomspace(g0, g1, gn).round(5),
        }

    def _grid_search(
        self,
        train=pd.DataFrame(
            columns=[DataFrameStrings.PROTEIN_IDS, DataFrameStrings.COMPARTMENT]
        ),
        parameters={"C": np.logspace(0, 1.7, 4), "gamma": np.logspace(0, 2, 4)},
        random_state=None,
    ):
        """
        This runs the specified gridsearch for each experiment and returns the best parameters per experiment and all accuracies.
        Arguments:
            train: training ids and labels
            parameters: search space for C and gamma
        """
        accuracies = dict()
        best_params = dict()

        if random_state == None:
            random_state = self.random_state

        labels = train[DataFrameStrings.COMPARTMENT]
        svc = svm.SVC(probability=False, random_state=random_state)

        for exp in self.experiments:
            # Only use full coverage data and retrieve training set
            data = (
                self.df.xs(exp, level=DataFrameStrings.EXPERIMENT, axis=0)
                .unstack(DataFrameStrings.MAP)
                .dropna()
            )
            data = data.reset_index(
                level=[DataFrameStrings.COMPARTMENT], drop=True
            ).loc[train[DataFrameStrings.PROTEIN_IDS], :]

            # COnfigure and run gridsearch
            clf = GridSearchCV(svc, parameters)
            clf.fit(data.values, labels)

            # Store results
            best_params[exp] = clf.best_params_
            accuracies[exp] = pd.DataFrame(
                [
                    [k["C"], k["gamma"], v]
                    for k, v in zip(
                        clf.cv_results_["params"], clf.cv_results_["mean_test_score"]
                    )
                ],
                columns=["C", "gamma", "accuracy"],
            ).set_index(["C", "gamma"])

        accuracies = pd.concat(accuracies, axis=1)
        accuracies.columns.names = [DataFrameStrings.EXPERIMENT, "Value"]
        return best_params, accuracies

    def _best_estimator(
        accuracies=pd.DataFrame(
            columns=pd.MultiIndex.from_tuples(
                [(".*", "accuracy")], names=[DataFrameStrings.EXPERIMENT, "Value"]
            ),
            index=pd.MultiIndex.from_tuples([], names=["C", "gamma"]),
        ),
        mode="maxsum",
    ):
        """
        Determine the best estimator from an accuracies table, as returned by _grid_search.
        Only valid mode right now is 'maxsum'.
        """

        if mode == "maxsum":
            sums = accuracies.sum(axis=1)
            return accuracies.index[np.where(sums == max(sums))].to_list()[0]

    def _plot_training(accuracies):
        accuracies = accuracies.stack([DataFrameStrings.EXPERIMENT]).reset_index()
        return px.scatter_3d(
            accuracies,
            x="C",
            y="gamma",
            color=DataFrameStrings.EXPERIMENT,
            z="accuracy",
            log_x=True,
            log_y=True,
        )

    def train_svm(
        self,
        train=pd.DataFrame(
            columns=[DataFrameStrings.PROTEIN_IDS, DataFrameStrings.COMPARTMENT]
        ),
        rounds=5,
        autoselect="maxsum",
        output="print",
        canvas=None,
        C0=1,
        C1=30,
        g0=1,
        g1=50,
    ):
        """
        To find the best SVM parameters C and gamma several rounds of a gridsearch are run, employing 5 fold training validation.
        In the first round the initially given range (default C[1,30], gamma[1,50]) is searched with 36 grid points.
        For every next round the grid is determined based on the minimum and maximum best parameters across the experiments.
        How the grid is expanded depends on the round and the agreement (agreed := (max-min)/max < 0.2) of the experiments.
        Rounds 2 and 3:
            agreement -> 5 points sampled from [0.7*min,max/0.7]
            disagreement -> 6 points sampled from [min-0.5*(max-min),max+0.5(max-min)]
            During these round the best parameters are only taken from the newly added grid points, to ensure sampling of many different points first.
        Round 4+:
            agreement -> 3 points sampled from [0.9*min,max/0.9]
            disagreement -> 4 points sampled from [min-0.3*(max-min),max+0.3(max-min)]
            The best parameters are now taken from all gridpoints sampled in any round.
            If the calculated grid has been sampled before (because the best parameters didn't change), the limits are slightly expanded and more poitns are sampled (up to 6).

        If all expriments agree after at least 4 rounds, training is stopped early.

        In a final step the average best and the overall best parameters (maxizing the summed accuracy) are compared.
        """
        if len(train) == 0:
            train = self._train_test_split(
                test_percent=self.test_split, random_state=self.random_state
            )[0]

        n = 0
        best_params = pd.Series()
        accuracies = pd.DataFrame()
        while rounds > n:
            if n == 0:
                grid = SVMComp._get_grid(C0=C0, C1=C1, Cn=6, g0=g0, g1=g1, gn=6)
            else:
                if n <= 2:
                    agreed_expansion = 0.7
                    disagreed_expansion = 0.5
                    n_agreed = 5
                    n_disagreed = 6
                if n > 2:
                    agreed_expansion = 0.9
                    disagreed_expansion = 0.3
                    n_agreed = 3
                    n_disagreed = 4

                c0 = min([el["C"] for el in best_params.values()])
                c1 = max([el["C"] for el in best_params.values()])
                if (c1 - c0) / c1 < 0.2:
                    c0 = agreed_expansion * c0
                    c1 = min(c1 / agreed_expansion, 95)
                    cn = n_agreed
                else:
                    c0 = max(c0 - disagreed_expansion * (c1 - c0), 0.01)
                    c1 = min(c1 + disagreed_expansion * (c1 - c0), 95)
                    cn = n_disagreed

                g0 = min([el["gamma"] for el in best_params.values()])
                g1 = max([el["gamma"] for el in best_params.values()])
                if (g1 - g0) / g1 < 0.2:
                    g0 = agreed_expansion * g0
                    g1 = min(g1 / agreed_expansion, 500)
                    gn = n_agreed
                else:
                    g0 = max(g0 - disagreed_expansion * (g1 - g0), 0.0001)
                    g1 = min(g1 + disagreed_expansion * (g1 - g0), 500)
                    gn = n_disagreed

                grid_new = SVMComp._get_grid(Cn=cn, gn=gn, C0=c0, C1=c1, g0=g0, g1=g1)
                ng = 1
                while all(
                    [
                        (ii, jj) in list(accuracies.index.values)
                        for i, j in zip(*np.meshgrid(*grid_new.values()))
                        for ii, jj in zip(i, j)
                    ]
                ):
                    if output == "print":
                        print(
                            "Calculated grid was already searched. Resampling from larger range..."
                        )
                    elif output == "canvas":
                        canvas.append(
                            "Calculated grid was already searched. Resampling from larger range..."
                        )
                    grid_new = SVMComp._get_grid(
                        Cn=min(cn + ng, 6),
                        gn=min(gn + ng, 6),
                        C0=c0 * (1 - ng * 0.1),
                        C1=c1 / (1 - ng * 0.1),
                        g0=g0 * (1 - ng * 0.1),
                        g1=g1 / (1 - ng * 0.1),
                    )
                    ng += 1

                grid = grid_new

            if n > 3:
                if gn == n_agreed and cn == n_agreed:
                    n = rounds

            # Actual training
            best_params, accuracies_n = self._grid_search(train=train, parameters=grid)

            # Store results
            if n == 0:
                accuracies = accuracies_n
                if output == "canvas":
                    canvas[0] = SVMComp._plot_training(accuracies)
            else:
                accuracies = pd.concat([accuracies, accuracies_n])
                accuracies = accuracies[~accuracies.index.duplicated(keep="first")]
                if output == "canvas":
                    canvas[0] = SVMComp._plot_training(accuracies)
                if n > 2:
                    for k in best_params.keys():
                        best_params[k] = {
                            n: v
                            for n, v in zip(
                                accuracies.index.names, accuracies[k].idxmax(axis=0)[0]
                            )
                        }
            best_scoring = SVMComp._best_estimator(accuracies, mode=autoselect)
            if output == "print" and n <= 2:
                print(
                    f"Round {str(n + 1)} gave the following optimal parameters:\n{pd.DataFrame(best_params)},\nfrom this grid:\nC    :{grid['C']}\ngamma:{grid['gamma']}"
                )
                print(
                    f"The best scoring parameter set now is:\n{pd.DataFrame(best_scoring)}\n-----\n"
                )
            elif output == "canvas" and n <= 2:
                canvas.append(
                    f"Round {str(n + 1)} gave the following optimal parameters:\n{pd.DataFrame(best_params)},\nfrom this grid:\nC    :{grid['C']}\ngamma:{grid['gamma']}"
                )
                canvas.append(
                    f"The best scoring parameter set now is:\n{pd.DataFrame(best_scoring)}\n-----\n"
                )
            if output == "print" and n > 2:
                print(
                    f"Round {str(n + 1)} gave the following optimal parameters:\n{pd.DataFrame(best_params)},\nfrom this grid and all previous iterations:\nC    :{grid['C']}\ngamma:{grid['gamma']}"
                )
                print(
                    f"The best scoring parameter set now is:\n{pd.DataFrame(best_scoring)}\n-----\n"
                )
            elif output == "canvas" and n > 2:
                canvas.append(
                    f"Round {str(n + 1)} gave the following optimal parameters:\n{pd.DataFrame(best_params)},\nfrom this grid and all previous iterations:\nC    :{grid['C']}\ngamma:{grid['gamma']}"
                )
                canvas.append(
                    f"The best scoring parameter set now is:\n{pd.DataFrame(best_scoring)}\n-----\n"
                )
            n += 1

        # Final comparison
        grid = {
            "C": [best_scoring[0], pd.DataFrame(best_params).mean(axis=1)[0]],
            "gamma": [best_scoring[1], pd.DataFrame(best_params).mean(axis=1)[1]],
        }
        best_params, accuracies_n = self._grid_search(train=train, parameters=grid)
        accuracies = pd.concat([accuracies, accuracies_n])
        accuracies = accuracies[~accuracies.index.duplicated(keep="first")]
        if output == "canvas":
            canvas[0] = SVMComp._plot_training(accuracies)
        best_scoring = SVMComp._best_estimator(accuracies, mode=autoselect)
        if output == "print":
            print(
                f"Comparison of mean and best scoring parameters gave the following optimal parameters:\n{pd.DataFrame(best_params)},\nfrom this grid:\nC    :{grid['C']}\ngamma:{grid['gamma']}"
            )
            print(
                f"The best scoring parameter set now is:\n{pd.DataFrame(best_scoring)}"
            )
        elif output == "canvas":
            canvas.append(
                f"Comparison of mean and best scoring parameters gave the following optimal parameters:\n{pd.DataFrame(best_params)},\nfrom this grid:\nC    :{grid['C']}\ngamma:{grid['gamma']}"
            )
            canvas.append(
                f"The best scoring parameter set now is:\n{pd.DataFrame(best_scoring)}"
            )

        # Store accuracies and set C and gamma
        self.accuracies = accuracies
        self.set_parameters(C=best_scoring[0], gamma=best_scoring[1])

    def _predict_svm(
        self,
        train=pd.DataFrame(
            columns=[DataFrameStrings.PROTEIN_IDS, DataFrameStrings.COMPARTMENT]
        ),
        C=-1,
        gamma=-1,
        random_state=None,
    ):
        """
        Apply C and gamma to the training set and get predictions for df using the predict_proba parameter of sklearn SVC.
        """
        if len(train) == 0:
            train = self._train_test_split(
                test_percent=self.test_split, random_state=self.random_state
            )[0]

        if C == -1:
            C = self.C
        if gamma == -1:
            gamma = self.gamma
        if random_state == None:
            random_state = self.random_state

        labels = train[DataFrameStrings.COMPARTMENT]

        results = pd.DataFrame()

        for exp in self.experiments:
            svc = svm.SVC(probability=True, C=C, gamma=gamma, random_state=random_state)
            data_all = (
                self.df.xs(exp, level=DataFrameStrings.EXPERIMENT, axis=0)
                .unstack(DataFrameStrings.MAP)
                .dropna()
            )
            data = data_all.reset_index(
                level=[DataFrameStrings.COMPARTMENT], drop=True
            ).loc[train[DataFrameStrings.PROTEIN_IDS], :]
            svc.fit(data.values, labels)
            proba = svc.predict_proba(data_all.values)
            proba = pd.DataFrame(
                proba,
                columns=pd.MultiIndex.from_tuples(
                    [(exp, el) for el in svc.classes_],
                    names=[DataFrameStrings.EXPERIMENT, DataFrameStrings.COMPARTMENT],
                ),
                index=data_all.index,
            )
            if len(results) == 0:
                results = proba.copy()
            else:
                results = results.join(proba, how="outer")
        self.probabilities = results.copy()

    def _score_predictions(
        self, probabilities=pd.DataFrame(), min_p=0.4, min_diff=0.15
    ):
        """
        To draw actual predictions from the probabilities, take the highest probability (Winner) and handle it according to the cutoffs.
        All probabilities > min_p are considered (Winners). If the difference between two such winners is at least min_diff only the larger one is the classification. If the difference is smaller, they are combined with an ' OR '. By the default values e.g. two predicitons with 0.45 probability will be '1 OR 2', two probabilities 0.57 and 0.41 will just be classified as '1', but winners '1,2'.
        """
        predictions = pd.DataFrame()
        if len(probabilities) == 0:
            probabilities = self.probabilities.copy()

        for exp in set(
            probabilities.columns.get_level_values(DataFrameStrings.EXPERIMENT)
        ):
            print(exp)
            p_e = probabilities[exp].dropna()
            p_max = p_e.max(axis=1)
            p_max.name = "Max probability"
            conf = pd.Series(
                [
                    "best guess"
                    if p <= min_p
                    else "very high confidence"
                    if p > 0.95
                    else "high confidence"
                    if p > 0.8
                    else "medium confidence"
                    if p > 0.65
                    else "low confidence"
                    if p > 0.4
                    else "best guess"
                    for p in p_max
                ],
                name="Confidence",
                index=p_max.index,
            )
            winner = p_e.idxmax(axis=1)
            winner.name = "Winner"
            winners = p_e.apply(
                lambda y: ";".join(y.index[np.where(y > min_p)]), axis=1
            )
            winners.name = "Winners"
            prediction = p_e.apply(
                lambda y: winners[y.name]
                if ";" not in winners[y.name]
                else winner[y.name]
                if y[winner[y.name]] - y.drop(winner[y.name]).max() > min_diff
                else winners[y.name].replace(";", " OR "),
                axis=1,
            )
            prediction.name = "Classification"
            if len(predictions) == 0:
                predictions = pd.DataFrame(
                    [p_max, winner, winners, prediction, conf],
                    index=pd.MultiIndex.from_arrays(
                        [
                            np.repeat(exp, 5),
                            [
                                "Max probability",
                                "Winner",
                                "Winners",
                                "Classification",
                                "Confidence",
                            ],
                        ],
                        names=[DataFrameStrings.EXPERIMENT, "Score"],
                    ),
                ).T
            else:
                predictions = predictions.join(
                    pd.DataFrame(
                        [p_max, winner, winners, prediction, conf],
                        index=pd.MultiIndex.from_arrays(
                            [
                                np.repeat(exp, 5),
                                [
                                    "Max probability",
                                    "Winner",
                                    "Winners",
                                    "Classification",
                                    "Confidence",
                                ],
                            ],
                            names=[DataFrameStrings.EXPERIMENT, "Score"],
                        ),
                    ).T,
                    how="outer",
                )
        self.predictions = predictions
        return predictions

    def _get_confusion(
        self,
        experiment: str,
        predictions=pd.DataFrame(),
        test=pd.DataFrame(
            columns=[DataFrameStrings.PROTEIN_IDS, DataFrameStrings.COMPARTMENT]
        ),
    ):
        """
        For whatever subset of df get the confusion matrix. This is absed on the wunner column, regardless of the probability cutoffs used for classification.
        """
        if len(predictions) == 0:
            predictions = self.predictions.copy()
        if len(test) == 0:
            test = self._train_test_split(test_percent=self.test_split)[1]

        test_predictions = predictions.loc[
            [
                el in test[DataFrameStrings.PROTEIN_IDS].values
                for el in predictions.index.get_level_values(
                    DataFrameStrings.PROTEIN_IDS
                )
            ],
            :,
        ]
        winner = test_predictions[(experiment, "Winner")].dropna()
        confusion = (
            winner.groupby(DataFrameStrings.COMPARTMENT)
            .value_counts()
            .unstack(DataFrameStrings.COMPARTMENT)
            .T
        )
        confusion_formatted = pd.DataFrame(
            index=list(set(test[DataFrameStrings.COMPARTMENT].values)),
            columns=list(set(test[DataFrameStrings.COMPARTMENT].values)),
        )
        for i in confusion_formatted.columns:
            for j in confusion_formatted.index:
                try:
                    confusion_formatted.loc[j, i] = confusion.loc[j, i]
                except:
                    pass

        return confusion_formatted


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

    # df_SVM = self.df_SVM.copy()
    # if hasattr(df_SVM, "keys") == True:
    try:
        df_SVM = pd.read_json(df_SVM["Misclassification Matrix"])
        df_SVM = df_SVM.set_index("T: True group")[::-1]

    except:
        pass

    y_axis_label = df_SVM.index
    x_axis_label = df_SVM.columns
    data_svm = df_SVM.values
    fig_SVMheatmap = go.Figure()

    fig_SVMheatmap.add_trace(
        go.Heatmap(
            z=data_svm,
            x=x_axis_label,
            y=y_axis_label if y_axis_label[0] != 0 else x_axis_label,
            colorscale=[[0.0, "green"], [0.01, "white"], [1.0, "red"]],
        )
    )

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
            if (
                data_type == "0/1 normalized data"
                and exp_name == list(json_dict.keys())[0]
            ):
                df_01_combined = pd.read_json(json_dict[exp_name][data_type])
                df_01_combined = df_01_combined.set_index(
                    [
                        DataFrameStrings.GENE_NAMES,
                        DataFrameStrings.PROTEIN_IDS,
                        DataFrameStrings.COMPARTMENT,
                    ]
                ).copy()
                df_01_combined.drop(
                    [
                        col
                        for col in df_01_combined.columns
                        if not col.startswith(DataFrameStrings.NORMALIZED_PROFILE)
                    ]
                )
                df_01_combined.columns = pd.MultiIndex.from_tuples(
                    [el.split("?") for el in df_01_combined.columns],
                    names=[
                        DataFrameStrings.SET,
                        DataFrameStrings.MAP,
                        DataFrameStrings.FRACTION,
                    ],
                )
                df_01_combined.rename(
                    columns={DataFrameStrings.NORMALIZED_PROFILE: exp_name},
                    inplace=True,
                )

            elif (
                data_type == "0/1 normalized data"
                and exp_name != list(json_dict.keys())[0]
            ):
                df_01_toadd = pd.read_json(json_dict[exp_name][data_type])
                df_01_toadd = df_01_toadd.set_index(
                    [
                        DataFrameStrings.GENE_NAMES,
                        DataFrameStrings.PROTEIN_IDS,
                        DataFrameStrings.COMPARTMENT,
                    ]
                ).copy()
                df_01_toadd.drop(
                    [
                        col
                        for col in df_01_toadd.columns
                        if not col.startswith(DataFrameStrings.NORMALIZED_PROFILE)
                    ]
                )
                df_01_toadd.columns = pd.MultiIndex.from_tuples(
                    [el.split("?") for el in df_01_toadd.columns],
                    names=[
                        DataFrameStrings.SET,
                        DataFrameStrings.MAP,
                        DataFrameStrings.FRACTION,
                    ],
                )
                df_01_toadd.rename(
                    columns={DataFrameStrings.NORMALIZED_PROFILE: exp_name},
                    inplace=True,
                )
                df_01_combined = pd.concat([df_01_combined, df_01_toadd], axis=1)

    df_01_combined.columns.names = [
        DataFrameStrings.EXPERIMENT,
        DataFrameStrings.MAP,
        DataFrameStrings.FRACTION,
    ]
    df = df_01_combined.stack(
        [DataFrameStrings.EXPERIMENT, DataFrameStrings.MAP]
    ).dropna(axis=0)
    df = df.div(df.sum(axis=1), axis=0)
    index_ExpMap = (
        df.index.get_level_values(DataFrameStrings.EXPERIMENT)
        + "_"
        + df.index.get_level_values(DataFrameStrings.MAP)
    )
    index_ExpMap.name = DataFrameStrings.EXP_MAP
    df.set_index(index_ExpMap, append=True, inplace=True)

    df.index = df.index.droplevel([DataFrameStrings.MAP, DataFrameStrings.EXPERIMENT])
    df = df.stack(DataFrameStrings.FRACTION).unstack(
        [DataFrameStrings.EXP_MAP, DataFrameStrings.FRACTION]
    )
    df.columns = ["_".join(col) for col in df.columns.values]

    return df


def parse_reannotation_source(mode, source):
    reannotation_source = {}
    reannotation_source["source"] = mode
    reannotate = bool(mode)
    if mode == "fasta_headers":
        if "\n" in source:
            idmapping = [el for el in source.split("\n")]
        else:
            idmapping = [
                el.decode("UTF-8").strip()
                for el in pkg_resources.resource_stream(
                    "domaps", f"annotations/idmapping/{source}.txt"
                ).readlines()
            ]
        reannotation_source["idmapping"] = idmapping
    elif mode == "tsv":
        if "\n" in source:
            idmapping = pd.read_csv(
                StringIO(source), sep="\t", usecols=["Entry", "Gene names  (primary )"]
            )
        else:
            idmapping = pd.read_csv(
                pkg_resources.resource_stream(
                    "domaps", f"annotations/idmapping/{source}.tab"
                ),
                sep="\t",
                usecols=["Entry", "Gene names  (primary )"],
            )
        reannotation_source["idmapping"] = idmapping
    return reannotation_source


def generate_usecols_regex(
    sets: dict,
    original_protein_ids: str,
    genes: str,
    samples: str = None,
    column_filters: dict = None,
    annotation_columns: list = None,
):
    columns = []

    # sets
    for v in sets.values():
        columns.append(v)

    # mandatory columns
    columns.append(original_protein_ids)
    columns.append(genes)
    if samples is not None:
        columns.append(samples)

    # column filters
    if column_filters is not None:
        for k in column_filters.keys():
            columns.append(k)

    # index columns
    if annotation_columns is not None:
        columns += annotation_columns

    columns = set(columns)
    return "^" + "$|^".join(columns) + "$"


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

        url = "https://www.uniprot.org/uploadlists/"

        params = {
            "from": "ACC+ID",
            "to": "GENENAME",
            "format": "tab",
            "query": "\t".join(protein_ids),
        }

        data = urllib.parse.urlencode(params)
        data = data.encode("utf-8")
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
            response = f.read()
        for l in response.decode("utf-8").split("\n")[1:-1]:
            l = l.split("\t")
            if l[0] in protein_gene.keys():
                protein_gene[l[0]] = protein_gene[l[0]] + "/" + l[1]
            else:
                protein_gene[l[0]] = l[1]

    elif source == "fasta_headers":
        for header in idmapping:
            p_g = re.findall("\|([^-\|]+)\|[^;]+GN=([^ ]+) ", header)
            for el in p_g:
                protein_gene[el[0]] = el[1]

    elif source == "tsv":
        for p, g in zip(idmapping["Entry"], idmapping["Gene names  (primary )"]):
            if type(g) == str:
                protein_gene[p] = g.replace("; ", "/")

    else:
        raise ValueError(f"Unkown source for gene reannotation: {source}")

    genes = list()
    for ids in split_ids:
        ids = [ids[el] for el in sorted([ids.index(i) for i in set(ids)])]
        genes.append(
            ";".join(
                [p if p not in protein_gene.keys() else protein_gene[p] for p in ids]
            )
        )

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

    src = np.array(
        [
            exp,
            [el.split(sep) for el in pgs],
            [el.count(sep) for el in pgs],
            [el.split(sep) for el in genes],
            [el.count(sep) for el in genes],
            [i for i in range(len(genes))],
        ],
        dtype=object,
    )

    src_single_gene = src[:, src[4] == 0]
    src_single_gene = src_single_gene[:, src_single_gene[2].argsort()[::-1]]
    src_multi_gene = src[:, src[4] != 0]
    src_multi_gene = src_multi_gene[
        :, np.lexsort((src_multi_gene[2], src_multi_gene[4]))[::-1]
    ]

    out = np.empty([5, 0])

    # Process multi gene entries
    while src_multi_gene.shape[1] > 0:
        group_max = src_multi_gene[3, 0]
        pg_max = src_multi_gene[1, 0]

        # find members from groups with multiple gene names
        members_multi = [0]
        for i, pg in enumerate(src_multi_gene[1]):
            if (
                all(el in pg_max for el in pg)
                and src_multi_gene[0, i] not in src_multi_gene[0, members_multi]
            ):
                members_multi.append(i)
        group = src_multi_gene[:, members_multi]
        src_multi_gene = np.delete(src_multi_gene, members_multi, axis=1)

        # find members from groups with single gene names
        members_single = []
        for i, pg in enumerate(src_single_gene[1]):
            if (
                all(el in pg_max for el in pg)
                and src_single_gene[0, i] not in group[0]
                and src_single_gene[0, i] not in src_single_gene[0, members_single]
            ):
                members_single.append(i)
        group = np.append(group, src_single_gene[:, members_single], axis=1)
        src_single_gene = np.delete(src_single_gene, members_single, axis=1)

        # sort group by number of genes, then by number of protein groups and lastly
        group[1] = [sep.join(el) for el in group[1]]
        group[3] = [sep.join(el) for el in group[3]]
        pg_counts = {k: v for k, v in zip(*np.unique(group[1], return_counts=True))}
        group = np.append(group, [[pg_counts[el] for el in group[1]]], axis=0)
        group = group[:, np.lexsort((group[6], group[2], group[4]))]
        # print(group_max, group)

        group_out = np.array(
            [
                group[0],
                np.repeat(group[1, -1], group.shape[1]),
                np.repeat(group[3, -1], group.shape[1]),
                np.repeat(
                    "multiple genes"
                    if len(np.unique(group[3])) == 1
                    else "gene level conflict",
                    group.shape[1],
                ),
                group[5],
            ]
        )
        out = np.append(out, group_out, axis=1)

    # Process leftover single gene entries
    single_out = np.array(
        [
            src_single_gene[0],
            [split_ids(sep.join(el)) for el in src_single_gene[1]],
            [el[0] for el in src_single_gene[3]],
            ["primary id" for el in src_single_gene[2]],
            src_single_gene[5],
        ]
    )
    out = np.append(out, single_out, axis=1)
    out = out[0:4, out[-1].argsort()].T
    return pd.DataFrame(
        out,
        columns=[
            DataFrameStrings.EXPERIMENT,
            DataFrameStrings.PROTEIN_IDS,
            DataFrameStrings.GENE_NAMES,
            "merge type",
        ],
    )


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

    assert all(
        [
            el in df.index.names
            for el in [
                DataFrameStrings.EXPERIMENT,
                DataFrameStrings.GENE_NAMES,
                DataFrameStrings.COMPARTMENT,
                DataFrameStrings.ORIGINAL_PROTEIN_IDS,
            ]
        ]
    )
    assert all(
        [
            el in df.columns.names
            for el in [DataFrameStrings.MAP, DataFrameStrings.FRACTION]
        ]
    )

    df_01 = df.copy()
    n_exp = len(set(df_01.index.get_level_values(DataFrameStrings.EXPERIMENT)))

    # 2. Retrieve data with simple full coverage alignment so it doesn't get torn during the more complicated matching
    df_01.set_index(
        pd.Index(
            [
                split_ids(el)
                for el in df_01.index.get_level_values(
                    DataFrameStrings.ORIGINAL_PROTEIN_IDS
                )
            ],
            name=DataFrameStrings.PROTEIN_IDS,
        ),
        append=True,
        inplace=True,
    )
    full_cov = df_01.groupby(DataFrameStrings.PROTEIN_IDS).filter(
        lambda x: x.shape[0] >= n_exp
    )

    # 3. Relabel more complicated cases
    compl = df_01.drop(full_cov.index)
    new_labels = relabel_groups(
        compl.index.get_level_values(DataFrameStrings.EXPERIMENT),
        compl.index.get_level_values(DataFrameStrings.ORIGINAL_PROTEIN_IDS),
        [str(el) for el in compl.index.get_level_values(DataFrameStrings.GENE_NAMES)],
    )
    compl_relabeled = compl.copy()
    for i in new_labels.columns:
        if i in compl_relabeled.index.names:
            compl_relabeled.index = compl_relabeled.index.droplevel(i)
        compl_relabeled.set_index(
            pd.Index(new_labels[i], name=i), append=True, inplace=True
        )
    full_cov.set_index(
        pd.Index(np.repeat("primary id", len(full_cov)), name="merge type"),
        append=True,
        inplace=True,
    )
    compl_relabeled.index = compl_relabeled.index.reorder_levels(full_cov.index.names)

    # 4. Get index mapping
    index_mapping = pd.concat(
        [
            full_cov.index.to_frame(index=False),
            compl_relabeled.index.to_frame(index=False),
        ],
        axis=0,
    )
    index_mapping.set_index(
        [
            DataFrameStrings.EXPERIMENT,
            DataFrameStrings.GENE_NAMES,
            DataFrameStrings.PROTEIN_IDS,
            "merge type",
        ],
        inplace=True,
    )
    index_mapping = index_mapping.unstack(DataFrameStrings.EXPERIMENT)

    # 5. Rejoin datasets and relabel Compartments
    df_01_aligned = pd.concat([full_cov, compl_relabeled], axis=0)
    df_01_aligned.index = df_01_aligned.index.droplevel(
        [DataFrameStrings.ORIGINAL_PROTEIN_IDS, DataFrameStrings.COMPARTMENT]
    )
    df_01_aligned = df_01_aligned.unstack(DataFrameStrings.EXPERIMENT)  # rejoin index
    c = []
    for el in df_01_aligned.index.get_level_values(DataFrameStrings.PROTEIN_IDS):
        cs = [
            el
            for el in index_mapping.xs(el, level=DataFrameStrings.PROTEIN_IDS, axis=0)
            .Compartment.iloc[0, :]
            .values
            if type(el) == str
        ]
        if len(set(cs)) == 1:
            c.append(cs[0])
        else:
            c.append("undefined")
    df_01_aligned.set_index(
        pd.Index(c, name=DataFrameStrings.COMPARTMENT), append=True, inplace=True
    )
    df_01_aligned = df_01_aligned.stack(
        [DataFrameStrings.EXPERIMENT, DataFrameStrings.MAP]
    )

    return df_01_aligned, index_mapping


def relabel_fractions(df: pd.DataFrame, fraction_mapping: dict = dict()):
    """
    Relabel (and remove) fractions according to a dictionary.

    Old fraction labels are keys, new fraction labels are values. If the value is None the fraction is dropped.

    """

    df_fractions = df.copy()

    ## Rename and potentially delete fractions
    # empty dictionary: leave as is
    if len(fraction_mapping) == 0:
        pass
    else:
        for k, v in fraction_mapping.items():
            if v == None:
                df_fractions.drop(
                    k, axis=1, inplace=True, level=DataFrameStrings.FRACTION
                )
            elif k != v:
                df_fractions.rename(columns={k: v}, inplace=True)
    return df_fractions


def format_data_long(
    df: pd.DataFrame,
    original_protein_ids: str,
    genes: str,
    samples: str,
    name_pattern: str,
    sets: dict,
    index_cols: list = [],
    fraction_mapping: dict = dict(),
):
    # fmt: off
    """
    This formats proteomic input data in long format to be compatible with the SpatialDataset class.

    >>> pd.set_option('display.max_columns', 500)
    >>> pd.set_option('display.width', 1000)
    >>> format_data_long(pd.DataFrame([["foo", "lorem", "rep1_F1", 100, 4],
    ...                                ["bar;bar-1", "ipsum", "rep1_F1", 200, 1],
    ...                                ["foo", "lorem", "rep2_F1", 50, 0],
    ...                                ["bar;bar-1", "ipsum", "rep2_F1", 30.1, 20],
    ...                                ["foo", "lorem", "rep1_F2", 0, np.nan],
    ...                                ["bar;bar-1", "ipsum", "rep2_F2", np.nan, np.nan],
    ...                                ["foo", "lorem", "rep2_F2", 2e10, 6]],
    ...                               columns=["PG.ProteinGroups", "PG.Genes", "R.Condition", "PG.Quantity", "PG.RunEvidenceCount"]
    ...                               ),
    ...                  original_protein_ids="PG.ProteinGroups", genes="PG.Genes",
    ...                  samples="R.Condition", name_pattern="(?P<rep>.*)_(?P<frac>.*)",
    ...                  sets={"LFQ intensity": "PG.Quantity", "MS/MS count": "PG.RunEvidenceCount"})
    Set                                         LFQ intensity                          MS/MS count               
    Map                                                  rep1       rep2                      rep1      rep2     
    Fraction                                               F1   F2    F1            F2          F1  F2    F1   F2
    Original Protein IDs Gene names Protein IDs                                                                  
    bar;bar-1            ipsum      bar                 200.0  NaN  30.1           NaN         1.0 NaN  20.0  NaN
    foo                  lorem      foo                 100.0  0.0  50.0  2.000000e+10         4.0 NaN   0.0  6.0


    Use fractions parameter to relabel fractions
    >>> format_data_long(pd.DataFrame([["foo", "lorem", "rep1_F1", 100],
    ...                                ["bar;bar-1", "ipsum", "rep1_F1", 200],
    ...                                ["foo", "lorem", "rep1_F3", 50],
    ...                                ["bar;bar-1", "ipsum", "rep1_F3", 30.1],
    ...                                ["foo", "lorem", "rep1_F2", 0],
    ...                                ["bar;bar-1", "ipsum", "rep1_F2", np.nan]],
    ...                               columns=["PG.ProteinGroups", "PG.Genes", "R.Condition", "PG.Quantity"]
    ...                               ),
    ...                  original_protein_ids="PG.ProteinGroups", genes="PG.Genes",
    ...                  samples="R.Condition", name_pattern="(?P<rep>.*)_(?P<frac>.*)",
    ...                  sets={"LFQ intensity": "PG.Quantity"},
    ...                  fraction_mapping={"F1": None, "F2": "F2", "F3": "F1"}) # Drop F1, relabel F3 as F1, leave F2 as is
    Set                                         LFQ intensity      
    Map                                                  rep1      
    Fraction                                               F2    F1
    Original Protein IDs Gene names Protein IDs                    
    bar;bar-1            ipsum      bar                   NaN  30.1
    foo                  lorem      foo                   0.0  50.0


    >>> format_data_long(pd.DataFrame([["foo", "rep1_F1", 100],
    ...                                ["bar;bar-1", "rep1_F1", 200],
    ...                                ["foo", "rep1_F3", 50],
    ...                                ["bar;bar-1", "rep1_F3", 30.1],
    ...                                ["foo", "rep1_F2", 0],
    ...                                ["bar;bar-1", "rep1_F2", np.nan]],
    ...                               columns=["PG.ProteinGroups", "R.Condition", "PG.Quantity"]
    ...                               ),
    ...                  original_protein_ids="PG.ProteinGroups", genes="PG.ProteinGroups",
    ...                  samples="R.Condition", name_pattern="(?P<rep>.*)_(?P<frac>.*)",
    ...                  sets={"LFQ intensity": "PG.Quantity"},
    ...                  fraction_mapping={"F1": None, "F2": "F2", "F3": "F1"}) # Drop F1, relabel F3 as F1, leave F2 as is
    Set                                         LFQ intensity      
    Map                                                  rep1      
    Fraction                                               F2    F1
    Original Protein IDs Gene names Protein IDs                    
    bar;bar-1            bar;bar-1  bar                   NaN  30.1
    foo                  foo        foo                   0.0  50.0
    """
    # fmt: on

    ## Rename columns
    df_index = df.rename(
        columns={
            original_protein_ids: DataFrameStrings.ORIGINAL_PROTEIN_IDS,
            genes: DataFrameStrings.GENE_NAMES,
            samples: "Samples",
            **{v: k for k, v in sets.items()},
        },
        errors="raise",
    )
    if genes == original_protein_ids:
        df_index.insert(
            0,
            DataFrameStrings.ORIGINAL_PROTEIN_IDS,
            df_index[DataFrameStrings.GENE_NAMES],
        )

    ## Set index columns and column name "Set"
    df_index.dropna(subset=[DataFrameStrings.ORIGINAL_PROTEIN_IDS], inplace=True)
    df_index.set_index(
        [DataFrameStrings.ORIGINAL_PROTEIN_IDS, DataFrameStrings.GENE_NAMES, "Samples"]
        + index_cols,
        inplace=True,
    )
    df_index.columns.names = [DataFrameStrings.SET]

    ## Catch any additional columns, which are not accounted for
    if any([col not in sets.keys() for col in df_index.columns]):
        raise ValueError(
            f"Additional unknown column(s) {', '.join([col for col in df_index.columns if col not in sets.keys()])} present. Please specify what these are or remove them from the uploaded file."
        )

    ## Convert to numeric type
    df_index = df_index.apply(pd.to_numeric, errors="coerce")

    ## Unstack samples and apply regex
    df_index = df_index.unstack("Samples")
    df_index.columns = pd.MultiIndex.from_arrays(
        [
            df_index.columns.get_level_values(DataFrameStrings.SET),
            [
                re.match(name_pattern, i).group("rep")
                if "<cond>" not in name_pattern
                else "_".join(re.match(name_pattern, i).group("cond", "rep"))
                for i in df_index.columns.get_level_values("Samples")
            ],
            [
                re.match(name_pattern, i).group("frac")
                for i in df_index.columns.get_level_values("Samples")
            ],
        ],
        names=[DataFrameStrings.SET, DataFrameStrings.MAP, DataFrameStrings.FRACTION],
    )

    ## Remnant from deprecated code in case this is ever needed: Fractionated data that is not actually aggregated by the search engine
    # In case fractionated data was used this needs to be catched and aggregated
    # try:
    #    df_index = df_index.unstack(["Map", "Fraction"])
    # except ValueError:
    #    df_index = df_index.groupby(by=df_index.index.names).agg(np.nansum, axis=0)
    #    df_index = df_index.unstack(["Map", "Fraction"])

    ## Shorten Protein IDs
    if DataFrameStrings.PROTEIN_IDS in df_index.index.names:
        df_index.index.rename(
            "Protein IDs loaded", level=DataFrameStrings.PROTEIN_IDS, inplace=True
        )
    df_index.set_index(
        pd.Index(
            [
                split_ids(el)
                for el in df_index.index.get_level_values(
                    DataFrameStrings.ORIGINAL_PROTEIN_IDS
                )
            ],
            name=DataFrameStrings.PROTEIN_IDS,
        ),
        append=True,
        inplace=True,
    )

    ## Relabel and remove fractions
    if len(fraction_mapping) == 0:
        pass
    else:
        df_index = relabel_fractions(df_index, fraction_mapping)

    return df_index


def format_data_pivot(
    df: pd.DataFrame,
    original_protein_ids: str,
    genes: str,
    name_pattern: str,
    sets: dict,
    index_cols: list = [],
    fraction_mapping: dict = dict(),
):
    # fmt: off
    """
    >>> pd.set_option('display.max_columns', 500)
    >>> pd.set_option('display.width', 1000)
    >>> format_data_pivot(pd.DataFrame([["foo", "lorem", 20,np.nan,20,10, 2,0,4,3, None],
    ...                                 ["bar;bar-1", "ipsum", 0,40,20,40, np.nan,10,5,3, None]],
    ...                                columns=["Majority Protein IDs", "Gene Names",
    ...                                         "LFQ intensity 1_F1","LFQ intensity 1_F2","LFQ intensity 2_F1","LFQ intensity 2_F2",
    ...                                         "MS/MS counts 1_F1","MS/MS counts 1_F2","MS/MS counts 2_F1","MS/MS counts 2_F2", "Reverse"]),
    ...                       original_protein_ids="Majority Protein IDs", genes="Gene Names",
    ...                       name_pattern=".* (?P<rep>.*)_(?P<frac>.*)",
    ...                       sets={"LFQ intensity": "LFQ intensity .*", "MS/MS count": "MS/MS counts .*"},
    ...                       index_cols=["Reverse"])
    Set                                                 LFQ intensity               MS/MS count            
    Map                                                             1         2               1        2   
    Fraction                                                       F1    F2  F1  F2          F1    F2 F1 F2
    Original Protein IDs Gene names Reverse Protein IDs                                                    
    foo                  lorem      NaN     foo                  20.0   NaN  20  10         2.0   NaN  4  3
    bar;bar-1            ipsum      NaN     bar                   NaN  40.0  20  40         NaN  10.0  5  3
    """
    # fmt: on

    ## Rename columns
    df_index = df.rename(
        columns={
            original_protein_ids: DataFrameStrings.ORIGINAL_PROTEIN_IDS,
            genes: DataFrameStrings.GENE_NAMES,
        },
        errors="raise",
    )
    if genes == original_protein_ids:
        df_index.insert(
            0,
            DataFrameStrings.ORIGINAL_PROTEIN_IDS,
            df_index[DataFrameStrings.GENE_NAMES],
        )

    ## Set index columns and column name "Set"
    df_index.dropna(
        subset=[DataFrameStrings.ORIGINAL_PROTEIN_IDS], inplace=True, axis=0
    )
    df_index.set_index(
        [DataFrameStrings.ORIGINAL_PROTEIN_IDS, DataFrameStrings.GENE_NAMES]
        + index_cols,
        inplace=True,
    )

    ## Catch any additional columns, which are not accounted for
    if any([not re.match("|".join(sets.values()), col) for col in df_index.columns]):
        raise ValueError(
            f"Additional unknown column(s) {', '.join([col for col in df_index.columns if not re.match('|'.join(sets.values()), col)])} present. Please specify what these are or remove them from the uploaded file."
        )

    ## Convert to numeric type
    df_index = df_index.apply(pd.to_numeric, errors="coerce")

    # multindex will be generated, by extracting the information about the Map, Fraction and Type from each individual column name
    multiindex = pd.MultiIndex.from_arrays(
        arrays=[
            [
                [k for k, v in sets.items() if re.match(v, col)][0]
                for col in df_index.columns
            ],
            [re.match(name_pattern, col).group("rep") for col in df_index.columns]
            if "<cond>" not in name_pattern
            else [
                "_".join(re.match(name_pattern, col).group("cond", "rep"))
                for col in df_index.columns
            ],
            [re.match(name_pattern, col).group("frac") for col in df_index.columns],
        ],
        names=[DataFrameStrings.SET, DataFrameStrings.MAP, DataFrameStrings.FRACTION],
    )

    df_index.columns = multiindex
    #### TODO change this in the old and in the new code to only affect the main column set. This mostly affects SILAC data, where ratio variabilities of 0 would be voided, thereby affecting the normalization
    df_index.replace({0: np.nan}, inplace=True)

    ## Shorten Protein IDs
    if DataFrameStrings.PROTEIN_IDS in df_index.index.names:
        df_index.index.rename(
            "Protein IDs loaded", level=DataFrameStrings.PROTEIN_IDS, inplace=True
        )
    df_index.set_index(
        pd.Index(
            [
                split_ids(el)
                for el in df_index.index.get_level_values(
                    DataFrameStrings.ORIGINAL_PROTEIN_IDS
                )
            ],
            name=DataFrameStrings.PROTEIN_IDS,
        ),
        append=True,
        inplace=True,
    )

    ## Relabel and remove fractions
    if len(fraction_mapping) == 0:
        pass
    else:
        df_index = relabel_fractions(df_index, fraction_mapping)

    return df_index


def filter_SILAC_countvar(
    df,
    RatioCount: int = 2,
    RatioVariability: float = 30,
    sets: dict = dict(
        count=DataFrameStrings.RATIO_COUNT,
        variability=DataFrameStrings.RATIO_VARIABILITY,
    ),
):
    # fmt: off
    """
    The multiindex dataframe is subjected to stringency filtering. Proteins were retained with 3 or more quantifications in each
    subfraction (=count). Furthermore, proteins with only 2 quantification events in one or more subfraction were retained, if their ratio variability for
    ratios obtained with 2 quantification events was below 30% (=var). SILAC ratios were linearly normalized by division through the fraction median.
    Subsequently normalization to SILAC loading was performed.Data is annotated based on specified marker set e.g. eLife.

    Args:
        df: pd.DataFrame with column names as defined by sets on column level "Set" and a column or index level "Fraction".
        RatioCount: int, 2; Minimum ratio count
        RatioVariability: float, 30; If ratio count is exactly its set value, additionaly check ratio variability <= this value.
        sets: dictionary, dict(ratio="Ratio", count="Ratio count", variability="Ratio variability")

    Returns:
        df_filtered: pd.DataFrame, same as original data frame, but only with data for complete profiles remaining after filtering

    Doctest for expected input
    >>> filter_SILAC_countvar(pd.DataFrame([[1,0, 3,2, 3,15], # no bad values
    ...                                     [2,0, 3,2, 5,40], # second value 2 counts, but high variability
    ...                                     [3,0, 1,2, 5,5]], # first value only one count
    ...                                    index=pd.Index(["a","b","c"],name="id"),
    ...                                    columns=pd.MultiIndex.from_arrays(
    ...                                        [["Ratio","Ratio","Ratio count","Ratio count","Ratio variability","Ratio variability"],
    ...                                         ["F1","F2","F1","F2","F1","F2"]],
    ...                                        names=["Set", "Fraction"])))
    Set      Ratio      Ratio count      Ratio variability      
    Fraction    F1   F2          F1   F2                F1    F2
    id                                                          
    a          1.0  0.0         3.0  2.0               3.0  15.0
    b          2.0  NaN         3.0  NaN               5.0   NaN
    c          NaN  0.0         NaN  2.0               NaN   5.0
    """
    # fmt: on

    ## Fetch column levels and assert presence of the required sets and fraction annotations
    column_levels = df.columns.names
    # All sets present
    if DataFrameStrings.SET not in column_levels:
        raise KeyError("Column level 'Set' is missing.")
    else:
        if any(
            [
                v not in df.columns.get_level_values(DataFrameStrings.SET)
                for v in sets.values()
            ]
        ):
            raise KeyError(
                f"One of {', '.join(sets.values())} is not in the sets required for filtering SILAC data by counts and variability."
            )
    # Fraction annotation present
    if DataFrameStrings.FRACTION not in column_levels:
        if DataFrameStrings.FRACTION not in df.index.names:
            raise KeyError("Column or index level 'Fraction' is missing.")

    ## Convert to long dataframe
    if len(column_levels) > 1:
        df_stack = df.stack([el for el in column_levels if el != DataFrameStrings.SET])
    else:
        df_stack = df.copy()

    ## FILTER DATA
    # count and variability
    df_stack_filtered = df_stack.loc[
        [
            count > RatioCount or (count == RatioCount and var < RatioVariability)
            for var, count in zip(
                df_stack[sets["variability"]], df_stack[sets["count"]]
            )
        ]
    ]

    # Return to input format
    df_filtered = (
        df_stack_filtered.unstack(
            [el for el in column_levels if el != DataFrameStrings.SET]
        )
        .reorder_levels(column_levels, axis=1)
        .sort_index(axis=1, key=natsort_index_keys)
    )

    return df_filtered


def filter_msms_count(
    df,
    average_MSMS_counts: int = 2,
    sets: dict = dict(msms=DataFrameStrings.MSMS_COUNT),
):
    # fmt: off
    """
    Filter proteomic profiling data for average ms/ms counts per profile.

    If summed MS/MS counts are fewer than (number of fraction) * (required average) profiles and all associated data are removed.

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

    Doctest for expected input:
    >>> filter_msms_count(pd.DataFrame([[1,2,1, 5,2,3], # msms sum = 10 > 6
    ...                                 [2,0,1, 5,np.nan,1], # msms sum = 6 = 6
    ...                                 [1,1,1, 1,1,1]], # msms sum = 3 < 6 -> remove
    ...                                index=pd.Index(["a","b","c"],name="id"),
    ...                                columns=pd.MultiIndex.from_arrays([["LFQ","LFQ","LFQ","MS/MS count","MS/MS count","MS/MS count"],
    ...                                                                  ["F1","F2","F3","F1","F2","F3"]], names=["Set", "Fraction"])))
    Set      LFQ       MS/MS count          
    Fraction  F1 F2 F3          F1   F2   F3
    id                                      
    a          1  2  1         5.0  2.0  3.0
    b          2  0  1         5.0  NaN  1.0
    """
    # fmt: on

    ## Fetch column levels and assert presence of the required sets and fraction annotations
    column_levels = df.columns.names
    # All sets present
    if DataFrameStrings.SET not in column_levels:
        raise KeyError("Column level 'Set' is missing.")
    else:
        if any(
            [
                v not in df.columns.get_level_values(DataFrameStrings.SET)
                for v in sets.values()
            ]
        ):
            raise KeyError(
                f"One of {', '.join(sets.values())} is not in the sets required for filtering abundances by MS/MS counts."
            )
    # Fraction annotation present
    if DataFrameStrings.FRACTION not in column_levels:
        if DataFrameStrings.FRACTION not in df.index.names:
            raise KeyError("Column or index level 'Fraction' is missing.")

    ## Convert to profile dataframe
    if len(column_levels) > 1:
        df_profiles = df.stack(
            [el for el in column_levels if el != DataFrameStrings.SET]
        ).unstack(DataFrameStrings.FRACTION)
    else:
        df_profiles = df.unstack(DataFrameStrings.FRACTION)

    # calculate minimal number of MS/MS counts per profile from number of fractions and required average
    minms = (
        len(set(df_profiles.columns.get_level_values(DataFrameStrings.FRACTION)))
        * average_MSMS_counts
    )

    # FILTER DATA
    if minms > 0:
        df_filtered = df_profiles.loc[
            df_profiles[sets["msms"]].apply(np.nansum, axis=1) >= minms
        ]
    else:
        df_filtered = df_profiles.copy()

    # Return to input format
    df_filtered = (
        df_filtered.stack(DataFrameStrings.FRACTION)
        .unstack([el for el in column_levels if el != DataFrameStrings.SET])
        .reorder_levels(column_levels, axis=1)
        .sort_index(axis=1, key=natsort_index_keys)
    )

    return df_filtered


def filter_consecutive(
    df,
    consecutive: int = 4,
    fraction_order: list = [],
    sets: dict = {DataFrameStrings.ABUNDANCE: DataFrameStrings.LFQ_INTENSITY},
):
    # fmt: off
    """
    >>> filter_consecutive(pd.DataFrame([[1,2,1,5,2,3], # complete profile
    ...                                  [0,0,1,5,3,1], # partial profile
    ...                                  [1,1,2,np.nan,1,1], # fragmented profile with nan
    ...                                  [1,1,0,0,1,1]], # fragmented profile
    ...                                 index=pd.Index(["a","b","c","d"],name="id"),
    ...                                 columns=pd.MultiIndex.from_arrays([["LFQ","LFQ","LFQ","LFQ","LFQ","LFQ"],
    ...                                                                   ["F1","F2","F3","F4","F5","F6"]], names=["Set", "Fraction"])), sets=dict(Abundance="LFQ"))
    Set       LFQ                         
    Fraction   F1   F2   F3   F4   F5   F6
    id                                    
    a         1.0  2.0  1.0  5.0  2.0  3.0
    b         NaN  NaN  1.0  5.0  3.0  1.0
    """
    # fmt: on

    ## Fetch column levels and assert presence of the required sets and fraction annotations
    column_levels = df.columns.names
    # All sets present
    if DataFrameStrings.SET not in column_levels:
        raise KeyError("Column level 'Set' is missing.")
    else:
        if any(
            [
                v not in df.columns.get_level_values(DataFrameStrings.SET)
                for v in sets.values()
            ]
        ):
            raise KeyError(
                f"One of {', '.join(sets.values())} is not in the sets required for filtering for consecutive values."
            )
    # Fraction annotation present
    if DataFrameStrings.FRACTION not in column_levels:
        if DataFrameStrings.FRACTION not in df.index.names:
            raise KeyError("Column or index level 'Fraction' is missing.")

    ## Convert to profile dataframe
    if len(column_levels) > 1:
        df_profiles = df.stack(
            [el for el in column_levels if el != DataFrameStrings.SET]
        ).unstack(DataFrameStrings.FRACTION)
    else:
        df_profiles = df.unstack(DataFrameStrings.FRACTION)

    ## Make sure fractions are in the correct order:
    if len(fraction_order) == 0:
        df_profiles.sort_index(
            level=DataFrameStrings.FRACTION,
            axis=1,
            key=natsort_index_keys,
            inplace=True,
        )
    else:
        df_profiles = df_profiles.reindex(
            fraction_order, level=DataFrameStrings.FRACTION, axis=1
        )

    ## Replace 0 with np.nan
    df_profiles = df_profiles.replace({0: np.nan})

    ## FILTER DATA
    df_filtered = df_profiles.loc[
        df_profiles[sets[DataFrameStrings.ABUNDANCE]]
        .apply(lambda x: np.isfinite(x), axis=0)
        .apply(
            lambda x: sum(x) >= consecutive
            and any(x.rolling(window=consecutive).sum() >= consecutive),
            axis=1,
        )
    ]

    # Return to input format
    df_filtered = (
        df_filtered.stack(DataFrameStrings.FRACTION)
        .unstack([el for el in column_levels if el != DataFrameStrings.SET])
        .reorder_levels(column_levels, axis=1)
        .sort_index(axis=1, key=natsort_index_keys)
    )

    return df_filtered


def filter_singlecolumn_keep(df, column: str, operator: str = "!=", value="'+'"):
    # fmt: off
    """
    >>> filter_singlecolumn_keep(pd.DataFrame([[1,2],[3,4],[5,6]],
    ...                                   columns=pd.MultiIndex.from_tuples([("Intensity", "F1"), ("Intensity", "F2")], names=["Set", "Fraction"]),
    ...                                   index=pd.MultiIndex.from_tuples([("foo", "lorem", np.nan), ("bar", "ipsum", ""), ("baz", "dolor", "+")], names=["Protein IDs", "Gene names", "Reverse"])),
    ...                      column="Reverse")
    Set                            Intensity   
    Fraction                              F1 F2
    Protein IDs Gene names Reverse             
    foo         lorem      NaN             1  2
    bar         ipsum                      3  4

    >>> filter_singlecolumn_keep(pd.DataFrame([[1,2],[3,4],[5,6]],
    ...                                   columns=pd.MultiIndex.from_tuples([("Intensity", "F1"), ("Intensity", "F2")], names=["Set", "Fraction"]),
    ...                                   index=pd.MultiIndex.from_tuples([("foo", "lorem", np.nan), ("bar", "ipsum", 0.5), ("baz", "dolor", 1)], names=["Protein IDs", "Gene names", "Score"])),
    ...                      column="Score", operator=">=", value=1)
    Set                    Intensity   
    Fraction                      F1 F2
    Protein IDs Gene names             
    baz         dolor              5  6
    """
    # fmt: on

    df_filtered = df.query(f"`{column}` {operator} {value}")
    if column in df_filtered.columns:
        if len(set(df_filtered[column])) == 1:
            df_filtered.drop(column, axis=1, inplace=True)
    else:
        if len(set(df_filtered.index.get_level_values(column))) == 1:
            df_filtered.reset_index(column, drop=True, inplace=True)

    return df_filtered


def transform_data(df, unlog=2, invert=False, log=False):
    """
    Perform basic data transformations: invert, logarithmize or unlogarithmize.

    Arguments:
        df: pd.DataFrame
        invert: boolean
        unlog: number or False
        log: one of [2, 10, e] or False

    Returns:
        df_transformed: pd.DataFrame

    Doctests for the different transformations:
    # log data
    >>> transform_data(pd.DataFrame([[0,1,np.nan,10,-2.5]]), unlog=False, invert=False, log="e")
         0    1   2         3   4
    0 -inf  0.0 NaN  2.302585 NaN

    # shift data to a different log
    >>> transform_data(pd.DataFrame([[0,1,np.nan,10,-2.5]]), unlog=2, invert=False, log=10)
         0        1   2       3         4
    0  0.0  0.30103 NaN  3.0103 -0.752575

    # invert and log data
    >>> transform_data(pd.DataFrame([[0,1,np.nan,10,-2.5]]), unlog=False, invert=True, log=2)
         0    1   2         3   4
    0  inf  0.0 NaN -3.321928 NaN
    """
    if unlog:
        df_transformed = df.transform(lambda x: unlog**x)
    else:
        df_transformed = df.copy()

    if invert:
        df_transformed = df_transformed.transform(lambda x: 1 / x)

    if log == 2:
        df_transformed = df_transformed.transform(lambda x: np.log2(x))
    elif log == 10:
        df_transformed = df_transformed.transform(lambda x: np.log10(x))
    elif log == "e":
        df_transformed = df_transformed.transform(lambda x: np.log(x))
    elif not log:
        pass
    else:
        raise ValueError(
            f"Only numpy logarithms are supported (and reasonable), not {log}."
        )

    return df_transformed


def normalize_samples(df, method="sum"):
    """
    Normalization of samples by different methods.

    sum: Normalize all samples to have same (average) sum (assuming equal sample loading). Useful for raw intensities.
    median: Normalize all samples to have a median of 1. Useful for SILAC samples assuming equal channel loading.

    Arguments:
        df: pd.DataFrame
        method: str, one of [sum, median]

    Returns:
        df_transformed: pd.DataFrame

    Doctests for each method:
    >>> normalize_samples(pd.DataFrame([[1,1,1,0],[1,2,3,4],[0,1,2,np.nan]]),method="sum")
         0    1         2    3
    0  2.0  1.0  0.666667  0.0
    1  2.0  2.0  2.000000  4.0
    2  0.0  1.0  1.333333  NaN
    >>> normalize_samples(pd.DataFrame([[1,1,1,0],[1,2,3,4],[0,1,2,np.nan]]),method="median")
         0    1    2    3
    0  1.0  1.0  0.5  0.0
    1  1.0  2.0  1.5  2.0
    2  0.0  1.0  1.0  NaN
    """

    if method == "sum":
        ## Normalize data by ensuring all columns have the same average sum.
        sums = df.sum(axis=0)
        targetsum = np.mean(sums)
        df_transformed = df / sums * targetsum
    elif method == "median":
        df_transformed = df / df.apply(np.nanmedian, axis=0)
    else:
        raise ValueError(f"Unknown method for sample normalization: {method}")

    return df_transformed


# Future
# def weigh_yields(df):
#    """
#
#    """
#    if len(df) == 0:
#        raise ValueError
#    return df_transformed


def normalize_sum1(
    df, sets: dict = {DataFrameStrings.ABUNDANCE: DataFrameStrings.LFQ_INTENSITY}
):
    """ """
    ## Fetch column levels and assert presence of the required sets and fraction annotations
    column_levels = df.columns.names
    # All sets present
    if DataFrameStrings.SET not in column_levels:
        raise KeyError("Column level 'Set' is missing.")
    else:
        if any(
            [
                v not in df.columns.get_level_values(DataFrameStrings.SET)
                for v in sets.values()
            ]
        ):
            raise KeyError(
                f"One of {', '.join(sets.values())} is not in the sets required for 0/1 normalization."
            )
    # Fraction annotation present
    if DataFrameStrings.FRACTION not in column_levels:
        if DataFrameStrings.FRACTION not in df.index.names:
            raise KeyError("Column or index level 'Fraction' is missing.")

    ## Convert to profile dataframe
    if len(column_levels) > 1:
        df_profiles = df.stack(
            [el for el in column_levels if el != DataFrameStrings.SET]
        ).unstack(DataFrameStrings.FRACTION)
    else:
        df_profiles = df.unstack(DataFrameStrings.FRACTION)

    df_01norm = df_profiles[[sets[DataFrameStrings.ABUNDANCE]]].div(
        df_profiles[sets[DataFrameStrings.ABUNDANCE]].sum(axis=1), axis=0
    )
    df_01norm = df_01norm.replace(np.NaN, 0)
    df_01norm.rename(
        {sets[DataFrameStrings.ABUNDANCE]: DataFrameStrings.NORMALIZED_PROFILE},
        axis=1,
        inplace=True,
    )

    df_01 = df_profiles.join(df_01norm)

    return df_01


## The following function is taken from https://doi.org/10.1002/pmic.202100103, and was authored by Schessner et al.
def plot_sample_correlations(
    df: pd.DataFrame,
    data_columns: str = "Intensity (.*)",
    correlation_function: callable = lambda x: np.corrcoef(x.T),
    mode: str = "scatter",
    log10: bool = True,
    binning: int = 10,
):
    """
    Generates either a grid of paired scatter plots, or a heatmap of pairwise correlation values.

    Parameters
    ----------
    df: pandas DataFrame
    data_columns: regular expression string, default = "Intensity (.*)"
        Regular expression matching to columns containing data and a capture group for sample labels.
    correlation_function: callable, default = lambda x: np.corrcoef(x.T)
        Callable function to calculate correlation between samples.
    mode: str, default = "scatter"
        One of 'scatter' and 'heatmap' to swtich modes.
    log10: bool = True
        If True, data is log10 transformed prior to analysis
    binning: int = 10
        Only relevant in scatter mode. Number of bins per full axis unit (e.g. 1e10 Intensity)
        used for datashader density rendering.

    Returns
    -------
    a plotly Figure
        either a subplotted, or a single heatmap depending on the mode
    """
    # pick and process data
    df_sub = df[[el for el in df.columns if re.match(data_columns, el)]].copy()
    if log10:
        df_sub = df_sub.apply(np.log10)
    df_sub = df_sub.replace([np.inf, -np.inf], np.nan)
    df_sub.columns = [re.findall(data_columns, el)[0] for el in df_sub.columns]

    if mode == "scatter":
        # setup subplots and axes
        fig = make_subplots(
            rows=len(df_sub.columns),
            cols=len(df_sub.columns),
            start_cell="bottom-left",
            shared_yaxes=True,
            shared_xaxes=True,
            horizontal_spacing=0.03,
            vertical_spacing=0.03,
        )
        i_range = (
            np.floor(np.nanmin(df_sub)),
            np.ceil(np.nanmax(df_sub)) + 1 / binning,
        )
        j_range = (
            np.floor(np.nanmin(df_sub)),
            np.ceil(np.nanmax(df_sub)) + 1 / binning,
        )
        i_width = int((i_range[1] - i_range[0] - 1 / binning) * binning + 1)
        j_width = int((j_range[1] - j_range[0] - 1 / binning) * binning + 1)

        # fill plots
        for i, ni in enumerate(df_sub.columns):
            for j, nj in enumerate(df_sub.columns):
                # apply datashader
                dc = ds.Canvas(
                    plot_width=i_width,
                    plot_height=j_width,
                    x_range=i_range,
                    y_range=j_range,
                )
                df_ij = (
                    df_sub[[ni, nj]].dropna()
                    if i != j
                    else pd.DataFrame(df_sub[ni].dropna())
                )
                da = dc.points(df_ij, x=ni, y=nj)
                zero_mask = da.values == 0
                da.values = da.values.astype(float)
                da.values[zero_mask] = np.nan

                # add trace
                fig.add_trace(
                    go.Heatmap(
                        z=da, coloraxis="coloraxis1" if i != j else "coloraxis2"
                    ),
                    row=j + 1,
                    col=i + 1,
                )

                # add annotations
                if j == 0:
                    fig.update_xaxes(
                        title_text=ni,
                        row=j + 1,
                        col=i + 1,
                        tickvals=list(range(0, i_width, binning)),
                        ticktext=np.round(da[nj].values[0:i_width:binning]),
                    )
                if i == 0:
                    fig.update_yaxes(
                        title_text=nj,
                        row=j + 1,
                        col=i + 1,
                        tickvals=list(range(0, j_width, binning)),
                        ticktext=np.round(da[ni].values[0:j_width:binning]),
                    )
                if i != j:
                    fig.add_annotation(
                        dict(
                            text=str(
                                np.round(
                                    np.min(
                                        correlation_function(df_sub[[ni, nj]].dropna())
                                    ),
                                    4,
                                )
                            ),
                            x=binning,
                            y=j_width,
                            showarrow=False,
                        ),
                        row=j + 1,
                        col=i + 1,
                    )

        # layout figure
        fig.update_layout(
            template="simple_white",
            coloraxis2=dict(showscale=False, colorscale=["black", "black"]),
            width=i * 200 + 100,
            height=j * 200 + 70,
            margin_t=0,
        )
    elif mode == "heatmap":
        da = np.ones((len(df_sub.columns), len(df_sub.columns)))
        for i, ni in enumerate(df_sub.columns):
            for j, nj in enumerate(df_sub.columns):
                # filter data and store correlation values
                df_ij = (
                    df_sub[[ni, nj]].dropna()
                    if i != j
                    else pd.DataFrame(df_sub[ni].dropna())
                )
                if i != j:
                    da[i, j] = np.round(
                        np.min(correlation_function(df_sub[[ni, nj]].dropna())), 4
                    )
        # create figure and label axes
        fig = go.Figure(data=go.Heatmap(z=da))
        fig.update_xaxes(
            tickvals=list(range(0, i + 1, 1)), ticktext=list(df_sub.columns)
        )
        fig.update_yaxes(
            tickvals=list(range(0, j + 1, 1)), ticktext=list(df_sub.columns)
        )
        fig.update_layout(
            template="simple_white", width=i * 40 + 150, height=j * 40 + 150
        )
    else:
        raise ValueError
    return fig


def calc_width_categorical(f, itemwidth=40, charwidth=7):
    legendlength = max([len(el.name) for el in f.data])
    try:
        items = set([i for el in f.data for i in el.x])
    except:
        items = [tup for el in f.data for tup in zip(*el.x)]
        items = [el for i, el in enumerate(items) if el not in items[:i]]
    xwidth = 80 + len(items) * itemwidth
    legendwidth = 50 + charwidth * legendlength
    width = xwidth + legendwidth
    return dict(
        width=width, xaxis_domain=[0, xwidth / width], legend_x=(xwidth + 5) / width
    )


def outlier_test(
    delta: np.array,
    support_fraction: float = 0.75,
    iterations: int = 10,
    record_medians: bool = True,
    stop_at_95_05: bool = True,
):
    mh = []
    medians = []
    for i in range(iterations):
        with warnings.catch_warnings():
            if i != 0:
                warnings.simplefilter("ignore")
            else:
                warnings.simplefilter(action="ignore", category=RuntimeWarning)
            cov = MinCovDet(support_fraction=support_fraction, random_state=i).fit(
                delta
            )
            mh.append(cov.mahalanobis(delta))
            if record_medians:
                medians.append(np.median(np.asarray(mh), axis=0))
                if stop_at_95_05 and i > 0:
                    q95 = np.quantile(
                        abs(medians[-1] - medians[-2]) / medians[-2], 0.95
                    )
                    if q95 < 0.005:
                        break
    if not record_medians:
        medians.append(np.median(np.asarray(mh), axis=0))
    pvals = chi2.sf(medians[-1], delta.shape[1])
    return pvals, medians


def rscores_multi(df: pd.DataFrame, pairs: list = []):
    paired_pvals = pd.DataFrame(index=df.index)

    if pairs == []:
        if DataFrameStrings.MAP in df.columns.names:
            pairs = [
                el
                for el in combinations(
                    np.unique(df.columns.get_level_values(DataFrameStrings.MAP)), 2
                )
            ]
    print([el for el in pairs])

    fractions = df[pairs[0][0]].columns.get_level_values(DataFrameStrings.FRACTION)
    delta = pd.DataFrame(
        df[[p[0] for p in pairs]].values - df[[p[1] for p in pairs]].values,
        columns=pd.MultiIndex.from_tuples(
            [
                (p, f)
                for p, f in zip(
                    np.repeat(["-".join(pair) for pair in pairs], len(fractions)),
                    cycle(fractions),
                )
            ],
            names=["Pair", DataFrameStrings.FRACTION],
        ),
    )
    print(delta.head())

    # Preparation of new dataframe
    reshapedRow = reshapeRow(
        delta.iloc[0, :], replicates="Pair", variables=DataFrameStrings.FRACTION
    )
    nameMatrix = np.array(
        [
            ["{}_{}".format(r1, r2) for r2 in np.unique(reshapedRow.index)]
            for r1 in np.unique(reshapedRow.index)
        ]
    )

    # Calculation of correlations
    correlations = delta.apply(
        lambda x: pd.Series(
            correlation_row(
                reshapeRow(x, replicates="Pair", variables=DataFrameStrings.FRACTION),
                method="cosine",
            )
        ),
        axis=1,
    )

    # Returning output table
    correlations.columns = nameMatrix[np.tril_indices(len(nameMatrix), k=-1)]
    correlations.index = df.index

    return correlations


def calculate_correlation(m, method="pearson"):
    methods = {
        "pearson": lambda x: x.T.corr(method="pearson").values,
        "spearman": lambda x: x.T.corr(method="spearman").values,
        "kendall": lambda x: x.T.corr(method="kendall").values,
        "cosine": lambda x: cosine_similarity(x),
    }
    return methods[method](m)


# get lower triangle from the correlation matrix
def correlation_row(m, method="pearson"):
    cm = calculate_correlation(m, method=method)
    cv = cm[np.tril_indices(len(cm), k=-1)]
    return cv


# pivot a single df line for correlation matrix calculation
def reshapeRow(
    x, replicates="Replicate", variables=["Condition", DataFrameStrings.FRACTION]
):
    x = x.reset_index()
    x.columns = np.append(x.columns.values[0:-1], ["value"])
    if type(variables) == list:
        x["variables"] = [i + "_" + j for i, j in zip(x[variables[0]], x[variables[1]])]
        variables = "variables"
    x = x.pivot(index=replicates, columns=variables, values="value")
    return x


def subsample_misclassifcation(mc, sample_size=0.75, random_state=42):
    mc = mc.copy()
    mc.index.name = "correct"
    mc.columns.name = "predicted"

    prots = pd.DataFrame(columns=["correct", "predicted"])
    for el in mc.stack("predicted").iteritems():
        for _ in range(int(el[1])):
            prots.loc[len(prots), :] = [el[0][0], el[0][1]]

    if sample_size != 1:
        _, sample = train_test_split(
            prots,
            stratify=prots["correct"],
            test_size=sample_size,
            random_state=random_state,
        )
    else:
        sample = prots.copy()
    mc1 = sample.groupby(["correct", "predicted"]).size().unstack("predicted")
    confusion_formatted = pd.DataFrame(index=mc.index.values, columns=mc.index.values)
    for i in confusion_formatted.columns:
        for j in confusion_formatted.index:
            try:
                confusion_formatted.loc[j, i] = mc1.loc[j, i]
            except:
                pass
    return confusion_formatted


def process_mc(confusion):
    true_positives = np.diag(confusion)
    class_sizes = confusion.sum(axis=1)
    mask_true = np.ones(confusion.shape)
    mask_true[np.diag_indices(confusion.shape[0])] = 0
    false_predicitons = confusion * mask_true
    false_positives = false_predicitons.sum(axis=0)
    false_negatives = false_predicitons.sum(axis=1)
    recall = true_positives / (true_positives + false_negatives)
    precision = true_positives / (true_positives + false_positives)
    F1_scores = pd.Series(
        [statistics.harmonic_mean([p, r]) for p, r in zip(precision, recall)],
        index=confusion.index,
    )
    overall = true_positives.sum() / (true_positives.sum() + false_negatives.sum())
    recall["Average all classes"] = recall.mean()
    precision["Average all classes"] = precision.mean()
    F1_scores["Average all classes"] = F1_scores.mean()
    class_sizes["Average all classes"] = class_sizes.mean()
    recall["Overall performance"] = overall
    precision["Overall performance"] = overall
    F1_scores["Overall performance"] = overall
    class_sizes["Overall performance"] = class_sizes.sum()
    # Average performance across membranous organelles
    organelles = [
        "Plasma membrane",
        "Peroxisome",
        "Mitochondrion",
        "Lysosome",
        "ER",
        "Endosome",
        "Golgi",
        "Ergic/cisGolgi",
        "ER_high_curvature",
        "PM",
        "endosome",
        "Nucleus",
    ]
    recall["Average membraneous organelles"] = recall[
        [el in organelles for el in recall.index]
    ].mean()
    precision["Average membraneous organelles"] = precision[
        [el in organelles for el in precision.index]
    ].mean()
    F1_scores["Average membraneous organelles"] = F1_scores[
        [el in organelles for el in F1_scores.index]
    ].mean()
    class_sizes["Average membraneous organelles"] = class_sizes[
        [el in organelles for el in class_sizes.index]
    ].mean()
    return pd.DataFrame(
        np.array([F1_scores, recall, precision, class_sizes]).T,
        columns=["F1 score", "Recall", "Precision", "Class size"],
        index=pd.Index(F1_scores.index, name=DataFrameStrings.COMPARTMENT),
    )


def process_mc_errorbars(mc, n=20, sample_size=0.75):
    score_iterations = pd.DataFrame(
        columns=["F1 score", "Recall", "Precision", "Class size"],
        index=pd.MultiIndex.from_tuples(
            [], names=[DataFrameStrings.COMPARTMENT, "iteration"]
        ),
    )

    for i in range(n):
        confusion = subsample_misclassifcation(
            mc=mc, random_state=i, sample_size=sample_size
        ).fillna(0)
        score_iterations = score_iterations.append(
            process_mc(confusion).set_index(
                pd.Index(np.repeat(i, len(confusion) + 3), name="iteration"),
                append=True,
            )
        )

    scores = (
        score_iterations.fillna(0)
        .groupby([DataFrameStrings.COMPARTMENT])
        .apply(
            lambda x: pd.DataFrame(
                [x.mean(), x.std()], index=pd.Index(["mean", "std"], name="value")
            )
        )
    )
    scores = scores.append(
        process_mc(mc).set_index(
            pd.Index(np.repeat("direct", len(confusion) + 3), name="value"), append=True
        )
    )
    return scores
