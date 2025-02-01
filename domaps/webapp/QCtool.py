#!/usr/bin/env python
# coding: utf-8

# ## Table of contents<a id="TOC"></a>
# 1. [Library imports](#libraries)
# 2. [panel and plotly settings and customization](#styling)
# 3. [Global variables](#globals)
# 4. [App interface skeleton](#skeleton)
# 5. [App serving](#serving)
# 6. [Cell structuring](#structuring)
# 7. [Home tab](#home)
# 8. [Analysis tab](#analysis)
# 9. [Comparison output tabs](#comparison)
# 10. [Benchmark tab upload section](#benchmarkupload)
# 11. [Code interactions](#interactions)

# [<div style="text-align: right; font-size: 8pt">back to top</div>](#TOC)
# ## Library imports<a id="libraries"></a>

# In[ ]:


import natsort
import numpy as np
import pandas as pd
import param
import re
import traceback
import panel as pn

pn.extension()
pn.extension("plotly")
import io
from io import BytesIO
from io import StringIO
from bokeh.models.widgets.tables import NumberFormatter
import plotly.express as px
import json
import os
from importlib import reload
import pkg_resources
import time
import copy
from scipy.stats import spearmanr

import domaps.gui as gui
from domaps.network import format_data

import domaps

from datetime import datetime


# [<div style="text-align: right; font-size: 8pt">back to top</div>](#TOC)
# ## panel and plotly settings and customization<a id="styling"></a>

# In[ ]:


css = """
.detail_menu .bk-headers-wrapper{
  border-bottom: 2px solid #0a5da8 !important;
  margin-bottom: 10px;
  min-width: 800px;
  max-width: 2000px !important;
}
.detail_menu .bk-tab{
  color: #0a5da8;
}
.detail_menu .bk-active{
  border-color: #0a5da8 !important;
  border-width: 2px 1px 0px 1px !important;
  color: #111111 !important;
}

.main_menu .bk-headers-wrapper{
  border-bottom: 2px solid #0a5da8 !important;
  margin-bottom: 10px;
  min-width: 800px;
  max-width: 2000px !important;
}
.main_menu .bk-tab{
  color: #0a5da8;
  font-size: 120%;
  font-weight: bold;
}
.main_menu .bk-active{
  border-color: #0a5da8 !important;
  border-width: 2px 1px 0px 1px !important;
  color: #111111 !important;
}

.bk-tabs-header{
  min-width: 1000px;
}

.card-title:first-child{
  font-size: 13px;
}

.button-main .bk-btn{
  font-size: 120%;
  font-weight: bold;
}

body{
  margin:0;
}
"""
pn.extension(raw_css=[css])

plotly_config = {
    "toImageButtonOptions": {
        "format": "svg",  # one of png, svg, jpeg, webp
        "filename": "figure",
    }
}


def resize(el):
    try:
        el.append(None)
        el.pop(-1)
    except:
        pass


# [<div style="text-align: right; font-size: 8pt">back to top</div>](#TOC)
# ## Global variables<a id="globals"></a>

# In[ ]:


#### Data that will be stored in memory during each session
## dataset analysed and displayed in the single analysis tab
mem_single_analysis = None
i_class = None
## dataset collection analysed and displayed in the benchmarking tab
mem_benchmark = None
i_class_comp = None
## currently available datasets to select for benchmarking
mem_available_datasets = dict()

with open("textfragments.json", "r") as file:
    textfragments = json.load(file)

DEBUG = False
MAX_SIZE_MB = 80
CONTENT_WIDTH = 1000
HOSTING = (
    "online"  # switch to local to remove limits on iterations and cross-validation
)


# [<div style="text-align: right; font-size: 8pt">back to top</div>](#TOC)
# ## App interface skeleton<a id="skeleton"></a>

# In[ ]:


#### Served panel object:
app = pn.GridSpec(sizing_mode="stretch_both", margin=0)
app[0, 0] = pn.Spacer(background="#DDDDDD", margin=0)  #
app[0, 16] = pn.Spacer(background="#DDDDDD", margin=0)  #

#### Insert main content container
app_center = pn.Column(
    pn.Row(
        pn.Pane("# QC tool for Spatial Proteomics", width=600),
        pn.layout.HSpacer(),
        pn.Pane(f"domaps version {domaps.__version__}", width=200),
        margin=10,
    ),
    pn.Row(name="main_content", margin=50),
)
app[0, 1:15] = app_center

#### Insert main menu tab object
app_tabs = pn.Tabs(
    margin=10, css_classes=["main_menu"], dynamic=True, sizing_mode="stretch_width"
)
app_center.objects[[i.name for i in app_center].index("main_content")] = app_tabs

#### Append individual dashboards
## Home
dashboard_home = pn.Column(
    "Interface loading ...", name="home", sizing_mode="stretch_width"
)
app_tabs.append(("Home", dashboard_home))

## Single analysis
dashboard_analysis = pn.Column(
    "Interface loading ...", name="analysis", sizing_mode="stretch_width"
)
app_tabs.append(("Analysis", dashboard_analysis))
analysis_tabs = pn.Tabs(
    margin=10, css_classes=["detail_menu"], dynamic=True, sizing_mode="stretch_width"
)

## Benchmark
dashboard_benchmark = pn.Column(
    "Interface loading ...", name="benchmark", sizing_mode="stretch_width"
)
app_tabs.append(("Benchmark", dashboard_benchmark))
comparison_tabs = pn.Tabs(
    margin=10, css_classes=["detail_menu"], dynamic=True, sizing_mode="stretch_width"
)

## About
app_tabs.append(
    (
        "About",
        pn.Row(pn.Pane(textfragments["about_intro"], sizing_mode="stretch_width")),
    )
)


# [<div style="text-align: right; font-size: 8pt">back to top</div>](#TOC)
# ## App serving<a id="serving"></a>
# Switch cells below between markup and code to set up for server hosting from the command line (app.servable) vs. local hosting from python.

# try:
#     server.stop()
# except Exception:
#     print("First server startup")
# server = app.show(port=5067, websocket_max_message_size=MAX_SIZE_MB*1024*1024, admin=True,
#                   http_server_kwargs={'max_buffer_size': MAX_SIZE_MB*1024*1024})

# In[ ]:


app.servable()


# [<div style="text-align: right; font-size: 8pt">back to top</div>](#TOC)
# ## Cell structuring<a id="structuring"></a>
# All cells below contain one or several sets of these points:
# - (Dashboard structure)
# - Layout and widget elements (outside in)
# - Layout assembly and appending
# - Callback definitions
# - Positioning of callback outputs

# In[ ]:


#### Dashboard structure
########################

#### Layout elements
####################

#### Append layout to dashboard
###############################

#### Callbacks
##############
# list
# of
# callbacks

#### Callback output positioning
################################


# [<div style="text-align: right; font-size: 8pt">back to top</div>](#TOC)
# ## Home tab<a id="home"></a>

# In[ ]:


#### Dashboard structure
########################
# already defined as single column
dashboard_home.objects = []

#### Layout elements
####################
lo_home_intro = pn.Pane(textfragments["home_intro"], width=CONTENT_WIDTH)
btn_home_analysesingle = pn.widgets.Button(
    name="Format and analyse single experiment",
    button_type="success",
    width=400,
    height=50,
    css_classes=["button-main"],
)
lo_home_singleinstructions = pn.Column(
    pn.Pane(textfragments["home_single_shortinstructions"], width=CONTENT_WIDTH),
    pn.Card(
        textfragments["quick_start_guide"],
        header="DOM-QC 1-min Quick Start Guide",
        width=CONTENT_WIDTH,
        name="tutorial_single",
        collapsed=True,
    ),
)
btn_home_benchmark = pn.widgets.Button(
    name="Benchmark multiple experiments",
    button_type="success",
    width=400,
    height=50,
    css_classes=["button-main"],
)
lo_home_benchmarkinstructions = pn.Column(
    pn.Pane(textfragments["home_benchmark_shortinstructions"], width=CONTENT_WIDTH),
    #    pn.Card("Add screenshot tutorial here.", header="Tutorial", width=CONTENT_WIDTH,
    #            name="tutorial_benchmark", collapsed=True)
)

#### Append layout to dashboard
###############################
for el in [
    lo_home_intro,
    btn_home_analysesingle,
    lo_home_singleinstructions,
    btn_home_benchmark,
    lo_home_benchmarkinstructions,
]:
    dashboard_home.append(el)

#### Callbacks
##############
# home_gotosingleanalysis
# home_gotobenchmark


def home_gotosingleanalysis(event):
    app_tabs.active = [el.name for el in app_tabs].index("analysis")


btn_home_analysesingle.on_click(home_gotosingleanalysis)


def home_gotobenchmark(event):
    app_tabs.active = [el.name for el in app_tabs].index("benchmark")


btn_home_benchmark.on_click(home_gotobenchmark)

#### Callback output positioning
################################
# None


# [<div style="text-align: right; font-size: 8pt">back to top</div>](#TOC)
# ## Analysis tab<a id="analysis"></a>
# - File config
# - Analysis output

# In[ ]:


#### Dashboard structure
########################
dashboard_analysis.objects = [
    pn.Row(name="file_config"),
    pn.Column(name="analysis_output", sizing_mode="stretch_width"),
]

#### Layout elements
#### File config
####################
# lo_read_file = pn.Card(header="### Upload configuration", min_width=400)
lo_config_instructions = pn.Card(
    pn.Pane(textfragments["upload_instructions"]), header="### Instructions", width=400
)
lo_config_details = pn.Card(
    pn.Pane(textfragments["upload_details"]),
    header="### Details on configuring your data",
    width=400,
    collapsed=True,
)
lo_config_error_messages = pn.Card(
    pn.Pane(textfragments["upload_error_messages"]),
    header="### Common error messages",
    width=400,
    collapsed=True,
)

loading_status = pn.Row()
idle = pn.indicators.LoadingSpinner(value=False, width=100, height=100, color="primary")
loading = pn.indicators.LoadingSpinner(
    value=True, width=100, height=100, color="primary"
)
analysis_status = pn.Pane("", width=300)
filereading_status = pn.Pane("No data import yet", width=300)
i_FileConfig = gui.ConfigureSingleFile(width=540)
# lo_read_file.append(i_FileConfig)

#### Append layout to dashboard
#### File config
###############################
dashboard_analysis.objects[
    [el.name for el in dashboard_analysis].index("file_config")
].objects = []
for el in [
    pn.Column(i_FileConfig, analysis_status, loading_status, width=600),
    pn.Column(
        lo_config_instructions,
        lo_config_details,
        # lo_config_error_messages # currently empty
    ),
]:
    dashboard_analysis.objects[
        [el.name for el in dashboard_analysis].index("file_config")
    ].append(el)

#### Callbacks
#### File config
################
# future_execution

cache_uploaded = pn.widgets.Checkbox(value=False)
cache_run = pn.widgets.Checkbox(value=False)
# define widgets that should be disbled after run==True
# wdgts = [i_FileConfig]


def execution(event):
    loading_status.objects = [loading]
    analysis_status.object = "Analysis in progress"
    lo_config_instructions.collapsed = True
    lo_config_details.collapsed = True
    # lo_config_error_messages.collapsed = True
    # lo_read_file.collapsed = True
    output_layoutpos = [el.name for el in dashboard_analysis].index("analysis_output")
    dashboard_analysis.objects[output_layoutpos].objects = []
    cache_run.value = False
    try:
        global i_class
        i_class = domaps.SpatialDataSet.from_settings(i_FileConfig.get_settings())
        i_class.run_pipeline(
            content=BytesIO(i_FileConfig._content.file.value),
            progressbar=analysis_status,
        )

        analysis_status.object = "Analysis finished!"
        update_object_selector()
        loading_status.objects = []
        dashboard_analysis.objects[output_layoutpos].append(analysis_tabs)
        mem_available_datasets[i_class.expname] = i_class.analysed_datasets_dict[
            i_class.expname
        ]
        try:
            i_dfs_available.options = i_dfs_available.options + [i_class.expname]
            i_dfs_available.value = i_dfs_available.value + [i_class.expname]
            coll_activatebuttons(i_dfs_available.value)
            resize(lo_dfs_available)
        except:
            pass
        cache_run.value = True
    except:
        loading_status.objects = [""]
        analysis_status.object = traceback.format_exc()
        cache_run.value = False


i_FileConfig._btn_run.on_click(execution)


#### Dashboard structure
#### Analysis output
########################

#### Layout elements
####################

i_logOR01_selection = pn.widgets.Select(
    options=[
        "0/1 normalized data",
        "log transformed data",
        "stringency filtered raw data",
        "Overview - cluster distances",
    ],
    name="Select type of data for download",
    width=300,
)

i_clusterwidget = pn.widgets.Select(
    options=["Proteasome", "Lysosome"], name="Cluster of interest", width=300
)
i_mapwidget = pn.widgets.Select(
    options=["Map1", "Map2"], name="Map of interest", width=300
)

i_collapse_maps_PCA = pn.widgets.Checkbox(value=False, name="Collapse maps")
i_pca_ncomp = pn.widgets.IntSlider(value=3, start=2, end=5)
btn_pca_ncomp = pn.widgets.Button(
    name="Recalculate PCA with different number of components."
)

#### Correlation plot in single analysis tab
i_corr_mode = pn.widgets.Select(
    options=["heatmap", "scatter"], name="Show correlations as"
)
i_corr_measure = pn.widgets.Select(
    options=["Pearson", "Spearman"], name="Correlation metric"
)
i_corr_selection = pn.widgets.Select(options=["all samples"], name="Select fraction(s)")
i_corr_level = pn.widgets.Select(
    options=["original data", "filtered and processed data", "normalized profiles"],
    value="filtered and processed data",
    name="Show dataset",
)
lo_corr_widgets = pn.WidgetBox(
    i_corr_mode, i_corr_measure, i_corr_selection, i_corr_level
)

#### Append layout to dashboard
###############################

#### Callbacks
##############
# update_corr_selection
# show_correlation_single
# update_object_selector
# recalculate_pca
# update_data_overview # tab content
# update_cluster_overview # tab content
# update_cluster_details
# update_quantity # tab content
# show_tabular_overview
# json_download
# df01_download_widget
# df01_download
# dflog_download
# df_filteredRawData_download
# table_download # tab content


@param.depends(cache_run.param.value, watch=True)
def update_corr_selection(run):
    if run == True:
        i_corr_selection.options = [i_corr_selection.options[0]] + i_class.fractions


@param.depends(
    i_corr_mode.param.value,
    i_corr_measure.param.value,
    i_corr_selection.param.value,
    i_corr_level.param.value,
    cache_run.param.value,
)
def show_correlation_single(mode, measure, selection, level, run):
    if run == True:
        try:
            if mode == "scatter" and selection == "all samples":
                return "Please select only one fraction in scatter mode."
            if level == "original data":
                df = i_class.df_index[i_class.mainset]
            elif level == "filtered and processed data":
                df = i_class.df_log_stacked["log profile"].unstack(["Map", "Fraction"])
            elif level == "normalized profiles":
                df = i_class.df_01_stacked["normalized profile"].unstack(
                    ["Map", "Fraction"]
                )
            else:
                raise ValueError(f"Unknown data level {level} for correlation plot.")
            df = df.sort_index(axis=1, level="Fraction")
            if selection != "all samples":
                df = df.xs(selection, level="Fraction", axis=1, drop_level=False)

            df.columns = ["_".join(el) for el in df.columns]

            measure_dict = {
                "Spearman": lambda x: spearmanr(x.values).correlation,
                "Pearson": lambda x: np.corrcoef(x.T),
            }

            plot = domaps.plot_sample_correlations(
                df,
                data_columns="(.*)",
                log10=False if level != "original data" else True,
                binning=10 if level != "normalized profiles" else 50,
                mode=mode,
                correlation_function=measure_dict[measure],
            )
            return plot

        except:
            return traceback.format_exc()


lo_corr = pn.Row(lo_corr_widgets, show_correlation_single)


def update_object_selector():
    i_mapwidget.options = list(i_class.map_names)
    i_clusterwidget.options = list(i_class.markerproteins.keys())
    i_pca_ncomp.end = len(i_class.fractions)


def recalculate_pca(event):
    cache_run.value = False
    i_class.perform_pca(n=i_pca_ncomp.value)
    cache_run.value = True


btn_pca_ncomp.on_click(recalculate_pca)


@pn.depends(cache_run.param.value, i_collapse_maps_PCA.param.value)
def update_data_overview(run, collapse_maps_PCA):
    try:
        if run == True:
            compartments = i_class.df_organellarMarkerSet["Compartment"].unique()
            compartment_color = dict(zip(compartments, i_class.css_color))
            compartment_color["undefined"] = "lightgrey"

            pca_plot = gui.pca_plot(
                df_pca=i_class.df_pca
                if not collapse_maps_PCA
                else i_class.df_pca_combined,
                df_var=i_class.df_pca_var
                if not collapse_maps_PCA
                else i_class.df_pca_combined_var,
                df_loadings=i_class.df_pca_loadings
                if not collapse_maps_PCA
                else i_class.df_pca_combined_loadings,
                color="Compartment",
                color_map=compartment_color,
                facet_col="Map" if not collapse_maps_PCA else None,
                title="PCA plot",
            )
            pca_plot.highlight_dict = i_class.markerproteins

            log_histogram = i_class.plot_log_data()
            visualization_map = pn.Column(
                pn.Pane(textfragments["analysis_overview_top"], width=600),
                pn.Row(i_collapse_maps_PCA, i_pca_ncomp, btn_pca_ncomp),
                pca_plot,
                pn.layout.VSpacer(height=3, background="#AAAAAA"),
                lo_corr,
                pn.layout.VSpacer(height=3, background="#AAAAAA"),
                pn.Row(pn.Pane(log_histogram, width=1000, config=plotly_config)),
            )
            app_tabs.active = 1
            return visualization_map
        else:
            visualization_map = "Run analysis first!"
            return visualization_map
    except Exception:
        update_status = traceback.format_exc()
        return pn.Column(pn.Row(i_collapse_maps_PCA), update_status)


@pn.depends(i_clusterwidget.param.value, i_mapwidget.param.value, cache_run.param.value)
def update_cluster_overview(clusterwidget, mapwidget, run):
    try:
        if run == True:
            list_genes = [
                goi
                for goi in i_class.genenames_sortedout_list
                if goi in i_class.markerproteins[clusterwidget]
            ]
            i_class.cache_cluster_quantified = True
            distance_boxplot = i_class.distance_boxplot(
                cluster_of_interest=clusterwidget
            )
            if i_class.cache_cluster_quantified == False:
                return pn.Column(
                    pn.Row(
                        pn.Pane(textfragments["analysis_intramap_top"], width=600),
                        pn.Column(width=30),
                        pn.Column(i_clusterwidget, i_mapwidget),
                    ),
                    "This protein cluster was not quantified",
                )

            else:
                df_quantification_overview = i_class.quantification_overview(
                    cluster_of_interest=clusterwidget
                )
                profiles_plot = i_class.profiles_plot(
                    map_of_interest=mapwidget, cluster_of_interest=clusterwidget
                )
                pca_plot = gui.pca_plot(
                    df_pca=i_class.df_pca_all_marker_cluster_maps.xs(
                        clusterwidget, level="Cluster", axis=0
                    ),
                    show_variability=False,
                    show_loadings=False,
                    enable_highlight=False,
                    facet_col=None,
                    color="Map",
                    color_map=dict(),
                    title=f"PCA plot of {clusterwidget}",
                )
                pca_plot._dimensions.value = "3D"
                cluster_overview = pn.Column(
                    pn.Row(
                        pn.Pane(textfragments["analysis_intramap_top"], width=600),
                        pn.Column(width=30),
                        pn.Column(i_clusterwidget, i_mapwidget),
                    ),
                    pn.Row(
                        pca_plot,
                        pn.Pane(profiles_plot, width=500, config=plotly_config),
                        pn.Pane(distance_boxplot, width=500, config=plotly_config),
                    ),
                    pn.layout.VSpacer(height=3, background="#AAAAAA"),
                    pn.Row(
                        pn.Pane(
                            "In total {} proteins across all maps were quantified, whereas the following proteins were not consistently quantified throughout all maps: {}".format(
                                i_class.proteins_quantified_across_all_maps,
                                ", ".join(list_genes),
                            )
                            if len(list_genes) != 0
                            else "All genes from this cluster are quantified in all maps."
                        ),
                        width=1000,
                    ),
                    pn.Row(
                        pn.widgets.DataFrame(
                            df_quantification_overview,
                            height=200,
                            width=500,
                            disabled=True,
                        )
                    ),
                    pn.layout.VSpacer(height=3, background="#AAAAAA"),
                    update_cluster_details,
                )
                return cluster_overview

        else:
            cluster_overview = "Run analysis first!"
            return cluster_overview
    except Exception:
        update_status = pn.Column(
            pn.Row(
                pn.Pane(textfragments["analysis_intramap_top"], width=600),
                pn.Column(width=30),
                pn.Column(i_clusterwidget, i_mapwidget),
            ),
            pn.Pane(traceback.format_exc(), width=600),
        )
        return update_status


@pn.depends(i_clusterwidget.param.value, cache_run.param.value)
def update_cluster_details(clusterwidget, run):
    try:
        if run == True:
            cluster_details = i_class.distance_to_median_boxplot(
                cluster_of_interest=clusterwidget
            )
            return pn.Pane(cluster_details, width=1000, config=plotly_config)
        else:
            cluster_details = "Run analysis first!"
            return pn.Pane(cluster_details, width=1000)
    except Exception:
        update_status = traceback.format_exc()
        return update_status


@pn.depends(cache_run.param.value)
def update_quantity(run):
    try:
        if run == True:
            fig_npg, fig_npr, fig_npr_dc, fig_npg_F, fig_npgf_F, fig_npg_F_dc = (
                i_class.plot_quantity_profiles_proteinGroups()
            )
            return pn.Column(
                pn.Pane(textfragments["analysis_depth_top"], width=600),
                pn.Row(
                    pn.Pane(fig_npg, config=plotly_config),
                    # pn.Pane(fig_npr, config=plotly_config),
                    pn.Pane(fig_npr_dc, config=plotly_config),
                ),
                pn.layout.VSpacer(background="#AAAAAA", height=3),
                pn.Row(
                    pn.Pane(fig_npg_F, config=plotly_config),
                    pn.Pane(fig_npgf_F, config=plotly_config),
                    pn.Pane(fig_npg_F_dc, config=plotly_config),
                ),
            )
        else:
            return "Run analysis first!"
    except Exception:
        update_status = traceback.format_exc()
        return update_status


@pn.depends(cache_run.param.value)
def show_tabular_overview(run):
    try:
        if run == True:
            content = pn.Column(
                pn.widgets.DataFrame(
                    pd.read_json(
                        i_class.analysed_datasets_dict[i_class.expname][
                            "Overview table"
                        ]
                    ),
                    height=200,
                    width=600,
                    disabled=True,
                ),
                i_logOR01_selection,
                df01_download_widget,
                pn.widgets.FileDownload(
                    callback=json_download, filename="AnalysedDatasets.json"
                ),
            )
            return pn.Row(content, width=1000)
        else:
            content = "Please, upload a file first and press ‘Analyse clusters’"
            return pn.Row(content, width=1000)
    except Exception:
        content = traceback.format_exc()
        return pn.Row(content, width=1000)


@pn.depends(cache_run.param.value)
def json_download(run):
    sio = StringIO()
    json.dump(i_class.analysed_datasets_dict, sio, indent=4, sort_keys=True)
    sio.seek(0)
    return sio


@pn.depends(cache_run.param.value, i_logOR01_selection.param.value)
def df01_download_widget(run, logOR01_selection):
    if logOR01_selection == "0/1 normalized data":
        return pn.Column(
            pn.widgets.FileDownload(
                callback=df01_download, filename="01_normalized_data.csv"
            ),
            width=650,
        )
    if logOR01_selection == "log transformed data":
        return pn.Column(
            pn.widgets.FileDownload(
                callback=dflog_download, filename="log_transformed_data.csv"
            ),
            width=650,
        )
    if logOR01_selection == "Overview - cluster distances":
        return pn.Column(
            pn.widgets.FileDownload(
                callback=table_download, filename="cluster_distances.csv"
            ),
            width=650,
        )
    else:
        return pn.Column(
            pn.widgets.FileDownload(
                callback=df_filteredRawData_download,
                filename="stringency_filtered_raw_data.csv",
            ),
            width=650,
        )


@pn.depends(cache_run.param.value)
def df01_download(run):
    df_01 = i_class.reframe_df_01ORlog_for_Perseus(i_class.df_01_stacked)
    sio = StringIO()
    df_01.to_csv(sio)
    sio.seek(0)
    return sio


@pn.depends(cache_run.param.value)
def dflog_download(run):
    df_log = i_class.reframe_df_01ORlog_for_Perseus(i_class.df_log_stacked)
    sio = StringIO()
    df_log.to_csv(sio)
    sio.seek(0)
    return sio


@pn.depends(cache_run.param.value)
def df_filteredRawData_download(run):
    df = i_class.reframe_df_01ORlog_for_Perseus(
        i_class.df_filtered.stack(["Map", "Fraction"])
    )
    sio = StringIO()
    df.to_csv(sio)
    sio.seek(0)
    return sio


@pn.depends(cache_run.param.value)
def table_download(run):
    df = i_class.results_overview_table()
    sio = StringIO()
    df.to_csv(sio)
    sio.seek(0)
    return sio


#### Dashboard structure
#### Movement analysis
########################

lo_movement_analysis = pn.Column(
    pn.panel(textfragments["mr_top"], width=600),
    name="movement_analysis",
    sizing_mode="stretch_width",
)
btn_calc_mr = pn.widgets.Button(
    name="Calculate/display MR-plot", width=300, button_type="success"
)


def calculate_mr(event):
    lo_movement_analysis.objects = lo_movement_analysis.objects[0:2]
    MRstatus = pn.Card(title="Temporary processing data")
    lo_movement_analysis.append(MRstatus)
    try:
        mr = i_class.run_outliertest(**i_MRConfig.get_settings(), canvas=MRstatus)
    except:
        MRstatus.append(traceback.format_exc())
    try:
        MRstatus.append("Loading M-R plot display ...")
        lo_movement_analysis.append(gui.MRPlot(data=mr))
        MRstatus.pop(-1)
    except:
        MRstatus.append(traceback.format_exc())


btn_calc_mr.on_click(calculate_mr)


@pn.depends(cache_run.param.value)
def update_MRConfig(run):
    if run:
        try:
            global i_MRConfig
            i_MRConfig = gui.ConfigureMR(
                conditions=i_class.conditions,
                maps={
                    c: [
                        el.split("_")[1] for el in i_class.map_names if el.startswith(c)
                    ]
                    for c in i_class.conditions
                },
            )
            return pn.Column(i_MRConfig, btn_calc_mr)
        except:
            if i_class.conditions == [""]:
                return pn.Column("**This analysis does not contain any conditions**")
            else:
                return pn.Column(
                    i_class.conditions,
                    pn.Str(traceback.format_exc(), width=600),
                )
    else:
        lo_movement_analysis.objects = lo_movement_analysis.objects[0:2]
        return "Run analysis first"


lo_movement_analysis.append(update_MRConfig)


#### Dashboard structure
#### Neighborhood
########################
lo_neighborhood_analysis = pn.Column(
    pn.panel(textfragments["neighborhood_top"], width=600),
    name="neighborhood_analysis",
    sizing_mode="stretch_width",
)

#### Layout elements
#### Neighborhood
####################
btn_format_neighborhood = pn.widgets.Button(
    name="Format data for neighborhood analysis", width=500, button_type="success"
)

#### Append layout to dashboard
#### Neighborhood
###############################

#### Callbacks
#### Neighborhood
##############
# format_neighborhood
# reset_neighborhood


def format_neighborhood(event):
    lo_neighborhood_analysis.objects = lo_neighborhood_analysis.objects[0:2]
    df_core, meta_dict = format_data(i_class)
    lo_neighborhood_analysis.append(
        gui.NeighborhoodAnalyzer(df_core=df_core, meta_dict=meta_dict)
    )
    btn_format_neighborhood.disabled = True


btn_format_neighborhood.on_click(format_neighborhood)


@pn.depends(cache_run.param.value)
def reset_neighborhood(run):
    if run:
        lo_neighborhood_analysis.objects = lo_neighborhood_analysis.objects[0:2]
        btn_format_neighborhood.disabled = False
        return btn_format_neighborhood
    else:
        return "Run analysis first"


#### Callback output positioning
#### Neighborhood
################################
lo_neighborhood_analysis.append(reset_neighborhood)

#### Callback output positioning
################################
analysis_tabs.clear()
analysis_tabs.append(("Data overview", update_data_overview))
analysis_tabs.append(("Depth and Coverage", update_quantity))
analysis_tabs.append(("Intramap Scatter", update_cluster_overview))
analysis_tabs.append(("Movement analysis", lo_movement_analysis))
analysis_tabs.append(("Neighborhood analysis", lo_neighborhood_analysis))
analysis_tabs.append(("Download", show_tabular_overview))


# [<div style="text-align: right; font-size: 8pt">back to top</div>](#TOC)
# ## Comparison output tabs<a id="comparison"></a>

# In[ ]:


#### Dashboard structure
#### Overview tab
########################
lo_benchmark_pca = pn.Column()

#### Layout elements
####################

loading_status_comparison = pn.Row()
idle_comparison = pn.indicators.LoadingSpinner(
    value=False, width=100, height=100, color="primary"
)
loading_comparison = pn.indicators.LoadingSpinner(
    value=True, width=100, height=100, color="primary"
)
cache_uploaded_json = pn.widgets.Checkbox(value=False)
cache_run_json = pn.widgets.Checkbox(value=False)
m_diverget_fractions = pn.Pane("")
comparison_status = pn.Pane("No datasets were compared yet")
i_ExpOverview = pn.Row(pn.Pane("", width=1000))

i_multi_choice = pn.widgets.CrossSelector(
    name="Select experiments for comparison",
    value=["a", "b"],
    options=["a", "b", "c"],
    definition_order=False,
)

i_markerset_or_cluster = pn.widgets.Checkbox(
    value=False, name="Display only protein clusters"
)
i_pca_comp_ncomp = pn.widgets.IntSlider(value=3, start=2, end=5)
btn_pca_comp_ncomp = pn.widgets.Button(
    name="Recalculate PCA with different number of components."
)

#### Append layout to dashboard
###############################
for el in [
    pn.Pane(textfragments["benchmark_pca_top"], width=600),
    pn.Row(i_markerset_or_cluster, i_pca_comp_ncomp, btn_pca_comp_ncomp),
]:
    lo_benchmark_pca.append(el)

#### Callbacks
##############
# recalculate_comp_pca
# update_visualization_map_comparison # tab content


def recalculate_comp_pca(event):
    cache_run_json.value = False
    i_class_comp.perform_pca_comparison(n=i_pca_comp_ncomp.value)
    cache_run_json.value = True


btn_pca_comp_ncomp.on_click(recalculate_comp_pca)


@pn.depends(
    i_multi_choice.param.value,
    cache_run_json.param.value,
    i_markerset_or_cluster.param.value,
)
def update_visualization_map_comparison(multi_choice, run_json, markerset_or_cluster):
    try:
        if run_json == True:
            if multi_choice == []:
                return "Please select experiments for comparison"
            else:
                pass

            df_pca = i_class_comp.df_pca.loc[
                i_class_comp.df_pca["Experiment"].isin(multi_choice)
            ].set_index(
                [col for col in i_class_comp.df_pca.columns if not col.startswith("PC")]
            )

            pca_global_comparison = gui.pca_plot(
                df_pca=df_pca,
                df_var=i_class_comp.df_pca_var,
                df_loadings=i_class_comp.df_pca_loadings,
                show_variability=True,
                show_loadings=True,
                color="Cluster" if markerset_or_cluster else "Compartment",
                color_map=i_class_comp.color_maps["Clusters"]
                if markerset_or_cluster
                else i_class_comp.color_maps["Compartments"],
                enable_highlight=False if markerset_or_cluster else True,
                facet_col="Experiment",
                title="PCA plot",
            )
            pca_global_comparison.highlight_dict = i_class_comp.markerproteins
            return pca_global_comparison
        else:
            pca_global_comparison = "Run analysis first!"
            return pca_global_comparison
    except Exception:
        update_status = traceback.format_exc()
        return update_status


#### Callback output positioning
################################
lo_benchmark_pca.append(update_visualization_map_comparison)

#### Dashboard structure
#### Depth tab
########################
lo_benchmark_depth = pn.Column()

#### Layout elements
####################

#### Append layout to dashboard
###############################
for el in [pn.Pane(textfragments["benchmark_depth_top"], width=600)]:
    lo_benchmark_depth.append(el)

#### Callbacks
##############
# update_npr_ngg_nprDc
# update_venn


@pn.depends(i_multi_choice.param.value, cache_run_json.param.value)
def update_npr_ngg_nprDc(multi_choice, run_json):
    try:
        if run_json == True:
            if multi_choice == []:
                return pn.Column(pn.Row("Please select experiments for comparison"))
            else:
                fig_quantity_pg, fig_quantity_pr = (
                    i_class_comp.quantity_pr_pg_barplot_comparison(
                        multi_choice=multi_choice
                    )
                )
                coverage_barplot = i_class_comp.coverage_comparison(
                    multi_choice=multi_choice
                )
                return pn.Row(
                    pn.Pane(fig_quantity_pg, config=plotly_config),
                    pn.Pane(coverage_barplot, config=plotly_config),
                )
        else:
            completeness_barplot = "Run analysis first!"
            return completeness_barplot
    except Exception:
        update_status = traceback.format_exc()
        return update_status


@pn.depends(i_multi_choice.param.value, cache_run_json.param.value)
def update_venn(multi_choice, run_json):
    try:
        if run_json == True:
            venn_plot = []
            if len(multi_choice) <= 1:
                return pn.Column(
                    pn.Row(
                        pn.Pane("Please select 2 or more experiments for comparison"),
                        width=1000,
                    )
                )
            else:
                (
                    venn_plot_total,
                    venn_plot_int,
                    figure_UpSetPlot_total,
                    figure_UpSetPlot_int,
                ) = i_class_comp.venn_sections(multi_choice_venn=multi_choice)
                return pn.Row(
                    pn.Column(
                        "Proteins quantified in at least one map",
                        pn.Pane(venn_plot_total),
                        pn.Row(figure_UpSetPlot_total, width=1000),
                    ),
                    pn.Column(
                        "Proteins quantified in all maps",
                        pn.Pane(venn_plot_int),
                        pn.Row(figure_UpSetPlot_int, width=1000),
                    ),
                )
        else:
            venn_plot = "Run analysis first!"
            return venn_plot
    except Exception:
        update_status = traceback.format_exc()
        return update_status


#### Callback output positioning
################################
lo_benchmark_depth.append(update_npr_ngg_nprDc)
lo_benchmark_depth.append(update_venn)

#### Dashboard structure
#### Intramap scatter tab
########################

#### Layout elements
####################
i_complexes_norm = pn.widgets.Select(
    name="Normalize complexes to", options=["Median across all experiments"]
)
i_complexes_plot = pn.widgets.Select(
    name="Plot type", options=["stacked", "strip", "box", "violin", "histogram"]
)
i_complexes_sync = pn.widgets.Checkbox(
    name="Highlight selected cluster (summary plots only)"
)
i_complexes_detail = pn.widgets.Select(
    options=[], name="Cluster of interest", width=300
)

lo_complexes = pn.WidgetBox(
    "**Compare experiments**", i_complexes_norm, i_complexes_plot, i_complexes_sync
)

i_clusters_for_ranking = pn.widgets.CrossSelector(
    name="Select clusters to be considered for ranking calculation", options=[], size=8
)
i_minn_proteins = pn.widgets.IntSlider(
    name="Minimum number of proteins per complex", start=3, end=13, step=1, value=5
)
i_collapse_maps = pn.widgets.Checkbox(value=False, name="Collapse maps")

#### Append layout to dashboard
###############################

#### Callbacks
##############
# update_comp_cluster_coverage
# update_comp_bp_global
# update_comp_bp_single


@pn.depends(
    i_multi_choice.param.options,
    i_minn_proteins.param.value,
    cache_run_json.param.value,
)
def update_comp_cluster_coverage(exp_names, min_n, run_json):
    try:
        if not run_json:
            return ""
        [f, p, n] = i_class_comp.get_complex_coverage(min_n)
        i_clusters_for_ranking.options = [
            el for el in i_class_comp.markerproteins.keys() if el not in n.keys()
        ]
        i_complexes_detail.options = [
            el for el in i_class_comp.markerproteins.keys() if el not in n.keys()
        ]
        i_clusters_for_ranking.value = [
            el for el in i_class_comp.markerproteins.keys() if el in f.keys()
        ]
        return pn.Row(
            "Coverage in all experiments \[>= n proteins]:<br>"
            + "<br>".join(["- {} ({})".format(k, v) for k, v in f.items()]),
            "Coverage in some experiments \[proteins/experiment]:<br>"
            + "<br>".join(["- {} \{}".format(k, str(v)) for k, v in p.items()]),
            "No sufficient coverage in any experiment \[proteins/experiment]:<br>"
            + "<br>".join(["- {} \{}".format(k, str(v)) for k, v in n.items()]),
        )
    except Exception:
        update_status = traceback.format_exc()
        return update_status


@pn.depends(
    i_multi_choice.param.value,
    i_clusters_for_ranking.param.value,
    i_minn_proteins.param.value,
    i_complexes_norm.param.value,
    i_complexes_plot.param.value,
    i_complexes_sync.param.value,
    i_complexes_detail.param.value,
    cache_run_json.param.value,
)
def update_comp_bp_global(
    multi_choice,
    clusters_for_ranking,
    min_n,
    norm,
    plot_type,
    sync,
    highlight,
    run_json,
):
    try:
        if not run_json:
            return ""
        if set(multi_choice) != set(i_class_comp.df_distance_comp.Experiment.values):
            i_class_comp.calc_biological_precision(multi_choice)
            i_complexes_norm.options = [i_complexes_norm.options[0]] + multi_choice
        if clusters_for_ranking == []:
            return pn.Row(lo_complexes, "Select at least one cluster")
        else:
            medians, plot = i_class_comp.plot_intramap_scatter(
                normalization=np.median
                if norm == "Median across all experiments"
                else norm,
                aggregate_proteins=True,
                aggregate_maps=False,
                plot_type=plot_type,
                highlight=None if sync == False else highlight,
                min_size=min_n,
                multi_choice=multi_choice,
                clusters_for_ranking=clusters_for_ranking,
            )
            return pn.Row(
                lo_complexes,
                medians.rename({"distance": "median distance"}, axis=1),
                pn.Pane(plot, config=plotly_config),
            )

    except Exception:
        update_status = traceback.format_exc()
        return pn.Row(lo_complexes, update_status)


@pn.depends(
    i_multi_choice.param.value,
    i_complexes_detail.param.value,
    i_collapse_maps.param.value,
    cache_run_json.param.value,
)
def update_comp_bp_single(
    multi_choice, clusterwidget_comparison, collapse_maps, run_json
):
    try:
        i_class_comp.cache_cluster_quantified = True
        distance_comparison = i_class_comp.plot_intramap_scatter_cluster(
            collapse_maps=collapse_maps,
            cluster_of_interest_comparison=clusterwidget_comparison,
            multi_choice=multi_choice,
        )
        if i_class_comp.cache_cluster_quantified == False:
            return "Cluster was not quantified in any experiment"
        else:
            df_pca = i_class_comp.df_cluster_pca.xs(
                clusterwidget_comparison, level="Cluster", axis=0
            )
            df_pca = df_pca.query(f"Experiment in {str(multi_choice)}")
            if collapse_maps:
                df_pca = df_pca.groupby(
                    [el for el in df_pca.index.names if "Map" not in el]
                ).mean()
            pca_comparison = gui.pca_plot(
                df_pca=df_pca,
                color="Experiment",
                color_map=dict(),
                facet_col=None,
                enable_highlight=False,
                show_variability=False,
                show_loadings=False,
                title="PCA plot for {}".format(clusterwidget_comparison),
            )
            return pn.Row(
                pca_comparison, pn.Pane(distance_comparison, config=plotly_config)
            )
    except Exception:
        update_status = traceback.format_exc()
        return update_status


#### Callback output positioning
################################
comparison_tab_bp = pn.Column(
    pn.Pane(textfragments["benchmark_intra_top"], width=600),
    pn.Row(
        pn.Column(i_clusters_for_ranking, i_minn_proteins), update_comp_cluster_coverage
    ),
    update_comp_bp_global,
    pn.Row(i_complexes_detail, i_collapse_maps),
    update_comp_bp_single,
)

#### Dashboard structure
#### SMV analysis
########################
lo_benchmark_SVMs_tabs = pn.Tabs()
lo_benchmark_SVMs = pn.Column(
    pn.Pane(textfragments["benchmark_SVM_top"], width=600),
    lo_benchmark_SVMs_tabs,
    pn.Row(name="svm_output"),
)
lo_benchmark_SVMs_addmcm = pn.Column(name="add_mcmatrix")
lo_benchmark_SVMs_runsvm = pn.Column(name="run_svm")
lo_benchmark_SVMs_tabs.append(("Upload misclassification", lo_benchmark_SVMs_addmcm))
lo_benchmark_SVMs_tabs.append(("Run SVMs", lo_benchmark_SVMs_runsvm))

#### Layout elements
#### SVM analysis
####################
SVM_status = pn.pane.Markdown(width=400)
lo_SVM_heatmap = pn.Column(SVM_status)
i_SVMmatrix = pn.widgets.input.TextAreaInput(
    name="Misclassification matrix", placeholder="Copy matrix here..."
)
i_SVMsource = pn.widgets.Select(
    options=["Perseus", "MetaMass", "direct"],
    value="Perseus",
    name="Select source of misclassification matrix",
)
i_SVMname = pn.widgets.TextInput(
    name="Set name",
    value="default",
    placeholder="Identify the set of misclassification matrices (e.g. rbf_SVM)",
)
lo_SVM_source = pn.Row(
    i_SVMsource,
    gui.help_icon(
        "For Perseus, copy whole matrix via Ctrl-A, Ctrl-C. For MetaMass copy matrix including column headings, but without row labels. Direct assumes true classes in rows and predicted classes in columns, with column headings only."
    ),
)
i_SVMcomment = pn.widgets.TextInput(name="Comment", placeholder="Add comments here...")
i_SVMexp = pn.widgets.Select(
    name="Select experiments for the assignment of a misclassification matrix",
    options=["a", "b", "c"],
)
btn_SVM_addmatrix = pn.widgets.Button(
    name="Update misclassification matrix",
    button_type="success",
    width=400,
    height=50,
    css_classes=["button-main"],
)
i_svm_set = pn.widgets.Select(options=["default"], name="Select set")
i_svm_score = pn.widgets.Select(
    options=["F1 score", "Precision", "Recall", "Class size"], value="F1 score"
)

#### Append layout to dashboard
#### SVM analysis
###############################
for el in [
    i_SVMexp,
    lo_SVM_source,
    i_SVMname,
    i_SVMcomment,
    i_SVMmatrix,
    lo_SVM_heatmap,
    btn_SVM_addmatrix,
]:
    lo_benchmark_SVMs_addmcm.append(el)

#### Callbacks
# update_SVMexp
# add_SVM_result
# update_SVM_Analysis # needs update, including the backend function in domaps.py


@pn.depends(i_multi_choice.param.value, watch=True)
def update_SVMexp(multi_choice):
    i_SVMexp.options = multi_choice


def add_SVM_result(event):
    """
    Adds SVM matrix uploaded in the tool to memory and to the current analysis
    """
    # 1. Get input from interface

    experiment, SVMname, SVMsource, SVMcomment = (
        i_SVMexp.value,
        i_SVMname.value,
        i_SVMsource.value,
        i_SVMcomment.value,
    )

    # change comment of a stored dataset
    try:
        if (
            i_SVMmatrix.value == ""
            and i_class_comp.svm_results[experiment][SVMname]["misclassification"]
            is not None
        ):
            SVMmatrix = i_class_comp.svm_results[experiment][SVMname][
                "misclassification"
            ]
        # return error if no misclassification  is uploaded or no comment is changed
        # if i_SVMmatrix.value == "":
        #    SVM_status.object = "No misclassification matrix is uploaded"
        else:
            SVMmatrix = pd.read_table(StringIO(i_SVMmatrix.value), sep="\t")
    except Exception:
        update_status = traceback.format_exc()
        SVM_status.object = update_status
    #
    # 2. Add to i_class_comp
    # defaults to name="default" and overwrite=True
    i_class_comp.add_svm_result(
        experiment, SVMmatrix, name=SVMname, source=SVMsource, comment=SVMcomment
    )
    if SVMname not in i_svm_set.options:
        i_svm_set.options = i_svm_set.options + [SVMname]
    # 3. Add to mem_available_datasets so it can be downloaded together with the data
    mem_available_datasets[experiment]["SVM results"] = copy.deepcopy(
        i_class_comp.svm_results[experiment]
    )
    for k in i_class_comp.svm_results[experiment].keys():
        mc_json = mem_available_datasets[experiment]["SVM results"][k][
            "misclassification"
        ].to_json()
        mem_available_datasets[experiment]["SVM results"][k]["misclassification"] = (
            mc_json
        )
        pr_json = mem_available_datasets[experiment]["SVM results"][k][
            "prediction"
        ].to_json()
        mem_available_datasets[experiment]["SVM results"][k]["prediction"] = pr_json

    # empty misclassification matrix of the textarea-input widget
    i_SVMmatrix.value = ""
    # 4. Trigger update of figure
    update_SVM_Analysis(i_multi_choice.value)


btn_SVM_addmatrix.on_click(add_SVM_result)


@pn.depends(
    i_SVMexp.param.value, i_SVMname.param.value, i_SVMmatrix.param.value, watch=True
)
def fill_svm_comment(SVMexp, SVMname, SVMmatrix):
    """
    Acess stored comment if possible and no new matrix was uploaded.
    """
    try:
        comment = i_class_comp.svm_results[SVMexp][SVMname]["comment"]
        if SVMmatrix != "" or comment == "":
            timestampStr = datetime.now().strftime("%d-%b-%Y (%H:%M:%S)")
            comment = "Matrix added to set {} on {}".format(SVMname, timestampStr)
        else:
            pass
    except Exception:
        timestampStr = datetime.now().strftime("%d-%b-%Y (%H:%M:%S)")
        comment = "Matrix added to set {} on {}".format(SVMname, timestampStr)
    i_SVMcomment.value = comment


@pn.depends(i_SVMexp.param.value, i_SVMname.param.value, i_SVMmatrix.param.value)
def show_misclassification(SVMexp, SVMname, SVMmatrix):
    """
    Display heatmap of freshly uploaded or reloaded SVM misclassification matrix
    """
    if SVMmatrix == "":
        try:
            df_SVM = i_class_comp.svm_results[SVMexp][SVMname]["misclassification"]
            SVMheatmap = domaps.svm_heatmap(df_SVM)
        except:
            SVMheatmap = "No misclassification matrix is uploaded"
    else:
        df_SVM = pd.read_table(StringIO(SVMmatrix), sep="\t")
        SVMheatmap = domaps.svm_heatmap(df_SVM)
    return SVMheatmap


@pn.depends(i_SVMexp.param.value, i_SVMname.param.value)
def load_matrix(SVMexp, SVMname):
    """
    Load data if new experiment is selected
    """
    try:
        SVMmatrix = i_class_comp.svm_results[SVMexp][SVMname]["misclassification"]
        SVM_status.object = "Misclassification Matrix from class"
    except:
        SVM_status.object = "Upload missclassificationmatrix first"


@pn.depends(i_multi_choice.param.value, i_svm_set.param.value, i_svm_score.param.value)
def show_svm_results(multi_choice, svmset, score):
    try:
        if len(i_class_comp.svm_results) == 0:
            return "Please upload SVM results first"
        plot_detail = i_class_comp.plot_svm_detail(
            multi_choice=multi_choice, score=score, svmset=svmset
        )
        plot_summary = i_class_comp.plot_svm_summary(
            multi_choice=multi_choice, score=score, svmset=svmset
        )
        return pn.Column(
            pn.Pane(plot_summary, config=plotly_config),
            pn.Pane(plot_detail, config=plotly_config),
        )
    except:
        return traceback.format_exc()


@pn.depends(i_multi_choice.param.value, watch=True)
def update_SVM_Analysis(multi_choice):
    # empty lo if available
    lo_benchmark_SVMs.objects[
        [i.name for i in lo_benchmark_SVMs].index("svm_output")
    ].objects = []
    try:
        if multi_choice == []:
            lo = pn.Column(pn.Row("Please select experiments for comparison"))
        else:
            lo = pn.Column(i_svm_set, i_svm_score, show_svm_results)
    except Exception:
        update_status = traceback.format_exc()
        lo = update_status
    lo_benchmark_SVMs.objects[
        [el.name for el in lo_benchmark_SVMs.objects].index("svm_output")
    ].append(lo)


#### Callback output positioning
#### SVM analysis
################################
lo_SVM_heatmap.append(show_misclassification)
# lo_benchmark_SVMs.objects[[el.name for el in lo_benchmark_SVMs.objects].index("svm_output")].append(update_SVM_Analysis())


i_svm_svmcomp = domaps.SVMComp("exp1;exp2__42__0.2__cl1;cl2")
i_svm_randomstate = pn.widgets.IntInput(
    value=42, width=150, name="Seed for randomization"
)
i_svm_testsplit = pn.widgets.FloatInput(
    start=0, end=1, step=0.05, value=0.2, width=150, name="Test set proportion"
)
i_svm_classes = gui.ConfigureSVMClasses()
i_svm_Crange = pn.widgets.IntRangeSlider(
    start=1, end=50, step=1, value=(1, 30), name="C (% misclassification)", width=250
)
i_svm_gammarange = pn.widgets.IntRangeSlider(
    start=1, end=100, step=1, value=(1, 50), name="gamma (RBF diameter)", width=250
)
i_svm_trainingrounds = pn.widgets.IntSlider(
    start=3, end=10, value=5, name="Training iterations", width=150
)
btn_svm_train = pn.widgets.Button(name="Run training", width=150, button_type="success")
i_svm_canvas = pn.Column(pn.Row())
i_svm_C = pn.widgets.FloatInput(name="C", background="salmon", width=150)
i_svm_gamma = pn.widgets.FloatInput(name="gamma", background="salmon", width=150)
i_svm_nameprediction = pn.widgets.TextInput(
    name="Name prediction set", width=150, value="SVM run"
)
i_svm_minp = pn.widgets.FloatInput(
    start=0.1, end=1.0, step=0.05, value=0.4, width=150, name="Minimum SVM probability"
)
i_svm_minpdiff = pn.widgets.FloatInput(
    start=0.00,
    end=0.9,
    step=0.05,
    value=0.15,
    width=150,
    name="Minimum difference to second",
)
btn_svm_predict = pn.widgets.Button(name="Run predictions", width=150)
i_svm_canvas_prediction = pn.Column()

for el in [
    pn.panel(textfragments["benchmark_SVM_internal"], width=600),
    "**Define classes and training/test split for SVMs**",
    pn.Row(
        i_svm_classes,
        pn.Column(i_svm_randomstate, i_svm_testsplit),
        gui.help_icon(
            "The seed for randomization is used for defining the test-split, as well as the folds for crossvalidation. By setting this the analysis can be exactly repeated at any time."
        ),
    ),
    "**Configure training for hyperparameter optimization**",
    pn.Row(
        i_svm_Crange,
        i_svm_gammarange,
        i_svm_trainingrounds,
        btn_svm_train,
        gui.help_icon(
            "The hyper parameter optimization is done in an iterative grid search. These paramters define where the search starts and how many iterations it should run."
        ),
    ),
    pn.Card(
        i_svm_canvas,
        header=pn.pane.Markdown("**Training output**", width=860),
        width=860,
    ),
    "**Select hyper parameters and run predictions**",
    pn.Row(
        i_svm_C,
        i_svm_gamma,
    ),
    pn.Row(
        i_svm_minp,
        i_svm_minpdiff,
        i_svm_nameprediction,
        btn_svm_predict,
        gui.help_icon(
            "The minimum probability is used to determine whether any organelle matches the protein well. The minimum difference determines the labelling if more than one class matches the minimum probability."
        ),
    ),
    i_svm_canvas_prediction,
]:
    lo_benchmark_SVMs_runsvm.append(el)


@pn.depends(i_multi_choice.param.value, watch=True)
def update_svm_class(multi_choice):
    try:
        i_svm_svmcomp.experiments = multi_choice
        available_classes = list(
            set(
                i_class_comp.df_01_filtered_combined.drop(
                    "undefined", axis=0, level="Compartment"
                ).index.get_level_values("Compartment")
            )
        )
        i_svm_svmcomp.classes = available_classes
        i_svm_svmcomp.set_df(i_class_comp.df_01_filtered_combined)
        class_counts = i_svm_svmcomp._get_markers("shared").value_counts("Compartment")
        i_svm_classes.class_counts = class_counts
        i_svm_classes.classes = list(class_counts.index)
        train, test = i_svm_svmcomp._train_test_split(
            test_percent=i_svm_testsplit.value, random_state=i_svm_randomstate.value
        )
        i_svm_classes.train = train
        i_svm_classes.test = test
        update_canvas_hash()
    except:
        lo_benchmark_SVMs_runsvm.append(traceback.format_exc())


@pn.depends(
    i_svm_randomstate.param.value,
    i_svm_testsplit.param.value,
    i_svm_classes.param.classes,
    watch=True,
)
def update_test_split(random_state, test_split, classes):
    try:
        i_svm_svmcomp.random_state = random_state
        i_svm_svmcomp.test_split = test_split
        if len(classes) > 0:
            i_svm_svmcomp.classes = classes
        train, test = i_svm_svmcomp._train_test_split(
            test_percent=test_split, random_state=random_state
        )
        i_svm_classes.train = train
        i_svm_classes.test = test
        update_canvas_hash()
    except:
        lo_benchmark_SVMs_runsvm.append(traceback.format_exc())


def update_canvas_hash():
    hash_data = domaps.SVMComp._construct_hash_data(
        experiments=i_svm_svmcomp.experiments,
        random_state=i_svm_svmcomp.random_state,
        test_split=i_svm_svmcomp.test_split,
        classes=i_svm_svmcomp.classes,
    )
    if hash_data in i_class_comp.svm_runs.keys():
        i_svm_canvas.objects = [
            pn.Row(
                domaps.SVMComp._plot_training(
                    i_class_comp.svm_runs[hash_data].accuracies
                )
            )
        ]
        i_svm_C.value = i_class_comp.svm_runs[hash_data].C
        i_svm_gamma.value = i_class_comp.svm_runs[hash_data].gamma
        sio = StringIO()
        i_class_comp.svm_runs[hash_data].accuracies.to_csv(sio)
        sio.seek(0)
        download_acc = pn.widgets.FileDownload(
            file=sio, filename="SVMtraining_accuracies.csv"
        )
        i_svm_canvas.append(download_acc)
    else:
        i_svm_canvas.objects = [pn.Row()]
        btn_svm_train.button_type = "success"


def svm_run_training(event):
    try:
        hash_data = domaps.SVMComp._construct_hash_data(
            experiments=i_multi_choice.value,
            random_state=i_svm_randomstate.value,
            test_split=i_svm_testsplit.value,
            classes=i_svm_classes.classes,
        )
        if (
            hash_data not in i_class_comp.svm_runs.keys()
            or btn_svm_train.button_type != "success"
        ):
            i_svm_canvas.objects = ["Starting training"]
        C, gamma = i_class_comp.train_svm(
            experiments=i_multi_choice.value,
            random_state=i_svm_randomstate.value,
            test_split=i_svm_testsplit.value,
            classes=i_svm_classes.classes,
            output="canvas",
            canvas=i_svm_canvas,
            rounds=i_svm_trainingrounds.value,
            C0=i_svm_Crange.value[0],
            C1=i_svm_Crange.value[1],
            g0=i_svm_gammarange.value[0],
            g1=i_svm_gammarange.value[1],
            overwrite=False if btn_svm_train.button_type == "success" else True,
        )
        btn_svm_train.button_type = "success"
        sio = StringIO()
        i_class_comp.svm_runs[hash_data].accuracies.to_csv(sio)
        sio.seek(0)
        download_acc = pn.widgets.FileDownload(
            file=sio, filename="SVMtraining_accuracies.csv"
        )
        i_svm_canvas.append(download_acc)
        i_svm_C.value = C
        i_svm_gamma.value = gamma
    except RuntimeError:
        btn_svm_train.button_type = "danger"
        i_svm_canvas[0] = "Please confirm that you want to retrain"
        pass
    except:
        i_svm_canvas[0] = traceback.format_exc()
        pass


btn_svm_train.on_click(svm_run_training)


def svm_run_prediction(event):
    try:
        i_svm_canvas_prediction.objects = ["Running predictions"]
        prediction = i_class_comp.predict_svm(
            experiments=i_multi_choice.value,
            random_state=i_svm_randomstate.value,
            test_split=i_svm_testsplit.value,
            classes=i_svm_classes.classes,
            C=i_svm_C.value,
            gamma=i_svm_gamma.value,
            min_p=i_svm_minp.value,
            min_diff=i_svm_minpdiff.value,
            svmset=i_svm_nameprediction.value,
        )
        hash_data = domaps.SVMComp._construct_hash_data(
            experiments=i_multi_choice.value,
            random_state=i_svm_randomstate.value,
            test_split=i_svm_testsplit.value,
            classes=i_svm_classes.classes,
        )
        prob = i_class_comp.svm_runs[hash_data].probabilities
        out = (
            i_class_comp.id_alignment.drop(i_class_comp.id_alignment.columns, axis=1)
            .join(prediction)
            .join(prob)
        )
        sio = StringIO()
        out.to_csv(sio)
        sio.seek(0)
        download_pred = pn.widgets.FileDownload(file=sio, filename="SVMprediction.csv")
        i_svm_canvas_prediction.append(download_pred)

        for exp in i_multi_choice.value:
            mem_available_datasets[exp]["SVM results"] = copy.deepcopy(
                i_class_comp.svm_results[exp]
            )
            for k in i_class_comp.svm_results[exp].keys():
                mc_json = mem_available_datasets[exp]["SVM results"][k][
                    "misclassification"
                ].to_json()
                mem_available_datasets[exp]["SVM results"][k]["misclassification"] = (
                    mc_json
                )
                pr_json = mem_available_datasets[exp]["SVM results"][k][
                    "prediction"
                ].to_json()
                mem_available_datasets[exp]["SVM results"][k]["prediction"] = pr_json

        i_svm_set.options = i_svm_set.options + [i_svm_nameprediction.value]
        i_svm_set.value = i_svm_nameprediction.value
    except:
        i_svm_canvas_prediction.objects = [traceback.format_exc()]


btn_svm_predict.on_click(svm_run_prediction)


#### Dashboard structure
#### Intermap scatter tab
########################
lo_benchmark_intermap = pn.Column(
    pn.Row(name="intermap_top"),
    pn.Row(name="intermap_scatter"),
    pn.layout.VSpacer(background="#AAAAAA", height=3),
    pn.Row(name="correlation_plots"),
)

#### Layout elements
####################

i_scatter_metric = pn.widgets.Select(
    name="Distance metric",
    options=[
        "manhattan distance to average profile",
        "manhattan distance to median profile",
        "euclidean distance",
        "manhattan distance",
        "1 - cosine correlation",
        "1 - pearson correlation",
    ],
)
i_scatter_consolidation = pn.widgets.Select(
    name="Consolidation of replicate distances", options=["average", "median", "sum"]
)
i_scatter_quantile = pn.widgets.FloatInput(name="Highlight quantile", value=0.5)
i_scatter_showrug = pn.widgets.Checkbox(
    name="Show rugplot in bottom margin", value=False
)
i_scatter_showfull = pn.widgets.Checkbox(
    name="Show additional traces for full datasets", value=False
)
i_scatter_x_cut = pn.widgets.FloatInput(name="Cut x-axis at", value=1)

i_comp_corr_mode = pn.widgets.Select(
    options=["heatmap", "scatter"], name="Show correlations as"
)
i_comp_corr_measure = pn.widgets.Select(
    options=["Pearson", "Spearman"], name="Correlation metric"
)
i_comp_corr_fraction = pn.widgets.Select(
    options=["all fractions"], name="Select fraction(s)"
)
i_comp_corr_experiment = pn.widgets.Select(
    options=["all experiments"], name="Select experiment(s)"
)

#### Append layout to dashboard
###############################
lo_benchmark_intermap.objects[0].append(
    pn.Pane(textfragments["benchmark_inter_top"], width=600)
)
lo_benchmark_inter_widgets = pn.WidgetBox(
    i_scatter_metric,
    i_scatter_consolidation,
    i_scatter_showrug,
    i_scatter_showfull,
    i_scatter_quantile,
    i_scatter_x_cut,
)
lo_benchmark_intermap.objects[1].append(lo_benchmark_inter_widgets)
lo_benchmark_corr_widgets = pn.WidgetBox(
    i_comp_corr_mode, i_comp_corr_measure, i_comp_corr_fraction, i_comp_corr_experiment
)
lo_benchmark_intermap.objects[3].append(lo_benchmark_corr_widgets)

#### Callbacks
##############
# update_corr_selection_comp
# show_correlation_comp
# update_global_scatter_comparison


@param.depends(cache_run_json.param.value, i_multi_choice.param.value, watch=True)
def update_corr_selection_comp(run, multi_choice):
    if run == True:
        i_comp_corr_fraction.options = [
            i_comp_corr_fraction.options[0]
        ] + i_class_comp.fractions
        i_comp_corr_experiment.options = [
            i_comp_corr_experiment.options[0]
        ] + multi_choice


@param.depends(
    i_comp_corr_mode.param.value,
    i_comp_corr_measure.param.value,
    i_comp_corr_fraction.param.value,
    i_comp_corr_experiment.param.value,
    i_multi_choice.param.value,
    cache_run_json.param.value,
)
def show_correlation_comp(mode, measure, fractions, experiments, multi_choice, run):
    if run == True:
        try:
            if mode == "scatter" and (
                fractions == "all fractions" or experiments == "all experiments"
            ):
                return "Please select only one fraction and experiment in scatter mode."
            df = i_class_comp.df_01_filtered_combined.loc[
                i_class_comp.df_01_filtered_combined.index.get_level_values(
                    "Experiment"
                ).isin(multi_choice),
                :,
            ].copy()
            df = (
                df.reset_index("Exp_Map", drop=True)
                .unstack(["Experiment", "Map"])
                .sort_index(axis=1, level="Fraction")
            )

            if fractions != "all fractions":
                df = df.xs(fractions, level="Fraction", axis=1, drop_level=False)
            if experiments != "all experiments":
                df = df.xs(experiments, level="Experiment", axis=1, drop_level=False)

            df.columns = ["_".join(el) for el in df.columns]

            measure_dict = {
                "Spearman": lambda x: spearmanr(x.values).correlation,
                "Pearson": lambda x: np.corrcoef(x.T),
            }

            plot = domaps.plot_sample_correlations(
                df,
                data_columns="(.*)",
                log10=False,
                binning=50,
                mode=mode,
                correlation_function=measure_dict[measure],
            )
            return plot

        except:
            return traceback.format_exc()
    else:
        return "Run analysis first"


@pn.depends(
    i_multi_choice.param.value,
    i_scatter_metric.param.value,
    i_scatter_consolidation.param.value,
    cache_run_json.param.value,
    i_scatter_showfull.param.value,
    i_scatter_quantile.param.value,
    i_scatter_showrug.param.value,
    i_scatter_x_cut.param.value,
)
def update_global_scatter_comparison(
    multi_choice, metric, consolidation, run_json, show_full, quantile, show_rug, x_cut
):
    try:
        if run_json == False:
            return "Run analysis first!"
        if multi_choice == []:
            return pn.Column(pn.Row("Please select experiments for comparison"))
        if (
            "reproducibility consolidation" not in i_class_comp.parameters.keys()
            or i_class_comp.parameters["reproducibility metric"] != metric
            or i_class_comp.parameters["reproducibility consolidation"] != consolidation
        ):
            i_class_comp.calculate_global_scatter(metric, consolidation)
        reproducibility_plot = i_class_comp.plot_reproducibility_distribution(
            multi_choice=multi_choice,
            x_cut=x_cut,
            q=quantile,
            show_rug=show_rug,
            show_full=show_full,
        )
        return pn.Pane(reproducibility_plot, config=plotly_config)
    except Exception:
        update_status = traceback.format_exc()
        return update_status


#### Callback output positioning
################################
lo_benchmark_intermap[1].append(update_global_scatter_comparison)
lo_benchmark_intermap[3].append(show_correlation_comp)


#### Dashboard structure
#### Profile comparison tab
########################

#### Layout elements
####################

i_compare_gene = pn.widgets.TextInput(
    value="PLEC", name="Enter gene name or protein ID to see profile."
)
i_compare_profile_style = pn.widgets.Select(
    options=[
        # "all profiles",
        "mean +- stdev",
        "mean +- SEM",
    ]
)
i_compare_compartment = pn.widgets.MultiSelect(
    options=[], name="Select compartments for which to show summary profiles."
)

#### Append layout to dashboard
###############################

#### Callbacks
##############
# update_profile_comparison


@pn.depends(
    i_multi_choice.param.value,
    i_compare_gene.param.value,
    i_compare_compartment.param.value,
    i_compare_profile_style.param.value,
    cache_run_json.param.value,
)
def update_profile_comparison(
    multi_choice, compare_gene, compare_compartment, compare_profile_style, run_json
):
    if not run_json:
        return "Run analysis first!"
    try:
        try:
            plotdata = (
                i_class_comp.df_01_filtered_combined.xs(
                    compare_gene, level="Gene names", axis=0, drop_level=False
                )
                .stack("Fraction")
                .reset_index()
                .rename({0: "Profile [% total signal]"}, axis=1)
            )
        except:
            try:
                plotdata = (
                    i_class_comp.df_01_filtered_combined.xs(
                        compare_gene, level="Protein IDs", axis=0, drop_level=False
                    )
                    .stack("Fraction")
                    .reset_index()
                    .rename({0: "Profile [% total signal]"}, axis=1)
                )
            except:
                try:
                    plotdata = (
                        i_class_comp.df_01_filtered_combined.loc[
                            [
                                compare_gene in " ".join(el)
                                for el in i_class_comp.df_01_filtered_combined.index.values
                            ],
                            :,
                        ]
                        .stack("Fraction")
                        .reset_index()
                        .rename({0: "Profile [% total signal]"}, axis=1)
                    )
                except:
                    plotdata = pd.DataFrame()
        if len(plotdata) > 0:
            plotdata.drop("Exp_Map", axis=1, inplace=True)
            plotdata.sort_values("Fraction", key=domaps.natsort_list_keys, inplace=True)
            experiments = [
                el for el in multi_choice if el in plotdata["Experiment"].values
            ]
            plotdata = (
                plotdata.set_index("Experiment").loc[experiments, :].reset_index()
            )
            plotprofile = px.line(
                plotdata,
                x="Fraction",
                y="Profile [% total signal]",
                line_group="Map",
                facet_col="Protein IDs",
                template="simple_white",
                line_dash="Sequence" if "Sequence" in plotdata.columns else None,
                color="Experiment",
                hover_data=list(plotdata.columns),
            )
        else:
            plotprofile = "No gene or protein ID matching {} found.".format(
                compare_gene
            )

        plotdata = pd.DataFrame()
        if compare_profile_style == "all profiles":
            for el in compare_compartment:
                el_df = (
                    i_class_comp.df_01_filtered_combined.xs(
                        el, level="Compartment", axis=0, drop_level=False
                    )
                    .stack("Fraction")
                    .reset_index()
                    .rename({0: "Profile [% total signal]"}, axis=1)
                )
                plotdata = plotdata.append(el_df)
            if len(plotdata) > 0:
                plotdata.sort_values(
                    "Fraction", key=domaps.natsort_list_keys, inplace=True
                )
                plotdata = (
                    plotdata.set_index("Experiment").loc[multi_choice, :].reset_index()
                )
                plotdata.insert(
                    0,
                    "PG_Map",
                    [
                        str(p) + "_" + str(m)
                        for p, m in zip(plotdata["Protein IDs"], plotdata["Map"])
                    ],
                )
                plotcompartments = px.box(
                    plotdata,
                    x="Fraction",
                    y="Profile [% total signal]",
                    color="Experiment",
                    facet_row="Compartment",
                    template="simple_white",
                )
            else:
                plotcompartments = "Please select at least one compartment"
        else:
            for el in compare_compartment:
                el_df = i_class_comp.df_01_filtered_combined.xs(
                    el, level="Compartment", axis=0, drop_level=False
                )
                plotdata = plotdata.append(
                    el_df.stack("Fraction")
                    .groupby(["Compartment", "Map", "Experiment", "Fraction"])
                    .apply(
                        lambda x: pd.Series(
                            {
                                "Profile [% total signal]": np.nanmean(x),
                                "std": np.nanstd(x),
                                "sem": np.nanstd(x) / np.sqrt(sum(np.isfinite(x))),
                            }
                        )
                    )
                    .reset_index()
                    .rename(columns={"level_4": "measure", 0: "value"})
                    .set_index(
                        ["Compartment", "Map", "Experiment", "Fraction", "measure"]
                    )
                    .unstack("measure")
                    .droplevel(0, axis=1)
                    .reset_index()
                )
            if len(plotdata) > 0:
                plotdata.sort_values(
                    "Fraction", key=domaps.natsort_list_keys, inplace=True
                )
                plotdata = (
                    plotdata.set_index("Experiment").loc[multi_choice, :].reset_index()
                )
                plotcompartments = px.line(
                    plotdata,
                    x="Fraction",
                    y="Profile [% total signal]",
                    color="Experiment",
                    line_group="Map",
                    line_dash="Compartment",
                    template="simple_white",
                    error_y="std" if "stdev" in compare_profile_style else "sem",
                )
            else:
                plotcompartments = "Please select at least one compartment"

        return pn.Column(
            pn.Pane(textfragments["benchmark_profile_top"], width=600),
            pn.Row(
                pn.Column(i_compare_gene, pn.Pane(plotprofile, config=plotly_config)),
                pn.Column(
                    pn.Row(i_compare_compartment, i_compare_profile_style),
                    pn.Pane(plotcompartments, config=plotly_config),
                ),
            ),
        )
    except Exception:
        return pn.Column(
            pn.Pane(textfragments["benchmark_profile_top"], width=600),
            pn.Row(
                pn.Column(i_compare_gene, traceback.format_exc()),
                pn.Column(pn.Row(i_compare_compartment, i_compare_profile_style)),
            ),
        )


#### Callback output positioning
################################


#### Dashboard structure
#### Download data
########################
lo_benchmark_download = pn.Column()

#### Layout elements
#### Download data
####################
i_benchmark_downloadselector = pn.widgets.Select(
    options=[
        "0-1 normalized data.csv",
        "pca coordinates.csv",
        "complex scatter.csv",
        "reproducibility.csv",
        "protein id alignment.csv",
        "benchmarking results collection.xlsx",
    ],
    value="0-1 normalized data.csv",
    width=200,
)
i_benchmark_download = pn.widgets.FileDownload(
    label="Download file", width=200, button_type="success"
)
o_benchmark_downloadpreview = pn.Row()

#### Append layout to dashboard
#### Download data
###############################
for el in [
    i_benchmark_downloadselector,
    i_benchmark_download,
    o_benchmark_downloadpreview,
]:
    lo_benchmark_download.append(el)

#### Callbacks
#### Download data
##################
# benchmark_download_getsheet
# benchmark_download_preview
# benchmark_download


def benchmark_download_getsheet(sheet, mode="csv"):
    if sheet == "0-1 normalized data":
        out = i_class_comp.df_01_filtered_combined.copy()
        out.index = out.index.droplevel("Exp_Map")
        out = out.unstack(["Experiment", "Map"])
        out.columns = [
            "_".join(el)
            for el in out.columns.reorder_levels(
                ["Experiment", "Map", "Fraction"]
            ).values
        ]
    elif sheet == "pca coordinates":
        out = i_class_comp.df_pca.copy()
    elif sheet == "complex scatter":
        out = i_class_comp.df_distance_comp.copy()
        if mode == "csv":
            out["Cluster"] = ['"' + el + '"' for el in out["Cluster"]]
        out = (
            out.set_index(
                [
                    "Cluster",
                    "Gene names",
                    "Protein IDs",
                    "Compartment",
                    "Experiment",
                    "Map",
                ]
            )
            .drop(["Exp_Map", "merge type"], axis=1)
            .unstack(["Experiment", "Map"])
            .copy()
        )
        out.columns = ["_".join(el) for el in out.columns.values]
    elif sheet == "reproducibility":
        out = i_class_comp.distances.copy()
    elif sheet == "protein id alignment":
        out = i_class_comp.id_alignment.copy()
        out.columns = [" ".join(el) for el in out.columns.values]
    else:
        raise KeyError(sheet)
    return out


@pn.depends(i_benchmark_downloadselector.param.value, cache_run_json.param.value)
def benchmark_download_preview(file_selection, run_json):
    if not run_json:
        return "Run analysis first."
    if file_selection.endswith(".csv"):
        out = benchmark_download_getsheet(file_selection.split(".csv")[0], mode="csv")
        return pn.Column(
            f"{out.shape[0]} rows x {out.head().reset_index().shape[1]} columns",
            pn.widgets.DataFrame(out.head(10), editable=False),
        )
    else:
        return "This will download a .xlsx file with all tables as individual sheets."


@pn.depends(i_benchmark_downloadselector.param.value)
def benchmark_download(file_selection):
    i_benchmark_downloadselector.disabled = True
    i_benchmark_download.loading = True
    # set up list of files to format
    if file_selection.endswith(".csv"):
        sheets = [file_selection.split(".csv")[0]]
        mode = "csv"
    elif file_selection.endswith(".xlsx"):
        sheets = [
            "0-1 normalized data",
            "pca coordinates",
            "complex scatter",
            "reproducibility",
            "protein id alignment",
        ]
        mode = "xlsx"

    # format outputs
    out_dict = dict()
    for sheet in sheets:
        out_dict[sheet] = benchmark_download_getsheet(sheet, mode=mode)

    # return file object
    i_benchmark_download.filename = file_selection
    if file_selection.endswith(".csv"):
        sio = StringIO()
        out_dict[sheets[0]].to_csv(sio)
        sio.seek(0)
        i_benchmark_download.loading = False
        i_benchmark_downloadselector.disabled = False
        return sio
    elif file_selection.endswith(".xlsx"):
        bio = BytesIO()
        excel = pd.ExcelWriter(bio, engine_kwargs={"data_only": True})
        for sheet, df in out_dict.items():
            df.reset_index().T.reset_index().T.to_excel(
                excel, sheet_name=sheet, merge_cells=False, index=False, header=False
            )
        excel.save()
        bio.seek(0)
        i_benchmark_download.loading = False
        i_benchmark_downloadselector.disabled = False
        return bio


i_benchmark_download.callback = benchmark_download

#### Callback output positioning
#### Download data
################################
o_benchmark_downloadpreview.append(benchmark_download_preview)

#### Dashboard structure
#### Overview tab
########################

#### Layout elements
####################

#### Append layout to dashboard
###############################

#### Callbacks
##############
# update_multi_choice
# update_ExpOverview
# update_benchmark_overview


def update_multi_choice():
    i_multi_choice.options = i_class_comp.exp_names
    i_complexes_norm.options = [i_complexes_norm.options[0]] + i_class_comp.exp_names
    i_clusterwidget.options = list(i_class_comp.markerproteins.keys())
    i_clusters_for_ranking.options = list(i_class_comp.markerproteins.keys())
    i_clusters_for_ranking.value = list(i_class_comp.markerproteins.keys())
    i_multi_choice.value = i_class_comp.exp_names
    i_svm_set.options = list(
        set(
            [
                k
                for exp in i_class_comp.exp_names
                if exp in i_class_comp.svm_results.keys()
                for k in i_class_comp.svm_results[exp].keys()
            ]
        )
    )


@pn.depends(i_multi_choice.param.value, watch=True)
def update_ExpOverview(multi_choice):
    dict_analysis_parameters = {}
    for exp_name in multi_choice:
        dict_analysis_parameters[exp_name] = i_class_comp.json_dict[exp_name][
            "Analysis parameters"
        ]
    i_ExpOverview[0] = pn.widgets.DataFrame(
        pd.DataFrame.from_dict(dict_analysis_parameters), height=300, disabled=True
    )


@pn.depends(
    i_multi_choice.param.value,
    i_clusters_for_ranking.param.value,
    i_scatter_metric.param.value,
    i_scatter_consolidation.param.value,
    i_scatter_quantile.param.value,
    cache_run_json.param.value,
)
def update_benchmark_overview(
    multi_choice, clusters, metric, consolidation, quantile, run_json
):
    try:
        if run_json == False:
            return "Run analysis first!"
        if multi_choice == []:
            return "Please select experiments for comparison"
        if (
            "reproducibility consolidation" not in i_class_comp.parameters.keys()
            or i_class_comp.parameters["reproducibility metric"] != metric
            or i_class_comp.parameters["reproducibility consolidation"] != consolidation
        ):
            i_class_comp.calculate_global_scatter(metric, consolidation)
        if len(clusters) == 0:
            return "Please select clusters for comparison in the intramap scatter tab"

        fig = i_class_comp.plot_overview(multi_choice, clusters, quantile)

        return pn.Column(
            pn.Pane(fig, config=plotly_config),
            "Note that this overview figure is affected by the settings you choose in the intra- and inter-map scatter tabs.",
        )
    except:
        return traceback.format_exc()


#### Callback output positioning
################################


comparison_tabs.clear()
comparison_tabs.append(
    (
        "Overview",
        pn.Column(
            "Once you have run the analysis you can find different benchmarking outputs here and dive into the data.",
            update_benchmark_overview,
        ),
    )
)
comparison_tabs.append(("PCA maps", lo_benchmark_pca))
comparison_tabs.append(("Depth & Coverage", lo_benchmark_depth))
comparison_tabs.append(("Intermap scatter", lo_benchmark_intermap))
comparison_tabs.append(("Intramap scatter", comparison_tab_bp))
comparison_tabs.append(("SVM Analysis", lo_benchmark_SVMs))
comparison_tabs.append(("Compare profiles", update_profile_comparison))
comparison_tabs.append(("Download data", lo_benchmark_download))


# [<div style="text-align: right; font-size: 8pt">back to top</div>](#TOC)
# ## Benchmark tab upload section<a id="benchmarkupload"></a>

# In[ ]:


dashboard_benchmark.objects = [
    pn.Card(
        objects=[], header="## Manage data", name="manage_data", height_policy="fit"
    ),
    pn.Row(objects=[], name="benchmark_output", height_policy="fit"),
]

#### Manage data Card layout
# This accesses mem_available_datasets and mem_benchmark.
####

## Adding datasets row
i_upload_collection = pn.widgets.FileInput(name="Upload collection")
btn_load_reference = pn.widgets.Button(name="Load", width=100)
i_load_reference = pn.widgets.Select(
    options=pkg_resources.resource_listdir("domaps", "referencedata"),
    value=None,
    width=200,
)
lo_add_datasets = pn.Row(
    objects=[
        "**Upload collection from file (.json):**",
        i_upload_collection,
        "**Add reference set:**",
        i_load_reference,
        btn_load_reference,
    ]
)

## Selection checkbox
i_dfs_available = pn.widgets.CheckBoxGroup(
    options=[], value=[], name="Datasets available"
)
lo_dfs_available = pn.WidgetBox(
    objects=["**Datasets available**", i_dfs_available], width=250
)

## Management button group
i_coll_download = pn.widgets.FileDownload(
    label="Download selected as collection (.json)",
    filename="AnalysedDatasets.json",
    disabled=True,
)
btn_coll_editnames = pn.widgets.Button(
    name="Edit names and comments", disabled=True
)  # move from management
btn_coll_reannotate = pn.widgets.Button(
    name="Reannotate genes/organelles/complexes", disabled=True
)  # new functionality
btn_coll_runmain = pn.widgets.Button(
    name="Align and analyse selected datasets", button_type="success", disabled=True
)  # change from main comparison
btn_coll_dropmem = pn.widgets.Button(
    name="Drop selected datasets from memory", button_type="danger", disabled=True
)  # move from top of page
lo_coll_buttons = pn.WidgetBox(
    objects=[
        i_coll_download,
        #    btn_coll_editnames,
        #    btn_coll_reannotate,
        btn_coll_runmain,
        btn_coll_dropmem,
    ],
    width=300,
)

## Interaction pane
lo_instructions_datamanagement = pn.Card(
    pn.pane.Markdown(textfragments["coll_status_default"]),
    header="**Explanation**",
    width=400,
    collapsed=True,
)

lo_instructions_error_messages = pn.Card(
    pn.pane.Markdown(textfragments["benchmark_error_messages"]),
    header="**Common error messages**",
    width=400,
    collapsed=True,
)


o_status_datamanagement = pn.pane.Markdown(width=400)


def set_status_datamanagement(x, append=False):
    if not append:
        o_status_datamanagement.object = x + "<br><br>"
    else:
        o_status_datamanagement.object += x + "<br><br>"
    if (
        x.startswith("Traceback")
        or o_status_datamanagement.object.count("<br><br>") > 1
    ):
        resize(
            dashboard_benchmark.objects[
                [i.name for i in dashboard_benchmark].index("manage_data")
            ]
        )
    if DEBUG:
        time.sleep(0.3)
    o_status_datamanagement.object = x


set_status_datamanagement("Step 1: Add datasets")
o_dynamic_collectionmanagement = pn.Row()
lo_coll_interactions = pn.Column(
    objects=[
        lo_instructions_datamanagement,
        # lo_instructions_error_messages, # currently empty
        o_status_datamanagement,
        o_dynamic_collectionmanagement,
    ]
)

## Assemble collection management row
lo_manage_collection = pn.Row(
    objects=[lo_dfs_available, lo_coll_buttons, lo_coll_interactions]
)

#### Append elements to manage data row
dashboard_benchmark.objects[
    [i.name for i in dashboard_benchmark].index("manage_data")
].objects = []
for el in [
    pn.Pane(textfragments["benchmark_management_top"], width=600),
    lo_add_datasets,
    lo_manage_collection,
]:
    dashboard_benchmark.objects[
        [i.name for i in dashboard_benchmark].index("manage_data")
    ].append(el)

#### Management callbacks
# upload_collection #Done
# load_reference
# coll_activatebuttons # Done
# coll_downloadjson
# coll_editnames
# coll_reannotate
# coll_runmain #Button change in place
# coll_dropmem #Done
####

lock_collection_change = [
    # btn_coll_editnames,
    # btn_coll_reannotate,
    i_coll_download,
    btn_coll_runmain,
    btn_coll_dropmem,
    i_dfs_available,
    i_upload_collection,
    btn_load_reference,
]


@pn.depends(i_upload_collection.param.value, watch=True)
def upload_collection(file):
    """
    This callback adds the datasets from a .json collection to the global memory object
    and updates interface elements accordingly.
    """
    if file == None:
        return
    # deactivate interface
    for el in lock_collection_change:
        el.disabled = True
    try:
        set_status_datamanagement("Loading data ...")

        status = ""
        json_loaded = json.load(BytesIO(file))

        # Check if experiment names are still free
        renamed_exps = []
        n_sets = 0
        keys = list(json_loaded.keys())
        for exp in keys:
            if exp in i_dfs_available.options:
                renamed_exps.append(exp)
                json_loaded[exp + i_upload_collection.filename] = json_loaded.pop(exp)
            n_sets += 1
        if len(renamed_exps) != 0:
            status += f"""These datasets were already available:<br><br>{", ".join(renamed_exps)}
            <br><br>Please rename them prior to analysis.<br><br><br><br>"""

        # Load datasets into memory
        keys = list(json_loaded.keys())
        for exp in keys:
            set_status_datamanagement(f"Loading dataset {exp} ...")
            mem_available_datasets[exp] = json_loaded[exp]

        # Adjust list of available dataset
        i_dfs_available.options = i_dfs_available.options + keys
        i_dfs_available.value = i_dfs_available.value + keys
        resize(lo_dfs_available)

        status += (
            f"Loaded **{n_sets}** datasets from file **{i_upload_collection.filename}**"
        )
        set_status_datamanagement(status)
    except Exception:
        set_status_datamanagement(traceback.format_exc())
    finally:
        # reactivate interface
        for el in lock_collection_change:
            if type(el) != pn.widgets.button.Button:
                el.disabled = False
        coll_activatebuttons(i_dfs_available.value)


def load_reference(event):
    upload_collection(
        pkg_resources.resource_stream(
            "domaps", "referencedata/" + i_load_reference.value
        ).read()
    )


btn_load_reference.on_click(load_reference)


@pn.depends(i_dfs_available.param.value, watch=True)
def coll_activatebuttons(v):
    """
    Activate/Deactivate buttons based on selection of available datasets.
    """
    try:
        btn_load_reference.disabled = False
        if len(v) == 0:
            i_coll_download.disabled = True
            btn_coll_reannotate.disabled = True
            btn_coll_runmain.disabled = True
            btn_coll_dropmem.disabled = True
        else:
            btn_coll_reannotate.disabled = False
            i_coll_download.disabled = False
            btn_coll_runmain.disabled = False
            btn_coll_dropmem.disabled = False
        if len(i_dfs_available.options) == 0:
            btn_coll_editnames.disabled = True
        else:
            btn_coll_editnames.disabled = False
            pass
    except Exception:
        set_status_datamanagement(traceback.format_exc())


@pn.depends(i_dfs_available.param.value)
def coll_downloadjson(dfs):
    sio = StringIO()
    json.dump(
        {k: mem_available_datasets[k] for k in dfs}, sio, indent=4, sort_keys=True
    )
    sio.seek(0)
    return sio


i_coll_download.callback = coll_downloadjson


def coll_runmain(event):
    """
    Runs main analysis
    """
    # deactivate interface
    for el in lock_collection_change:
        el.disabled = True
    try:
        if btn_coll_runmain.button_type == "success":
            set_status_datamanagement("Aligning and analysing data ...")
            cache_run_json.value = False
            #### Main execution of the comparison
            loading_status_comparison.objects = [loading_comparison]
            selection = i_dfs_available.value
            global i_class_comp
            i_class_comp = domaps.SpatialDataSetComparison(
                ref_exp=selection[0]
            )  # , clusters_for_ranking=protein_cluster, organism=i_organism_comparison.value)
            i_class_comp.json_dict = {k: mem_available_datasets[k] for k in selection}
            set_status_datamanagement("Aligning data ...", append=True)
            i_class_comp.read_jsonFile()
            set_status_datamanagement("Analysing intra-map scatter ...", append=True)
            i_class_comp.calc_biological_precision()
            i_class_comp.get_complex_coverage()
            update_multi_choice()
            set_status_datamanagement("Running PCA ...", append=True)
            i_class_comp.perform_pca_comparison()
            i_pca_comp_ncomp.end = len(i_class_comp.df_01_filtered_combined.columns)
            i_compare_compartment.options = list(
                set(
                    i_class_comp.df_01_filtered_combined.index.get_level_values(
                        "Compartment"
                    )
                )
            )
            loading_status_comparison.objects = []
            m_diverget_fractions.object = (
                ""
                if not i_class_comp.mixed
                else "**Caution: You are comparing experiments with differently labelled fractions. This does not affect distance metrics, but the PCA and profile plots.**"
            )
            m_diverget_fractions.background = (
                None if not i_class_comp.mixed else "salmon"
            )
            cache_run_json.value = True
            set_status_datamanagement("Comparison finished!", append=True)

            #### Switch button mode
            btn_coll_runmain.button_type = "danger"
            btn_coll_runmain.name = "Reset analysis to make new selection"
            set_status_datamanagement(
                "Next step: Use interface below to evaluate and download benchmark results.",
                append=True,
            )
        elif btn_coll_runmain.button_type == "danger":
            set_status_datamanagement("Resetting data analysis ...")
            cache_run_json.value = False

            btn_coll_runmain.button_type = "success"
            btn_coll_runmain.name = "Align and analyse selected datasets"
            set_status_datamanagement("Analysis results have been reset.")
        else:
            pass

    except Exception:
        set_status_datamanagement(traceback.format_exc())
        cache_run_json.value = False

    finally:
        # reactivate interface
        if btn_coll_runmain.button_type == "danger":
            for el in [i_coll_download, btn_coll_runmain, btn_coll_editnames]:
                el.disabled = False
        elif btn_coll_runmain.button_type == "success":
            for el in lock_collection_change:
                if type(el) != pn.widgets.button.Button:
                    el.disabled = False
            coll_activatebuttons(i_dfs_available.value)
        else:
            for el in lock_collection_change:
                if type(el) != pn.widgets.button.Button:
                    el.disabled = False
            coll_activatebuttons(i_dfs_available.value)


btn_coll_runmain.on_click(coll_runmain)


def coll_dropmem(event):
    """
    Drops selected datasets from collection stored in RAM.
    """
    # deactivate interface
    for el in lock_collection_change:
        el.disabled = True
    try:
        set_status_datamanagement("Deleting data ...")

        keys = list(i_dfs_available.value)
        for exp in keys:
            set_status_datamanagement(f"Deleting dataset {exp} ...")
            del mem_available_datasets[exp]

        # Adjust list of available dataset
        i_dfs_available.options = [
            el for el in i_dfs_available.options if el not in keys
        ]
        i_dfs_available.value = [el for el in i_dfs_available.value if el not in keys]
        resize(lo_dfs_available)
        i_dfs_available.disabled = True

        set_status_datamanagement(
            f"Deleted **{len(keys)}** datasets **{', '.join(keys)}**"
        )

    except Exception:
        set_status_datamanagement(traceback.format_exc())
    finally:
        # reactivate interface
        for el in lock_collection_change:
            if type(el) != pn.widgets.button.Button:
                el.disabled = False
        coll_activatebuttons(i_dfs_available.value)


btn_coll_dropmem.on_click(coll_dropmem)


#### Benchmark output
@pn.depends(cache_run_json.param.value)
def display_benchmark_output(run_json):
    if run_json:
        return pn.Column(
            "## Benchmark results",
            "Select dataset overlap to plot:",
            i_multi_choice,
            comparison_tabs,
            sizing_mode="stretch_width",
        )
    else:
        return "Select data and run analysis."


#### Append callback output to benchmark output
dashboard_benchmark.objects[
    [i.name for i in dashboard_benchmark].index("benchmark_output")
].objects = []
for el in [display_benchmark_output]:
    dashboard_benchmark.objects[
        [i.name for i in dashboard_benchmark].index("benchmark_output")
    ].append(el)


# [<div style="text-align: right; font-size: 8pt">back to top</div>](#TOC)
# ## Code interactions<a id="interactions"></a>

# In[ ]:


# In case of loading a json comparison larger than 80 MB
# with open(r"C:\Documents\AnalysedDatasets.json", "br") as file:
#    i_upload_collection.value = file.read()
#    i_upload_collection.filename = "bla"


# In[ ]:


# In case of loading a dingle file larger than 80 MB
# i_FileConfig._content.file.filename = "proteinGroups.txt"
# with open(r"path\proteinGroups.txt", "br") as file:
#    i_FileConfig._content.file.value = file.read()


# In[ ]:


# Set the order of the multi choice widget manually
# i_multi_choice.value=["21 min", "44 min", "100 min", "DDA"]
