#!/usr/bin/env python
# coding: utf-8

# ## Library imports

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

try:
    type(domaps)
    print("reloading library")
    domaps = reload(domaps)
except Exception as ex:
    print("loading library first time")
    import domaps
    
from datetime import datetime


# ## panel and plotly settings and customization

# In[ ]:


css = """
.detail_menu .bk-headers-wrapper{
  border-bottom: 2px solid #0a5da8 !important;
  margin-bottom: 10px;
  min-width: 1300px;
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
  min-width: 1300px;
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
.content-width{
  min-width: 1300px;
}

.bk-tabs-header{
  min-width: 1300px;
}

.card-title:first-child{
  font-size: 13px;
}

.button-main .bk-btn{
  font-size: 120%;
  font-weight: bold;
}
"""
pn.extension(raw_css=[css])

plotly_config={
      'toImageButtonOptions': {
            'format': 'svg', # one of png, svg, jpeg, webp
            'filename': 'figure'
      }
}

def resize(el):
    try:
        el.append(None)
        el.pop(-1)
    except:
        pass


# ## Global variables

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

DEBUG = True
MAX_SIZE_MB = 80
CONTENT_WIDTH = 1000


# ## App interface

# In[ ]:


#### Served panel object:
app = pn.GridSpec()#sizing_mode="stretch_both", margin=0)
app[0,0] = pn.Spacer(background="white", margin=0) #"#DDDDDD"
app[0,9] = pn.Spacer(background="white", margin=0) #"#DDDDDD"

#### Insert main content container
app_center = pn.Column(pn.Row(pn.Pane("# QC tool for Spatial Proteomics", width = 600),
                              pn.layout.HSpacer(),
                              margin=10),
                       pn.Row(name="main_content"),
                       pn.Spacer(background="#DDDDDD", height=100, margin=0)
                      )
app[0,1:8] = app_center

#### Insert main menu tab object
app_tabs = pn.Tabs(margin=10, css_classes=["content-width", "main_menu"], dynamic=True)
app_center.objects[[i.name for i in app_center].index("main_content")] = app_tabs

#### Append individual dashboards
## Home
dashboard_home = pn.Column(
    "Interface loading ...",
    name="home", css_classes=["content-width"])
app_tabs.append(("Home", dashboard_home))

## Single analysis
dashboard_analysis = pn.Column(
    "Interface loading ...",
    name="analysis", css_classes=["content-width"])
app_tabs.append(("Analysis", dashboard_analysis))
analysis_tabs = pn.Tabs(margin=10, css_classes=["content-width", "detail_menu"], dynamic=True)

## Benchmark
dashboard_benchmark = pn.Column(
    "Interface loading ...",
    name="benchmark", css_classes=["content-width"])
app_tabs.append(("Benchmark", dashboard_benchmark))
comparison_tabs = pn.Tabs(margin=10, css_classes=["content-width", "detail_menu"], dynamic=True)

## Manage datasets
#dashboard_MissclassificationMatrix = pn.Column(
#    "Please, upload a file first and press 'Analyse clusters'",
#    name="SVM Analysis", css_classes=["content-width"])
#dashboard_amendment = pn.Column(
#    "Please, upload a json file first",
#    name="Renaming", css_classes=["content-width"])
dashboard_manageDatasets = pn.Column(
    "Interface loading ...",
    name="Manage datasets", css_classes=["content-width"])
#app_tabs.append(("Manage Datasets", dashboard_manageDatasets))

#amendment_tabs = pn.Tabs(margin=10, css_classes=["content-width", "detail_menu"], dynamic=True)
#amendment_tabs.append(("Change Experiment name and comment", dashboard_amendment))
#amendment_tabs.append(("SVM Upload", dashboard_MissclassificationMatrix))

## About
app_tabs.append(("About", pn.Row(pn.Pane(textfragments["about_intro"], width=1000))))


# ## App serving
# Switch cells below between markup and code to set up for server hosting from the command line (app.servable) vs. local hosting from python.

# try:
#     server.stop()
# except Exception:
#     print("First server startup")
# server = app.show(port=5067, websocket_max_message_size=MAX_SIZE_MB*1024*1024,
#                   http_server_kwargs={'max_buffer_size': MAX_SIZE_MB*1024*1024})

# In[ ]:


app.servable()


# ## Cell structuring
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


# ## Home tab

# In[ ]:


#### Dashboard structure
########################
# already defined as single column
dashboard_home.objects = []

#### Layout elements
####################
lo_home_intro = pn.Pane(textfragments["home_intro"], width=CONTENT_WIDTH)
btn_home_analysesingle = pn.widgets.Button(name="Format and analyse single experiment",
                                           button_type="success", width=400, height=50,
                                           css_classes=["button-main"])
lo_home_singleinstructions = pn.Column(
    pn.Pane(textfragments["home_single_shortinstructions"], width=CONTENT_WIDTH),
    pn.Card("Add screenshot tutorial here.", header="Tutorial", width=CONTENT_WIDTH,
            name="tutorial_single", collapsed=True)
)
btn_home_benchmark = pn.widgets.Button(name="Benchmark multiple experiments",
                                       button_type="success", width=400, height=50,
                                       css_classes=["button-main"])
lo_home_benchmarkinstructions = pn.Column(
    pn.Pane(textfragments["home_benchmark_shortinstructions"], width=CONTENT_WIDTH),
    pn.Card("Add screenshot tutorial here.", header="Tutorial", width=CONTENT_WIDTH,
            name="tutorial_benchmark", collapsed=True)
)

#### Append layout to dashboard
###############################
for el in [lo_home_intro,
           btn_home_analysesingle, lo_home_singleinstructions,
           btn_home_benchmark, lo_home_benchmarkinstructions]:
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


# ## Analysis tab
# - File config
# - Analysis output

# In[ ]:


#### Dashboard structure
########################
dashboard_analysis.objects = [
    pn.Column(name="file_config"),
    pn.Column(name="analysis_output")
]

#### Layout elements
#### File config
####################
lo_read_file = pn.Card(header="### Upload configuration", min_width=400)
lo_config_instructions = pn.Card(
    pn.Pane(textfragments["upload_details"]), header="### Details on configuring your data", width=400)

i_file = pn.widgets.FileInput(name="Upload file")
loading_status = pn.Row()
idle = pn.indicators.LoadingSpinner(value=False, width=100, height=100, color="primary")
loading = pn.indicators.LoadingSpinner(value=True, width=100, height=100, color="primary")
i_acquisition = pn.widgets.Select(options=["LFQ6 - MQ", "LFQ5 - MQ", "LFQ6 - Spectronaut", "LFQ5 - Spectronaut",
                                           "SILAC - MQ", "Custom"], name="Acquisition", width=300)
i_organism = pn.widgets.Select(options=[el.split(".csv")[0] for el in
                                        pkg_resources.resource_listdir("domaps", "annotations/complexes")],
                               name="Organism", width=300, value="Homo sapiens - Uniprot")
i_comment = pn.widgets.input.TextAreaInput(
    name="Additional Comments",
    placeholder="Write any kind of information associated with this dataset here...")
analysis_status = pn.Pane("", width=300)
filereading_status = pn.Pane("No data import yet", width=300)
i_expname = pn.widgets.TextInput(name="Experiment Name", placeholder="Enter your experiment name here here...")
i_custom_namepattern = pn.widgets.TextInput(name="Customized Name Pattern",
                                            placeholder="Enter a string here...e.g.: .* (?P<rep>.*)_(?P<frac>.*)")
regex_pattern = {
    "(?P<rep>.*)_(?P<frac>.*)" : ["Spectronaut MAP1_03K"],
    ".* (?P<rep>.*)_(?P<frac>.*)" : ["MAP1_03K","MAP3_03K"],
    ".* (?P<cond>.*)_(?P<rep>.*)_(?P<frac>.*)" : ["EGF_rep1_06K","EGF_rep3_06K"],
    ".* (?P<cond>.*)_(?P<frac>.*)_(?P<rep>.*)" : ["Control_Mem_1", "Control_Cyt_1"],
    ".* (?P<cond>.*)_(?P<frac>.*_.*)_(?P<rep>.*)" : ["4h_mem_frac3_1", "co_mem_frac2_2"]
    }
i_name_pattern = pn.widgets.Select(name="Name pattern",
                                   options=[".* (?P<rep>.*)_(?P<frac>.*)", ".* (?P<frac>.*)_(?P<rep>.*)",
                                            ".* (?P<cond>.*)_(?P<rep>.*)_(?P<frac>.*)",
                                            ".* (?P<cond>.*)_(?P<frac>.*)_(?P<rep>.*)",
                                            ".* (?P<cond>.*)_(?P<frac>.*_.*)_(?P<rep>.*)",
                                            "(?P<rep>.*)_(?P<frac>.*)", "Custom"])
i_pattern_examples = pn.widgets.Select(name = "Examples", options=regex_pattern[i_name_pattern.value])

button_analysis = pn.widgets.Button(name="Analyse dataset", width=50)

# gene reannotation
i_reannotate_genes = pn.widgets.Select(
    options=["don't reannotate", "from uniprot fasta headers", "from uniprot.org", "from uniprot tab download"],
    name="Reannotate genes", value="from uniprot fasta headers")
i_reannotate_genes_source = pn.widgets.Select(options=["upload custom file"], value="select last",
                                              name="Annotation source:")
i_reannotate_genes_file = pn.widgets.FileInput(name="Upload gene annotation file")

# custom acqusition
i_custom_cids = pn.widgets.Select(name="column containing protein ids")
i_custom_cgenes = pn.widgets.Select(name="column containing gene names")
i_custom_cdata = pn.widgets.MultiSelect(name="columns containing data")
i_custom_normalized = pn.widgets.Checkbox(name="profiles are already 0-1 normalized (profile sum = 1)")

# LFQ acquisition
i_consecutiveLFQi = pn.widgets.IntSlider(name="Consecutive LFQ intensities", start=1, end=10, step=1, value=4)
i_summed_MSMS_counts = pn.widgets.IntSlider(name="Summed MS/MS counts ≥ 2n x number of fractions",
                                            start=0, end=10, step=1, value=2)

# SILAC acquisition
i_RatioHLcount = pn.widgets.IntSlider(name="Quantification events (Ratio H/L Count)", start=1, end=10, step=1, value=2)
i_RatioVariability = pn.widgets.IntSlider(name="Ratio H/L variability [%]", start=0, end=100, step=5, value=30)
#### Append layout to dashboard
#### File config
###############################
dashboard_analysis.objects[[el.name for el in dashboard_analysis].index("file_config")].objects = []
for el in [i_file,
    pn.Row(
        pn.Column(lo_read_file, analysis_status, loading_status),
        lo_config_instructions
    )
]:
    dashboard_analysis.objects[[el.name for el in dashboard_analysis].index("file_config")].append(el)

#### Callbacks
#### File config
################
# response_pattern
# response_acquisition
# response_reannotation
# read_file
# execution

cache_uploaded = pn.widgets.Checkbox(value=False)
cache_run = pn.widgets.Checkbox(value=False)
#define widgets that should be disbled after run==True
wdgts = [i_acquisition,i_name_pattern,i_expname, i_pattern_examples, button_analysis, i_expname, i_organism,
         i_consecutiveLFQi, i_summed_MSMS_counts, i_RatioHLcount, i_RatioVariability, i_comment,
         i_custom_cids, i_custom_cgenes, i_custom_cdata, i_custom_normalized,
         i_reannotate_genes, i_reannotate_genes_file, i_reannotate_genes_source]


@pn.depends(i_file.param.value)
def read_file(file):
    if file is None:
        filereading_status = "No file is uploaded"
        cache_uploaded.value = False
        return filereading_status
    else:
        cache_uploaded.value = False
        try:
            if i_file.filename[-3:] == "xls" or i_file.filename[-3:] == "txt":
                df_original = pd.read_csv(BytesIO(file), sep="\t", comment="#", nrows=5, usecols=lambda x: bool(re.match(domaps.SpatialDataSet.regex["imported_columns"], x)), low_memory=True)
                if len(df_original.columns) < 5:
                    df_original = pd.read_csv(BytesIO(file), sep="\t", comment="#", nrows=5)
            elif i_file.filename[-3:] == "csv":
                df_original = pd.read_csv(BytesIO(file), sep=",", comment="#", nrows=5, usecols=lambda x: bool(re.match(domaps.SpatialDataSet.regex["imported_columns"], x)), low_memory=True)
                if len(df_original.columns) < 5:
                    df_original = pd.read_csv(BytesIO(file), sep=",", comment="#", nrows=5)
            else:
                return pn.Column("Upload either csv, xls, or txt formatted files.")
            cache_uploaded.value = True
            analysis_status.object = "No analysis run yet"
            for wdgt in wdgts:
                wdgt.disabled = False

            return pn.Column(pn.Row(pn.widgets.DataFrame(df_original, height=200, width=600, disabled=True)),#i_class.
                             pn.Row(i_expname), 
                             pn.Row(i_organism, i_acquisition),
                             pn.Row(i_name_pattern, response_pattern),
                             pn.Row(i_reannotate_genes, response_reannotation),
                             pn.Row(response_acquisition), 
                             pn.Row(i_comment),
                             pn.Row(button_analysis))

        except Exception: 
            filereading_status = traceback.format_exc()
            cache_uploaded.value = False
            return filereading_status

@pn.depends(i_name_pattern.param.value)
def response_pattern(name_pattern):
    if name_pattern == "Custom":
        return i_custom_namepattern
    else:
        example_for_name_pattern = regex_pattern[name_pattern]
        i_pattern_examples.options = example_for_name_pattern
        return i_pattern_examples


@pn.depends(i_acquisition.param.value)
def response_acquisition(acquisition):
    if acquisition == "SILAC - MQ":
        return pn.Column(pn.Row(pn.Pane("Stringency filtering")), pn.Row(i_RatioHLcount, i_RatioVariability))
    elif acquisition == "Custom":
        try:
            if i_file.filename[-3:] == "xls" or i_file.filename[-3:] == "txt":
                columns = list(pd.read_csv(BytesIO(i_file.value), sep="\t", comment="#", nrows=2).columns)
            elif i_file.filename[-3:] == "csv":
                columns = list(pd.read_csv(BytesIO(i_file.value), sep=",", comment="#", nrows=2).columns)
            else:
                return "Failed to read file"
            i_custom_cids.options = columns
            i_custom_cgenes.options = columns
            i_custom_cdata.options = columns
            return pn.Column("Please note that the custom upload is still under development and currently assumes 0-1 normalized data.",
                             pn.Row(pn.Column(i_custom_cids, i_custom_cgenes), i_custom_cdata),
                             #i_custom_normalized
                            )
        except Exception:
            return traceback.format_exc()
    else:
        return pn.Column(pn.Row(pn.Pane("Stringency filtering")), pn.Row(i_consecutiveLFQi, i_summed_MSMS_counts))

@pn.depends(i_reannotate_genes.param.value, i_reannotate_genes_source.param.value)
def response_reannotation(mode, source):
    if mode in ["don't reannotate", "from uniprot.org"]:
        i_reannotate_genes_source.options = i_reannotate_genes_source.options[0:1]
        if source not in i_reannotate_genes_source.options:
            i_reannotate_genes_source.value = i_reannotate_genes_source.options[-1]
        return None
    elif mode == "from uniprot fasta headers":
        stored = [el for el in pkg_resources.resource_listdir("domaps", "annotations/idmapping") if el.endswith(".txt")]
        i_reannotate_genes_source.options = i_reannotate_genes_source.options[0:1]+stored
        if source not in i_reannotate_genes_source.options:
            i_reannotate_genes_source.value = i_reannotate_genes_source.options[-1]
    elif mode == "from uniprot tab download":
        stored = [el for el in pkg_resources.resource_listdir("domaps", "annotations/idmapping") if el.endswith(".tab")]
        i_reannotate_genes_source.options = i_reannotate_genes_source.options[0:1]+stored
        if source not in i_reannotate_genes_source.options:
            i_reannotate_genes_source.value = i_reannotate_genes_source.options[-1]
    else:
        return f"Reannotation mode {mode} is not implemented."
    if source == i_reannotate_genes_source.options[0]:
        return pn.Column(i_reannotate_genes_source, i_reannotate_genes_file)
    else:
        return i_reannotate_genes_source

def execution(event):
    if cache_uploaded.value == False:
        analysis_status.object = "Please upload a file first"
    elif i_expname.value == "":
        analysis_status.object = "Please enter an experiment name first"
    elif i_reannotate_genes_file.value == "upload custom file" and i_reannotate_genes_source.value is None:
        analysis_status.object = "Please upload your custom file for gene reannotation"
    else:
        lo_config_instructions.collapsed = True
        lo_read_file.collapsed = True
        output_layoutpos = [el.name for el in dashboard_analysis].index("analysis_output")
        dashboard_analysis.objects[output_layoutpos].objects = []
        cache_run.value = False
        for wdgt in wdgts:
            wdgt.disabled = True
        #if you did already your comparison, but add another experiment afterwards - without reloading your
        #AnylsedDatasets.json
        try:
            loading_status.objects = [loading]
            analysis_status.object = "Analysis in progress"
            
            # get naming pattern
            if i_name_pattern.value == "Custom":
                namePattern = i_custom_namepattern.value
            else:
                namePattern = i_name_pattern.value
            
            # get gene reannotation
            reannotation_source = {}
            if i_reannotate_genes.value == "don't reannotate":
                reannotate = False
            else:
                reannotate = True
                if i_reannotate_genes.value == "from uniprot.org":
                    reannotation_source["source"] = "uniprot"
                elif i_reannotate_genes.value == "from uniprot fasta headers":
                    reannotation_source["source"] = "fasta_headers"
                    if i_reannotate_genes_source.value == "upload custom file":
                        idmapping = [el.decode('UTF-8').strip()
                                     for el in BytesIO(i_reannotate_genes_file.value).readlines()]
                    else:
                        idmapping = [el.decode('UTF-8').strip()
                                     for el in pkg_resources.resource_stream(
                            "domaps", 'annotations/idmapping/{}'.format(i_reannotate_genes_source.value)).readlines()]
                    reannotation_source["idmapping"] = idmapping
                elif i_reannotate_genes.value == "from uniprot tab download":
                    reannotation_source["source"] = "tsv"
                    if i_reannotate_genes_source.value == "upload custom file":
                        idmapping = pd.read_csv(BytesIO(i_reannotate_genes_file.value), sep='\t',
                                                usecols=["Entry", "Gene names  (primary )"])
                    else:
                        idmapping = pd.read_csv(pkg_resources.resource_stream(
                            "domaps", 'annotations/idmapping/{}'.format(i_reannotate_genes_source.value)), sep='\t',
                                                usecols=["Entry", "Gene names  (primary )"])
                    reannotation_source["idmapping"] = idmapping
                else:
                    raise ValueError(i_reannotate_genes.value)
            
            global i_class
            if i_acquisition.value == "SILAC":
                i_class = domaps.SpatialDataSet(i_file.filename, i_expname.value, i_acquisition.value,
                                                comment=i_comment.value, name_pattern=namePattern,
                                                organism=i_organism.value,
                                                reannotate_genes=reannotate, reannotation_source=reannotation_source,
                                                RatioHLcount=i_RatioHLcount.value,
                                                RatioVariability=i_RatioVariability.value)
            if i_acquisition.value == "Custom":
                i_class = domaps.SpatialDataSet(i_file.filename, i_expname.value, i_acquisition.value,
                                                comment=i_comment.value, name_pattern=namePattern,
                                                organism=i_organism.value,
                                                reannotate_genes=reannotate, reannotation_source=reannotation_source,
                                                custom_columns = {"ids": i_custom_cids.value,
                                                                  "genes": i_custom_cgenes.value,
                                                                  "data": i_custom_cdata.value},
                                                custom_normalized = i_custom_normalized.value)
            else:
                i_class = domaps.SpatialDataSet(i_file.filename, i_expname.value, i_acquisition.value,
                                                comment=i_comment.value, name_pattern=namePattern,
                                                organism=i_organism.value,
                                                reannotate_genes=reannotate, reannotation_source=reannotation_source,
                                                consecutiveLFQi=i_consecutiveLFQi.value,
                                                summed_MSMS_counts=i_summed_MSMS_counts.value)
            #analysis_status.object = "Data Reading"
            #i_class.data_reading(content=BytesIO(i_file.value))
            #analysis_status.object = "Data Processing"
            #i_class.processingdf()
            #update_object_selector(i_mapwidget, i_clusterwidget)
            #i_class.quantity_profiles_proteinGroups()
            #analysis_status.object = "PCA"
            #i_class.perform_pca()
            #analysis_status.object = "Calculating Manhattan Distance"
            #i_class.calc_biological_precision()
            #analysis_status.object = "Writing Overview Table"
            #i_class.results_overview_table()
            #analysis_status.object = "Writing Analysed Dataset Dictionary"
            #loading_status.objects = []
            i_class.run_pipeline(content=BytesIO(i_file.value), progressbar=analysis_status)
            analysis_status.object = "Analysis finished!"
            update_object_selector(i_mapwidget, i_clusterwidget)
            loading_status.objects = []
            
            dashboard_analysis.objects[output_layoutpos].append(i_clusterwidget)
            dashboard_analysis.objects[output_layoutpos].append(i_mapwidget)
            dashboard_analysis.objects[output_layoutpos].append(analysis_tabs)
            mem_available_datasets[i_class.expname] = i_class.analysed_datasets_dict[i_class.expname]
            try:
                i_dfs_available.options.append(i_class.expname)
                i_dfs_available.value.append(i_class.expname)
                coll_activatebuttons(i_dfs_available.value)
            except:
                pass
            cache_run.value = True
        except:
            for wdgt in wdgts:
                wdgt.disabled = False
            loading_status.objects = [""]
            analysis_status.object = traceback.format_exc()
            cache_run.value = False
button_analysis.on_click(execution)

#### Callback output positioning
#### File config
################################
lo_read_file.append(read_file)



i_logOR01_selection = pn.widgets.Select(options=["0/1 normalized data", "log transformed data",
                                                 "stringency filtered raw data", "Overview - cluster distances"],
                                        name="Select type of data for download", width=300)

i_clusterwidget = pn.widgets.Select(options=["Proteasome", "Lysosome"], name="Cluster of interest", width=300)
i_mapwidget = pn.widgets.Select(options=["Map1", "Map2"], name="Map of interest", width=300)

i_collapse_maps_PCA = pn.widgets.Checkbox(value=False, name="Collapse maps")

i_x_vs_yAxis_PCA = {
    "PC1" : ["PC3", "PC2"],
    "PC2" : ["PC1", "PC3"],
    "PC3" : ["PC1", "PC2"],
    }

i_xAxis_PCA = pn.widgets.Select(name="X-Axis", options=["PC1", "PC2","PC3"])
i_yAxis_PCA = pn.widgets.Select(name="Y-Axis", options=i_x_vs_yAxis_PCA[i_xAxis_PCA.value])

i_xAxis_PCA_comp = pn.widgets.Select(name="X-Axis", options=["PC1", "PC2","PC3"])
i_yAxis_PCA_comp = pn.widgets.Select(name="Y-Axis", options=i_x_vs_yAxis_PCA[i_xAxis_PCA_comp.value])

@pn.depends(i_xAxis_PCA_comp.param.value, watch=True)
def custimization_PCA_comp(xAxis_PCA_comp):
    yAxis_PCA_comp = i_x_vs_yAxis_PCA[xAxis_PCA_comp]
    i_yAxis_PCA_comp.options = yAxis_PCA_comp
    #return i_yAxis_PCA_comp

@pn.depends(i_xAxis_PCA.param.value, watch=True)
def custimization_PCA(xAxis_PCA):
    yAxis_PCA = i_x_vs_yAxis_PCA[xAxis_PCA]
    i_yAxis_PCA.options = yAxis_PCA
    #return i_yAxis_PCA
        


## Analysis tab

def update_object_selector(i_mapwidget, i_clusterwidget):
    i_mapwidget.options = list(i_class.map_names)
    i_clusterwidget.options = list(i_class.markerproteins.keys())


@pn.depends(i_mapwidget.param.value, cache_run.param.value, i_collapse_maps_PCA.param.value, i_clusterwidget.param.value, i_xAxis_PCA.param.value, i_yAxis_PCA.param.value)
def update_data_overview(mapwidget, run, collapse_maps_PCA, clusterwidget, xAxis_PCA, yAxis_PCA):
    try:
        if run == True:
            pca_plot = i_class.plot_global_pca(map_of_interest=mapwidget , cluster_of_interest=clusterwidget, x_PCA=xAxis_PCA, y_PCA=yAxis_PCA, collapse_maps=collapse_maps_PCA)
            if i_acquisition.value != "Custom":
                log_histogram = i_class.plot_log_data()
            else:
                log_histogram = ""
            visualization_map = pn.Column(
                pn.Row(i_collapse_maps_PCA),
                pn.Row(pn.Pane(pca_plot, width=1000, config=plotly_config)),
                pn.Row(i_xAxis_PCA, i_yAxis_PCA),
                pn.Row(pn.Pane(log_histogram, width=1000, config=plotly_config))
            )
            app_tabs.active = 1
            return visualization_map
        else:
            visualization_map = "Run analysis first!"
            return visualization_map
    except Exception:
        update_status = traceback.format_exc()
        return pn.Column(
            pn.Row(i_collapse_maps_PCA),
            pn.Row(i_xAxis_PCA, i_yAxis_PCA),             
            update_status
        )
            

@pn.depends(i_clusterwidget.param.value,i_mapwidget.param.value, cache_run.param.value)
def update_cluster_overview(clusterwidget, mapwidget, run):
    try:
        if run == True:
            list_genes = [goi for goi in i_class.genenames_sortedout_list if goi in i_class.markerproteins[clusterwidget]]
            i_class.cache_cluster_quantified = True
            distance_boxplot = i_class.distance_boxplot(cluster_of_interest=clusterwidget)
            if i_class.cache_cluster_quantified == False:
                return "This protein cluster was not quantified"
            
            else:
                df_quantification_overview = i_class.quantification_overview(cluster_of_interest=clusterwidget)
                profiles_plot = i_class.profiles_plot(map_of_interest = mapwidget, cluster_of_interest=clusterwidget)
                pca_plot = i_class.plot_cluster_pca(cluster_of_interest=clusterwidget)
                cluster_overview = pn.Column(
                        pn.Row(pn.Pane(pca_plot, width=500, config=plotly_config),
                               pn.Pane(profiles_plot, width=500, config=plotly_config),
                               pn.Pane(distance_boxplot, width=500, config=plotly_config),
                              ),
                        pn.Row(pn.Pane("In total {} proteins across all maps were quantified, whereas the following proteins were not consistently quantified throughout all maps: {}".format(
                                i_class.proteins_quantified_across_all_maps, ", ".join(list_genes)) if len(list_genes) != 0 else
                            "All genes from this cluster are quantified in all maps."), width=1000),
                        pn.Row(pn.widgets.DataFrame(df_quantification_overview, height=200, width=500, disabled=True))  
                        )
                return cluster_overview
        
        else:
            cluster_overview = "Run analysis first!"
            return cluster_overview
    except Exception:
        update_status = traceback.format_exc()
        return update_status
    
    
@pn.depends(i_clusterwidget.param.value, cache_run.param.value)
def update_cluster_details(clusterwidget, run):
    try:
        if run == True:
            cluster_details = i_class.distance_to_median_boxplot(cluster_of_interest = clusterwidget)
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
            fig_npg, fig_npr, fig_npr_dc, fig_npg_F, fig_npgf_F, fig_npg_F_dc = i_class.plot_quantity_profiles_proteinGroups()
            return pn.Column(
                    pn.Row(pn.Pane(fig_npg, config=plotly_config), pn.Pane(fig_npr, config=plotly_config), pn.Pane(fig_npr_dc, config=plotly_config)) ,
                pn.Row(pn.Pane(fig_npg_F, config=plotly_config), pn.Pane(fig_npgf_F, config=plotly_config), pn.Pane(fig_npg_F_dc, config=plotly_config)))
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
                pn.widgets.DataFrame(pd.read_json(i_class.analysed_datasets_dict[i_expname.value]["Overview table"]), height=200, width=600, disabled=True),

                i_logOR01_selection,
                df01_download_widget,
                pn.widgets.FileDownload(callback=json_download, filename="AnalysedDatasets.json")
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
    json.dump(
        domaps.SpatialDataSetComparison.analysed_datasets_dict, 
        sio, 
        indent=4, 
        sort_keys=True
    )
    sio.seek(0)
    return sio


@pn.depends(cache_run.param.value, i_logOR01_selection.param.value)
def df01_download_widget(run, logOR01_selection):
    if logOR01_selection == "0/1 normalized data":
        return pn.Column(pn.widgets.FileDownload(callback=df01_download, filename = "01_normalized_data.csv"), width=650) 
    if logOR01_selection == "log transformed data":
        return pn.Column(pn.widgets.FileDownload(callback=dflog_download, filename = "log_transformed_data.csv"), width=650)
    if logOR01_selection == "Overview - cluster distances":
        return pn.Column(pn.widgets.FileDownload(callback=table_download, filename = "cluster_distances.csv"), width=650)  
    else:
        return pn.Column(pn.widgets.FileDownload(callback=df_filteredRawData_download, filename = "stringency_filtered_raw_data.csv"), width=650)

    
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
    df = i_class.reframe_df_01ORlog_for_Perseus(i_class.df_stringencyFiltered)
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

analysis_tabs.clear()
analysis_tabs.append(("Data overview", update_data_overview))
analysis_tabs.append(("Depth and Coverage", update_quantity))
analysis_tabs.append(("Cluster Overview", update_cluster_overview))
analysis_tabs.append(("Cluster Details", update_cluster_details))
analysis_tabs.append(("Download", show_tabular_overview))


# In[ ]:


## Comparison tab

loading_status_comparison = pn.Row()
idle_comparison = pn.indicators.LoadingSpinner(value=False, width=100, height=100, color="primary")
loading_comparison = pn.indicators.LoadingSpinner(value=True, width=100, height=100, color="primary")

cache_uploaded_json = pn.widgets.Checkbox(value=False)
cache_run_json = pn.widgets.Checkbox(value=False)
#button_comparison = pn.widgets.Button(name="Compare experiments", width=50)
#i_jsonFile = pn.widgets.FileInput(name="Upload JSON file for comparison")
i_multi_choice = pn.widgets.CrossSelector(name="Select experiments for comparison", value=["a", "b"],
                                          options=["a", "b", "c"], definition_order=False)
i_ref_exp = pn.widgets.Select(name="Select experiments as reference", options=["a", "b", "c"])
i_collapse_cluster = pn.widgets.Checkbox(value=True, name="Collapse cluster")
comparison_status = pn.Pane("No datasets were compared yet")
i_markerset_or_cluster = pn.widgets.Checkbox(value=False, name="Display only protein clusters")
#i_ranking_boxPlot = pn.widgets.Checkbox(value=False, name="Display box plot")
i_ranking_boxPlot = pn.widgets.RadioBoxGroup(name="Types of ranking", options=["Box plot", "Bar plot - median", "Bar plot - sum"], inline=True)
#i_toggle_sumORmedian = pn.widgets.Toggle(name="Sum or Median", button_type="primary")
i_clusterwidget_comparison = pn.widgets.Select(options=[], name="Cluster of interest", width=300)
i_ExpOverview = pn.Row(pn.Pane("", width=1000))
i_include_dataset = pn.widgets.Checkbox(value=False, name="Include data analysed under 'Analysis' tab")
i_compare_gene = pn.widgets.TextInput(value="PLEC", name="Enter gene name or protein ID to see profile.")
i_compare_profile_style = pn.widgets.Select(options=[
    #"all profiles",
    "mean +- stdev", "mean +- SEM"])
i_compare_compartment = pn.widgets.MultiSelect(options=[], name="Select compartments for which to show summary profiles.")
#wdgts_comparison = [button_comparison]#,i_organism_comparison]#, i_include_dataset]
json_exp_name_cache = []


@pn.depends(i_multi_choice.param.value, i_clusterwidget_comparison.param.value, cache_run_json.param.value, i_xAxis_PCA_comp.param.value, i_yAxis_PCA_comp.param.value, 
            i_markerset_or_cluster.param.value)
def update_visualization_map_comparison(multi_choice, clusterwidget_comparison, run_json, xAxis_PCA_comp, yAxis_PCA_comp, markerset_or_cluster):
    try:
        if run_json == True:
            if multi_choice == []:
                return pn.Column(#pn.Row(i_multi_choice),
                                 pn.Row("Please select experiments for comparison")
                                )
            else:
                pass
            pca_global_comparison = i_class_comp.plot_global_pca_comparison(
                cluster_of_interest_comparison=clusterwidget_comparison,
                x_PCA=xAxis_PCA_comp, y_PCA=yAxis_PCA_comp,
                markerset_or_cluster=markerset_or_cluster, multi_choice=multi_choice
            )
            if markerset_or_cluster == False:
                return pn.Column(
                    pn.Row(i_clusterwidget_comparison),
                    pn.Row(i_markerset_or_cluster),
                    pn.Pane(pca_global_comparison, config=plotly_config),
                    pn.Row(i_xAxis_PCA_comp, i_yAxis_PCA_comp)
                )
            else:
                return pn.Column(
                    pn.Row(i_markerset_or_cluster),
                    pn.Pane(pca_global_comparison, config=plotly_config),
                    pn.Row(i_xAxis_PCA_comp, i_yAxis_PCA_comp)
                )
        else:
            pca_global_comparison = "Run analysis first!"
            return pca_global_comparison
    except Exception:
        update_status = traceback.format_exc()
        return pn.Column(
            pn.Row(i_markerset_or_cluster),
            pn.Row(i_xAxis_PCA_comp, i_yAxis_PCA_comp),
            update_status
        )

#### Biological precision tab
## Widgets
i_clusters_for_ranking = pn.widgets.CrossSelector(name="Select clusters to be considered for ranking calculation",
                                                  options=[], size=8)
i_minn_proteins = pn.widgets.IntSlider(name="Minimum number of proteins per complex", start=3, end=13, step=1, value=5)
i_collapse_maps = pn.widgets.Checkbox(value=False, name="Collapse maps")
i_reference_map = pn.widgets.Select(options=[], value="", name="Select reference map")

## Callbacks
@pn.depends(i_multi_choice.param.options, i_minn_proteins.param.value, cache_run_json.param.value)
def update_comp_cluster_coverage(exp_names, min_n, run_json):
    try:
        if not run_json:
            return ""
        [f,p,n] = i_class_comp.get_complex_coverage(min_n)
        i_clusters_for_ranking.options = [el for el in i_class_comp.markerproteins.keys() if el not in n.keys()]
        i_clusterwidget_comparison.options = [el for el in i_class_comp.markerproteins.keys() if el not in n.keys()]
        i_clusters_for_ranking.value = [el for el in i_class_comp.markerproteins.keys() if el in f.keys()]
        return pn.Row(
            "Coverage in all experiments \[>= n proteins]:<br>"+"<br>".join(["- {} ({})".format(k,v) for k,v in f.items()]),
            "Coverage in some experiments \[proteins/experiment]:<br>"+"<br>".join(["- {} \{}".format(k,str(v)) for k,v in p.items()]),
            "No sufficient coverage in any experiment \[proteins/experiment]:<br>"+"<br>".join(["- {} \{}".format(k,str(v)) for k,v in n.items()])
        )
    except Exception:
        update_status = traceback.format_exc()
        return update_status

@pn.depends(i_multi_choice.param.value, i_clusters_for_ranking.param.value, i_reference_map.param.value,
            i_minn_proteins.param.value, cache_run_json.param.value)
def update_comp_bp_global(multi_choice, clusters_for_ranking, reference_map, min_n, run_json):
    try:
        if not run_json:
            return ""
        if set(multi_choice) != set(i_class_comp.df_distance_comp.Experiment.values):
            i_class_comp.calc_biological_precision(multi_choice)
            i_reference_map.options = multi_choice
            if reference_map not in multi_choice:
                i_reference_map.value = multi_choice[0]
                reference_map = multi_choice[0]
        if clusters_for_ranking == []:
            return "Select at least one cluster"
        else:
            bp_bargraph, bp_boxplot_abs, bp_boxplot_rel = i_class_comp.plot_biological_precision(
                multi_choice, clusters_for_ranking, min_members=min_n, reference = reference_map)
            return pn.Column(
                pn.Row(i_reference_map),
                pn.Row(pn.Pane(bp_bargraph, config=plotly_config),
                       pn.Pane(bp_boxplot_abs, config=plotly_config),
                       pn.Pane(bp_boxplot_rel, config=plotly_config)))
        
    except Exception:
        update_status = traceback.format_exc()
        return pn.Column(
            pn.Row(i_reference_map),
            update_status
        )

@pn.depends(i_multi_choice.param.value, i_clusterwidget_comparison.param.value, i_collapse_maps.param.value,
            cache_run_json.param.value)
def update_comp_bp_single(multi_choice, clusterwidget_comparison, collapse_maps, run_json):
    try:
        i_class_comp.cache_cluster_quantified = True
        distance_comparison = i_class_comp.distance_boxplot_comparison(collapse_maps=collapse_maps, cluster_of_interest_comparison=clusterwidget_comparison, multi_choice=multi_choice)
        if i_class_comp.cache_cluster_quantified == False:
            return "Cluster was not quantified in any experiment"
        else:
            pca_comparison = i_class_comp.plot_pca_comparison(cluster_of_interest_comparison=clusterwidget_comparison, multi_choice=multi_choice)
            return pn.Row(pn.Pane(pca_comparison, config=plotly_config),
                          pn.Pane(distance_comparison, config=plotly_config))
    except Exception:
        update_status = traceback.format_exc()
        return update_status

## Tab assembly
comparison_tab_bp = pn.Column(
    pn.Row(pn.Column(i_clusters_for_ranking,i_minn_proteins),
           update_comp_cluster_coverage),
    update_comp_bp_global,
    pn.Row(i_clusterwidget_comparison,i_collapse_maps),
    update_comp_bp_single
)
    
@pn.depends(i_multi_choice.param.value, cache_run_json.param.value)
def update_npr_ngg_nprDc(multi_choice, run_json):
    try:
        if run_json == True: 
            if multi_choice == []:
                return pn.Column(
                                 pn.Row("Please select experiments for comparison"))
            else:
                fig_quantity_pg, fig_quantity_pr = i_class_comp.quantity_pr_pg_barplot_comparison(multi_choice=multi_choice)
                coverage_barplot = i_class_comp.coverage_comparison(multi_choice=multi_choice)
                return pn.Row(pn.Column(
                                 pn.Pane(fig_quantity_pg, config=plotly_config), 
                                 #pn.Pane(fig_quantity_pr, config=plotly_config), #removed for now, as the added information is not a lot
                                 pn.Pane(coverage_barplot, config=plotly_config)
                                ))
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
            if len(multi_choice)<=1:
                return pn.Column(
                    pn.Row(pn.Pane("Please select 2 or more experiments for comparison"), width=1000))
            else:
                venn_plot_total, venn_plot_int, figure_UpSetPlot_total, figure_UpSetPlot_int = i_class_comp.venn_sections(multi_choice_venn = multi_choice)
                return pn.Row(
                    pn.Column(
                        "Proteins quantified in at least one map",
                        pn.Pane(venn_plot_total),
                        pn.Row(figure_UpSetPlot_total,width=1000)
                    ),
                    pn.Column(
                        "Proteins quantified in all maps",
                        pn.Pane(venn_plot_int),
                        pn.Row(figure_UpSetPlot_int,width=1000)
                    )
                )
        else:
            venn_plot = "Run analysis first!"
            return venn_plot
    except Exception:
        update_status = traceback.format_exc()
        return update_status    


#### Dashboard structure
#### SMV analysis
########################
lo_benchmark_SVMs = pn.Column(
    pn.Card(header="###Add misclassification matrix", name="add_mcmatrix"),
    pn.Row(name="svm_output")
)

#### Layout elements
#### SVM analysis
####################
SVM_status = pn.pane.Markdown(width=400)
lo_SVM_heatmap = pn.Column(SVM_status)
i_SVMmatrix = pn.widgets.input.TextAreaInput(name="Misclassification matrix from Perseus", placeholder="Copy matrix here...")
i_SVMcomment = pn.widgets.TextInput(name="Comment", placeholder="Add comments here...")
i_SVMexp = pn.widgets.Select(name="Select experiments for the assignment of a misclassification matrix", options=["a", "b", "c"])
btn_SVM_addmatrix = pn.widgets.Button(name="Update misclassification matrix",
                                      button_type="success", width=400, height=50,
                                      css_classes=["button-main"])
#### Append layout to dashboard
#### SVM analysis
###############################
for el in [i_SVMexp, i_SVMcomment, i_SVMmatrix, lo_SVM_heatmap, btn_SVM_addmatrix]:
    lo_benchmark_SVMs.objects[[el.name for el in lo_benchmark_SVMs.objects].index("add_mcmatrix")].append(el)

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
    
    experiment, SVMcomment = i_SVMexp.value, i_SVMcomment.value
    
    # change comment of a stored dataset
    try:
        if i_SVMmatrix.value == "" and i_class_comp.svm_results[experiment]["default"]["misclassification"] is not None:
            SVMmatrix = i_class_comp.svm_results[experiment]["default"]["misclassification"]
        # return error if no misclassification  is uploaded or no comment is changed
        #if i_SVMmatrix.value == "":
        #    SVM_status.object = "No misclassification matrix is uploaded"
        else:
            SVMmatrix = pd.read_table(StringIO(i_SVMmatrix.value), sep="\t")
    except Exception:
        update_status = traceback.format_exc()
        SVM_status.object = update_status
# 
    # 2. Add to i_class_comp
    # defaults to name="default" and overwrite=True
    i_class_comp.add_svm_result(experiment, SVMmatrix, comment=SVMcomment)
    # 3. Add to mem_available_datasets so it can be downloaded together with the data
    mem_available_datasets[experiment]["SVM results"] = copy.deepcopy(i_class_comp.svm_results[experiment])
    for k in i_class_comp.svm_results[experiment].keys():
        mc_json = mem_available_datasets[experiment]["SVM results"][k]["misclassification"].to_json()
        mem_available_datasets[experiment]["SVM results"][k]["misclassification"] = mc_json
        pr_json = mem_available_datasets[experiment]["SVM results"][k]["prediction"].to_json()
        mem_available_datasets[experiment]["SVM results"][k]["prediction"] = pr_json
        
    #empty misclassification matrix of the textarea-input widget
    i_SVMmatrix.value = ""
    # 4. Trigger update of figure
    update_SVM_Analysis(i_multi_choice.value)


btn_SVM_addmatrix.on_click(add_SVM_result)
    
@pn.depends(i_SVMexp.param.value, i_SVMmatrix.param.value, watch=True)
def fill_svm_comment(SVMexp, SVMmatrix):
    """
    Acess stored comment if possible and no new matrix was uploaded.
    """
    try: 
        comment = i_class_comp.svm_results[SVMexp]["default"]["comment"]
        if SVMmatrix != "" or comment == "":
            timestampStr = datetime.now().strftime("%d-%b-%Y (%H:%M:%S)")
            comment = "Matrix added on {date}".format(date=timestampStr)
        else:
            pass
    except Exception:
            timestampStr = datetime.now().strftime("%d-%b-%Y (%H:%M:%S)")
            comment = "Matrix added on {date}".format(date=timestampStr)
    i_SVMcomment.value = comment
    #except Exception:
    #    update_status = traceback.format_exc()
    #    SVM_status.object = update_status
        
@pn.depends(i_SVMexp.param.value, i_SVMmatrix.param.value)
def show_misclassification(SVMexp, SVMmatrix):
    """
    Display heatmap of freshly uploaded or reloaded SVM misclassification matrix
    """
    if SVMmatrix == "":
        try: 
            df_SVM = i_class_comp.svm_results[SVMexp]["default"]["misclassification"]
            SVMheatmap = domaps.svm_heatmap(df_SVM)
        except: 
            SVMheatmap = "No misclassification matrix is uploaded"
    else: 
        df_SVM = pd.read_table(StringIO(SVMmatrix), sep="\t")
        SVMheatmap = domaps.svm_heatmap(df_SVM)
    return SVMheatmap

@pn.depends(i_SVMexp.param.value)
def load_matrix(SVMexp):
    """
    Load data if new experiment is selected 
    """
    try:
        SVMmatrix = i_class_comp.svm_results[SVMexp]["default"]["misclassification"]
        SVM_status.object = "Misclassification Matrix from class"
    except:
        SVM_status.object = "Upload missclassificationmatrix first"

@pn.depends(i_multi_choice.param.value, watch=True)        
def update_SVM_Analysis(multi_choice):     
    #empty lo if available
    lo_benchmark_SVMs.objects[[i.name for i in lo_benchmark_SVMs].index("svm_output")].objects = []
    try:
        if multi_choice == []:
            lo = pn.Column(pn.Row("Please select experiments for comparison"))
        else:
            for exp in multi_choice:
                # catch exception, if misclassification matrix is missing
                try:
                    i_class_comp.svm_results[exp]
                    i_class_comp.svm_processing()
                    fig_markerPredictionAccuracy, fig_clusterPerformance = i_class_comp.svm_plotting(multi_choice)
                    lo = pn.Column(pn.Pane(fig_markerPredictionAccuracy, config=plotly_config),
                                pn.Pane(fig_clusterPerformance, config=plotly_config),
                                pn.Row(pn.Pane("", width=1000))
                               )
                except KeyError:
                #if i_class_comp.svm_results[exp]["default"]["misclassification"] is None:
                    lo = pn.Column(pn.Row("Please upload misclassfication matrix for {}".format(exp)))
    except Exception:
        update_status = traceback.format_exc()
        lo = update_status
    lo_benchmark_SVMs.objects[[el.name for el in lo_benchmark_SVMs.objects].index("svm_output")].append(lo)
    
#### Callback output positioning
#### SVM analysis
################################
lo_SVM_heatmap.append(show_misclassification)
#lo_benchmark_SVMs.objects[[el.name for el in lo_benchmark_SVMs.objects].index("svm_output")].append(update_SVM_Analysis())

i_scatter_metric = pn.widgets.Select(name="Distance metric",
                                     options=["manhattan distance to average profile",
                                              "manhattan distance to median profile",
                                              "euclidean distance", "manhattan distance",
                                              "1 - cosine correlation", "1 - pearson correlation",])
i_scatter_consolidation = pn.widgets.Select(name="Consolidation of replicate distances",
                                            options=["average","median","sum"])
i_scatter_quantile = pn.widgets.FloatInput(name="Highlight quantile", value=0.5)
i_scatter_showrug = pn.widgets.Checkbox(name="Show rugplot in bottom margin", value=False)
i_scatter_showfull = pn.widgets.Checkbox(name="Show additional traces for full datasets", value=False)
i_scatter_x_cut = pn.widgets.FloatInput(name="Cut x-axis at", value=1)

@pn.depends(i_multi_choice.param.value, i_scatter_metric.param.value, i_scatter_consolidation.param.value,
            cache_run_json.param.value, i_scatter_showfull.param.value,
            i_scatter_quantile.param.value, i_scatter_showrug.param.value, i_scatter_x_cut.param.value)
def update_global_scatter_comparison(multi_choice, metric, consolidation, run_json,
                                     show_full, quantile, show_rug, x_cut):
    try:
        if run_json == False:
            return "Run analysis first!"
        if multi_choice == []:
            return pn.Column(pn.Row("Please select experiments for comparison"))
        if "reproducibility consolidation" not in i_class_comp.parameters.keys() or         i_class_comp.parameters["reproducibility metric"] != metric or         i_class_comp.parameters["reproducibility consolidation"] != consolidation:
            i_class_comp.calculate_global_scatter(metric, consolidation)
        reproducibility_plot = i_class_comp.plot_reproducibility_distribution(
            multi_choice=multi_choice,
            x_cut=x_cut, q=quantile,
            show_rug=show_rug, show_full=show_full)
        return pn.Column(pn.Row(pn.Column(i_scatter_metric,i_scatter_consolidation),
                                pn.Column("Customize:",i_scatter_showrug,i_scatter_showfull,
                                          i_scatter_quantile,i_scatter_x_cut)),
                         pn.Pane(reproducibility_plot, config=plotly_config)
                        )
    except Exception:
        update_status = traceback.format_exc()
        return pn.Column(pn.Row(pn.Column(i_scatter_metric,i_scatter_consolidation),
                                pn.Column("Customize:",i_scatter_showrug,i_scatter_showfull,
                                          i_scatter_quantile,i_scatter_x_cut)),
                         update_status
                        )

@pn.depends(i_multi_choice.param.value, i_compare_gene.param.value, i_compare_compartment.param.value,
            i_compare_profile_style.param.value, cache_run_json.param.value)
def update_profile_comparison(multi_choice, compare_gene, compare_compartment, compare_profile_style, run_json):
    if not run_json:
        return "Run analysis first!"
    try:
        try:
            plotdata = i_class_comp.df_01_filtered_combined.xs(compare_gene,level="Gene names",
                                                               axis=0, drop_level=False)\
                                   .stack("Fraction").reset_index().rename({0:"Profile [% total signal]"}, axis=1)
        except:
            try:
                plotdata = i_class_comp.df_01_filtered_combined.xs(compare_gene,level="Protein IDs",
                                                                   axis=0, drop_level=False)\
                                       .stack("Fraction").reset_index().rename({0:"Profile [% total signal]"}, axis=1)
            except:
                try:
                    plotdata = i_class_comp.df_01_filtered_combined.loc[[
                        compare_gene in " ".join(el) for el in i_class_comp.df_01_filtered_combined.index.values],:]\
                                           .stack("Fraction").reset_index().rename({0:"Profile [% total signal]"}, axis=1)
                except:
                    plotdata = pd.DataFrame()
        if len(plotdata) > 0:
            plotdata.drop("Exp_Map", axis=1, inplace=True)
            plotdata.sort_values("Fraction", key=domaps.natsort_list_keys, inplace=True)
            experiments = [el for el in multi_choice if el in plotdata["Experiment"].values]
            plotdata = plotdata.set_index("Experiment").loc[experiments,:].reset_index()
            plotprofile = px.line(plotdata, x="Fraction", y="Profile [% total signal]", line_group="Map",
                                  facet_col="Protein IDs", template="simple_white",
                                  line_dash = "Sequence" if "Sequence" in plotdata.columns else None,
                                  color="Experiment", hover_data=list(plotdata.columns))
        else:
            plotprofile = "No gene or protein ID matching {} found.".format(compare_gene)
        
        plotdata = pd.DataFrame()
        if compare_profile_style == "all profiles":
            for el in compare_compartment:
                el_df = i_class_comp.df_01_filtered_combined.xs(el, level="Compartment", axis=0, drop_level=False)                .stack("Fraction").reset_index().rename({0:"Profile [% total signal]"}, axis=1)
                plotdata = plotdata.append(el_df)
            if len(plotdata) > 0:
                plotdata.sort_values("Fraction", key=domaps.natsort_list_keys, inplace=True)
                plotdata = plotdata.set_index("Experiment").loc[multi_choice,:].reset_index()
                plotdata.insert(0, "PG_Map", [str(p)+"_"+str(m) for p,m in zip(plotdata["Protein IDs"], plotdata["Map"])])
                plotcompartments = px.box(plotdata, x="Fraction", y="Profile [% total signal]", color="Experiment",
                                          facet_row="Compartment", template="simple_white")
            else:
                plotcompartments = "Please select at least one compartment"
        else:
            for el in compare_compartment:
                el_df = i_class_comp.df_01_filtered_combined.xs(el, level="Compartment", axis=0, drop_level=False)
                plotdata = plotdata.append(el_df.stack("Fraction").groupby(["Compartment", "Map", "Experiment", "Fraction"])                    .apply(lambda x: pd.Series({"Profile [% total signal]": np.nanmean(x), "std":np.nanstd(x),
                                                "sem":np.nanstd(x)/np.sqrt(sum(np.isfinite(x)))}))\
                    .reset_index().rename(columns={"level_4": "measure", 0: "value"})\
                    .set_index(["Compartment", "Map", "Experiment", "Fraction", "measure"]).unstack("measure")\
                    .droplevel(0, axis=1).reset_index())
            if len(plotdata) > 0:
                plotdata.sort_values("Fraction", key=domaps.natsort_list_keys, inplace=True)
                plotdata = plotdata.set_index("Experiment").loc[multi_choice,:].reset_index()
                plotcompartments = px.line(plotdata, x="Fraction", y="Profile [% total signal]", color="Experiment",
                                           line_group="Map", line_dash="Compartment", template="simple_white",
                                           error_y="std" if "stdev" in compare_profile_style else "sem")
            else:
                plotcompartments = "Please select at least one compartment"
        
        return pn.Row(pn.Column(i_compare_gene,
                                pn.Pane(plotprofile, config=plotly_config)),
                      pn.Column(pn.Row(i_compare_compartment, i_compare_profile_style),
                                pn.Pane(plotcompartments, config=plotly_config)))
    except Exception:
        return pn.Row(pn.Column(i_compare_gene,
                                traceback.format_exc()),
                      pn.Column(pn.Row(i_compare_compartment, i_compare_profile_style)))


def update_multi_choice():
    i_multi_choice.options = i_class_comp.exp_names
    i_reference_map.options = i_class_comp.exp_names
    i_clusterwidget.options = list(i_class_comp.markerproteins.keys())
    i_clusters_for_ranking.options = list(i_class_comp.markerproteins.keys())
    i_clusters_for_ranking.value = list(i_class_comp.markerproteins.keys())
    i_multi_choice.value = i_class_comp.exp_names
    i_reference_map.value = i_class_comp.exp_names[0]

    
@pn.depends(i_multi_choice.param.value, watch=True)
def update_ref_exp(multi_choice):
    i_ref_exp.options = i_multi_choice.value
    #return i_ref_exp


@pn.depends(i_multi_choice.param.value, watch=True)
def update_ExpOverview(multi_choice):
    dict_analysis_parameters={}
    for exp_name in multi_choice:
        dict_analysis_parameters[exp_name] = i_class_comp.json_dict[exp_name]["Analysis parameters"]
    i_ExpOverview[0] = pn.widgets.DataFrame(pd.DataFrame.from_dict(dict_analysis_parameters))
    i_ExpOverview.value = pd.DataFrame.from_dict(dict_analysis_parameters)
    i_ExpOverview.disabled = True
    i_ExpOverview.height = 300


#### Dashboard structure
#### Download data
########################
lo_benchmark_download = pn.Column()

#### Layout elements
#### Download data
####################
i_benchmark_downloadselector = pn.widgets.Select(options=[
    "0-1 normalized data.csv",
    "pca coordinates.csv",
    "complex scatter.csv",
    "reproducibility.csv",
    "protein id alignment.csv",
    "benchmarking results collection.xlsx"
], value="0-1 normalized data.csv", width=200)
i_benchmark_download = pn.widgets.FileDownload(label="Download file", width=200, button_type="success")
o_benchmark_downloadpreview = pn.Row()

#### Append layout to dashboard
#### Download data
###############################
for el in [i_benchmark_downloadselector,i_benchmark_download,o_benchmark_downloadpreview]:
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
        out.columns = ["_".join(el) for el in out.columns.reorder_levels(["Experiment", "Map", "Fraction"]).values]
    elif sheet == "pca coordinates":
        out = i_class_comp.df_pca.copy()
        out.index = out.index.droplevel("merge type")
    elif sheet == "complex scatter":
        out = i_class_comp.df_distance_comp.copy()
        if mode=="csv":
            out["Cluster"] = ['"'+el+'"' for el in out["Cluster"]]
        out = out.set_index(
            ["Cluster", "Gene names", "Protein IDs", "Compartment", "Experiment", "Map"]).drop(
            ["Exp_Map", "merge type"], axis=1).unstack(["Experiment", "Map"]).copy()
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
        return pn.Column(f"{out.shape[0]} rows x {out.head().reset_index().shape[1]} columns",
                          pn.widgets.DataFrame(
            out.head(10),
            editable=False)
                         )
    else:
        return "This will download a .xlsx file with all tables as individual sheets."

@pn.depends(i_benchmark_downloadselector.param.value)
def benchmark_download(file_selection):
    i_benchmark_downloadselector.disabled=True
    i_benchmark_download.loading=True
    # set up list of files to format
    if file_selection.endswith(".csv"):
        sheets = [file_selection.split(".csv")[0]]
        mode="csv"
    elif file_selection.endswith(".xlsx"):
        sheets = [
            "0-1 normalized data",
            "pca coordinates",
            "complex scatter",
            "reproducibility",
            "protein id alignment",
        ]
        mode="xlsx"
    
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
        i_benchmark_download.loading=False
        i_benchmark_downloadselector.disabled=False
        return sio
    elif file_selection.endswith(".xlsx"):
        bio = BytesIO()
        excel = pd.ExcelWriter(bio, engine_kwargs = {"data_only": True})
        for sheet, df in out_dict.items():
            df.reset_index().T.reset_index().T.to_excel(
                excel, sheet_name=sheet,
                merge_cells=False, index=False, header=False)
        excel.save()
        bio.seek(0)
        i_benchmark_download.loading=False
        i_benchmark_downloadselector.disabled=False
        return bio
i_benchmark_download.callback = benchmark_download

#### Callback output positioning
#### Download data
################################
o_benchmark_downloadpreview.append(benchmark_download_preview)
    
#def execution_comparison(event):
#    if cache_uploaded_json.value == False:
#        comparison_status.object = "Please upload a JSON-file first"
#    else:        
#        #dashboard_comparison.objects[2:] = []
#        cache_run_json.value = False
#        for wdgt in wdgts_comparison:
#            wdgt.disabled = True
#        try:
#            loading_status_comparison.objects = [loading_comparison]
#            comparison_status.object = "Analysis in progress"
#            #protein_cluster = SpatialDataSet.markerproteins_set[i_organism_comparison.value].keys()
#            update_ref_exp(i_ref_exp)
#            global i_class_comp
#            comparison_status.object = "Initialization"
#            i_class_comp = domaps.SpatialDataSetComparison(ref_exp=i_ref_exp.value)#, clusters_for_ranking=protein_cluster, organism=i_organism_comparison.value)
#            i_class_comp.json_dict = json_dict
#            comparison_status.object = "Reading"
#            i_class_comp.read_jsonFile()
#            i_class_comp.calc_biological_precision()
#            i_class_comp.get_complex_coverage()
#            update_multi_choice(i_multi_choice, i_clusterwidget_comparison, i_clusters_for_ranking)
#            i_compare_compartment.options = list(set(
#                i_class_comp.df_01_filtered_combined.index.get_level_values("Compartment")))
#            comparison_status.object = "SVM Processing"
#            i_class_comp.svm_processing()
#            comparison_status.object = "PCA"
#            i_class_comp.perform_pca_comparison()
#            dashboard_comparison.append(pn.Row(i_multi_choice, i_ExpOverview))
#            dashboard_comparison.append(comparison_tabs)
#            loading_status_comparison.objects = []
#            comparison_status.object = "Comparison finished!"
#            cache_run_json.value = True
#        except Exception:
#            loading_status_comparison.objects = [""]
#            for wdgt in wdgts_comparison:
#                wdgt.disabled = False
#            comparison_status.object = traceback.format_exc()
#            cache_run_json.value = False
#button_comparison.on_click(execution_comparison)

@pn.depends(i_multi_choice.param.value, i_clusters_for_ranking.param.value,
            i_scatter_metric.param.value, i_scatter_consolidation.param.value, i_scatter_quantile.param.value,
            cache_run_json.param.value)
def update_benchmark_overview(multi_choice, clusters,
                              metric, consolidation, quantile,
                              run_json):
    try:
        if run_json == False:
            return "Run analysis first!"
        if multi_choice == []:
            return "Please select experiments for comparison"
        if "reproducibility consolidation" not in i_class_comp.parameters.keys() or         i_class_comp.parameters["reproducibility metric"] != metric or         i_class_comp.parameters["reproducibility consolidation"] != consolidation:
            i_class_comp.calculate_global_scatter(metric, consolidation)
        if len(clusters) == 0:
            return "Please select clusters for comparison in the intramap scatter tab"
        
        fig = i_class_comp.plot_overview(multi_choice, clusters, quantile)
        
        return pn.Column(
            pn.Pane(fig, config=plotly_config),
            "Note that this overview figure is affected by the settings you choose in the intra- and inter-map scatter tabs."
        )
    except:
        return traceback.format_exc()


comparison_tabs.clear()
comparison_tabs.append(("Overview", pn.Column(
    "Once you have run the analysis you can find different benchmarking outputs here and dive into the data.",
    update_benchmark_overview
)))
comparison_tabs.append(("PCA maps", update_visualization_map_comparison))
comparison_tabs.append(("Depth & Coverage", pn.Column(update_npr_ngg_nprDc, update_venn)))
comparison_tabs.append(("Intermap scatter", update_global_scatter_comparison))
comparison_tabs.append(("Intramap scatter", comparison_tab_bp))
comparison_tabs.append(("SVM Analysis", lo_benchmark_SVMs))
comparison_tabs.append(("Compare profiles", update_profile_comparison))
comparison_tabs.append(("Download data", lo_benchmark_download))


# In[ ]:


dashboard_benchmark.objects = [
    pn.Card(objects=[], header="## Manage data", name="manage_data", height_policy="fit"),
    pn.Row(objects=[], name="benchmark_output", height_policy="fit")
]

#### Manage data Card layout
# This accesses mem_available_datasets and mem_benchmark.
####

## Adding datasets row
i_upload_collection = pn.widgets.FileInput(name="Upload collection")
btn_load_reference = pn.widgets.Button(name="Load", width=100)
i_load_reference = pn.widgets.Select(options=pkg_resources.resource_listdir("domaps", "referencedata"),
                                     value=None, width=200)
lo_add_datasets = pn.Row(objects=[
    "**Upload collection from file (.json):**", i_upload_collection,
    "**Add reference set:**", i_load_reference, btn_load_reference
])

## Selection checkbox
i_dfs_available = pn.widgets.CheckBoxGroup(options=[], value=[],
                                           name="Datasets available")
lo_dfs_available = pn.WidgetBox(objects=["**Datasets available**", i_dfs_available], width=250)

## Management button group
i_coll_download = pn.widgets.FileDownload(label="Download selected as collection (.json)", filename="AnalysedDatasets.json", disabled=True)
btn_coll_editnames = pn.widgets.Button(name="Edit names and comments", disabled=True) # move from management
btn_coll_reannotate = pn.widgets.Button(name="Reannotate genes/organelles/complexes", disabled=True) # new functionality
btn_coll_runmain = pn.widgets.Button(name="Align and analyse selected datasets",
                                     button_type="success", disabled=True) # change from main comparison
btn_coll_dropmem = pn.widgets.Button(name="Drop selected datasets from memory",
                                     button_type="danger", disabled=True) # move from top of page
lo_coll_buttons = pn.WidgetBox(objects=[
    i_coll_download,
    btn_coll_editnames,
    btn_coll_reannotate,
    btn_coll_runmain,
    btn_coll_dropmem
], width=300)

## Interaction pane
lo_instructions_datamanagement = pn.Card(pn.pane.Markdown(textfragments["coll_status_default"]),
                                         header="**Explanation**", width=400)
o_status_datamanagement = pn.pane.Markdown(width=400)
def set_status_datamanagement(x, append=False):
    if not append:
        o_status_datamanagement.object = x+"<br><br>"
    else:
        o_status_datamanagement.object += x+"<br><br>"
    if x.startswith("Traceback") or o_status_datamanagement.object.count("<br><br>") > 1:
        resize(dashboard_benchmark.objects[[i.name for i in dashboard_benchmark].index("manage_data")])
    if DEBUG:
        time.sleep(0.3)
    o_status_datamanagement.object = x

set_status_datamanagement("Step 1: Add datasets")
o_dynamic_collectionmanagement = pn.Row()
lo_coll_interactions = pn.Column(objects=[lo_instructions_datamanagement,
                                          o_status_datamanagement,
                                          o_dynamic_collectionmanagement])

## Assemble collection management row
lo_manage_collection = pn.Row(objects=[lo_dfs_available, lo_coll_buttons, lo_coll_interactions])

#### Append elements to manage data row
dashboard_benchmark.objects[[i.name for i in dashboard_benchmark].index("manage_data")].objects = []
for el in [lo_add_datasets, lo_manage_collection]:
    dashboard_benchmark.objects[[i.name for i in dashboard_benchmark].index("manage_data")].append(el)

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
    #btn_coll_editnames,
    #btn_coll_reannotate,
    i_coll_download,
    btn_coll_runmain,
    btn_coll_dropmem,
    i_dfs_available,
    i_upload_collection,
    btn_load_reference
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
                json_loaded[exp+i_upload_collection.filename] = json_loaded.pop(exp)
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
        
        status += f"Loaded **{n_sets}** datasets from file **{i_upload_collection.filename}**"
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
    upload_collection(pkg_resources.resource_stream("domaps", "referencedata/"+i_load_reference.value).read())
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
        {k: mem_available_datasets[k] for k in dfs}, 
        sio,
        indent=4, 
        sort_keys=True
    )
    sio.seek(0)
    return sio
i_coll_download.callback = coll_downloadjson


def coll_runmain(event):
    """
    Drops selected datasets from collection stored in RAM.
    """
    # deactivate interface
    for el in lock_collection_change:
        el.disabled = True
    try:
        if btn_coll_runmain.button_type == "success":
            set_status_datamanagement("Aligning and analysing data ...")
            cache_run_json.value=False
            #### Main execution of the comparison
            loading_status_comparison.objects = [loading_comparison]
            selection = i_dfs_available.value
            global i_class_comp
            i_class_comp = domaps.SpatialDataSetComparison(ref_exp=selection[0])#, clusters_for_ranking=protein_cluster, organism=i_organism_comparison.value)
            i_class_comp.json_dict = {k: mem_available_datasets[k] for k in selection}
            set_status_datamanagement("Aligning data ...", append=True)
            i_class_comp.read_jsonFile()
            set_status_datamanagement("Analysing intra-map scatter ...", append=True)
            i_class_comp.calc_biological_precision()
            i_class_comp.get_complex_coverage()
            update_multi_choice()
            set_status_datamanagement("Running PCA ...", append=True)
            i_class_comp.perform_pca_comparison()
            i_compare_compartment.options = list(set(
                i_class_comp.df_01_filtered_combined.index.get_level_values("Compartment")))
            loading_status_comparison.objects = []
            cache_run_json.value=True
            set_status_datamanagement("Comparison finished!", append=True)
            
            #### Switch button mode
            btn_coll_runmain.button_type = "danger"
            btn_coll_runmain.name = "Reset analysis to make new selection"
            set_status_datamanagement(
                "Next step: Use interface below to evaluate and download benchmark results.", append=True)
        elif btn_coll_runmain.button_type == "danger":
            set_status_datamanagement("Resetting data analysis ...")
            cache_run_json.value=False
            
            btn_coll_runmain.button_type = "success"
            btn_coll_runmain.name = "Align and analyse selected datasets"
            set_status_datamanagement("Analysis results have been reset.")
        else:
            pass
        
    except Exception:
        set_status_datamanagement(traceback.format_exc())
        cache_run_json.value=False
        
    finally:
        # reactivate interface
        if btn_coll_runmain.button_type == "danger":
            for el in [
                i_coll_download,
                btn_coll_runmain,
                btn_coll_editnames
            ]:
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
        i_dfs_available.options = [el for el in i_dfs_available.options if el not in keys]
        i_dfs_available.value = [el for el in i_dfs_available.value if el not in keys]
        resize(lo_dfs_available)
        i_dfs_available.disabled = True
        
        set_status_datamanagement(f"Deleted **{len(keys)}** datasets **{', '.join(keys)}**")
        
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
            pn.Row(i_multi_choice, i_ExpOverview),
            comparison_tabs
        )
    else:
        return "Select data and run analysis."

#### Append callback output to benchmark output
dashboard_benchmark.objects[[i.name for i in dashboard_benchmark].index("benchmark_output")].objects = []
for el in [display_benchmark_output]:
    dashboard_benchmark.objects[[i.name for i in dashboard_benchmark].index("benchmark_output")].append(el)


# ## Management tab

# i_jsonFile_amendments_intended = pn.widgets.FileInput(name="Upload JSON file to be amended")
# i_json_ExpSelector = pn.widgets.CrossSelector(name="Select experiments, that will be removed from JSON file", width=1000)
# cache_uploaded_json_amendment = pn.widgets.Checkbox(value=False)
# cache_run_json_amendment = pn.widgets.Checkbox(value=False)
# button_reset = pn.widgets.Button(name="Reset", width=630)
# download_status = pn.Pane("Upload a JSON file first", width=1000)
# i_df_ExpComment = pn.widgets.DataFrame()
# wdgt_json = [button_reset]
# json_dict_amendments_intended = {}
# #make a cache, and say, if this hasnt been executed so far, please reset it
# dict_new_expNames = {}
# dict_new_comments = {}
# # i_exp_SVM = pn.widgets.Select(name="Select experiments as reference", options=["a", "b", "c"])
# button_SVM_analysis = pn.widgets.Button(name="Analyse misclassification matrix", width=50)
# # i_SVM_table = pn.widgets.input.TextAreaInput(name="Misclassification matrix from Perseus", placeholder="Copy matrix here...")
# cache_uploaded_SVM = pn.widgets.Checkbox(value=False)
# analysis_status_SVM = pn.Row(pn.Pane("No SVM analysis run yet", width=1000))
# i_all_or_marker = pn.widgets.Select(name="Select type of data for download", options=["Modified AnalysedDataset.json file", "0/1 normalized data, all experiments", 
#                                                                                "0/1 normalized data, markerset only, all experiments"], width=300)
#     
# @pn.depends(i_jsonFile_amendments_intended.param.value)#cache_run.param.value
# def open_jsonFile_amendment(jsonFile_amendments):#run
#     cache_run_json_amendment.value = False
#     if jsonFile_amendments is None:
#         cache_uploaded_json_amendment.value = False
#     else:
#         dashboard_manageDatasets.objects[2:] = []
#         #dashboard_manageDatasets.objects = []
#         dashboard_MissclassificationMatrix.objects = []
#         dashboard_amendment.objects = []
#         cache_uploaded_json_amendment.value = False
#         try:
#             
#             json_dict_cache = json.load(BytesIO(i_jsonFile_amendments_intended.value))
#             if hasattr(json_dict_cache, "keys") == False: 
#                     return "Your json-File does not fulfill the requirements"
#             else:
#                 global json_dict_amendments_intended
#                 try:
#                     json_dict_amendments_intended.update(json_dict_cache)
#                 except Exception:
#                     json_dict_amendments_intended = json_dict_cache
#                 i_json_ExpSelector.options = list(json_dict_amendments_intended.keys())
#                 cache_uploaded_json_amendment.value = True
#                 for wdgt in wdgt_json:
#                     wdgt.disabled = False
#                 download_status.object = "Upload successful! Select experiments now."
#                 dashboard_manageDatasets.append(amendment_tabs)
#                 dashboard_manageDatasets.append(button_reset)
#                 dashboard_amendment.append(pn.Column(i_df_ExpComment, download_status))
#                 analysis_status_SVM[0] = pn.Pane("Upload successful! Select experiments now.")
#                 dashboard_MissclassificationMatrix.append(analysis_status_SVM)
#                 return pn.Column(
#                                  i_json_ExpSelector,
#                                  i_all_or_marker,           
#                                  df01_json_download_widget,#pn.Column(pn.widgets.FileDownload(callback=df01_fromJson_download, filename = "all_01_normalized_data.csv", width=650)),
#                                  #button_reset,
#                                  )
#         except Exception: 
#             filereading_status_json = traceback.format_exc()
#             cache_uploaded_json_amendment.value = False
#             return filereading_status_json
# 
#     
# #@pn.depends(i_exp_SVM.param.value, watch=True)
# #def update_SVM_Matrix(exp_SVM):
# #    try:
# #        i_SVM_table.value = json_dict_amendments_intended[exp_SVM]["Misclassification Matrix"]
# #        
# #    except:
# #        i_SVM_table.value = ''
#     
# 
# # @pn.depends(i_json_ExpSelector.param.value, watch=True)
# # def update_exp_for_SVM(json_ExpSelector):
# #     if json_ExpSelector == []:
# #         dashboard_MissclassificationMatrix.objects = []
# #         analysis_status_SVM[0] = pn.Pane("Select experiments first")
# #         dashboard_MissclassificationMatrix.append(analysis_status_SVM)
# #     else:
# #         #i_SVM_table.value = ""       
# #         dashboard_MissclassificationMatrix.objects = []
# #         i_exp_SVM.options = json_ExpSelector
# #         analysis_status_SVM[0] = pn.Pane("Please paste a SVM Matrix first")
# #         dashboard_MissclassificationMatrix.append(pn.Row(i_exp_SVM, i_SVM_table))
# #         dashboard_MissclassificationMatrix.append(read_SVM_matrix)
# #         dashboard_MissclassificationMatrix.append(analysis_status_SVM)
#     
#         
# # @pn.depends(i_SVM_table.param.value)
# # def read_SVM_matrix(SVM_table):   
# #     if SVM_table == "":
# #         SVM_reading_status = "No misclassification matrix is uploaded"
# #         cache_uploaded_SVM.value = False
# #         analysis_status_SVM[0] = pn.Pane("Please paste a SVM Matrix first")
# #     else:
# #         cache_uploaded_SVM.value = False
# #         try:
# #             try:
# #                 df_SVM = pd.DataFrame(json.loads(SVM_table))
# #             except json.JSONDecodeError:
# #                 df_SVM = pd.read_table(StringIO(SVM_table), sep="\t")
# #             json_dict_amendments_intended[i_exp_SVM.value]["Misclassification Matrix"] = df_SVM.to_json()
# #             SVM_reading_status = domaps.svm_heatmap(df_SVM)
# #             cache_uploaded_SVM.value = True
# #             #button_SVM_analysis.disabled = False
# #             return pn.Column(#pn.Row(button_SVM_analysis),
# #                              pn.Row(SVM_reading_status, height=600),
# #                             )
# #         except Exception: 
# #             SVM_reading_status = "Paste the SVM matrix from Perseus only!"
# #             #SVM_reading_status = traceback.format_exc()
# #             cache_uploaded_SVM.value = False
# #             analysis_status_SVM[0] = pn.Pane(traceback.format_exc())    #pn.Pane("")
# #             return SVM_reading_status 
# 
# 
# @pn.depends(cache_run_json_amendment.param.value)
# def df01_fromJson_download(run):
#     if i_json_ExpSelector.value == []:
#         download_status.object = "No experiments are selected"
#         return
#     else:
#         df = domaps.reframe_df_01_fromJson_for_Perseus(json_dict_amendments_intended)
#         df = df.reset_index("Compartment")
#         df["Compartment"] = [el if el != "undefined" else "" for el in df.Compartment.values]
#         df = df.set_index("Compartment", append=True)
#         sio = StringIO()
#         df.to_csv(sio)
#         sio.seek(0)
#         return sio 
# 
#     
# @pn.depends(cache_run_json_amendment.param.value)
# def df01_marker_fromJson_download(run):
#     if i_json_ExpSelector.value == []:
#         download_status.object = "No experiments are selected"
#         return
#     else:
#         df = domaps.reframe_df_01_fromJson_for_Perseus(json_dict_amendments_intended)
#         df = df.loc[df.index.get_level_values("Compartment")!= "undefined"]
#         sio = StringIO()
#         df.to_csv(sio)
#         sio.seek(0)
#         return sio 
#     
#     
# @pn.depends(cache_run_json_amendment.param.value, i_all_or_marker.param.value)
# def df01_json_download_widget(run, all_or_marker):
#     if all_or_marker == "Modified AnalysedDataset.json file":
#         return pn.Column(pn.widgets.FileDownload(callback=json_amendment_download, filename = "AnalysedDatasets2.0.json"), width=650) 
#     elif all_or_marker == "0/1 normalized data, all experiments":
#         return pn.Column(pn.widgets.FileDownload(callback=df01_fromJson_download, filename = "all_01_normalized_data.csv"), width=650) 
#     else:
#         return pn.Column(pn.widgets.FileDownload(callback=df01_marker_fromJson_download, filename = "marker_01_normalized_data.csv"), width=650)
#     
#     
# def json_amendment_download():
#     if i_json_ExpSelector.value == []:
#         download_status.object = "No experiments are selected"
#         return
#     else:
#         json_new = json_dict_amendments_intended.copy()
#         exp_names_del = [elem for elem in i_json_ExpSelector.options if not elem in i_json_ExpSelector.value]        
#         for key in exp_names_del:
#             del json_new[key]
#         checked_exp = set()
#         redundant_expNames = set(new_Exp for new_Exp in dict_new_expNames.values() if new_Exp in checked_exp or checked_exp.add(new_Exp))
#         if redundant_expNames != set():
#             download_status.object = "Experiments are not allowed to be labelled identically"
#             return
#         else:
#             for exp_name in json_new:
#                 json_new[exp_name]["Analysis parameters"]["comment"] = dict_new_comments[exp_name] #dict_new_expNames
#             json_new = {dict_new_expNames[oldK]: value for oldK, value in json_new.items()}
#             sio = StringIO()
#             json.dump(
#                 json_new, 
#                 sio, 
#                 indent=4, 
#                 sort_keys=True
#             )
#             sio.seek(0)
#             download_status.object = "Download sucessful"
#             return sio
# 
#         
# def reset_json_amendment(event):
#     global json_dict_amendments_intended
#     json_dict_amendments_intended = {}
#     i_json_ExpSelector.options = []
#     i_json_ExpSelector.value = []
#     i_df_ExpComment.value = pd.DataFrame()
#     for wdgt in wdgt_json:
#         wdgt.disabled = True
#     download_status.object = "Reset sucessful"
# button_reset.on_click(reset_json_amendment)
# 
# 
# @pn.depends(i_json_ExpSelector.param.value, watch=True)
# def update_renameExp(json_ExpSelector):
#     dict_ExpComments = {}
#     for exp_name in json_ExpSelector:
#         dict_ExpComments[exp_name] = json_dict_amendments_intended[exp_name]["Analysis parameters"]["comment"]
#     df_ExpComments = pd.DataFrame(dict_ExpComments.items(), columns=["Experiment name - old", "Comment"])#pd.DataFrame.from_dict(dict_ExpComments)
#     df_ExpComments.insert(0, "Experiment name - new", df_ExpComments["Experiment name - old"])
#     df_ExpComments.set_index("Experiment name - old", inplace=True)
#     df_ExpComments.replace({"Experiment name - new": dict_new_expNames}, inplace=True)
#     exp_previous = list(dict_new_comments.keys())
#     for exp in exp_previous:
#         if exp not in json_ExpSelector:
#             del dict_new_comments[exp]
#     df_ExpComments.loc[dict_new_comments.keys(),"Comment"] = list(dict_new_comments.values())
#     i_df_ExpComment.value = df_ExpComments
#     i_df_ExpComment.height=len(json_ExpSelector)*50
#     #return i_df_ExpComment
# 
# @pn.depends(i_df_ExpComment.param.value, watch=True)
# def update_newExpNames(df_ExpComment):
#     try:        
#         global dict_new_expNames 
#         changed_expName = set(list(dict_new_expNames.values())+list(df_ExpComment["Experiment name - new"]))-set(dict_new_expNames.values())
#         dict_new_expNames = dict(zip(df_ExpComment.index, df_ExpComment["Experiment name - new"]))
#         global dict_new_comments
#         dict_new_comments_cache = dict_new_comments.copy()
#         changed_comment = set(list(dict_new_comments.values())+list(df_ExpComment["Comment"]))-set(dict_new_comments.values())
#         dict_new_comments = dict(zip(df_ExpComment.index, df_ExpComment["Comment"]))
#         if changed_expName == set() and changed_comment == set():
#             tracked_change = "No changes saved"
#         elif changed_expName != set() and changed_comment != set():
#             tracked_change = "Upload sucessful"
#         elif changed_expName != set():
#             new_exp = list(changed_expName)[0]
#             old_exp = list(dict_new_expNames.keys())[list(dict_new_expNames.values()).index(new_exp)]
#             tracked_change = "Experiment name was changed from {} to {}".format(old_exp, new_exp)
#         else:# changed_comment != set():
#             new_comment = list(changed_comment)[0]
#             exp = list(dict_new_comments.keys())[list(dict_new_comments.values()).index(new_comment)]
#             old_comment = dict_new_comments_cache[exp]
#             tracked_change = "Comment of the experiment {} was changed from {} to {}".format(exp, old_comment, new_comment)
#         download_status.object = tracked_change
#     except Exception:
#         pass
# 
# 
# def save_parameters():
#     parameters = dict()
#     parameters["multi_choice"] = i_multi_choice.value
#     
#     # PCA plot
#     parameters["xAxis_PCA_comp"] = i_xAxis_PCA_comp.value
#     parameters["yAxis_PCA_comp"] = i_yAxis_PCA_comp.value
#     
#     # biological precision
#     parameters["clusters_for_ranking"] = i_clusters_for_ranking.value
#     parameters["reference_map"] = i_reference_map.value
#     parameters["minn_proteins"] = i_minn_proteins.value
#     
#     # scatter
#     parameters["scatter_metric"] = i_scatter_metric.value
#     parameters["scatter_consolidation"] = i_scatter_consolidation.value
#     
#     # profile comparison
#     parameters["compare_compartment"] = i_compare_compartment.value
#     parameters["compare_profile_style"] = i_compare_profile_style.value
#     
#     sio = StringIO()
#     json.dump(parameters, sio, indent=4)
#     sio.seek(0)
#     return sio
# o_download_parameters = pn.widgets.FileDownload(callback=save_parameters, filename="parameters.json", width=200)

# In[ ]:


#In case of loading a json comparison larger than 80 MB
#with open("G:\_DIA manuscript\Figure panes and data\Figure 1\Figure1_v2_peptides.json", "br") as file:
#    i_jsonFile.value = file.read()


# In[ ]:


#In case of loading a dingle file larger than 80 MB
#i_file.filename = "filename.txt"
#with open("pathtofile.txt", "br") as file:
#    i_file.value = file.read()


# In[ ]:


#Set the order of the multi choice widget manually
#i_multi_choice.value=["21 min", "44 min", "100 min", "DDA"]

