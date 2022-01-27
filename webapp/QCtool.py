#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import natsort
import numpy as np
import pandas as pd
import param
import re
import traceback
import panel as pn
pn.extension("plotly")
import io
from io import BytesIO
from io import StringIO
from bokeh.models.widgets.tables import NumberFormatter
import json
import os
from importlib import reload
import pkg_resources

css = """
.detail_menu .bk-headers-wrapper{
  border-bottom: 2px solid #0a5da8 !important;
  margin-bottom: 10px;
  min-width: 1000px;
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
  min-width: 1000px;
}
.main_menu .bk-tab{
  color: #0a5da8;
}
.main_menu .bk-active{
  border-color: #0a5da8 !important;
  border-width: 2px 1px 0px 1px !important;
  color: #111111 !important;
}
.content-width{
  min-width: 1000px;
}

.bk-tabs-header{
  min-width: 1000px;
}

.card-title:first-child{
  font-size: 13px;
}
"""
pn.extension(raw_css=[css])

plotly_config={
      'toImageButtonOptions': {
            'format': 'svg', # one of png, svg, jpeg, webp
            'filename': 'figure'
      }
}


# In[ ]:


try:
    type(domaps)
    domaps = reload(domaps)
except Exception as ex:
    print(ex)
    import domaps

with open("textfragments.json", "r") as file:
    textfragments = json.load(file)

i_file = pn.widgets.FileInput(name="Upload file")

loading_status = pn.Row()
idle = pn.indicators.LoadingSpinner(value=False, width=100, height=100, color="primary")
loading = pn.indicators.LoadingSpinner(value=True, width=100, height=100, color="primary")

button_analysis = pn.widgets.Button(name="Analyse dataset", width=50)

i_logOR01_selection = pn.widgets.Select(options=["0/1 normalized data", "log transformed data",
                                                 "stringency filtered raw data", "Overview - cluster distances"],
                                        name="Select type of data for download", width=300)

i_acquisition = pn.widgets.Select(options=["LFQ6 - MQ", "LFQ5 - MQ", "LFQ6 - Spectronaut", "LFQ5 - Spectronaut",
                                           "SILAC - MQ", "Custom"], name="Acquisition", width=300)
i_organism = pn.widgets.Select(options=[el.split(".csv")[0] for el in
                                        pkg_resources.resource_listdir("domaps", "annotations/complexes")],
                               name="Organism", width=300, value="Homo sapiens - Uniprot")
i_reannotate_genes = pn.widgets.RadioButtonGroup(
    options=["skip", "from fasta column", "from Uniprot (requires an internet connection)"],
    name="Reannotate genes", value="skip")
i_custom_cids = pn.widgets.Select(name="column containing protein ids")
i_custom_cgenes = pn.widgets.Select(name="column containing gene names")
i_custom_cdata = pn.widgets.MultiSelect(name="columns containing data")
i_custom_normalized = pn.widgets.Checkbox(name="profiles are already 0-1 normalized (profile sum = 1)")
i_comment = pn.widgets.input.TextAreaInput(
    name="Additional Comments",
    placeholder="Write any kind of information associated with this dataset here...")

i_clusterwidget = pn.widgets.Select(options=["Proteasome", "Lysosome"], name="Cluster of interest", width=300)
i_mapwidget = pn.widgets.Select(options=["Map1", "Map2"], name="Map of interest", width=300)

i_collapse_maps_PCA = pn.widgets.Checkbox(value=False, name="Collapse maps")

cache_uploaded = pn.widgets.Checkbox(value=False)


cache_run = pn.widgets.Checkbox(value=False)

analysis_status = pn.Pane("", width=300)
filereading_status = pn.Pane("No data import yet", width=300)

i_expname = pn.widgets.TextInput(name="Experiment Name", placeholder="Enter your experiment name here here...")

i_consecutiveLFQi = pn.widgets.IntSlider(name="Consecutive LFQ intensities", start=1, end=10, step=1, value=4)
i_summed_MSMS_counts = pn.widgets.IntSlider(name="Summed MS/MS counts ≥ 2n x number of fractions",
                                            start=0, end=10, step=1, value=2)

i_RatioHLcount = pn.widgets.IntSlider(name="Quantification events (Ratio H/L Count)", start=1, end=10, step=1, value=2)
i_RatioVariability = pn.widgets.IntSlider(name="Ratio H/L variability [%]", start=0, end=100, step=5, value=30)

i_name_pattern = pn.widgets.Select(name="Name pattern",
                                   options=[".* (?P<rep>.*)_(?P<frac>.*)", ".* (?P<cond>.*)_(?P<rep>.*)_(?P<frac>.*)",
                                            ".* (?P<cond>.*)_(?P<frac>.*)_(?P<rep>.*)",
                                            ".* (?P<cond>.*)_(?P<frac>.*_.*)_(?P<rep>.*)",
                                            "(?P<rep>.*)_(?P<frac>.*)", "Custom"])
i_custom_namepattern = pn.widgets.TextInput(name="Customized Name Pattern",
                                            placeholder="Enter a string here...e.g.: .* (?P<rep>.*)_(?P<frac>.*)")
regex_pattern = {
    "(?P<rep>.*)_(?P<frac>.*)" : ["Spectronaut MAP1_03K"],
    ".* (?P<rep>.*)_(?P<frac>.*)" : ["MAP1_03K","MAP3_03K"],
    ".* (?P<cond>.*)_(?P<rep>.*)_(?P<frac>.*)" : ["EGF_rep1_06K","EGF_rep3_06K"],
    ".* (?P<cond>.*)_(?P<frac>.*)_(?P<rep>.*)" : ["Control_Mem_1", "Control_Cyt_1"],
    ".* (?P<cond>.*)_(?P<frac>.*_.*)_(?P<rep>.*)" : ["4h_mem_frac3_1", "co_mem_frac2_2"]
    }
i_pattern_examples = pn.widgets.Select(name = "Examples", options=regex_pattern[i_name_pattern.value])

@pn.depends(i_name_pattern.param.value, i_custom_namepattern, i_pattern_examples)
def custimization(name_pattern, custom_namepattern, pattern_examples):
    if name_pattern == "Custom":
        return i_custom_namepattern
    else:
        example_for_name_pattern = regex_pattern[name_pattern]
        i_pattern_examples.options = example_for_name_pattern
        return i_pattern_examples

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

@pn.depends(i_acquisition.param.value, i_consecutiveLFQi, i_summed_MSMS_counts)
def acquisition_response(acquisition, consecutiveLFQi, summed_MSMS_counts):
    if acquisition == "SILAC - MQ":
        return pn.Column(pn.Row(i_name_pattern, custimization),
                         pn.Row(pn.Pane("Stringency filtering")), pn.Row(i_RatioHLcount, i_RatioVariability))
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
                             #i_custom_normalized,
                             pn.Row(i_name_pattern, custimization))
        except Exception:
            return traceback.format_exc()
    else:
        return pn.Column(pn.Row(i_name_pattern, custimization),
                         pn.Row(pn.Pane("Stringency filtering")), pn.Row(i_consecutiveLFQi, i_summed_MSMS_counts))

#define widgets that should be disbled after run==True
wdgts = [i_acquisition,i_name_pattern,i_expname, i_pattern_examples, button_analysis, i_expname, i_organism,
         i_consecutiveLFQi, i_summed_MSMS_counts, i_RatioHLcount, i_RatioVariability, i_comment,
         i_custom_cids, i_custom_cgenes, i_custom_cdata, i_custom_normalized, i_reannotate_genes]            

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
                             pn.Row(acquisition_response), 
                             pn.Row(i_comment),
                             pn.Row(button_analysis))

        except Exception: 
            filereading_status = traceback.format_exc()
            cache_uploaded.value = False
            return filereading_status   
        
        
def execution(event):
    if cache_uploaded.value == False:
        analysis_status.object = "Please upload a file first"
    elif i_expname.value == "":
        analysis_status.object = "Please enter an experiment name first"
    else:        
        dashboard_analysis.objects = []
        cache_run.value = False
        for wdgt in wdgts:
            wdgt.disabled = True
        #if you did already your comparison, but add another experiment afterwards - without reloading your
        #AnylsedDatasets.json
        try:
            dashboard_analysis.append(i_clusterwidget)
            dashboard_analysis.append(i_mapwidget)
            dashboard_analysis.append(analysis_tabs)
            loading_status.objects = [loading]
            analysis_status.object = "Analysis in progress"
            if i_name_pattern.value == "Custom":
                namePattern = i_custom_namepattern.value
            else:
                namePattern = i_name_pattern.value
            global i_class
            if i_acquisition.value == "SILAC":
                i_class = domaps.SpatialDataSet(i_file.filename, i_expname.value, i_acquisition.value,
                                                comment=i_comment.value, name_pattern=namePattern,
                                                organism=i_organism.value, reannotate_genes=i_reannotate_genes.value,
                                                RatioHLcount=i_RatioHLcount.value,
                                                RatioVariability=i_RatioVariability.value)
            if i_acquisition.value == "Custom":
                i_class = domaps.SpatialDataSet(i_file.filename, i_expname.value, i_acquisition.value,
                                                comment=i_comment.value, name_pattern=namePattern,
                                                organism=i_organism.value, reannotate_genes=i_reannotate_genes.value,
                                                custom_columns = {"ids": i_custom_cids.value,
                                                                  "genes": i_custom_cgenes.value,
                                                                  "data": i_custom_cdata.value},
                                                custom_normalized = i_custom_normalized.value)
            else:
                i_class = domaps.SpatialDataSet(i_file.filename, i_expname.value, i_acquisition.value,
                                                comment=i_comment.value, name_pattern=namePattern,
                                                organism=i_organism.value, reannotate_genes=i_reannotate_genes.value,
                                                consecutiveLFQi=i_consecutiveLFQi.value,
                                                summed_MSMS_counts=i_summed_MSMS_counts.value)
            analysis_status.object = "Data Reading"
            i_class.data_reading(content=BytesIO(i_file.value))
            analysis_status.object = "Data Processing"
            i_class.processingdf()
            update_object_selector(i_mapwidget, i_clusterwidget)
            i_class.quantity_profiles_proteinGroups()
            analysis_status.object = "PCA"
            i_class.perform_pca()
            analysis_status.object = "Calculating Manhattan Distance"
            i_class.calc_biological_precision()
            analysis_status.object = "Calculating Dynamic Range"
            i_class.dynamic_range()
            analysis_status.object = "Writing Overview Table"
            i_class.results_overview_table()
            analysis_status.object = "Writing Analysed Dataset Dictionary"
            domaps.SpatialDataSetComparison.analysed_datasets_dict.update(i_class.analysed_datasets_dict)
            loading_status.objects = []
            analysis_status.object = "Analysis finished! Please open the 'Analysis' tab!"
            cache_run.value = True
        except:
            for wdgt in wdgts:
                wdgt.disabled = False
            loading_status.objects = [""]
            #The traceback gives no traceback, so out of that there will be still the output: Analysis in progress, although it is not possible. Out of that i removed the traceback
            analysis_status.object = traceback.format_exc()
            cache_run.value = False

button_analysis.on_click(execution)  

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
        return update_status
            

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
def update_dynamic_range(run):
    try:
        if run == True:
            fig_dynamic_range = i_class.plot_dynamic_range()
            return pn.Pane(fig_dynamic_range, config=plotly_config)
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

loading_status_comparison = pn.Row()
idle_comparison = pn.indicators.LoadingSpinner(value=False, width=100, height=100, color="primary")
loading_comparison = pn.indicators.LoadingSpinner(value=True, width=100, height=100, color="primary")

cache_uploaded_json = pn.widgets.Checkbox(value=False)
cache_run_json = pn.widgets.Checkbox(value=False)
button_comparison = pn.widgets.Button(name="Compare experiments", width=50)
i_jsonFile = pn.widgets.FileInput(name="Upload JSON file for comparison")
i_multi_choice = pn.widgets.CrossSelector(name="Select experiments for comparison", value=["a", "b"],
                                          options=["a", "b", "c"], definition_order=False)
i_scatter_metric = pn.widgets.Select(name="Distance metric",
                                     options=["manhattan distance to average profile",
                                              "manhattan distance to median profile",
                                              "euclidean distance", "manhattan distance",
                                              "1 - cosine correlation", "1 - pearson correlation",])
i_scatter_consolidation = pn.widgets.Select(name="Consolidation of replicate distances",
                                            options=["average","median","sum"])
i_ref_exp = pn.widgets.Select(name="Select experiments as reference", options=["a", "b", "c"])
i_collapse_cluster = pn.widgets.Checkbox(value=True, name="Collapse cluster")
dashboard_json = pn.Column("Please, upload a file first and press 'Compare clusters'", name="Comparison", css_classes=["content-width"])
comparison_status = pn.Pane("No datasets were compared yet")
i_markerset_or_cluster = pn.widgets.Checkbox(value=False, name="Display only protein clusters")
#i_ranking_boxPlot = pn.widgets.Checkbox(value=False, name="Display box plot")
i_ranking_boxPlot = pn.widgets.RadioBoxGroup(name="Types of ranking", options=["Box plot", "Bar plot - median", "Bar plot - sum"], inline=True)
#i_toggle_sumORmedian = pn.widgets.Toggle(name="Sum or Median", button_type="primary")
i_clusterwidget_comparison = pn.widgets.Select(options=[], name="Cluster of interest", width=300)
i_ExpOverview = pn.Row(pn.Pane("", width=1000))
i_include_dataset = pn.widgets.Checkbox(value=False, name="Include data analysed under 'Analysis' tab")
wdgts_comparison = [button_comparison]#,i_organism_comparison]#, i_include_dataset] 
json_dict = {}
json_exp_name_cache = []



@pn.depends(cache_run.param.value, i_jsonFile.param.value)#cache_run.param.value
def open_jsonFile(run, jsonFile):#run
    cache_run_json.value = False
    if run == False and jsonFile is None:
        filereading_status_json = "No file is uploaded"
        cache_uploaded_json.value = False
        return filereading_status_json
    else:
        try:
            if run == True:
                global json_dict
                json_dict.update(domaps.SpatialDataSet.analysed_datasets_dict)
            if jsonFile != None:
                global json_exp_name_cache
                if cache_uploaded_json.value==False:
                    json_dict.update(json.load(BytesIO(jsonFile)))
                    json_exp_name_cache = list(json.load(BytesIO(jsonFile)).keys())
                else:
                    for key in json_exp_name_cache:
                        del json_dict[key]
                    json_dict.update(json.load(BytesIO(jsonFile)))
                    json_exp_name_cache = list(json.load(BytesIO(jsonFile)).keys())
                    
            #    except:
            #        pass
            #elif jsonFile is not None:
            #    if i_include_dataset.value == False:
            #        json_dict = json.load(BytesIO(jsonFile))
            #    else:
            #        #global json_dict
            #        json_dict.update(json.load(BytesIO(jsonFile))) #i_class.
            try:
                dashboard_comparison.objects = dashboard_comparison.objects[0:4]
            except:
                pass
            if hasattr(json_dict, "keys") == False: #i_class.
                return "Your json-File does not fulfill the requirements"
            else:
                i_multi_choice.options = []
                filereading_status_json = list(json_dict.keys())# list(set(list(SpatialDataSet.analysed_datasets_dict.keys()) + )) #i_class.
                cache_uploaded_json.value = True
                for wdgt in wdgts_comparison:
                    wdgt.disabled = False
                return pn.Column(pn.Row("Experiments for comparison: {}".format(", ".join(filereading_status_json[:]))),
                                 #pn.Row(i_include_dataset),
                                 #pn.Row(i_organism_comparison),
                                 pn.Row(button_comparison),
                                 )
        
        except Exception: 
            filereading_status_json = traceback.format_exc()
            cache_uploaded_json.value = False
            return filereading_status_json                 

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
            pca_global_comparison = i_class_comp.plot_global_pca_comparison(cluster_of_interest_comparison=clusterwidget_comparison, x_PCA=xAxis_PCA_comp, y_PCA=yAxis_PCA_comp, 
                                                                       markerset_or_cluster=markerset_or_cluster, multi_choice=multi_choice)
            if markerset_or_cluster == False:
                return pn.Column(#pn.Row(i_multi_choice),
                                 pn.Row(i_clusterwidget_comparison),
                                 pn.Row(i_markerset_or_cluster),
                                 pn.Pane(pca_global_comparison, config=plotly_config),
                                 pn.Row(i_xAxis_PCA_comp, i_yAxis_PCA_comp)  
                                )
            else:
                return pn.Column(#pn.Row(i_multi_choice),
                                 pn.Row(i_markerset_or_cluster),
                                 pn.Pane(pca_global_comparison, config=plotly_config),
                                 pn.Row(i_xAxis_PCA_comp, i_yAxis_PCA_comp)   
                                )
        else:
            pca_global_comparison = "Run analysis first!"
            return pca_global_comparison
    except Exception:
        update_status = traceback.format_exc()
        return update_status

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
        return update_status

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

@pn.depends(i_multi_choice.param.value, cache_run_json.param.value)
def update_SVM_Analysis(multi_choice, run_json):
    try:
        if run_json == True: 
            if multi_choice == []:
                return pn.Column(
                                 pn.Row("Please select experiments for comparison"))
            else:
                if i_class_comp.cache_stored_SVM == False:
                    return pn.Column(
                                 pn.Row("No Missclassifiaction Matrix is stored"))
                else:
                    fig_markerPredictionAccuracy, fig_clusterPerformance = i_class_comp.svm_plotting(multi_choice)
                    return pn.Column(
                                 pn.Pane(fig_markerPredictionAccuracy, config=plotly_config),
                                 pn.Pane(fig_clusterPerformance, config=plotly_config), 
                                 pn.Row(pn.Pane("", width=1000)),
                    )
                
    except Exception:
        update_status = traceback.format_exc()
        return update_status 
                
@pn.depends(i_multi_choice.param.value, i_ref_exp.param.value, i_collapse_cluster.param.value, cache_run_json.param.value)
def update_dynamic_range_comparison(multi_choice, ref_exp, collapse_cluster, run_json):
    try:
        if run_json == True: 
            if multi_choice == []:
                return pn.Column(
                                 pn.Row("Please select experiments for comparison"))
            else:
                dynamic_range_barplot = i_class_comp.dynamic_range_comparison(collapse_cluster=collapse_cluster, multi_choice=multi_choice, ref_exp=ref_exp)
                return pn.Column(
                                 pn.Pane(dynamic_range_barplot, config=plotly_config),
                                 pn.Row(i_collapse_cluster),
                                 pn.Row(i_ref_exp)
                                )
        else:
            dynamic_range_barplot = "Run analysis first!"
            return dynamic_range_barplot
    except Exception:
        update_status = traceback.format_exc()
        return update_status

@pn.depends(i_multi_choice.param.value, i_scatter_metric.param.value, i_scatter_consolidation.param.value, cache_run_json.param.value)
def update_global_scatter_comparison(multi_choice, metric, consolidation, run_json):
    try:
        if run_json == True: 
            if multi_choice == []:
                return pn.Column(#pn.Row(i_multi_choice),
                                 pn.Row("Please select experiments for comparison"))
            else:
                scatter_histogram = i_class_comp.calculate_global_scatter(multi_choice=multi_choice,
                                                                              metric=metric, consolidation=consolidation)
                return pn.Column(pn.Row(i_scatter_metric),
                                 pn.Row(i_scatter_consolidation),
                                 pn.Pane(scatter_histogram, config=plotly_config)
                                )
        else:
            return "Run analysis first!"
    except Exception:
        update_status = traceback.format_exc()
        return update_status


def update_multi_choice(i_multi_choice, i_clusterwidget, i_clusters_for_ranking):
    i_multi_choice.options = list(json_dict.keys())
    i_reference_map.options = list(json_dict.keys())
    i_clusterwidget.options = list(i_class_comp.markerproteins.keys())
    i_clusters_for_ranking.options = list(i_class_comp.markerproteins.keys())
    i_clusters_for_ranking.value = list(i_class_comp.markerproteins.keys())
    i_multi_choice.value = list(json_dict.keys())
    i_reference_map.value = list(json_dict.keys())[0]

    
@pn.depends(i_multi_choice.param.value, watch=True)
def update_ref_exp(multi_choice):
    i_ref_exp.options = i_multi_choice.value
    #return i_ref_exp


@pn.depends(i_multi_choice.param.value, watch=True)
def update_ExpOverview(multi_choice):
    dict_analysis_parameters={}
    for exp_name in multi_choice:
        dict_analysis_parameters[exp_name] = json_dict[exp_name]["Analysis parameters"]
    i_ExpOverview[0] = pn.widgets.DataFrame(pd.DataFrame.from_dict(dict_analysis_parameters))
    i_ExpOverview.value = pd.DataFrame.from_dict(dict_analysis_parameters)
    i_ExpOverview.disabled = True
    i_ExpOverview.height = 300

    
def execution_comparison(event):
    if cache_uploaded_json.value == False:
        comparison_status.object = "Please upload a JSON-file first"
    else:        
        #dashboard_comparison.objects[2:] = []
        cache_run_json.value = False
        for wdgt in wdgts_comparison:
            wdgt.disabled = True
        try:
            loading_status_comparison.objects = [loading_comparison]
            comparison_status.object = "Analysis in progress"
            #protein_cluster = SpatialDataSet.markerproteins_set[i_organism_comparison.value].keys()
            update_ref_exp(i_ref_exp)
            global i_class_comp
            comparison_status.object = "Initialization"
            i_class_comp = domaps.SpatialDataSetComparison(ref_exp=i_ref_exp.value)#, clusters_for_ranking=protein_cluster, organism=i_organism_comparison.value)
            i_class_comp.json_dict = json_dict
            comparison_status.object = "Reading"
            i_class_comp.read_jsonFile()
            i_class_comp.calc_biological_precision()
            i_class_comp.get_complex_coverage()
            update_multi_choice(i_multi_choice, i_clusterwidget_comparison, i_clusters_for_ranking)
            comparison_status.object = "SVM Processing"
            i_class_comp.svm_processing()
            comparison_status.object = "PCA"
            i_class_comp.perform_pca_comparison()
            dashboard_comparison.append(pn.Row(i_multi_choice, i_ExpOverview))
            dashboard_comparison.append(comparison_tabs)
            loading_status_comparison.objects = []
            comparison_status.object = "Comparison finished!"
            cache_run_json.value = True
        except Exception:
            loading_status_comparison.objects = [""]
            for wdgt in wdgts_comparison:
                wdgt.disabled = False
            comparison_status.object = traceback.format_exc()
            cache_run_json.value = False
button_comparison.on_click(execution_comparison)


i_jsonFile_amendments_intended = pn.widgets.FileInput(name="Upload JSON file to be amended")
i_json_ExpSelector = pn.widgets.CrossSelector(name="Select experiments, that will be removed from JSON file", width=1000)
cache_uploaded_json_amendment = pn.widgets.Checkbox(value=False)
cache_run_json_amendment = pn.widgets.Checkbox(value=False)
button_reset = pn.widgets.Button(name="Reset", width=630)
download_status = pn.Pane("Upload a JSON file first", width=1000)
i_df_ExpComment = pn.widgets.DataFrame()
wdgt_json = [button_reset]
json_dict_amendments_intended = {}
#make a cache, and say, if this hasnt been executed so far, please reset it
dict_new_expNames = {}
dict_new_comments = {}
i_exp_SVM = pn.widgets.Select(name="Select experiments as reference", options=["a", "b", "c"])
button_SVM_analysis = pn.widgets.Button(name="Analyse misclassification matrix", width=50)
i_SVM_table = pn.widgets.input.TextAreaInput(name="Misclassification matrix from Perseus", placeholder="Copy matrix here...")
cache_uploaded_SVM = pn.widgets.Checkbox(value=False)
analysis_status_SVM = pn.Row(pn.Pane("No SVM analysis run yet", width=1000))
i_all_or_marker = pn.widgets.Select(name="Select type of data for download", options=["Modified AnalysedDataset.json file", "0/1 normalized data, all experiments", 
                                                                               "0/1 normalized data, markerset only, all experiments"], width=300)
    
@pn.depends(i_jsonFile_amendments_intended.param.value)#cache_run.param.value
def open_jsonFile_amendment(jsonFile_amendments):#run
    cache_run_json_amendment.value = False
    if jsonFile_amendments is None:
        cache_uploaded_json_amendment.value = False
    else:
        dashboard_manageDatasets.objects[2:] = []
        #dashboard_manageDatasets.objects = []
        dashboard_MissclassificationMatrix.objects = []
        dashboard_amendment.objects = []
        cache_uploaded_json_amendment.value = False
        try:
            
            json_dict_cache = json.load(BytesIO(i_jsonFile_amendments_intended.value))
            if hasattr(json_dict_cache, "keys") == False: 
                    return "Your json-File does not fulfill the requirements"
            else:
                global json_dict_amendments_intended
                try:
                    json_dict_amendments_intended.update(json_dict_cache)
                except Exception:
                    json_dict_amendments_intended = json_dict_cache
                i_json_ExpSelector.options = list(json_dict_amendments_intended.keys())
                cache_uploaded_json_amendment.value = True
                for wdgt in wdgt_json:
                    wdgt.disabled = False
                download_status.object = "Upload successful! Select experiments now."
                dashboard_manageDatasets.append(amendment_tabs)
                dashboard_manageDatasets.append(button_reset)
                dashboard_amendment.append(pn.Column(i_df_ExpComment, download_status))
                analysis_status_SVM[0] = pn.Pane("Upload successful! Select experiments now.")
                dashboard_MissclassificationMatrix.append(analysis_status_SVM)
                return pn.Column(
                                 i_json_ExpSelector,
                                 i_all_or_marker,           
                                 df01_json_download_widget,#pn.Column(pn.widgets.FileDownload(callback=df01_fromJson_download, filename = "all_01_normalized_data.csv", width=650)),
                                 #button_reset,
                                 )
        except Exception: 
            filereading_status_json = traceback.format_exc()
            cache_uploaded_json_amendment.value = False
            return filereading_status_json

    
@pn.depends(i_exp_SVM.param.value, watch=True)
def update_SVM_Matrix(exp_SVM):
    try:
        i_SVM_table.value = json_dict_amendments_intended[exp_SVM]["Misclassification Matrix"]
        
    except:
        i_SVM_table.value = ''
    

@pn.depends(i_json_ExpSelector.param.value, watch=True)
def update_exp_for_SVM(json_ExpSelector):
    if json_ExpSelector == []:
        dashboard_MissclassificationMatrix.objects = []
        analysis_status_SVM[0] = pn.Pane("Select experiments first")
        dashboard_MissclassificationMatrix.append(analysis_status_SVM)
    else:
        #i_SVM_table.value = ""       
        dashboard_MissclassificationMatrix.objects = []
        i_exp_SVM.options = json_ExpSelector
        analysis_status_SVM[0] = pn.Pane("Please paste a SVM Matrix first")
        dashboard_MissclassificationMatrix.append(pn.Row(i_exp_SVM, i_SVM_table))
        dashboard_MissclassificationMatrix.append(read_SVM_matrix)
        dashboard_MissclassificationMatrix.append(analysis_status_SVM)
    
        
@pn.depends(i_SVM_table.param.value)
def read_SVM_matrix(SVM_table):   
    if SVM_table == "":
        SVM_reading_status = "No misclassification matrix is uploaded"
        cache_uploaded_SVM.value = False
        analysis_status_SVM[0] = pn.Pane("Please paste a SVM Matrix first")
    else:
        cache_uploaded_SVM.value = False
        try:
            try:
                df_SVM = pd.DataFrame(json.loads(SVM_table))
            except json.JSONDecodeError:
                df_SVM = pd.read_table(StringIO(SVM_table), sep="\t")
            json_dict_amendments_intended[i_exp_SVM.value]["Misclassification Matrix"] = df_SVM.to_json()
            SVM_reading_status = svm_heatmap(df_SVM)
            cache_uploaded_SVM.value = True
            #button_SVM_analysis.disabled = False
            return pn.Column(#pn.Row(button_SVM_analysis),
                             pn.Row(SVM_reading_status, height=600),
                            )
        except Exception: 
            SVM_reading_status = "Paste the SVM matrix from Perseus only!"
            #SVM_reading_status = traceback.format_exc()
            cache_uploaded_SVM.value = False
            analysis_status_SVM[0] = pn.Pane(traceback.format_exc())    #pn.Pane("")
            return SVM_reading_status 


@pn.depends(cache_run_json_amendment.param.value)
def df01_fromJson_download(run):
    if i_json_ExpSelector.value == []:
        download_status.object = "No experiments are selected"
        return
    else:
        df = reframe_df_01_fromJson_for_Perseus(json_dict_amendments_intended)
        sio = StringIO()
        df.to_csv(sio)
        sio.seek(0)
        return sio 

    
@pn.depends(cache_run_json_amendment.param.value)
def df01_marker_fromJson_download(run):
    if i_json_ExpSelector.value == []:
        download_status.object = "No experiments are selected"
        return
    else:
        df = reframe_df_01_fromJson_for_Perseus(json_dict_amendments_intended)
        df = df.loc[df.index.get_level_values("Compartment")!= "undefined"]
        sio = StringIO()
        df.to_csv(sio)
        sio.seek(0)
        return sio 
    
    
@pn.depends(cache_run_json_amendment.param.value, i_all_or_marker.param.value)
def df01_json_download_widget(run, all_or_marker):
    if all_or_marker == "Modified AnalysedDataset.json file":
        return pn.Column(pn.widgets.FileDownload(callback=json_amendment_download, filename = "AnalysedDatasets2.0.json"), width=650) 
    if all_or_marker == "0/1 normalized data, all experiments":
        return pn.Column(pn.widgets.FileDownload(callback=df01_fromJson_download, filename = "all_01_normalized_data.csv"), width=650) 
    else:
        return pn.Column(pn.widgets.FileDownload(callback=df01_marker_fromJson_download, filename = "marker_01_normalized_data.csv"), width=650)
    
    
def json_amendment_download():
    if i_json_ExpSelector.value == []:
        download_status.object = "No experiments are selected"
        return
    else:
        json_new = json_dict_amendments_intended.copy()
        exp_names_del = [elem for elem in i_json_ExpSelector.options if not elem in i_json_ExpSelector.value]        
        for key in exp_names_del:
            del json_new[key]
        checked_exp = set()
        redundant_expNames = set(new_Exp for new_Exp in dict_new_expNames.values() if new_Exp in checked_exp or checked_exp.add(new_Exp))
        if redundant_expNames != set():
            download_status.object = "Experiments are not allowed to be labelled identically"
            return
        else:
            for exp_name in json_new:
                json_new[exp_name]["Analysis parameters"]["comment"] = dict_new_comments[exp_name] #dict_new_expNames
            json_new = {dict_new_expNames[oldK]: value for oldK, value in json_new.items()}
            sio = StringIO()
            json.dump(
                json_new, 
                sio, 
                indent=4, 
                sort_keys=True
            )
            sio.seek(0)
            download_status.object = "Download sucessful"
            return sio

        
def reset_json_amendment(event):
    global json_dict_amendments_intended
    json_dict_amendments_intended = {}
    i_json_ExpSelector.options = []
    i_json_ExpSelector.value = []
    i_df_ExpComment.value = pd.DataFrame()
    for wdgt in wdgt_json:
        wdgt.disabled = True
    download_status.object = "Reset sucessful"
button_reset.on_click(reset_json_amendment)


@pn.depends(i_json_ExpSelector.param.value, watch=True)
def update_renameExp(json_ExpSelector):
    dict_ExpComments = {}
    for exp_name in json_ExpSelector:
        dict_ExpComments[exp_name] = json_dict_amendments_intended[exp_name]["Analysis parameters"]["comment"]
    df_ExpComments = pd.DataFrame(dict_ExpComments.items(), columns=["Experiment name - old", "Comment"])#pd.DataFrame.from_dict(dict_ExpComments)
    df_ExpComments.insert(0, "Experiment name - new", df_ExpComments["Experiment name - old"])
    df_ExpComments.set_index("Experiment name - old", inplace=True)
    df_ExpComments.replace({"Experiment name - new": dict_new_expNames}, inplace=True)
    exp_previous = list(dict_new_comments.keys())
    for exp in exp_previous:
        if exp not in json_ExpSelector:
            del dict_new_comments[exp]
    df_ExpComments.loc[dict_new_comments.keys(),"Comment"] = list(dict_new_comments.values())
    i_df_ExpComment.value = df_ExpComments
    i_df_ExpComment.height=len(json_ExpSelector)*50
    #return i_df_ExpComment

@pn.depends(i_df_ExpComment.param.value, watch=True)
def update_newExpNames(df_ExpComment):
    try:        
        global dict_new_expNames 
        changed_expName = set(list(dict_new_expNames.values())+list(df_ExpComment["Experiment name - new"]))-set(dict_new_expNames.values())
        dict_new_expNames = dict(zip(df_ExpComment.index, df_ExpComment["Experiment name - new"]))
        global dict_new_comments
        dict_new_comments_cache = dict_new_comments.copy()
        changed_comment = set(list(dict_new_comments.values())+list(df_ExpComment["Comment"]))-set(dict_new_comments.values())
        dict_new_comments = dict(zip(df_ExpComment.index, df_ExpComment["Comment"]))
        if changed_expName == set() and changed_comment == set():
            tracked_change = "No changes saved"
        elif changed_expName != set() and changed_comment != set():
            tracked_change = "Upload sucessful"
        elif changed_expName != set():
            new_exp = list(changed_expName)[0]
            old_exp = list(dict_new_expNames.keys())[list(dict_new_expNames.values()).index(new_exp)]
            tracked_change = "Experiment name was changed from {} to {}".format(old_exp, new_exp)
        else:# changed_comment != set():
            new_comment = list(changed_comment)[0]
            exp = list(dict_new_comments.keys())[list(dict_new_comments.values()).index(new_comment)]
            old_comment = dict_new_comments_cache[exp]
            tracked_change = "Comment of the experiment {} was changed from {} to {}".format(exp, old_comment, new_comment)
        download_status.object = tracked_change
    except Exception:
        pass

dasboard_home = pn.Column(pn.Pane(textfragments["home_intro"], width=800), pn.Row(
    pn.Column(i_file, read_file,analysis_status, loading_status),
    pn.Accordion(pn.Pane(textfragments["upload_details"], name="Details on configuring your data"), width=500)), name="Home", css_classes=["content-width"])
dashboard_analysis = pn.Column("Please, upload a file first and press 'Analyse clusters'", name="Analysis", css_classes=["content-width"])
dashboard_MissclassificationMatrix = pn.Column("Please, upload a file first and press 'Analyse clusters'", name="SVM Analysis", css_classes=["content-width"])
dashboard_amendment = pn.Column("Please, upload a json file first", name="Renaming", css_classes=["content-width"])
dashboard_comparison = pn.Column(i_jsonFile, open_jsonFile, comparison_status, loading_status_comparison)
dashboard_manageDatasets = pn.Column(i_jsonFile_amendments_intended, open_jsonFile_amendment)
    
analysis_tabs = pn.Tabs(margin=10, css_classes=["content-width", "detail_menu"], dynamic=True)
analysis_tabs.append(("Data overview", update_data_overview))
analysis_tabs.append(("Depth and Coverage", update_quantity))
analysis_tabs.append(("Cluster Overview", update_cluster_overview))
analysis_tabs.append(("Cluster Details", update_cluster_details))
analysis_tabs.append(("Download", show_tabular_overview))
#analysis_tabs.append(("Dynamic Range", update_dynamic_range))

comparison_tabs = pn.Tabs(margin=10, css_classes=["content-width", "detail_menu"], dynamic=True, width=1000)
comparison_tabs.append(("Data Overview", update_visualization_map_comparison))
comparison_tabs.append(("Depth and Coverage", update_npr_ngg_nprDc))
comparison_tabs.append(("Unique and shared protein groups", update_venn))
comparison_tabs.append(("Global Scatter", update_global_scatter_comparison))
comparison_tabs.append(("Biological Precision", comparison_tab_bp))
comparison_tabs.append(("SVM Analysis", update_SVM_Analysis))
#comparison_tabs.append(("Dynamic Range", update_dynamic_range_comparison))

amendment_tabs = pn.Tabs(margin=10, css_classes=["content-width", "detail_menu"], dynamic=True)
amendment_tabs.append(("Change Experiment name and comment", dashboard_amendment))
amendment_tabs.append(("SVM Upload", dashboard_MissclassificationMatrix))

app_tabs = pn.Tabs(margin=10, css_classes=["content-width", "main_menu"], dynamic=True)
app_tabs.append(("Home", dasboard_home))
app_tabs.append(("Analysis", dashboard_analysis))
app_tabs.append(("Data comparison", dashboard_comparison))
app_tabs.append(("Manage Datasets", dashboard_manageDatasets))

app_tabs.append(("About", pn.Row(pn.Pane(textfragments["about_intro"], width=1000))))

#i_search = pn.widgets.TextInput(name="Search")
app_center = pn.Column(pn.Row(pn.Pane("# QC tool for Spatial Proteomics", width = 600),
                              pn.layout.HSpacer(),
                              #i_search,
                              #width=1600, 
                              margin=10),
                       app_tabs,
                       #pn.Spacer(background="#DDDDDD", height=100, margin=0)
                      )
app = pn.GridSpec()#sizing_mode="stretch_both", margin=0)
app[0,0] = pn.Spacer(background="white", margin=0) #"#DDDDDD"
app[0,9] = pn.Spacer(background="white", margin=0) #"#DDDDDD"
app[0,1:8] = app_center

pwd = pn.widgets.PasswordInput(name="Please enter password for access.")
app_container = pn.Column(pwd)

def check_pwd(event, app=app):
    pwd = event.new
    if pwd == "pwd":
        app_container[0]=app
pwd.param.watch(check_pwd, "value")

#try:
#    server.stop()
#except Exception:
#    print("First server startup")
#server = app.show(port=5066, websocket_max_message_size=2000000000)
app.servable()


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

