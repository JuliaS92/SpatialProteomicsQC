import panel as pn
pn.extension()
from panel.viewable import Viewer
import param
import pkg_resources
import traceback
from domaps import SpatialDataSet
import re
import pandas as pd
from io import BytesIO, StringIO
import json
import natsort

class ValidationError(Exception):
    pass

def help_icon(text="help text displayed on hover", sizing_mode="stretch_both", **kwargs):
    return pn.pane.HTML(f'<img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABUAAAAVCAYAAACpF6WWAAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAALoAAAC6ABUywAuAAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoAAANRSURBVDiNnZRNaF1VEMd/M+e+vCdJm7RNq0JNaZJ2U8QUIQuFriIWPyjFutGCuNGS6Eohbly50bgQJGkrrqUFK8VKUdOFS9FC1YgLm0+agiWJkMQ8SN69d8bFue/58lFBDwyHe878/zP3zPxHuM86MrbWZ+IvismAmT8CdAJLqjLvwW+YyOez59p+2QkrWw+6R1ePCjLiOac8dzAHJ5oUpoIEQQJXPbHhmXPtk/cl7RlbO+m5X/bM28kczx1yx5tIRYAQSUkESeQvTfTs1GDrtW2k3eerz5Lal55a8NTxzCFz3BysKbKCaIMQKQlS0lxK+vz0YOvXDdLu0dWjmNz0mu32muOp4alzojvhicMlDnYopUSYWcz57Id15lesmRBpEaRFV7ScPz71ese0Rmb50DPf7WlBWHO8Znx6dhcv9Vc4sEt5IMCrT1b47u0O+g8GvFb4FQl45u1eK40ASO+F6nGr2S3fMGzD8A1rAE4eKzF+a500Bdw5dCBh/N1OJu9lnPp4uZ4hUla0rEhZ3ctyPHG3M54XRcniW3oWM7j+4zqexmJhztzdlIm5lGNdCZ4aqII6EhxPHHIXzXkhEZMBzGKVrah2gzg+g6cO7lTKSs9Dgck/sljIUBCaI/X2M3lKzfwQEQPWtOfEAAVxxZ3zQ3vY3x5479JqvLctuGhdCbCn0dzUL4ooFvcTj5Z5/7W97N8beOWDJW7+voFUtPBrxoFDZwL8ifBwo2OF2OEioHH/6M19LC7nnH5jgYWqoxX9x2cTDhAWVVXmG0rRpj0UyilF4K9zKQtVj99Ropv9C1KB+cTVx1HpJwiigoeoFHJBLEpo9HqVOwtZzLAIJIWi6jhC/DNRvpUjY2t9eeY/7dSnnlrsgrx4Z4kZRjXpzn1akj6dHGr7WQJXt2g5AqIjt0cf5NJb+9CKNs6kpfArso54vTI72DahAJ7YsCSyskXLSEvM4PvplN/u5UilyKilbpuIl4OEdzZPqYvVZ9iwa/97SiXy3PRQ2zfb5mnvhbWnLfXLnnnHf5qnQV+eGmr9aschDdD7yXKP10ojbnbac5d/mfxO0C8SDcO3ByszzRzbSOvr8MW1x9T8DDkDQJc7nSIsAXdQbpjIldnBtomdsH8DjIPo5+eBi9cAAAAASUVORK5CYII=" title="{text}" />', sizing_mode=sizing_mode, **kwargs)


def wrapcallback_return(callback):
    def wrapper(*args, **kwargs):
        try:
            return callback(*args, **kwargs)
        except ValidationError as ex:
            return pn.pane.Markdown(str(ex), sizing_mode="stretch_width")
        except:
            return pn.pane.Markdown(traceback.format_exc(), sizing_mode="stretch_width")
    return wrapper


def get_settings(objects):
    settings = dict()
    for o in objects:
        if 'get_settings' not in dir(o):
            raise ValidationError(f"{o} has no method 'get_settings'.")
        s = o.get_settings()
        if s == None:
            continue
        for k,v in s.items():
            if k not in settings.keys():
                settings[k] = v
            else:
                if type(settings[k]) == list:
                    settings[k] += v
                elif type(settings[k]) == dict:
                    for vk in v.keys():
                        settings[k][vk] = v[vk]
    return settings


class FileInput(Viewer):
    
    value = param.ClassSelector(class_=bytes)
    filename = param.String(default="")
    accept = param.String(default=".*")
    width = param.Integer(default=200)
    disabled = param.Boolean(default=False)
    
    def __init__(self, **params):
        super().__init__(**params)
        self.layout = pn.Column(
            pn.Row(),
            pn.Row()
        )
        self._construct()
    
    def __panel__(self):
        return self.layout
    
    @param.depends('width', 'accept', watch=True)
    def _construct(self):
        try:
            file = pn.widgets.FileInput(width=self.width, accept=self.accept, disabled=self.disabled)
            self.layout[0] = file
        except Exception as ex:
            self.layout[1] = ex
        
        @pn.depends(file.param.filename, watch=True)
        def _sync_and_reconstruct(_):
            file = self.layout[0]
            if file.value != None:
                self.value = file.value
                self.filename = self.layout[0].filename
                self._construct()
    
    @param.depends('disabled', watch=True)
    def _disable(self):
        state = self.disabled
        self.layout[0].disabled = state


class SelectOrFile(Viewer):
    
    value = param.Selector(doc="Value selected from dropdown list", default=None)
    options = param.List(doc="List of options", default=["upload custom file"])
    custom_option = param.String(default="upload custom file")
    custom_file = param.String(default=None)
    accept = param.String(default="*")
    title = param.String(doc="Will be placed above selector widget", default="Select option or upload file")
    width = param.Integer(default=265)
    setting = param.String(doc="Name of the setting being selected", default="SelectOrUpload")
    disabled = param.Boolean(default=False)
    
    def __init__(self, **params):
        self._value = pn.widgets.Select()
        self._file = pn.widgets.FileInput()
        super().__init__(**params)
        self._layout = pn.Column(self._value, self._custom)
        self._sync_layout()
        self._sync_options()
        self._sync_value()
        
    def __panel__(self):
        return self._layout
    
    @param.depends('width', 'title', 'accept', watch=True)
    def _sync_layout(self):
        self._layout.width = self.width
        self._value.name = self.title
        self._file.accept = self.accept
    
    @param.depends('options', 'custom_option', watch=True)
    def _sync_options(self):
        value = self.value
        self._value.options = self.options if self.custom_option in self.options else self.options+[self.custom_option]
        if value not in self._value.options:
            self.custom_file = value
            self._file.filename = "loaded"
            self._value.value = self.custom_option
        else:
            self._value.value = value
    
    @param.depends('value', watch=True)
    def _sync_value(self):
        value = self.value
        if value not in self._value.options:
            self.custom_file = value
            self._file.filename = "loaded"
            self._value.value = self.custom_option
        else:
            self._value.value = value
    
    @param.depends('_value.value')
    def _custom(self):
        self.value = self._value.value
        if self._value.value == self.custom_option:
            return self._file
        else:
            return pn.Row()
    
    @param.depends('_file.value', watch=True)
    def _sync_custom_file(self):
        if self._file.value is not None:
            self.custom_file = self._file.value.decode("utf-8")
            self._file.value = None # delete reduntant copy from memory
    
    @param.depends('disabled', watch=True)
    def _disable(self):
        state = self.disabled
        self._value.disabled = state
        self._file.disabled = state
    
    def get_settings(self):
        if self.value == self.custom_option:
            return {self.setting: self.custom_file}
        else:
            return {self.setting: self.value}


class ConfigureSILACFilter(Viewer):
    
    toggle = param.Boolean(doc="Is this filter enabled?", default=True)
    value_count = param.Integer(doc="Minimal ratio count", default=2)
    value_variability = param.Integer(doc="Maximum variability", default=30)
    regex_count = param.String(doc="Regular expression uniquely identifying columns containing the ratio counts", default="Ratio H/L count .+")
    regex_variability = param.String(doc="Regular expression uniquely identifying columns containing the ratio variabilities", default="Ratio H/L variability.... .+")
    width = param.Integer(default=500)
    title = param.String(default="**SILAC filter** (count > x or count = x and variability <= y)")
    disabled = param.Boolean(default=False)
    regex_disabled = param.Boolean(default=False)
    
    def __init__(self, **params):
        self._toggle = pn.widgets.Checkbox(name="Apply filter")
        self._count = pn.widgets.IntSlider(
            name="Quantification events (Ratio H/L Count)",
            start=1, end=10, step=1
        )
        self._variability = pn.widgets.IntSlider(
            name="Ratio H/L variability [%]",
            start=0, end=100, step=1)
        self._regex_count = pn.widgets.TextInput(name="Regex/column name for ratio counts")
        self._regex_variability = pn.widgets.TextInput(name="Regex/column name for ratio variability")
        super().__init__(**params)
        self._layout = pn.Card(
            header=self.title,
            collapsible=False,
            objects=[pn.Row(self._toggle, help_icon("This is a filter intended for SILAC or other ratios, for which a variability and a count of quant events is available (e.g. MaxQuant output). Any ratios with a count > the specification are kept regardless of the variability, while ratios with exactly = the specified count are only kept if the variability is below the specified threshold.")),
                     pn.Row(self._count, self._variability),
                     pn.Row(self._regex_count, self._regex_variability)
                     ]
        )
        self._sync_layout()
        self._sync_widgets()
    
    def __panel__(self):
        return self._layout
    
    @param.depends('width', 'title', watch=True)
    def _sync_layout(self):
        self._layout.width = self.width
        self._layout.header = pn.pane.Markdown(self.title, width=self.width)
        self._toggle.width = self.width-20-45
        self._count.width = self.width//2-20
        self._variability.width = self.width//2-20
        self._regex_count.width = self.width//2-20
        self._regex_variability.width = self.width//2-20
    
    @param.depends('toggle',
                   'value_count', 'value_variability',
                   'regex_count', 'regex_variability',
                   watch=True)
    def _sync_widgets(self):
        self._toggle.value = self.toggle
        self._count.value = self.value_count
        self._variability.value = self.value_variability
        self._regex_count.value = self.regex_count
        self._regex_variability.value = self.regex_variability
    
    @param.depends('_toggle.value',
                   watch=True)
    def _sync_toggle(self):
        self.toggle = self._toggle.value
        status = not self._toggle.value
        self._count.disabled = status
        self._variability.disabled = status
        self._regex_count.disabled = any([status, self.regex_disabled])
        self._regex_variability.disabled = any([status, self.regex_disabled])
    
    @param.depends('_count.value', watch=True)
    def _sync_value_count(self):
        self.value_count = self._count.value
    
    @param.depends('_variability.value', watch=True)
    def _sync_value_variability(self):
        self.value_variability = self._variability.value
    
    @param.depends('_regex_count.value', watch=True)
    def _sync_regex_count(self):
        self.regex_count = self._regex_count.value
    
    @param.depends('_regex_variability.value',watch=True)
    def _sync_regex_variability(self):
        self.regex_variability = self._regex_variability.value
    
    @param.depends('disabled',
                   watch=True)
    def _disable(self):
        status = self.disabled
        self._toggle.disabled = status
        self._count.disabled = any([status, not self._toggle.value])
        self._variability.disabled = any([status, not self._toggle.value])
        self._regex_count.disabled = any([status, not self._toggle.value, self.regex_disabled])
        self._regex_variability.disabled = any([status, not self._toggle.value, self.regex_disabled])
    
    @param.depends('regex_disabled',
                   watch=True)
    def _disable_regex(self):
        status = any([self.regex_disabled, not self.toggle])
        self._regex_count.disabled = status
        self._regex_variability.disabled = status
    
    def get_settings(self):
        if self.toggle:
            settings = dict(
                quality_filter=["SILAC_countvar"],
                RatioCount=self.value_count,
                RatioVariability=self.value_variability,
                sets={
                    "Ratio count": self.regex_count,
                    "Ratio variability": self.regex_variability
                }
            )
            return settings
        return None


class ConfigureMSMScountFilter(Viewer):
    
    toggle = param.Boolean(doc="Is this filter enabled?", default=True)
    value_count = param.Integer(doc="Minimal average MS/MS count", default=2)
    regex_count = param.String(doc="Regular expression uniquely identifying columns containing the MS/MS counts", default="MS/MS count .+")
    width = param.Integer(default=240)
    title = param.String(default="**MS/MS count filter** (average >= x)")
    disabled = param.Boolean(default=False)
    regex_disabled = param.Boolean(default=False)
    
    def __init__(self, **params):
        self._toggle = pn.widgets.Checkbox(name="Apply filter")
        self._count = pn.widgets.IntSlider(
            name="Average MS/MS counts per profile",
            start=0, end=20, step=1
        )
        self._regex_count = pn.widgets.TextInput(name="Regex/column name for MS/MS counts")
        super().__init__(**params)
        self._layout = pn.Card(
            header=self.title,
            collapsible=False,
            objects=[pn.Row(self._toggle, help_icon("This filter will replace any profiles with NaN, if the average MS/MS count per profile point is < the specified value. The regular expression specifies which columns contain those counts. If for example peptide counts are preferred, simply change the regular expression.")),
                     pn.Row(self._count),
                     pn.Row(self._regex_count)
                     ]
        )
        self._sync_layout()
        self._sync_widgets()
    
    def __panel__(self):
        return self._layout
    
    @param.depends('width', 'title', watch=True)
    def _sync_layout(self):
        self._layout.width = self.width
        self._layout.header = pn.pane.Markdown(self.title, width=self.width)
        self._toggle.width = self.width-20-45
        self._count.width = self.width-20
        self._regex_count.width = self.width-20
        
    @param.depends('toggle', 'value_count', 'regex_count', watch=True)
    def _sync_widgets(self):
        self._toggle.value = self.toggle
        self._count.value = self.value_count
        self._regex_count.value = self.regex_count
    
    @param.depends('_toggle.value',
                   watch=True)
    def _sync_toggle(self):
        self.toggle = self._toggle.value
        status = not self._toggle.value
        self._count.disabled = status
        self._regex_count.disabled = any([status, self.regex_disabled])
    
    @param.depends('_count.value', watch=True)
    def _sync_value_count(self):
        self.value_count = self._count.value
    
    @param.depends('_regex_count.value', watch=True)
    def _sync_regex_count(self):
        self.regex_count = self._regex_count.value
    
    @param.depends('disabled',
                   watch=True)
    def _disable(self):
        status = self.disabled
        self._toggle.disabled = status
        self._count.disabled = any([status, not self._toggle.value])
        self._regex_count.disabled = any([status, not self._toggle.value, self.regex_disabled])
    
    @param.depends('regex_disabled',
                   watch=True)
    def _disable_regex(self):
        status = any([self.regex_disabled, not self.toggle])
        self._regex_count.disabled = status
    
    def get_settings(self):
        if self.toggle:
            settings = dict(
                quality_filter=["msms_count"],
                average_MSMS_counts=self.value_count,
                sets={
                    "MS/MS count": self.regex_count
                }
            )
            return settings
        return None


class ConfigureConsecutiveFilter(Viewer):
    
    toggle = param.Boolean(doc="Is this filter enabled?", default=True)
    value_count = param.Integer(doc="Minimal consecutive values", default=4)
    width = param.Integer(default=240)
    title = param.String(default="**Consecutive value filter**")
    disabled = param.Boolean(default=False)
    
    def __init__(self, **params):
        self._toggle = pn.widgets.Checkbox(name="Apply filter")
        self._count = pn.widgets.IntSlider(
            name="Minimum consecutive values per profile",
            start=0, end=20, step=1
        )
        super().__init__(**params)
        self._layout = pn.Card(
            header=self.title,
            collapsible=False,
            objects=[pn.Row(self._toggle, help_icon("As flexible filter for data completeness, specify how many consecutive non-zero values a profile must contain to be considered reasonably quantified. Any profiles with less will be removed and not included in any quality control measures.")),
                     pn.Row(self._count)
                    ]
        )
        self._sync_layout()
        self._sync_widgets()
    
    def __panel__(self):
        return self._layout
    
    @param.depends('width', 'title', watch=True)
    def _sync_layout(self):
        self._layout.width = self.width
        self._layout.header = pn.pane.Markdown(self.title, width=self.width)
        self._toggle.width = self.width-20-45
        self._count.width = self.width-20
        
    @param.depends('toggle',
                   'value_count', watch=True)
    def _sync_widgets(self):
        self._toggle.value = self.toggle
        self._count.value = self.value_count
    
    @param.depends('_toggle.value',
                   watch=True)
    def _sync_toggle(self):
        self.toggle = self._toggle.value
        status = not self._toggle.value
        self._count.disabled = status
    
    @param.depends('_count.value', watch=True)
    def _sync_value_count(self):
        self.value_count = self._count.value
    
    @param.depends('disabled',
                   watch=True)
    def _disable(self):
        status = self.disabled
        self._toggle.disabled = status
        self._count.disabled = any([status, not self._toggle.value])
    
    def get_settings(self):
        if self.toggle:
            settings = dict(
                quality_filter=["consecutive"],
                consecutive=self.value_count
            )
            return settings
        return None


class ConfigureFileContent(Viewer):
    
    source = param.Selector(default="MaxQuant")
    acquisition = param.Selector(default="LFQ")
    level = param.Selector(default="proteins")
    orientation = param.Selector(default="pivot")
    column_ids = param.String(default="")
    column_genes = param.String(default="")
    column_samples = param.String(default="")
    column_mainset = param.String(default="")
    name_pattern = param.String(default=".* (?P<cond>.*)_(?P<rep>.*)_(?P<frac>.*)")
    columns_annotation = param.List(default=[])
    columns = param.List(default=[])
    width = param.Integer(default=540)
    disabled = param.Boolean(default=False)
    
    def __init__(self, **params):
        self.file = FileInput()
        self.experiment_name = pn.widgets.TextInput(name="Experiment name")
        self._source = pn.widgets.Select(
            options=["MaxQuant", "Spectronaut", "custom"],
            name="Data source"
        )
        self._acquisition = pn.widgets.Select(
            options=["SILAC", "LFQ", "Intensity", "custom"],
            name="Acquisition mode"
        )
        self._level = pn.widgets.Select(
            options=["proteins", "peptides"], name="Data type"
        )
        self._orientation = pn.widgets.Select(
            options=["long", "pivot"], name="Data orientation"
        )
        self._column_custom_ids = pn.widgets.Select(name="Protein ids", value="")
        self._column_custom_genes = pn.widgets.Select(name="Gene symbols", value="")
        self._column_custom_samples = pn.widgets.Select(name="Sample names", value="")
        self._column_custom_mainset_long = pn.widgets.Select(name="Main dataset", value="")
        self._column_ids = pn.widgets.TextInput(name="Protein ids", disabled=True)
        self._column_genes = pn.widgets.TextInput(name="Gene symbols", disabled=True, value="")
        self._column_samples = pn.widgets.TextInput(name="Sample names", disabled=True, value="")
        self._column_mainset_long = pn.widgets.TextInput(name="Main dataset", disabled=True, value="")
        self._column_mainset_pivot = pn.widgets.TextInput(name="Main dataset regex", disabled=True, value="")
        self._columns_annotation = pn.widgets.MultiSelect(name="Additional columns to load", value=[])
        self._name_pattern = pn.widgets.TextInput(name="Sample name pattern", value="")
        self._pattern_presets = pn.widgets.Select(name="preset patterns",
                                   options=[".* (?P<rep>.*)_(?P<frac>.*)",
                                            ".* (?P<frac>.*)_(?P<rep>.*)",
                                            ".* (?P<cond>.*)_(?P<rep>.*)_(?P<frac>.*)",
                                            ".* (?P<cond>.*)_(?P<frac>.*)_(?P<rep>.*)",
                                            ".* (?P<cond>.*)_(?P<frac>.*_.*)_(?P<rep>.*)",
                                            "(?P<rep>.*)_(?P<frac>.*)"])
        self.fraction_ordering = OrderedMapper(name="Detected fractions")
        super().__init__(**params)
        self._layout = pn.Column(
            pn.Row(self.file, self.experiment_name),
            pn.Row(name="preview"),
            pn.Row(
                objects=[
                    pn.Column(pn.Row(self._source,
                                     self._acquisition,
                                     self._orientation,
                                     #self._level,
                                     help_icon("These selections determine, which options you have available for filtering and how the data is handled internally. For many combinations settings like column names and regular expressions are fixed by default. If you want to modify them you need to change the source to 'custom' instead. Orientation refers to how the data is arranged in the uploaded file.")
                                    )
                             )
                ],
                name="format"
            ),
            pn.Row(name="dependent configuration")
        )
        self._dependent_layout_long = pn.Column(
            pn.Card(
                objects=[
                    pn.Row(self._column_ids, self._column_genes, self._column_mainset_long),
                    pn.Row(self._name_pattern, self._pattern_presets, self._column_samples, help_icon("Correct specification of this regular expression is crucial to ensure correct processing of the profiles. Open the fractions card below to control that fractions are labelled and ordered correctly.")),
                    self._columns_annotation
                ],
                header=pn.pane.Markdown("**Column configuration**", width=self.width)),
            pn.Card(self.fraction_ordering, header=pn.pane.Markdown("**Fractions**"), collapsed=True)
        )
        self._dependent_layout_pivot = pn.Column(
            pn.Card(
                objects=[
                    pn.Row(self._column_ids, self._column_genes, self._column_mainset_pivot),
                    pn.Row(self._name_pattern, self._pattern_presets, help_icon("Correct specification of this regular expression is crucial to ensure correct processing of the profiles. Open the fractions card below to control that fractions are labelled and ordered correctly.")),
                    self._columns_annotation
                ],
                header=pn.pane.Markdown("**Column configuration**", width=self.width)),
            pn.Card(self.fraction_ordering, header=pn.pane.Markdown("**Fractions**", width=self.width), collapsed=True)
        )
        self._sync_width()
        self._sync_layout()
        self._sync_fileupload()
        self._sync_widgets_values()
        
    def __panel__(self):
        return self._layout
    
    @param.depends('width', watch=True)
    def _sync_width(self):
        self._layout.width = self.width
        self.file.width = self.width//2-20
        self.experiment_name.width = self.width//2-20
        self._source.width = self.width//4+20-11
        self._acquisition.width = self.width//4-20-11
        self._level.width = self.width//4-40-11
        self._orientation.width = self.width//4-40-11
        self._dependent_layout_long.width = self.width
        self._dependent_layout_pivot.width = self.width
        self._dependent_layout_long[0].width = self.width
        self._dependent_layout_pivot[0].width = self.width
        self._dependent_layout_long[1].width = self.width
        self._dependent_layout_pivot[1].width = self.width
        self._column_custom_ids.width = self.width//3-20
        self._column_custom_genes.width = self.width//3-20
        self._column_custom_mainset_long.width = self.width//3-20
        self._column_custom_samples.width = self.width//3-20
        self._column_ids.width = self.width//3-20
        self._column_genes.width = self.width//3-20
        self._column_mainset_long.width = self.width//3-20
        self._column_mainset_pivot.width = self.width//3-20
        self._column_samples.width = self.width//3-20-15
        self._name_pattern.width = self.width//3-20-15
        self._pattern_presets.width = self.width//3-20-15
        self._columns_annotation.width= self.width//2-20
    
    @param.depends('acquisition', 'level', 'orientation',
                   'source', watch=True)
    def _sync_layout(self):
        if self.name_pattern is None:
            self._pattern_presets.param.trigger('value')
        else:
            self._name_pattern.value = self.name_pattern
        set_source_default = False
        set_acquisition_default = False
        if self._source.value != self.source:
            self._source.value = self.source
            set_source_default = True
        if self._level.value != self.level:
            self._level.value = self.level
            set_source_default = True
        if self._orientation.value != self.orientation:
            self._orientation.value = self.orientation
            set_source_default = True
        if self._acquisition.value != self.acquisition:
            self._acquisition.value = self.acquisition
            set_acquisition_default = True
        if set_source_default:
            self._source.param.trigger('value')
        elif set_acquisition_default:
            self._acquisition.param.trigger('value')
    
    def _sync_dependent_layout(self):
        self._sync_pattern()
        if self.file.value is not None:
            self._layout[3] = self._dependent_layout_long if self.orientation == "long"\
            else self._dependent_layout_pivot
        
    @param.depends('file.value', watch=True)
    def _sync_fileupload(self):
        if self.file.value == None:
            self._layout[1] = pn.Row()
            self._layout[3] = pn.Row()
        else:
            df_head = pd.read_csv(BytesIO(self.file.value), nrows=5, sep="\t")
            if len(df_head.columns) < 2:
                df_head = pd.read_csv(BytesIO(self.file.value), nrows=5)
            if len(df_head.columns) < 2:
                df_head = pd.read_excel(BytesIO(self.file.value), nrows=5)                
            self._layout[1] = pn.Column("Data preview:", pn.widgets.DataFrame(df_head, width=self.width-10, disabled=True))
            self.columns = [el for el in df_head.columns]
            self._layout[3] = self._dependent_layout_long if self.orientation == "long"\
                              else self._dependent_layout_pivot
            self._set_defaults()
    
    @param.depends('columns', watch=True)
    def _sync_availablecolumns(self):
        self._column_custom_ids.options = self.columns
        self._column_custom_genes.options = self.columns
        self._column_custom_mainset_long.options = self.columns
        self._column_custom_samples.options = self.columns
        self._columns_annotation.options = self.columns
    
    @param.depends('_source.value', watch=True)
    def _sync_source(self):
        self.source = self._source.value
        if self.source.startswith("MaxQuant"):
            self.file.accept = ".txt"
            self._orientation.value = "pivot"
            self._orientation.disabled = True
        else:
            self.file.accept = ".tsv,.txt,.csv,.xls"
            self._orientation.disabled = False
        self._sync_fileupload()
    
    @param.depends('_orientation.value', watch=True)
    def _sync_orientation(self):
        self.orientation = self._orientation.value
        self._set_defaults()
    
    @param.depends('_level.value', watch=True)
    def _sync_level(self):
        self.level = self._level.value
        self._set_defaults()
    
    @param.depends('_acquisition.value', watch=True)
    def _sync_acquisition(self):
        self.acquisition = self._acquisition.value
        self._set_defaults_acquisition()
    
    def _set_defaults(self):
        source = "_".join([self.source, self.level, self.orientation])
        if source in SpatialDataSet.defaultsettings["sources"].keys():
            source = SpatialDataSet.defaultsettings["sources"][source]
            try:
                self._column_ids.value = source["original_protein_ids"]
                self._dependent_layout_long[0][0][0] = self._column_ids
                self._dependent_layout_pivot[0][0][0] = self._column_ids
            except:
                self._dependent_layout_long[0][0][0] = self._column_custom_ids
                self._dependent_layout_pivot[0][0][0] = self._column_custom_ids
                self._column_custom_ids.param.trigger('value')
            try:
                self._column_genes.value = source["genes"]
                self._dependent_layout_long[0][0][1] = self._column_genes
                self._dependent_layout_pivot[0][0][1] = self._column_genes
            except:
                self._dependent_layout_long[0][0][1] = self._column_custom_genes
                self._dependent_layout_pivot[0][0][1] = self._column_custom_genes
                self._column_custom_genes.param.trigger('value')
            if self.orientation == "long":
                try:
                    self._column_samples.value = source["samples"]
                    self._dependent_layout_long[0][1][2] = self._column_samples
                except:
                    self._dependent_layout_long[0][1][2] = self._column_custom_samples
                    self._column_custom_samples.param.trigger('value')
            try:
                self._columns_annotation.value = source["annotation_columns"]
                self._dependent_layout_long[0][2] = self._columns_annotation
                self._dependent_layout_pivot[0][2] = self._columns_annotation
            except:
                self._dependent_layout_long[0][2] = self._columns_annotation
                self._dependent_layout_pivot[0][2] = self._columns_annotation
        
        else:
            self._dependent_layout_long[0][0][0] = self._column_custom_ids
            self._dependent_layout_pivot[0][0][0] = self._column_custom_ids
            self._column_custom_ids.param.trigger('value')
            self._dependent_layout_long[0][0][1] = self._column_custom_genes
            self._dependent_layout_pivot[0][0][1] = self._column_custom_genes
            self._column_custom_genes.param.trigger('value')
            self._column_mainset_pivot.disabled = False
            self._dependent_layout_long[0][1][2] = self._column_custom_samples
            self._column_custom_samples.param.trigger('value')
            self._dependent_layout_long[0][2] = self._columns_annotation
            self._dependent_layout_pivot[0][2] = self._columns_annotation
        
        self._set_defaults_acquisition()
    
    
    def _set_defaults_acquisition(self):
        source = "_".join([self.source, self.level, self.orientation])
        if source in SpatialDataSet.defaultsettings["sources"].keys():
            source = SpatialDataSet.defaultsettings["sources"][source]
            if self.orientation == "pivot":
                try:
                    self._column_mainset_pivot.value = source["sets"]\
                    [SpatialDataSet.defaultsettings["acquisition_modes"][self.acquisition]["sets"][0]]
                    self._column_mainset_pivot.disabled = True
                except:
                    self._column_mainset_pivot.disabled = False
                    self._column_mainset_pivot.param.trigger('value')
            elif self.orientation == "long":
                try:
                    self._column_mainset_long.value = source["sets"]\
                    [SpatialDataSet.defaultsettings["acquisition_modes"][self.acquisition]["sets"][0]]
                    self._dependent_layout_long[0][0][2] = self._column_mainset_long
                except:
                    self._dependent_layout_long[0][0][2] = self._column_custom_mainset_long
                    self._column_custom_mainset_long.param.trigger('value')
        else:
            if self.orientation == "pivot":
                self._column_mainset_pivot.disabled = False
                self._column_mainset_pivot.param.trigger('value')
            elif self.orientation == "long":
                self._dependent_layout_long[0][0][2] = self._column_custom_mainset_long
                self._column_custom_mainset_long.param.trigger('value')
        
        self._sync_dependent_layout()
    
    @param.depends('columns_annotation', 'column_ids', 'column_genes', 'column_samples', 'column_mainset', 'name_pattern', watch=True)
    def _sync_widgets_values(self):
        if self._columns_annotation.value != self.columns_annotation:
            self._columns_annotation.value = self.columns_annotation
        try:
            self._column_custom_ids.value = self.column_ids
        except:
            pass
        try:
            self._column_custom_genes.value = self.column_genes
        except:
            pass
        try:
            self._column_custom_samples.value = self.column_samples
        except:
            pass
        if self.orientation == "long":
            try:
                self._column_custom_mainset_long.value = self.column_mainset
            except:
                pass
        else:
            self._column_mainset_pivot.value = self.column_mainset
        self._name_pattern.value = self.name_pattern
        
    @param.depends('_column_ids.value', watch=True)
    def _sync_column_ids(self):
        self.column_ids = self._column_ids.value
        self._column_ids.background = 'salmon' if self.column_ids not in self.columns else 'lightgreen'
    
    @param.depends('_column_custom_ids.value', watch=True)
    def _sync_column_custom_ids(self):
        self._column_ids.value = self._column_custom_ids.value
        
    @param.depends('_column_genes.value', watch=True)
    def _sync_column_genes(self):
        self.column_genes = self._column_genes.value
        self._column_genes.background = 'salmon' if self.column_genes not in self.columns else 'lightgreen'
    
    @param.depends('_column_custom_genes.value', watch=True)
    def _sync_column_custom_genes(self):
        self._column_genes.value = self._column_custom_genes.value
        
    @param.depends('_column_mainset_long.value', watch=True)
    def _sync_column_mainset_long(self):
        if self.orientation == "long":
            self.column_mainset = self._column_mainset_long.value
        self._column_mainset_long.background = 'salmon' if self.column_mainset not in self.columns else 'lightgreen'
    
    @param.depends('_column_custom_mainset_long.value', watch=True)
    def _sync_column_custom_mainset_long(self):
        self._column_mainset_long.value = self._column_custom_mainset_long.value
    
    @param.depends('_column_mainset_pivot.value', watch=True)
    def _sync_column_mainset_pivot(self):
        if self.orientation == "pivot":
            self.column_mainset = self._column_mainset_pivot.value
    
    @param.depends('_column_samples.value', watch=True)
    def _sync_column_samples(self):
        self.column_samples = self._column_samples.value
        self._column_samples.background = 'salmon' if self.column_samples not in self.columns else 'lightgreen'
    
    @param.depends('_column_custom_samples.value', watch=True)
    def _sync_column_custom_samples(self):
        self._column_samples.value = self._column_custom_samples.value
    
    @param.depends('_columns_annotation.value', watch=True)
    def _sync_columns_annotation(self):
        self.columns_annotation = self._columns_annotation.value
    
    @param.depends('experiment_name.value', watch=True)
    def _sync_name_color(self):
        self.experiment_name.background = 'salmon' if self.experiment_name.value == "" else None
    
    @param.depends('_pattern_presets.value', watch=True)
    def _sync_preset(self):
        self._name_pattern.value = self._pattern_presets.value
    
    @param.depends('_name_pattern.value', 'column_mainset', watch=True)
    def _sync_pattern(self):
        self.name_pattern = self._name_pattern.value
        if self.file.value is not None:
            try:
                if "(?P<frac>" not in self.name_pattern:
                    raise ValidationError("Name pattern must contain a capture group name 'frac' (see presets).")
                if self.orientation == "long":
                    if len(pd.read_csv(BytesIO(self.file.value), nrows=5, sep="\t").columns) > 2:
                        samples = pd.read_csv(BytesIO(self.file.value),
                                              sep="\t", usecols=[self.column_samples],
                                              dtype=str)[self.column_samples].unique()
                    elif len(pd.read_csv(BytesIO(self.file.value), nrows=5).columns) > 2:
                        samples = pd.read_csv(BytesIO(self.file.value),
                                              usecols=[self.column_samples],
                                              dtype=str)[self.column_samples].unique()
                    elif len(pd.read_excel(BytesIO(self.file.value), nrows=5, sep="\t").columns) > 2:
                        samples = pd.read_excel(BytesIO(self.file.value),
                                                usecols=[self.column_samples],
                                                dtype=str)[self.column_samples].unique()
                    else:
                        samples = []
                    if len(samples) == 0:
                        raise ValidationError("Please specify the column containing samples names correctly")
                else:
                    samples = [i for i in self.columns if re.match(self.column_mainset, i)]
                    if len(samples) == 0:
                        raise ValidationError("Please specify the regular expression for the main data column correctly")
                fractions = [re.match(self.name_pattern, i).group("frac") for i in samples if re.match(self.name_pattern, str(i))]
                fractions = natsort.natsorted(set(fractions))
                if len(fractions) == 0:
                    raise ValidationError("The name pattern didn't detect any fractions.")
                self.fraction_ordering.options = fractions
                self._dependent_layout_long[1][0] = self.fraction_ordering
                self._dependent_layout_pivot[1][0] = self.fraction_ordering
            except ValidationError as ex:
                self._dependent_layout_long[1][0] = pn.pane.Markdown(str(ex), sizing_mode="stretch_width")
                self._dependent_layout_pivot[1][0] = pn.pane.Markdown(str(ex), sizing_mode="stretch_width")
            except:
                self._dependent_layout_long[1][0] = pn.pane.Markdown(traceback.format_exc(), sizing_mode="stretch_width")
                self._dependent_layout_pivot[1][0] = pn.pane.Markdown(traceback.format_exc(), sizing_mode="stretch_width")
    
    @param.depends('disabled', watch=True)
    def _disable(self):
        state = self.disabled
        for el in [self.file, self.experiment_name,
                   self._source, self._acquisition, self._level, self._orientation,
                   self._column_custom_ids, self._column_custom_genes,
                   self._column_custom_samples, self._column_custom_mainset_long,
                   self._column_mainset_pivot, self._name_pattern, self._pattern_presets,
                   self.fraction_ordering, self._columns_annotation]:
            el.disabled = state
    
    def get_settings(self):
        settings = dict()
        if self.experiment_name.value == "":
            raise ValueError("Please give the experiment a name")
        settings["filename"] = self.file.filename
        settings["expname"] = self.experiment_name.value
        settings["source"] = self.source
        settings["acquisition"] = self.acquisition
        settings["level"] = self.level
        settings["orientation"] = self.orientation
        settings["original_protein_ids"] = self.column_ids
        settings["genes"] = self.column_genes
        settings["sets"] = {
            SpatialDataSet.defaultsettings["acquisition_modes"][self.acquisition]["sets"][0]: self.column_mainset
        }
        
        settings["name_pattern"] = self.name_pattern
        settings["fractions"] = self.fraction_ordering.get_ordered_labels()
        settings["fraction_mapping"] = self.fraction_ordering.get_mapping()
        settings["columns_annotation"] = self.columns_annotation
        if self.orientation == "long":
            settings["samples"] = self.column_samples
        return settings


class OrderedMapper(Viewer):
    
    options = param.List(doc="List of input values")
    labels = param.List(default=[])
    order = param.List(default=[])
    name = param.String(doc="Name shown above input list", default="Input list")
    disabled = param.Boolean(default=False)
    
    def __init__(self, **params):
        super().__init__(**params)
        self.layout = pn.Row(None, help_icon("Reorder and rename as required. Every order number can only be present once, emptying the label will lead to deletion."))
        self._sync_options()
    
    def __panel__(self):
        return self.layout
    
    @param.depends('options', 'name', 'labels', 'order', watch=True)
    def _sync_options(self):
        if len(self.labels) != len(self.options):
            self.labels = self.options
        if len(self.order) != len(self.options):
            self.order = [el+1 for el in range(len(self.options))]
        gs = pn.GridSpec(sizing_mode="fixed", height=30*(len(self.options)+1), width=400)
        gs[0,0] = pn.pane.Markdown("Order", margin=0)
        gs[0,1:4] = pn.pane.Markdown(self.name, margin=0)
        gs[0,4:7] = pn.pane.Markdown("Label (empty to delete)", margin=0)
        for i,el in enumerate(self.options):
            gs[i+1,0] = pn.widgets.Select(options=[str(j+1) for j in range(len(self.order))],
                                          value=str(self.order[i]), margin=0)
            gs[i+1,1:4] = pn.pane.Markdown(str(el), margin=0)
            gs[i+1,4:7] = pn.widgets.TextInput(value=self.labels[i], width=100, margin=0)
        self.layout[0] = gs
    
    @param.depends('disabled', watch=True)
    def _disable(self):
        state = self.disabled
        for i in range(len(self.options)):
            self.layout[0][i+1,0].disabled = state
            self.layout[0][i+1,4].disabled = state
    
    def get_ordered_labels(self):
        order = [int(el.value) for el in self.layout[0][1:,0]]
        if len(set(order)) != len(order):
            raise ValidationError("Please select each order number exactly once.")
        labels = [el.value for el in self.layout[0][1:,4]]
        return [labels[order.index(i+1)] for i in range(len(order)) if labels[order.index(i+1)] != ""]
    
    def get_mapping(self):
        return {o: el.value if el.value != "" else None for o,el in zip(self.options, self.layout[0][1:,4])}
    
    def set_ordered_mapping(self, ordered, mapping):
        labels = ["" if k not in mapping.keys() else mapping[k] if mapping[k] != None else "" for k in self.options]
        self.labels = labels
        order = [0 if l == "" else ordered.index(l)+1 for l in labels]
        for i,el in enumerate(order):
            if el == 0:
                order[i] = max(order)+1
        self.order = order


    
class ConfigureAnnotations(Viewer):
    organism = param.String(doc="Organism of uploaded data", default="Homo sapiens")
    organelles = param.String(doc="Filename from package ressources or string input from file",
                              default="")
    complexes = param.String(doc="Filename from package ressources or string input from file",
                           default="")
    gene_mode = param.String(default="don't reannotate")
    gene_source = param.String(doc="Filename from package ressources or string input from file",
                               default="")
    custom_option = param.String(default=SelectOrFile.custom_option)
    title = param.String(default="**Protein Annotations**")
    disabled = param.Boolean(default=False)
    width = param.Integer(default=540)
    
    def __init__(self, **params):
        self._organism = pn.widgets.Select(
            options=["Homo sapiens",
                     "Mus musculus",
                     "Saccharomyces cerevisiae",
                     "Arabidopsis thaliana",
                     "other"
                    ],
            name="Select organism"
        )
        self._organelles_options = [
            el.split(".")[0]
            for el in pkg_resources.resource_listdir("domaps", "annotations/organellemarkers")
            if el.endswith(".csv")]
        self._complexes_options = [
            el.split(".")[0]
            for el in pkg_resources.resource_listdir("domaps", "annotations/complexes")
            if el.endswith(".csv")]
        self._fasta_options = [
            el.split(".")[0]
            for el in pkg_resources.resource_listdir("domaps", "annotations/idmapping")
            if el.endswith(".txt")]
        self._tab_options = [
            el.split(".")[0]
            for el in pkg_resources.resource_listdir("domaps", "annotations/idmapping")
            if el.endswith(".tab")]
        self._organelles = SelectOrFile(accept=".csv", setting="organelles",
                                        title="Organelle markers", options=self._organelles_options)
        self._complexes = SelectOrFile(accept=".csv", setting="complexes",
                                       title="Protein complexes", options=self._complexes_options)
        self._gene_mode = pn.widgets.Select(
            options=["don't reannotate",
                     "from uniprot fasta headers",
                     "from uniprot.org",
                     "from uniprot tab download"],
            name="Reannotation of gene names")
        self._gene_source = SelectOrFile(
            setting="reannotation_source",
            title="Annotation file"
        )
        super().__init__(**params)
        self._layout = pn.Card(
            objects=[
                pn.Row(self._organism, help_icon("The tool uses marker proteins for organelles/subcellular locations and members of stable protein complexes for benchmarking. You can either use predefined lists curated by the authors, or upload your own .csv files (check the abouts page for details). Additionally you can decide to reannotate gene symbols from a selected source. If you want to use the gene symbols from the uniprot fasta headers, you need to reduce the file to headers only.")),
                pn.Row(self._organelles, self._complexes),
                pn.Row(self._gene_mode, self._update_gene_source)
            ],
            collapsed=False
        )
        self._sync_layout()
        self._organism.param.trigger('value')
        self._sync_values()
        
        
        
    def __panel__(self):
        return self._layout
    
    @param.depends('custom_option', 'width', 'title', watch=True)
    def _sync_layout(self):
        self._layout.width = self.width
        self._layout.header = pn.pane.Markdown(self.title, width=self.width)
        self._organism.width = self.width-20-45
        self._organelles.width = self.width//2
        self._complexes.width = self.width//2
        self._gene_mode.width = self.width//2-20
        self._gene_source.width = self.width//2
        self._organelles.custom_option = self.custom_option
        self._complexes.custom_option = self.custom_option
        self._gene_source.custom_option = self.custom_option
        
    
    @param.depends('organism', watch=True)
    def _sync_widgets(self):
        self._organism.value = self.organism
    
    @param.depends('organelles', 'complexes', 'gene_mode', 'gene_source', watch=True)
    def _sync_values(self):
        self._organelles.value = self.organelles
        self._complexes.value = self.complexes
        self._gene_mode.value = self.gene_mode
        self._gene_source.value = self.gene_source
    
    @param.depends('_organism.value', watch=True)
    def _sync_organism(self):
        self.organism = self._organism.value
        ov, cv = self.organelles, self.complexes
        self._organelles.options = [el for el in self._organelles_options if el.startswith(self.organism)]
        self._complexes.options = [el for el in self._complexes_options if el.startswith(self.organism)]
        if ov == "" or (ov not in self._organelles.options and ov in self._organelles_options):
            try:
                self._organelles.value = self._organelles.options[0]
            except:
                pass
        else:
            self._organelles.value = ov
        if cv == "" or (cv not in self._complexes.options and cv in self._complexes_options):
            try:
                self._complexes.value = self._complexes.options[0]
            except:
                pass
        else:
            self._complexes.value = cv
        self._gene_mode.param.trigger("value")
    
    @param.depends('_gene_mode.value')
    def _update_gene_source(self):
        self.gene_mode = self._gene_mode.value
        if self._gene_mode.value == "from uniprot fasta headers": #fasta_headers
            mode = "fasta"
        elif self._gene_mode.value == "from uniprot tab download": #tsv
            mode = "tab"
        else:
            return pn.Column()
        
        gv = self.gene_source
        opts = [el for el in self._fasta_options if el.startswith(self.organism)] if mode == "fasta" else\
               [el for el in self._tab_options if el.startswith(self.organism)]
        self._gene_source.options = opts
        self._gene_source.accept = ".txt" if mode == "fasta" else ".tab"
        
        if gv == "" or (gv not in opts and (gv in self._fasta_options or gv in self._tab_options)):
            try:
                self._gene_source.value = self._gene_source.options[0]
            except:
                pass
        else:
            self._gene_source.value = gv
        return self._gene_source
    
    @param.depends('_organelles.value', '_organelles.custom_file', watch=True)
    def _sync_organelles(self):
        if self._organelles.value in self._organelles.options:
            self.organelles = self._organelles.value
        else:
            if self._organelles.custom_file is None:
                self.organelles = ""
            else:
                self.organelles = self._organelles.custom_file
        self._layout[1][0] = self._organelles
    
    @param.depends('_complexes.value', '_complexes.custom_file', watch=True)
    def _sync_complexes(self):
        if self._complexes.value in self._complexes.options:
            self.complexes = self._complexes.value
        else:
            if self._complexes.custom_file is None:
                self.complexes = ""
            else:
                self.complexes = self._complexes.custom_file
        self._layout[1][1] = self._complexes
    
    @param.depends('_gene_source.value', '_gene_source.custom_file', watch=True)
    def _sync_genes(self):
        if self._gene_source.value in self._fasta_options or self._gene_source.value in self._tab_options:
            self.gene_source = self._gene_source.value
        else:
            if self._gene_source.custom_file is None:
                self.gene_source = ""
            else:
                self.gene_source = self._gene_source.custom_file
        #self._layout[2][1] = self._gene_source
    
    @param.depends('disabled',
                   watch=True)
    def _disable(self):
        status = self.disabled
        self._organism.disabled = status
        self._organelles.disabled = status
        self._complexes.disabled = status
        self._gene_mode.disabled = status
        self._gene_source.disabled = status
    
    def get_settings(self):
        reannotate = {
            "don't reannotate": False,
            "from uniprot fasta headers": "fasta_headers",
            "from uniprot.org": "uniprot",
            "from uniprot tab download": "tsv"
        }[self.gene_mode]
        return dict(
            organism=self.organism,
            reannotate=reannotate,
            reannotation_source=self.gene_source,
            organelles=self.organelles,
            complexes=self.complexes
        )


class ConfigureDataTransformations(Viewer):
    
    unlog = param.Selector(default = False)
    invert = param.Boolean(default = False)
    samplenormalization = param.Selector(default = None)
    yields = param.List(default=[])
    fractions = param.List(default=[])
    title = param.String(default="**Data Transformations**")
    width = param.Integer(default=540)
    disabled = param.Boolean(default=False)
    
    def __init__(self, **params):
        self._unlog_chk = pn.widgets.Checkbox(name="Data was previously log-transformed with base:")
        self._unlog_base = pn.widgets.Select(options=[2,10,"e"], margin=0)
        self._invert_chk = pn.widgets.Checkbox(name="Data needs to be inverted to represent biological profiles.")
        self._weigh_yields = pn.widgets.Checkbox(name="Weigh fractions (will be applied before 0-1 normalization)")
        self._sample_norm = pn.widgets.Select(options=[None, "sum", "median"], margin=0)
        super().__init__(**params)
        self._layout = pn.Card(
            pn.Row(self._invert_chk, help_icon("These transformations will be done after indexing and filtering. Inverting the data is meant e.g. for SILAC ratios. Undoing a log transformation is only required if data has been previously processed with a different tool.")),
            pn.Row(self._unlog_chk, self._unlog_base),
            pn.Row("Sample normalization:", self._sample_norm),
            #self._weigh_yields, TODO: Implement yield weighing and uncomment this line afterwards
            pn.Row(),
            collapsed=True
        )
        if len(self.yields) != 0:
            if len(self.fractions) != len(self.yields):
                raise ValidationError("Fractions and yields must have the same length.")
            self._weigh_yields.value = True
        self._sync_layout()
        self._sync_widgets()
    
    def __panel__(self):
        return self._layout
    
    @param.depends('title', 'width', watch=True)
    def _sync_layout(self):
        self._layout.width = self.width
        self._layout.header = pn.pane.Markdown(self.title, width=self.width)
        self._unlog_chk.width = 300
        self._unlog_base.width = 50
        self._invert_chk.width = self.width-20-45
        self._sample_norm.width = 120
        self._weigh_yields.width = self.width-20
    
    @param.depends('unlog', 'invert', 'samplenormalization', 'yields', 'fractions', watch=True)
    def _sync_widgets(self):
        self._invert_chk.value = self.invert
        self._unlog_chk.value = True if self.unlog != False else False
        self._unlog_base.value = self.unlog if self.unlog != False else None
        self._sample_norm.value = self.samplenormalization
    
    @param.depends('_invert_chk.value', watch=True)
    def _sync_invert(self):
        self.invert = self._invert_chk.value
    
    @param.depends('_unlog_chk.value', '_unlog_base.value', watch=True)
    def _sync_unlog(self):
        self.unlog = self._unlog_base.value if self._unlog_chk.value == True else False
        self._unlog_base.disabled = not self._unlog_chk.value
    
    @param.depends('_sample_norm.value', watch=True)
    def _sync_samplenorm(self):
        self.samplenormalization = self._sample_norm.value
    
    @param.depends('yields', 'fractions', 'width', '_weigh_yields.value', watch=True)
    def _layout_yields(self):
        if self._weigh_yields.value == False:
            try:
                self.get_yields()
            except:
                pass
            self._layout[-1] =  pn.Row()
        elif len(self.fractions) == 0:
            self._layout[-1] = "No fractions set"
        else:
            layout = pn.Column()
            n = self.width//80
            for i,el in enumerate(self.fractions):
                if i%n == 0:
                    layout.objects.append(pn.Row())
                layout[-1].objects.append(pn.widgets.FloatInput(name=el,
                                                                value=1 if len(self.yields) == 0 else self.yields[i],
                                                                width=self.width//n-10, margin=5))
            self._layout[-1] = layout
    
    def get_yields(self):
        yields = []
        for r in self._layout[-1][0:]:
            for f in r:
                yields.append(f.value)
        self.yields = yields
        return yields
    
    @param.depends('disabled', watch=True)
    def _disable(self):
        state = self.disabled
        for r in self._layout[-1][0:]:
            for f in r:
                f.disabled = state
        self._unlog_chk.disabled = state
        self._unlog_base.disabled = any([state, not self.unlog])
        self._invert_chk.disabled = state
        self._weigh_yields.disabled = state
        self._sample_norm.disabled = state
        
    
    def get_settings(self):
        settings = dict()
        settings["input_invert"] = self.invert
        settings["input_logged"] = self.unlog
        settings["input_samplenormalization"] = self.samplenormalization
        if self._weigh_yields.value == True:
            self.get_yields()
            settings["yields"] = self.yields
        return settings


class ConfigureColumnFilter(Viewer):
    
    column = param.String(default="")
    operator = param.String(default="!=")
    value = param.String(default="")
    columns = param.List(default=[])
    width = param.Integer(default=540)
    settings = param.ClassSelector(class_=dict)
    disabled = param.Boolean(default=False)
    
    def __init__(self, **params):
        self._column = pn.widgets.Select()
        self._operator = pn.widgets.Select(options=[">",">=","<","<=", "==","!="], width=60)
        self._value = pn.widgets.TextInput()
        super().__init__(**params)
        self._layout = pn.Row(self._column, self._operator, self._value)
        self._sync_layout()
        
    def __panel__(self):
        return self._layout
    
    @param.depends('column', 'operator', 'value', 'columns', 'width', watch=True)
    def _sync_layout(self):
        self._layout.width = self.width
        self._column.width = (self.width-80)//2-20
        self._value.width = (self.width-80)//2-20
        #cv = self.column
        #self._column.options = self.columns
        #if cv is not "":
        #    self._column.value = cv
        self._value.param.set_param(value=self.value)
        self._operator.value = self.operator
        self._column.param.set_param(options=self.columns, value=self.column)
    
    @param.depends('_column.value', watch=True)
    def _sync_column(self):
        if self.column != self._column.value:
            self.column = self._column.value
            self.get_settings()
    
    @param.depends('_value.value', watch=True)
    def _sync_value(self):
        if self.value != self._value.value:
            self.value = self._value.value
            self.get_settings()
    
    @param.depends('_operator.value', watch=True)
    def _sync_operator(self):
        if self.operator != self._operator.value:
            self.operator = self._operator.value
            self.get_settings()
    
    @param.depends('disabled', watch=True)
    def _disable(self):
        state = self.disabled
        self._column.disabled = state
        self._value.disabled = state
        self._operator.disabled = state
    
    def get_settings(self):
        self.settings = {self.column: [self.operator, self.value]}
        return self.settings


class ConfigureColumnFilters(Viewer):
    
    value = param.ClassSelector(class_=dict, default=dict())
    columns = param.List(default=[])
    width = param.Integer(default=450)
    title = param.String(default="**Simple column filters**")
    disabled = param.Boolean(default=False)
    
    def __init__(self, **params):
        self.btn_add = pn.widgets.Button(name="+", width=32, button_type="success")
        self.btn_add.on_click(self.add_filter)
        super().__init__(**params)
        self._layout = pn.Card(pn.Row(
            "Add filters to only keep entries matching the conditions.\nOnly one filter per column", help_icon("This is intended to e.g. remove contaminants, that are flagged in a separate column, or to filter for scores/aggregated quality measures on a global level. If the column only contains a single level after filtering, it is removed from the output.", width=30), self.btn_add),
            collapsed=True)
        self._sync_layout()
        self._sync_filters()
    
    def __panel__(self):
        return self._layout
    
    @param.depends('width', 'title', watch=True)
    def _sync_layout(self):
        self._layout.width = self.width
        self._layout[0][0].width = self.width-92
        self._layout.header = pn.pane.Markdown(self.title, width=self.width)
    
    @param.depends('value', 'columns', watch=True)
    def _sync_filters(self):
        self._layout.objects = [pn.Row(self._layout[0])]
        for i,k in enumerate(self.value.keys()):
            x = pn.widgets.Button(name="X", width=32, button_type="danger")
            def delete(event, index=i):
                self._layout[index+1] = pn.Row()
            y = self.value[k][1]
            if y.startswith("'") and y.endswith("'"):
                y = y[1:-1]
            self._layout.append(pn.Row(
                ConfigureColumnFilter(columns=self.columns, operator=self.value[k][0],
                                      column=k, value=y, width=self.width-50),
                x))
            x.on_click(delete)
    
    def add_filter(self, event):
        x = pn.widgets.Button(name="X", width=30, button_type="danger")
        n = len(self._layout.objects)
        def delete(event):
            self._layout[n] = pn.Row()
        self._layout.append(pn.Row(
            ConfigureColumnFilter(columns=self.columns, width=self.width-50),
            x))
        x.on_click(delete)
    
    @param.depends('disabled', watch=True)
    def disable(self):
        state = self.disabled
        self.btn_add.disabled = state
        for r in self._layout:
            for el in r:
                if type(el) == pn.layout.base.Row:
                    for widget in el:
                        widget.disabled = state
                else:
                    el.disabled = state
    
    def get_settings(self):
        value = dict()
        for i,r in enumerate(self._layout):
            if len(r) == 0 or i == 0:
                continue
            x = r[0][2].value
            if len(x) == 0:
                raise ValidationError(f"Please specify a value or delete filter for column {r[0][0].value}.")
            try:
                float(x)
            except:
                x = "'"+x+"'"
            if r[0][0].value in value.keys():
                raise ValidationError(f"Two or more column filters for {r[0][0].value} were defined.")
            value[r[0][0].value] = [r[0][1].value, 
                                    x]
        self.value = value
        return dict(column_filters = self.value)


class ConfigureSingleFile(Viewer):
    
    settings = param.ClassSelector(class_=dict, default={})
    width = param.Integer(default=540)
    disabled = param.Boolean(default=False)
    
    def __init__(self, **params):
        self._content = ConfigureFileContent(width=540)
        self._ann = ConfigureAnnotations(organism="Homo sapiens")
        self._transform = ConfigureDataTransformations(width=540)
        self._SILAC = ConfigureSILACFilter(width=540)
        self._MSMS = ConfigureMSMScountFilter(width=265)
        self._cons = ConfigureConsecutiveFilter(width=265)
        self._col = ConfigureColumnFilters(width=540)
        self._comment = pn.widgets.TextInput(name="Comment", value="", width=540, height=70)
        self._error_output = pn.Row(pn.pane.Markdown())
        self._btn_save = pn.widgets.FileDownload(callback=self.save_settings, filename="domqc_settings.json",
                                                 name="Save settings", button_type="primary", width=150, disabled=True, align="end")
        self._btn_load = FileInput(width=110, disabled=True, align="end", height=30)
        self._btn_run = pn.widgets.Button(name="Run processing", button_type="success", width=150, disabled=True, align="end")
        super().__init__(**params)
        self._layout = pn.Column(
            self._content, 
            self._file_loaded,
            pn.Row(self._btn_save, pn.Column(pn.pane.Markdown("Load settings", height=40), self._btn_load), self._btn_run),
            self._error_output
            )
        self._layout_loaded = pn.Column(
            self._transform, 
            self._ann, 
            self._select_acquisition,
            self._col,
            self._comment
        )
        self._sync_layout()
        self._set_defaults()
    
    def __panel__(self):
        return self._layout
    
    @param.depends('width', watch=True)
    def _sync_layout(self):
        self._layout.width = self.width
        self._content.width = self.width
        self._ann.width = self.width
        self._transform.width = self.width
        self._SILAC.width = self.width
        self._MSMS.width = self.width//2-5
        self._cons.width = self.width//2-5
        self._col.width = self.width
        self._comment.width = self.width-10
        self._error_output.width = self.width
        self._btn_save.width = self.width//2-75
        self._btn_load.width = 90
        self._btn_run.width = self.width//2-75
    
    @param.depends('_content.file.value')
    def _file_loaded(self):
        if self._content.file.value != None:
            try:
                self._set_defaults()
            except Exception as ex:
                return str(ex)
            self._btn_load.disabled=False
            self._btn_save.disabled=False
            self._btn_run.disabled=False
            return self._layout_loaded
        else:
            return "Load file for further configuration"
    
    @param.depends('_content.source', '_content.level', '_content.orientation', watch=True)
    def _sync_default_settings(self):
        self._set_defaults()
    
    @param.depends('_content.acquisition')
    def _select_acquisition(self):
        self._set_defaults()
        filters = SpatialDataSet.defaultsettings["acquisition_modes"][self._content.acquisition]["quality_filter"]
        if "SILAC_countvar" in filters:
            layout = self._SILAC
            self._SILAC.toggle = True
        else:
            self._SILAC.toggle = False
            layout = pn.Row()
        if "msms_count" in filters:
            layout.append(self._MSMS)
            self._MSMS.toggle = True
        else:
            self._MSMS.toggle = False
        if "consecutive" in filters:
            layout.append(self._cons)
            self._cons.toggle = True
        else:
            self._cons.toggle = False
        return layout
    
    def _set_defaults(self):
        source = "_".join([self._content.source, self._content.level, self._content.orientation])
        if source in SpatialDataSet.defaultsettings["sources"].keys():
            source = SpatialDataSet.defaultsettings["sources"][source]
            
        # Column filters:
        try:
            self._col.columns = self._content.columns
        except Exception as ex:
            self._error_output[0] = str(ex)
        try:
            self._col.value = source["column_filters"]
        except:
            self._col.value = {}
        
        acquisition = SpatialDataSet.defaultsettings["acquisition_modes"][self._content.acquisition]
        # tranformations
        try:
            self._transform.invert = acquisition["input_invert"]
        except:
            pass
        try:
            self._transform.unlog = acquisition["input_unlog"]
        except:
            pass
        try:
            self._transform.samplenormalization = acquisition["input_samplenormalization"]
        except:
            pass
        
        # SILAC Filter
        if self._content.acquisition == "SILAC":
            try:
                self._SILAC.param.set_param(
                    value_count = acquisition["RatioCount"],
                    value_variability = acquisition["RatioVariability"]
                )
            except:
                pass
            try:
                self._SILAC.param.set_param(
                    regex_count = source["sets"][acquisition["sets"][1]],
                    regex_variability = source["sets"][acquisition["sets"][2]],
                    regex_disabled = True
                )
            except:
                self._SILAC.regex_disabled = False
        
        # MSMS Filter
        if self._content.acquisition != "SILAC":
            try:
                self._MSMS.value_count = acquisition["average_MSMS_counts"]
            except:
                pass
            try:
                self._MSMS.regex_count = source["sets"][acquisition["sets"][1]]
                self._MSMS.regex_disabled = True
            except:
                self._MSMS.regex_disabled = False
        
        # Consecutive Filter
        if self._content.acquisition != "SILAC":
            try:
                self._cons.value_count = acquisition["consecutive"]
            except:
                pass
    
    def get_settings(self):
        settings = get_settings([self._content, self._ann, self._transform,
                                 self._SILAC, self._MSMS, self._cons, self._col])
        settings["comment"] = self._comment.value
        if "quality_filter" not in settings.keys():
            settings["quality_filter"] = []
        # the following line will update the interface once with the settings.
        # It should not influence the widgets besides the order sof the fractions, which will still not change the outcome.
        self.settings = settings
        return settings
    
    @param.depends('settings', watch=True)
    def set_settings(self):
        settings = self.settings.copy()
        if len(settings) == 0:
            return
        try:
            # content
            self._content.experiment_name.value = settings["expname"]
            self._content.param.set_param(source=settings['source'],
                                          acquisition=settings['acquisition'],
                                          level=settings['level'],
                                          orientation=settings['orientation'])
            self._content.param.set_param(column_ids= settings['original_protein_ids'],
                                          column_genes= settings['genes'],
                                          column_mainset= settings['sets'][
                                              SpatialDataSet.defaultsettings["acquisition_modes"][self._content.acquisition]["sets"][0]
                                          ],
                                          columns_annotation = [] if "columns_annotation" not in settings.keys() else settings['columns_annotation']
                                         )
            if settings['orientation'] == "long":
                self._content.param.set_param(column_samples= settings['samples'])
            self._content.param.set_param(name_pattern= settings['name_pattern'])
        except Exception as ex:
            self._error_output[0] = pn.Column(settings, traceback.format_exc())
            return
            
        try:
            # fractions
            self._content.fraction_ordering.set_ordered_mapping(settings["fractions"],
                                                                settings["fraction_mapping"])
        except Exception as ex:
            self._error_output[0] = pn.Column(settings, traceback.format_exc())
            return
        
        try:
            # transformations
            self._transform.param.set_param(invert=settings['input_invert'],
                                            unlog=settings['input_logged'],
                                            samplenormalization=settings['input_samplenormalization'])
        except Exception as ex:
            self._error_output[0] = pn.Column(settings, traceback.format_exc())
            return
        
        try:
            # annotations
            self._ann.param.set_param(organism=settings["organism"],
                                      organelles=settings["organelles"],
                                      complexes=settings["complexes"],
                                      gene_mode={False: "don't reannotate",
                                                 "fasta_headers": "from uniprot fasta headers",
                                                 "uniprot": "from uniprot.org",
                                                 "tsv": "from uniprot tab download"}[settings["reannotate"]],
                                      gene_source=settings["reannotation_source"])
        except Exception as ex:
            self._error_output[0] = pn.Column(settings, traceback.format_exc())
            return
        
        try:
            # filters
            if "quality_filter" in settings.keys() and "SILAC_countvar" in settings["quality_filter"]:
                self._SILAC.param.set_param(toggle=True,
                                            value_count=settings["RatioCount"],
                                            value_variability=settings["RatioVariability"],
                                            regex_count=settings["sets"]['Ratio count'],
                                            regex_variability=settings["sets"]['Ratio variability'])
            else:
                self._SILAC.param.set_param(toggle=False)
            if "quality_filter" in settings.keys() and "msms_count" in settings["quality_filter"]:
                self._MSMS.param.set_param(toggle=True,
                                           value_count=settings["average_MSMS_counts"],
                                           regex_count=settings["sets"]['MS/MS count'])
            else:
                self._MSMS.param.set_param(toggle=False)
            if "quality_filter" in settings.keys() and "consecutive" in settings["quality_filter"]:
                self._cons.param.set_param(toggle=True,
                                           value_count=settings["consecutive"])
            else:
                self._cons.param.set_param(toggle=False)
            self._col.param.set_param(value=settings["column_filters"])
        except Exception as ex:
            self._error_output[0] = pn.Column(settings, traceback.format_exc())
            return
        
        try:
            # comment
            self._comment.value = settings['comment']
        except Exception as ex:
            self._error_output[0] = pn.Column(settings, traceback.format_exc())
            return
        
        self._error_output[0] = None
        self.settings = dict()
    
    #@param.depends('_btn_save._clicks', watch=True)
    #def _set_settingsfilename(self):
    #    self._btn_save.filename = "domqc_settings_"+str(self._btn_save._clicks)+".json"
    
    def save_settings(self):
        try:
            settings = self.get_settings()
            settings_2 = self.get_settings()
        except:
            self._error_output[0] = traceback.format_exc()
            return StringIO()
        if settings != settings_2:
            self._error_output[0] = "The settings you are saving will not be loadable exactly as are. Please check them when you reuse them and ideally send them to us for debugging."
        else:
            self._error_output[0] = "Downloading settings."
        sio = StringIO()
        json.dump(
            settings, 
            sio, 
            indent=4
        )
        sio.seek(0)
        self.settings = {}
        return sio
    
    @param.depends('_btn_load.value', watch=True)
    def _load_settings(self):
        if self._btn_load.value != None:
            self._error_output[0] = "Loading settings."
            self.settings = json.load(BytesIO(self._btn_load.value))
            self._btn_load.value = None
    
    @param.depends('disabled', watch=True)
    def _disable(self):
        state = self.disabled
        for el in [self._content,self._ann,self._transform,
                   self._SILAC,self._MSMS,self._cons,self._col,self._comment,
                   self._btn_save,self._btn_load,self._btn_run]:
            el.disabled=state

