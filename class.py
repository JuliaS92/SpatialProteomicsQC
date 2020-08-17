class SpatialDataSet:


    def __init__(self, filename):
        """reading the data considering the seperation-mode and comments"""
        self.filename = filename
        self.df_long = pd.read_csv(filename, sep="\t", comment="#")
        self.df = pd.read_csv(filename, sep="\t", comment="#")
        self.regex = {
            "index_col_silac": "[Pp]rotein.[Ii][Dd]s|[Mm]ajority.[Pp]rotein.[Ii][Dd]s|[Pp]rotein.[Nn]ames|[Gg]ene.[Nn]ames|[Ii][Dd]",
            "Map_silac": ".*([Mm][Aa][pP]\d+).*",
            "Fraction_silac": ".*_(\d+[Kk])",
            "Type_silac": "(.*).[Mm][Aa][pP].*",
            "filter1": "[Pp]otential.[cC]ontaminant",
            "filter2": "[Oo]nly.[iI]dentified.[bB]y.[sS]ite",
            "filter3": "[Rr]everse",

            "col_shortened" : "id$|[Mm][Ss].*[cC]ount.+$|[Ll][Ff][Qq].*|.*[nN]ames.*|.*[Pp][rR].*[Ii][Dd]s.*|[Pp]otential.[cC]ontaminant|[Oo]nly.[iI]dentified.[bB]y.[sS]ite|[Rr]everse|[Ss]core|[Qq]-[Vv]alue",
            "index_col_lfq" : ".*[Pp][rR].*[Ii][Dd]s.*|.*[nN]ames.*|[Ii][Dd]|[Ss]core|[Qq]-[Vv]alue|MS/MS.count$",
            "LFQ_nan": "[Ll][Ff][Qq].*",
            "Map_lfq": ".*[nNTt]{1}[yYtT]{1}[ .](.*)_\d+[Kk]$",
            "Fraction_lfq": ".*_(\d+[Kk])",
            "Type_lfq": "(.*[nNTt]{1}[yYtT]{1})[ .].*_\d+[Kk]$"        }
        self.markerproteins = {"Proteasome": ["PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1",
                                              "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7"]}
        self.lfq_filter_param={"summed MS/MS counts" : 2, "consecutive LFQ_I" : 4}

    def analyze(self, acquisition):
        """
        Args:
            acquisition mode: SILAC, LFQ
        """


        def filterdf_silac(df, regex):
            """"filtering the data in terms of removing entries caracterized by potential contaminants, etc.

            Args:
                df: dataframe
                regex: dictionary, in which the values represent the regular expressions
            """
            df_filt = df.loc[df[[col for col in df.columns if bool(re.match(regex["filter1"], col))][0]] != "+"]
            df_filt = df_filt.loc[df_filt[[col for col in df_filt.columns if
                                           bool(re.match(regex["filter2"], col))][0]] != "+"]
            df_filt = df_filt.loc[df_filt[[col for col in df_filt.columns if
                                           bool(re.match(regex["filter3"], col))][0]] != "+"]
            return df_filt

        def indexingdf_silac(df_filt, regex):
            """ a multiindex will be generated, containing Map, Fraction and Type,
            allowing the stacking and unstacking of the dataframe;
            the dataset is then subjected to 0:1-normalization and filtering processes
            considering the Ratio H/L count and Ratio H/L variability;

            Args:
                filtered dataset
                regular expressions can be accessed by:  SpatialDataSet().regex["key"]="value"

            Returns:
                dataframe, which contains 3 columns (Ratio H/L count,Ratio H/L variability [%], Ratio H/L),
                rest of the information is stored in the index
            """

            df_index = df_filt.copy()
            # es wird eine Kopie des Dataframes angelegt
            col_to_index = [col for col in df_filt.columns if re.match(regex["index_col_silac"], col)]
            # pro column, fuer jede re: schaut ob spaltenname zu re passt: list out of 1*true and x*false,
            df_index = df_index.set_index(col_to_index)
            # es werden die colum names als index festgelegt
            df_index = df_index.drop([column for column in df_index.columns if not column.startswith("Ratio")], axis=1)

            columns = pd.MultiIndex.from_arrays([
                [re.findall(regex["Map_silac"], colum)[0] for colum in df_index.columns],
                [re.findall(regex["Fraction_silac"], colum)[0] for colum in df_index.columns],
                [re.findall(regex["Type_silac"], colum)[0] for colum in df_index.columns]],
                names=["Map", "Fraction", "Type"])
            df_index.columns = columns

            df_stack = df_index.stack(["Fraction", 'Map'])
            df_filt2 = df_stack.loc[[count >= 3 or (count >= 2 and var < 30) for var, count in
                                     zip(df_stack["Ratio H/L variability [%]"], df_stack['Ratio H/L count'])]]

            df_normsilac = df_filt2["Ratio H/L"].unstack(["Fraction", "Map"]).apply(
                lambda x: x / np.median([el for el in x if not np.isnan(el)]), axis=0).stack(["Map", "Fraction"])
            df_filt3 = df_filt2[["Ratio H/L count", "Ratio H/L variability [%]"]].join(pd.DataFrame(df_normsilac,
                                                                                                    columns=[
                                                                                                        "Ratio H/L"]))
            df_filt3 = df_filt3.groupby(["Map", "id"]).filter(lambda x: len(x) >= 5)
            df_unstack = df_filt3["Ratio H/L"].unstack("Fraction").transform(lambda x: 1 / x)

            df_norm2 = df_unstack.div(df_unstack.sum(axis=1), axis=0)

            df_final = df_filt3[["Ratio H/L count", "Ratio H/L variability [%]"]].join(pd.DataFrame
                                                                                       (df_norm2.stack("Fraction"),
                                                                                        columns=["Ratio H/L"]))
            df_final.columns = [col if col != "Ratio H/L" else "normalized profile" for col in df_final.columns]
            return df_final


        def filterdf_lfq(df, regex):
            """"filtering the data in terms of removing entries caracterized by potential contaminants, etc.

            Args:
                df: dataframe
                regex: dictionary, in which the values represent the regular expressions
            """
            df_filt = df.loc[df[[col for col in df.columns if bool(re.match(regex["filter1"], col))][0]] != "+"]
            df_filt = df_filt.loc[df_filt[[col for col in df_filt.columns if
                                           bool(re.match(regex["filter2"], col))][0]] != "+"]
            df_filt = df_filt.loc[df_filt[[col for col in df_filt.columns if
                                           bool(re.match(regex["filter3"], col))][0]] != "+"]

            return df_filt

        def readingdf_lfq(filename, df_long, regex):
            """Creation of a data frame, which contains only desired information

            Args:
                filename: Data, which is supposed to be analysed
                df_long: original data frame, containing all information ("Protein groups")
                regex: regular expression as dictionary: Key (col_shortened) refers to desired column names

            Returns:
                df_shortened: shortened data frame, contains only the information, which was determined via the Values of the Key col_shortened)
                """

            col_list_in_list = []
            col_list_in_list = [e1 for e1 in [re.findall(regex["col_shortened"], colum) for colum in df_long.columns] if
                                e1]
            col_list = [item for sublist in col_list_in_list for item in sublist]
            df = pd.read_csv(filename, sep="\t", comment="#", usecols=col_list)
            return df

        def filterdf_lfq(df, regex):
            """"filtering the data in terms of removing entries caracterized by potential contaminants, etc.

            Args:
                df: raw dataframe as protein groups file, produced by MaxQuant
                regex: dictionary, in which the values represent the regular expressions
            """
            df_filt = df.loc[df[[col for col in df.columns if bool(re.match(regex["filter1"], col))][0]] != "+"]
            df_filt = df_filt.loc[df_filt[[col for col in df_filt.columns if
                                           bool(re.match(regex["filter2"], col))][0]] != "+"]
            df_filt = df_filt.loc[df_filt[[col for col in df_filt.columns if
                                           bool(re.match(regex["filter3"], col))][0]] != "+"]

            return df_filt

        def indexingdf_lfq(df_filt, regex):
            """

            """
            df_index = df_filt.copy()
            col_to_index = [col for col in df_filt.columns if re.match(regex["index_col_lfq"], col)]
            df_index = df_index.set_index(col_to_index)
            df_index = df_index.drop([column for column in df_index.columns if not column.startswith(("LFQ", "MS"))],
                                     axis=1)
                                    #df_indexx.columns was there before
            df_index[[col for col in df_index.columns if re.match(regex["LFQ_nan"], col)]] = df_index[[col
                                for col in df_index.columns if re.match(regex["LFQ_nan"], col)]].replace(0, np.nan)

            columns = pd.MultiIndex.from_arrays([
                [re.findall(regex["Map_lfq"], colum)[0] for colum in df_index.columns],
                [re.findall(regex["Fraction_lfq"], colum)[0] for colum in df_index.columns],
                [re.findall(regex["Type_lfq"], colum)[0] for colum in df_index.columns]],
                names=["Map", "Fraction", "Type"])
            df_index.columns = columns

            df_test = df_index.copy().stack("Fraction")
            number_fractions = len(df_test.index.get_level_values("Fraction").unique())
            # self.number_fractions=number_fractions

            df_index = df_index.stack("Map")
            df_index.columns = df_index.columns.swaplevel(0, 1)
            df_index.sort_index(axis=1, level=0, inplace=True)

            df_ms = df_index.loc[df_index[('MS/MS count')].apply(np.sum,axis=1) >= (
                                             number_fractions * self.lfq_filter_param["summed MS/MS counts"])]
            # for the column MS/MS count, you have to generate the sum, and if the sum is larger than n[fraction]*2, you take it into new f
            df_lfq = df_ms.copy()
            df_lfq = df_lfq.loc[df_lfq[("LFQ intensity")].apply(lambda x: any(
                np.invert(np.isnan(x)).rolling(window=self.lfq_filter_param["consecutive LFQ_I"]).sum() >=
                self.lfq_filter_param["consecutive LFQ_I"]), axis=1)]
            # series no DF; if there are at least i.e. 4 consecutive non-NANs
            df_norm = df_lfq["LFQ intensity"].div(df_lfq["LFQ intensity"].sum(axis=1), axis=0)
            # 0-1 normalization but only for LFQ intensity
            df_ms_lfq = df_lfq.copy()
            df_ms_lfq["LFQ intensity"] = df_norm
            # you overwrite the column ["LFQ intensity"] with df_norm, which contains only the data for LFQ intensity
            df_ms_lfq = df_ms_lfq.fillna(0).stack("Fraction")
            # NAN will be replaced with 0
            df_ms_lfq.columns = [col if col != "LFQ intensity" else "normalized profile" for col in df_ms_lfq.columns]
            # rename columns: "LFQ intensity" into "normalized profile"
            return df_ms_lfq


        if acquisition == "SILAC":
            df_filter = filterdf_silac(self.df, self.regex)
            df_final = indexingdf_silac(df_filter, self.regex)
            self.df_final = df_final
            fractions = df_final.index.get_level_values("Fraction").unique()
            self.fractions = fractions
            # fractions als return value
            return df_final


        elif acquisition == "LFQ":
            df = readingdf_lfq(self.filename, self.df_long, self.regex)
            df_filter = filterdf_lfq(df, self.regex)
            df_ms_lfq = indexingdf_lfq(df_filter, self.regex)
            self.df_final = df_ms_lfq
            return df_ms_lfq

        else:
            return "I don't know this"


    def plottingdf(self):
        """
        The function allows the plotting of filtered spatial proteomic data using plotly.express

        Args:
            requires processed data
            f.e. SILAC-data:
                 Dataframe with multiindex and the column names: Ratio H/L count | Ratio H/L variability [%] | Ratio H/L

        Returns:
            line chart of proteasomal genes of MAP1
        """

        df_complex = pd.DataFrame()
        df_plot = self.df_final
        markerproteins = self.markerproteins
        for marker in markerproteins["Proteasome"]:
            plot_try = df_plot.xs((marker, "EGF_rep1"), level=["Gene names", "Map"], drop_level=False)
     #       plot_try = df_plot.xs((marker, "MAP1"), level=["Gene names", "Map"], drop_level=False)
            plot_try = plot_try.reset_index()
            df_complex = df_complex.append(plot_try)
        fig = px.line(df_complex, x="Fraction", y="normalized profile", color="Gene names", title='Relative abundance profile')
        return fig


    def __repr__(self):
        return "This is a spatial dataset with {} lines.".format(len(self.df))

    # try catch system
    # try: code block wird ausgefuehrt
    # except: wenn try Fehler erzeugt hat