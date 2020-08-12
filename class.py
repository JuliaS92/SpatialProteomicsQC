class SpatialDataSet:

    def __init__(self, filename):
        """reading the data considering the seperation-mode and comments"""
        self.df = pd.read_csv(filename, sep="\t", comment="#")
        self.regex = {
            "index_col": "[Pp]rotein.[Ii][Dd]s|[Mm]ajority.[Pp]rotein.[Ii][Dd]s|[Pp]rotein.[Nn]ames|[Gg]ene.[Nn]ames|[Ii][Dd]",
            "Map": ".*([Mm][Aa][pP]\d+).*",
            "Fraction": ".*_(\d+[Kk])",
            "Type": "(.*).[Mm][Aa][pP].*",
            "filter1": "[Pp]otential.[cC]ontaminant",
            "filter2": "[Oo]nly.[iI]dentified.[bB]y.[sS]ite",
            "filter3": "[Rr]everse"}
        self.markerproteins = {"Proteasome": ["PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1",
                                              "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7"]}


    def analyze(self, acquisition):
        """
        Args:
            acquisition mode: SILAC, LFQ
        """


        def filterdf_SILAC(df, regex):
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

        def indexingdf_SILAC(df_filt, regex):
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
            col_to_index = [col for col in df_filt.columns if re.match(regex["index_col"], col)]
            # pro column, für jede re: schaut ob spaltenname zu re passt: list out of 1*true and x*false,
            df_index = df_index.set_index(col_to_index)
            # es werden die colum names als index festgelegt
            df_index = df_index.drop([column for column in df_index.columns if not column.startswith("Ratio")], axis=1)

            columns = pd.MultiIndex.from_arrays([
                [re.findall(regex["Map"], colum)[0] for colum in df_index.columns],
                [re.findall(regex["Fraction"], colum)[0] for colum in df_index.columns],
                [re.findall(regex["Type"], colum)[0] for colum in df_index.columns]],
                names=["Map", "Fraction", "Type"])
            df_index.columns = columns

            df_stack = df_index.stack(["Fraction", 'Map'])
            df_filt2 = df_stack.loc[[count >= 3 or (count >= 2 and var < 30) for var, count in
                                     zip(df_stack["Ratio H/L variability [%]"], df_stack['Ratio H/L count'])]]

            df_normSILAC = df_filt2["Ratio H/L"].unstack(["Fraction", "Map"]).apply(
                lambda x: x / np.median([el for el in x if not np.isnan(el)]), axis=0).stack(["Map", "Fraction"])
            df_filt3 = df_filt2[["Ratio H/L count", "Ratio H/L variability [%]"]].join(pd.DataFrame(df_normSILAC,
                                                                                                    columns=[
                                                                                                        "Ratio H/L"]))
            df_filt3 = df_filt3.groupby(["Map", "id"]).filter(lambda x: len(x) >= 5)
            df_unstack = df_filt3["Ratio H/L"].unstack("Fraction").transform(lambda x: 1 / x)

            df_norm2 = df_unstack.div(df_unstack.sum(axis=1), axis=0)

            df_final = df_filt3[["Ratio H/L count", "Ratio H/L variability [%]"]].join(pd.DataFrame
                                                                                       (df_norm2.stack("Fraction"),
                                                                                        columns=["Ratio H/L"]))
            return df_final

        if acquisition == "SILAC":
            df_filt = filterdf_SILAC(self.df, self.regex)
            df_final = indexingdf_SILAC(df_filt, self.regex)
            self.df_final = df_final
            fractions = df_final.index.get_level_values("Fraction").unique()
            self.fractions = fractions
            # fractions als return value

            return df_final.head()

            # return df_processed
        elif acquisition == "LFQ":
            return self.df.tail()
        else:
            return "I don't know this"

    def plottingdf_SILAC(self):

        df_plot = self.df_final
        df_plot = df_plot.reset_index()
        plot_try = df_plot.xs("PSMA1", "MAP1"), level("Gene names", "Map")

        # goi=loc[df_plot["Gene names"] == self.markerproteins["Proteasome"][0]
        # g = sns.catplot(data = df_plot[df_plot["Gene names"] == gene], x = "Fraction", y = "Ratio H/L")
        # fig = px.line(plot_try, x="", y="Ratio H/L", title='Relative abundance profile')

        # fig = px.line(df_plot['Ratio H/L count','MAP1',df_plot["Gene names"] == markerproteins["Proteasome"],
        #                     x="Fraction", y="Ratio H/L", title='Relative abundance profile')

        #  df_long=pd.melt(df_plot, id_vars=['Ratio H/L count'], value_vars=[markerproteins["Proteasome"]])
        return plot_try

    def __repr__(self):
        return "This is a spatial dataset with {} lines.".format(len(self.df))

    # try catch system
    # try: code block wird ausgeführt
    # except: wenn try Fehler erzeugt hat