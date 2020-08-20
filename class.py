class SpatialDataSet:

    def __init__(self, filename):
        """Import of the raw file of interest. Dataframe will be generated, which contains only the data of desired
        column names, specified by the dictionary entry regex["col_shortened"]

        Dictionaries are created, that ate used for filtering and plotting, respectively.

        Args:
            filename: raw file obtained by the LFQ/SILAC approach (Protein Groups Files), processed by MaxQuant

        Returns:
            df_original: shortened data frame, contains only the information, which was determined due to
            the values of the key "col_shortened"
        """
        # df_original contains all information of the raw file; tab separated file is imported,
        # without considering comments, marked with #
        self.filename = filename
        self.regex = {
            "imported_columns": "^[Rr]atio H/L (?!normalized|type|is.*).+|id$|[Mm][Ss].*[cC]ount.+$|[Ll][Ff][Qq].*|.*[nN]ames.*|.*[Pp][rR].*[Ii][Dd]s.*|[Pp]otential.[cC]ontaminant|[Oo]nly.[iI]dentified.[bB]y.[sS]ite|[Rr]everse|[Ss]core|[Qq]-[Vv]alue",
            "index_col_silac": "[Pp]rotein.[Ii][Dd]s|[Mm]ajority.[Pp]rotein.[Ii][Dd]s|[Pp]rotein.[Nn]ames|[Gg]ene.[Nn]ames|[Ii][Dd]|[Ss]core|[Qq]-[Vv]alue",
            "index_col_lfq": ".*[Pp][rR].*[Ii][Dd]s.*|.*[nN]ames.*|[Ii][Dd]|[Ss]core|[Qq]-[Vv]alue|MS/MS.count$",

            "map": ".*[ .](.*)_.*K",

            "fraction_silac": ".*_(\d+[Kk])",
            "fraction_lfq": ".*_(\d+[Kk])",

            "type_count_silac": "([Rr]atio.[Hh]/[Ll].[cC]ount)[ .].*_.*K",
            "type_var_silac": "([Rr]atio.[Hh]/[Ll].[Vv]ariability....)[ .].*_.*K",
            "type_ratio_silac": "([rR]atio.[Hh]/[Ll])[ .](?![cC]ount|[Vv]ariability).*_.*K",

            "type_count_lfq": "([Rr]atio.[Hh]/[Ll].[cC]ount)[ .].*_.*K",
            "type_var_lfq": "([Rr]atio.[Hh]/[Ll].[Vv]ariability....)[ .].*_.*K",
            "type_ratio_silac": "([rR]atio.[Hh]/[Ll])[ .](?![cC]ount|[Vv]ariability).*_.*K",

            "type_msms_lfq": "([Mm][Ss]/[Mm][Ss].[cC]ount)[ .].*_.*K",
            "type_intensity_lfq": "([Ll][Ff][Qq].[Ii]ntensity)[ .].*_.*K",

            # "type_lfq": "(.*[nNTt]{1}[yYtT]{1})[ .].*_\d+[Kk]$",

            "lfq_nan": "[Ll][Ff][Qq].*",

            "contaminants": "[Pp]otential.[cC]ontaminant",
            "sites": "[Oo]nly.[iI]dentified.[bB]y.[sS]ite",
            "reverse": "[Rr]everse"
        }
        self.lfq_filter_param = {"summed MS/MS counts": 2, "consecutive LFQ_I": 4}
        self.markerproteins = {"Proteasome": ["PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1",
                                              "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7"]}

        self.df_original = pd.read_csv(filename, sep="\t", comment="#",
                                       usecols=lambda x: bool(re.match(self.regex["imported_columns"], x)))

    def analyze(self, acquisition):
        """Analysis of the SILAC/LFQ data will be performed.

        The dataframe will be filtered, normalized and converted into a dataframe, characterized by a flat column index,
        which can be used for plotting

        Args:
            acquisition mode: "SILAC" or "LFQ", which is referring to the acquisition method

        Returns:
            A dataframe, in which "Fraction" and "Map" are stacked, containing "normalized profile" as column,
            additionally "Ratio H/L count", "Ratio H/L variability [%]" is found for SILAC data and "MS/MS count"
            for LFQ data; represented as a flat column index
        """

        def filterdf(df_original, regex):
            """"The dataframe will be filtered by removing matches to the reverse database, matches only identified by
            site, and potential contaminants.

            Args:
                df_original: raw dataframe, MaxQuant returns as the primary output the "protein groups" file,
                used as the basis for the analysis
                regex: dictionary, unique keys, that correspond to regular expressions, that allow the identification
                of the dataframe entries "Only identified by site", "Reverse" and "Potential contaminant", respectively
            """

            # f.e. matches to the reverse database are depicted with "+". Therefore only entries without the "+" (!=+)
            # are taken over into the new dataframe
            df_filt = df_original.loc[df_original[[col for col in df_original.columns if
                                                   bool(re.match(regex["contaminants"], col))][0]] != "+"]
            df_filt = df_filt.loc[df_filt[[col for col in df_filt.columns if
                                           bool(re.match(regex["sites"], col))][0]] != "+"]
            df_filt = df_filt.loc[df_filt[[col for col in df_filt.columns if
                                           bool(re.match(regex["reverse"], col))][0]] != "+"]

            return df_filt

        def indexingdf_silac(df_filt, regex):
            """ A multiindex will be generated, characterized by Map, Fraction and Type as level labels,
            allowing the stacking and unstacking of the dataframe;

            Args:
                df_filt: dataframe, that was filtered for matches to the reverse database, matches only identified by
                site, and potential contaminants.

            Returns:
                df_index: multiindex dataframe, which contains 3 level labels: Map (e.g. MAP1, MAP2, ...)
                Fraction (03K, 06K, 12K, 24K, 80K), Type (Ratio H/L count,Ratio H/L variability [%], Ratio H/L),
                rest of the information is stored in the index (Protein IDs, Majority protein IDs,
                Protein names, Gene names, id)
            """

            # deep copy of the dataframe
            df_index = df_filt.copy()

            # iteration through the column names, and searches for column names that are mentioned in the dictionary
            # entry regex["index_col_silac"]
            col_to_index = [col for col in df_filt.columns if re.match(regex["index_col_silac"], col)]

            # column names are defined as index
            df_index = df_index.set_index(col_to_index)

            # all other columns, that are not needed for data processing (they do not start with "Ratio") are removed

            list_drop_col = [column for column in df_index.columns if not column.startswith("Ratio ")]
            list_endcount = [column for column in df_index.columns if column.endswith("count")]
            list_endpct = [column for column in df_index.columns if column.endswith("[%]")]

            list_drop_col.extend(list_endcount)
            list_drop_col.extend(list_endpct)

            df_index = df_index.drop(list_drop_col, axis=1)
            # df_index = df_index.drop([column for column in df_index.columns if not column.startswith("Ratio")], axis=1)

            # multindex will be generated, by isolating the information about the Map, Fraction and Type from each
            # individual column name
            # names=["Map", "Fraction", "Type"] defines the label of the multiindex

            desired_columns_types = [e1 for e1 in [re.findall(regex["type_count_silac"], column) or re.findall(regex["type_var_silac"],
                                column) or re.findall(regex["type_ratio_silac"], column) for column in df_index.columns] if e1]
            col_list_type_silac = [item for sublist in desired_columns_types for item in sublist]

            columns = pd.MultiIndex.from_arrays([
                [re.findall(regex["map"], column)[0] for column in df_index.columns],
                [re.findall(regex["fraction_silac"], column)[0] for column in df_index.columns],
                col_list_type_silac],
                # [re.findall(regex["type_silac"], column)[0] for column in df_index.columns],
                # [re.findall(regex["type_count_silac"], column)[0] or re.findall(regex["type_var_silac"], column)[0] or re.findall(regex["type_ratio_silac"], column)[0] for column in df_index.columns]],
                names=["Map", "Fraction", "Type"])
            df_index.columns = columns

            return df_index

        def stringency_silac(df_index):
            """The multiindex dataframe is subjected to stringency filtering. Only Proteins with complete profiles are
            considered (a set of f.e. 5 SILAC ratios in case you have 5 fractions / any proteins with missing values
            were rejected). Proteins were retained with 3 or more quantifications in each subfraction (=count). Furthermore,
            proteins with only 2 quantification events in one or more subfraction were retained, if their ratio
            variability for ratios obtained with 2 quantification events was below 30% (=var). Subsequently normalization
            to SILAC loading was performed.

            Args:
                df_index: multiindex dataframe, which contains 3 level labels: MAP, Fraction, Type

            Returns:
                df_filteredjoint_stacked: dataframe, in which "MAP" and "Fraction" are stacked;
                the columns "Ratio H/L count", "Ratio H/L variability [%]", and "Ratio H/L" stored as single level indices
            """

            # Fraction and Map will be stacked
            df_stack = df_index.stack(["Fraction", 'Map'])

            # filtering for sufficient number of quantifications (count in 'Ratio H/L count'), taken
            # variability (var in Ratio H/L variability [%]) into account
            # zip: allows direct comparison of count and var
            # only if the filtering parameters are fulfilled the data will be introduced into df_countvarfiltered_stacked
            df_countvarfiltered_stacked = df_stack.loc[[count >= 3 or (count >= 2 and var < 30) for var, count in
                                                        zip(df_stack["Ratio H/L variability [%]"],
                                                            df_stack['Ratio H/L count'])]]

            # "Ratio H/L":normalization to SILAC loading, each individual experiment (FractionXMap) will be divided by its median
            # np.median([...]): only entries, that are not NANs are considered
            df_normsilac_stacked = df_countvarfiltered_stacked["Ratio H/L"].unstack(["Fraction", "Map"]).apply(
                lambda x: x / np.median([el for el in x if not np.isnan(el)]), axis=0).stack(["Map", "Fraction"])

            df_filteredjoint_stacked = df_countvarfiltered_stacked[["Ratio H/L count",
                                                                    "Ratio H/L variability [%]"]].join(
                pd.DataFrame(df_normsilac_stacked, columns=["Ratio H/L"]))

            # dataframe is grouped (Map, id), that allows the filtering for complete profiles
            df_filteredjoint_stacked = df_filteredjoint_stacked.groupby(["Map", "id"]).filter(lambda x: len(x) >= 5)

            # Ratio H/L is converted into Ratio L/H
            df_filteredjoint_stacked["Ratio H/L"] = df_filteredjoint_stacked["Ratio H/L"].transform(lambda x: 1 / x)

            return df_filteredjoint_stacked

        def normalization_01_silac(df_filteredjoint_stacked):
            """The multiindex dataframe, that was subjected to stringency filtering, is 0-1 normalized ("Ratio H/L").

            Args:
                df_filteredjoint_stacked: dataframe, in which "MAP" and "Fraction" are stacked;
                the columns "Ratio H/L count", "Ratio H/L variability [%]", and "Ratio H/L" stored as single level indices

            Returns:
                df_01_stacked: dataframe, in which "MAP" and "Fraction" are stacked; data in the column
                "Ratio H/L" is 0-1 normalized; the columns "Ratio H/L count", "Ratio H/L variability [%]",
                and "Ratio H/L" stored as single level indices; plotting is possible now
            """

            df_01norm_unstacked = df_filteredjoint_stacked["Ratio H/L"].unstack("Fraction")

            # 0:1 normalization of Ratio L/H
            df_01norm_unstacked = df_01norm_unstacked.div(df_01norm_unstacked.sum(axis=1), axis=0)

            df_01_stacked = df_filteredjoint_stacked[["Ratio H/L count", "Ratio H/L variability [%]"]].join(pd.DataFrame
                                                                                                            (
                                                                                                                df_01norm_unstacked.stack(
                                                                                                                    "Fraction"),
                                                                                                                columns=[
                                                                                                                    "Ratio H/L"]))

            # "Ratio H/L" will be renamed to "normalized profile"
            df_01_stacked.columns = [col if col != "Ratio H/L" else "normalized profile" for col in
                                     df_01_stacked.columns]

            return df_01_stacked

        def conversion_to_log_silac(df_filteredjoint_stacked):
            """The multiindex dataframe, that was subjected to stringency filtering, is logarithmized ("Ratio H/L").

            Args:
                df_filteredjoint_stacked: dataframe, in which "MAP" and "Fraction" are stacked;
                the columns "Ratio H/L count", "Ratio H/L variability [%]", and "Ratio H/L" stored as single level indices

            Returns:
                df_log_stacked: dataframe, in which "MAP" and "Fraction" are stacked; data in the column
                "log profile" originates from logarithmized "Ratio H/L" data; the columns "Ratio H/L count",
                "Ratio H/L variability [%]" and  "log profile" are stored as single level indices; PCA is possible now

            """
            # logarithmizing, basis of 2
            df_lognorm_ratio_stacked = df_filteredjoint_stacked["Ratio H/L"].transform(np.log2)
            df_log_stacked = df_filteredjoint_stacked[["Ratio H/L count", "Ratio H/L variability [%]"]].join(
                pd.DataFrame
                (df_lognorm_ratio_stacked, columns=["Ratio H/L"]))

            # "Ratio H/L" will be renamed to "normalized profile"
            df_log_stacked.columns = [col if col != "Ratio H/L" else "log profile" for col in df_log_stacked.columns]

            return df_log_stacked

        def indexingdf_lfq(df_filt, regex):
            """ A multiindex will be generated, characterized by Map, Fraction and Type as level labels,
            allowing the stacking and unstacking of the dataframe;

            Args:
                df_filt: dataframe, that was filtered for matches to the reverse database, matches only identified by
                site, and potential contaminants.

            Returns:
                df_index: multiindex dataframe, which contains 3 level labels: Map (f.e. EGF_rep1, nt_rep1)
                Fraction (03K, 06K, 12K, 24K, 80K), Type (LFQ intensity, MS/MS count),
                rest of the information is stored in the index (Protein IDs, Majority protein IDs, Protein names,
                Gene names, Q-value, Score, id)
            ##same for SILAC
            """
            # deep copy of the dataframe
            df_index = df_filt.copy()

            # iteration through the column names, and searches for column names that are mentioned in the dictionary
            # entry regex["index_col_lfq"]
            col_to_index = [col for col in df_filt.columns if re.match(regex["index_col_lfq"], col)]

            # column names are defined as index
            df_index = df_index.set_index(col_to_index)

            # all other columns, that do not start with "LFQ" and "MS" (essentially "Only identified by site",
            # "Reverse", "Potential contaminant") are removed
            df_index = df_index.drop([column for column in df_index.columns
                                      if not column.startswith(("LFQ", "MS"))], axis=1)

            # "LFQ intensity"-values: converting 0 into NAN
            # iteration through the column names, and searches for column names (LFQ intensity, deposited as regular
            # expression in the dictionary entry regex["lfq_nan"]
            df_index[[col for col in df_index.columns if re.match(regex["lfq_nan"], col)]] = \
                df_index[[col for col in df_index.columns if re.match(regex["lfq_nan"], col)]].replace(0, np.nan)

            # multindex will be generated, by isolating the information about the Map, Fraction and Type from each
            # individual column name
            # names=["Map", "Fraction", "Type"] defines the label of the multiindex

            desired_columns_types = [e1 for e1 in [
                re.findall(regex["type_msms_lfq"], column) or re.findall(regex["type_intensity_lfq"], column) for column
                in df_index.columns] if e1]
            col_list_type_lfq = [item for sublist in desired_columns_types for item in sublist]

            columns = pd.MultiIndex.from_arrays([
                [re.findall(regex["map"], column)[0] for column in df_index.columns],
                [re.findall(regex["fraction_lfq"], column)[0] for column in df_index.columns],
                col_list_type_lfq],
                names=["Map", "Fraction", "Type"])
            df_index.columns = columns

            return df_index

        def stringency_lfq(df_index):
            """The multiindex dataframe is subjected to stringency filtering. Only Proteins which were identified with
            at least [4] consecutive data points regarding the "LFQ intensity", and if summed MS/MS counts >= n(fractions)*[2]
            (LFQ5: min 10 and LFQ6: min 12, respectively; coverage filtering) were included.

            Args:
                df_index: multiindex dataframe, which contains 3 level labels: MAP, Fraction, Typ

            Returns:
                df_consecutive_mapstacked: dataframe, in which "MAP" and "Fraction" are stacked;
                the columns "normalized profile" and "MS/MS count" stored as single level indices
                """

            # retrieve number of fractions that are present in the dataset
            df_fractionnumber_stacked = df_index.copy().stack("Fraction")
            number_fractions = len(df_fractionnumber_stacked.index.get_level_values("Fraction").unique())

            df_index = df_index.stack("Map")

            # level 0 = Fraction, level 1 = Type is converted into level 0 = Type, level 1 = Fraction
            df_index.columns = df_index.columns.swaplevel(0, 1)

            # sorting the level 0, in order to have LFQ intensity -	MS/MS count instead of continuous alternation
            df_index.sort_index(axis=1, level=0, inplace=True)

            # "MS/MS count"-column: take the sum over the fractions;
            # if the sum is larger than n[fraction]*2, it will be stored in the new dataframe
            df_mscount_mapstacked = df_index.loc[df_index[('MS/MS count')].apply(np.sum, axis=1) >= (
                    number_fractions * self.lfq_filter_param["summed MS/MS counts"])]

            df_consecutive_mapstacked = df_mscount_mapstacked.copy()

            # series no dataframe is generated; if there are at least i.e. 4 consecutive non-NANs, data will be retained
            df_consecutive_mapstacked = df_consecutive_mapstacked.loc[
                df_consecutive_mapstacked[("LFQ intensity")].apply(lambda x: any(
                    np.invert(np.isnan(x)).rolling(window=self.lfq_filter_param["consecutive LFQ_I"]).sum() >=
                    self.lfq_filter_param["consecutive LFQ_I"]), axis=1)]

            return df_consecutive_mapstacked

        def normalization_01_lfq(df_consecutive_mapstacked):
            """The multiindex dataframe, that was subjected to stringency filtering, is 0-1 normalized ("LFQ intensity").

            Args:
                df_consecutive_mapstacked: dataframe, in which "Map" is stacked; the column names are stored as
                2 level labels: Fraction (03K, 06K, 12K, 24K, 80K), Type (LFQ intensity, MS/MS count)


            Returns:
                df_01_stacked: dataframe, in which "MAP" and "Fraction" are stacked; data in the column
                "LFQ intensity" is 0-1 normalized; the columns "LFQ intensity" and "MS/MS count" are stored as
                single level indices; plotting is possible now

            """

            # 0-1 normalization but only for LFQ intensity
            df_01norm_mapstacked = df_consecutive_mapstacked["LFQ intensity"].div(
                df_consecutive_mapstacked["LFQ intensity"].sum(axis=1), axis=0)

            df_01_stacked = df_consecutive_mapstacked.copy()

            # you overwrite the column ["LFQ intensity"] with df_01norm_mapstacked, which contains only the data for LFQ intensity
            df_01_stacked["LFQ intensity"] = df_01norm_mapstacked

            # Replace missing values in the remaining maps with 0
            df_01_stacked = df_01_stacked.fillna(0).stack("Fraction")

            # rename columns: "LFQ intensity" into "normalized profile"
            df_01_stacked.columns = [col if col != "LFQ intensity" else "normalized profile" for col in
                                     df_01_stacked.columns]

            return df_01_stacked

        def conversion_to_log_lfq(df_consecutive_mapstacked):
            """The multiindex dataframe, that was subjected to stringency filtering, is logarithmized ("LFQ intensity").

            Args:
                df_01_stacked: dataframe, in which "MAP" and "Fraction" are stacked;
                the columns "Ratio H/L count", "Ratio H/L variability [%]", and "Ratio H/L" stored as single level indices

            Returns:
                df_log_stacked: dataframe, in which "MAP" and "Fraction" are stacked; data in the column "log profile"
                originates from logarithmized  "LFQ intensity"; the columns "log profile" and "MS/MS count" are
                stored as single level indices; PCA is possible now
            """

            df_lognorm_mapstacked = df_consecutive_mapstacked["LFQ intensity"].transform(np.log2)
            df_log_stacked = df_consecutive_mapstacked.copy()
            df_log_stacked["LFQ intensity"] = df_lognorm_mapstacked
            df_log_stacked = df_log_stacked.fillna(0).stack("Fraction")
            df_log_stacked.columns = [col if col != "LFQ intensity" else "log profile" for col in
                                      df_log_stacked.columns]
            return df_log_stacked

        if acquisition == "SILAC":
            df_filter = filterdf(self.df_original, self.regex)
            df_index = indexingdf_silac(df_filter, self.regex)
            df_filteredjoint_stacked = stringency_silac(df_index)
            df_01_stacked = normalization_01_silac(df_filteredjoint_stacked)
            df_log_stacked = conversion_to_log_silac(df_filteredjoint_stacked)
            self.df_log_stacked = df_log_stacked
            self.df_01_stacked = df_01_stacked
            fractions = df_01_stacked.index.get_level_values("Fraction").unique()
            self.fractions = fractions

            return df_log_stacked
            # return df_index

        elif acquisition == "LFQ":
            df_filter = filterdf(self.df_original, self.regex)
            df_index = indexingdf_lfq(df_filter, self.regex)
            df_consecutive_mapstacked = stringency_lfq(df_index)
            df_01_stacked = normalization_01_lfq(df_consecutive_mapstacked)
            df_log_stacked = conversion_to_log_lfq(df_consecutive_mapstacked)
            self.df_log_stacked = df_log_stacked
            self.df_01_stacked = df_01_stacked
            fractions = df_01_stacked.index.get_level_values("Fraction").unique()
            self.fractions = fractions

            return df_log_stacked

        else:
            return "I don't know this"


    def analysis_01_df(self, map):
        """
        The function will generate a new dataframe called "df_setofproteins", which contains exclusively these data
        of proteins, specified in  the dictionary markerproteins

        Args:
            map: data from the specified map (e.g. "MAP1" or ""EGF_rep1") is used for plotting
            requires dataframe (df_01_stacked, single level indices), stored as attribute (self.df_01_stacked),
            in which "MAP" and "Fraction" are stacked; the data in the column "normalized profile" is used for plotting.
            Additionally the columns "MS/MS count" and "Ratio H/L count | Ratio H/L variability [%] | Ratio H/L" are
            found in LFQ and SILAC data respectively
            markerproteins, stored as attribute (self.markerproteins), that will be used for plotting

        Returns:
            self.df_setofproteins: multiindex dataframe, containing the "normalized" of desired markerproteins (e.g. proteasomal genes)
        """

        df_setofproteins = pd.DataFrame()
        df_01_stacked = self.df_01_stacked
        markerproteins = self.markerproteins

        # datapoints of each individual markerprotein is written into plot_try and appended to df_setofproteins
        for marker in markerproteins["Proteasome"]:
            plot_try = df_01_stacked.xs((marker, map), level=["Gene names", "Map"], drop_level=False)
            df_setofproteins = df_setofproteins.append(plot_try)

        self.df_setofproteins = df_setofproteins


    def plottingdf(self):
        """
        The function allows the plotting of filtered spatial proteomic data using plotly.express

        """
        df_setofproteins = self.df_setofproteins.copy()
        df_proteinset_median = df_setofproteins["normalized profile"].unstack("Fraction").median()
        df_setofproteins.reset_index(inplace=True)
        fig_abundance_profiles = px.line(df_setofproteins, x="Fraction", y="normalized profile",
                                         color="Gene names", title='Relative abundance profile')

        # make it available for plotting
        df_proteinset_median.name = "normalized profile"
        df_proteinset_median = df_proteinset_median.reset_index()
        df_proteinset_median.insert(0, "Gene names", np.repeat("Median profile", len(df_proteinset_median)))
        # return df_proteinset_median
        fig_abundance_profiles_and_median = fig_abundance_profiles.add_scatter(x=df_proteinset_median["Fraction"],
                                                                               y=df_proteinset_median[
                                                                                   "normalized profile"],
                                                                               name="Median profile")

        return fig_abundance_profiles_and_median


    def cluster_analysis(self):
        """


        """
        df_setofproteins = self.df_setofproteins.copy()
        # calculate the median of the "normalized profile" of individual proteins over the fractions
        df_proteinset_median = df_setofproteins["normalized profile"].unstack("Fraction").median()

        # substract the median individually from each entry (e.g. 03K_rotein_xy - 03K_median)
        df_onlynorm_fracunstacked = df_setofproteins["normalized profile"].unstack("Fraction")
        df_onlynorm_fracunstacked = df_onlynorm_fracunstacked.apply(lambda x: x - df_proteinset_median, axis=1)

        self.df_onlynorm_fracunstacked = df_onlynorm_fracunstacked

    def distance_boxplot(self):
        """

        """
        markerproteins = self.markerproteins
        df_onlynorm_fracunstacked = self.df_onlynorm_fracunstacked.copy()
        # np.linalg.norm requires array; ord=1: Manhattan distance will be calculated over 5 dimensions (Fractions)
        distance_array = np.linalg.norm(df_onlynorm_fracunstacked.to_numpy(), axis=1, ord=1)

        df_distance_boxplot = pd.DataFrame(distance_array, columns=["Manhattan distance"],
                                           index=markerproteins["Proteasome"]).sort_index()

        distance_boxplot_figure = px.box(df_distance_boxplot, points="all")

        # calculating range/mean/stdv
        statistic_table = {"max-min": df_distance_boxplot.max(axis=0) - df_distance_boxplot.min(axis=0),
                           "mean": df_distance_boxplot.mean(axis=0),
                           "standardeviation": df_distance_boxplot.std(axis=0)}
        df_statistic_distance_boxplot = pd.DataFrame(data=statistic_table)
        self.df_statistic_distance_boxplot = df_statistic_distance_boxplot
        return distance_boxplot_figure


    def violinplot_distance_to_median(self):
        """

        """
        df_onlynorm_fracunstacked = self.df_onlynorm_fracunstacked.copy()
        # convert dataframe, such that the "Fractions" become to column name
        df_onlynorm_fracstacked = df_onlynorm_fracunstacked.sort_index(level="Gene names").stack("Fraction")
        # df_onlynorm_fracstackedseries  contains only one "column", which is renamed to "distance to median"
        df_onlynorm_fracstacked.name = "distance to median"
        df_distance_to_median_singleindex = df_onlynorm_fracstacked.reset_index()

        violinplot_distance_to_median = px.violin(df_distance_to_median_singleindex, x="Fraction",
                                                  y="distance to median")

        return violinplot_distance_to_median


    def perform_pca(self, map):
        """PCA will be performed, using logarithmized data.

        Args:
            self.df_log_stacked: dataframe, in which "MAP" and "Fraction" are stacked; data in the column "log profile"
                originates from logarithmized  "LFQ intensity" and "Ratio H/L", respectively; additionally the columns
                "MS/MS count" and "Ratio H/L count|Ratio H/L variability [%]" are stored as single level indices

        Returns:
                df_pca: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3"
        """

        markerproteins = self.markerproteins

        # isolate only logarithmized profile, and unstack "Fraction"
        df_log_stacked = self.df_log_stacked
        df_log_fracunstacked = df_log_stacked["log profile"].unstack("Fraction")

        pca = PCA(n_components=3)
        df_pca = pd.DataFrame(pca.fit_transform(df_log_fracunstacked))

        df_pca.columns = ["PC1", "PC2", "PC3"]
        df_pca.index = df_log_fracunstacked.index

        df_setofproteins = pd.DataFrame()
        for marker in markerproteins["Proteasome"]:
            plot_try_pca = df_pca.xs((marker, map), level=["Gene names", "Map"], drop_level=False)
            # plot_try_pca = plot_try_pca.reset_index()
            df_setofproteins = df_setofproteins.append(plot_try_pca)

        df_setofproteins.reset_index(inplace=True)
        fig = go.Figure(data=[go.Scatter3d(x=df_setofproteins.PC1, y=df_setofproteins.PC2, z=df_setofproteins.PC3,
                                           hovertext=df_setofproteins["Gene names"])])

        ##, mode="markers", marker=dict(color=[f'rgb({np.random.randint(0,256)}, {np.random.randint(0,256)}, {np.random.randint(0,256)})' for _   in range(25)])
        return fig


    def __repr__(self):
        return "This is a spatial dataset with {} lines.".format(len(self.df))