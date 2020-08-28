class SpatialDataSet:

    def __init__(self, **kwargs):
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

        self.filename = "6_deep_maps.txt" if "filename" not in kwargs.keys() else kwargs["filename"]

        self.map_of_interest = "MAP1" if "map_of_interest" not in kwargs.keys() else kwargs["map_of_interest"]
        self.cluster_of_interest = "Proteasome" if "cluster_of_interest" not in kwargs.keys() else kwargs[
            "cluster_of_interest"]

        self.map_index = ".*[ .](.*)_.*K" if "map_index" not in kwargs.keys() else kwargs["map_index"]
        self.fraction_index = ".*_(\d+[Kk])" if "fraction_index" not in kwargs.keys() else kwargs["fraction_index"]

        self.summed_MSMS_counts = 2 if "summed_MSMS_counts" not in kwargs.keys() else kwargs["summed_MSMS_counts"]
        self.consecutive_LFQ_I = 4 if "consecutive_LFQ_I" not in kwargs.keys() else kwargs["consecutive_LFQ_I"]

        # self.lfq_filter_param = {"summed MS/MS counts": 2, "consecutive LFQ_I": 4}

        self.regex = {
            "imported_columns": "^[Rr]atio H/L (?!normalized|type|is.*).+|id$|[Mm][Ss].*[cC]ount.+$|[Ll][Ff][Qq].*|.*[nN]ames.*|.*[Pp][rR].*[Ii][Dd]s.*|[Pp]otential.[cC]ontaminant|[Oo]nly.[iI]dentified.[bB]y.[sS]ite|[Rr]everse|[Ss]core|[Qq]-[Vv]alue",
            "index_col_silac": "[Pp]rotein.[Ii][Dd]s|[Mm]ajority.[Pp]rotein.[Ii][Dd]s|[Pp]rotein.[Nn]ames|[Gg]ene.[Nn]ames|[Ii][Dd]|[Ss]core|[Qq]-[Vv]alue",
            "index_col_lfq": ".*[Pp][rR].*[Ii][Dd]s.*|.*[nN]ames.*|[Ii][Dd]|[Ss]core|[Qq]-[Vv]alue|MS/MS.count$",

            # "map": ".*[ .](.*)_.*K",
            # "fraction_silac": ".*_(\d+[Kk])",
            # "fraction_lfq": ".*_(\d+[Kk])",

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

        self.markerproteins = {
            "Proteasome": ["PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1", "PSMB2", "PSMB3",
                           "PSMB4", "PSMB5", "PSMB6", "PSMB7"],
            "CCT complex": ["CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8","CCT6B", "TCP1"],
            "V-type proton ATPase": ["ATP6AP1", "ATP6V0A1", "ATP6V0A2", "ATP6V0A4", "ATP6V0D1", "ATP6V1A", "ATP6V1B2",
                                     "ATP6V1C1", "ATP6V1E1", "ATP6V1G1", "ATP6V1H"],
            "EMC": ["EMC1", "EMC2", "EMC3", "EMC4", "EMC7", "EMC8", "EMC10","EMC6","EMC9"],
            "Lysosome" : ["LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5", "LAMP1", "LAMP2", "CTSA", "CTSB", "CTSC", "CTSD", "CTSL", "CTSZ"]
            }


        self.acquisition = "SILAC" if "acquisition" not in kwargs.keys() else kwargs["acquisition"]


    def data_reading(self):
        """Data import"""

        self.df_original = pd.read_csv(self.filename, sep="\t", comment="#",
                                       usecols=lambda x: bool(re.match(self.regex["imported_columns"], x)))


    def analyze(self):
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

        def indexingdf_silac(df_filt, regex, map_index, fraction_index):
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

            desired_columns_types = [e1 for e1 in [
                re.findall(regex["type_count_silac"], column) or re.findall(regex["type_var_silac"],
                                                                            column) or re.findall(
                    regex["type_ratio_silac"], column) for column in df_index.columns] if e1]
            col_list_type_silac = [item for sublist in desired_columns_types for item in sublist]

            columns = pd.MultiIndex.from_arrays([
                [re.findall(map_index, column)[0] for column in df_index.columns],
                [re.findall(fraction_index, column)[0] for column in df_index.columns],
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

            len_fractions = len(df_stack.index.get_level_values("Fraction").unique())

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
            df_filteredjoint_stacked = df_filteredjoint_stacked.groupby(["Map", "id"]).filter(
                lambda x: len(x) >= len_fractions)

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

        def indexingdf_lfq(df_filt, regex, map_index, fraction_index):
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
                [re.findall(map_index, column)[0] for column in df_index.columns],
                [re.findall(fraction_index, column)[0] for column in df_index.columns],
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
                    number_fractions * self.summed_MSMS_counts)]

            df_consecutive_mapstacked = df_mscount_mapstacked.copy()

            # series no dataframe is generated; if there are at least i.e. 4 consecutive non-NANs, data will be retained
            df_consecutive_mapstacked = df_consecutive_mapstacked.loc[
                df_consecutive_mapstacked[("LFQ intensity")].apply(lambda x: any(
                    np.invert(np.isnan(x)).rolling(window=self.consecutive_LFQ_I).sum() >=
                    self.consecutive_LFQ_I), axis=1)]

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

        if self.acquisition == "SILAC":
            df_filter = filterdf(self.df_original, self.regex)
            df_index = indexingdf_silac(df_filter, self.regex, self.map_index, self.fraction_index)
            df_filteredjoint_stacked = stringency_silac(df_index)
            df_01_stacked = normalization_01_silac(df_filteredjoint_stacked)
            df_log_stacked = conversion_to_log_silac(df_filteredjoint_stacked)
            self.df_log_stacked = df_log_stacked
            self.df_01_stacked = df_01_stacked
            fractions = df_01_stacked.index.get_level_values("Fraction").unique()
            self.fractions = fractions
            map_names = self.df_01_stacked.index.get_level_values("Map").unique()
            self.map_names = map_names
            return df_01_stacked
            # return df_index


        elif self.acquisition == "LFQ":
            df_filter = filterdf(self.df_original, self.regex)
            df_index = indexingdf_lfq(df_filter, self.regex, self.map_index, self.fraction_index)
            df_consecutive_mapstacked = stringency_lfq(df_index)
            df_01_stacked = normalization_01_lfq(df_consecutive_mapstacked)
            df_log_stacked = conversion_to_log_lfq(df_consecutive_mapstacked)
            self.df_log_stacked = df_log_stacked
            self.df_01_stacked = df_01_stacked
            fractions = df_01_stacked.index.get_level_values("Fraction").unique()
            self.fractions = fractions

            return df_01_stacked


        else:
            return "I don't know this"


    def perform_pca(self):
        """PCA will be performed, using logarithmized data.

        Args:
            self.df_log_stacked: dataframe, in which "MAP" and "Fraction" are stacked; data in the column "log profile"
                originates from logarithmized  "LFQ intensity" and "Ratio H/L", respectively; additionally the columns
                "MS/MS count" and "Ratio H/L count|Ratio H/L variability [%]" are stored as single level indices

        Returns:
                df_pca_all_marker_cluster_maps: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3",
                filtered for marker genes, that are consistent throughout all maps / coverage filtering.
        """

        map_names = self.map_names
        markerproteins = self.markerproteins

        # isolate only logarithmized profile, and unstack "Fraction"
        df_log_stacked = self.df_log_stacked
        df_log_fracunstacked = df_log_stacked["log profile"].unstack("Fraction")

        pca = PCA(n_components=3)

        # df_pca: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3"
        df_pca = pd.DataFrame(pca.fit_transform(df_log_fracunstacked))
        df_pca.columns = ["PC1", "PC2", "PC3"]
        df_pca.index = df_log_fracunstacked.index

        df_pca_all_marker_cluster_maps_unfiltered = pd.DataFrame()

        for maps in map_names:
            for clusters in markerproteins:
                for marker in markerproteins[clusters]:
                    if marker not in df_pca.index.get_level_values("Gene names"):
                        continue
                    plot_try_pca = df_pca.xs((marker, maps), level=["Gene names", "Map"], drop_level=False)
                    df_pca_all_marker_cluster_maps_unfiltered = df_pca_all_marker_cluster_maps_unfiltered.append(
                        plot_try_pca)

        # genes are droped, if they are not present in all maps
        df_pca_all_marker_cluster_maps = df_pca_all_marker_cluster_maps_unfiltered.groupby(["Gene names"]).filter(
            lambda x: len(x) >= len(map_names))
        self.df_pca_all_marker_cluster_maps = df_pca_all_marker_cluster_maps


    def plot_pca(self):
        """
        PCA plot will be generated

        Args:
            df_pca_all_marker_cluster_maps: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3",
                filtered for marker genes, that are consistent throughout all maps / coverage filtering.

        Returns:
            fig_pca: PCA plot, for one protein cluster all maps are plotted
        """

        df_pca_all_marker_cluster_maps = self.df_pca_all_marker_cluster_maps
        map_names = self.map_names
        markerproteins = self.markerproteins

        for maps in map_names:
            df_setofproteins = pd.DataFrame()
            for marker in markerproteins[self.cluster_of_interest]:
                if marker not in df_pca_all_marker_cluster_maps.index.get_level_values("Gene names"):
                    continue
                plot_try_pca = df_pca_all_marker_cluster_maps.xs((marker, maps), level=["Gene names", "Map"],
                                                                 drop_level=False)
                df_setofproteins = df_setofproteins.append(plot_try_pca)

            df_setofproteins.reset_index(inplace=True)
            if maps == map_names[0]:
                fig_pca = go.Figure(
                    data=[go.Scatter3d(x=df_setofproteins.PC1, y=df_setofproteins.PC2, z=df_setofproteins.PC3,
                                       hovertext=df_setofproteins["Gene names"], mode="markers", name=maps
                                       # marker=dict(color=[f'rgb({np.random.randint(0,256)}, {np.random.randint(0,256)}, {np.random.randint(0,256)})' for _   in range(25)])
                                       )])
            else:
                fig_pca.add_trace(go.Scatter3d(x=df_setofproteins.PC1, y=df_setofproteins.PC2, z=df_setofproteins.PC3,
                                               hovertext=df_setofproteins["Gene names"], mode="markers", name=maps
                                               ))

        fig_pca.update_layout(autosize=False, width=500, height=500,
                              title="PCA plot for <br>the protein cluster: {}".format(self.cluster_of_interest))

        return fig_pca


    def plottingdf(self):
        """
        The function allows the plotting of filtered and normalized spatial proteomic data using plotly.express.
        The median profile is also calculated and displayed

        Args:
            df_setofproteins: multiindex dataframe, that contains data about the desired protein cluster,
            stored in the columns "Ratio H/L count", "Ratio H/L variability [%]" and "normalized profile"

        Returns:
            fig_abundance_profiles_and_median: Line plot: relative abundan

        """

        df_setofproteins = self.cluster_isolation_df(self.map_of_interest, self.cluster_of_interest)

        df_setofproteins = df_setofproteins.copy()
        df_proteinset_median = df_setofproteins["normalized profile"].unstack("Fraction").median()

        df_setofproteins = df_setofproteins.reindex(index=natsort.natsorted(df_setofproteins.index))

        df_setofproteins.reset_index(inplace=True)
        fig_abundance_profiles = px.line(df_setofproteins, x="Fraction", y="normalized profile",
                                         color="Gene names",
                                         title="Relative abundance profile for {} of <br>the protein cluster: {}".format(
                                             self.map_of_interest, self.cluster_of_interest))

        # make it available for plotting
        df_proteinset_median.name = "normalized profile"

        df_proteinset_median = df_proteinset_median.reindex(index=natsort.natsorted(df_proteinset_median.index))

        df_proteinset_median = df_proteinset_median.reset_index()
        df_proteinset_median.insert(0, "Gene names", np.repeat("Median profile", len(df_proteinset_median)))
        # return df_proteinset_median


        fig_abundance_profiles_and_median = fig_abundance_profiles.add_scatter(x=df_proteinset_median["Fraction"],
                                                                               y=df_proteinset_median[
                                                                                   "normalized profile"],
                                                                               name="Median profile")

        return fig_abundance_profiles_and_median


    def multiple_iterations(self):
        """
        For each individual protein cluster and for each indiviual fraction the distance to the median will be calculated
        and stored in a multiindex dataframe.

        Args:
            self:
             requires dataframe (df_01_stacked, single level column), stored as attribute (self.df_01_stacked),
             in which "MAP" and "Fraction" are stacked. Additionally the columns "MS/MS count" and
             "Ratio H/L count | Ratio H/L variability [%] | Ratio H/L" are found in LFQ and SILAC data respectively.
             markerproteins, stored as attribute (self.markerproteins), that will be used for plotting

        Returns:
            df_allclusters_onlynorm_fracunstacked: multiindex dataframe, that contains only the normalized data of individual
            protein clusters substracted by the median of the respective protein cluster for each fraction. The fractions are
            unstacked.
        """

        markerproteins = self.markerproteins
        df_01_stacked = self.df_01_stacked

        df_allclusters_onlynorm_fracunstacked_unfiltered = pd.DataFrame()

        map_names = self.map_names

        # for each individual map, and each individual protein cluster, defined in the dictionary markerproteins,
        # the functions "cluster_isolation_df" and "distance_to_median_calculation" will be performed
        for maps in map_names:
            for clusters in markerproteins:
                df_setofproteins = self.cluster_isolation_df(maps, clusters)
                df_distance_to_median_fracunstacked = self.distance_to_median_calculation(df_setofproteins)
                # new column is introduced: Column name = "Cluster"; values: clustername, to wich the individual protein belongs to
                df_distance_to_median_fracunstacked["Cluster"] = clusters
                df_distance_to_median_fracunstacked.set_index("Cluster", inplace=True, append=True)
                # the isolated and processed
                df_allclusters_onlynorm_fracunstacked_unfiltered = df_allclusters_onlynorm_fracunstacked_unfiltered.append(
                    df_distance_to_median_fracunstacked)

        # genes are droped, if they are not present in all maps
        df_allclusters_onlynorm_fracunstacked = df_allclusters_onlynorm_fracunstacked_unfiltered.groupby(["Gene names"]).filter(lambda x: len(x) >= len(map_names))
        self.df_allclusters_onlynorm_fracunstacked = df_allclusters_onlynorm_fracunstacked
        # concatenate both original and filtered dataframe, and drop duplicates
        dfs_dictionary = {'DF1': df_allclusters_onlynorm_fracunstacked_unfiltered,
                          'DF2': df_allclusters_onlynorm_fracunstacked}
        df = pd.concat(dfs_dictionary)
        df_genenames_sortedout = df.drop_duplicates(keep=False)
        # retrieve all gene names, that are sorted out
        genenames_sortedout_index = df_genenames_sortedout.index.get_level_values("Gene names").unique()
        genenames_sortedout_list = genenames_sortedout_index.tolist()

        self.genenames_sortedout_list = genenames_sortedout_list

        return df_allclusters_onlynorm_fracunstacked


    def cluster_isolation_df(self, maps, clusters):
        """
        The function returns a dataframe for a specified proteincluster and map, that will be hand over to the
        definition "distance_to_median_calculation"

        Args:
            maps: specified map (e.g. "MAP1")
            clusters: protein cluster, specified in self.markerproteins (e.g. "Proteasome")
            self:
                requires dataframe (df_01_stacked, single level column), stored as attribute (self.df_01_stacked),
            in which "MAP" and "Fraction" are stacked; the data in the column "normalized profile" is used for plotting.
            Additionally the columns "MS/MS count" and "Ratio H/L count | Ratio H/L variability [%] | Ratio H/L" are
            found in LFQ and SILAC data respectively
                markerproteins, stored as attribute (self.markerproteins), that will be used for plotting

        Returns:
            df_setofproteins: multiindex dataframe, that contains data about the desired protein cluster,
            stored in the columns "Ratio H/L count", "Ratio H/L variability [%]" and "normalized profile"
        """

        df_setofproteins = pd.DataFrame()
        df_01_stacked = self.df_01_stacked
        markerproteins = self.markerproteins

        # datapoints of each individual markerprotein is written into plot_try and appended to df_setofproteins
        for marker in markerproteins[clusters]:
            if marker not in df_01_stacked.index.get_level_values("Gene names"):
                continue
            plot_try = df_01_stacked.xs((marker, maps), level=["Gene names", "Map"], drop_level=False)
            df_setofproteins = df_setofproteins.append(plot_try)

        return df_setofproteins


    def distance_to_median_calculation(self, df_setofproteins):
        """
        A dataframe will be generated, in which the distances to the median profile will be stored.

        Args:
            df_setofproteins: multiindex dataframe, that contains data about the desired protein cluster,
            stored in the columns "Ratio H/L count", "Ratio H/L variability [%]" and "normalized profile"

        Returns:
            df_distance_to_median_fracunstacked: multiindex dataframe, for each individual protein and fraction the distances
            are claculated. The data is stored in single level columns (Column names: e.g. "03K", "06K", "12K", "24K", "80K")
        """
        # calculate the median of the "normalized profile" of individual proteins over the fractions
        df_proteinset_median = df_setofproteins["normalized profile"].unstack("Fraction").median()

        # substract the median individually from each entry (e.g. 03K_protein_xy - 03K_median)
        df_distance_to_median_fracunstacked = df_setofproteins["normalized profile"].unstack("Fraction")
        df_distance_to_median_fracunstacked = df_distance_to_median_fracunstacked.apply(
            lambda x: x - df_proteinset_median, axis=1)

        return df_distance_to_median_fracunstacked


    def boxplot_distance_to_median_manymaps(self):
        """
        A box plot for 1 desired cluster, 1 fraction, and across all maps is generated displaying the
        distribution of the distance to the median. For each fraction, one box plot will be displayed.

        Args:
            self:
            df_allclusters_onlynorm_fracunstacked, dataframe with single level column, stored as attribute
            (self.df_allclusters_onlynorm_fracunstacked), in which "Fraction" is unstacked. It contains only the
            normalized data of individual protein clusters substracted by the median of the respective protein cluster
            for each fraction.
            map_names: individual map names are stored as an index

        Returns:
            bx_median_distance_manymaps: Box plot. Along the x-axis, the maps are shown, along the y-axis
            the distances is plotted
        """

        map_names = self.map_names
        df_allclusters_onlynorm_fracunstacked = self.df_allclusters_onlynorm_fracunstacked

        df_boxplot_manymaps = pd.DataFrame()

        # for each individual map and a defined cluster data will be extracted from the dataframe
        # "df_allclusters_onlynorm_fracunstacked" and appended to the new dataframe df_boxplot_manymaps
        for maps in map_names:
            plot_try = df_allclusters_onlynorm_fracunstacked.xs((self.cluster_of_interest, maps),
                                                                level=["Cluster", "Map"], drop_level=False)
            df_boxplot_manymaps = df_boxplot_manymaps.append(plot_try)

        # index will be reset, required by px.violin
        df_boxplot_manymaps = abs(df_boxplot_manymaps.stack("Fraction"))

        df_boxplot_manymaps.name = "distance"

        df_boxplot_manymaps = df_boxplot_manymaps.reindex(index=natsort.natsorted(df_boxplot_manymaps.index))

        df_boxplot_manymaps = df_boxplot_manymaps.reset_index()

        # box plot will be generated, every fraction will be displayed in a single plot
        bx_median_distance_manymaps = px.box(df_boxplot_manymaps, x="Map", y="distance", facet_col="Cluster",
                                             facet_row="Fraction",
                                             boxmode="overlay", height=300 * 5, width=250 * 4, points="all",
                                             hover_name="Gene names",
                                             title="Distribution of the distance to the median for <br>the protein cluster: {}".format(
                                                 self.cluster_of_interest))

        # determine the layout
    #    bx_median_distance_manymaps.update_layout(
#
 #           # setting a black boy around the graph
  #          xaxis=go.layout.XAxis(linecolor='black',
   #                               linewidth=1,
    #                              title="Map",
     #                             mirror=True),
#
 #           yaxis=go.layout.YAxis(linecolor='black',
  #                                linewidth=1,
    #                              mirror=True),
     #   )

        return bx_median_distance_manymaps


    def distance_calculation(self):
        """
        The distances to the median profile of the individual proteins, belonging to specified clusters, across the
        fractions will be used to calculate the distance profile (e.g. Manhattan distance).

        Args:
            self:
            requires dataframe (df_allclusters_onlynorm_fracunstacked, single level column), stored as attribute
            (self.df_allclusters_onlynorm_fracunstacked), in which "Fraction" is unstacked. It contains only the
            normalized data of individual protein clusters substracted by the median of the respective protein cluster
            for each fraction.
            markerproteins: dictionary, key: protein cluster, value: set of gene names, corresponding to a specific cluster

        Returns:
            df_distance_noindex: dataframe, index is reset. It contains the column name "distance", in which the e.g.
            Manhattan distances for each individual protein of the specified clusters (see self.markerproteins) are stored
        """

        markerproteins = self.markerproteins
        df_allclusters_onlynorm_fracunstacked = self.df_allclusters_onlynorm_fracunstacked.copy()

        # np.linalg.norm requires array; ord=1: Manhattan distance will be calculated over 5 dimensions (Fractions)
        distance_array = np.linalg.norm(df_allclusters_onlynorm_fracunstacked.to_numpy(), axis=1, ord=1)

        # series is created, that allows finally the formation of a dataframe again
        distance_series = pd.Series(distance_array)

        # new "Column" name = "Distance"
        distance_series.name = "distance"
        distance_series.index = df_allclusters_onlynorm_fracunstacked.index

        # df is created
        df_distance_noindex = distance_series.reset_index()

        self.df_distance_noindex = df_distance_noindex

        return df_distance_noindex


    def distance_boxplot(self):
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

        # cluster_of_interest = "Proteasome"
        # self.cluster_of_interest=cluster_of_interest

        df_cluster_xmaps_distance_with_index = pd.DataFrame()

        # for each individual map and a defined cluster data will be extracted from the dataframe
        # "df_distance_map_cluster_gene_in_index" and appended to the new dataframe df_cluster_xmaps_distance_with_index
        for maps in map_names:
            plot_try = df_distance_map_cluster_gene_in_index.xs((self.cluster_of_interest, maps),
                                                                level=["Cluster", "Map"], drop_level=False)
            df_cluster_xmaps_distance_with_index = df_cluster_xmaps_distance_with_index.append(plot_try)

        # index will be reset, required by px.box
        df_cluster_xmaps_distance = df_cluster_xmaps_distance_with_index.reset_index()

        # optinal: points=all (= next to each indiviual boxplot the corresponding datapoints are displayed)
        distance_boxplot_figure = px.box(df_cluster_xmaps_distance, x="Map", y="distance",points="all",
                                         hover_name="Gene names",
                                         title="Manhattan distance distribution for <br>the protein cluster: {}".format(
                                             self.cluster_of_interest)
                                         )

        # determine the layout
        distance_boxplot_figure.update_layout(
            autosize=False,
            width=500,
            height=500,

            # setting a black boy around the graph
            xaxis=go.layout.XAxis(linecolor='black',
                                  linewidth=1,
                                  title="Map",
                                  mirror=True),

            yaxis=go.layout.YAxis(linecolor='black',
                                  linewidth=1,
                                  title="distance",
                                  mirror=True),
        )

        return distance_boxplot_figure


    def results_overview_table(self):
        """
        Dataframe will be created, that provides information about "range", "mean" and "standardeviation",
        given as the column names, based on the data given in df_distance_noindex

        Args:
            df_distance_noindex: stored as attribute (self.df_distance_noindex),index is reset.
            It contains the column name "distance", in which the e.g. Manhattan distances for each individual protein
            of the specified clusters (see self.markerproteins) are stored
            markerproteins
        """
        df_distance_noindex = self.df_distance_noindex
        df_distance_map_cluster_gene_in_index = df_distance_noindex.set_index(["Gene names", "Map", "Cluster"])
        map_names = self.map_names
        markerproteins = self.markerproteins

        df_overview = pd.DataFrame()

        for maps in map_names:
            for clusters in markerproteins:
                plot_try = df_distance_map_cluster_gene_in_index.xs((clusters, maps), level=["Cluster", "Map"],
                                                                    drop_level=False)
                statistic_table = {"range": (plot_try["distance"].max(axis=0)) - (plot_try["distance"].min(axis=0)),
                                   "mean": plot_try["distance"].median(axis=0),
                                   "standardeviation": plot_try["distance"].std(axis=0),
                                   "Cluster": clusters,
                                   "Map": maps}
                statistic_series = pd.Series(data=statistic_table)
                df_statistic_table_individual_cluster = pd.DataFrame(statistic_series).T
                df_overview = df_overview.append(df_statistic_table_individual_cluster)

        df_overview.set_index(["Cluster", "Map"], inplace=True)
        df_overview.sort_index(axis=0, level=0, inplace=True)

        return df_overview


    def __repr__(self):
        return "This is a spatial dataset with {} lines.".format(len(self.df_original))