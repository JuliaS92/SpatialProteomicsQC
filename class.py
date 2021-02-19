class SpatialDataSet:
    
    regex = {
        "imported_columns": "^[Rr]atio H/L (?!normalized|type|is.*).+|id$|[Mm][Ss].*[cC]ount.+$|[Ll][Ff][Qq].*|.*[nN]ames.*|.*[Pp][rR].*[Ii][Dd]s.*|[Pp]otential.[cC]ontaminant|[Oo]nly.[iI]dentified.[bB]y.[sS]ite|[Rr]everse|[Ss]core|[Qq]-[Vv]alue|R.Condition|PG.Genes|PG.ProteinGroups|PG.Cscore|PG.Qvalue|PG.RunEvidenceCount|PG.Quantity"
    }
    
    acquisition_set_dict = {
        "LFQ" : ["[Ll][Ff][Qq].[Ii]ntensity", "[Mm][Ss]/[Mm][Ss].[cC]ount", "[Ii]ntensity"],
        "LFQ Spectronaut" : ["LFQ intensity", "MS/MS count"],
        "SILAC" : [ "[Rr]atio.[Hh]/[Ll](?!.[Vv]aria|.[Cc]ount)","[Rr]atio.[Hh]/[Ll].[Vv]ariability.\[%\]", "[Rr]atio.[Hh]/[Ll].[cC]ount"]
    }
    
    fraction_dict = {
        "1K": "01K","3K": "03K", "6K": "06K", "12K": "12K", "24K": "24K", "80K": "80K", "01K": "01K","03K": "03K", "06K": "06K", "012K": "12K", "024K": "24K",
        "080K": "80K", "Cyt": "Cyt", "Mem": "Mem", "Nuc": "Nuc", "Prot": "Prot", "cyt": "Cyt", "mem": "Mem", "nuc": "Nuc", "Prot": "Prot", "prot": "Prot"
    }
    
    Spectronaut_columnRenaming = {
        "R.Condition": "Map", "PG.Genes" : "Gene names", "PG.Qvalue": "Q-value", "PG.Cscore":"C-Score", 
        "PG.ProteinGroups" : "Protein IDs", "PG.RunEvidenceCount" : "MS/MS count", "PG.Quantity" : "LFQ intensity"
    }
    
    css_color = ["#b2df8a", "#6a3d9a", "#e31a1c", "#b15928", "#fdbf6f", "#ff7f00", "#cab2d6", "#fb9a99", "#1f78b4", "#ffff99", "#a6cee3", 
                      "#33a02c", "blue", "orange", "goldenrod", "lightcoral", "magenta", "brown", "lightpink", "red", "turquoise",
                      "khaki", "darkgoldenrod","darkturquoise", "darkviolet", "greenyellow", "darksalmon", "hotpink", "indianred", "indigo","darkolivegreen", 
                      "coral", "aqua", "beige", "bisque", "black", "blanchedalmond", "blueviolet", "burlywood", "cadetblue", "yellowgreen", "chartreuse",
                      "chocolate", "cornflowerblue", "cornsilk", "darkblue", "darkcyan", "darkgray", "darkgrey", "darkgreen", "darkkhaki", "darkmagenta", 
                      "darkorange", "darkorchid", "darkred", "darkseagreen", "darkslateblue", "snow", "springgreen", "darkslategrey", "mediumpurple", "oldlace", 
                      "olive", "lightseagreen", "deeppink", "deepskyblue", "dimgray", "dimgrey", "dodgerblue", "firebrick", "floralwhite", "forestgreen", 
                      "fuchsia", "gainsboro", "ghostwhite", "gold", "gray", "ivory", "lavenderblush", "lawngreen", "lemonchiffon", "lightblue", "lightcyan",
                      "lightgoldenrodyellow", "lightgray", "lightgrey", "lightgreen", "lightsalmon", "lightskyblue", "lightslategray", "lightslategrey",
                      "lightsteelblue", "lightyellow", "lime", "limegreen", "linen", "maroon", "mediumaquamarine", "mediumblue", "mediumseagreen",
                      "mediumslateblue", "mediumspringgreen", "mediumturquoise", "mediumvioletred", "midnightblue", "mintcream", "mistyrose", "moccasin",
                      "olivedrab", "orangered", "orchid", "palegoldenrod", "palegreen", "paleturquoise", "palevioletred", "papayawhip", "peachpuff", "peru",
                      "pink", "plum", "powderblue", "rosybrown", "royalblue", "saddlebrown", "salmon", "sandybrown", "seagreen", "seashell", "sienna", "silver",
                      "skyblue", "slateblue", "steelblue", "teal", "thistle", "tomato", "violet", "wheat", "white", "whitesmoke", "slategray", "slategrey",
                      "aquamarine", "azure","crimson", "cyan", "darkslategray", "grey","mediumorchid","navajowhite", "navy"]
    
    markerproteins = {
            "Human - Swissprot" :
            {
                "Proteasome" : ["PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7"], 
                        #       "PSMC1", "PSMC2", "PSMC3"],
                #"CCT complex" : ["CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8","CCT6B", "TCP1"],
                #"V-type proton ATPase": ["ATP6AP1", "ATP6V0A1", "ATP6V0A2", "ATP6V0A4", "ATP6V0D1", "ATP6V1A", "ATP6V1B2", "ATP6V1E1", "ATP6V1G1", "ATP6V1H"],
                "EMC" : ["EMC1", "EMC2", "EMC3", "EMC4", "EMC7", "EMC8", "EMC10","EMC6","EMC9"],
                "Lysosome" : ["LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5", "LAMP1", "LAMP2", "CTSA", "CTSB", "CTSC", "CTSD", "CTSL", "CTSZ"],
                #"MCM complex" : ["MCM2", "MCM3", "MCM4", "MCM5", "MCM7"],
                
                "Arp2/3 protein complex" : ["ACTR2", "ACTR3", "ARPC1B", "ARPC2", "ARPC3", "ARPC4", "ARPC5"], 
                #"Prefoldin complex" : [ "PFDN1", "PFDN2", "PFDN4", "PFDN5", "PFDN6", "VBP1"],
                #"AP1 adaptor complex" : ["AP1B1", "AP1G1", "AP1M1", "AP1S1", "AP1S2", "AP1S3"],
                "AP2 adaptor complex" : ["AP2A1", "AP2A2", "AP2B1", "AP2M1",  "AP2S1", ],
                #"AP3 adaptor / AP3-BLOC1 complex" : ["AP3B1", "AP3D1", "AP3M1", "AP3M2", "AP3S1", "AP3S2"],
                #"AP4 adaptor complex" : ["AP4B1", "AP4E1","AP4M1",  "AP4S1"],
                #"Anaphas,e-promoting complex" : ["ANAPC1", "ANAPC10", "ANAPC16", "ANAPC2", "ANAPC4","ANAPC5", "ANAPC7", "CDC16", "CDC23","CDC27"] ,
                #"Rnase/Mrp complex" : ["POP1", "POP4", "POP5", "RPP14","RPP25", "RPP30", "RPP38", "RPP40"],
                "Class C, Vps complex" : ["VPS11","VPS16", "VPS18", "VPS33A"],
                #"Dynactin complex" : ["DCTN1", "DCTN2", "DCTN3", "DCTN4", "DCTN6", "ACTR1A", "CAPZA1"],
                #"CTLH complex" : ["ARMC8", "MAEA", "MKLN1", "RANBP9", "RMND5A"],
                #"Coatomer complex" : ["ARCN1", "COPA", "COPB1", "COPB2", "COPE", "COPG1", "COPZ1"],
                "Wave-2 complex": ["NCKAP1", "CYFIP1", "ABI1", "BRK1", "WASF2"],
                
                "TREX complex / THO complex": ["ALYREF", "THOC5", "THOC2", "THOC1", "THOC3", "DDX39B"], #["THOC5", "THOC2", "THOC1"]
                "Exon junction complex #TREX: ALYREF,DDX39B": ["SRRM1", "RBM8A", "RNPS1", "EIF4A3", "UPF3B", "UPF2"],
                #"TNF-alpha/NF-kappa B signaling complex 5": ["POLR2H", "POLR1A", "POLR1B", "CUL1", "KPNA2"],
                #"Septin complex": ["SEPT7", "SEPT9", "SEPT11", "SEPT8", "SEPT2"],
                #"Sec6/8 exocyst complex": ["EXOC4", "EXOC2", "EXOC1", "EXOC7", "EXOC5", "EXOC3", "EXOC8", "EXOC6"],
                #"SNW1 complex": ["EFTUD2", "SNRNP200", "PRPF8", "MSH2", "DDX23", "SNW1", "PFKL"],
                #"SF3b complex #SF3B4": ["SF3B1", "SF3B3", "SF3B5", "SF3B6", "PHF5A"],
                "CDC5L complex #SF3b: SF3B4": ["SNRPD1", "SNRPD3", "SNRPD2", "PRPF19", "SRSF1", "SF3B2", "SNRPA1", "SF3B4"],
                #"Retromer complex": ["VPS29", "VPS35", "VPS26A"],
                "Respiratory chain complex I": ["NDUFB6", "NDUFB10", "NDUFA10", "NDUFA8", "NDUFA6", "NDUFB11", "NDUFB3", "NDUFB5", "NDUFAB1", "NDUFA4", "NDUFB9",
                                                "NDUFB7", "NDUFA9", "NDUFA5", "NDUFV3", "NDUFA11", "NDUFV1", "NDUFA12", "NDUFV2", "NDUFA7", "NDUFS6", "NDUFS2",
                                                "NDUFA2", "NDUFS8", "NDUFS1"],
                #"RFC complex": ["RFC4", "RFC2", "RFC5", "RFC3", "RFC1"],
                "Nup 107-160 subcomplex": ["NUP85", "NUP37", "NUP160", "NUP98", "NUP107", "NUP133"],
                "Multisynthetase complex": ["EEF1E1", "IARS", "DARS", "EPRS", "AIMP1", "KARS", "LARS", "RARS", "AIMP2", "MARS"],
                #"MCM complex": ["MCM4", "MCM6", "MCM7", "MCM3", "MCM2", "MCM5"],
                "GAA1-GPI8-PIGT-PIG-PIGS complex": ["PIGT", "PIGS", "PIGU", "PIGK", "GPAA1"],
                "Frataxin complex /f1f0: ATP5L": ["SDHA", "HSPD1", "HSPA9", "AFG3L2"],
                "F1F0-ATP synthase": ["ATP5O", "ATP5I", "ATPIF1", "ATP5A1", "ATP5F1", "ATP5B", "ATP5H", "ATP5L", "ATP5J"],
                "Exosome": ["EXOSC1", "EXOSC3", "EXOSC8", "EXOSC4", "EXOSC2", "EXOSC10"],
                #"Large Drosha complex, DGCR8: FUS,HNRNPH1, DDX17, DDX5": ["HNRNPDL", "RALY", "TARDBP", "HNRNPM", "DDX3X", "EWSR1"],
                #"DGCR8 multiprotein complex": ["HNRNPR", "HNRNPH1", "DDX17", "DDX5", "DHX9", "FUS", "NCL"],
                #"COP9 signalosome complex / CNS-P53 complex": ["COPS2", "COPS3", "COPS4", "COPS5", "COPS6", "COPS8", "GPS1"],
                "Arp2/3 protein complex": ["ACTR2", "ARPC1B", "ARPC2", "ARPC3", "ARPC5", "ARPC4", "ACTR3"],
                "60S ribosomal subunit, cytoplasmic": ["RPL10", "RPL10A", "RPL27", "RPL37A", "RPL7A", "RPL23A", "RPL23", "RPL31", "RPL15", "RPL26", "RPL18A", "RPL11",
                                                    "RPL38", "RPL24","RPL36A", "RPL36", "RPL19", "RPL18", "RPL32", "RPL14", "RPL35A", "RPL29", "RPL34", "RPLP0",
                                                    "RPL7", "RPL17", "RPL13", "RPL12", "RPL9", "RPL22", "RPLP1", "RPLP2", "RPL3", "RPL13A", "RPL35", "RPL27A",
                                                    "RPL5", "RPL21", "RPL28", "RPL30", "RPL8", "RPL6", "RPL4"],
                "40S ribosomal subunit, cytoplasmic": ["RPS9", "RPS18", "RPS29", "RPS4X", "RPS6", "RPS15", "FAU", "RPS26", "RPS28", "RPS21", "RPS23", "RPS25", "RPS14",
                                                    "RPS16", "RPS3","RPSA", "RPS2", "RPS12", "RPS19", "RPS27", "RPS17", "RPS5", "RPS20", "RPS3A", "RPS7", "RPS8",
                                                    "RPS10", "RPS15A", "RPS11", "RPS13"],#RPS24 ###RPS17;RPS17L; RPS26;RPS26P11
                "39S ribosomal subunit, mitochondrial": ["MRPL37", "MRPL20", "MRPL9", "MRPL46", "MRPL4", "MRPL44", "MRPL17", "MRPL22", "MRPL39", "MRPL11", "MRPL47",
                                                        "MRPL32", "MRPL48", "MRPL45", "MRPL3", "MRPL19", "MRPL28", "MRPL49", "MRPL2", "MRPL14", "MRPL12", "MRPL55",
                                                        "MRPL10", "MRPL50", "MRPL43", "MRPL24", "MRPL53"],
                "28S ribosomal subunit, mitochondrial": ["MRPS31", "MRPS24", "MRPS30", "MRPS28", "MRPS17", "MRPS7", "MRPS16", "MRPS23", "MRPS18B", "MRPS27", "MRPS9",
                                                        "MRPS34", "DAP3", "MRPS15", "MRPS11", "MRPS36", "MRPS5", "MRPS35", "MRPS10", "MRPS22", "MRPS12", "MRPS6"],
 ###"OLDmitochondrial ribosomal subunits" : ["MRPL1", "MRPL13", "MRPL15", "MRPL16", "MRPL18", "MRPL23", "MRPL30", "MRPL38", "MRPL40",#39S proteins;MRPL12;SLC25A10
            }, 
            "Arabidopsis - Araport" :
            {
                "CCT complex": ["AT3G11830", "AT5G16070", "AT5G20890", "AT5G26360", "AT1G24510", "AT3G18190", "AT3G20050", "AT3G03960"],
                "Coatomer": ["AT1G62020", "AT4G34450", "AT4G31490", "AT1G79990", "AT1G30630", "AT4G08520"],
                "SAGA complex": ["AT5G25150", "AT3G54610", "AT1G54360", "AT1G54140", "AT4G38130", "AT4G31720"],
                "AP1/2": ["AT1G60070","AT2G17380", "AT4G23460", "AT5G22780", "AT1G47830", "AT5G46630", "AT1G10730"],
                "20S proteasome": ["AT1G13060",  "AT1G21720", "AT3G22110", "AT3G22630", "AT5G40580", "AT1G16470", "AT1G47250", "AT1G53850",
                                  "AT1G56450","AT2G27020", "AT3G60820", "AT4G31300", "AT5G35590", "AT3G51260", "AT3G53230"],
                "cis Golgi proteins": ["AT1G05720", "AT1G07230", "AT1G10950", "AT1G15020", "AT1G18580", "AT1G20270", "AT1G29060", "AT1G29310",
                                      "AT1G51590", "AT1G52420", "AT1G53710", "AT1G62330", "AT1G65820", "AT1G76270", "AT1G77370", "AT1G78920",
                                      "AT2G01070", "AT2G14740", "AT2G17720", "AT2G20130", "AT2G20810", "AT2G40280", "AT2G43080", "AT2G47320",
                                      "AT3G06300", "AT3G09090", "AT3G21160", "AT3G24160", "AT3G28480", "AT3G48280", "AT4G01210", "AT4G24530",
                                      "AT5G04480", "AT5G06660", "AT5G14430", "AT5G14950", "AT5G18900", "AT5G27330", "AT5G47780", "AT5G65470",
                                      "AT5G66060"],
                "photosystem": ["AT2G33040", "AT5G13450", "ATCG00280", "AT1G31330", "AT1G29920", "ATCG00340", "ATCG00350", "AT1G61520",
                               "ATCG00580", "ATCG00710", "AT4G12800", "AT4G10340", "AT3G08940", "AT3G54890", "ATCG00270", "ATCG00020",
                               "AT5G13440", "AT1G55670", "AT4G22890", "AT3G47470", "AT1G45474"]
            },
            "Mouse - Swissprot" :
            {
                "STH" : ["STH"]
            },
        }

    def __init__(self, **kwargs):
        
        # df_original contains all information of the raw file; tab separated file is imported,
        self.filename = "6_deep_maps.txt" if "filename" not in kwargs.keys() else kwargs["filename"]
        
        self.name_pattern = ".* (?P<cond>.*)_(?P<rep>.*)_(?P<frac>.*)" if "name_pattern" not in kwargs.keys() else kwargs["name_pattern"]
        
        #self.map_of_interest = "MAP1" if "map_of_interest" not in kwargs.keys() else kwargs["map_of_interest"]
        
        self.cluster_of_interest = "Proteasome" if "cluster_of_interest" not in kwargs.keys() else kwargs["cluster_of_interest"]
        
        self.cluster_of_interest_comparison = "Proteasome" if "cluster_of_interest_comparison" not in kwargs.keys() else kwargs["cluster_of_interest_comparison"]

        self.summed_MSMS_counts = 2 if "summed_MSMS_counts" not in kwargs.keys() else kwargs["summed_MSMS_counts"]
        self.consecutiveLFQi = 4 if "consecutiveLFQi" not in kwargs.keys() else kwargs["consecutiveLFQi"]
        
        self.RatioHLcount = 2 if "RatioHLcount" not in kwargs.keys() else kwargs["RatioHLcount"]
        self.RatioVariability = 30 if "RatioVariability" not in kwargs.keys() else kwargs["RatioVariability"]
            
        self.collapse_maps = False if "collapse_maps" not in kwargs.keys() else kwargs["collapse_maps"]
        self.collapse_cluster = False if "collapse_cluster" not in kwargs.keys() else kwargs["collapse_cluster"]
        #self.collapse_maps_PCA = False if "collapse_maps_PCA" not in kwargs.keys() else kwargs["collapse_maps_PCA"]
        self.markerset_or_cluster = False if "markerset_or_cluster" not in kwargs.keys() else kwargs["markerset_or_cluster"]
        
        self.clusters_for_ranking = ["x", "y"] if "clusters_for_ranking" not in kwargs.keys() else kwargs["clusters_for_ranking"]
        
        self.multi_choice = ["x", "y"] if "multi_choice" not in kwargs.keys() else kwargs["multi_choice"]
        self.multi_choice_venn = ["x", "y"] if "multi_choice_venn" not in kwargs.keys() else kwargs["multi_choice_venn"]
        
        self.expname = "Experiment_name" if "expname" not in kwargs.keys() else kwargs["expname"]      
        
        #self.x_PCA = "PC1" if "x_PCA" not in kwargs.keys() else kwargs["x_PCA"]
        #self.y_PCA = "PC3" if "y_PCA" not in kwargs.keys() else kwargs["y_PCA"]
        
        self.x_PCA_comp = "PC1" if "x_PCA_comp" not in kwargs.keys() else kwargs["x_PCA_comp"]
        self.y_PCA_comp = "PC3" if "y_PCA_comp" not in kwargs.keys() else kwargs["y_PCA_comp"]
                
        self.acquisition = "SILAC" if "acquisition" not in kwargs.keys() else kwargs["acquisition"]
        self.ref_exp = "Exp_name" if "ref_exp" not in kwargs.keys() else kwargs["ref_exp"]
        
        self.analysed_datasets_dict = {}
        self.analysis_summary_dict = {}
        
        self.markerproteins = self.markerproteins["Human - Swissprot"] if "organism" not in kwargs.keys() else self.markerproteins[kwargs["organism"]]

        
    def data_reading(self, filename=None, content=None):
        """
        Data import. Can read the df_original from a file or buffer.

        Args:
            self:
                filename: string
                regex["imported_columns"] : dictionry; columns that correspond to this regular expression will be imported
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

        self.df_original = pd.read_csv(content, sep="\t", comment="#", usecols=lambda x: bool(re.match(self.regex["imported_columns"], x)))
        
        self.filename = filename

        return self.df_original
    

    def processingdf(self):
        """
        Analysis of the SILAC/LFQ-MQ/LFQ-Spectronaut data will be performed. The dataframe will be filtered, normalized, and converted into a dataframe, 
        characterized by a flat column index. These tasks is performed by following functions:
            def indexingdf(df_original, acquisition_set_dict, acquisition, fraction_dict, name_pattern)
            def spectronaut_LFQ_indexingdf(df_original, Spectronaut_columnRenaming, acquisition_set_dict, acquisition, fraction_dict, name_pattern)
            def stringency_silac(df_index)
            def normalization_01_silac(df_stringency_mapfracstacked):
            def logarithmization_silac(df_stringency_mapfracstacked):
            def stringency_lfq(df_index):
            def normalization_01_lfq(df_stringency_mapfracstacked):
            def logarithmization_lfq(df_stringency_mapfracstacked):

        Args:
            self.acquisition: string, "SILAC", "LFQ", or "LFQ Spectronaut"

        Returns:
            self:
                map_names: list of Map names
                df_01_stacked: df; 0-1 normalized data with "normalized profile" as column name
                df_log_stacked: df; log transformed data
                analysis_summary_dict["0/1 normalized data - mean"] : 0/1 normalized data across all maps by calculating the mean
                                     ["changes in shape after filtering"]
                                     ["Unique Proteins"] : unique proteins, derived from the first entry of Protein IDs, seperated by a ";"
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
        
        shape_dict = {}
        
        def indexingdf(self):
            """
            For data output from MaxQuant, all columns - except of "MS/MS count" and "LFQ intensity" (LFQ) | "Ratio H/L count", "Ratio H/L variability [%]" 
            (SILAC) - will be set as index. A multiindex will be generated, containing "Set" ("MS/MS count", "LFQ intensity"|  "Ratio H/L count", "Ratio H/L
            variability [%]"), "Fraction" (= defined via "name_pattern") and "Map" (= defined via "name_pattern") as level labels, allowing the stacking and 
            unstacking of the dataframe. The dataframe will be filtered by removing matches to the reverse database, matches only identified by site, and 
            potential contaminants.
            
            Args:
                self:
                    df_original: dataframe, columns defined through self.regex["imported_columns"]
                    acquisition_set_dict: dictionary, all columns will be set as index, except of those that are listed in acquisition_set_dict
                    acquisition: string, "SILAC" or "LFQ"
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
            df_original = df_original.set_index([col for col in df_original.columns if any([re.match(s, col) for s in self.acquisition_set_dict[self.acquisition]]) == False])
    
            # multindex will be generated, by isolating the information about the Map, Fraction and Type from each individual column name
            multiindex = pd.MultiIndex.from_arrays(
                    arrays=[
                        [item for sublist in [[re.findall(s, col)[0] for s in self.acquisition_set_dict[self.acquisition] if re.match(s,col)] 
                                              for col in df_original.columns] for item in sublist],
                        [re.match(self.name_pattern, col).group("rep") for col in df_original.columns] if not "<cond>" in self.name_pattern 
                                              else ["_".join(re.match(self.name_pattern, col).group("cond", "rep")) for col in df_original.columns],
                        [self.fraction_dict[re.match(self.name_pattern, col).group("frac")] for col in df_original.columns],
                    ],
                    names=["Set", "Map", "Fraction"]
            )
            
            df_original.columns = multiindex
            df_original.sort_index(1, inplace=True)
            
            shape_dict["Original size"] = df_original.shape
            
            df_index = df_original.xs(
                    np.nan, 0, "Reverse").xs(
                    np.nan, 0, "Potential contaminant").xs(
                    np.nan, 0, "Only identified by site"
            )
            
            df_index.replace(0, np.nan, inplace=True)
            shape_dict["Shape after categorical filtering"] = df_index.shape

            fraction_wCyt = list(df_index.columns.get_level_values("Fraction").unique())
            
            ##############Cyt should get only be removed if it is not an NMC split
            if "Cyt" in fraction_wCyt and len(fraction_wCyt) >= 4:
                df_index.drop("Cyt", axis=1, level="Fraction", inplace=True)
            else:
                pass
            
            self.fractions = list(df_index.columns.get_level_values("Fraction").unique())
            self.df_index = df_index
            
            return df_index


        def spectronaut_LFQ_indexingdf(self):
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
                    df_original: dataframe, columns defined through self.regex["imported_columns"]
                    Spectronaut_columnRenaming
                    acquisition_set_dict: dictionary, all columns will be set as index, except of those that are listed in acquisition_set_dict
                    acquisition: string, "LFQ Spectronaut":
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
            
            df_renamed["Fraction"] = [re.match(self.name_pattern, i).group("frac") for i in df_renamed["Map"]]
            df_renamed["Map"] = [re.match(self.name_pattern, i).group("rep") for i in df_renamed["Map"]] if not "<cond>" in self.name_pattern else ["_".join(
                        re.match(self.name_pattern, i).group("cond", "rep")) for i in df_renamed["Map"]]
            
            df_index = df_renamed.set_index([col for col in df_renamed.columns if any([re.match(s, col) for s in self.acquisition_set_dict[self.acquisition]])==False])
            
            df_index.columns.names = ["Set"]
            
            df_index = df_index.unstack(["Map", "Fraction"])
            df_index.replace(0, np.nan, inplace=True)
            df_index.rename(columns=self.fraction_dict, inplace=True)
            shape_dict["Original size"]=df_index.shape
            
            fraction_wCyt = list(df_index.columns.get_level_values("Fraction").unique())
            #Cyt is removed only if it is not an NMC split
            if "Cyt" in fraction_wCyt and len(fraction_wCyt) >= 4:
                df_index.drop("Cyt", axis=1, level="Fraction", inplace=True)
            else:
                pass
            
            self.fractions = list(df_index.columns.get_level_values("Fraction").unique())
            
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
                df_eLifeMarkers: df, columns: "Gene names", "Compartment", no index 
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
            #defasult setting: RatioHLcount_1 = 3 ; RatioHLcount_2 = 2 ; RatioVariability = 30
            
            df_countvarfiltered_stacked = df_stack.loc[[count>self.RatioHLcount or (count>=self.RatioHLcount and var<self.RatioVariability) 
                                            for var, count in zip(df_stack["Ratio H/L variability [%]"], df_stack["Ratio H/L count"])]]
            
            shape_dict["Shape after Ratio H/L count (>=3)/var (count>=2, var<30) filtering"] = df_countvarfiltered_stacked.shape

            # "Ratio H/L":normalization to SILAC loading, each individual experiment (FractionXMap) will be divided by its median
            # np.median([...]): only entries, that are not NANs are considered
            df_normsilac_stacked = df_countvarfiltered_stacked["Ratio H/L"].unstack(["Fraction", "Map"]).apply(
                lambda x: x/np.median([el for el in x if not np.isnan(el)]), axis=0).stack(["Map", "Fraction"])

            df_stringency_mapfracstacked = df_countvarfiltered_stacked[["Ratio H/L count", "Ratio H/L variability [%]"]].join(
                pd.DataFrame(df_normsilac_stacked, columns=["Ratio H/L"]))

            # dataframe is grouped (Map, id), that allows the filtering for complete profiles
            df_stringency_mapfracstacked = df_stringency_mapfracstacked.groupby(["Map", "id"]).filter(lambda x: len(x)>=len(self.fractions))
            
            shape_dict["Shape after filtering for complete profiles"]=df_stringency_mapfracstacked.shape
            
            # Ratio H/L is converted into Ratio L/H
            df_stringency_mapfracstacked["Ratio H/L"] = df_stringency_mapfracstacked["Ratio H/L"].transform(lambda x: 1/x)
            
            #Annotation with marker genes
            df_eLifeMarkers = self.df_eLifeMarkers
            
            df_stringency_mapfracstacked.reset_index(inplace=True)
            df_stringency_mapfracstacked = df_stringency_mapfracstacked.merge(df_eLifeMarkers, how="outer", on="Gene names", indicator=True)
            df_stringency_mapfracstacked = df_stringency_mapfracstacked.loc[df_stringency_mapfracstacked["_merge"].isin(["both", 
                                                                                                                         "left_only"])].drop("_merge", axis=1)
            df_stringency_mapfracstacked.set_index([c for c in df_stringency_mapfracstacked.columns
                                                    if c != "Ratio H/L count" and c != "Ratio H/L variability [%]" and c!="Ratio H/L"], inplace=True)
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
                    df_eLifeMarkers: df, columns: "Gene names", "Compartment", no index
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

            # retrieve number of fractions that are present in the dataset
            df_fractionnumber_stacked = df_index.copy().stack("Fraction")

            df_index = df_index.stack("Map")

            # sorting the level 0, in order to have LFQ intensity -    MS/MS count instead of continuous alternation
            df_index.sort_index(axis=1, level=0, inplace=True)
            
            # "MS/MS count"-column: take the sum over the fractions; if the sum is larger than n[fraction]*2, it will be stored in the new dataframe
            df_mscount_mapstacked = df_index.loc[df_index[("MS/MS count")].apply(np.sum, axis=1) >= (len(self.fractions) * self.summed_MSMS_counts)]

            shape_dict["Shape after MS/MS value filtering"]=df_mscount_mapstacked.shape
            
            df_stringency_mapfracstacked = df_mscount_mapstacked.copy()

            # series no dataframe is generated; if there are at least i.e. 4 consecutive non-NANs, data will be retained
            df_stringency_mapfracstacked = df_stringency_mapfracstacked.loc[
                df_stringency_mapfracstacked[("LFQ intensity")].apply(lambda x: any(
                    np.invert(np.isnan(x)).rolling(window=self.consecutiveLFQi).sum() >=
                    self.consecutiveLFQi), axis=1)]
            
            shape_dict["Shape after consecutive value filtering"]=df_stringency_mapfracstacked.shape

            df_stringency_mapfracstacked = df_stringency_mapfracstacked.copy().stack("Fraction")
            
            #Annotation with marker genes
            df_eLifeMarkers = self.df_eLifeMarkers
            
            df_stringency_mapfracstacked.reset_index(inplace=True)
            df_stringency_mapfracstacked = df_stringency_mapfracstacked.merge(df_eLifeMarkers, 
                                                                              how="outer", on="Gene names", indicator=True)
            df_stringency_mapfracstacked = df_stringency_mapfracstacked.loc[df_stringency_mapfracstacked["_merge"].isin(["both", 
                                                                                                                         "left_only"])].drop("_merge", axis=1)
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
                
                self.df_01_stacked_for_quantity: see df_01_stacked, w/o imputation

            """

            df_01norm_mapstacked = df_stringency_mapfracstacked["LFQ intensity"].unstack("Fraction")

            # 0:1 normalization of Ratio L/H
            df_01norm_unstacked = df_01norm_mapstacked.div(df_01norm_mapstacked.sum(axis=1), axis=0)

            df_01_stacked = df_stringency_mapfracstacked[["MS/MS count"]].join(pd.DataFrame(df_01norm_unstacked.stack(
                   "Fraction"),columns=["LFQ intensity"]))

            # rename columns: "LFQ intensity" into "normalized profile"
            df_01_stacked.columns = [col if col!="LFQ intensity" else "normalized profile" for col in
                                     df_01_stacked.columns]
            
            self.df_01_stacked_for_quantity = df_01_stacked
            
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

            df_log_stacked = df_stringency_mapfracstacked[["MS/MS count"]].join(pd.DataFrame(df_lognorm_ratio_stacked, columns=["LFQ intensity"]))
            
            # "LFQ intensity" will be renamed to "log profile"
            df_log_stacked.columns = [col if col!="LFQ intensity" else "log profile" for col in df_log_stacked.columns]

            return df_log_stacked


        if self.acquisition == "SILAC":
            df_index = indexingdf(self)#self.df_original, self.acquisition_set_dict, self.acquisition, self.fraction_dict, self.name_pattern)
            
            map_names = df_index.columns.get_level_values("Map").unique()
            self.map_names = map_names
            
            df_stringency_mapfracstacked = stringency_silac(df_index)
            self.df_01_stacked = normalization_01_silac(df_stringency_mapfracstacked)
            self.df_log_stacked = logarithmization_silac(df_stringency_mapfracstacked)
            
            self.analysis_summary_dict["0/1 normalized data - mean"] = self.df_01_stacked["normalized profile"].unstack("Map").dropna().mean(axis=1).to_frame(
                name="normalized profile - mean").reset_index().to_json()

            unique_proteins = list(dict.fromkeys([i.split(";")[0] for i in self.df_01_stacked.reset_index()["Protein IDs"]]))
            self.analysis_summary_dict["Unique Proteins"] = unique_proteins
            self.analysis_summary_dict["changes in shape after filtering"] = shape_dict.copy() 
            analysis_parameters = {"acquisition" : self.acquisition, 
                                   "filename" : self.filename,
                                   "Ratio H/L count 1 (>=X)" : self.RatioHLcount_1,
                                   "Ratio H/L count 2 (>=Y, var<Z)" : self.RatioHLcount_2,
                                   "Ratio variability (<Z, count>=Y)" : self.RatioVariability
                                  }
            self.analysis_summary_dict["Analysis parameters"] = analysis_parameters.copy() 
            self.analysed_datasets_dict[self.expname] = self.analysis_summary_dict.copy()
            #return self.df_01_stacked


        elif self.acquisition == "LFQ" or self.acquisition == "LFQ Spectronaut":
            if self.acquisition == "LFQ":
                df_index = indexingdf(self)
            elif self.acquisition == "LFQ Spectronaut":
                df_index = spectronaut_LFQ_indexingdf(self)
            
            map_names = df_index.columns.get_level_values("Map").unique()
            self.map_names = map_names
            
            df_stringency_mapfracstacked = stringency_lfq(df_index)
            self.df_log_stacked = logarithmization_lfq(df_stringency_mapfracstacked)
            self.df_01_stacked = normalization_01_lfq(df_stringency_mapfracstacked)
            
            self.analysis_summary_dict["0/1 normalized data - mean"] = self.df_01_stacked["normalized profile"].unstack("Map").dropna().mean(axis=1).to_frame(
                name="normalized profile - mean").reset_index().to_json() 
            unique_proteins = list(dict.fromkeys([i.split(";")[0] for i in self.df_01_stacked.reset_index()["Protein IDs"]]))
            self.analysis_summary_dict["Unique Proteins"] = unique_proteins
            self.analysis_summary_dict["changes in shape after filtering"] = shape_dict.copy() 
            analysis_parameters = {"acquisition" : self.acquisition, 
                                   "filename" : self.filename,
                                   "consecutive data points" : self.consecutiveLFQi,
                                   "summed MS/MS counts" : self.summed_MSMS_counts
                                  }
            self.analysis_summary_dict["Analysis parameters"] = analysis_parameters.copy() 
            self.analysed_datasets_dict[self.expname] = self.analysis_summary_dict.copy()
            #return self.df_01_stacked

        else:
            return "I do not know this"    

            
    def quantity_profiles_proteinGroups(self):
        """
        Number of profiles, protein groups per experiment, and the data completness of profiles (total quantity, intersection) is calculated.
        
        Args:
            self:
                acquisition: string, "SILAC", "LFQ", or "LFQ Spectronaut"
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
    
        if self.acquisition == "SILAC":
            df_index = self.df_index["Ratio H/L"]
            df_01_stacked = self.df_01_stacked["normalized profile"]
        elif self.acquisition == "LFQ" or self.acquisition == "LFQ Spectronaut":
            df_index = self.df_index["LFQ intensity"]
            df_01_stacked = self.df_01_stacked_for_quantity["normalized profile"]
        
        #unfiltered
        npg_t = df_index.shape[0]
        df_index_MapStacked = df_index.stack("Map")
        npr_t = df_index_MapStacked.shape[0]
        npr_t_dc = 1-df_index_MapStacked.isna().sum().sum()/(df_index_MapStacked.shape[0]*df_index_MapStacked.shape[1])
        
        #filtered
        df_01_stacked.unstack(["Map", "Fraction"]).sort_index(axis=1)
        npgf_t = df_01_stacked.unstack(["Map", "Fraction"]).shape[0]
        df_01_MapStacked = df_01_stacked.unstack("Fraction")
        nprf_t = df_01_MapStacked.shape[0]
        nprf_t_dc = 1-df_01_MapStacked.isna().sum().sum()/(df_01_MapStacked.shape[0]*df_01_MapStacked.shape[1])
        
        #unfiltered
        df_index_intersection = df_index_MapStacked.groupby(level="Protein IDs").filter(lambda x : len(x)==len(self.map_names))
        npr_i = df_index_intersection.shape[0]
        npr_i_dc = 1-df_index_intersection.isna().sum().sum()/(df_index_intersection.shape[0]*df_index_intersection.shape[1])
        npg_i = df_index_intersection.unstack("Map").shape[0]
        
        #filtered
        df_01_intersection = df_01_MapStacked.groupby(level = "Protein IDs").filter(lambda x : len(x)==len(self.map_names))
        nprf_i = df_01_intersection.shape[0]
        nprf_i_dc = 1-df_01_intersection.isna().sum().sum()/(df_01_intersection.shape[0]*df_01_intersection.shape[1])
        npgf_i = df_01_intersection.unstack("Map").shape[0]
        
                
        dict_npgf = {}
        dict_npgf_dc = {}
        dict_npg = {}
        dict_npg_dc = {}
        for df_intersection in [df_index_intersection, df_01_intersection]:
            for fraction in i_class.fractions:
                df_intersection_frac = df_intersection[fraction]
                npgF_f_dc = 1-df_intersection_frac.isna().sum()/len(df_intersection_frac)
                npgF_f = df_intersection_frac.unstack("Map").isnull().sum(axis=1).value_counts()
                if df_intersection.shape==df_index_intersection.shape:
                    dict_npg[fraction] = npgF_f
                    dict_npg_dc[fraction] = [npgF_f_dc]
                else:
                    dict_npgf[fraction] = npgF_f
                    dict_npgf_dc[fraction] = [npgF_f_dc]
                    
                    
        df_npg = pd.DataFrame(dict_npg)
        df_npg.index.name =  "Number of NaN"
        df_npg.rename_axis("Fraction", axis=1, inplace=True)
        df_npg = df_npg.stack("Fraction").reset_index()
        self.df_npg = df_npg.rename({0: "Protein Groups"}, axis=1)
        
        df_npg_dc = pd.DataFrame(dict_npgf_dc)
        df_npg_dc.rename_axis("Fraction", axis=1, inplace=True)
        df_npg_dc = df_npg_dc.rename({0:"Data completeness"}, axis="index").T
        self.df_npg_dc = df_npg_dc.reset_index()
                    
        df_npgf = pd.DataFrame(dict_npgf)
        df_npgf.index.name =  "Number of NaN"
        df_npgf.rename_axis("Fraction", axis=1, inplace=True)
        df_npgf = df_npgf.stack("Fraction").reset_index()
        self.df_npgf = df_npgf.rename({0: "Protein Groups"}, axis=1)
        
        df_npgf_dc = pd.DataFrame(dict_npgf_dc)
        df_npgf_dc.rename_axis("Fraction", axis=1, inplace=True)
        df_npgf_dc = df_npgf_dc.rename({0:"Data completeness"}, axis="index").T
        self.df_npgf_dc = df_npgf_dc.reset_index()
        
        
        df_quantity_pr_pg = pd.DataFrame(np.array([["before filtering", "total", npg_t, npr_t/len(self.map_names), npr_t_dc],
                                                   ["before filtering", "intersection", npg_i, npr_i/len(self.map_names), npr_i_dc],
                                                   ["after filtering", "total", npgf_t, nprf_t/len(self.map_names), nprf_t_dc],
                                                   ["after filtering", "intersection", npgf_i, nprf_i/len(self.map_names), nprf_i_dc]]), 
                                         columns=["filtering", "type", "number of protein groups", "number of profiles", "data completeness of profiles"])
        
        self.df_quantity_pr_pg = df_quantity_pr_pg.reset_index()
        self.analysis_summary_dict["quantity: profiles/protein groups"] = self.df_quantity_pr_pg.to_json() 
            
            
    def plot_quantity_profiles_proteinGroups(self):
        """
        
        Args:
            self:
                df_quantity_pr_pg: df; no index, columns: "filtering", "type", "npg", "npr", "npr_dc"; further information: see above
                
        Returns:
            
        
        """
        df_quantity_pr_pg = self.df_quantity_pr_pg
        
        fig_npg = go.Figure()
        for t in df_quantity_pr_pg["type"].unique():
            plot_df = df_quantity_pr_pg[df_quantity_pr_pg["type"] == t]
            fig_npg.add_trace(go.Bar(
                x=plot_df["filtering"],
                y=plot_df["number of protein groups"],
                name=t))
        fig_npg.update_layout(barmode="overlay", title="Number of Protein Groups", autosize=False, width=300, height=500)
                
        fig_npr = go.Figure()
        for t in df_quantity_pr_pg["type"].unique():
            plot_df = df_quantity_pr_pg[df_quantity_pr_pg["type"] == t]
            fig_npr.add_trace(go.Bar(
                x=plot_df["filtering"],
                y=plot_df["number of profiles"],
                name=t))
        fig_npr.update_layout(barmode="overlay", title="Number of Profiles", autosize=False, width=300, height=500)
        
        df_quantity_pr_pg = df_quantity_pr_pg.sort_values("filtering")
        fig_npr_dc = go.Figure()
        for t in df_quantity_pr_pg["filtering"].unique():
            plot_df = df_quantity_pr_pg[df_quantity_pr_pg["filtering"] == t]
            fig_npr_dc.add_trace(go.Bar(
                x=plot_df["type"],
                y=plot_df["data completeness of profiles"],
                name=t))
        fig_npr_dc.update_layout(barmode="overlay", title="Coverage", autosize=False, width=300, height=500)
        
        return pn.Row(pn.Column(fig_npg), pn.Column(fig_npr), pn.Column(fig_npr_dc)) 
                                                                              
                                                                                
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

        if self.acquisition == "SILAC":
            df_01orlog_fracunstacked = self.df_log_stacked["log profile"].unstack("Fraction").dropna()
            df_01orlog_MapFracUnstacked = self.df_log_stacked["log profile"].unstack(["Fraction", "Map"]).dropna()  
            
        elif self.acquisition == "LFQ" or self.acquisition == "LFQ Spectronaut":
            df_01orlog_fracunstacked = self.df_01_stacked["normalized profile"].unstack("Fraction").dropna()
            df_01orlog_MapFracUnstacked = self.df_01_stacked["normalized profile"].unstack(["Fraction", "Map"]).dropna()  
            
        pca = PCA(n_components=3)

        # df_pca: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3"
        df_pca = pd.DataFrame(pca.fit_transform(df_01orlog_fracunstacked))
        df_pca.columns = ["PC1", "PC2", "PC3"]
        df_pca.index = df_01orlog_fracunstacked.index
        self.df_pca = df_pca
        
        # df_pca: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3"
        df_pca_combined = pd.DataFrame(pca.fit_transform(df_01orlog_MapFracUnstacked))
        df_pca_combined.columns = ["PC1", "PC2", "PC3"]
        df_pca_combined.index = df_01orlog_MapFracUnstacked.index
        self.df_pca_combined = df_pca_combined
        
        map_names = self.map_names
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
        df_pca_all_marker_cluster_maps = df_pca_all_marker_cluster_maps_unfiltered.groupby(["Gene names"]).filter(lambda x: len(x) >= len(map_names))
        self.df_pca_all_marker_cluster_maps = df_pca_all_marker_cluster_maps

        
    def global_pca_plot(self, map_of_interest, cluster_of_interest, x_PCA="PC1", y_PCA="PC3", collapse_maps=False):
        """"
        PCA plot will be generated

        Args:
            self:
                df_eLifeMarkers: df, columns: "Gene names", "Compartment", no index
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
            df_global_pca.loc[df_global_pca["Gene names"] == i, "Compartment"] = "Selection"

        compartments = self.df_eLifeMarkers["Compartment"].unique()
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
                                    opacity=0.9
                                    )
        return fig_global_pca                         
            
            
    def pca_plot(self, cluster_of_interest):
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
                    if marker not in df_pca_all_marker_cluster_maps.index.get_level_values("Gene names"):
                        continue
                    plot_try_pca = df_pca_all_marker_cluster_maps.xs((marker, maps), level=["Gene names", "Map"],
                                                                     drop_level=False)
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
                                  title="PCA plot for <br>the protein cluster: {}".format(cluster_of_interest))
            return pca_figure
        
        except:
            return "This protein cluster was not quantified"
            

    def profiles_plot(self, map_of_interest, cluster_of_interest):
        """
        The function allows the plotting of filtered and normalized spatial proteomic data using plotly.express.
        The median profile is also calculated and displayed

        Args:
            df_ setofproteins: multiindex dataframe, that contains data about the desired protein cluster,
            stored in the columns "Ratio H/L count", "Ratio H/L variability [%]" and "normalized profile"

        Returns:
            abundance_profiles_and_median_figure: Line plot, displaying the relative abundance profiles.
        """
        
        try:
            df_setofproteins = self.cluster_isolation_df(map_of_interest, cluster_of_interest)
            
            df_setofproteins = df_setofproteins.copy()
    
            # fractions get sorted
            df_setofproteins = df_setofproteins.reindex(index=natsort.natsorted(df_setofproteins.index))
    
            df_setofproteins_median = df_setofproteins["normalized profile"].unstack("Fraction").median()
    
            # make it available for plotting
            df_setofproteins.reset_index(inplace=True)
            abundance_profiles_figure = px.line(df_setofproteins, 
                                                x="Fraction", 
                                                y="normalized profile",
                                                color="Gene names",
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
    
            return abundance_profiles_and_median_figure
        
        except:
            return "This protein cluster was not quantified"
            


    def multiple_iterations(self):
        """
        For each individual protein cluster and for each indiviual fraction the distance to the median will be calculated
        and stored in a multiindex dataframe.

        Args:
            self:
             requires dataframe (df_01_stacked, single level column), stored as attribute (self.df_01_stacked),
             in which "MAP" and "Fraction" are stacked. Additionally the columns "MS/MS count" and
             "Ratio H/L count | Ratio H/L variability [%] | Ratio H/L" are found in LFQ and SILAC data respectively.
             markerproteins: dictionary, key: cluster name, value: gene names (e.g. {"Proteasome" : ["PSMA1", "PSMA2",...], ...}


        Returns:
            df_allclusters_onlynorm_fracunstacked: multiindex dataframe, that contains only the normalized data of individual
            protein clusters substracted by the median of the respective protein cluster for each fraction. The fractions are
            unstacked.
        """

        df_allclusters_onlynorm_fracunstacked_unfiltered = pd.DataFrame()
        df_allclusters_01 = pd.DataFrame()
        
        # for each individual map, and each individual protein cluster, defined in the dictionary markerproteins,
        # the functions "cluster_isolation_df" and "distance_to_median_calculation" will be performed
        for maps in self.map_names:
            for clusters in self.markerproteins:
                #if a certain cluster is not available in the dataset at all
                try:
                    df_setofproteins = self.cluster_isolation_df(maps, clusters)
                    df_allclusters_01 = df_allclusters_01.append(df_setofproteins)
                    df_distance_to_median_fracunstacked = self.distance_to_median_calculation(df_setofproteins)
                    # new column is introduced: Column name = "Cluster"; values: clustername, to wich the individual protein belongs to
                    df_distance_to_median_fracunstacked["Cluster"] = clusters
                    df_distance_to_median_fracunstacked.set_index("Cluster", inplace=True, append=True)
                    
                except:
                    continue
                df_allclusters_onlynorm_fracunstacked_unfiltered = df_allclusters_onlynorm_fracunstacked_unfiltered.append(
                    df_distance_to_median_fracunstacked)
        
        self.df_allclusters_01_test = df_allclusters_01
        #storage of 0/1 normalized data in global dictionary
        self.analysis_summary_dict["0/1 normalized data"] = df_allclusters_01.reset_index().to_json() 
        # genes are droped, if they are not present in all maps
        
        self.df_allclusters_onlynorm_fracunstacked_unfiltered = df_allclusters_onlynorm_fracunstacked_unfiltered
        
        df_allclusters_onlynorm_fracunstacked = df_allclusters_onlynorm_fracunstacked_unfiltered.groupby(["Gene names"]).filter(lambda x: len(x) >= len(self.map_names))
        self.df_allclusters_onlynorm_fracunstacked = df_allclusters_onlynorm_fracunstacked
        
        # concatenate both original and filtered dataframe, and drop duplicates
        dfs_dictionary = {"DF1": df_allclusters_onlynorm_fracunstacked_unfiltered,
                          "DF2": df_allclusters_onlynorm_fracunstacked}
        df = pd.concat(dfs_dictionary)
        df_genenames_sortedout = df.drop_duplicates(keep=False)
        
        # retrieve all gene names, that are sorted out
        genenames_sortedout_index = df_genenames_sortedout.index.get_level_values("Gene names").unique()
        genenames_sortedout_list = genenames_sortedout_index.tolist()

        self.genenames_sortedout_list = genenames_sortedout_list


    def quantification_overview(self, cluster_of_interest):
        """
        
        Args:
            self.df_allclusters_onlynorm_fracunstacked_unfiltered 
                columns: 01K, 03K, 06K, 12K, 24K, 80K
                index: Gene names, Protein IDs, C-Score, Q-value, Map, Compartment, Cluster
            
        Returns:
            df
        """
        
        df_quantification_overview = self.df_allclusters_onlynorm_fracunstacked_unfiltered.xs(cluster_of_interest, 
                                                                 level="Cluster").unstack("Map")[self.df_allclusters_onlynorm_fracunstacked_unfiltered.columns[0]]
        df_quantification_overview = df_quantification_overview.droplevel([i for i in df_quantification_overview.index.names if not i=="Gene names"])
        df_quantification_overview = df_quantification_overview.notnull().replace({True: "x", False: "-"})
        
        return df_quantification_overview
        
        
    def cluster_isolation_df(self, maps, clusters):
        """
        The function returns a dataframe for a specified proteincluster and map, that will be hand over to the
        definition "distance_to_median_calculation"

        Args:
            maps: specified map (e.g. "MAP1")
            clusters: protein cluster, specified in self.markerproteins (e.g. "Proteasome")
            self:
                markerproteins: dictionary, key: cluster name, value: gene names (e.g. {"Proteasome" : ["PSMA1", "PSMA2",...], ...}
                
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
        df_setofproteins_allMaps = pd.DataFrame()
        df_01_stacked = self.df_01_stacked

        # datapoints of each individual markerprotein is written into plot_try and appended to df_setofproteins
        for marker in self.markerproteins[clusters]:
            if marker not in df_01_stacked.index.get_level_values("Gene names"):
                continue
            df_marker = df_01_stacked.xs((marker, maps), level=["Gene names", "Map"], drop_level=False)
            df_setofproteins = df_setofproteins.append(df_marker)
            
            df_marker_allMaps = df_01_stacked.xs(marker, level="Gene names", drop_level=False)
            df_setofproteins_allMaps = df_setofproteins_allMaps.append(df_marker_allMaps)
            

        return df_setofproteins


    def distance_to_median_calculation(self, df_setofproteins):
        """
        A dataframe will be generated, in which the distances to the median profile will be stored.

        Args:
            df_setofproteins: multiindex dataframe, that contains data about the desired protein cluster,
            stored in the columns "Ratio H/L count", "Ratio H/L variability [%]" and "normalized profile"

        Returns:
            df_distance_to_median_fracunstacked: multiindex dataframe, for each individual protein and fraction the distances
            are calculated. The data is stored in single level columns (Column names: e.g. "03K", "06K", "12K", "24K", "80K")
        """
        # calculate the median of the "normalized profile" of individual proteins over the fractions
        df_setofproteins_median = df_setofproteins["normalized profile"].unstack("Fraction").median()

        # substract the median individually from each entry (e.g. 03K_protein_xy - 03K_median)
        df_distance_to_median_fracunstacked = df_setofproteins["normalized profile"].unstack("Fraction")
        df_distance_to_median_fracunstacked = df_distance_to_median_fracunstacked.apply(lambda x: x - df_setofproteins_median, axis=1)

        return df_distance_to_median_fracunstacked

    
    def distance_boxplot(self, cluster_of_interest):
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
            self.proteins_qunatified_across_all_maps = df_cluster_xmaps_distance_with_index.unstack("Map").shape[0]
        
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
            )
    
            return distance_boxplot_figure
    
        except:
            self.cache_cluster_quantified = False
            

        
    def distance_to_median_boxplot(self, cluster_of_interest):
        """
        A box plot for 1 desired cluster, across all maps and fractions is generated displaying the
        distribution of the distance to the median. For each fraction, one box plot will be displayed.

        Args:
            self:
            df_allclusters_onlynorm_fracunstacked, dataframe with single level column, stored as attribute
            (self.df_allclusters_onlynorm_fracunstacked), in which "Fraction" is unstacked. It contains only the
            normalized data of individual protein clusters substracted by the median of the respective protein cluster
            for each fraction.
            map_names: individual map names are stored as an index

        Returns:
            distance_to_median_boxplot_figure: Box plot. Along the x-axis, the maps are shown, along the y-axis
            the distances is plotted
        """

        map_names = self.map_names
        df_allclusters_onlynorm_fracunstacked = self.df_allclusters_onlynorm_fracunstacked
        
        #storage of "Distances to the median profile" in global dictionary
        df_dist_to_median = abs(df_allclusters_onlynorm_fracunstacked.stack("Fraction"))
        df_dist_to_median.name = "distance"
        df_dist_to_median = df_dist_to_median.reindex(index=natsort.natsorted(df_dist_to_median.index))
        self.analysis_summary_dict["Distances to the median profile"] = df_dist_to_median.reset_index().to_json() 

        
        df_boxplot_manymaps = pd.DataFrame()

        try:
        # for each individual map and a defined cluster data will be extracted from the dataframe
        # "df_allclusters_onlynorm_fracunstacked" and appended to the new dataframe df_boxplot_manymaps
            for maps in map_names:
                plot_try = df_allclusters_onlynorm_fracunstacked.xs((cluster_of_interest, maps), level=["Cluster", "Map"], drop_level=False)
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
                                                       facet_col="Cluster",
                                                       facet_row="Fraction",
                                                       boxmode="overlay", height=300 * 5, width=250 * 4, points="all",
                                                       hover_name="Gene names",
                                                       title="Distribution of the distance to the median for <br>the protein cluster: {}".format(cluster_of_interest))
    
            return distance_to_median_boxplot_figure
            
        except:
            return "This protein cluster was not quantified"


    def distance_calculation(self):
        """
        The distances to the median profile of the individual proteins, belonging to specified clusters, across the
        fractions will be used to calculate the distance profile (e.g. Manhattan distance).

        Args:
            self:
                df_allclusters_onlynorm_fracunstacked: single level column, "Fraction" is unstacked. It contains only the normalized data of individual protein
                                                       clusters substracted by the median of the respective protein cluster for each fraction.

        Returns:
            df_distance_noindex: dataframe, no index. It contains the column name "distance", in which the e.g. Manhattan distances for each individual
                                 protein of the specified clusters (see self.markerproteins) are stored
        """
        
        df_allclusters_onlynorm_fracunstacked = self.df_allclusters_onlynorm_fracunstacked.copy()

        # np.linalg.norm requires array; ord=1: Manhattan distance will be calculated over 5 dimensions (Fractions)
        distance_array = np.linalg.norm(df_allclusters_onlynorm_fracunstacked.to_numpy(), axis=1, ord=1)

        distance_series = pd.Series(distance_array)
        distance_series.name = "distance"
        distance_series.index = df_allclusters_onlynorm_fracunstacked.index
        df_distance_noindex = distance_series.reset_index()

        self.df_distance_noindex = df_distance_noindex
        
        self.analysis_summary_dict["Manhattan distances"] = df_distance_noindex.to_json() 

        
    def dynamic_range(self):
        """
        Dynamic range of each individual protein clusters (of the median profile) across all maps is calculated and displayed"

        Args:
            self:
                markerproteins: dictionary, key: cluster name, value: gene names (e.g. {"Proteasome" : ["PSMA1", "PSMA2",...], ...}
                df_01_stacked: "MAP" and "Fraction" are stacked; the data in the column "normalized profile" is used for plotting. Additionally the columns 
                               "MS/MS count" and "Ratio H/L count | Ratio H/L variability [%] | Ratio H/L" are found in LFQ and SILAC data respectively

        Returns:
            fig_dynamicRange: Bar plot, displaying the dynamic range for each protein cluster
            self.df_dynamicRange: df, no index, columns: "Max", "Min", "Dynamic Range", "Cluster"
        """

        df_setofproteins_allMaps = pd.DataFrame()
        df_dynamicRange = pd.DataFrame()
        df_01_stacked = self.df_01_stacked

        for clusters in self.markerproteins:
            try:
                df_setofproteins_allMaps = pd.DataFrame()
                for marker in self.markerproteins[clusters]:
                    if marker not in df_01_stacked.index.get_level_values("Gene names"):
                        continue
                    df_marker_allMaps = df_01_stacked.xs(marker, level="Gene names", drop_level=False)
                    df_setofproteins_allMaps = df_setofproteins_allMaps.append(df_marker_allMaps)
                df_setofproteins_allMaps_median = df_setofproteins_allMaps["normalized profile"].unstack("Fraction").median()
                
                df_dynamicRange = df_dynamicRange.append(pd.DataFrame(np.array([[max(df_setofproteins_allMaps_median), 
                                                                                 min(df_setofproteins_allMaps_median), 
                                                                                 max(df_setofproteins_allMaps_median)-min(df_setofproteins_allMaps_median),
                                                                                 clusters]]), 
                                                                      columns=["Max", "Min", "Dynamic Range", "Cluster"]))
            except:
                continue
        
        self.analysis_summary_dict["Dynamic Range"] = df_dynamicRange.reset_index(drop=True).to_json()
        
        fig_dynamicRange = px.bar(df_dynamicRange, x="Cluster", y="Dynamic Range", base="Min", width=1200, height=800).update_xaxes(categoryorder="total ascending")
        return fig_dynamicRange
        
    
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
            
        df_overview.set_index(["Cluster", "Map"], inplace=True)
        df_overview.sort_index(axis=0, level=0, inplace=True)

        self.analysis_summary_dict["Overview table"] = df_overview.reset_index().to_json()
        self.analysed_datasets_dict[self.expname] = self.analysis_summary_dict.copy() 
        self.analysis_summary_dict.clear()
        
        return df_overview


    def reframe_df_01ORlog_for_svm(self, df_01ORlog):
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
        df_01ORlog_svm.columns = ["_".join(col) for col in df_01ORlog_svm.columns.values]
        df_01ORlog_svm.rename(index={"undefined" : np.nan}, level="Compartment", inplace=True)
        
        return df_01ORlog_svm

        
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
                {"Experiment name" : {see def read_jsonFile(self) [below]}
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
        
        df_SVM = self.df_SVM
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
            if df_SVM["T: True group"][i]!="True: Nuclear pore complex" and df_SVM["T: True group"][i]!="True: Large Protein Complex" and df_SVM["T: True group"][i]!="True: Actin binding proteins" :
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
                    
        SVM_dict["Summary"]={"Total - Recall": total_recall, 
                           "Membrane - Recall" : membrane_recall, 
                           "Av per organelle - Recall": av_per_organelle_recall,
                           "Median per organelle - Recall" : median_per_organelle_recall,
                           "Av precision organelles" : av_per_organelle_precision,
                           "Av F1 organelles" : avg_organelle_f1,
                           "Av F1 all clusters" :  avg_F1_all_cluster,
                          }
        self.analysed_datasets_dict[self.expname]["Misclassification Analysis"] = SVM_dict.copy()
        SVM_dict.clear()


    def read_jsonFile(self):
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
                "Unique Proteins": list,
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
                "Dynamic Range": df - individual cluster,
                "Overview table" : df - individual cluster,

               ##if user dis the Misclassification Analysis befor downloading the dictionary AnalysedDatasets.json##
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
                df_dynamicRange_combined: df, no index, column names: "Max", "Min", "Dynamic Range", "Cluster", "Experiment"
                unique_proteins_total: dict, key: Experiment name, value: unique protein (groups)
                exp_map_names: list of unique Exp_Map - fusions e.g. LFQ_Map1
                exp_names: list of unique Experiment names - e.g. LFQ
        """
        
        json_dict = self.json_dict
        
        #add experiments that are not stored in AnalysedDAtasets.json for comparison
        try:
            json_dict.update(self.analysed_datasets_dict)
        except:
            pass
    
        acquisition_loaded = []
        unique_proteins_total = {}
        
        for exp_name in json_dict.keys():
            for data_type in json_dict[exp_name].keys():
                if data_type == "0/1 normalized data" and exp_name == list(json_dict.keys())[0]:
                    #convert into dataframe
                    df_01_combined = pd.read_json(json_dict[exp_name][data_type])
                    #get only 01 normalized data 
                    df_01_combined = df_01_combined.set_index(["Fraction", "Map", "Gene names", "Protein IDs", 
                                                               "Compartment"])[["normalized profile"]].unstack(["Fraction", "Map"])
                    df_01_combined.rename(columns = {"normalized profile":exp_name}, inplace=True)
        
                elif data_type == "0/1 normalized data" and exp_name != list(json_dict.keys())[0]:
                    df_01_toadd = pd.read_json(json_dict[exp_name][data_type])
                    df_01_toadd = df_01_toadd.set_index(["Fraction", "Map", "Gene names", "Protein IDs", 
                                                         "Compartment"])[["normalized profile"]].unstack(["Fraction", "Map"])
                    df_01_toadd.rename(columns = {"normalized profile":exp_name}, inplace=True)
                    #dataframes will be concatenated, all proteins/Profiles (also one sided) will be retained
                    df_01_combined = pd.concat([df_01_combined, df_01_toadd], axis=1)#, join="inner")
                    
                elif data_type == "0/1 normalized data - mean" and exp_name == list(json_dict.keys())[0]:
                    df_01_mean_combined = pd.read_json(json_dict[exp_name][data_type])
                    df_01_mean_combined = df_01_mean_combined.set_index(["Fraction", "Gene names", "Protein IDs", 
                                                                         "Compartment"])[["normalized profile - mean"]].unstack(["Fraction"])
                    df_01_mean_combined.rename(columns = {"normalized profile - mean":exp_name}, inplace=True)
        
                elif data_type == "0/1 normalized data - mean" and exp_name != list(json_dict.keys())[0]:
                    df_01_mean_toadd = pd.read_json(json_dict[exp_name][data_type])
                    df_01_mean_toadd = df_01_mean_toadd.set_index(["Fraction", "Gene names", "Protein IDs", 
                                                                   "Compartment"])[["normalized profile - mean"]].unstack(["Fraction"])
                    df_01_mean_toadd.rename(columns = {"normalized profile - mean":exp_name}, inplace=True)
                    df_01_mean_combined = pd.concat([df_01_mean_combined, df_01_mean_toadd], axis=1)#, join="inner")  
                    
                elif data_type == "quantity: profiles/protein groups" and exp_name == list(json_dict.keys())[0]:
                    df_quantity_pr_pg_combined = pd.read_json(json_dict[exp_name][data_type])
                    df_quantity_pr_pg_combined["Experiment"] = exp_name
                    
                elif data_type == "quantity: profiles/protein groups" and exp_name != list(json_dict.keys())[0]:
                    df_quantity_pr_pg_toadd = pd.read_json(json_dict[exp_name][data_type])
                    df_quantity_pr_pg_toadd["Experiment"] = exp_name
                    df_quantity_pr_pg_combined = pd.concat([df_quantity_pr_pg_combined, df_quantity_pr_pg_toadd])  
                    
                elif data_type == "Manhattan distances" and exp_name == list(json_dict.keys())[0]:
                    df_distances_combined = pd.read_json(json_dict[exp_name][data_type])
                    df_distances_combined = df_distances_combined.set_index(["Map", "Gene names", "Cluster", "Protein IDs", 
                                                                             "Compartment"])[["distance"]].unstack(["Map"])
                    df_distances_combined.rename(columns = {"distance":exp_name}, inplace=True)
        
                elif data_type == "Manhattan distances" and exp_name != list(json_dict.keys())[0]:
                    df_distances_toadd = pd.read_json(json_dict[exp_name][data_type])
                    df_distances_toadd = df_distances_toadd.set_index(["Map", "Gene names", "Cluster", "Protein IDs", "Compartment"])[["distance"]].unstack(["Map"])
                    df_distances_toadd.rename(columns = {"distance":exp_name}, inplace=True)
                    df_distances_combined = pd.concat([df_distances_combined, df_distances_toadd], axis=1)#, join="inner")
                
                elif data_type == "Dynamic Range" and exp_name == list(json_dict.keys())[0]:
                    df_dynamicRange_combined = pd.read_json(json_dict[exp_name][data_type])
                    df_dynamicRange_combined["Experiment"] = exp_name
                    
                elif data_type == "Dynamic Range" and exp_name != list(json_dict.keys())[0]:
                    df_dynamicRange_toadd = pd.read_json(json_dict[exp_name][data_type])
                    df_dynamicRange_toadd["Experiment"] = exp_name
                    df_dynamicRange_combined = pd.concat([df_dynamicRange_combined, df_dynamicRange_toadd]) 
                
 #               if data_type == "Overview table" and exp_name == list(json_dict.keys())[0]:
 #                   #convert into dataframe
 #                   df_distanceOverview_combined = pd.read_json(json_dict[exp_name][data_type])
 #                   df_distanceOverview_combined["Experiment"] = exp_name
 #                   df_distanceOverview_combined = df_distanceOverview_combined.set_index(["Map", "Cluster", "Experiment"]).unstack(["Cluster"])
 #       
 #               elif data_type == "Overview table" and exp_name != list(json_dict.keys())[0]:
 #                   df_distanceOverview_toadd = pd.read_json(json_dict[exp_name][data_type])
 #                   df_distanceOverview_toadd["Experiment"] = exp_name
 #                   df_distanceOverview_toadd = df_distanceOverview_toadd.set_index(["Map", "Cluster", "Experiment"]).unstack(["Cluster"])
 #                   #dataframes will be concatenated, only proteins/Profiles that are in both df will be retained
 #                   df_distanceOverview_combined = pd.concat([df_distanceOverview_combined, df_distanceOverview_toadd])
                
                elif data_type == "Unique Proteins":
                    unique_proteins_total[exp_name] = json_dict[exp_name][data_type]
                    
                try:
                    for paramters in json_dict[exp_name][data_type].keys():
                        if paramters=="acquisition":
                            acquisition_loaded.append(json_dict[exp_name][data_type][paramters])
                        #elif parameters=="Non valid profiles":
                except:
                    continue

        #filter for consistently quantified proteins (they have to be in all fractions and all maps)
        df_01_filtered_combined = df_01_combined.dropna()
        df_01_filtered_combined.columns.names = ["Experiment", "Fraction", "Map"]
        #reframe it to make it ready for PCA | dropna: to make sure, that you do consider only fractions that are in all experiments
        df_01_filtered_combined = df_01_filtered_combined.stack(["Experiment", "Map"]).swaplevel(0,1, axis=0).dropna(axis=1)
        index_ExpMap = df_01_filtered_combined.index.get_level_values("Experiment")+"_"+df_01_filtered_combined.index.get_level_values("Map")
        index_ExpMap.name = "Exp_Map"
        df_01_filtered_combined.set_index(index_ExpMap, append=True, inplace=True)
        df_01_filtered_combined = df_01_filtered_combined.div(df_01_filtered_combined.sum(axis=1), axis=0)
        
        #filter for consistently quantified proteins (they have to be in all fractions and all maps)
        #df_01_mean_filtered_combined = df_01_mean_combined.dropna()    
        df_01_mean_combined.columns.names = ["Experiment", "Fraction"]
        #reframe it to make it ready for PCA
        df_01_mean_filtered_combined = df_01_mean_combined.stack(["Experiment"]).dropna(axis=1)
        df_01_mean_filtered_combined = df_01_mean_filtered_combined.div(df_01_mean_filtered_combined.sum(axis=1), axis=0)
        
        df_distances_combined.columns.names = ["Experiment", "Map"]
        series = df_distances_combined.stack(["Experiment", "Map"])
        series.name = "distance"
        
        df_distance_comp = series.to_frame()
        #fuse Experiment and Map into one column = "Exp_Map"
        index_dist_ExpMap = df_distance_comp.index.get_level_values("Experiment")+"_"+df_distance_comp.index.get_level_values("Map")
        index_dist_ExpMap.name = "Exp_Map"
        df_distance_comp.set_index(index_dist_ExpMap, append=True, inplace=True)
        df_distance_comp.reset_index(inplace=True)
        
        self.unique_proteins_total = unique_proteins_total
        self.exp_names = list(df_01_filtered_combined.index.get_level_values("Experiment").unique())
        self.exp_map_names = list(index_dist_ExpMap.unique())
        
        self.df_01_filtered_combined = df_01_filtered_combined 
        self.df_01_mean_filtered_combined = df_01_mean_filtered_combined
        
        self.df_quantity_pr_pg_combined = df_quantity_pr_pg_combined
        self.df_dynamicRange_combined = df_dynamicRange_combined
        
        self.df_distance_comp = df_distance_comp

    
    def perform_pca_comparison(self):
        """
        PCA will be performed, using logarithmized data.

        Args:
            self:
                df_01_filtered_combined: df, which contains 0/1 normalized data across all maps (mean) - for all experiments and for the specified protein clusters
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
                df_global_pca_for_plotting: PCA processed dataframe
                    index: "Gene names", "Protein IDs", "Compartment", "Experiment", 
                    columns: "PC1", "PC2", "PC3"
                    contains all protein IDs, that are consistent throughout all experiments
        """

        markerproteins = self.markerproteins.copy()
        
        df_01_filtered_combined = self.df_01_filtered_combined
        df_01_mean_filtered_combined = self.df_01_mean_filtered_combined 
        
        pca = PCA(n_components=3)

        df_pca = pd.DataFrame(pca.fit_transform(df_01_filtered_combined))
        df_pca.columns = ["PC1", "PC2", "PC3"]
        df_pca.index = df_01_filtered_combined.index
        
        df_global_pca = pd.DataFrame(pca.fit_transform(df_01_mean_filtered_combined))
        df_global_pca.columns = ["PC1", "PC2", "PC3"]
        df_global_pca.index = df_01_mean_filtered_combined.index
        
        try:
            markerproteins["PSMA subunits"] = [item for sublist in [re.findall("PSMA.*",p) for p in markerproteins["Proteasome"]] for item in sublist]
            markerproteins["PSMB subunits"] = [item for sublist in [re.findall("PSMB.*",p) for p in markerproteins["Proteasome"]] for item in sublist]
            del markerproteins["Proteasome"]
        except:
            pass 
        
        
        df_cluster = pd.DataFrame([(k, i) for k, l in markerproteins.items() for i in l], columns=["Cluster", "Gene names"])
        df_global_pca_annotated = df_global_pca.reset_index().merge(df_cluster, how="left", on="Gene names")
        df_global_pca_annotated.Cluster.replace(np.NaN, "Undefined", inplace=True)
        
        self.markerproteins_splitProteasome = markerproteins
        self.df_pca_for_plotting = df_pca
        self.df_global_pca_for_plotting = df_global_pca_annotated
            
            
    def plot_pca_comparison(self):
        """
        A PCA plot for desired experiments (multi_choice) and 1 desired cluster is generated.
        Either the maps for every single experiment are displayed individually or in a combined manner
    
        Args:
            self:
                markerproteins: dictionary, key: cluster name, value: gene names (e.g. {"Proteasome" : ["PSMA1", "PSMA2",...], ...}
                multi_choice: list of experiment names
                collapse_maps: boolean
                cluster_of_interest_comparison: string, protein cluster (key in markerproteins, e.g. "Proteasome")
                df_pca_for_plotting: PCA processed dataframe
                    index: "Experiment", "Gene names", "Map", "Exp_Map"
                    columns: "PC1", "PC2", "PC3"
                    contains only marker genes, that are consistent throughout all maps / experiments
    
        Returns:
            pca_figure: PCA plot for a specified protein cluster.
        """
    
        df_pca_for_plotting = self.df_pca_for_plotting.copy()
        markerproteins = self.markerproteins
        multi_choice = self.multi_choice        

        if self.collapse_maps == False:
            map_or_exp_names= df_pca_for_plotting[df_pca_for_plotting.index.get_level_values("Experiment").isin(multi_choice)].index.get_level_values("Exp_Map").unique()
            level_of_interest = "Exp_Map"
            symbol_pca = "Experiment"     
        else:
            map_or_exp_names = multi_choice
            level_of_interest = "Experiment"
            symbol_pca = None
        try:
            df_setofproteins_PCA = pd.DataFrame()
            for map_or_exp in map_or_exp_names:
                    for marker in markerproteins[self.cluster_of_interest_comparison]:
                        if marker not in df_pca_for_plotting.index.get_level_values("Gene names"):
                            continue
                        plot_try_pca = df_pca_for_plotting.xs((marker, map_or_exp), level=["Gene names", level_of_interest], drop_level=False)
                        df_setofproteins_PCA = df_setofproteins_PCA.append(plot_try_pca)
            df_setofproteins_PCA.reset_index(inplace=True)
            
            df_setofproteins_PCA = df_setofproteins_PCA.assign(Experiment_lexicographic_sort=pd.Categorical(df_setofproteins_PCA["Experiment"], categories=self.sorting_list,
                                                                                                              ordered=True))
            df_setofproteins_PCA.sort_values("Experiment_lexicographic_sort", inplace=True)
                    
            pca_figure = px.scatter_3d(df_setofproteins_PCA, 
                                       x="PC1", 
                                       y="PC2", 
                                       z="PC3", 
                                       color=level_of_interest, 
                                       symbol=symbol_pca, 
                                       hover_data=["Gene names"]
                                      )
                    
            pca_figure.update_layout(autosize=False, 
                                     width=1400, 
                                     height=500,
                                     title="PCA plot for <br>the protein cluster: {}".format(self.cluster_of_interest_comparison)
                                    )
    
            return pca_figure
        except:
            return "This protein cluster was not identified in across all experiments"
    
                
    def plot_global_pca_comparison(self):
        """"
        PCA plot will be generated
    
        Args:
            self:
                df_eLifeMarkers: df, columns: "Gene names", "Compartment", no index
                multi_choice: list of experiment names
                css_color: list of colors
                df_global_pca_for_plotting: PCA processed dataframe
                    index: "Gene names", "Protein IDs", "Compartment", "Experiment", 
                    columns: "PC1", "PC2", "PC3"
                    contains all protein IDs, that are consistent throughout all experiments    
    
        Returns:
            pca_figure: global PCA plot, clusters based on the markerset based (df_eLifeMarkers) are color coded. 
        """
        
        multi_choice = self.multi_choice
        
        df_global_pca_exp = self.df_global_pca_for_plotting.loc[self.df_global_pca_for_plotting["Experiment"].isin(self.multi_choice)]
        df_global_pca_exp.reset_index(inplace=True)

        compartments = list(self.df_eLifeMarkers["Compartment"].unique())
        compartment_color = dict(zip(compartments, self.css_color))
        compartment_color["Selection"] = "black"
        compartment_color["undefined"] = "lightgrey"
        compartments.insert(0, "undefined")
        compartments.insert(len(compartments), "Selection")
            
        cluster = self.markerproteins_splitProteasome.keys()
        cluster_color = dict(zip(cluster, self.css_color))
        cluster_color["Undefined"] = "lightgrey"
                
        
        if self.markerset_or_cluster == True:
            df_global_pca = df_global_pca_exp[df_global_pca_exp.Cluster!="Undefined"].sort_values(by="Cluster")
            df_global_pca = df_global_pca_exp[df_global_pca_exp.Cluster=="Undefined"].append(df_global_pca)
        else:
            for i in self.markerproteins[self.cluster_of_interest_comparison]:
                df_global_pca_exp.loc[df_global_pca_exp["Gene names"] == i, "Compartment"] = "Selection"
            df_global_pca = df_global_pca_exp.assign(Compartment_lexicographic_sort = pd.Categorical(df_global_pca_exp["Compartment"], 
                                                                                                     categories=[x for x in compartments], 
                                                                                                     ordered=True))
            df_global_pca.sort_values(["Compartment_lexicographic_sort", "Experiment"], inplace=True)
            
        fig_global_pca = px.scatter(data_frame=df_global_pca,
                                    x=self.x_PCA_comp,
                                    y=self.y_PCA_comp,
                                    color="Compartment" if self.markerset_or_cluster == False else "Cluster",
                                    color_discrete_map=compartment_color if self.markerset_or_cluster == False else cluster_color,
                                    title="Protein subcellular localization by PCA",
                                    hover_data=["Protein IDs", "Gene names", "Compartment"],
                                    facet_col="Experiment",
                                    facet_col_wrap=2,
                                    opacity=0.9
                                    )
        
#            px.scatter(data_frame=df_marker,
#           x=self.x_PCA_comp,
#           y=i_class.y_PCA_comp,
#           color="Cluster",
#           color_discrete_map=compartment_color,
#           title="Protein subcellular localization by PCA",
#           hover_data=["Protein IDs", "Gene names", "Compartment"],
#           facet_col="Experiment",
#           facet_col_wrap=2,
#           opacity=0.9,
#           height=2000
#           )
        
        
        fig_global_pca.update_layout(autosize=False, 
                                     width=1500 if self.markerset_or_cluster == False else 1600, 
                                     height=500*(int(len(multi_choice) / 2) + (len(multi_choice) % 2 > 0))
                                    )
        
        return fig_global_pca 
    
    
    def distance_boxplot_comparison(self):
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
        
        multi_choice = self.multi_choice
        
        #an error massage, if no Experiments are selected, will be displayed already, that is why: return ""
        if len(multi_choice)>=1:
            pass
        else:
            return ("")
        
        df_distance_comp = self.df_distance_comp.copy()
        #set categroical column, allowing lexicographic sorting
        df_distance_comp["Experiment_lexicographic_sort"] = pd.Categorical(df_distance_comp["Experiment"], categories=self.sorting_list, ordered=True)
        df_distance_comp.sort_values(["Experiment_lexicographic_sort", "Map"], inplace=True)
        
        if self.collapse_maps == False:
            #get only values form experiment of interest
            df_distance_selectedExp = df_distance_comp.loc[df_distance_comp["Experiment"].isin(multi_choice)]
            #get only values form cluster of interest
            df_distance_selectedExp = df_distance_selectedExp.loc[df_distance_selectedExp["Cluster"]==self.cluster_of_interest_comparison]
            
            if df_distance_selectedExp.shape[0] == 0:
                self.cache_cluster_quantified = False
                
            else:
                individual_distance_boxplot_figure=go.Figure()
                for i, exp in enumerate(df_distance_selectedExp["Experiment"].unique()):
                    df_plot=df_distance_selectedExp[df_distance_selectedExp["Experiment"]==exp]
                    individual_distance_boxplot_figure.add_trace(go.Box(
                        x=[df_plot["Experiment"], df_plot["Map"]],
                        y=df_plot["distance"],
                        line=dict(color=px.colors.qualitative.Plotly[i]),
                        boxpoints="all",
                        whiskerwidth=0.2,
                        marker_size=2,
                        name=exp, 
                        hovertext=df_plot["Gene names"]
                    ))
                    
                individual_distance_boxplot_figure.update_layout(boxmode="group", 
                                                                 xaxis_tickangle=0, 
                                                                 title="Manhattan distance distribution for <br>the protein cluster: {}".format(self.cluster_of_interest_comparison),
                                                                 autosize=False,
                                                                 width=350*len(multi_choice),
                                                                 height=500,
                                                                 xaxis=go.layout.XAxis(linecolor="black",
                                                                                       linewidth=1,
                                                                                       title="Map",
                                                                                       mirror=True),
                                                                 yaxis=go.layout.YAxis(linecolor="black",
                                                                                       linewidth=1,
                                                                                       title="distance",
                                                                                       mirror=True))
                
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
                plot_try = df_distance_selectedExp_global.xs((self.cluster_of_interest_comparison, map_or_exp), level=["Cluster", 
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
                                             title="Global Manhattan distance distribution for <br>the protein cluster: {}".format(self.cluster_of_interest_comparison)
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
                                                                        mirror=True)
                                                 )
            
            return distance_boxplot_figure
        
        
    def distance_ranking_barplot_comparison(self):
        """
            For each cluster and for each experiment the median ditance is calculated. The smallest median is set to 1, while the other medians will be 
            normalized accordingly. For the "collapsed view" - investigating the global ranking of an experiment across the cluster, the sum of the normlaized
            medians of each clusters for each experimetn will be calculated and displayed. 
            
        Args:
            self:
                df_distance_comp: no index, column names: "Gene names", "Cluster", "Protein IDs", "Compartment", "Experiment", "Map", "Exp_Map", "distance"
                                 "distance": Manhattan distances for each individual protein of the specified clusters (see self.markerproteins) are stored
                markerproteins: dictionary, key: cluster name, value: gene names (e.g. {"Proteasome" : ["PSMA1", "PSMA2",...], ...}
                clusters_for_ranking: list of clusters, that will be used to calculate the ranking
                multi_choice: list of experiment names
                collapse_maps: boolean
                ref_exp: stirng, reference experiment for normalization
        
        Returns:
            fig_globalRanking: figure, barplot. 
            self:
                df_quantified_cluster: df, cluster as column index e.g. "Proteasome", "Lysosome", "Experiment" as row index. dataframe contains x 
                                       (=quantified in certain experiment) and - (=not quantified in certain experiment)
                sorting_list: list of experiment names. They are sorted in ascending order - from the lowest to the highest "Normalized Median - Sum". 
                              It is used to sort the other dataframes, such that the color coding for each experiment is in all plots the same 
            
        """
        
        multi_choice = self.multi_choice
        #an error massage, if no Experiments are selected, will be displayed already, that is why: return ""
        if len(multi_choice)>=1:
            pass
        else:
            return ("")
        
        df_distance_comp = self.df_distance_comp.copy()
        df_distance_comp = df_distance_comp[df_distance_comp["Experiment"].isin(multi_choice)]
        df_distance_comp = df_distance_comp[df_distance_comp["Cluster"].isin(self.clusters_for_ranking)]

        df_quantified_cluster = df_distance_comp.reset_index()
        df_quantified_cluster = df_quantified_cluster.drop_duplicates(subset=["Cluster", "Experiment"]).set_index(["Cluster", 
                                                                                                                   "Experiment"])["distance"].unstack("Cluster")
        self.df_quantified_cluster = df_quantified_cluster.notnull().replace({True: "x", False: "-"})
        
        #dict_cluster_normalizedMedian = {}
        dict_cluster_normalizedMedian_ref = {}
        dict_median_distance_ranking = {}
        for cluster in self.markerproteins.keys():
            try:
                df_cluster = df_distance_comp[df_distance_comp["Cluster"]==cluster]
                all_median_one_cluster_several_exp = {}
                for exp in multi_choice:
                    median = df_cluster[df_cluster["Experiment"]==exp].median()
                    all_median_one_cluster_several_exp[exp] = float(median)
                    #new
                    if exp == self.ref_exp:
                        ref = median
                dict_median_distance_ranking[cluster] = all_median_one_cluster_several_exp
                min_median = min(all_median_one_cluster_several_exp.items(), key=lambda x: x[1])[1]
                #median_ranking = {exp: median/min_median for exp, median in all_median_one_cluster_several_exp.items()}
                #dict_cluster_normalizedMedian[cluster] = median_ranking
                #new
                median_ranking_ref = {exp: median/ref[0] for exp, median in all_median_one_cluster_several_exp.items()}
                dict_cluster_normalizedMedian_ref[cluster] = median_ranking_ref
            except:
                continue
        
        #df_cluster_normalizedMedian = pd.DataFrame(dict_cluster_normalizedMedian)
        #df_cluster_normalizedMedian.index.name="Experiment"
        #df_cluster_normalizedMedian.rename_axis("Cluster", axis=1, inplace=True)
        #df_ranking = df_cluster_normalizedMedian.stack("Cluster")
        #df_ranking.name="Normalized Median"
        #df_ranking = df_ranking.reset_index()
        #ranking_sum = df_cluster_normalizedMedian.sum(axis=1).round(2)
        #ranking_sum.name = "Normalized Median - Sum"
        #ranking_product = df_cluster_normalizedMedian.product(axis=1).round(2)
        #ranking_product.name = "Normalized Median - Product"
        #df_globalRanking = pd.concat([pd.DataFrame(ranking_sum), pd.DataFrame(ranking_product)], axis=1).reset_index()
        #self.sorting_list = list(df_globalRanking.sort_values("Normalized Median - Sum")["Experiment"])
        #set categroical column, allowing lexicographic sorting
        #df_ranking = df_ranking.assign(Experiment_lexicographic_sort = pd.Categorical(df_ranking["Experiment"], categories=self.sorting_list, ordered=True))
        #df_ranking.sort_values("Experiment_lexicographic_sort", inplace=True)
        #df_globalRanking = df_globalRanking.assign(Experiment_lexicographic_sort = pd.Categorical(df_globalRanking["Experiment"], categories=self.sorting_list,
        #                                                                                          ordered=True))
        #df_globalRanking.sort_values("Experiment_lexicographic_sort", inplace=True)
        #df_ranking["Experiment_lexicographic_sort"] = pd.Categorical(df_ranking["Experiment"], categories=self.sorting_list, ordered=True)
        #df_ranking.sort_values("Experiment_lexicographic_sort", inplace=True)
        #df_globalRanking["Experiment_lexicographic_sort"] = pd.Categorical(df_globalRanking["Experiment"], categories=self.sorting_list, ordered=True)
        #df_globalRanking.sort_values("Experiment_lexicographic_sort", inplace=True)
        
        df_cluster_normalizedMedian_ref = pd.DataFrame(dict_cluster_normalizedMedian_ref)
        df_cluster_normalizedMedian_ref.index.name="Experiment"
        df_cluster_normalizedMedian_ref.rename_axis("Cluster", axis=1, inplace=True)
        
        #median makes a huge differnece, improves result of DIA, MQ, libary
        df_RelDistanceRanking = pd.concat([df_cluster_normalizedMedian_ref.median(axis=1), df_cluster_normalizedMedian_ref.sem(axis=1)], axis=1, 
                                          keys=["Distance Ranking (rel, median)", "SEM"]).reset_index().sort_values("Distance Ranking (rel, median)")
        self.sorting_list = list(df_RelDistanceRanking["Experiment"])
        
        df_cluster_normalizedMedian_ref = df_cluster_normalizedMedian_ref.stack("Cluster")
        df_cluster_normalizedMedian_ref.name="Normalized Median"
        df_cluster_normalizedMedian_ref = df_cluster_normalizedMedian_ref.reset_index()
        
        df_cluster_normalizedMedian_ref = df_cluster_normalizedMedian_ref.assign(Experiment_lexicographic_sort = pd.Categorical(df_cluster_normalizedMedian_ref["Experiment"], categories=self.sorting_list, ordered=True))
        df_cluster_normalizedMedian_ref.sort_values("Experiment_lexicographic_sort", inplace=True)

        
        if self.collapse_cluster == False:
            fig_ranking = px.bar(df_cluster_normalizedMedian_ref, 
                                 x="Cluster", 
                                 y="Normalized Median", 
                                 color="Experiment", 
                                 barmode="group", 
                                 title="Ranking - normalization to reference experiment: {}".format(self.ref_exp)
                                )
            
            fig_ranking.update_xaxes(categoryorder="total ascending")
            fig_ranking.update_layout(autosize=False,
                                      width=1200 if len(multi_choice)<=3 else 300*len(multi_choice),
                                      height=500
                                     )
            
            #fig_ranking = px.bar(df_ranking, 
            #                     x="Cluster", 
            #                     y="Normalized Median", 
            #                     color="Experiment", 
            #                     barmode="group", 
            #                     title="Ranking - normalization to smallest median (=1)"
            #                    )
            #
            #fig_ranking.update_xaxes(categoryorder="total ascending")
            #fig_ranking.update_layout(autosize=False,
            #                          width=1200 if len(multi_choice)<=3 else 300*len(multi_choice),
            #                          height=500
            #                         )
            #
            #return pn.Column(pn.Row(fig_ranking), pn.Row(fig_ranking2))
            return fig_ranking
        
        else:
            #fig_globalRanking = px.bar(df_globalRanking, 
            #                           x="Experiment", 
            #                           y="Normalized Median - Sum", 
            #                           color="Experiment", 
            #                           title="Ranking - sum of all individual normalized medians"
            #                          )
            #fig_globalRanking.update_layout(autosize=False,
            #                                width=250*len(multi_choice),
            #                                height=500,
            #                               )
            fig_globalRanking = px.bar(df_RelDistanceRanking.sort_values("Distance Ranking (rel, median)"), 
                                        x="Experiment",
                                        y="Distance Ranking (rel, median)", 
                                        title="Ranking - median of all individual normalized medians - reference experiment: {}".format(self.ref_exp),
                                        error_x="SEM", error_y="SEM", 
                                        color="Experiment")
            
            fig_globalRanking.update_layout(autosize=False,
                                            width=250*len(multi_choice),
                                            height=500,
                                           )
            #return pn.Column(pn.Row(fig_globalRanking), pn.Row(fig_globalRanking2))
            return fig_globalRanking
    
    def quantity_pr_pg_barplot_comparison(self):
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
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined[df_quantity_pr_pg_combined["Experiment"].isin(self.multi_choice)]
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined.assign(Experiment_lexicographic_sort = pd.Categorical(df_quantity_pr_pg_combined["Experiment"],
                                                                                                                        categories=self.sorting_list, ordered=True))
        #df_quantity_pr_pg_combined.sort_values("Experiment_lexicographic_sort", inplace=True)
        df_quantity_pr_pg_combined.sort_values(["Experiment_lexicographic_sort", "type"], ascending=[True, False], inplace=True)
        #df_quantity_pr_pg_combined.sort_values("type", ascending=False, inplace=True)
        
        layout = go.Layout(barmode="overlay", 
          xaxis_tickangle=0, 
          autosize=False,
          width=350*len(self.multi_choice),
          height=500,
          xaxis=go.layout.XAxis(linecolor="black",
                                linewidth=1,
                                title="Map",
                                mirror=True),
          yaxis=go.layout.YAxis(linecolor="black",
                                linewidth=1,
                                title="distance",
                                mirror=True))
        
        fig_quantity_pg = go.Figure()
        for t in df_quantity_pr_pg_combined["type"].unique():
            plot_df = df_quantity_pr_pg_combined[df_quantity_pr_pg_combined["type"] == t]
            fig_quantity_pg.add_trace(go.Bar(
                x=[plot_df["Experiment"], plot_df["filtering"]],
                y=plot_df["number of protein groups"],
                name=t))
        fig_quantity_pg.update_layout(title="Number of Protein Groups")
        fig_quantity_pg.update_layout(layout)
         
        
        fig_quantity_pr = go.Figure()
        for t in df_quantity_pr_pg_combined["type"].unique():
            plot_df = df_quantity_pr_pg_combined[df_quantity_pr_pg_combined["type"] == t]
            fig_quantity_pr.add_trace(go.Bar(
                x=[plot_df["Experiment"], plot_df["filtering"]],
                y=plot_df["number of profiles"],
                name=t))
        fig_quantity_pr.update_layout(title="Number of Profiles")
        fig_quantity_pr.update_layout(layout) 
        
        return pn.Column(pn.Row(fig_quantity_pg), pn.Row(fig_quantity_pr)) 
    
    
    def coverage_comparison(self):
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
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined[df_quantity_pr_pg_combined["Experiment"].isin(self.multi_choice)].sort_values("filtering")
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined.assign(Experiment_lexicographic_sort = pd.Categorical(df_quantity_pr_pg_combined["Experiment"],
                                                                                                                        categories=self.sorting_list, ordered=True))
        #df_quantity_pr_pg_combined.sort_values("Experiment_lexicographic_sort", inplace=True)
        df_quantity_pr_pg_combined.sort_values(["Experiment_lexicographic_sort", "filtering"], inplace=True)

        
        fig_pr_dc = go.Figure()
        for t in df_quantity_pr_pg_combined["filtering"].unique():
            plot_df = df_quantity_pr_pg_combined[df_quantity_pr_pg_combined["filtering"] == t]
            fig_pr_dc.add_trace(go.Bar(
                x=[plot_df["Experiment"], plot_df["type"]],
                y=plot_df["data completeness of profiles"],
                name=t))
            
        fig_pr_dc.update_layout(barmode="overlay", 
                                         xaxis_tickangle=0, 
                                         title="Data Completeness of Profiles",
                                         autosize=False,
                                         width=350*len(self.multi_choice),
                                         height=500,
                                         xaxis=go.layout.XAxis(linecolor="black",
                                                               linewidth=1,
                                                               title="Map",
                                                               mirror=True),
                                         yaxis=go.layout.YAxis(linecolor="black",
                                                               linewidth=1,
                                                               title="distance",
                                                               mirror=True))
        
        return fig_pr_dc
    
    
    def venn_diagram(self):
        """
        Venn diagram based on 2/3 experiments is created using matplotlib. Figure has to be transforme dfrom matplotlib object to jpg, to make it available for the 
        webinterface via panel/holoviz
        
        Args:
            self:
                multi_choice_venn: list of experiment names, max 3 (slelect widget can only hold max 3 values) 
                unique_proteins_total: dict, key: Experiment name, value: unique protein (groups)
        
        Returns:
            im: Venn diagram, made availabe flor plotly/webinterface
        """
        
        
        multi_choice_venn = self.multi_choice_venn
        
        if len(multi_choice_venn)==2:
            vd = venn2([set(self.unique_proteins_total[i]) for i in multi_choice_venn], 
                       set_labels=([i for i in multi_choice_venn]),
                       set_colors=("darksalmon", "darkgrey")
                      )
        elif len(multi_choice_venn)==3:
            vd = venn3([set(self.unique_proteins_total[i]) for i in multi_choice_venn],
                       set_labels=([i for i in multi_choice_venn]),
                       set_colors=("darksalmon", "darkgrey","rosybrown"),
                       alpha=0.8
                      )
        
        else:
            return ("Please select 2 or more experiments for comparison")

        vd = plt.title("Unique Protein Groups - Venn Diagram",
                       pad=30,
                       fontsize=20
                      )
        
        #make matplot figure available for plotly
        vd = vd.figure
        out_img = io.BytesIO()
        plt.savefig(out_img, bbox_inches="tight",format="jpg", dpi=72)
        out_img.seek(0)  # rewind file
        im = Image.open(out_img)
        plt.clf()
        
        return im
    
    
    def dynamic_range_comparison(self):
        """
        A box plot for desired experiments (multi_choice) and all protein clusters is generated displaying the dynamic range
        
        Args:
            self:
                multi_choice: list of experiment names 
                df_dynamicRange_combined: df, no index, column names: "Max", "Min", "Dynamic Range", "Cluster", "Experiment"
        
        Returns:
            fig_dynamic_range: bar plot, dynamic range of each protein cluster for desired experiments is displayed.
        """
        multi_choice = self.multi_choice
        df_dynamicRange_combined = self.df_dynamicRange_combined.copy()
        df_dynamicRange_combined = df_dynamicRange_combined[df_dynamicRange_combined["Experiment"].isin(multi_choice)]
        df_dynamicRange_combined = df_dynamicRange_combined.assign(Experiment_lexicographic_sort = pd.Categorical(df_dynamicRange_combined["Experiment"],
                                                                                                                categories=self.sorting_list, ordered=True))
        
        df_dynamicRange_combined.sort_values(["Experiment_lexicographic_sort", "Dynamic Range"], inplace=True)
        
        fig_dynamic_range = px.bar(df_dynamicRange_combined, x="Cluster", y="Dynamic Range", base="Min", facet_row="Experiment", height=400*len(multi_choice),
                                   width=1200)
        
        df_dynamicRange_combined_ref = df_dynamicRange_combined.drop(["Experiment_lexicographic_sort"], axis=1)
        df_dynamicRange_combined_ref = df_dynamicRange_combined.set_index(["Cluster", "Experiment"], drop=False).unstack("Cluster")["Dynamic Range"]
        df_dynamicRange_combined_ref = df_dynamicRange_combined_ref.div(df_dynamicRange_combined_ref.xs(self.ref_exp))
        df_RelDynamicRange = pd.concat([df_dynamicRange_combined_ref.median(axis=1), df_dynamicRange_combined_ref.sem(axis=1)], axis=1, 
                                       keys=["Dynamic Range (rel, median)", "SEM"]).reset_index()
        
        if self.collapse_cluster == False:
            df_dynamicRange_combined_ref = df_dynamicRange_combined_ref.stack("Cluster")
            df_dynamicRange_combined_ref.name="Normalized Dynamic Range"
            df_dynamicRange_combined_ref = df_dynamicRange_combined_ref.reset_index()
            
            fig_RelDynamicRange = px.bar(df_dynamicRange_combined_ref, 
                                         x="Cluster", 
                                         y="Normalized Dynamic Range", 
                                         title="Dynamic Range - normalization to reference experiment: {}".format(self.ref_exp),
                                         barmode="group", 
                                         color="Experiment")
            fig_RelDynamicRange.update_xaxes(categoryorder="total ascending")
            fig_RelDynamicRange.update_layout(autosize=False,
                                              width=1200 if len(multi_choice)<=3 else 300*len(multi_choice),
                                              height=500,
                                               )
        else:
            fig_RelDynamicRange = px.bar(df_RelDynamicRange.sort_values("Dynamic Range (rel, median)"), 
                                         x="Experiment", 
                                         y="Dynamic Range (rel, median)", 
                                         error_x="SEM", error_y="SEM",
                                         title="Dynamic Range - median of all individual normalized medians - reference experiment: {}".format(self.ref_exp),
                                         color="Experiment")
            fig_RelDynamicRange.update_layout(autosize=False,
                                                width=250*len(multi_choice),
                                                height=500,
                                               )
            
            
        return pn.Column(pn.Row(fig_dynamic_range), pn.Row(fig_RelDynamicRange)) 
    
    
        
    def __repr__(self):
        return "This is a spatial dataset with {} lines.".format(len(self.df_original))