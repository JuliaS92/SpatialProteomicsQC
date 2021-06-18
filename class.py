class SpatialDataSet:
    
    regex = {
        "imported_columns": "^[Rr]atio H/L (?!normalized|type|is.*).+|id$|[Mm][Ss].*[cC]ount.+$|[Ll][Ff][Qq].*|.*[nN]ames.*|.*[Pp][rR]otein.[Ii][Dd]s.*|[Pp]otential.[cC]ontaminant|[Oo]nly.[iI]dentified.[bB]y.[sS]ite|[Rr]everse|[Ss]core|[Qq]-[Vv]alue|R.Condition|PG.Genes|PG.ProteinGroups|PG.Cscore|PG.Qvalue|PG.RunEvidenceCount|PG.Quantity|^Proteins$|^Sequence$"
    }
    
    acquisition_set_dict = {
        "LFQ6 - Spectronaut" : ["LFQ intensity", "MS/MS count"],
        "LFQ5 - Spectronaut" : ["LFQ intensity", "MS/MS count"],
        "LFQ5 - MQ" : ["[Ll][Ff][Qq].[Ii]ntensity", "[Mm][Ss]/[Mm][Ss].[cC]ount", "[Ii]ntensity"],
        "LFQ6 - MQ" : ["[Ll][Ff][Qq].[Ii]ntensity", "[Mm][Ss]/[Mm][Ss].[cC]ount", "[Ii]ntensity"],
        "SILAC - MQ" : [ "[Rr]atio.[Hh]/[Ll](?!.[Vv]aria|.[Cc]ount)","[Rr]atio.[Hh]/[Ll].[Vv]ariability.\[%\]", "[Rr]atio.[Hh]/[Ll].[cC]ount"]
    }
    
    fraction_dict = {
        "mem_frac1": "03K", "mem_frac2": "06K", "mem_frac3": "12K", "mem_frac4": "24K", "mem_frac5": "80K",
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
                      "fuchsia", "gainsboro", "ghostwhite", "gold", "gray", "ivory", "lavenderblush", "lawngreen", "lemonchiffon", "lightblue", "lightcyan",
                      "lightgoldenrodyellow", "lightgray", "lightgrey", "lightgreen", "lightsalmon", "lightskyblue", "lightslategray", "lightslategrey",
                      "lightsteelblue", "lightyellow", "lime", "limegreen", "linen", "maroon", "mediumaquamarine", "mediumblue", "mediumseagreen",
                      "mediumslateblue", "mediumspringgreen", "mediumturquoise", "mediumvioletred", "midnightblue", "mintcream", "mistyrose", "moccasin",
                      "olivedrab", "orangered", "orchid", "palegoldenrod", "palegreen", "paleturquoise", "palevioletred", "papayawhip", "peachpuff", "peru",
                      "pink", "plum", "powderblue", "rosybrown", "royalblue", "saddlebrown", "salmon", "sandybrown", "seagreen", "seashell", "sienna", "silver",
                      "skyblue", "slateblue", "steelblue", "teal", "thistle", "tomato", "violet", "wheat", "white", "whitesmoke", "slategray", "slategrey",
                      "aquamarine", "azure","crimson", "cyan", "darkslategray", "grey","mediumorchid","navajowhite", "navy"]
    
    markerproteins_set = {
            "Human - Swissprot" :
            {
                "Proteasome" : ["PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7"], 
                        #       "PSMC1", "PSMC2", "PSMC3"],
                "EMC" : ["EMC1", "EMC2", "EMC3", "EMC4", "EMC7", "EMC8", "EMC10","EMC6","EMC9"],
                "Lysosome" : ["LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5", "LAMP1", "LAMP2", "CTSA", "CTSB", "CTSC", "CTSD", "CTSL", "CTSZ"],
                "Arp2/3 protein complex" : ["ACTR2", "ACTR3", "ARPC1B", "ARPC2", "ARPC3", "ARPC4", "ARPC5"],               
                "AP2 adaptor complex" : ["AP2A1", "AP2A2", "AP2B1", "AP2M1",  "AP2S1", ],              
                "Class C, Vps complex" : ["VPS11","VPS16", "VPS18", "VPS33A"],
                "Wave-2 complex": ["NCKAP1", "CYFIP1", "ABI1", "BRK1", "WASF2"],
                "TREX complex / THO complex": ["ALYREF", "THOC5", "THOC2", "THOC1", "THOC3", "DDX39B"], #["THOC5", "THOC2", "THOC1"]
                #"Exon junction complex #TREX: ALYREF,DDX39B": ["SRRM1", "RBM8A", "RNPS1", "EIF4A3", "UPF3B", "UPF2"],
                #"CDC5L complex #SF3b: SF3B4": ["SNRPD1", "SNRPD3", "SNRPD2", "PRPF19", "SRSF1", "SF3B2", "SNRPA1", "SF3B4"],
                "Respiratory chain complex I": ["NDUFB6", "NDUFB10", "NDUFA10", "NDUFA8", "NDUFA6", "NDUFB11", "NDUFB3", "NDUFB5", "NDUFAB1", "NDUFA4", "NDUFB9",
                                                "NDUFB7", "NDUFA9", "NDUFA5", "NDUFV3", "NDUFA11", "NDUFV1", "NDUFA12", "NDUFV2", "NDUFA7", "NDUFS6", "NDUFS2",
                                                "NDUFA2", "NDUFS8", "NDUFS1"],
                "Nup 107-160 subcomplex": ["NUP85", "NUP37", "NUP160", "NUP98", "NUP107", "NUP133"],
                "Multisynthetase complex": ["EEF1E1", "IARS", "DARS", "EPRS", "AIMP1", "KARS", "LARS", "RARS", "AIMP2", "MARS"],       
                "GAA1-GPI8-PIGT-PIG-PIGS complex": ["PIGT", "PIGS", "PIGU", "PIGK", "GPAA1"],
                "Frataxin complex /f1f0: ATP5L": ["SDHA", "HSPD1", "HSPA9", "AFG3L2"],
                "F1F0-ATP synthase": ["ATP5O", "ATP5I", "ATPIF1", "ATP5A1", "ATP5F1", "ATP5B", "ATP5H", "ATP5L", "ATP5J"],
                "Exosome": ["EXOSC1", "EXOSC3", "EXOSC8", "EXOSC4", "EXOSC2", "EXOSC10"],
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
 #CCT and V type were not in set
                "CCT complex" : ["CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8","CCT6B", "TCP1"],
                "V-type proton ATPase": ["ATP6AP1", "ATP6V0A1", "ATP6V0A2", "ATP6V0A4", "ATP6V0D1", "ATP6V1A", "ATP6V1B2", "ATP6V1E1", "ATP6V1G1", "ATP6V1H"],
                "MCM complex" : ["MCM2", "MCM3", "MCM4", "MCM5", "MCM7"],
                "Prefoldin complex" : [ "PFDN1", "PFDN2", "PFDN4", "PFDN5", "PFDN6", "VBP1"],
                "AP1 adaptor complex" : ["AP1B1", "AP1G1", "AP1M1", "AP1S1", "AP1S2", "AP1S3"],
                "AP3 adaptor / AP3-BLOC1 complex" : ["AP3B1", "AP3D1", "AP3M1", "AP3M2", "AP3S1", "AP3S2"],
                "AP4 adaptor complex" : ["AP4B1", "AP4E1","AP4M1",  "AP4S1"],
                "Anaphas,e-promoting complex" : ["ANAPC1", "ANAPC10", "ANAPC16", "ANAPC2", "ANAPC4","ANAPC5", "ANAPC7", "CDC16", "CDC23","CDC27"] ,
                "Rnase/Mrp complex" : ["POP1", "POP4", "POP5", "RPP14","RPP25", "RPP30", "RPP38", "RPP40"],
                "Dynactin complex" : ["DCTN1", "DCTN2", "DCTN3", "DCTN4", "DCTN6", "ACTR1A", "CAPZA1"],
                "CTLH complex" : ["ARMC8", "MAEA", "MKLN1", "RANBP9", "RMND5A"],
                "Coatomer complex" : ["ARCN1", "COPA", "COPB1", "COPB2", "COPE", "COPG1", "COPZ1"],
                "TNF-alpha/NF-kappa B signaling complex 5": ["POLR2H", "POLR1A", "POLR1B", "CUL1", "KPNA2"],
                "Septin complex": ["SEPT7", "SEPT9", "SEPT11", "SEPT8", "SEPT2"],
                "Sec6/8 exocyst complex": ["EXOC4", "EXOC2", "EXOC1", "EXOC7", "EXOC5", "EXOC3", "EXOC8", "EXOC6"],
                "SNW1 complex": ["EFTUD2", "SNRNP200", "PRPF8", "MSH2", "DDX23", "SNW1", "PFKL"],
                #"SF3b complex #SF3B4": ["SF3B1", "SF3B3", "SF3B5", "SF3B6", "PHF5A"],
                "RFC complex": ["RFC4", "RFC2", "RFC5", "RFC3", "RFC1"],
                "MCM complex": ["MCM4", "MCM6", "MCM7", "MCM3", "MCM2", "MCM5"],
                #"Large Drosha complex, DGCR8: FUS,HNRNPH1, DDX17, DDX5": ["HNRNPDL", "RALY", "TARDBP", "HNRNPM", "DDX3X", "EWSR1"],
                "DGCR8 multiprotein complex": ["HNRNPR", "HNRNPH1", "DDX17", "DDX5", "DHX9", "FUS", "NCL"],
                "COP9 signalosome complex / CNS-P53 complex": ["COPS2", "COPS3", "COPS4", "COPS5", "COPS6", "COPS8", "GPS1"],
                #"Autophagy test" : ["ULK1", "ULK2", "ATG13", "FIP200", "ATG101", "VPS34", "VPS15", "BECN1", "ATG14L", "NRBF2", "ATG9", "Rab1", 
                #"TRAPPC8", "MAP1LC3A", "MAP1LC3C", "MAP1LC3B"]
 
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
                "20 SProteasome" : ["Psmb1","Psmb5","Psma3","Psma2","Psmb7","Psmb4","Psmb6","Psma6","Psma4","Psmb3","Psmb2","Psma1","Psma7","Psma5"], 
                "CCT complex": ["Tcp1","Cct8","Cct7","Cct2","Cct4","Cct5","Cct6a","Cct3"],
                "MCM complex": ["Mcm3","Mcm4","Mcm5","Mcm2","Mcm6","Mcm7"],
                "COP9 signalosome complex": ["Cops5","Cops3","Cops4","Cops6","Cops2","Cops8","Gps1","Cops7a"],
                "Anaphase-promoting complex":["Anapc1","Ccdc27","Cdc23","Anapc5","Anapc2","Anapc10","Cdc16","Anapc4","Anapc7"],
                "Caveolar macromolecular signaling complex": ["Adrb2","Prkar2b","Cav3","Gnas","Adcy8","Cacna1c","Ppp2r1a"],
                "(ER)-localized multiprotein complex": ["Pdia4","Hsp90b1","P4hb","Hspa5","Ppib","Uggt1","Dnajb11","Sdf2l1","Hyou1","Cabp1"],
                "Cytochrome bc1-complex, mitochondrial": ["Mt-Cyb,Uqcrh,Uqcr10,Uqcr11,Uqcrq,Uqcrfs1,Uqcrc1,Cyc1,Uqcrb,Uqcrc2"],
                "PYR complex": ["Smarce1","Actb","Hmgb1","Hdac2","Smarcc1","Ikzf1","Rbbp7","Smarcd1","Smarcc2","Chd4","Smarcb1","Actl6a"],
                "Ikaros complex": ["Ikzf3","Hdac1","Hdac2","Ikzf2","Smarcc1","Ikzf1","Smarca4","Rbbp4","Smarcd1","Smarcd3","Chd4","Smarcd2"],
                "Ubiquitin E3 ligase": ["Tceb2","Rbx1","Tceb1","Neurl2","Cul5"],
                "Gata1-Fog1-MeCP1 complex": ["Hdac1","Zfpm1","Gata1","Hdac2","Rbbp4","Rbbp7","Chd4","Mta1","Gatad2b","Mta3","Mta2","Mbd3","Mbd2"],
                "BLOC-1 complex": ["Bloc1s1","Bloc1s3","Bloc1s5","Bloc1s4","Dtnbp1","Bloc1s2","Bloc1s6","Snapin"],

            },
        }

    analysed_datasets_dict = {}
    
    df_organellarMarkerSet = pd.read_csv("eLife_markers.txt", sep="\t", comment="#",
                                       usecols=lambda x: bool(re.match("Gene name|Compartment", x)))
    df_organellarMarkerSet = df_organellarMarkerSet.rename(columns={"Gene name":"Gene names"})
    df_organellarMarkerSet = df_organellarMarkerSet.astype({"Gene names": "str"})

    def __init__(self, filename, expname, acquisition, comment, name_pattern="e.g.:.* (?P<cond>.*)_(?P<rep>.*)_(?P<frac>.*)", **kwargs):
        
        self.filename = filename
        self.expname = expname
        self.acquisition = acquisition
        self.name_pattern = name_pattern
        self.comment =  comment
        
        
        
        self.fractions, self.map_names = [], []
        self.df_01_stacked, self.df_log_stacked = pd.DataFrame(), pd.DataFrame()
        
        if acquisition == "SILAC - MQ":
            if "RatioHLcount" not in kwargs.keys():
                self.RatioHLcount = 2
            else:
                self.RatioHLcount = kwargs["RatioHLcount"]
                del kwargs["RatioHLcount"]
            if "RatioVariability" not in kwargs.keys():
                self.RatioVariability = 30
            else:
                self.RatioVariability = kwargs["RatioVariability"]
                del kwargs["RatioVariability"]
                       
        #elif acquisition == "LFQ5 - MQ" or acquisition == "LFQ6 - MQ" or acquisition == "LFQ6 - Spectronaut" or acquisition == "LFQ5 - Spectronaut":
        else:
            if "summed_MSMS_counts" not in kwargs.keys():
                self.summed_MSMS_counts = 2
            else:
                self.summed_MSMS_counts = kwargs["summed_MSMS_counts"]
                del kwargs["summed_MSMS_counts"]
            if "consecutiveLFQi" not in kwargs.keys():
                self.consecutiveLFQi = 4
            else:
                self.consecutiveLFQi = kwargs["consecutiveLFQi"]
                del kwargs["consecutiveLFQi"]
        
        #self.markerset_or_cluster = False if "markerset_or_cluster" not in kwargs.keys() else kwargs["markerset_or_cluster"]
        if "organism" not in kwargs.keys():
            self.markerproteins = self.markerproteins_set["Human - Swissprot"]
        else:
            assert kwargs["organism"] in self.markerproteins_set.keys()
            self.markerproteins = self.markerproteins_set[kwargs["organism"]]
            del kwargs["organism"]
        
        self.analysed_datasets_dict = {}
        self.analysis_summary_dict = {}
    
    
    def data_reading(self, filename=None, content=None):
        """
        Data import. Can read the df_original from a file or buffer. 
        df_original contains all information of the raw file; tab separated file is imported,

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

        if filename.endswith("xls") or filename.endswith("txt"):
            self.df_original = pd.read_csv(content, sep="\t", comment="#", usecols=lambda x: bool(re.match(self.regex["imported_columns"], x)), low_memory = True)
        else: #assuming csv file
            self.df_original = pd.read_csv(content, sep=",", comment="#", usecols=lambda x: bool(re.match(self.regex["imported_columns"], x)), low_memory = True)
        assert self.df_original.shape[0]>10 and self.df_original.shape[1]>5
        
        self.filename = filename

        return self.df_original
    

    def processingdf(self, name_pattern=None, summed_MSMS_counts=None, consecutiveLFQi=None, RatioHLcount=None, RatioVariability=None):
        """
        Analysis of the SILAC/LFQ-MQ/LFQ-Spectronaut data will be performed. The dataframe will be filtered, normalized, and converted into a dataframe, 
        characterized by a flat column index. These tasks is performed by following functions:
            indexingdf(df_original, acquisition_set_dict, acquisition, fraction_dict, name_pattern)
            spectronaut_LFQ_indexingdf(df_original, Spectronaut_columnRenaming, acquisition_set_dict, acquisition, fraction_dict, name_pattern)
            stringency_silac(df_index)
            normalization_01_silac(df_stringency_mapfracstacked):
            logarithmization_silac(df_stringency_mapfracstacked):
            stringency_lfq(df_index):
            normalization_01_lfq(df_stringency_mapfracstacked):
            logarithmization_lfq(df_stringency_mapfracstacked):

        Args:
            self.acquisition: string, "LFQ6 - Spectronaut", "LFQ5 - Spectronaut", "LFQ5 - MQ", "LFQ6 - MQ", "SILAC - MQ"
            additional arguments can be used to override the value set by the class init function

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
        
        if name_pattern is None:
            name_pattern = self.name_pattern
        
        if self.acquisition == "SILAC - MQ":
            if RatioHLcount is None:
                RatioHLcount = self.RatioHLcount
            if RatioVariability is None:
                RatioVariability = self.RatioVariability
        else:
            if summed_MSMS_counts is None:
                summed_MSMS_counts = self.summed_MSMS_counts
            if consecutiveLFQi is None:
                consecutiveLFQi = self.consecutiveLFQi
            
        shape_dict = {}
        
        
        def indexingdf():
            """
            For data output from MaxQuant, all columns - except of "MS/MS count" and "LFQ intensity" (LFQ) | "Ratio H/L count", "Ratio H/L variability [%]" 
            (SILAC) - will be set as index. A multiindex will be generated, containing "Set" ("MS/MS count", "LFQ intensity"|  "Ratio H/L count", "Ratio H/L
            variability [%]"), "Fraction" (= defined via "name_pattern") and "Map" (= defined via "name_pattern") as level names, allowing the stacking and 
            unstacking of the dataframe. The dataframe will be filtered by removing matches to the reverse database, matches only identified by site, and 
            potential contaminants.
            
            Args:
                self:
                    df_original: dataframe, columns defined through self.regex["imported_columns"]
                    acquisition_set_dict: dictionary, all columns will be set as index, except of those that are listed in acquisition_set_dict
                    acquisition: string, one of "LFQ6 - Spectronaut", "LFQ5 - Spectronaut", "LFQ5 - MQ", "LFQ6 - MQ", "SILAC - MQ"
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
            df_original.rename({"Proteins": "Protein IDs"}, axis=1, inplace=True)
            df_original = df_original.set_index([col for col in df_original.columns
                                                 if any([re.match(s, col) for s in self.acquisition_set_dict[self.acquisition]]) == False])
    
            # multindex will be generated, by extracting the information about the Map, Fraction and Type from each individual column name
            multiindex = pd.MultiIndex.from_arrays(
                arrays=[
                    [[re.findall(s, col)[0] for s in self.acquisition_set_dict[self.acquisition] if re.match(s,col)][0]
                    for col in df_original.columns],
                    [re.match(self.name_pattern, col).group("rep") for col in df_original.columns] if not "<cond>" in self.name_pattern
                     else ["_".join(re.match(self.name_pattern, col).group("cond", "rep")) for col in df_original.columns],
                    [self.fraction_dict[re.match(self.name_pattern, col).group("frac")] for col in df_original.columns],
                ],
                names=["Set", "Map", "Fraction"]
            )
            
            df_original.columns = multiindex
            df_original.sort_index(1, inplace=True)
            
            shape_dict["Original size"] = df_original.shape
            
            try:
                df_index = df_original.xs(
                    np.nan, 0, "Reverse")
            except:
                pass
            
            try:
                df_index = df_index.xs(
                np.nan, 0, "Potential contaminant")
            except:
                pass
            
            try:
                df_index = df_index.xs(
                np.nan, 0, "Only identified by site")
            except:
                pass
            
            df_index.replace(0, np.nan, inplace=True)
            shape_dict["Shape after categorical filtering"] = df_index.shape
            
            df_index.rename(columns={"MS/MS Count":"MS/MS count"}, inplace=True)

            fraction_wCyt = list(df_index.columns.get_level_values("Fraction").unique())
            
            ##############Cyt should get only be removed if it is not an NMC split
            if "Cyt" in fraction_wCyt and len(fraction_wCyt) >= 4:
                df_index.drop("Cyt", axis=1, level="Fraction", inplace=True)
            try:
                if self.acquisition == "LFQ5 - MQ":
                    df_index.drop("01K", axis=1, level="Fraction", inplace=True)
            except:
                pass
            
            self.fractions = list(df_index.columns.get_level_values("Fraction").unique())
            self.df_index = df_index
            
            return df_index
        
        
        def spectronaut_LFQ_indexingdf():
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
                    acquisition: string, "LFQ6 - Spectronaut", "LFQ5 - Spectronaut"
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
            
            # In case fractionated data was used this needs to be catched and aggregated
            try:
                df_index = df_index.unstack(["Map", "Fraction"])
            except ValueError:
                df_index = df_index.groupby(by=df_index.index.names).agg(np.nansum, axis=0)
                df_index = df_index.unstack(["Map", "Fraction"])
            
            df_index.replace(0, np.nan, inplace=True)
            df_index.rename(columns=self.fraction_dict, inplace=True)
            shape_dict["Original size"]=df_index.shape
            
            fraction_wCyt = list(df_index.columns.get_level_values("Fraction").unique())
            #Cyt is removed only if it is not an NMC split
            if "Cyt" in fraction_wCyt and len(fraction_wCyt) >= 4:
                df_index.drop("Cyt", axis=1, level="Fraction", inplace=True)
            try:
                if self.acquisition == "LFQ5 - Spectronaut":
                    df_index.drop("01K", axis=1, level="Fraction", inplace=True)
            except:
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
                df_organellarMarkerSet: df, columns: "Gene names", "Compartment", no index 
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
            #default setting: RatioHLcount = 2 ; RatioVariability = 30
            
            df_countvarfiltered_stacked = df_stack.loc[[count>RatioHLcount or (count==RatioHLcount and var<RatioVariability) 
                                            for var, count in zip(df_stack["Ratio H/L variability [%]"], df_stack["Ratio H/L count"])]]
            
            shape_dict["Shape after Ratio H/L count (>=3)/var (count==2, var<30) filtering"] = df_countvarfiltered_stacked.unstack(["Fraction", "Map"]).shape

            # "Ratio H/L":normalization to SILAC loading, each individual experiment (FractionXMap) will be divided by its median
            # np.median([...]): only entries, that are not NANs are considered
            df_normsilac_stacked = df_countvarfiltered_stacked["Ratio H/L"]\
                .unstack(["Fraction", "Map"])\
                .apply(lambda x: x/np.nanmedian(x), axis=0)\
                .stack(["Map", "Fraction"])

            df_stringency_mapfracstacked = df_countvarfiltered_stacked[["Ratio H/L count", "Ratio H/L variability [%]"]].join(
                pd.DataFrame(df_normsilac_stacked, columns=["Ratio H/L"]))

            # dataframe is grouped (Map, id), that allows the filtering for complete profiles
            df_stringency_mapfracstacked = df_stringency_mapfracstacked.groupby(["Map", "id"]).filter(lambda x: len(x)>=len(self.fractions))
            
            shape_dict["Shape after filtering for complete profiles"]=df_stringency_mapfracstacked.unstack(["Fraction", "Map"]).shape
            
            # Ratio H/L is converted into Ratio L/H
            df_stringency_mapfracstacked["Ratio H/L"] = df_stringency_mapfracstacked["Ratio H/L"].transform(lambda x: 1/x)
            
            #Annotation with marker genes
            df_organellarMarkerSet = self.df_organellarMarkerSet
            
            df_stringency_mapfracstacked.reset_index(inplace=True)
            df_stringency_mapfracstacked = df_stringency_mapfracstacked.merge(df_organellarMarkerSet, how="left", on="Gene names")
            df_stringency_mapfracstacked.set_index([c for c in df_stringency_mapfracstacked.columns
                                                    if c not in ["Ratio H/L count","Ratio H/L variability [%]","Ratio H/L"]], inplace=True)
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
                    df_organellarMarkerSet: df, columns: "Gene names", "Compartment", no index
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

            df_index = df_index.stack("Map")

            # sorting the level 0, in order to have LFQ intensity -    MS/MS count instead of continuous alternation
            df_index.sort_index(axis=1, level=0, inplace=True)
            
            # "MS/MS count"-column: take the sum over the fractions; if the sum is larger than n[fraction]*2, it will be stored in the new dataframe
            minms = (len(self.fractions) * self.summed_MSMS_counts)
            
            if minms > 0:
                df_mscount_mapstacked = df_index.loc[df_index[("MS/MS count")].apply(np.sum, axis=1) >= minms]
    
                shape_dict["Shape after MS/MS value filtering"]=df_mscount_mapstacked.unstack("Map").shape
                df_stringency_mapfracstacked = df_mscount_mapstacked.copy()
            else:
                df_stringency_mapfracstacked = df_index.copy()

            # series no dataframe is generated; if there are at least i.e. 4 consecutive non-NANs, data will be retained
            df_stringency_mapfracstacked = df_stringency_mapfracstacked.loc[
                df_stringency_mapfracstacked[("LFQ intensity")]\
                .apply(lambda x: np.isfinite(x), axis=0)\
                .apply(lambda x: sum(x) >= self.consecutiveLFQi and any(x.rolling(window=self.consecutiveLFQi).sum() >= self.consecutiveLFQi), axis=1)]
            
            shape_dict["Shape after consecutive value filtering"]=df_stringency_mapfracstacked.unstack("Map").shape

            df_stringency_mapfracstacked = df_stringency_mapfracstacked.copy().stack("Fraction")
            
            #Annotation with marker genes
            df_organellarMarkerSet = self.df_organellarMarkerSet
            
            df_stringency_mapfracstacked.reset_index(inplace=True)
            df_stringency_mapfracstacked = df_stringency_mapfracstacked.merge(df_organellarMarkerSet, how="left", on="Gene names")
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

            """
            df_01norm_mapstacked = df_stringency_mapfracstacked["LFQ intensity"].unstack("Fraction")
            
            # 0:1 normalization of Ratio L/H
            df_01norm_unstacked = df_01norm_mapstacked.div(df_01norm_mapstacked.sum(axis=1), axis=0)
            
            df_rest = df_stringency_mapfracstacked.drop("LFQ intensity", axis=1)
            df_01_stacked = df_rest.join(pd.DataFrame(df_01norm_unstacked.stack(
                   "Fraction"),columns=["LFQ intensity"]))
            
            # rename columns: "LFQ intensity" into "normalized profile"
            df_01_stacked.columns = [col if col!="LFQ intensity" else "normalized profile" for col in
                                     df_01_stacked.columns]
            
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
            
            df_rest = df_stringency_mapfracstacked.drop("LFQ intensity", axis=1)
            df_log_stacked = df_rest.join(pd.DataFrame(df_lognorm_ratio_stacked, columns=["LFQ intensity"]))
            
            # "LFQ intensity" will be renamed to "log profile"
            df_log_stacked.columns = [col if col!="LFQ intensity" else "log profile" for col in df_log_stacked.columns]

            return df_log_stacked
        
        def split_ids_uniprot(el):
            """
            This finds the primary canoncial protein ID in the protein group. If no canonical ID is present it selects the first isoform ID.
            """
            p1 = el.split(";")[0]
            if "-" not in p1:
                return p1
            else:
                p = p1.split("-")[0]
                if p in el.split(";"):
                    return p
                else:
                    return p1
        
        if self.acquisition == "SILAC - MQ":
            
            # Index data
            df_index = indexingdf()
            map_names = df_index.columns.get_level_values("Map").unique()
            self.map_names = map_names
            
            # Run stringency filtering and normalization
            df_stringency_mapfracstacked = stringency_silac(df_index)
            self.df_stringencyFiltered = df_stringency_mapfracstacked
            self.df_01_stacked = normalization_01_silac(df_stringency_mapfracstacked)
            self.df_log_stacked = logarithmization_silac(df_stringency_mapfracstacked)
            
            # format and reduce 0-1 normalized data for comparison with other experiments
            df_01_comparison = self.df_01_stacked.copy()
            comp_ids = pd.Series([split_ids_uniprot(el) for el in df_01_comparison.index.get_level_values("Protein IDs")], name="Protein IDs")
            df_01_comparison.index = df_01_comparison.index.droplevel("Protein IDs")
            df_01_comparison.set_index(comp_ids, append=True, inplace=True)
            df_01_comparison.drop(["Ratio H/L count", "Ratio H/L variability [%]"], inplace=True, axis=1)
            df_01_comparison = df_01_comparison.unstack(["Map", "Fraction"])
            df_01_comparison.columns = ["?".join(el) for el in df_01_comparison.columns.values]
            df_01_comparison = df_01_comparison.copy().reset_index().drop(["C-Score", "Q-value", "Score", "Majority protein IDs", "Protein names", "id"], axis=1, errors="ignore")
            
            # poopulate analysis summary dictionary with (meta)data
            unique_proteins = [split_ids_uniprot(i) for i in set(self.df_01_stacked.index.get_level_values("Protein IDs"))]
            unique_proteins.sort()
            self.analysis_summary_dict["0/1 normalized data"] = df_01_comparison.to_json()
            self.analysis_summary_dict["Unique Proteins"] = unique_proteins
            self.analysis_summary_dict["changes in shape after filtering"] = shape_dict.copy() 
            analysis_parameters = {"acquisition" : self.acquisition, 
                                   "filename" : self.filename,
                                   "comment" : self.comment,
                                   "Ratio H/L count" : self.RatioHLcount,
                                   "Ratio variability" : self.RatioVariability
                                  }
            self.analysis_summary_dict["Analysis parameters"] = analysis_parameters.copy()
            
            # TODO this line needs to be removed.
            self.analysed_datasets_dict[self.expname] = self.analysis_summary_dict.copy()


        elif self.acquisition == "LFQ5 - MQ" or self.acquisition == "LFQ6 - MQ" or self.acquisition == "LFQ5 - Spectronaut" or self.acquisition == "LFQ6 - Spectronaut":
        
            #if not summed_MS_counts:
            #    summed_MS_counts = self.summed_MS_counts
            #if not consecutiveLFQi:
            #    consecutiveLFQi = self.consecutiveLFQi
        
            if self.acquisition == "LFQ5 - MQ" or self.acquisition == "LFQ6 - MQ":
                df_index = indexingdf()
            elif self.acquisition == "LFQ5 - Spectronaut" or self.acquisition == "LFQ6 - Spectronaut":
                df_index = spectronaut_LFQ_indexingdf()
            
            map_names = df_index.columns.get_level_values("Map").unique()
            self.map_names = map_names
            
            df_stringency_mapfracstacked = stringency_lfq(df_index)
            self.df_stringencyFiltered = df_stringency_mapfracstacked
            self.df_log_stacked = logarithmization_lfq(df_stringency_mapfracstacked)
            self.df_01_stacked = normalization_01_lfq(df_stringency_mapfracstacked)
            
            df_01_comparison = self.df_01_stacked.copy()
            comp_ids = pd.Series([split_ids_uniprot(el) for el in df_01_comparison.index.get_level_values("Protein IDs")], name="Protein IDs")
            df_01_comparison.index = df_01_comparison.index.droplevel("Protein IDs")
            df_01_comparison.set_index(comp_ids, append=True, inplace=True)
            df_01_comparison.drop("MS/MS count", inplace=True, axis=1, errors="ignore")
            df_01_comparison = df_01_comparison.unstack(["Map", "Fraction"])
            df_01_comparison.columns = ["?".join(el) for el in df_01_comparison.columns.values]
            df_01_comparison = df_01_comparison.copy().reset_index().drop(["C-Score", "Q-value", "Score", "Majority protein IDs", "Protein names", "id"], axis=1, errors="ignore")
            self.analysis_summary_dict["0/1 normalized data"] = df_01_comparison.to_json()#double_precision=4) #.reset_index()
            
            unique_proteins = [split_ids_uniprot(i) for i in set(self.df_01_stacked.index.get_level_values("Protein IDs"))]
            unique_proteins.sort()
            self.analysis_summary_dict["Unique Proteins"] = unique_proteins
            self.analysis_summary_dict["changes in shape after filtering"] = shape_dict.copy() 
            analysis_parameters = {"acquisition" : self.acquisition, 
                                   "filename" : self.filename,
                                   "comment" : self.comment,
                                   "consecutive data points" : self.consecutiveLFQi,
                                   "summed MS/MS counts" : self.summed_MSMS_counts
                                  }
            self.analysis_summary_dict["Analysis parameters"] = analysis_parameters.copy() 
            self.analysed_datasets_dict[self.expname] = self.analysis_summary_dict.copy()
            #return self.df_01_stacked

        else:
            return "I do not know this"    
    
    
    def plot_log_data(self):     
        """
        
        Args:
            self.df_log_stacked
        
        
        Returns:
            log_histogram: Histogram of log transformed data
        
        """  
      
        log_histogram = px.histogram(self.df_log_stacked.reset_index().sort_values(["Map"]), 
                                     x="log profile",
                                     facet_col="Fraction",
                                     facet_row="Map",
                                     template="simple_white",
                                     labels={"log profile": "log tranformed data ({})".format("LFQ intenisty" if self.acquisition != "SILAC - MQ" else "Ratio H/L")}
                                    )
        
        log_histogram.for_each_xaxis(lambda axis: axis.update(title={"text":""}))
        log_histogram.for_each_yaxis(lambda axis: axis.update(title={"text":""}))
        log_histogram.add_annotation(x=0.5, y=0, yshift=-50, xref="paper",showarrow=False, yref="paper",
                    text="log2(LFQ intensity)")
        log_histogram.add_annotation(x=0, y=0.5, textangle=270, xref="paper",showarrow=False, yref="paper", xshift=-50,
                    text="count")
        log_histogram.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        
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
    
        if self.acquisition == "SILAC - MQ":
            df_index = self.df_index["Ratio H/L"]
            df_01_stacked = self.df_01_stacked["normalized profile"]
        elif self.acquisition == "LFQ5 - MQ" or self.acquisition == "LFQ6 - MQ" or self.acquisition == "LFQ5 - Spectronaut" or self.acquisition == "LFQ6 - Spectronaut":
            df_index = self.df_index["LFQ intensity"]
            df_01_stacked = self.df_01_stacked["normalized profile"].replace(0, np.nan)
        
        #unfiltered
        npg_t = df_index.shape[0]
        df_index_MapStacked = df_index.stack("Map")
        npr_t = df_index_MapStacked.shape[0]/len(self.map_names)
        npr_t_dc = 1-df_index_MapStacked.isna().sum().sum()/np.prod(df_index_MapStacked.shape)
        
        #filtered
        npgf_t = df_01_stacked.unstack(["Map", "Fraction"]).shape[0]
        df_01_MapStacked = df_01_stacked.unstack("Fraction")
        nprf_t = df_01_MapStacked.shape[0]/len(self.map_names)
        nprf_t_dc = 1-df_01_MapStacked.isna().sum().sum()/np.prod(df_01_MapStacked.shape)
        
        #unfiltered intersection
        try:
            df_index_intersection = df_index_MapStacked.groupby(level="Sequence").filter(lambda x : len(x)==len(self.map_names))
        except:
            df_index_intersection = df_index_MapStacked.groupby(level="Protein IDs").filter(lambda x : len(x)==len(self.map_names))
        npr_i = df_index_intersection.shape[0]/len(self.map_names)
        npr_i_dc = 1-df_index_intersection.isna().sum().sum()/np.prod(df_index_intersection.shape)
        npg_i = df_index_intersection.unstack("Map").shape[0]
        
        #filtered intersection
        try:
            df_01_intersection = df_01_MapStacked.groupby(level = "Sequence").filter(lambda x : len(x)==len(self.map_names))
        except:
            df_01_intersection = df_01_MapStacked.groupby(level = "Protein IDs").filter(lambda x : len(x)==len(self.map_names))
        nprf_i = df_01_intersection.shape[0]/len(self.map_names)
        nprf_i_dc = 1-df_01_intersection.isna().sum().sum()/np.prod(df_01_intersection.shape)
        npgf_i = df_01_intersection.unstack("Map").shape[0]
        
        # summarize in dataframe and save to attribute
        df_quantity_pr_pg = pd.DataFrame(
            {
            "filtering": pd.Series(["before filtering", "before filtering", "after filtering", "after filtering"], dtype=np.dtype("O")),
            "type": pd.Series(["total", "intersection", "total", "intersection"], dtype=np.dtype("O")),
            "number of protein groups": pd.Series([npg_t, npg_i, npgf_t, npgf_i], dtype=np.dtype("float")),
            "number of profiles": pd.Series([npr_t, npr_i, nprf_t, nprf_i], dtype=np.dtype("float")),
            "data completeness of profiles": pd.Series([npr_t_dc, npr_i_dc, nprf_t_dc, nprf_i_dc], dtype=np.dtype("float"))})
        
        self.df_quantity_pr_pg = df_quantity_pr_pg.reset_index()
        self.analysis_summary_dict["quantity: profiles/protein groups"] = self.df_quantity_pr_pg.to_json() 
        
        #additional depth assessment per fraction
        dict_npgf = {}
        dict_npg = {}
        list_npg_dc = []
        list_npgf_dc = []
        
        for df_intersection in [df_index_intersection, df_01_intersection]:
            for fraction in self.fractions:
                df_intersection_frac = df_intersection[fraction]
                npgF_f_dc = 1-df_intersection_frac.isna().sum()/len(df_intersection_frac)
                npgF_f = df_intersection_frac.unstack("Map").isnull().sum(axis=1).value_counts()
                if df_intersection.shape==df_index_intersection.shape:
                    dict_npg[fraction] = npgF_f
                    list_npg_dc.append(npgF_f_dc)
                else:
                    dict_npgf[fraction] = npgF_f
                    list_npgf_dc.append(npgF_f_dc)
        
        df_npg = pd.DataFrame(dict_npg)
        df_npg.index.name =  "Protein Groups present in:"
        df_npg.rename_axis("Fraction", axis=1, inplace=True)
        df_npg = df_npg.stack("Fraction").reset_index()
        df_npg = df_npg.rename({0: "Protein Groups"}, axis=1)
        df_npg.sort_values(["Fraction", "Protein Groups present in:"], inplace=True)     
        
        df_npgf = pd.DataFrame(dict_npgf)
        df_npgf.index.name =  "Protein Groups present in:"
        df_npgf.rename_axis("Fraction", axis=1, inplace=True)
        df_npgf = df_npgf.stack("Fraction").reset_index()
        df_npgf = df_npgf.rename({0: "Protein Groups"}, axis=1)
        df_npgf.sort_values(["Fraction", "Protein Groups present in:"], inplace=True)
        
        max_df_npg = df_npg["Protein Groups present in:"].max()
        min_df_npg = df_npg["Protein Groups present in:"].min()
        rename_numOFnans =  {}
        for x, y in zip(range(max_df_npg,min_df_npg-1, -1), range(max_df_npg+1)):
            if y == 1:
                rename_numOFnans[x] = "{} Map".format(y)
            elif y == 0:
                rename_numOFnans[x] = "PG not identified".format(y)
            else:
                rename_numOFnans[x] = "{} Maps".format(y)
        for keys in rename_numOFnans.keys():
            df_npg.loc[df_npg["Protein Groups present in:"] ==keys, "Protein Groups present in:"] = rename_numOFnans[keys]
            df_npgf.loc[df_npgf["Protein Groups present in:"] ==keys, "Protein Groups present in:"] = rename_numOFnans[keys]
        
        # summarize in dataframe and save to attributes
        self.df_npg_dc = pd.DataFrame(
            {
            "Fraction" : pd.Series(self.fractions),
            "Data completeness before filtering": pd.Series(list_npg_dc),
            "Data completeness after filtering": pd.Series(list_npgf_dc),
            }).sort_values(["Fraction"])
        
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
        
        layout = go.Layout(barmode="overlay", 
                           xaxis_tickangle=90, 
                           autosize=False,
                           width=300,
                           height=500,
                           xaxis=go.layout.XAxis(linecolor="black",
                                                 linewidth=1,
                                                 #title="Map",
                                                 mirror=True),
                           yaxis=go.layout.YAxis(linecolor="black",
                                                 linewidth=1,
                                                 mirror=True),
                           template="simple_white")
        
        fig_npg = go.Figure()
        for t in df_quantity_pr_pg["type"].unique():
            plot_df = df_quantity_pr_pg[df_quantity_pr_pg["type"] == t]
            fig_npg.add_trace(go.Bar(
                x=plot_df["filtering"],
                y=plot_df["number of protein groups"],
                name=t))
        fig_npg.update_layout(layout, title="Number of Protein Groups", yaxis=go.layout.YAxis(title="Protein Groups"))
        
                
        fig_npr = go.Figure()
        for t in df_quantity_pr_pg["type"].unique():
            plot_df = df_quantity_pr_pg[df_quantity_pr_pg["type"] == t]
            fig_npr.add_trace(go.Bar(
                x=plot_df["filtering"],
                y=plot_df["number of profiles"],
                name=t))
        fig_npr.update_layout(layout, title="Number of Profiles")
        
        df_quantity_pr_pg = df_quantity_pr_pg.sort_values("filtering")
        fig_npr_dc = go.Figure()
        for t in df_quantity_pr_pg["filtering"].unique():
            plot_df = df_quantity_pr_pg[df_quantity_pr_pg["filtering"] == t]
            fig_npr_dc.add_trace(go.Bar(
                x=plot_df["type"],
                y=plot_df["data completeness of profiles"],
                name=t))
        fig_npr_dc.update_layout(layout, title="Coverage", yaxis=go.layout.YAxis(title="Data completness"))
        #fig_npr_dc.update_xaxes(tickangle=30)
        
        fig_npg_F = px.bar(self.df_npg,
                          x="Fraction", 
                          y="Protein Groups", 
                          color="Protein Groups present in:",
                          template="simple_white",
                          title = "Protein groups per fraction - before filtering",
                          width=500)
        
        fig_npgf_F = px.bar(self.df_npgf,
                  x="Fraction", 
                  y="Protein Groups", 
                  color="Protein Groups present in:",
                  template="simple_white",
                  title = "Protein groups per fraction - after filtering",
                  width=500)
        
        fig_npg_F_dc = go.Figure()
        for data_type in ["Data completeness after filtering", "Data completeness before filtering"]:
            fig_npg_F_dc.add_trace(go.Bar(
                    x=self.df_npg_dc["Fraction"],
                    y=self.df_npg_dc[data_type],
                    name=data_type))
        fig_npg_F_dc.update_layout(layout, barmode="overlay", title="Data completeness per fraction", yaxis=go.layout.YAxis(title=""), height=450, width=600)
        
        
        return fig_npg, fig_npr, fig_npr_dc, fig_npg_F, fig_npgf_F, fig_npg_F_dc
                                                  
                                                                                
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

        if self.acquisition == "SILAC - MQ":
            df_01orlog_fracunstacked = self.df_log_stacked["log profile"].unstack("Fraction").dropna()
            df_01orlog_MapFracUnstacked = self.df_log_stacked["log profile"].unstack(["Fraction", "Map"]).dropna()  
            
        elif self.acquisition.startswith("LFQ"):
            df_01orlog_fracunstacked = self.df_01_stacked["normalized profile"].unstack("Fraction").dropna()
            df_01orlog_MapFracUnstacked = self.df_01_stacked["normalized profile"].unstack(["Fraction", "Map"]).dropna()  
            
        pca = PCA(n_components=3)

        # df_pca: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3"
        df_pca = pd.DataFrame(pca.fit_transform(df_01orlog_fracunstacked))
        df_pca.columns = ["PC1", "PC2", "PC3"]
        df_pca.index = df_01orlog_fracunstacked.index
        self.df_pca = df_pca.sort_index(level=["Gene names", "Compartment"])
        
        # df_pca: PCA processed dataframe, containing the columns "PC1", "PC2", "PC3"
        df_pca_combined = pd.DataFrame(pca.fit_transform(df_01orlog_MapFracUnstacked))
        df_pca_combined.columns = ["PC1", "PC2", "PC3"]
        df_pca_combined.index = df_01orlog_MapFracUnstacked.index
        self.df_pca_combined = df_pca_combined.sort_index(level=["Gene names", "Compartment"])
        
        map_names = self.map_names
        df_pca_all_marker_cluster_maps = pd.DataFrame()
        df_pca_filtered = df_pca.unstack("Map").dropna()
        for clusters in markerproteins:
            for marker in markerproteins[clusters]:
                try:
                    plot_try_pca = df_pca_filtered.xs(marker, level="Gene names", drop_level=False)
                except KeyError:
                    continue
                df_pca_all_marker_cluster_maps = df_pca_all_marker_cluster_maps.append(
                    plot_try_pca)
        df_pca_all_marker_cluster_maps = df_pca_all_marker_cluster_maps.stack("Map")
        self.df_pca_all_marker_cluster_maps = df_pca_all_marker_cluster_maps.sort_index(level=["Gene names", "Compartment"])

        
    def plot_global_pca(self, map_of_interest="Map1", cluster_of_interest="Proteasome", x_PCA="PC1", y_PCA="PC3", collapse_maps=False):
        """"
        PCA plot will be generated

        Args:
            self:
                df_organellarMarkerSet: df, columns: "Gene names", "Compartment", no index
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

        compartments = self.df_organellarMarkerSet["Compartment"].unique()
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
                                    template="simple_white",
                                    opacity=0.9
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
                        plot_try_pca = df_pca_all_marker_cluster_maps.xs((marker, maps), level=["Gene names", "Map"],
                                                                         drop_level=False)
                    except KeyError:
                        continue
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
                                  title="PCA plot for <br>the protein cluster: {}".format(cluster_of_interest), 
                                  template="simple_white")
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
            df_cluster_unfiltered = self.get_marker_proteins_unfiltered(cluster)
            df_allclusters_01_unfiltered_mapfracunstacked = df_allclusters_01_unfiltered_mapfracunstacked.append(df_cluster_unfiltered)
            
            # filter for coverage and calculate distances
            df_cluster = df_cluster_unfiltered.dropna()
            if len(df_cluster) == 0:
                continue
            df_distances_aggregated, df_distances_individual = self.calc_cluster_distances(df_cluster)
            df_alldistances_individual_mapfracunstacked = df_alldistances_individual_mapfracunstacked.append(df_distances_individual)
            df_alldistances_aggregated_mapunstacked = df_alldistances_aggregated_mapunstacked.append(df_distances_aggregated)
        
        ## Get compatibility with plotting functions, by mimicking assignment of old functions:
        # old output of distance_calculation
        self.df_distance_noindex = df_alldistances_aggregated_mapunstacked.stack("Map").reset_index().rename({0: "distance"}, axis=1) 
        self.analysis_summary_dict["Manhattan distances"] = self.df_distance_noindex.to_json()
        # old output of multiple_iterations
        # self.df_allclusters_clusterdist_fracunstacked_unfiltered --> this won't exist anymore, replaced by:
        self.df_allclusters_01_unfiltered_mapfracunstacked = df_allclusters_01_unfiltered_mapfracunstacked
        # kept for testing of quantification table:
        self.df_allclusters_clusterdist_fracunstacked_unfiltered = df_allclusters_01_unfiltered_mapfracunstacked.stack("Map")
        # same as before, but now already abs
        self.df_allclusters_clusterdist_fracunstacked = df_alldistances_individual_mapfracunstacked.stack("Map")
        df_dist_to_median = self.df_allclusters_clusterdist_fracunstacked.stack("Fraction")
        df_dist_to_median.name = "distance"
        df_dist_to_median = df_dist_to_median.reindex(index=natsort.natsorted(df_dist_to_median.index))
        self.analysis_summary_dict["Distances to the median profile"] = df_dist_to_median.reset_index().to_json()
        self.genenames_sortedout_list = [el for el in df_allclusters_01_unfiltered_mapfracunstacked.index.get_level_values("Gene names")
                                         if el not in df_alldistances_individual_mapfracunstacked.index.get_level_values("Gene names")]
                                         
        return df_alldistances_individual_mapfracunstacked, df_alldistances_aggregated_mapunstacked, df_allclusters_01_unfiltered_mapfracunstacked
    
    
    def get_marker_proteins_unfiltered(self, cluster):
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
        df_in = self.df_01_stacked["normalized profile"].unstack("Fraction")
        markers = self.markerproteins[cluster]
        
        # retrieve marker proteins
        df_cluster_unfiltered = pd.DataFrame()
        for marker in markers:
            try:
                df_p = df_in.xs(marker, level="Gene names", axis=0, drop_level=False)
            except:
                continue
            df_cluster_unfiltered = df_cluster_unfiltered.append(df_p)
        if len(df_cluster_unfiltered) == 0:
            return df_cluster_unfiltered
        
        # Unstack maps and add Cluster to index
        df_cluster_unfiltered = df_cluster_unfiltered.unstack("Map")
        df_cluster_unfiltered.set_index(pd.Index(np.repeat(cluster, len(df_cluster_unfiltered)), name="Cluster"), append=True, inplace=True)
        
        return df_cluster_unfiltered
    
    
    def calc_cluster_distances(self, df_cluster, complex_profile=np.median, distance_measure="manhattan"):
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
        
        ref_profile = pd.DataFrame(df_cluster.apply(complex_profile, axis=0, result_type="expand")).T
        df_distances_individual = df_cluster.apply(lambda x: np.abs(x-ref_profile.iloc[0,:]), axis=1)
            
        # loop over maps
        maps = set(df_cluster.columns.get_level_values("Map"))
        for m in maps:
            if distance_measure == "manhattan":
                d_m = pw.manhattan_distances(df_cluster.xs(m, level="Map", axis=1), ref_profile.xs(m, level="Map", axis=1))
            else:
                raise ValueError(distance_measure)
            d_m = pd.DataFrame(d_m, columns=[m], index=df_cluster.index)
            df_distances_aggregated = pd.concat([df_distances_aggregated, d_m], axis=1)
        
        df_distances_aggregated.columns.set_names(names="Map", inplace=True)
        return df_distances_aggregated, df_distances_individual
    
    
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
            df_setofproteins = self.df_allclusters_01_unfiltered_mapfracunstacked.xs(cluster_of_interest, level="Cluster", axis=0)
            df_setofproteins_median = df_setofproteins.dropna().xs(map_of_interest, level="Map", axis=1).median(axis=0)
    
            # fractions get sorted
            df_setofproteins = df_setofproteins.xs(map_of_interest, level="Map", axis=1).stack("Fraction")
            df_setofproteins = df_setofproteins.reindex(index=natsort.natsorted(df_setofproteins.index))
            df_setofproteins.name = "normalized profile"
    
            # make it available for plotting
            df_setofproteins = df_setofproteins.reset_index()
            abundance_profiles_figure = px.line(df_setofproteins, 
                                                x="Fraction", 
                                                y="normalized profile",
                                                color="Gene names",
                                                line_group="Sequence" if "Sequence" in df_setofproteins.columns else "Gene names",
                                                template="simple_white",
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
            # dash lines for proteins that have insufficient coverage across maps
            abundance_profiles_and_median_figure.for_each_trace(lambda x: x.update(line={"dash":"dash"}),
                                                                selector=lambda x: x.name in self.genenames_sortedout_list)
    
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
        
        df_quantification_overview = self.df_allclusters_clusterdist_fracunstacked_unfiltered.xs(cluster_of_interest, level="Cluster", axis=0)\
                                                                                             [self.fractions[0]].unstack("Map")
        if "Sequence" in df_quantification_overview.index.names:
            df_quantification_overview = df_quantification_overview.droplevel([i for i in df_quantification_overview.index.names if not i in ["Sequence","Gene names"]])
        else:
            df_quantification_overview = df_quantification_overview.droplevel([i for i in df_quantification_overview.index.names if not i=="Gene names"])
        df_quantification_overview = df_quantification_overview.notnull().replace({True: "x", False: "-"})
        
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
        df_distance_map_cluster_gene_in_index = df_distance_noindex.set_index(["Gene names", "Map", "Cluster"])
        if "Sequence" in df_distance_map_cluster_gene_in_index.columns:
            df_distance_map_cluster_gene_in_index.set_index("Sequence", append=True, inplace=True)

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
            self.proteins_quantified_across_all_maps = df_cluster_xmaps_distance_with_index.unstack("Map").shape[0]
        
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
                template="simple_white"
            )
    
            return distance_boxplot_figure
    
        except:
            self.cache_cluster_quantified = False
            
    
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
                plot_try = self.df_allclusters_clusterdist_fracunstacked.xs((cluster_of_interest, maps), level=["Cluster", "Map"], drop_level=False)
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
                                                       facet_col="Fraction",
                                                       facet_col_wrap=2,
                                                       boxmode="overlay", height=900, width=700, points="all",
                                                       hover_name="Gene names",
                                                       template="simple_white",
                                                       title="Distribution of the distance to the median for <br>the protein cluster: {}".format(cluster_of_interest))
    
            return distance_to_median_boxplot_figure
            
        except:
            return "This protein cluster was not quantified" 

        
    def dynamic_range(self):
        """
        Dynamic range of each individual protein clusters (of the median profile) across all maps is calculated"

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
                    try:
                        df_marker_allMaps = df_01_stacked.xs(marker, level="Gene names", drop_level=False)
                    except KeyError:
                        continue
                    df_setofproteins_allMaps = df_setofproteins_allMaps.append(df_marker_allMaps)
                df_setofproteins_allMaps_median = df_setofproteins_allMaps["normalized profile"].unstack("Fraction").median()
                
                df_dynamicRange = df_dynamicRange.append(pd.DataFrame(np.array([[max(df_setofproteins_allMaps_median), 
                                                                                 min(df_setofproteins_allMaps_median), 
                                                                                 max(df_setofproteins_allMaps_median)-min(df_setofproteins_allMaps_median),
                                                                                 clusters]]), 
                                                                      columns=["Max", "Min", "Dynamic Range", "Cluster"]),
                                                        ignore_index=True)
            except:
                continue
        
        self.analysis_summary_dict["Dynamic Range"] = df_dynamicRange.to_json()
        
    
    def plot_dynamic_range(self):
        """
        Dynamic range of each individual protein clusters (of the median profile) across all maps is displayed"

        Args:
            self:
                markerproteins: dictionary, key: cluster name, value: gene names (e.g. {"Proteasome" : ["PSMA1", "PSMA2",...], ...}
                df_01_stacked: "MAP" and "Fraction" are stacked; the data in the column "normalized profile" is used for plotting. Additionally the columns 
                               "MS/MS count" and "Ratio H/L count | Ratio H/L variability [%] | Ratio H/L" are found in LFQ and SILAC data respectively

        Returns:
            fig_dynamicRange: Bar plot, displaying the dynamic range for each protein cluster
            self.df_dynamicRange: df, no index, columns: "Max", "Min", "Dynamic Range", "Cluster"
        """
        
        fig_dynamicRange = px.bar(pd.read_json(self.analysis_summary_dict["Dynamic Range"]), 
                                  x="Cluster", 
                                  y="Dynamic Range", 
                                  base="Min", 
                                  template="simple_white",
                                  width=1000, 
                                  height=500).update_xaxes(categoryorder="total ascending")
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
#self.analysis_summary_dict.clear()
        
        return df_overview


    def reframe_df_01ORlog_for_Perseus(self, df_01ORlog):
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
        #df_01ORlog_svm = df_01ORlog_svm.dropna(axis=0, subset=df_01ORlog_svm.loc[[], ["normalized profile"]].columns)
        df_01ORlog_svm.columns = ["_".join(col) for col in df_01ORlog_svm.columns.values]
        df_01ORlog_svm.rename(index={"undefined" : np.nan}, level="Compartment", inplace=True)
        
        return df_01ORlog_svm

        

    

class SpatialDataSetComparison:
    
        
    analysed_datasets_dict = SpatialDataSet.analysed_datasets_dict
    css_color = SpatialDataSet.css_color
    markerproteins_set = SpatialDataSet.markerproteins_set
    cache_stored_SVM = True


    def __init__(self, clusters_for_ranking=["Proteasome", "Lysosome"], ref_exp="Exp2", **kwargs):
        
        self.clusters_for_ranking = clusters_for_ranking
        self.ref_exp = ref_exp
        self.json_dict = {}
        #self.fractions, self.map_names = [], []  #self.df_01_stacked, self.df_log_stacked = pd.DataFrame(), pd.DataFrame()
        #collapse_maps,collapse_cluster,  cluster_of_interest_comparison, multi_choice, multi_choice_venn, x_PCA_comp, y_PCA_comp
        
        
        if "organism" not in kwargs.keys():
            self.markerproteins = self.markerproteins_set["Human - Swissprot"]
        else:
            assert kwargs["organism"] in self.markerproteins_set.keys()
            self.markerproteins = self.markerproteins_set[kwargs["organism"]]
            del kwargs["organism"]
        
        #self.unique_proteins_total = unique_proteins_total
        
        self.exp_names, self.exp_map_names = [], []
        
        self.df_01_filtered_combined, self.df_distance_comp = pd.DataFrame(), pd.DataFrame()
        self.df_quantity_pr_pg_combined, self.df_dynamicRange_combined = pd.DataFrame(), pd.DataFrame()
        

    def read_jsonFile(self): #, content=None
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
                df_dynamicRange_combined: df, no index, column names: "Max", "Min", "Dynamic Range", "Cluster", "Experiment"
                unique_proteins_total: dict, key: Experiment name, value: unique protein (groups)
                exp_map_names: list of unique Exp_Map - fusions e.g. LFQ_Map1
                exp_names: list of unique Experiment names - e.g. LFQ
        """
        
        json_dict = self.json_dict
        
        #add experiments that are not stored in AnalysedDAtasets.json for comparison
        #try:
        #if len(SpatialDataSet.analysed_datasets_dict.keys())>=1:
        #    json_dict.update(SpatialDataSet.analysed_datasets_dict)
        ##except:
        #else:
        #    pass
    
        self.analysis_parameters_total = {}
        unique_proteins_total = {}
        
        for exp_name in json_dict.keys():
            for data_type in json_dict[exp_name].keys():                   
                if data_type == "0/1 normalized data" and exp_name == list(json_dict.keys())[0]:
                    df_01_combined = pd.read_json(json_dict[exp_name][data_type])
                    df_01_combined = df_01_combined.set_index(["Gene names", "Protein IDs", "Compartment"]).copy()
                    if "Sequence" in df_01_combined.columns:
                        df_01_combined.set_index(["Sequence"], inplace=True, append=True)
                    df_01_combined.drop([col for col in df_01_combined.columns if not col.startswith("normalized profile")])
                    df_01_combined.columns = pd.MultiIndex.from_tuples([el.split("?") for el in df_01_combined.columns], names=["Set", "Map", "Fraction"])
                    df_01_combined.rename(columns = {"normalized profile":exp_name}, inplace=True)
        
                elif data_type == "0/1 normalized data" and exp_name != list(json_dict.keys())[0]:
                    df_01_toadd = pd.read_json(json_dict[exp_name][data_type])
                    df_01_toadd = df_01_toadd.set_index(["Gene names", "Protein IDs", "Compartment"]).copy()
                    if "Sequence" in df_01_toadd.columns:
                        df_01_toadd.set_index(["Sequence"], inplace=True, append=True)
                    df_01_toadd.drop([col for col in df_01_toadd.columns if not col.startswith("normalized profile")])
                    df_01_toadd.columns = pd.MultiIndex.from_tuples([el.split("?") for el in df_01_toadd.columns], names=["Set", "Map", "Fraction"])
                    df_01_toadd.rename(columns = {"normalized profile":exp_name}, inplace=True)
                    df_01_combined = pd.concat([df_01_combined, df_01_toadd], axis=1)#, join="inner")  
                        
                elif data_type == "quantity: profiles/protein groups" and exp_name == list(json_dict.keys())[0]:
                    df_quantity_pr_pg_combined = pd.read_json(json_dict[exp_name][data_type])
                    df_quantity_pr_pg_combined["Experiment"] = exp_name
                    
                elif data_type == "quantity: profiles/protein groups" and exp_name != list(json_dict.keys())[0]:
                    df_quantity_pr_pg_toadd = pd.read_json(json_dict[exp_name][data_type])
                    df_quantity_pr_pg_toadd["Experiment"] = exp_name
                    df_quantity_pr_pg_combined = pd.concat([df_quantity_pr_pg_combined, df_quantity_pr_pg_toadd])  
                    
                elif data_type == "Manhattan distances" and exp_name == list(json_dict.keys())[0]:
                    df_distances_combined = pd.read_json(json_dict[exp_name][data_type])
                    df_distances_combined = df_distances_combined.set_index(["Map", "Gene names", "Cluster", "Protein IDs", "Compartment"]).copy()
                    if "Sequence" in df_distances_combined.columns:
                        df_distances_combined.set_index(["Sequence"], inplace=True, append=True)
                    df_distances_combined = df_distances_combined[["distance"]].unstack(["Map"])
                    df_distances_combined.rename(columns = {"distance":exp_name}, inplace=True)
        
                elif data_type == "Manhattan distances" and exp_name != list(json_dict.keys())[0]:
                    df_distances_toadd = pd.read_json(json_dict[exp_name][data_type])
                    df_distances_toadd = df_distances_toadd.set_index(["Map", "Gene names", "Cluster", "Protein IDs", "Compartment"]).copy()
                    if "Sequence" in df_distances_toadd.columns:
                        df_distances_toadd.set_index(["Sequence"], inplace=True, append=True)
                    df_distances_toadd = df_distances_toadd[["distance"]].unstack(["Map"])
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
                
                elif data_type == "Analysis parameters":
                    self.analysis_parameters_total[exp_name] = json_dict[exp_name][data_type]
                    
                #try:
                #    for paramters in json_dict[exp_name][data_type].keys():
                #        if paramters=="acquisition":
                #            acquisition_loaded.append(json_dict[exp_name][data_type][paramters])
                #        #elif parameters=="Non valid profiles":
                #except:
                #    continue
                #
                
        #filter for consistently quantified proteins (they have to be in all fractions and all maps)
        #df_01_filtered_combined = df_01_mean_combined.dropna()    
        df_01_combined.columns.names = ["Experiment", "Map", "Fraction"]
        #reframe it to make it ready for PCA
        df_01_filtered_combined = df_01_combined.stack(["Experiment", "Map"]).dropna(axis=0)
        #df_01_filtered_combined = df_01_combined.stack(["Experiment"]).dropna(axis=1)
        df_01_filtered_combined = df_01_filtered_combined.div(df_01_filtered_combined.sum(axis=1), axis=0)
        #df_01_filtered_combined = df_01_combined.copy()
        #df_01_filtered_combined.columns.names = ["Experiment", "Fraction", "Map"]
        ## Replace protein IDs by the unifying protein ID across experiments
        #comparison_IDs = pd.Series([split_ids_uniprot(el) for el in df_01_filtered_combined.index.get_level_values("Protein IDs")],
        #                           name="Protein IDs")
        #df_01_filtered_combined.index = df_01_filtered_combined.index.droplevel("Protein IDs")
        #df_01_filtered_combined.set_index(comparison_IDs, append=True, inplace=True)
        ##reframe it to make it ready for PCA | dropna: to make sure, that you do consider only fractions that are in all experiments
        #df_01_filtered_combined = df_01_filtered_combined.stack(["Experiment", "Map"]).swaplevel(0,1, axis=0).dropna(axis=1)
        index_ExpMap = df_01_filtered_combined.index.get_level_values("Experiment")+"_"+df_01_filtered_combined.index.get_level_values("Map")
        index_ExpMap.name = "Exp_Map"
        df_01_filtered_combined.set_index(index_ExpMap, append=True, inplace=True)
        
        
        df_distances_combined.columns.names = ["Experiment", "Map"]
        series = df_distances_combined.stack(["Experiment", "Map"])
        series.name = "distance"
        
        df_distance_comp = series.to_frame()
        #fuse Experiment and Map into one column = "Exp_Map"
        index_dist_ExpMap = df_distance_comp.index.get_level_values("Experiment")+"_"+df_distance_comp.index.get_level_values("Map")
        index_dist_ExpMap.name = "Exp_Map"
        df_distance_comp.set_index(index_dist_ExpMap, append=True, inplace=True)
        #new
        #self.df_distance_comp2 = df_distance_comp.copy()
        df_distance_comp.reset_index(level=['Protein IDs'], inplace=True)
        df_distance_comp["Protein IDs"] = df_distance_comp["Protein IDs"].str.split(";", expand=True)[0]
        df_distance_comp = df_distance_comp.set_index("Protein IDs", append=True).unstack(["Experiment", "Exp_Map", "Map"]).dropna().stack(["Experiment", "Exp_Map", "Map"]).reset_index()
        #df_distance_comp.reset_index(inplace=True)
        
        self.unique_proteins_total = unique_proteins_total
        self.exp_names = list(df_01_filtered_combined.index.get_level_values("Experiment").unique())
        self.exp_map_names = list(index_dist_ExpMap.unique())
        
        self.df_01_filtered_combined = df_01_filtered_combined 
        #self.df_01_mean_filtered_combined = df_01_mean_filtered_combined
        
        self.df_quantity_pr_pg_combined = df_quantity_pr_pg_combined
        self.df_dynamicRange_combined = df_dynamicRange_combined
        
        self.df_distance_comp = df_distance_comp

    
    def perform_pca_comparison(self):
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
                df_global_pca: PCA processed dataframe
                    index: "Gene names", "Protein IDs", "Compartment", "Experiment", 
                    columns: "PC1", "PC2", "PC3"
                    contains all protein IDs, that are consistent throughout all experiments
        """

        markerproteins = self.markerproteins.copy()
        
        #df_01_filtered_combined = self.df_01_filtered_combined
        #df_01_filtered_combined = self.df_01_filtered_combined 
        
        
        df_mean = pd.DataFrame()
        for exp in self.exp_names:
            df_exp = self.df_01_filtered_combined.stack("Fraction").unstack(["Experiment", "Map","Exp_Map"])[exp].mean(axis=1).to_frame(name=exp)
            df_mean = pd.concat([df_mean, df_exp], axis=1)
        df_mean = df_mean.rename_axis("Experiment", axis="columns").stack("Experiment").unstack("Fraction")
        
        markerproteins = i_class_comp.markerproteins.copy()
        pca = PCA(n_components=3)
        
        df_pca = pd.DataFrame(pca.fit_transform(df_mean))
        df_pca.columns = ["PC1", "PC2", "PC3"]
        df_pca.index = df_mean.index
        
        
        try:
            markerproteins["PSMA subunits"] = [item for sublist in [re.findall("PSMA.*",p) for p in markerproteins["Proteasome"]] for item in sublist]
            markerproteins["PSMB subunits"] = [item for sublist in [re.findall("PSMB.*",p) for p in markerproteins["Proteasome"]] for item in sublist]
            del markerproteins["Proteasome"]
        except:
            pass 
        
        ###only one df, make annotation at that time
        df_cluster = pd.DataFrame([(k, i) for k, l in markerproteins.items() for i in l], columns=["Cluster", "Gene names"])
        df_global_pca = df_pca.reset_index().merge(df_cluster, how="left", on="Gene names")
        df_global_pca.Cluster.replace(np.NaN, "Undefined", inplace=True)
        
        self.markerproteins_splitProteasome = markerproteins
        self.df_pca = df_pca 
        self.df_global_pca = df_global_pca
            
            
    def plot_pca_comparison(self, cluster_of_interest_comparison="Proteasome", multi_choice=["Exp1", "Exp2"]):
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
            df_setofproteins_PCA = pd.DataFrame()
            for map_or_exp in multi_choice:
                    for marker in markerproteins[cluster_of_interest_comparison]:
                        try:
                            plot_try_pca = df_pca.xs((marker, map_or_exp), level=["Gene names", "Experiment"], drop_level=False)
                        except KeyError:
                            continue
                        df_setofproteins_PCA = df_setofproteins_PCA.append(plot_try_pca)
            df_setofproteins_PCA.reset_index(inplace=True)
            
            df_setofproteins_PCA = df_setofproteins_PCA.assign(Experiment_lexicographic_sort=pd.Categorical(df_setofproteins_PCA["Experiment"], categories=self.sorting_list,
                                                                                                              ordered=True))
            df_setofproteins_PCA.sort_values("Experiment_lexicographic_sort", inplace=True)
                    
            pca_figure = px.scatter_3d(df_setofproteins_PCA, 
                                       x="PC1", 
                                       y="PC2", 
                                       z="PC3", 
                                       color="Experiment",
                                       template="simple_white",
                                       hover_data=["Gene names"]
                                      )
                    
            pca_figure.update_layout(autosize=False, 
                                     width=1400, 
                                     height=500,
                                     title="PCA plot for <br>the protein cluster: {}".format(cluster_of_interest_comparison),
                                     template="simple_white"
                                    )
    
            return pca_figure
        except:
            return "This protein cluster was not identified in all experiments"
            
                
    def plot_global_pca_comparison(self, cluster_of_interest_comparison="Proteasome", x_PCA="PC1", y_PCA="PC3", 
        markerset_or_cluster=False, multi_choice=["Exp1", "Exp2"]):
        """"
        PCA plot will be generated
    
        Args:
            self:
                df_organellarMarkerSet: df, columns: "Gene names", "Compartment", no index
                multi_choice: list of experiment names
                css_color: list of colors
                df_global_pca: PCA processed dataframe
                    index: "Gene names", "Protein IDs", "Compartment", "Experiment", 
                    columns: "PC1", "PC2", "PC3"
                    contains all protein IDs, that are consistent throughout all experiments    
    
        Returns:
            pca_figure: global PCA plot, clusters based on the markerset based (df_organellarMarkerSet) are color coded. 
        """
        
        
        df_global_pca_exp = self.df_global_pca.loc[self.df_global_pca["Experiment"].isin(multi_choice)]
        df_global_pca_exp.reset_index(inplace=True)

        compartments = list(SpatialDataSet.df_organellarMarkerSet["Compartment"].unique())
        compartment_color = dict(zip(compartments, self.css_color))
        compartment_color["Selection"] = "black"
        compartment_color["undefined"] = "lightgrey"
        compartments.insert(0, "undefined")
        compartments.insert(len(compartments), "Selection")
            
        cluster = self.markerproteins_splitProteasome.keys()
        cluster_color = dict(zip(cluster, self.css_color))
        cluster_color["Undefined"] = "lightgrey"
                
        
        if markerset_or_cluster == True:
            df_global_pca = df_global_pca_exp[df_global_pca_exp.Cluster!="Undefined"].sort_values(by="Cluster")
            df_global_pca = df_global_pca_exp[df_global_pca_exp.Cluster=="Undefined"].append(df_global_pca)
        else:
            for i in self.markerproteins[cluster_of_interest_comparison]:
                df_global_pca_exp.loc[df_global_pca_exp["Gene names"] == i, "Compartment"] = "Selection"
            df_global_pca = df_global_pca_exp.assign(Compartment_lexicographic_sort = pd.Categorical(df_global_pca_exp["Compartment"], 
                                                                                                     categories=[x for x in compartments], 
                                                                                                     ordered=True))
            df_global_pca.sort_values(["Compartment_lexicographic_sort", "Experiment"], inplace=True)
            
        fig_global_pca = px.scatter(data_frame=df_global_pca,
                                    x=x_PCA,
                                    y=y_PCA,
                                    color="Compartment" if markerset_or_cluster == False else "Cluster",
                                    color_discrete_map=compartment_color if markerset_or_cluster == False else cluster_color,
                                    title="Protein subcellular localization by PCA",
                                    hover_data=["Protein IDs", "Gene names", "Compartment"],
                                    facet_col="Experiment",
                                    facet_col_wrap=2,
                                    opacity=0.9, 
                                    template="simple_white"
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
                                     width=1800 if markerset_or_cluster == False else 1600, 
                                     height=400*(int(len(multi_choice) / 2) + (len(multi_choice) % 2 > 0)),
                                     template="simple_white"
                                    )
        
        return fig_global_pca
    
    
    def get_marker_proteins(self, experiments, cluster):
        df_in = self.df_01_filtered_combined.copy()
        markers = self.markerproteins[cluster]
        
        # retrieve marker proteins
        df_cluster = pd.DataFrame()
        for marker in markers:
            try:
                df_p = df_in.xs(marker, level="Gene names", axis=0, drop_level=False)
            except:
                continue
            df_cluster = df_cluster.append(df_p)
        if len(df_cluster) == 0:
            return df_cluster
        
        # filter for all selected experiments
        df_cluster = df_cluster.droplevel("Exp_Map", axis=0)
        df_cluster = df_cluster.unstack(["Experiment", "Map"])
        drop_experiments = [el for el in df_cluster.columns.get_level_values("Experiment") if el not in experiments]
        if len(drop_experiments) > 0:
            df_cluster.drop([el for el in df_cluster.columns.get_level_values("Experiment") if el not in experiments],
                            level="Experiment", axis=1, inplace=True)
        df_cluster.dropna(inplace=True)
        if len(df_cluster) == 0:
            return df_cluster
        df_cluster.set_index(pd.Index(np.repeat(cluster, len(df_cluster)), name="Cluster"), append=True, inplace=True)
        
        return df_cluster
    
    
    def calc_cluster_distances(self, df_cluster, complex_profile=np.median, distance_measure="manhattan"):
        df_distances = pd.DataFrame()
        
        # loop over experiments
        experiments = set(df_cluster.columns.get_level_values("Experiment"))
        for exp in experiments:
            df_exp = df_cluster.xs(exp, level="Experiment", axis=1)
            ref_profile = pd.DataFrame(df_exp.apply(complex_profile, axis=0, result_type="expand")).T
            
            # loop over maps
            maps = set(df_exp.columns.get_level_values("Map"))
            for m in maps:
                if distance_measure == "manhattan":
                    d_m = pw.manhattan_distances(df_exp.xs(m, level="Map", axis=1), ref_profile.xs(m, level="Map", axis=1))
                else:
                    raise ValueError(distance_measure)
                d_m = pd.DataFrame(d_m, columns=[(exp, m)], index=df_exp.index)
                df_distances = pd.concat([df_distances, d_m], axis=1)
        
        df_distances.columns = pd.MultiIndex.from_tuples(df_distances.columns, names=["Experiment", "Map"])
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
            df_cluster = self.get_marker_proteins(experiments, cluster)
            if len(df_cluster) == 0:
                continue
            dists_cluster = self.calc_cluster_distances(df_cluster)
            df_distances = df_distances.append(dists_cluster)
        df_distances = df_distances.stack(["Experiment", "Map"]).reset_index()\
            .sort_values(["Experiment","Gene names"]).rename({0: "distance"}, axis=1)
        df_distances.insert(0, "Exp_Map", ["_".join([e,m]) for e,m in zip(df_distances["Experiment"], df_distances["Map"])])
        
        self.df_distance_comp = df_distances
        
        return df_distances
    
    
    def distance_boxplot_comparison(self, cluster_of_interest_comparison="Proteasome", collapse_maps=False, multi_choice=["Exp1", "Exp2"]):
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
        
        #an error massage, if no Experiments are selected, will be displayed already, that is why: return ""
        if len(multi_choice)>=1:
            pass
        else:
            return ("")
        
        df_distance_comp = self.df_distance_comp.copy()
        #set categroical column, allowing lexicographic sorting
        df_distance_comp["Experiment_lexicographic_sort"] = pd.Categorical(df_distance_comp["Experiment"], categories=self.sorting_list, ordered=True)
        df_distance_comp.sort_values(["Experiment_lexicographic_sort", "Map"], inplace=True)
        
        if collapse_maps == False:
            #get only values form experiment of interest
            df_distance_selectedExp = df_distance_comp.loc[df_distance_comp["Experiment"].isin(multi_choice)]
            #get only values form cluster of interest
            df_distance_selectedExp = df_distance_selectedExp.loc[df_distance_selectedExp["Cluster"]==cluster_of_interest_comparison]
            
            if df_distance_selectedExp.shape[0] == 0:
                self.cache_cluster_quantified = False
                
            else:
                individual_distance_boxplot_figure=go.Figure()
                for i, exp in enumerate(df_distance_selectedExp["Experiment"].unique()):
                    df_plot=df_distance_selectedExp[df_distance_selectedExp["Experiment"]==exp]
                    individual_distance_boxplot_figure.add_trace(go.Box(
                        x=[df_plot["Experiment"], df_plot["Map"]],
                        y=df_plot["distance"],
                        #line=dict(color=pio.templates["simple_white"].layout["colorway"][i]),
                        boxpoints="all",
                        whiskerwidth=0.2,
                        marker_size=2,
                        name=exp, 
                        hovertext=df_plot["Gene names"]
                    ))
                    
                individual_distance_boxplot_figure.update_layout(boxmode="group", 
                                                                 xaxis_tickangle=90, 
                                                                 title="Manhattan distance distribution for <br>the protein cluster: {}".format(cluster_of_interest_comparison),
                                                                 autosize=False,
                                                                 width=350*len(multi_choice),
                                                                 height=500,
                                                                 xaxis=go.layout.XAxis(linecolor="black",
                                                                                       linewidth=1,
                                                                                       title="Experiment",
                                                                                       mirror=True),
                                                                 yaxis=go.layout.YAxis(linecolor="black",
                                                                                       linewidth=1,
                                                                                       title="Distance",
                                                                                       mirror=True),
                                                                 template="simple_white")
                
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
                plot_try = df_distance_selectedExp_global.xs((cluster_of_interest_comparison, map_or_exp), level=["Cluster", 
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
                                             template="simple_white",
                                             title="Global Manhattan distance distribution for the protein cluster: {}".format(cluster_of_interest_comparison)
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
                                                                        mirror=True),
                                                  template="simple_white"
                                                 )
            
            return distance_boxplot_figure
        
        
  
    def distance_ranking_barplot_comparison(self, collapse_cluster=False, multi_choice=["Exp1", "Exp2"], clusters_for_ranking=None, ranking_boxPlot="Box plot"):#, toggle_sumORmedian=False):
    #ref_exp="Exp1", 
        if clusters_for_ranking is None:
                clusters_for_ranking = self.clusters_for_ranking
            
            #an error massage, if no Experiments are selected, will be displayed already, that is why: return ""
        if len(multi_choice)>=1:
            pass
        else:
            return ("")
        
    #dict_cluster_normalizedMedian = {}
        #multi_choice = i_multi_choice.value
        #clusters_for_ranking =  i_clusters_for_ranking.value
        df_distance_comp = self.df_distance_comp.copy()
        df_distance_comp = df_distance_comp[df_distance_comp["Experiment"].isin(multi_choice)]
        df_distance_comp = df_distance_comp[df_distance_comp["Cluster"].isin(clusters_for_ranking)]
        
        df_quantified_cluster = df_distance_comp.reset_index()
        df_quantified_cluster = df_distance_comp.drop_duplicates(subset=["Cluster", "Experiment"]).set_index(["Cluster", 
                                                                                                                "Experiment"])["distance"].unstack("Cluster")
        self.df_quantified_cluster = df_quantified_cluster.notnull().replace({True: "x", False: "-"})
        
        
        dict_quantified_cluster = {}
        dict_cluster_normalizedMedian_ref = {}
        dict_median_distance_ranking = {}
        for cluster in clusters_for_ranking:
            try:
                df_cluster = df_distance_comp[df_distance_comp["Cluster"]==cluster]
                cluster_quantitity = df_cluster["Gene names"].unique().size
                if  cluster_quantitity>= 5:
                    dict_quantified_cluster[cluster] = cluster_quantitity
                    all_median_one_cluster_several_exp = {}
                    #ref = df_cluster["distance"].median()
                    for exp in multi_choice:
                        median = df_cluster[df_cluster["Experiment"]==exp]["distance"].median()
                        all_median_one_cluster_several_exp[exp] = float(median)
                        #new
                        #if exp == ref_exp:
                        #    ref = median
                    ref = np.median(list(all_median_one_cluster_several_exp.values()))
                    dict_median_distance_ranking[cluster] = all_median_one_cluster_several_exp
                    
                    median_ranking_ref = {exp: median/ref for exp, median in all_median_one_cluster_several_exp.items()}
                    dict_cluster_normalizedMedian_ref[cluster] = median_ranking_ref
                else:
                    continue
            except:
                continue
        
        self.cluster_above_treshold = dict_quantified_cluster.keys()
        self.df_quantified_cluster2 = pd.DataFrame.from_dict({"Number of PG per Cluster":dict_quantified_cluster}).T
        
        df_cluster_normalizedMedian_ref = pd.DataFrame(dict_cluster_normalizedMedian_ref)
        df_cluster_normalizedMedian_ref.index.name="Experiment"
        df_cluster_normalizedMedian_ref.rename_axis("Cluster", axis=1, inplace=True)
        
        #median makes a huge differnece, improves result of DIA, MQ, libary
        df_RelDistanceRanking = pd.concat([df_cluster_normalizedMedian_ref.median(axis=1), df_cluster_normalizedMedian_ref.sem(axis=1)], axis=1, 
                                        keys=["Distance Ranking (rel, median)", "SEM"]).reset_index().sort_values("Distance Ranking (rel, median)")
        self.sorting_list = multi_choice#list(df_RelDistanceRanking["Experiment"])
        
        ranking_sum = df_cluster_normalizedMedian_ref.sum(axis=1).round(2)
        ranking_sum.name = "Normalized Median - Sum"
        df_ranking_sum = ranking_sum.reset_index()
        
        #ranking_product = df_cluster_normalizedMedian.product(axis=1).round(2)
        #ranking_product.name = "Normalized Median - Product"
        #df_globalRanking = pd.concat([pd.DataFrame(ranking_sum), pd.DataFrame(ranking_product)], axis=1).reset_index()
        
        df_cluster_normalizedMedian_ref = df_cluster_normalizedMedian_ref.stack("Cluster")
        df_cluster_normalizedMedian_ref.name="Normalized Median"
        df_cluster_normalizedMedian_ref = df_cluster_normalizedMedian_ref.reset_index()
        self.df_cluster_normalizedMedian_ref = df_cluster_normalizedMedian_ref
        df_cluster_normalizedMedian_ref = df_cluster_normalizedMedian_ref.assign(Experiment_lexicographic_sort = pd.Categorical(df_cluster_normalizedMedian_ref["Experiment"], categories=self.sorting_list, ordered=True))
        df_cluster_normalizedMedian_ref.sort_values("Experiment_lexicographic_sort", inplace=True)
        
        
        if collapse_cluster == False:
            
            fig_ranking = px.bar(df_cluster_normalizedMedian_ref, 
                                x="Cluster", 
                                y="Normalized Median", 
                                color="Experiment", 
                                barmode="group", 
                                title="Ranking - normalization to reference experiments the median across all experiments for each cluster",
                                template="simple_white"
                                )
            
            fig_ranking.update_xaxes(categoryorder="total ascending")
                
            fig_ranking.update_layout(autosize=False,
                                    width=1200 if len(multi_choice)<=3 else 300*len(multi_choice),
                                    height=500,
                                    template="simple_white"
                                    )
            return fig_ranking
 
        else:
            if ranking_boxPlot == "Bar plot - median":
                fig_globalRanking = px.bar(df_RelDistanceRanking.sort_values("Distance Ranking (rel, median)"), 
                                            x="Experiment",
                                            y="Distance Ranking (rel, median)", 
                                            title="Median manhattan distance distribution for all protein clusters (n>=5 per cluster)",# - median of all individual normalized medians - reference experiment is the median across all experiments for each cluster",
                                            error_x="SEM", error_y="SEM", 
                                            color="Experiment", 
                                            template="simple_white")
                    
        
                                                
            if ranking_boxPlot == "Box plot": 
                fig_globalRanking = px.box(df_cluster_normalizedMedian_ref,
                                        x="Experiment",
                                        y="Normalized Median", 
                                        title="Median manhattan distance distribution for all protein clusters (n>=5 per cluster)",# "Ranking - median of all individual normalized medians - reference is the median across all experiments for each cluster",
                                        color="Experiment",
                                        points="all",
                                        template="simple_white",
                                        hover_name="Cluster")
                #return pn.Column(pn.Row(fig_globalRanking), pn.Row(fig_globalRanking2))
            else:
                fig_globalRanking = px.bar(df_ranking_sum.sort_values("Normalized Median - Sum"), 
                    x="Experiment",
                    template="simple_white",
                    y="Normalized Median - Sum", 
                    title="Ranking - median of all individual normalized medians - reference is the median across all experiments for each cluster",
                    color="Experiment")
            
            fig_globalRanking.update_layout(autosize=False,
                                            width=250*len(multi_choice),
                                            height=500,
                                            template="simple_white"
                                            )
            
            return fig_globalRanking
        
    
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
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined[df_quantity_pr_pg_combined["Experiment"].isin(multi_choice)]
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined.assign(Experiment_lexicographic_sort = pd.Categorical(df_quantity_pr_pg_combined["Experiment"],
                                                                                                                        categories=self.sorting_list, ordered=True))
        #df_quantity_pr_pg_combined.sort_values("Experiment_lexicographic_sort", inplace=True)
        df_quantity_pr_pg_combined.sort_values(["Experiment_lexicographic_sort", "type"], ascending=[True, False], inplace=True)
        #df_quantity_pr_pg_combined.sort_values("type", ascending=False, inplace=True)
        
        layout = go.Layout(barmode="overlay", 
          #xaxis_tickangle=90, 
          autosize=False,
          width=300*len(multi_choice),
          height=400,
          #xaxis=go.layout.XAxis(linecolor="black",
          #                      linewidth=1,
          #                      title="Experiment",
          #                      mirror=True),
          #yaxis=go.layout.YAxis(linecolor="black",
          #                      linewidth=1,
          #                      title="#",
          #                      mirror=True),
          #legend = dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
          template="simple_white")
        
        fig_quantity_pg = px.bar(df_quantity_pr_pg_combined, x="filtering", y="number of protein groups", color="type", barmode="overlay", labels={"Experiment":"", "filtering":""}, 
                                 facet_col="Experiment",template="simple_white", opacity=1).for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig_quantity_pg.update_layout(layout, title="Number of Protein Groups")

        #fig_quantity_pg = go.Figure()
        #for t in df_quantity_pr_pg_combined["type"].unique():
        #    plot_df = df_quantity_pr_pg_combined[df_quantity_pr_pg_combined["type"] == t]
        #    fig_quantity_pg.add_trace(go.Bar(
        #        x=[plot_df["Experiment"], plot_df["filtering"]],
        #        y=plot_df["number of protein groups"],
        #        name=t))
        #fig_quantity_pg.update_layout(title="Number of Protein Groups")
        #fig_quantity_pg.update_layout(layout)
         
        fig_quantity_pr = px.bar(df_quantity_pr_pg_combined, x="filtering", y="number of profiles", color="type", barmode="overlay", labels={"Experiment":"", "filtering":""}, 
                                 facet_col="Experiment",template="simple_white", opacity=1).for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig_quantity_pr.update_layout(layout, title="Number of Profiles" )

        
        #fig_quantity_pr = go.Figure()
        #for t in df_quantity_pr_pg_combined["type"].unique():
        #    plot_df = df_quantity_pr_pg_combined[df_quantity_pr_pg_combined["type"] == t]
        #    fig_quantity_pr.add_trace(go.Bar(
        #        x=[plot_df["Experiment"], plot_df["filtering"]],
        #        y=plot_df["number of profiles"],
        #        name=t))
        #fig_quantity_pr.update_layout(title="Number of Profiles")
        #fig_quantity_pr.update_layout(layout) 
        #
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
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined[df_quantity_pr_pg_combined["Experiment"].isin(multi_choice)].sort_values("filtering")
        df_quantity_pr_pg_combined = df_quantity_pr_pg_combined.assign(Experiment_lexicographic_sort = pd.Categorical(df_quantity_pr_pg_combined["Experiment"],
                                                                                                                        categories=self.sorting_list, ordered=True))
        #df_quantity_pr_pg_combined.sort_values("Experiment_lexicographic_sort", inplace=True)
        df_quantity_pr_pg_combined.sort_values(["Experiment_lexicographic_sort", "filtering"], inplace=True)

        fig_pr_dc = px.bar(df_quantity_pr_pg_combined, x="type", y="data completeness of profiles", color="filtering", barmode="overlay", labels={"Experiment":"", "type":""}, 
                                 facet_col="Experiment",template="simple_white", opacity=1).for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        #fig_pr_dc.update_layout(layout)
        #fig_pr_dc = go.Figure()
        #for t in df_quantity_pr_pg_combined["filtering"].unique():
        #    plot_df = df_quantity_pr_pg_combined[df_quantity_pr_pg_combined["filtering"] == t]
        #    fig_pr_dc.add_trace(go.Bar(
        #        x=[plot_df["Experiment"], plot_df["type"]],
        #        y=plot_df["data completeness of profiles"],
        #        name=t))
            
        fig_pr_dc.update_layout(#barmode="overlay", 
                                         #xaxis_tickangle=90, 
                                         title="Data Completeness of Profiles",
                                         autosize=False,
                                         width=300*len(multi_choice),
                                         height=400,
                                         #xaxis=go.layout.XAxis(linecolor="black",
                                         #                      linewidth=1,
                                         #                      title="Experiment",
                                         #                      mirror=True),
                                         #yaxis=go.layout.YAxis(linecolor="black",
                                         #                      linewidth=1,
                                         #                      #title="",
                                         #                      mirror=True),
                                         #legend = dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
                                         template="simple_white")
        
        return fig_pr_dc
    
    
    def venn_sections(self, multi_choice_venn=["Exp1"]):
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
                unique_proteins_total: dict, key: Experiment name, value: unique protein (groups)
        
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
        
        sets_uniqueProteins = [set(self.unique_proteins_total[i]) for i in multi_choice_venn]
        num_combinations = 2 ** len(sets_uniqueProteins)
        bit_flags = [2 ** n for n in range(len(sets_uniqueProteins))]
        flags_zip_sets = [z for z in zip(bit_flags, sets_uniqueProteins)]
    
        combo_sets = []
        overlapping_ids = []
        experiments = []
        #dictio = {}
        for bits in range(num_combinations - 1, 0, -1):
            include_sets = [s for flag, s in flags_zip_sets if bits & flag]
            exclude_sets = [s for flag, s in flags_zip_sets if not bits & flag]
            combo = set.intersection(*include_sets)
            combo = set.difference(combo, *exclude_sets)
            tag = "".join([str(int((bits & flag) > 0)) for flag in bit_flags])
            
            experiment_decoded = []
            for digit, exp in zip(list(tag), multi_choice_venn):
                if digit=="0":
                    continue
                else:
                    experiment_decoded.append(exp)
            #dictio[len(combo)] = experiment_decoded
            if len(multi_choice_venn)>3:
                if len(combo)>300:
                    overlapping_ids.append(len(combo))
                    experiments.append(experiment_decoded)
                else:
                    continue
            else:
                overlapping_ids.append(len(combo))
                experiments.append(experiment_decoded)
            #combo_sets.append((tag, len(combo)))
            
        figure_UpSetPlot = plt.Figure()
        series_UpSetPlot = from_memberships(experiments, data=overlapping_ids)
        plot(series_UpSetPlot, fig=figure_UpSetPlot, show_counts="%d")


        
        if len(multi_choice_venn) == 2:
            vd = venn2([set(self.unique_proteins_total[i]) for i in multi_choice_venn], 
                       set_labels=([i for i in multi_choice_venn]),
                       set_colors=("darksalmon", "darkgrey")
                      )
        elif len(multi_choice_venn) == 3:
            vd = venn3([set(self.unique_proteins_total[i]) for i in multi_choice_venn],
                       set_labels=([i for i in multi_choice_venn]),
                       set_colors=("darksalmon", "darkgrey","rosybrown"),
                       alpha=0.8
                      )
        
        else:
            im = "Venn diagram can be displayed for 3 Experiments or less"
            return im, figure_UpSetPlot

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
        
        return im, figure_UpSetPlot
    
    
    def dynamic_range_comparison(self, collapse_cluster=False, multi_choice=["Exp1", "Exp2"], ref_exp="Exp1"):
        """
        A box plot for desired experiments (multi_choice) and all protein clusters is generated displaying the dynamic range
        
        Args:
            self:
                multi_choice: list of experiment names 
                df_dynamicRange_combined: df, no index, column names: "Max", "Min", "Dynamic Range", "Cluster", "Experiment"
        
        Returns:
            fig_dynamic_range: bar plot, dynamic range of each protein cluster for desired experiments is displayed.
        """
        df_dynamicRange_combined = self.df_dynamicRange_combined.copy()
        df_dynamicRange_combined = df_dynamicRange_combined[df_dynamicRange_combined["Experiment"].isin(multi_choice)]
        df_dynamicRange_combined = df_dynamicRange_combined.assign(Experiment_lexicographic_sort = pd.Categorical(df_dynamicRange_combined["Experiment"],
                                                                                                                categories=self.sorting_list, ordered=True))
        
        df_dynamicRange_combined.sort_values(["Experiment_lexicographic_sort", "Dynamic Range"], inplace=True)
        
        fig_dynamic_range = px.bar(df_dynamicRange_combined,  
                                   x="Cluster", 
                                   y="Dynamic Range", 
                                   base="Min", 
                                   facet_row="Experiment", 
                                   template="simple_white",
                                   height=400*len(multi_choice),
                                   width=1200)
        
        df_dynamicRange_combined_ref = df_dynamicRange_combined.drop(["Experiment_lexicographic_sort"], axis=1)
        df_dynamicRange_combined_ref = df_dynamicRange_combined.set_index(["Cluster", "Experiment"], drop=False).unstack("Cluster")["Dynamic Range"]
        df_dynamicRange_combined_ref = df_dynamicRange_combined_ref.div(df_dynamicRange_combined_ref.xs(ref_exp))
        df_RelDynamicRange = pd.concat([df_dynamicRange_combined_ref.median(axis=1), df_dynamicRange_combined_ref.sem(axis=1)], axis=1, 
                                       keys=["Dynamic Range (rel, median)", "SEM"]).reset_index()
        
        if collapse_cluster == False:
            df_dynamicRange_combined_ref = df_dynamicRange_combined_ref.stack("Cluster")
            df_dynamicRange_combined_ref.name="Normalized Dynamic Range"
            df_dynamicRange_combined_ref = df_dynamicRange_combined_ref.reset_index()
            
            fig_RelDynamicRange = px.bar(df_dynamicRange_combined_ref, 
                                         x="Cluster", 
                                         y="Normalized Dynamic Range", 
                                         title="Dynamic Range - normalization to reference experiment: {}".format(ref_exp),
                                         barmode="group", 
                                         template="simple_white",
                                         color="Experiment")
            fig_RelDynamicRange.update_xaxes(categoryorder="total ascending")
            fig_RelDynamicRange.update_layout(autosize=False,
                                              width=1200 if len(multi_choice)<=3 else 300*len(multi_choice),
                                              height=500,
                                              template="simple_white"
                                               )
        else:
            fig_RelDynamicRange = px.bar(df_RelDynamicRange.sort_values("Dynamic Range (rel, median)"), 
                                         x="Experiment", 
                                         y="Dynamic Range (rel, median)", 
                                         error_x="SEM", error_y="SEM",
                                         template="simple_white",
                                         title="Dynamic Range - median of all individual normalized medians - reference experiment: {}".format(ref_exp),
                                         color="Experiment")
            fig_RelDynamicRange.update_layout(autosize=False,
                                                width=250*len(multi_choice),
                                                height=500,
                                                template="simple_white"
                                               )
            
            
        return pn.Column(pn.Row(fig_dynamic_range), pn.Row(fig_RelDynamicRange))
    
    
    def calculate_global_scatter(self, multi_choice, metric, consolidation):
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
        cons_functions = {
            "median": np.median,
            "average": np.mean,
            "sum": np.sum
        }
        metrics = {
            "euclidean distance": "euclidean",
            "manhattan distance": "manhattan",
            "1 - cosine correlation": "cosine",
            "1 - pearson correlation": lambda x,y: 1-np.corrcoef(x,y)[0][1]
        }
        
        # Option assertion
        assert consolidation in cons_functions.keys()
        assert metric in metrics.keys()
        
        # Filter experiments and intersection of proteins
        df = self.df_01_filtered_combined.loc[
            self.df_01_filtered_combined.index.get_level_values("Experiment").isin(multi_choice)].copy()
        
        df.index = df.index.droplevel("Exp_Map")
        df_across = df.unstack(["Experiment", "Map"]).dropna().stack(["Experiment", "Map"])
        nPG = df_across.unstack(["Experiment", "Map"]).shape[0]

        
        # Calculate and consolidate distances
        distances = pd.DataFrame()
        for exp in set(df_across.index.get_level_values("Experiment")):
            df_m = df_across.xs(exp, level="Experiment", axis=0)
            maps = list(set(df_m.index.get_level_values("Map")))
            distances_m = pd.DataFrame()
            for i,mapi in enumerate(maps):
                for j,mapj in enumerate(maps):
                    # only look at each comparison once
                    if j <= i:
                        continue
                    dist = pw.paired_distances(df_m.xs(mapi, level="Map", axis=0).values,
                                               df_m.xs(mapj, level="Map", axis=0).values,
                                               metric = metrics[metric])
                    dist = pd.Series(dist, name="_".join([mapi,mapj]))
                    distances_m = pd.concat([distances_m, dist], axis=1)
            distances_m.index = df_m.xs(maps[0], level="Map", axis=0).index
            distances = pd.concat([distances, pd.Series(distances_m.apply(cons_functions[consolidation], axis=1), name=exp)], axis=1)
        distances.index = distances_m.index
        
        self.distances = distances
        # Sort by median
        medians = distances.apply(np.median)
        
        distances = distances[multi_choice]
        #distances = distances[medians.rank().sort_values(ascending=False).index]
        #self.sorting_list = list(medians.rank().sort_values(ascending=False).index)
        
        #df_cluster_normalizedMedian_ref = df_cluster_normalizedMedian_ref.assign(Experiment_lexicographic_sort = pd.Categorical(df_cluster_normalizedMedian_ref["Experiment"], categories=self.sorting_list, ordered=True))
        #df_cluster_normalizedMedian_ref.sort_values("Experiment_lexicographic_sort", inplace=True)
        
        # Create and return plot
        plot = ff.create_distplot(distances.T.values, distances.columns, show_hist=False)
        plot.update_layout(title="Distribution of {} {}s, n = {}".format(metric, consolidation, nPG),
                           width=1500, height=600, template="simple_white")
        return plot

        
    
      
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
        
        global_SVM_dict_total = {}
        global_SVM_dict = {}
        for exp in self.json_dict.keys():
            try:
                df_SVM = pd.read_json(self.json_dict[exp]["Misclassification Matrix"])
                df_SVM["T: True group"] = df_SVM["T: True group"].str.replace(r'True: ', '')
            except KeyError:
                continue
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
                if df_SVM["T: True group"][i]!="Nuclear pore complex" and df_SVM["T: True group"][i]!="Large Protein Complex" and df_SVM["T: True group"][i]!="Actin binding proteins" :
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
            
            SVM_dict_total = {}
            SVM_dict_total["Avg. all clusters"] = {"Recall": total_recall, "F1": avg_F1_all_cluster} #total recall = marker prediction accuracy
            SVM_dict_total["Avg. all organelles"] = {"Recall": av_per_organelle_recall, "F1": avg_organelle_f1, "Precision": av_per_organelle_precision}
            SVM_dict_total["Membrane"] = {"Recall": membrane_recall}
            SVM_dict_total["Median. per organelle"] = {"Recall": median_per_organelle_recall}
            
            global_SVM_dict[exp] = SVM_dict
            global_SVM_dict_total[exp] = SVM_dict_total
            self.global_SVM_dict = global_SVM_dict
            self.global_SVM_dict_total = global_SVM_dict_total
            
        if global_SVM_dict=={}:
            self.cache_stored_SVM = False
            return
        else:
            df_clusterPerformance_global = pd.DataFrame.from_dict({(i,j): global_SVM_dict[i][j] 
                                    for i in global_SVM_dict.keys() 
                                    for j in global_SVM_dict[i].keys()},
                                orient='index')
            df_clusterPerformance_global.index.names = ["Experiment", "Type"]
            self.df_clusterPerformance_global = df_clusterPerformance_global.T 
            
            df_AvgClusterPerformance_global = pd.DataFrame.from_dict({(i,j): global_SVM_dict_total[i][j] 
                                    for i in global_SVM_dict_total.keys() 
                                    for j in global_SVM_dict_total[i].keys()},
                                orient='index')
            df_AvgClusterPerformance_global.index.names = ["Experiment", "Type"]
            self.df_AvgClusterPerformance_global = df_AvgClusterPerformance_global.T
            self.cache_stored_SVM = True
            return 
            
            
    def svm_plotting(self, multi_choice):
        """
        The markerperformance (line/scatter plot) as well as marker prediction accuracy (bar plot) is visuaized.
        
        Args:
            self: df_AvgClusterPerformance_global 
                  df_clusterPerformance_global
            multi_choice: list of experiment names
        """
        df_clusterPerformance_global = self.df_clusterPerformance_global
        df_AvgClusterPerformance_global = self.df_AvgClusterPerformance_global
        
        df_AvgAllCluster = df_AvgClusterPerformance_global.xs("Avg. all clusters", level='Type', axis=1)
        fig_markerPredictionAccuracy = go.Figure()#data=[go.Bar(x=df_test.columns, y=df_test.loc["Recall"])])
        for exp in multi_choice:
            fig_markerPredictionAccuracy.add_trace(go.Bar(x=[exp], y=[df_AvgAllCluster[exp].loc["Recall"]], name=exp))
        fig_markerPredictionAccuracy.update_layout(template="simple_white", #showlegend=False, 
                    title="Cluster Performance",
                    xaxis=go.layout.XAxis(linecolor="black",
                                        linewidth=1,
                                        mirror=True),
                    yaxis=go.layout.YAxis(linecolor="black",
                                        linewidth=1,
                                        title="F1 score - harmonic mean of recall and precision",
                                        mirror=True),
                    )
        
        fig_clusterPerformance = go.Figure()
        list_data_type = ["Avg. all clusters", "Avg. all organelles"]
        for i,exp in enumerate(multi_choice):
            df_clusterPerformance = df_clusterPerformance_global.xs(exp, level='Experiment', axis=1)
            df_AvgClusterPerformance = df_AvgClusterPerformance_global.xs(exp, level='Experiment', axis=1)
            fig_clusterPerformance.add_trace(go.Scatter(x=df_clusterPerformance.columns, y=df_clusterPerformance.loc["F1"], 
                                                    marker=dict(color=pio.templates["simple_white"].layout["colorway"][i]), name=exp))
            for  data_type in list_data_type:
                fig_clusterPerformance.add_trace(go.Scatter(x=[data_type], y=[df_AvgClusterPerformance[data_type].loc["F1"]],
                            mode="markers",
                            showlegend=False,
                            marker=dict(color=pio.templates["simple_white"].layout["colorway"][i])
                                ))
            fig_clusterPerformance.update_layout(template="simple_white", #showlegend=False, 
                            title="Cluster Performance",
                            xaxis=go.layout.XAxis(linecolor="black",
                                                linewidth=1,
                                                mirror=True),
                            yaxis=go.layout.YAxis(linecolor="black",
                                                linewidth=1,
                                                title="F1 score - harmonic mean of recall and precision",
                                                mirror=True),
                            )
        
        return fig_markerPredictionAccuracy, fig_clusterPerformance
        

        
    def __repr__(self):
        return str(self.__dict__)
        #return "This is a spatial dataset with {} lines.".format(len(self.df_original))

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
    
    #df_SVM = self.df_SVM.copy()
    #if hasattr(df_SVM, "keys") == True: 
    try:
        df_SVM = pd.read_json(df_SVM["Misclassification Matrix"])
        df_SVM = df_SVM.set_index("T: True group")[::-1]
    
    except:
        df_SVM = df_SVM.set_index("T: True group")[::-1]
    
    y_axis_label = df_SVM.index
    x_axis_label = df_SVM.columns
    data_svm = df_SVM.values
    fig_SVMheatmap = go.Figure()
    
    fig_SVMheatmap.add_trace(go.Heatmap(
                             z=data_svm,
                             x = x_axis_label,
                             y = y_axis_label,
                             colorscale=[
                                 [0.0, "green"],
                                 [0.01, "white"],
                                 [1.0, "red"]
                             ],
                         ))
               
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
            if data_type == "0/1 normalized data" and exp_name == list(json_dict.keys())[0]:
                df_01_combined = pd.read_json(json_dict[exp_name][data_type])
                df_01_combined = df_01_combined.set_index(["Gene names", "Protein IDs", "Compartment"]).copy()
                df_01_combined.drop([col for col in df_01_combined.columns if not col.startswith("normalized profile")])
                df_01_combined.columns = pd.MultiIndex.from_tuples([el.split("?") for el in df_01_combined.columns], names=["Set", "Map", "Fraction"])
                df_01_combined.rename(columns = {"normalized profile":exp_name}, inplace=True)
    
            elif data_type == "0/1 normalized data" and exp_name != list(json_dict.keys())[0]:
                df_01_toadd = pd.read_json(json_dict[exp_name][data_type])
                df_01_toadd = df_01_toadd.set_index(["Gene names", "Protein IDs", "Compartment"]).copy()
                df_01_toadd.drop([col for col in df_01_toadd.columns if not col.startswith("normalized profile")])
                df_01_toadd.columns = pd.MultiIndex.from_tuples([el.split("?") for el in df_01_toadd.columns], names=["Set", "Map", "Fraction"])
                df_01_toadd.rename(columns = {"normalized profile":exp_name}, inplace=True)
                df_01_combined = pd.concat([df_01_combined, df_01_toadd], axis=1)
                                   
    df_01_combined.columns.names = ["Experiment", "Map", "Fraction"]
    df = df_01_combined.stack(["Experiment", "Map"]).dropna(axis=0)
    df = df.div(df.sum(axis=1), axis=0)
    index_ExpMap = df.index.get_level_values("Experiment")+"_"+df.index.get_level_values("Map")
    index_ExpMap.name = "Exp_Map"
    df.set_index(index_ExpMap, append=True, inplace=True)
    
    df.index = df.index.droplevel(["Map", "Experiment"])
    df = df.stack("Fraction").unstack(["Exp_Map", "Fraction"])
    df.columns = ["_".join(col) for col in df.columns.values]
    
    return df
   
