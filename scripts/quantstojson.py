def benchmark_download_getsheet(comp, sheet, mode="csv"):
    if sheet == "0-1 normalized data":
        out = comp.df_01_filtered_combined.copy()
        out.index = out.index.droplevel("Exp_Map")
        out = out.unstack(["Experiment", "Map"])
        out.columns = ["_".join(el) for el in out.columns.reorder_levels(["Experiment", "Map", "Fraction"]).values]
    elif sheet == "pca coordinates":
        out = comp.df_pca.copy()
        out.index = out.index.droplevel("merge type")
    elif sheet == "complex scatter":
        out = comp.df_distance_comp.copy()
        if mode=="csv":
            out["Cluster"] = ['"'+el+'"' for el in out["Cluster"]]
        out = out.set_index(
            ["Cluster", "Gene names", "Protein IDs", "Compartment", "Experiment", "Map"]).drop(
            ["Exp_Map", "merge type"], axis=1).unstack(["Experiment", "Map"]).copy()
        out.columns = ["_".join(el) for el in out.columns.values]
    elif sheet == "reproducibility":
        out = comp.distances.copy()
    elif sheet == "protein id alignment":
        out = comp.id_alignment.copy()
        out.columns = [" ".join(el) for el in out.columns.values]
    else:
        raise KeyError(sheet)
    return out

if __name__ == "__main__":
    
    import argparse
    
    # Script arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='Processing and quality control of profiling proteomics data files (protein groups to json).',
                                     epilog='''
Simple usage examples:

Get 0-1 normalized data for single MQ output file with reannotated genes fresh from uniprot:
python quantstojson.py "LFQ6 - MQ" proteinGroups.txt . -p -r -n --rsource uniprot

Process a folder of files for benchmarking generating one json file per input:
python quantstojson.py "LFQ6 - MQ" inputfiles jsonfiles -f -i
                                 ''')
    
    parser.add_argument('acquisition',
                        help="Acquisition mode of the files (must all be identical), default for LFQ is 'LFQ6 - MQ'")
    parser.add_argument('files', nargs="+",
                        help='Input files, or directory if -d is specified')
    parser.add_argument('outputfolder',
                        help='Destination folder for output file(s)')
    parser.add_argument('--pattern', default=".* (?P<rep>.*)_(?P<frac>.*)",
                        help='Regular expression for parsing experiment names into fractions and maps. Must contain capture groups named "rep" and "frac", if applicable also "cond". default=".* (?P<rep>.*)_(?P<frac>.*)"')
    parser.add_argument('--organism', default="Homo sapiens - Uniprot",
                        help='Organism specifying complex and organelle annotations. Has to match a filename from the package ressources. default="Homo sapiens - Uniprot"')
    parser.add_argument('--name', '--names', nargs="+", default=["filename"],
                        help="List of experiment names, defaults to using the file names.")
    parser.add_argument('--comment', '--comments', nargs="+", default=["filename"],
                        help="List of comments, defaults to using the file names.")
    parser.add_argument('--directory', '-d', '-D', action="store_true",
                        help="Interpret input files as a folder and process all compatible files in that folder")
    parser.add_argument('--benchmark', '-b', '-B', action="store_true",
                        help='Set this flag to run the benchmarking analysis of all selected files. See argument group below for configuration.')
    parser.add_argument('--individualoutput', '-i', '-I', action="store_true",
                        help='Set this flag to create a single json file per input file instead of one collection.')
    parser.add_argument('--nojson', '-n', '-N', action="store_true",
                        help='Set this flag to create no json file at all.')
    parser.add_argument('--summary', '-s', '-S', action="store_true",
                        help='NOT IMPLEMENTED Set this flag to write a quality summary file to the output destination.')
    parser.add_argument('--tableoutput', '-t', '-T', action="store_true",
                        help='NOT IMPLEMENTED Set this flag to write various tabular output along the way to the output destination.')
    parser.add_argument('--overwrite', '-o', '-O', action="store_true",
                        help='Set this flag to overwrite existing files in the output destination.')
    parser.add_argument('--verbose', '-v', '-V', action="store_true",
                        help='Set this flag to write progress to STDOUT.')
    
    group_lfq = parser.add_argument_group('Quality filters for LFQ input')
    group_lfq.add_argument('--cons', default=4, help='Consecutive values requried per profile, default=4')
    group_lfq.add_argument('--msmsavg', default=2, help='Average number of MS/MS counts per profile point, default=2')
    
    group_SILAC = parser.add_argument_group('Quality filters for SILAC input')
    group_SILAC.add_argument('--count', default=2, help='Minimal ratio count per quant, default=2')
    group_SILAC.add_argument('--variability', default=30, help='Minimal ratio variability per quant, default=30')
    
    group_custom = parser.add_argument_group('Column specification for custom data input')
    group_custom.add_argument('--pgcolumn', default="Protein IDs", help='Name of column containing group groups, default="Protein IDs"')
    group_custom.add_argument('--genecolumn', default="Gene names", help='Name of column containing gene names, default="Gene names"')
    group_custom.add_argument('--datacolumns', nargs="+", default=".*",
                            help='Names or regular expression to match columns containing data (pgcolumn and genecolumn are excluded automatically), default=".*"')
    
    group_reannotate = parser.add_argument_group("Reannotate genes from annotation files or from uniprot.")
    group_reannotate.add_argument('--reannotate', '-r', '-R', action="store_true", help="Include this to reannotate genes")
    group_reannotate.add_argument('--rsource', default="fasta_headers", choices=["fasta_headers", "uniprot", "tsv"],
                                help='Type of reannotation source, default="fasta_headers"')
    group_reannotate.add_argument('--rfile', default="Human - uniprot reference_swissprot.txt",
                                help='Unless source is uniprot, path to the annotation file. This can be a filename of annotations stored in the dompas package, default="Human - uniprot reference_swissprot.txt"')
    
    group_benchmark = parser.add_argument_group("Reannotate genes from annotation files or from uniprot.")
    group_benchmark.add_argument('--bname',
                                 help="Name of this benchmark, used in filenames. Defaults to current date.")
    group_benchmark.add_argument('--excel', '-e', '-E', action="store_true",
                                 help='Generate excel output for the benchmark (this is a timeconsuming step).')
    group_benchmark.add_argument('--figures', '-f', '-F', action="store_true",
                                 help='Generate plotly html figures for the benchmark.')
    
    args = parser.parse_args()
    
    if args.verbose:
        print("\n----\nThis is the configuration:\n----")
        print(args)
    
    if args.verbose:
        print("\n----\nLoading base libraries and validating input\n----")
    import os
    from warnings import warn
    import traceback
    
    def uwarn(message):
        warn(message)
    
    for f in args.files:
        if not os.path.exists(f):
            raise FileNotFoundError(f)
    
    if args.directory and len(args.files) > 1:
        raise ValueError('Specify either -f and a single value for files or multiple filenames, not both')
    
    if args.directory:
        folder = args.files[0]
        args.files = [os.path.join(folder,el) for el in os.listdir(args.files[0]) if any([el.endswith(ending) for ending in ['.txt', '.xls', '.csv']])]
    
    if args.verbose:
        print(f'''Will try to process these files:
    {",   ".join(args.files)}''')
    
    if os.path.exists(args.outputfolder):
        if not args.overwrite:
            if not args.individualoutput and not args.nojson and os.path.exists(os.path.join(args.outputfolder, "collection.json")) and not args.benchmark:
                raise FileExistsError("Exiting because an output collection already exists in the destination folder. Either choose a different outputfolder, set -I to generate individual files per input, or -O to overwrite the existing file.")
            uwarn("Warning: The specified output folder already exists and any conflicting files will not be replaced by updated versions")
            
        if args.overwrite and args.verbose:
            uwarn("Warning: The specified output folder already exists and any conflicting files will be overwritten")
    
    if not os.path.exists(args.outputfolder):
        os.mkdir(args.outputfolder)
    
    if args.name == ['filename']:
        args.name = [os.path.basename(el).split(".")[0] for el in args.files]
    assert len(args.name) == len(args.files)
    if args.comment == ['filename']:
        args.comment = [os.path.basename(el).split(".")[0] for el in args.files]
    assert len(args.comment) == len(args.files)
    
    if args.verbose:
        print("\n----\nLoading main library and setting up the analysis\n----")
    
    import domaps
    import pandas as pd
    from io import BytesIO
    import json
    reannotation_source = {}
    
    if args.reannotate:
        import pkg_resources
        reannotation_source["source"] = args.rsource
        if args.rsource == "fasta_headers":
            if os.path.exists(args.rfile):
                with open(args.rfile, 'r') as rfile:
                    idmapping = [el.strip() for el in rfile.readlines()]
            else:
                idmapping = [el.decode('UTF-8').strip()
                             for el in pkg_resources.resource_stream(
                    "domaps", 'annotations/idmapping/{}'.format(args.rfile)).readlines()]
            reannotation_source["idmapping"] = idmapping
        elif args.rsource == "tsv":
            if os.path.exists(args.rfile):
                idmapping = pd.read_csv(args.rfile, sep='\t', usecols=["Entry", "Gene names  (primary )"])
            else:
                idmapping = pd.read_csv(pkg_resources.resource_stream(
                    "domaps", 'annotations/idmapping/{}'.format(args.rfile)), sep='\t',
                                        usecols=["Entry", "Gene names  (primary )"])
            reannotation_source["idmapping"] = idmapping
        
        if args.verbose:
            print("This will be used for reannotation:")
            print(reannotation_source["source"])
            if reannotation_source["source"] != "uniprot":
                print(reannotation_source["idmapping"][0:5])
    
    if args.verbose:
        print("\n----\nStarting file processing\n----")
    
    if args.benchmark or (not args.individualoutput and not args.nojson):
        collection = {}
    
    failed, success = [], []
    for file, name, comment in zip(args.files, args.name, args.comment):
        if args.verbose:
            print(f'\n Processing {file} ...')
        
        # check if output file already exists
        if args.individualoutput and os.path.exists(os.path.join(args.outputfolder, name+".json")):
            if args.overwrite and not args.nojson:
                uwarn(f"Will overwrite file {name}.json")
            elif args.benchmark:
                if args.verbose:
                    print(f"Loading previous analysis for {name} to be included in the benchmark.")
                with open(os.path.join(args.outputfolder, name+".json"), 'r') as ifile:
                    qc_json = json.load(ifile)
                collection[name] = qc_json[name]
                success.append(file)
                continue
            else:
                if args.verbose:
                    print(f"Analysis output for {name} already exists, thus it will be skipped.")
                success.append(file)
                continue
        
        try:
            qc = domaps.SpatialDataSet(
                file, name, args.acquisition,
                comment=comment, name_pattern=args.pattern, organism=args.organism,
                reannotate_genes=args.reannotate, reannotation_source=reannotation_source,
                consecutiveLFQi=args.cons, summed_MSMS_counts=args.msmsavg,
                RatioHLcount=args.count, RatioVariability=args.variability,
                custom_columns={"ids": args.pgcolumn,
                                "genes": args.genecolumn,
                                "data": args.datacolumns}
            )
        except:
            if args.verbose:
                uwarn(f'File {file} could not be initialized, continuing with next file.')
                print(traceback.format_exc())
            failed.append(file)
            continue
        try:
            if args.verbose:
                qc.run_pipeline(progressbar="STDOUT")
            else:
                qc.run_pipeline(progressbar=None)
        except:
            if args.verbose:
                uwarn(f'File {file} could not be processed, continuing with next file.')
                print(traceback.format_exc())
            failed.append(file)
            continue
        
        # Store output
        if args.benchmark or (not args.individualoutput and not args.nojson):
            collection[name] = qc.analysed_datasets_dict[name]
        if args.individualoutput and not args.nojson:
            with open(os.path.join(args.outputfolder, name+".json"), 'w') as ofile:
                json.dump(
                    {name: qc.analysed_datasets_dict[name]},
                    ofile,
                    indent=4,
                    sort_keys=True
                )
        
        success.append(file)
    
    if not args.individualoutput and not args.nojson:
        if args.verbose:
            print("Writing output collection ...")
        with open(os.path.join(args.outputfolder, "collection.json"), 'w') as ofile:
            json.dump(
                collection,
                ofile,
                indent=4,
                sort_keys=True
            )
    
    if args.verbose:
        if len(failed) > 0:
            print("\n----\nThese files failed to be processed:\n----")
            print(f'{",    ".join(failed)}')
        
        if len(success) > 0:
            print("\n----\nThese files were successfully processed:\n----")
            print(f'{",    ".join(success)}')
        else:
            print("\n----\nSomething went utterly wrong so no files were processed successfully\n----")
    
    if args.benchmark:
    
        if (
            (
             os.path.exists(os.path.join(args.outputfolder, args.bname+"_depth.html")) and args.figures
             ) or (
             os.path.exists(os.path.join(args.outputfolder, args.bname+"_data.xlsx")) and args.excel
             )
            ) and not args.overwrite:
            uwarn("The benchmark will not be run, because some files already exist and overwrite was set to False.")
        
        else:
            experiments = list(collection.keys())
            if args.verbose:
                print("\n----\nStarting benchmarking analysis with default parameters\n----")
            comp = domaps.SpatialDataSetComparison(ref_exp=experiments[0])
            comp.json_dict = {k: collection[k] for k in experiments}
            if args.verbose:
                print("Aligning protein groups ...")
            comp.read_jsonFile()
            if args.verbose:
                print("Calculating inter-map scatter ...")
            comp.calculate_global_scatter(metric="manhattan distance to average profile", consolidation="average")
            if args.verbose:
                print("Calculating intra-complex scatter ...")
            comp.calc_biological_precision()
            comp.get_complex_coverage()
            if args.verbose:
                print("Performing PCA analysis ...")
            comp.perform_pca_comparison()
            
            
            if args.figures:
                plotly_config={
                    'toImageButtonOptions': {
                            'format': 'svg', # one of png, svg, jpeg, webp
                            'filename': 'figure'
                    }
                }
                if args.verbose:
                    print("\n Writing figures ...")
                comp.quantity_pr_pg_barplot_comparison(experiments)[0].write_html(os.path.join(args.outputfolder, args.bname+"_depth.html"), config=plotly_config)
                comp.plot_reproducibility_distribution(show_full=True, q=0.5, x_cut=0.8).write_html(os.path.join(args.outputfolder, args.bname+"_intermap scatter.html"), config=plotly_config)
                comp.plot_biological_precision(experiments, reference=experiments[0],
                                            clusters_for_ranking=list(comp.coverage_lists[0].keys()))[0].write_html(os.path.join(args.outputfolder, args.bname+"_intracomplex scatter.html"), config=plotly_config)
                comp.venn_sections(experiments, omit_size=0)[3].savefig(os.path.join(args.outputfolder, args.bname+"_upsetplot.svg"))
                comp.plot_global_pca_comparison(multi_choice=experiments).write_html(os.path.join(args.outputfolder, args.bname+"_PCA.html"), config=plotly_config)
            
            if args.excel:
                if args.verbose:
                    print("\n Writing excel ...")
                sheets = [
                    "0-1 normalized data",
                    "pca coordinates",
                    "complex scatter",
                    "reproducibility",
                    "protein id alignment",
                ]
                mode="xlsx"
                out_dict = dict()
                for sheet in sheets:
                    out_dict[sheet] = benchmark_download_getsheet(comp, sheet, mode=mode)
                excel = pd.ExcelWriter(os.path.join(args.outputfolder, args.bname+"_data.xlsx"), engine_kwargs = {"data_only": True})
                for sheet, df in out_dict.items():
                    df.reset_index().T.reset_index().T.to_excel(
                        excel, sheet_name=sheet,
                        merge_cells=False, index=False, header=False)
                excel.save()
            
    if args.verbose:
        print("\n----\nDONE")