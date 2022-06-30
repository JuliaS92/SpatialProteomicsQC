def getcounts(pr_matrix):
    return pr_matrix[[col for col in pr_matrix.columns if ":" in col or col =="Protein.Group"]].melt("Protein.Group").dropna().groupby(["variable", "Protein.Group"]).size().unstack("variable")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Convert DIANN output to MQ style output for loading into domaps tools')
    parser.add_argument('directory', help='Directory containing DIANN output. Must include .*pg_matrix.tsv and *.pr_matrix.tsv')
    parser.add_argument('--type', default="pg", help='Type of data to convert, defaults to protein group quantities "pg". For upwards compatibility only.')
    parser.add_argument('--regex', default="(.*)", help='Regular expression to apply to basename of raw files to convert them to experiment names using re.findall()[0]')
    parser.add_argument('--verbose', '-v', '-V', action="store_true", help='Talk to me!')
    
    args = parser.parse_args()
    
    import os
    import pandas as pd
    import re
    
    files = os.listdir(args.directory)
    
    if args.type == "pg":
        pg_files = [f for f in files if re.match("^.*(?<!first-pass).pg_matrix.tsv$", f)]
        pr_files = [f for f in files if re.match("^.*(?<!first-pass).pr_matrix.tsv$", f)]
        
        for pg, pr in zip(pg_files, pr_files):
            
            if args.verbose:
                print(f'Converting {pg} ...')
                
            name = pg.split("pg_matrix")[0]
            pg = pd.read_csv(os.path.join(args.directory, pg), sep='\t')
            pr = pd.read_csv(os.path.join(args.directory, pr), sep='\t')
            
            pr_count = getcounts(pr).reset_index()
            
            pg.columns = [col if ":" not in col else "LFQ intensity "+re.findall(args.regex, os.path.basename(col).split(".")[0])[0] for col in pg.columns]
            pr_count.columns = [col if ":" not in col else "MS/MS count "+re.findall(args.regex, os.path.basename(col).split(".")[0])[0] for col in pr_count.columns]
            
            out = pg.merge(pr_count, on="Protein.Group", how="left")
            
            out.rename({"Protein.Group": "Majority protein IDs", "Genes": "Gene names", "Protein.Ids": "Protein IDs"}, axis=1, inplace=True)
            out.insert(loc=out.shape[1], column="Potential contaminant", value="")
            out.insert(loc=out.shape[1], column="Reverse", value="")
            out.insert(loc=out.shape[1], column="Only identified by site", value="")
            
            outpath = os.path.join(args.directory, name+"MQfordomaps.txt")
            
            if args.verbose:
                print(f'Writing {outpath} with columns {", ".join(out.columns)} ...')
            
            out.to_csv(outpath, sep="\t", index=False)