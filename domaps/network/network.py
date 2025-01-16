import pandas as pd
import numpy as np
import re
import networkx as nx
from sklearn.metrics import pairwise as pw

pd.set_option("mode.chained_assignment", None)
import panel as pn
import hvplot.pandas
import hvplot.networkx as hvnx
import holoviews as hv
from holoviews import opts
import matplotlib.pyplot as plt
import time


def format_data(sp):
    df = sp.df_01_stacked["normalized profile"].unstack("Fraction").copy()
    min_correl = pd.Series(sp.correlations.min(axis=1), name="Min correl")
    min_correl.reset_index(
        [el for el in min_correl.index.names if el != "Protein IDs"],
        drop=True,
        inplace=True,
    )
    df = (
        df.join(min_correl).dropna().set_index("Min correl", append=True).unstack("Map")
    )
    df_core = df.reset_index(
        [el for el in df.index.names if el not in ["Protein IDs", "Min correl"]],
        drop=True,
    )
    meta_dict = {}
    meta_dict["gene_id"] = {}
    for g, i in zip(
        df.index.get_level_values("Gene names"),
        df.index.get_level_values("Protein IDs"),
    ):
        if g not in meta_dict["gene_id"].keys():
            meta_dict["gene_id"][g] = i
        else:
            if type(meta_dict["gene_id"][g]) == list:
                meta_dict["gene_id"][g].append(i)
            else:
                meta_dict["gene_id"][g] = [meta_dict["gene_id"][g], i]
    meta_dict["id_meta"] = (
        df.index.to_frame()
        .set_index("Protein IDs")[["id", "Gene names", "Compartment"]]
        .to_dict("index")
    )

    return df_core, meta_dict


## single query dataset
def hv_query_df(
    df_core,
    meta_dict,
    gene,
    protein="average",
    dist="Manhattan",
    min_corr=0,
    size=50,
    rank_tolerance=25,
    perc_area=250,
):
    size += 1
    if gene.upper() not in meta_dict["gene_id"].keys():
        return "Gene {} not in dataset".format(gene)

    # Filter data
    df = df_core.loc[
        [corr >= min_corr for corr in df_core.index.get_level_values("Min correl")]
    ].copy()
    df.reset_index(["Min correl"], drop=True, inplace=True)

    # Retrieve query
    query_id = meta_dict["gene_id"][gene]
    msg_append = ""
    query_profile = df.loc[query_id, :]
    if type(query_id) == list:
        if protein == "average":
            query_profile = query_profile.mean(axis=0)
        else:
            query_profile = df.loc[protein, :]
    elif query_id not in df.index:
        return "The query protein {} was filtered out.".format(gene)

    # Calculate difference profiles and sum within replicates
    df_diff = df - query_profile
    df_diff = df_diff.transform(np.abs)
    df_diff_rep = df_diff.groupby("Map", axis=1).sum()

    # Score distances
    if dist == "Manhattan":
        dists = df_diff.apply(sum, axis=1)
    dists.name = "Distance measure"
    dists = dists.sort_values()

    # Calculate rank ranges across replicates
    df_diff_rep_rank = df_diff_rep.rank()
    df_diff_rep_rank.columns = [el + " rank" for el in df_diff_rep_rank.columns]
    df_diff_rep_rank["Rank range"] = df_diff_rep_rank.apply(
        lambda x: max(x) - min(x), axis=1
    )
    df_diff_rep_rank["Common lists"] = df_diff_rep_rank[
        [el for el in df_diff_rep_rank.columns if el.endswith(" rank")]
    ].apply(lambda x: sum([el <= size + rank_tolerance for el in x]), axis=1)

    # Calculate percentiles in local neighborhood
    area_members = dists.index[0:perc_area]
    area_members = df.loc[area_members, :]
    area_dists = pw.distance.pdist(np.array(area_members), "cityblock")
    qsteps = [0, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1]
    q = pd.Series(area_dists).quantile(q=qsteps)

    # Calculate z-scoring of distances to query
    z_dists_in = dists[1:perc_area]
    median = np.median(z_dists_in)
    dists_from_median = np.abs(z_dists_in - median)
    mad = np.median(dists_from_median)
    rob_sd = 1.4826 * mad
    z_dists = -(z_dists_in - median) / rob_sd
    z_dists = pd.DataFrame(z_dists)
    z_dists.columns = ["z-scored closeness"]
    z_dists["z-scoring based category"] = pd.cut(
        z_dists["z-scored closeness"],
        bins=[-1000, 1.28, 1.64, 2.33, 3.09, 1000],
        labels=["-", "B", "*", "**", "***"],
    )

    # assemble output table
    out_ids = dists.index[0:size]
    out = pd.DataFrame(dists[out_ids])
    out = out.join(df_diff_rep, how="left")
    out = out.join(df_diff_rep_rank, how="left")
    out = out.join(
        pd.DataFrame.from_dict(
            {el: meta_dict["id_meta"][el] for el in out_ids}, orient="index"
        ),
        how="left",
    )
    out = out.join(df_core, how="left")
    out["Local distance percentile"] = pd.cut(
        out["Distance measure"], bins=[-1] + list(q.values)[1:], labels=qsteps[1:]
    )
    out = out.join(z_dists, how="left")
    out = out.sort_values("Distance measure").reset_index()

    msg = "Neighborhood recalculated for {}.\n".format(gene) + msg_append

    return out, q, msg


## query multiple genes
## query multiple genes
def multi_query(
    df_core,
    meta_dict,
    genes,
    dist="Manhattan",
    min_corr=0,
    size=50,
    rank_tolerance=25,
    perc_area=250,
    pb=None,
):
    hits = []
    qs = []
    msgs = ""
    for gene in genes:
        output = hv_query_df(
            df_core,
            meta_dict,
            gene=gene,
            size=size,
            rank_tolerance=rank_tolerance,
            min_corr=min_corr,
            perc_area=perc_area,
        )
        if len(output) == 3:
            (h, q, m) = output
            h["Query"] = np.repeat(gene, len(h))
            hits.append(h)
            qs.append(q)
        else:
            m = output
        msgs = msgs + m + "\n"
        if pb:
            pb.value = pb.value + 1

    query_result = pd.concat(hits)
    q = pd.DataFrame(qs).apply(np.mean, axis=0)

    return query_result.reset_index(drop=True), q, msgs


## barplot function single query
def sq_barplot(df, var):
    q_m = df[
        [
            "Distance measure",
            "Rank range",
            "Min correl",
            "Local distance percentile",
            "Gene names",
            "Common lists",
            "z-scoring based category",
        ]
    ]
    q_m["Gene names"] = [
        str(n) + ", rank " + str(r) for n, r in zip(q_m["Gene names"], q_m.index)
    ]
    label_var_dict = {
        "Distance": "Distance measure",
        "z-scoring based category": "z-scoring based category",
        "Rank range between replicates": "Rank range",
        "Worst replicate correlation": "Min correl",
        "Local distance percentile": "Local distance percentile",
    }
    q_m = q_m.iloc[::-1, :]
    var = label_var_dict[var]
    if var == "z-scoring based category":
        z_dict = {"-": 5, "B": 4, "*": 3, "**": 2, "***": 1, np.nan: 0, None: 0}
        q_m["z-scoring based category"] = [
            z_dict[k] for k in q_m["z-scoring based category"]
        ]
        z_dict = {"-": 5, "B": 4, "*": 3, "**": 2, "***": 1, "query": 0}
        p = q_m.hvplot.barh(
            x="Gene names",
            y=var,
            color="Common lists",
            height=60 + 11 * len(q_m),
            xticks=[(z_dict[k], k) for k in z_dict.keys()],
        )
    else:
        p = q_m.hvplot.barh(
            x="Gene names", y=var, color="Common lists", height=60 + 11 * len(q_m)
        )
    p.opts(
        cmap=["orange", "darkgrey", "white", "lightskyblue"][
            0 : len(np.unique(q_m["Common lists"]))
        ][::-1]
    )
    p.opts(legend_position="top")
    return p


## Network construction functions
def filter_nwk_members(query_result, min_rep, max_z):
    nwk_members_1 = query_result[query_result["Common lists"] >= min_rep][
        "Protein IDs"
    ].values
    z_dict = {"-": 5, "B": 4, "*": 3, "**": 2, "***": 1, np.nan: 0, None: 0}
    query_result["z-scoring based category"] = [
        z_dict[k] for k in query_result["z-scoring based category"]
    ]
    nwk_members_2 = query_result[query_result["z-scoring based category"] <= max_z][
        "Protein IDs"
    ].values
    nwk_members = [el for el in nwk_members_1 if el in nwk_members_2]
    nwk_members = np.unique(nwk_members)
    if len(nwk_members) == 1:
        return None
    else:
        return nwk_members


def layout_network(
    df_core, meta_dict, query_result, nwk_members, max_q, q, layout_method="Spring"
):
    df = df_core.reset_index(["Min correl"], drop=True)
    nwk_profiles = df.loc[nwk_members, :]
    dists = pw.distance.pdist(np.array(nwk_profiles), "cityblock")
    dists_pd = pd.DataFrame(
        [
            [
                meta_dict["id_meta"][nwk_members[i1]]["Gene names"],
                meta_dict["id_meta"][nwk_members[i2]]["Gene names"],
            ]
            for i1 in range(len(nwk_members))
            for i2 in range(len(nwk_members))
            if i1 < i2
        ]
    )
    dists_pd.columns = ["Gene 1", "Gene 2"]
    dists_pd["dist"] = [1 / el for el in dists]
    dists_pd["quantile"] = pd.cut(
        pd.Series(dists),
        bins=[-1] + list(q.values)[1:],
        labels=[float(el) for el in q.index[1:]],
    )
    dists_pd = dists_pd[dists_pd["quantile"] <= float(max_q)]
    if len(dists_pd) == 0:
        raise ValueError("No data at this quantile.")

    # Genrate network graph
    G = nx.Graph()
    for idx in range(len(dists_pd)):
        G.add_edge(
            dists_pd.iloc[idx, 0],
            dists_pd.iloc[idx, 1],
            weight=dists_pd.iloc[idx, 2],
            quantile=dists_pd.iloc[idx, 3],
        )
    if layout_method == "Spectral":
        GP = nx.spectral_layout(G)
    elif layout_method == "Spring":
        GP = nx.spring_layout(G, k=5 / np.sqrt(len(dists_pd)))

    return dists_pd, GP


## Functions for drawing a network
def draw_network_figure(dists_pd, GP, query_result, highlight, q):
    # Regerate network graph from distances
    G = nx.Graph()
    for idx in range(len(dists_pd)):
        G.add_edge(
            dists_pd.iloc[idx, 0],
            dists_pd.iloc[idx, 1],
            weight=dists_pd.iloc[idx, 2],
            quantile=dists_pd.iloc[idx, 3],
        )

    fig = plt.figure(figsize=(7, 7))
    fig.set_dpi(300)

    # edges
    edge_styles = [
        {"width": 5, "edge_color": "darkred", "style": "solid", "alpha": 1},  # 0.1%
        {"width": 5, "edge_color": "darkred", "style": "solid", "alpha": 1},  # 0.5%
        {"width": 5, "edge_color": "darkred", "style": "solid", "alpha": 1},  # 1%
        {"width": 4, "edge_color": "red", "style": "solid", "alpha": 1},  # 2.5%
        {"width": 4, "edge_color": "red", "style": "solid", "alpha": 1},  # 5%
        {"width": 3, "edge_color": "darkorange", "style": "solid", "alpha": 1},  # 10%
        {"width": 2, "edge_color": "#ababab", "style": "solid", "alpha": 1},  # 25%
        {"width": 2, "edge_color": "lightgrey", "style": "solid", "alpha": 1},  # 50%
        {"width": 2, "edge_color": "lightgrey", "style": "dotted", "alpha": 1},
    ]  # 100%
    edges = []
    for q_cat, style in zip(q.index[-2::-1], edge_styles[-2::-1]):
        edgelist = [
            (u, v) for (u, v, d) in G.edges(data=True) if d["quantile"] == float(q_cat)
        ]
        if len(edgelist) == 0:
            continue
        Gi = nx.Graph(edgelist)
        nx.draw_networkx_edges(Gi, GP, **style)

    # nodes
    # highlighted gene
    if highlight is not None:
        Gi = nx.Graph([(hit, hit) for hit in [highlight] if hit in G.nodes()])
        nx.draw_networkx_nodes(Gi, GP, node_color="blue", node_size=800)
    # replicate 1
    Gi = nx.Graph(
        [
            (hit, hit)
            for hit, rep in zip(
                query_result["Gene names"], query_result["Common lists"]
            )
            if rep == 1 and hit in G.nodes()
        ]
    )
    nx.draw_networkx_nodes(Gi, GP, node_color="white", node_size=600)
    nx.draw_networkx_labels(Gi, GP, font_size=7)
    # replicate 2
    Gi = nx.Graph(
        [
            (hit, hit)
            for hit, rep in zip(
                query_result["Gene names"], query_result["Common lists"]
            )
            if rep == 2 and hit in G.nodes()
        ]
    )
    nx.draw_networkx_nodes(Gi, GP, node_color="lightgrey", node_size=600)
    nx.draw_networkx_labels(Gi, GP, font_size=7)
    # replicate 3
    Gi = nx.Graph(
        [
            (hit, hit)
            for hit, rep in zip(
                query_result["Gene names"], query_result["Common lists"]
            )
            if rep == 3 and hit in G.nodes()
        ]
    )
    nx.draw_networkx_nodes(Gi, GP, node_color="orange", node_size=600)
    nx.draw_networkx_labels(Gi, GP, font_size=7)
    # query
    if "Query" in query_result.columns:
        Gi = nx.Graph(
            [(hit, hit) for hit in np.unique(query_result["Query"]) if hit in G.nodes()]
        )
        nx.draw_networkx_nodes(Gi, GP, node_color="red", node_size=600)
    else:
        Gi = nx.Graph(
            [(hit, hit) for hit in [query_result["Gene names"][0]] if hit in G.nodes()]
        )
        nx.draw_networkx_nodes(Gi, GP, node_color="red", node_size=600)
    plt.close()

    return fig


def draw_network_interactive(dists_pd, GP, query_result, highlight, q):
    G = nx.Graph()
    for idx in range(len(dists_pd)):
        G.add_edge(
            dists_pd.iloc[idx, 0],
            dists_pd.iloc[idx, 1],
            weight=dists_pd.iloc[idx, 2],
            quantile=dists_pd.iloc[idx, 3],
        )

    width = min(900, int(np.sqrt(len(G.nodes))) * 150)
    if width < 400:
        width = 400

    Gi = nx.Graph([(hit, hit) for hit in [highlight] if hit in G.nodes()])
    highlight = hvnx.draw(
        Gi,
        GP,
        node_color="blue",
        node_size=3000,
        with_labels=False,
        width=width,
        height=width,
    )
    Gi = nx.Graph(
        [
            (hit, hit)
            for hit, rep in zip(
                query_result["Gene names"], query_result["Common lists"]
            )
            if rep == 1 and hit in G.nodes()
        ]
    )
    labels0 = hvnx.draw(
        Gi,
        GP,
        node_color="white",
        node_size=2500,
        with_labels=True,
        width=width,
        height=width,
    )
    Gi = nx.Graph(
        [
            (hit, hit)
            for hit, rep in zip(
                query_result["Gene names"], query_result["Common lists"]
            )
            if rep == 2 and hit in G.nodes()
        ]
    )
    labels1 = hvnx.draw(
        Gi,
        GP,
        node_color="lightgrey",
        node_size=2500,
        with_labels=True,
        width=width,
        height=width,
    )
    Gi = nx.Graph(
        [
            (hit, hit)
            for hit, rep in zip(
                query_result["Gene names"], query_result["Common lists"]
            )
            if rep == 3 and hit in G.nodes()
        ]
    )
    labels2 = hvnx.draw(
        Gi,
        GP,
        node_color="orange",
        node_size=2500,
        with_labels=True,
        width=width,
        height=width,
    )
    if "Query" in query_result.columns:
        Gi = nx.Graph(
            [(hit, hit) for hit in np.unique(query_result["Query"]) if hit in G.nodes()]
        )
        labels3 = hvnx.draw(
            Gi,
            GP,
            node_color="red",
            node_size=2500,
            with_labels=True,
            width=width,
            height=width,
        )
    else:
        Gi = nx.Graph(
            [(hit, hit) for hit in [query_result["Gene names"][0]] if hit in G.nodes()]
        )
        labels3 = hvnx.draw(
            Gi,
            GP,
            node_color="red",
            node_size=2500,
            with_labels=True,
            width=width,
            height=width,
        )
    # .opts(tools=[HoverTool(tooltips=[('index', '@index_hover')])])

    # edges
    edge_styles = [
        {
            "edge_width": 8,
            "edge_color": "darkred",
            "style": "solid",
            "alpha": 1,
        },  # 0.1%
        {
            "edge_width": 8,
            "edge_color": "darkred",
            "style": "solid",
            "alpha": 1,
        },  # 0.5%
        {"edge_width": 8, "edge_color": "darkred", "style": "solid", "alpha": 1},  # 1%
        {"edge_width": 6, "edge_color": "red", "style": "solid", "alpha": 1},  # 2.5%
        {"edge_width": 6, "edge_color": "red", "style": "solid", "alpha": 1},  # 5%
        {
            "edge_width": 4,
            "edge_color": "darkorange",
            "style": "solid",
            "alpha": 1,
        },  # 10%
        {"edge_width": 3, "edge_color": "#ababab", "style": "solid", "alpha": 1},  # 25%
        {
            "edge_width": 2,
            "edge_color": "lightgrey",
            "style": "solid",
            "alpha": 1,
        },  # 50%
        {"edge_width": 2, "edge_color": "lightgrey", "style": "dotted", "alpha": 1},
    ]  # 100%
    edges = []
    for q_cat, style in zip(q.index[1:], edge_styles):
        edgelist = [
            (u, v) for (u, v, d) in G.edges(data=True) if d["quantile"] == float(q_cat)
        ]
        if len(edgelist) == 0:
            continue
        Gi = nx.Graph(edgelist)
        edge = hvnx.draw(
            Gi,
            GP,
            node_size=2000,
            with_labels=False,
            width=width,
            height=width,
            **style,
        )
        edges.append(edge)

    nw = hv.Overlay(edges[::-1]) * highlight * labels0 * labels1 * labels2 * labels3
    return nw
