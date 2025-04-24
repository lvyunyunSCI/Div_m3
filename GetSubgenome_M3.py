#!/usr/bin/env python3
"""
Mash Chromosome Comparison and Visualization Tool (Polyploid Version, need seqkit,mash,python)

2025/04/23 feautures:
- Supports any ploidy level (diploid, triploid, tetraploid, etc.)
- Configurable number of subgenomes
- Dynamic color scheme generation
- Enhanced visualization for multiple subgenomes

2025/04/23 feautures:   
1. Split reference and query genomes
2. Build Mash reference database
3. Calculate inter-chromosomal distances
4. Process result data
5. Generate visualization plots 
    
Subcommands:
- all: Run complete pipeline
- calculate: Only perform calculations without visualization 
- plot: Only perform visualization using pre-calculated results
    
Dependencies:
- seqkit (must be pre-installed and in PATH)
- mash (must be pre-installed and in PATH)
- Python packages: pandas, matplotlib, seaborn
    
Author:Yunyun Lv;
Date:2025-4-24
"""

import os
import sys
import re
import argparse
import subprocess
import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import CubicSpline
from itertools import cycle
from colorsys import hsv_to_rgb

matplotlib.use('Agg')
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

# Default configuration
DEFAULT_SUBGENOMES = 2
DEFAULT_THREADS = 4
MAX_SUBGENOMES = 10  # Practical upper limit for visualization

def split_fasta(input_fasta, output_dir, prefix):
    """Split FASTA file using seqkit"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    cmd = [
        "seqkit", "split",
        "-f", "-i",
        "--by-id-prefix", "",
        "--out-dir", output_dir,
        input_fasta
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        print(f"[INFO] Successfully split {input_fasta} to {output_dir}")
    except subprocess.CalledProcessError as e:
        sys.exit(f"[ERROR] seqkit split failed: {e}")
def build_mash_db(chr_list, db_name, threads):
    """Build Mash database"""
    cmd = [
        "mash", "sketch",
        "-p", str(threads),
        "-k", "31",
        "-s", "5000000000",
        "-l", chr_list,
        "-o", db_name
    ]
    try:
        subprocess.run(cmd, check=True)
        print(f"[INFO] Successfully built Mash database: {db_name}.msh")
    except subprocess.CalledProcessError as e:
        sys.exit(f"[ERROR] Mash database construction failed: {e}")

def calculate_distances(db_msh, query_list, output_file, threads):
    """Calculate Mash distances"""
    cmd = [
        "mash", "dist",
        db_msh,
        "-p", str(threads),
        "-s", "5000000000",
        "-l", query_list
    ]
    try:
        with open(output_file, "w") as f:
            subprocess.run(cmd, check=True, stdout=f)
        print(f"[INFO] Distance calculation completed: {output_file}")
    except subprocess.CalledProcessError as e:
        sys.exit(f"[ERROR] Distance calculation failed: {e}")

def generate_color_scheme(n_subgenomes):
    """Generate a distinct color scheme for N subgenomes"""
    colors = []
    for i in range(n_subgenomes):
        hue = i / n_subgenomes  # Evenly distribute hues
        sat = 0.7 + 0.3 * (i % 2)  # Alternate saturation
        val = 0.8 - 0.2 * (i % 3)  # Vary brightness
        r, g, b = hsv_to_rgb(hue, sat, val)
        main_color = "#{:02x}{:02x}{:02x}".format(int(r*255), int(g*255), int(b*255))
        dark_color = "#{:02x}{:02x}{:02x}".format(int(r*200), int(g*200), int(b*200))
        colors.append({
            "main": main_color,
            "secondary": dark_color,
            "marker": ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h'][i % 10]
        })
    return {f"SG{i+1}": colors[i] for i in range(n_subgenomes)}

def process_results(input_file, output_prefix, n_subgenomes=DEFAULT_SUBGENOMES):
    """Process result data with configurable subgenome count"""
    # Read raw data
    df = pd.read_csv(input_file, sep="\t", header=None, 
                    names=["Reference", "Query", "Distance", "PValue", "MatchingHashes"])
    
    # Process filenames
    df["RefChr"] = df["Reference"].apply(lambda x: Path(x).stem)
    df["QryChr"] = df["Query"].apply(lambda x: Path(x).stem)
    
    # Filter top N closest results for each reference chromosome
    result = []
    for (ref_chr), group in df.groupby("RefChr"):
        top_n = group.nsmallest(n_subgenomes, "Distance")
        for i, (idx, row) in enumerate(top_n.iterrows(), 1):
            result.append({
                "Rchr": row["RefChr"],
                "Qchr": row["QryChr"],
                "Subg": f"SG{i}",
                "MashD": row["Distance"]
            })
    
    # Generate final data
    final_df = pd.DataFrame(result)
    final_df = final_df.sort_values(by=["Rchr", "MashD"])
    
    # Save results
    filter_file = f"{output_prefix}.filter.Gadd"
    final_df.to_csv(filter_file, sep="\t", index=False)
    print(f"[INFO] Processed results with {n_subgenomes} subgenomes saved to: {filter_file}")
    return final_df

def visualize(data, output_pdf,n_subgenomes=DEFAULT_SUBGENOMES):
    """Enhanced visualization supporting multiple subgenomes"""
    # Generate dynamic color scheme
    COLOR_SCHEME = generate_color_scheme(n_subgenomes)
    
    plt.figure(figsize=(20, 12))
    sns.set_theme(style="whitegrid", font_scale=1.2)
    
    # Data preprocessing
    data = data.copy()
    data["Rchr"] = data["Rchr"].astype(str)
    data["Qchr"] = data["Qchr"].astype(str)
    
    def natural_sort_key(s):
        s = str(s)
        return tuple(int(text) if text.isdigit() else text.lower()
                    for text in re.split('([0-9]+)', s))
    
    data["SortKey"] = data["Rchr"].apply(natural_sort_key)
    data = data.sort_values(["SortKey", "Subg"]).reset_index(drop=True)
    chr_order = sorted(data["Rchr"].unique(), key=natural_sort_key)
    chr_to_xpos = {chr_name: i for i, chr_name in enumerate(chr_order)}
    
    fig, ax = plt.subplots(figsize=(20, 12 + n_subgenomes))  # Adjust height based on ploidy
    
    # Plot trend lines for each subgenome
    for group in [f"SG{i+1}" for i in range(n_subgenomes)]:
        group_data = data[data["Subg"] == group]
        if len(group_data) > 1:
            x = [chr_to_xpos[chr_name] for chr_name in group_data["Rchr"]]
            y = group_data["MashD"].values
            color = COLOR_SCHEME[group]["main"] + "80"
            
            cs = CubicSpline(x, y)
            x_smooth = np.linspace(min(x), max(x), 300)
            ax.plot(
                x_smooth, cs(x_smooth),
                color=color,
                linestyle="-" if int(group[2:]) % 2 else "--",
                linewidth=3,
                alpha=0.5,
                zorder=1
            )
    
    # Connection lines between subgenomes
    for chr_name in chr_order:
        chr_data = data[data["Rchr"] == chr_name]
        if len(chr_data) > 1:
            y_values = chr_data["MashD"].tolist()
            x_pos = chr_to_xpos[chr_name]
            ax.plot(
                [x_pos, x_pos],
                [min(y_values), max(y_values)],
                color="#666666",
                linestyle=":",
                alpha=0.5,
                lw=1.0,
                zorder=2
            )
    
    # Scatter plot with dynamic markers
    plot_data = data.copy()
    plot_data["x_pos"] = plot_data["Rchr"].map(chr_to_xpos)
    
    markers = {f"SG{i+1}": COLOR_SCHEME[f"SG{i+1}"]["marker"] for i in range(n_subgenomes)}
    
    scatter = sns.scatterplot(
        data=plot_data,
        x="x_pos",
        y="MashD",
        hue="Subg",
        style="Subg",
        markers=markers,
        s=200,
        edgecolor=[COLOR_SCHEME[g]["secondary"] for g in plot_data["Subg"]],
        linewidth=2.0,
        palette={g: COLOR_SCHEME[g]["main"] for g in [f"SG{i+1}" for i in range(n_subgenomes)]},
        ax=ax,
        zorder=3,
        legend=False
    )
    
    # Data labels
    for index, row in plot_data.iterrows():
        color = COLOR_SCHEME[row["Subg"]]["main"]
        ax.annotate(
            text=row["Qchr"],
            xy=(row["x_pos"], row["MashD"]),
            xytext=(0, 15 + 5 * int(row["Subg"][2:])),  # Stagger labels
            textcoords="offset points",
            ha="center",
            va="bottom",
            rotation=45,
            fontsize=10,
            color=color,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor=color, alpha=0.8)
        )
    
    # Chart decoration
    ref_prefix = data['Rchr'].iloc[0].split('_')[0] if '_' in data['Rchr'].iloc[0] else data['Rchr'].iloc[0]
    qry_prefix = data['Qchr'].iloc[0].split('_')[0] if '_' in data['Qchr'].iloc[0] else data['Qchr'].iloc[0]
    
    ax.set_title(
        f"Chromosome Comparison ({n_subgenomes} Subgenomes)\n"
        f"Reference | Query",
        fontsize=16,
        pad=25
    )
    
    ax.set_xlabel("Reference Chromosome", labelpad=15)
    ax.set_ylabel("Mash Distance", labelpad=15)
    ax.set_xticks(np.arange(len(chr_order)))
    ax.set_xticklabels(chr_order, rotation=45, ha="right", fontsize=11)
    plt.yticks(fontsize=11)
    ax.yaxis.grid(True, linestyle="--", alpha=0.3)
    
    # Dynamic legend
    legend_elements = [
        Line2D([0], [0], 
               marker=COLOR_SCHEME[f"SG{i+1}"]["marker"],
               color='w',
               markerfacecolor=COLOR_SCHEME[f"SG{i+1}"]["main"],
               markeredgecolor=COLOR_SCHEME[f"SG{i+1}"]["secondary"],
               markersize=12,
               label=f'SG{i+1}') 
        for i in range(n_subgenomes)
    ]
    
    ax.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        title=f"Subgenomes (n={n_subgenomes})",
        frameon=True,
        title_fontsize=12
    )
    
    plt.tight_layout()
    plt.subplots_adjust(right=0.85 - 0.02*n_subgenomes)  # Adjust for legend width
    plt.savefig(output_pdf, dpi=300, bbox_inches="tight")
    print(f"[INFO] Visualization with {n_subgenomes} subgenomes saved to: {output_pdf}")
    plt.close()

def run_calculation(args):
    """Execute calculation pipeline with configurable subgenomes"""
    split_fasta(args.ref_fasta, f"{args.ref_abb}_split", "")
    split_fasta(args.qry_fasta, f"{args.qry_abb}_split", "")
    
    ref_chrlist = f"{args.ref_abb}.chrList"
    with open(ref_chrlist, "w") as f:
        for path in Path(f"{args.ref_abb}_split").glob("*"):
            f.write(f"{path.resolve()}\n")
    
    build_mash_db(ref_chrlist, args.ref_abb, args.threads)
    
    qry_chrlist = f"{args.qry_abb}.chrList"
    with open(qry_chrlist, "w") as f:
        for path in Path(f"{args.qry_abb}_split").glob("*"):
            f.write(f"{path.resolve()}\n")
    
    dist_file = f"{args.ref_abb}_{args.qry_abb}_mashDistance"
    calculate_distances(f"{args.ref_abb}.msh", qry_chrlist, dist_file, args.threads)
    
    return process_results(dist_file, f"{args.ref_abb}_{args.qry_abb}_mashDistance", args.subgenomes)

def run_visualization(args):
    """Execute visualization with configurable subgenomes"""
    try:
        processed_df = pd.read_csv(args.data_file, sep="\t")
        # Infer number of subgenomes from data if not specified
        n_subgenomes = args.subgenomes if args.subgenomes else processed_df["Subg"].nunique()
    except Exception as e:
        sys.exit(f"[ERROR] Unable to read data file: {e}") 
    output_pdf = args.output_pdf or f"{Path(args.data_file).stem}.pdf"
    visualize(processed_df,output_pdf, n_subgenomes)
def main():
    parser = argparse.ArgumentParser(
        description="Mash Chromosome Comparison Tool for Polyploid Genomes",
        formatter_class=argparse.RawTextHelpFormatter)
    
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    # Shared arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("-s", "--subgenomes", type=int, default=DEFAULT_SUBGENOMES,
                             choices=range(1, MAX_SUBGENOMES+1),
                             help=f"Number of subgenomes to consider (1-{MAX_SUBGENOMES}, default: {DEFAULT_SUBGENOMES})")
    parent_parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS,
                             help=f"Number of threads (default: {DEFAULT_THREADS})")
    
    # all command
    parser_all = subparsers.add_parser('all', parents=[parent_parser],
                                     help='Run complete pipeline')
    parser_all.add_argument("ref_abb", help="Reference genome abbreviation")
    parser_all.add_argument("ref_fasta", help="Path to reference genome FASTA")
    parser_all.add_argument("qry_abb", help="Query genome abbreviation")
    parser_all.add_argument("qry_fasta", help="Path to query genome FASTA")
    parser_all.add_argument("-o", "--output-pdf", help="Custom output PDF filename")
    parser_all.set_defaults(func=run_all)
    
    # calculate command
    parser_calc = subparsers.add_parser('calculate', parents=[parent_parser],
                                       help='Only perform calculations')
    parser_calc.add_argument("ref_abb", help="Reference genome abbreviation")
    parser_calc.add_argument("ref_fasta", help="Path to reference genome FASTA")
    parser_calc.add_argument("qry_abb", help="Query genome abbreviation")
    parser_calc.add_argument("qry_fasta", help="Path to query genome FASTA")
    parser_calc.set_defaults(func=run_calculation)
    
    # plot command
    parser_plot = subparsers.add_parser('plot', help='Only perform visualization')
    parser_plot.add_argument("data_file", help="Result data file (.filter.Gadd)")
    parser_plot.add_argument("-s", "--subgenomes", type=int,
                           help="Number of subgenomes (default: auto-detect from data)")
    parser_plot.add_argument("-o", "--output-pdf", help="Custom output PDF filename")
    parser_plot.set_defaults(func=run_visualization)
    
    args = parser.parse_args()
    args.func(args)

def run_all(args):
    processed_df = run_calculation(args)
    output_pdf = args.output_pdf or f"{args.ref_abb}_{args.qry_abb}_mashDistance.filter.Gadd.pdf"
    visualize(processed_df, output_pdf, args.subgenomes)

if __name__ == "__main__":
    main()
