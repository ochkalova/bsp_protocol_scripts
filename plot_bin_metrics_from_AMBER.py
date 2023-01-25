#!/usr/bin/env python

import argparse
from collections import defaultdict
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os

def read_amber_tables(input_folder, methods):
    pooled_data = defaultdict(list)
    n_methods = 0
    input_folder = os.path.join(input_folder, "genome")
    for method in methods:
        method_dir = os.path.join(input_folder, method)
        input_file = os.path.join(method_dir, "metrics_per_bin.tsv")
        metrics = pd.read_csv(input_file, sep='\t')
        metrics['Purity (bp)'] = metrics['Purity (bp)'] * 100
        metrics['Completeness (bp)'] = metrics['Completeness (bp)'] * 100
        metrics['F1'] = 2 * (metrics['Purity (bp)'] * metrics['Completeness (bp)']) / (
            metrics['Purity (bp)'] + metrics['Completeness (bp)'])

        pooled_data['Completeness'].extend(metrics['Completeness (bp)'].to_list())
        pooled_data['Purity'].extend(metrics['Purity (bp)'].to_list())
        pooled_data['F1'].extend(metrics['F1'].to_list())
        pooled_data['binner'].extend([method]*metrics['F1'].size)
        n_methods += 1
    return pooled_data, n_methods


def create_swarmplot(pooled_data, n_colors, output_folder, height, width):
    df = pd.DataFrame(pooled_data)
    for metric in ["Completeness", "Purity", "F1"]:
        plt.figure(figsize=(width, height))
        sns.color_palette("Set2")
        ax = sns.swarmplot(x=metric, y="binner", data=df[[metric, 'binner']], size=6,
                           palette=sns.color_palette("hls", n_colors))
        ax = sns.boxplot(x=metric, y="binner", data=df[[metric, 'binner']], linewidth=0.7,
                         palette=sns.color_palette("hls", n_colors))
        for patch in ax.artists:
            r, g, b, a = patch.get_facecolor()
            patch.set_facecolor((r, g, b, .2))
        ax.tick_params(labelsize=18)
        ax.set_ylabel(None)
        ax.set_xlim(-5, 105)
        ax.invert_xaxis()
        ax.set_xlabel("{} (%)".format(metric), fontsize=25)
        
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
            
        fig_path = os.path.join(output_folder, "AMBER_{}_swarmplot.png".format(metric))
        fig = ax.get_figure()
        fig.savefig(fig_path, dpi=300, bbox_inches='tight')


def main():
    parser = argparse.ArgumentParser(
        description="Create 3 swarmplots with Purity, Completeness and F1 with per bin from input AMBER tsv tables")
    parser.add_argument("-i", "--input_folder", required=True, help="AMBER results folder")
    parser.add_argument("-o", "--output_folder", required=True, help="Output directory to save generated plots")
    parser.add_argument("-l", "--labels", required=True, metavar="Binning_method_1,Binning_method_2,...",
                        help="Labels that were listed in AMBER run, comma-separated, without spaces")
    parser.add_argument("-W", "--width", type=int, metavar="N", default=6, required=False, help="Output picture width (default: 6)")
    parser.add_argument("-H", "--height", type=int, metavar="N", default=14, required=False, help="Output picture height (default: 14)")
    args = parser.parse_args()
    
    binnings_metrics, n_methods = read_amber_tables(args.input_folder, args.labels.split(","))
    create_swarmplot(binnings_metrics, n_methods, args.output_folder, args.height, args.width)


if __name__ == "__main__":
    main()
    


