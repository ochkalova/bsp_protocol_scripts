#!/usr/bin/env python

import argparse
from collections import defaultdict
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os

def read_checkm_tables(opened_file):
    metrics = pd.read_csv(opened_file, sep='\t')
    metrics['Contamination'] = metrics['Contamination'] / 100
    metrics['Completeness'] = metrics['Completeness'] / 100
    # We define purity as Purity = 1 / (1 + CheckM_Contamination)
    metrics['Purity'] =  1 / (1 + metrics['Contamination']) 
    metrics['F1'] = 2 * (metrics['Purity'] * metrics['Completeness']) / (
                metrics['Purity'] + metrics['Completeness']) * 100
    metrics['Completeness'] = metrics['Completeness'] * 100
    metrics['Purity'] = metrics['Purity'] * 100
    return metrics[['Completeness', 'Purity', 'F1']]


def collect_data(files: list, labels: dict):
    pooled_data = defaultdict(list)
    for binning in files:
        label = labels[binning]
        with open(binning) as opened_file:
            binning_data = read_checkm_tables(opened_file)
        pooled_data['Completeness'].extend(binning_data['Completeness'].to_list())
        pooled_data['Purity'].extend(binning_data['Purity'].to_list())
        pooled_data['F1'].extend(binning_data['F1'].to_list())
        pooled_data['binner'].extend([label] * binning_data['F1'].size)
    return pooled_data


def create_swarmplot(pooled_data, files, output_dir, height, width):
    df = pd.DataFrame(pooled_data)
    for metric in ["Completeness", "Purity", "F1"]:
        plt.figure(figsize=(width, height))
        sns.color_palette("Set2")
        n_colors = len(files)
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
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        fig_path = os.path.join(output_dir, "CheckM_{}_swarmplot.png".format(metric))
        fig = ax.get_figure()
        fig.savefig(fig_path, dpi=300, bbox_inches='tight')


def main():
    parser = argparse.ArgumentParser(
        description="Create 3 swarmplots with Purity, Completeness and F1 with per bin from input CheckM tsv tables")
    parser.add_argument("input_files", nargs='+', help="CheckM tab-separated files including full paths")
    parser.add_argument("-l", "--labels", required=True, metavar="Binning_method_1,Binning_method_2,...",
                        help="Names of binning methods to use in plots, comma-separated, without spaces")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory to save generated plots")
    parser.add_argument("-W", "--width", type=int, metavar="N", default=6, required=False, help="Output picture width (default: 6)")
    parser.add_argument("-H", "--height", type=int, metavar="N", default=14, required=False, help="Output picture height (default: 14)")
    args = parser.parse_args()

    files = args.input_files
    label2binning = dict(zip(files, args.labels.split(",")))
    binnings_metrics = collect_data(files, label2binning)
    create_swarmplot(binnings_metrics, files, args.output_dir, args.height, args.width)


if __name__ == "__main__":
    main()
    

