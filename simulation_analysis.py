"""
Script to plot SBS distribution over 96 types for Tracksig-generated simulation
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.ticker as mtick
from matplotlib import pyplot as plt


def plot_mut_spec(mut_dataframe, title: str, plot_prefix: str,
                  plot_dir: str,
                  y_lim=None, y_label='Percentage of SBS',):
    """
    Method to plot mutation spectrum from counts/frequency input.
    Inputs: 
        mut_dataframe: Counts or frequencies of mutations across 96 types (shape=(96,))
        y_lim: upper limit of y axis scale
        y_label: label for y axis
        title: plot title
    """
    # check if mut_dataframe vec sum to 1
    print(mut_dataframe.to_numpy().sum())
    # get 96 types
    sigs = pd.read_csv('annotation/alexSignatures_w_header.csv')
    mut_types_96 = sigs.iloc[:, 0]

    # get 6 types and shorten 96 types
    mut_types_6, mut_types_96_short = [], []
    for entry in mut_types_96:
        new_entry = entry[0]+'>'+entry[2]
        mut_types_6.append(new_entry)
        mut_types_96_short.append(entry[4:])

    # make color group num
    color_groups = mut_dataframe.index // 16

    plt.figure(figsize=(15, 5))
    # colors = ['red', 'blue', 'green', 'orange', 'purple', 'cyan']
    colors = ['deepskyblue', 'black', 'red', 'lightgray', 'yellowgreen', 'pink']

    for i in range(96):
        plt.bar(mut_types_96[i],
                mut_dataframe.iloc[i] * 100,  # convert to percentage
                color=colors[color_groups[i]],
                )
    plt.xlabel('96 types')
    plt.ylabel(y_label)
    plt.title(title)
    # remove white spacing on left and right of plot areas
    plt.margins(x=0.005)
    plt.xticks(mut_types_96, mut_types_96_short, rotation=75, fontsize=8)
    
    if y_lim is not None:
        plt.ylim(0, y_lim * 100)  # convert to percentage
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter())

    # Figure legend (grouping the 96 types into 6 main groups)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[i]) for i in range(6)]
    # labels = set(mut_types_6)  # get unique item from list --> set is NOT ordered
    labels = list(dict.fromkeys(mut_types_6))
    plt.legend(handles, labels, title='Groups')

    # savefig
    plt.savefig(os.path.join(
        plot_dir, f'{plot_prefix}.png'), dpi=150, bbox_inches='tight')


def plot_all_three(type_to_plot: str, filename: str, data: tuple, plot_dir: str,
                   y_lim, bin_num: int = None, bin_range: tuple = None):
    """Method to plot trios fo plots for sampled data, and simulated vs inferred freq."""

    if type_to_plot == 'sum':  # summing up all bins
        num_bins = 50  # 100 mutations per bin for 50 time points bin
        to_plot = [data_i.sum(axis=1) / num_bins for data_i in data]
        plot_suffix = '50_bins'
    elif type_to_plot == 'bin' and bin_num is not None:
        # already a single bin so no need to normalize
        to_plot = [data_i.iloc[:, bin_num] for data_i in data]
        plot_suffix = f'bin_{bin_num}'
    elif type_to_plot == 'range' and bin_range is not None:
        num_bins = bin_range[1] - bin_range[0]
        to_plot = [data_i.iloc[:, bin_range[0]:bin_range[1]] /
                   num_bins for data_i in data]
        plot_suffix = f'bins_{bin_range[0]}_{bin_range[1]}'
    else:
        print("Plotting range not specified")

    plot_mut_spec(to_plot[0], f'Simulated data mutation spectrum {filename} {plot_suffix}',
                  f'{filename}_{plot_suffix}_sample', plot_dir,  y_lim)
    plot_mut_spec(to_plot[1], f'Simulated loading mutation spectrum {filename} {plot_suffix}',
                  f'{filename}_{plot_suffix}_sim_loading', plot_dir,  y_lim)
    plot_mut_spec(to_plot[2], f'Inferred loading mutation spectrum {filename} {plot_suffix}',
                  f'{filename}_{plot_suffix}_inferred_loading', plot_dir,  y_lim)


if __name__ == "__main__":
    sample_name = sys.argv[1]
    # inputs from command line
    # type of plots: sum, range, bin
    # if range or bin need to specify range or bin num

    y_lim = float(sys.argv[2])

    # get the two signatures chosen (for alex table filtering later)
    sig1, sig2 = sample_name.split('_')[0], sample_name.split('_')[1]

    # make results dir
    result_dir = 'sim_analysis_plots'
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    # load files
    sim_sample = pd.read_csv(f'simulated_data/{sample_name}.csv', header=None)

    cosmic_sigs = pd.read_csv('annotation/alexSignatures_w_header.csv')
    cosmic_sigs_filter = cosmic_sigs.loc[:, ['S1', 'S5', sig1, sig2]]

    sim_exposure = pd.read_csv(
        f'simulated_data/{sample_name}.exposure.csv', header=None)
    inferred_exposure = pd.read_csv(
        f'results_signature_trajectories/{sample_name}/{sample_name}/mixtures.csv')

    # extract just the counts from data (ignore labels)
    cosmic_sigs_sample = cosmic_sigs_filter.to_numpy()
    sim_exposure_arr = sim_exposure.iloc[:, 1:].to_numpy()
    inferred_exposure_arr = np.c_[
        np.zeros(4), inferred_exposure.iloc[:, 1:].to_numpy()]

    # matrix multiplication
    sim_data = cosmic_sigs_sample @ sim_exposure_arr
    inferred_data = cosmic_sigs_sample @ inferred_exposure_arr

    data_to_plot = sim_sample.iloc[:, 1:] / \
        100, pd.DataFrame(sim_data), pd.DataFrame(inferred_data)
    # divide sim_sample (raw snp scounts) by 100 muts to normalize for each bin
    # to sum to 1, matching the other two dataframes

    # plot various scenarios
    bin_num, bin_range = None, None
    if sys.argv[3] == 'bin' and sys.argv[4] is not None:
        bin_num = int(sys.argv[4])
    if sys.argv[3] == 'range' and sys.argv[4] is not None:
        bin_range = int(sys.argv[4]), int(sys.argv[5])
    plot_all_three(sys.argv[3], sys.argv[1], data_to_plot, result_dir,
                   y_lim, bin_num, bin_range)
